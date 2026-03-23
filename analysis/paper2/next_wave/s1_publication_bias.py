#!/usr/bin/env python3
"""
Study 1: Publication & Excavation Density Bias Test
====================================================
Tests whether on-corridor archaeological sites (near the Alison Great Circle)
have more archaeological attention than off-corridor sites, creating detection bias.

Phases:
  1. Radiocarbon Intensity Test (p3k14c + XRONOS)
  2. Pleiades Attention Proxy
  3. Monument vs Settlement Decomposition
  4. Distance-Decile Analysis
"""

import sys, os, math, json, warnings
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

sys.stdout = os.fdopen(sys.stdout.fileno(), "w", buffering=1)
warnings.filterwarnings("ignore", category=FutureWarning)

# ── Shared Constants ──────────────────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
THRESHOLD_KM = 50
MC_TRIALS = 1000

np.random.seed(42)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = Path("/Users/elliotallan/megalith_site_research")
P3K_PATH = BASE / "data/p3k14c/p3k14c_data.csv"
PLEIADES_PATH = BASE / "data/pleiades/pleiades-places-latest.csv"
XRONOS_PATH = BASE / "data/xronos/xronos_sites.csv"
OUT_DIR = BASE / "outputs/next_wave/publication_bias"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Geometry ──────────────────────────────────────────────────────────────────

def _to_rad(d):
    return d * math.pi / 180.0


def angular_distance_deg(lat1, lon1, lat2, lon2):
    """Great-circle angular distance in degrees."""
    rlat1, rlon1 = _to_rad(lat1), _to_rad(lon1)
    rlat2, rlon2 = _to_rad(lat2), _to_rad(lon2)
    dlon = rlon2 - rlon1
    cos_d = (math.sin(rlat1) * math.sin(rlat2)
             + math.cos(rlat1) * math.cos(rlat2) * math.cos(dlon))
    cos_d = max(-1.0, min(1.0, cos_d))
    return math.degrees(math.acos(cos_d))


def dist_from_gc(lat, lon):
    """
    Minimum distance (km) from a point to the great circle whose pole
    is (POLE_LAT, POLE_LON).  Distance = |angular_dist_to_pole - 90deg|
    converted to km.
    """
    ang = angular_distance_deg(lat, lon, POLE_LAT, POLE_LON)
    offset_deg = abs(ang - 90.0)
    return offset_deg * (math.pi / 180.0) * EARTH_R_KM


def dist_from_gc_vec(lats, lons):
    """Vectorised wrapper."""
    return np.array([dist_from_gc(la, lo) for la, lo in zip(lats, lons)])


# ── Stats helpers ─────────────────────────────────────────────────────────────

def compare_groups(on_vals, off_vals, label=""):
    """Mann-Whitney U, medians, ratio, bootstrap CI."""
    on = np.asarray(on_vals, dtype=float)
    off = np.asarray(off_vals, dtype=float)
    if len(on) < 5 or len(off) < 5:
        return {"label": label, "n_on": int(len(on)), "n_off": int(len(off)),
                "skip": True, "reason": "too few sites"}

    u_stat, u_p = stats.mannwhitneyu(on, off, alternative="two-sided")
    med_on, med_off = float(np.median(on)), float(np.median(off))
    ratio = med_on / med_off if med_off > 0 else float("inf")

    # Bootstrap 95% CI on median ratio
    ratios = []
    for _ in range(MC_TRIALS):
        s_on = np.random.choice(on, size=len(on), replace=True)
        s_off = np.random.choice(off, size=len(off), replace=True)
        m_off = np.median(s_off)
        if m_off > 0:
            ratios.append(np.median(s_on) / m_off)
    if len(ratios) >= 20:
        ci_lo, ci_hi = np.percentile(ratios, [2.5, 97.5])
    else:
        ci_lo, ci_hi = np.nan, np.nan

    return {
        "label": label,
        "n_on": int(len(on)), "n_off": int(len(off)),
        "median_on": round(med_on, 3), "median_off": round(med_off, 3),
        "median_ratio": round(ratio, 3),
        "mann_whitney_U": float(u_stat), "mann_whitney_p": float(u_p),
        "bootstrap_ci_95": [round(float(ci_lo), 3), round(float(ci_hi), 3)],
        "mean_on": round(float(np.mean(on)), 3),
        "mean_off": round(float(np.mean(off)), 3),
    }


# ── Plotting ──────────────────────────────────────────────────────────────────

def _setup_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def save_histogram(on_vals, off_vals, title, fname):
    plt = _setup_mpl()
    fig, ax = plt.subplots(figsize=(8, 5))
    combined = np.concatenate([on_vals, off_vals])
    upper = min(np.percentile(combined, 99), 100)
    bins = np.linspace(0, upper, 40)
    ax.hist(off_vals, bins=bins, alpha=0.55,
            label=f"Off-corridor (n={len(off_vals)})", color="#3b82f6")
    ax.hist(on_vals, bins=bins, alpha=0.55,
            label=f"On-corridor (n={len(on_vals)})", color="#ef4444")
    ax.set_xlabel("Dates per site")
    ax.set_ylabel("Number of sites")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / fname, dpi=150)
    plt.close(fig)
    print(f"  -> saved {fname}")


def save_decile_plot(decile_df, title, fname):
    plt = _setup_mpl()
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(decile_df["decile"], decile_df["mean_dates"],
           color="#6366f1", edgecolor="white")
    ax.set_xlabel("Distance decile (1 = closest to GC)")
    ax.set_ylabel("Mean dates per site")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(OUT_DIR / fname, dpi=150)
    plt.close(fig)
    print(f"  -> saved {fname}")


def _fmt(res, key, fmt_str=".4g"):
    v = res.get(key, "N/A")
    if v == "N/A" or (isinstance(v, float) and np.isnan(v)):
        return "N/A"
    return f"{v:{fmt_str}}"


# ===========================================================================
#  PHASE 1: Radiocarbon Intensity Test
# ===========================================================================
print("=" * 70)
print("PHASE 1: Radiocarbon Intensity Test")
print("=" * 70)

# ── 1a. p3k14c ────────────────────────────────────────────────────────────
print("\n  --- p3k14c ---")
p3k = pd.read_csv(P3K_PATH, low_memory=False)
print(f"  Loaded {len(p3k)} radiocarbon records")

p3k = p3k.dropna(subset=["Lat", "Long", "SiteID"])
print(f"  After dropping missing lat/lon/SiteID: {len(p3k)} records")

p3k_sites = p3k.groupby("SiteID").agg(
    dates_per_site=("LabID", "count"),
    Lat=("Lat", "first"),
    Long=("Long", "first"),
    Continent=("Continent", "first"),
).reset_index()
print(f"  Unique sites: {len(p3k_sites)}")

p3k_sites["dist_km"] = dist_from_gc_vec(
    p3k_sites["Lat"].values, p3k_sites["Long"].values)

p3k_on = p3k_sites[p3k_sites["dist_km"] < THRESHOLD_KM]
p3k_off = p3k_sites[p3k_sites["dist_km"] > 200]
print(f"  On-corridor (<50 km): {len(p3k_on)} sites")
print(f"  Off-corridor (>200 km): {len(p3k_off)} sites")

if len(p3k_on) < 5:
    print("  WARNING: p3k14c has insufficient on-corridor coverage.")
    print("  The great circle passes through regions (North Africa, South America,")
    print("  South Pacific) with very sparse p3k14c sampling. This database is")
    print("  geographically biased toward North America and Europe.")
    phase1_p3k = {
        "label": "p3k14c_overall", "n_on": int(len(p3k_on)),
        "n_off": int(len(p3k_off)), "skip": True,
        "reason": "insufficient on-corridor coverage in p3k14c"
    }
    phase1_p3k_regional = {}
else:
    phase1_p3k = compare_groups(p3k_on["dates_per_site"].values,
                                p3k_off["dates_per_site"].values,
                                label="p3k14c_overall")
    print(f"  Mann-Whitney p = {_fmt(phase1_p3k, 'mann_whitney_p')}")
    print(f"  Median ratio (on/off) = {_fmt(phase1_p3k, 'median_ratio', '.3f')}")
    save_histogram(p3k_on["dates_per_site"].values,
                   p3k_off["dates_per_site"].values,
                   "p3k14c: Dates per Site (On- vs Off-Corridor)",
                   "phase1_p3k_histogram.png")

    # Regional control
    print("\n  Regional breakdown by continent (p3k14c):")
    phase1_p3k_regional = {}
    p3k_on_mask = p3k_sites["dist_km"] < THRESHOLD_KM
    p3k_off_mask = p3k_sites["dist_km"] > 200
    for cont in sorted(p3k_sites["Continent"].dropna().unique()):
        on_c = p3k_sites[p3k_on_mask & (p3k_sites["Continent"] == cont)]
        off_c = p3k_sites[p3k_off_mask & (p3k_sites["Continent"] == cont)]
        res = compare_groups(on_c["dates_per_site"].values,
                             off_c["dates_per_site"].values,
                             label=f"p3k14c_{cont}")
        phase1_p3k_regional[cont] = res
        if res.get("skip"):
            print(f"    {cont}: skipped ({res['reason']}, on={res['n_on']}, off={res['n_off']})")
        else:
            print(f"    {cont}: ratio={res['median_ratio']:.3f}, "
                  f"p={res['mann_whitney_p']:.4g} (on={res['n_on']}, off={res['n_off']})")

# ── 1b. XRONOS (primary radiocarbon source for this corridor) ─────────────
print("\n  --- XRONOS ---")
xronos = pd.read_csv(XRONOS_PATH, low_memory=False)
xronos = xronos.dropna(subset=["lat", "lon"])
print(f"  Loaded {len(xronos)} XRONOS sites with coordinates")

xronos["dist_km"] = dist_from_gc_vec(xronos["lat"].values, xronos["lon"].values)
xr_on = xronos[xronos["dist_km"] < THRESHOLD_KM]
xr_off = xronos[xronos["dist_km"] > 200]
print(f"  On-corridor (<50 km): {len(xr_on)} sites")
print(f"  Off-corridor (>200 km): {len(xr_off)} sites")

phase1_xronos = compare_groups(
    xr_on["n_dates"].dropna().values,
    xr_off["n_dates"].dropna().values,
    label="xronos_overall")

if not phase1_xronos.get("skip"):
    print(f"  Mann-Whitney p = {phase1_xronos['mann_whitney_p']:.4g}")
    print(f"  Median ratio (on/off) = {phase1_xronos['median_ratio']:.3f}")
    print(f"  Bootstrap 95% CI = {phase1_xronos['bootstrap_ci_95']}")
    save_histogram(xr_on["n_dates"].dropna().values,
                   xr_off["n_dates"].dropna().values,
                   "XRONOS: Dates per Site (On- vs Off-Corridor)",
                   "phase1_xronos_histogram.png")
else:
    print(f"  Skipped: {phase1_xronos.get('reason')}")

# XRONOS regional breakdown by country (group into broad regions)
print("\n  XRONOS regional breakdown by country (top represented):")
xr_on_mask = xronos["dist_km"] < THRESHOLD_KM
xr_off_mask = xronos["dist_km"] > 200
phase1_xronos_regional = {}
for country in xronos["country"].dropna().unique():
    on_c = xronos[xr_on_mask & (xronos["country"] == country)]
    off_c = xronos[xr_off_mask & (xronos["country"] == country)]
    if len(on_c) >= 5 and len(off_c) >= 5:
        res = compare_groups(on_c["n_dates"].dropna().values,
                             off_c["n_dates"].dropna().values,
                             label=f"xronos_{country}")
        phase1_xronos_regional[country] = res
        print(f"    {country}: ratio={res['median_ratio']:.3f}, "
              f"p={res['mann_whitney_p']:.4g} (on={res['n_on']}, off={res['n_off']})")


# ===========================================================================
#  PHASE 2: Pleiades Attention Proxy
# ===========================================================================
print("\n" + "=" * 70)
print("PHASE 2: Pleiades Attention Proxy")
print("=" * 70)

pleiades = pd.read_csv(PLEIADES_PATH, low_memory=False)
print(f"  Loaded {len(pleiades)} Pleiades places")

pleiades = pleiades.dropna(subset=["reprLat", "reprLong"])
print(f"  With coordinates: {len(pleiades)}")


def count_connections(val):
    if pd.isna(val) or str(val).strip() == "":
        return 0
    return len(str(val).split(","))


pleiades["n_connections"] = pleiades["hasConnectionsWith"].apply(count_connections)

pleiades["dist_km"] = dist_from_gc_vec(
    pleiades["reprLat"].values, pleiades["reprLong"].values)

pl_on = pleiades[pleiades["dist_km"] < THRESHOLD_KM]
pl_off = pleiades[pleiades["dist_km"] > 200]
print(f"  On-corridor: {len(pl_on)}, Off-corridor: {len(pl_off)}")

phase2_connections = compare_groups(
    pl_on["n_connections"].values, pl_off["n_connections"].values,
    label="pleiades_connections")

if not phase2_connections.get("skip"):
    print(f"  Connections — median ratio={phase2_connections['median_ratio']:.3f}, "
          f"p={phase2_connections['mann_whitney_p']:.4g}")
    print(f"  Bootstrap 95% CI = {phase2_connections['bootstrap_ci_95']}")
else:
    print(f"  Skipped: {phase2_connections.get('reason')}")

# Check for references / bibliography columns
ref_cols = [c for c in pleiades.columns if "ref" in c.lower() or "biblio" in c.lower()]
if ref_cols:
    print(f"  Found reference-like columns: {ref_cols}")
else:
    print("  No references/bibliography columns found in Pleiades data.")

save_histogram(pl_on["n_connections"].values, pl_off["n_connections"].values,
               "Pleiades: Connections per Site (On- vs Off-Corridor)",
               "phase2_pleiades_histogram.png")


# ===========================================================================
#  PHASE 3: Monument vs Settlement Decomposition
# ===========================================================================
print("\n" + "=" * 70)
print("PHASE 3: Monument vs Settlement Decomposition")
print("=" * 70)

MONUMENTAL_TYPES = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus", "shrine"
}
SETTLEMENT_TYPES = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement", "townhouse", "production"
}


def classify_feature_types(ft_str):
    """Parse comma-separated featureTypes and classify."""
    if pd.isna(ft_str) or str(ft_str).strip() == "":
        return "other"
    types = {t.strip().lower() for t in str(ft_str).split(",")}
    is_monument = bool(types & MONUMENTAL_TYPES)
    is_settlement = bool(types & SETTLEMENT_TYPES)
    if is_monument and not is_settlement:
        return "monument"
    if is_settlement and not is_monument:
        return "settlement"
    if is_monument and is_settlement:
        return "mixed"
    return "other"


pleiades["site_class"] = pleiades["featureTypes"].apply(classify_feature_types)

class_counts = pleiades["site_class"].value_counts()
print(f"  Site classifications: {dict(class_counts)}")

phase3_pleiades = {}
for cls in ["monument", "settlement"]:
    subset = pleiades[pleiades["site_class"] == cls]
    s_on = subset[subset["dist_km"] < THRESHOLD_KM]
    s_off = subset[subset["dist_km"] > 200]
    res = compare_groups(s_on["n_connections"].values,
                         s_off["n_connections"].values,
                         label=f"pleiades_{cls}")
    phase3_pleiades[cls] = res
    if res.get("skip"):
        print(f"  {cls}: skipped ({res['reason']}, on={res['n_on']}, off={res['n_off']})")
    else:
        print(f"  {cls}: median ratio={res['median_ratio']:.3f}, "
              f"p={res['mann_whitney_p']:.4g} (on={res['n_on']}, off={res['n_off']})")

# XRONOS monument/settlement decomposition
print("\n  XRONOS site_type decomposition:")
xronos_monument_kw = {"tomb", "temple", "monument", "megalith", "megalithic",
                      "tumulus", "cairn", "dolmen", "menhir", "stone circle",
                      "barrow", "pyramid", "sanctuary", "shrine", "church"}
xronos_settlement_kw = {"settlement", "village", "town", "farmstead", "cave",
                        "shelter", "camp", "habitation", "dwelling", "house",
                        "pit", "occupation"}


def xr_classify(row):
    st = str(row.get("site_type", "")).lower()
    ft = str(row.get("feature_type", "")).lower()
    combined = st + " " + ft
    if any(k in combined for k in xronos_monument_kw):
        return "monument"
    if any(k in combined for k in xronos_settlement_kw):
        return "settlement"
    return "other"


xronos["site_class"] = xronos.apply(xr_classify, axis=1)
xr_class_counts = xronos["site_class"].value_counts()
print(f"  XRONOS classifications: {dict(xr_class_counts)}")

phase3_xronos = {}
for cls in ["monument", "settlement"]:
    subset = xronos[xronos["site_class"] == cls]
    s_on = subset[subset["dist_km"] < THRESHOLD_KM]
    s_off = subset[subset["dist_km"] > 200]
    res = compare_groups(s_on["n_dates"].dropna().values,
                         s_off["n_dates"].dropna().values,
                         label=f"xronos_{cls}")
    phase3_xronos[cls] = res
    if res.get("skip"):
        print(f"  {cls}: skipped ({res['reason']}, on={res['n_on']}, off={res['n_off']})")
    else:
        print(f"  {cls}: median ratio={res['median_ratio']:.3f}, "
              f"p={res['mann_whitney_p']:.4g} (on={res['n_on']}, off={res['n_off']})")


# ===========================================================================
#  PHASE 4: Distance-Decile Analysis
# ===========================================================================
print("\n" + "=" * 70)
print("PHASE 4: Distance-Decile Analysis")
print("=" * 70)

# Use XRONOS for decile analysis (better corridor coverage)
print("\n  --- XRONOS decile analysis ---")
xronos_valid = xronos.dropna(subset=["n_dates"]).copy()
xronos_valid["decile"] = pd.qcut(xronos_valid["dist_km"], 10, labels=False) + 1
xr_decile = xronos_valid.groupby("decile").agg(
    mean_dates=("n_dates", "mean"),
    median_dates=("n_dates", "median"),
    n_sites=("n_dates", "count"),
    mean_dist_km=("dist_km", "mean"),
).reset_index()

print("  Decile | Mean dist (km) | Mean dates | Median dates | N sites")
print("  " + "-" * 65)
for _, row in xr_decile.iterrows():
    print(f"    {int(row['decile']):>2}   | {row['mean_dist_km']:>14.0f} | "
          f"{row['mean_dates']:>10.2f} | {row['median_dates']:>12.1f} | "
          f"{int(row['n_sites']):>7}")

spearman_r_xr, spearman_p_xr = stats.spearmanr(
    xr_decile["decile"], xr_decile["mean_dates"])
print(f"\n  Spearman rho = {spearman_r_xr:.4f}, p = {spearman_p_xr:.4g}")

save_decile_plot(xr_decile, "XRONOS: Dates per Site by Distance Decile",
                 "phase4_xronos_decile.png")

# Also do Pleiades decile analysis
print("\n  --- Pleiades decile analysis ---")
pleiades_copy = pleiades.copy()
pleiades_copy["decile"] = pd.qcut(pleiades_copy["dist_km"], 10, labels=False) + 1
pl_decile = pleiades_copy.groupby("decile").agg(
    mean_dates=("n_connections", "mean"),
    median_dates=("n_connections", "median"),
    n_sites=("n_connections", "count"),
    mean_dist_km=("dist_km", "mean"),
).reset_index()

print("  Decile | Mean dist (km) | Mean conn. | Median conn. | N sites")
print("  " + "-" * 65)
for _, row in pl_decile.iterrows():
    print(f"    {int(row['decile']):>2}   | {row['mean_dist_km']:>14.0f} | "
          f"{row['mean_dates']:>10.2f} | {row['median_dates']:>12.1f} | "
          f"{int(row['n_sites']):>7}")

spearman_r_pl, spearman_p_pl = stats.spearmanr(
    pl_decile["decile"], pl_decile["mean_dates"])
print(f"\n  Spearman rho = {spearman_r_pl:.4f}, p = {spearman_p_pl:.4g}")

save_decile_plot(pl_decile, "Pleiades: Connections per Site by Distance Decile",
                 "phase4_pleiades_decile.png")


# ===========================================================================
#  VERDICT
# ===========================================================================
print("\n" + "=" * 70)
print("VERDICT")
print("=" * 70)

# Primary verdict based on XRONOS (has adequate corridor coverage)
# and Pleiades (independent attention proxy)
xr_p = phase1_xronos.get("mann_whitney_p", 1.0)
xr_ratio = phase1_xronos.get("median_ratio", 1.0)
xr_skip = phase1_xronos.get("skip", False)

pl_p = phase2_connections.get("mann_whitney_p", 1.0)
pl_ratio = phase2_connections.get("median_ratio", 1.0)
pl_skip = phase2_connections.get("skip", False)

# Monument ratios for Phase 3 check
mon_ratio_pl = phase3_pleiades.get("monument", {}).get("median_ratio", 1.0)
set_ratio_pl = phase3_pleiades.get("settlement", {}).get("median_ratio", 1.0)
mon_skip_pl = phase3_pleiades.get("monument", {}).get("skip", True)
set_skip_pl = phase3_pleiades.get("settlement", {}).get("skip", True)

mon_ratio_xr = phase3_xronos.get("monument", {}).get("median_ratio", 1.0)
mon_skip_xr = phase3_xronos.get("monument", {}).get("skip", True)

# Use best available overall ratio
if not xr_skip:
    primary_ratio = xr_ratio
    primary_p = xr_p
    primary_source = "XRONOS"
elif not pl_skip:
    primary_ratio = pl_ratio
    primary_p = pl_p
    primary_source = "Pleiades"
else:
    primary_ratio = 1.0
    primary_p = 1.0
    primary_source = "none"

print(f"  Primary source: {primary_source}")
print(f"  Overall ratio: {primary_ratio:.3f}, p: {primary_p:.4g}")

if primary_p > 0.05 and primary_ratio < 1.5:
    verdict = "CLEAR"
    explanation = (
        f"{primary_source} Mann-Whitney p={primary_p:.4g} (>0.05) and "
        f"median ratio={primary_ratio:.3f} (<1.5). "
        "No significant publication density bias detected."
    )
elif primary_ratio > 3.0:
    # Check if monument-specific
    if (not mon_skip_xr and mon_ratio_xr > 3.0) or \
       (not mon_skip_pl and mon_ratio_pl > 3.0):
        verdict = "FATAL"
        explanation = (
            f"Overall median ratio={primary_ratio:.3f} and monument-specific "
            f"ratio exceeds 3x. Fatal publication bias detected."
        )
    else:
        verdict = "CONCERN"
        explanation = (
            f"Overall median ratio={primary_ratio:.3f} exceeds 3x, but "
            "monument-specific ratio is lower or unavailable. Flagged as concern."
        )
elif 1.5 <= primary_ratio <= 3.0:
    if not mon_skip_pl and not set_skip_pl and \
       abs(mon_ratio_pl - set_ratio_pl) < 0.5:
        verdict = "CONCERN"
        explanation = (
            f"Median ratio={primary_ratio:.3f} is elevated (1.5-3x), but monument "
            f"({mon_ratio_pl:.3f}) and settlement ({set_ratio_pl:.3f}) ratios are "
            "similar, suggesting general attention bias rather than monument-specific."
        )
    elif not mon_skip_pl and mon_ratio_pl > 3.0:
        verdict = "FATAL"
        explanation = (
            f"Monument-specific median ratio={mon_ratio_pl:.3f} exceeds 3x. "
            "Significant monument-targeted publication bias detected."
        )
    else:
        verdict = "CONCERN"
        explanation = (
            f"Median ratio={primary_ratio:.3f} is elevated (1.5-3x). "
            "Phase 3 decomposition partially inconclusive."
        )
else:
    verdict = "CLEAR"
    explanation = "No significant bias detected."

# Additional note about p3k14c coverage gap
p3k_note = ""
if phase1_p3k.get("skip"):
    p3k_note = (
        " Note: p3k14c has zero on-corridor sites (<50 km) due to geographic "
        "sampling bias toward North America and Europe; the great circle passes "
        "through North Africa, South America, and the South Pacific where p3k14c "
        "has minimal coverage. XRONOS serves as the primary radiocarbon proxy."
    )
    explanation += p3k_note

print(f"\n  Result: {verdict}")
print(f"  {explanation}")


# ===========================================================================
#  SAVE OUTPUTS
# ===========================================================================
print("\n" + "=" * 70)
print("SAVING OUTPUTS")
print("=" * 70)

results = {
    "study": "S1: Publication & Excavation Density Bias Test",
    "verdict": verdict,
    "explanation": explanation,
    "phase1": {
        "p3k14c_overall": phase1_p3k,
        "p3k14c_regional": phase1_p3k_regional if not phase1_p3k.get("skip") else {},
        "xronos_overall": phase1_xronos,
        "xronos_regional": phase1_xronos_regional,
    },
    "phase2_pleiades_connections": phase2_connections,
    "phase3_monument_vs_settlement": {
        "pleiades": phase3_pleiades,
        "xronos": phase3_xronos,
    },
    "phase4_decile": {
        "xronos_spearman": {"rho": round(spearman_r_xr, 4),
                            "p": round(spearman_p_xr, 6)},
        "xronos_table": xr_decile.to_dict(orient="records"),
        "pleiades_spearman": {"rho": round(spearman_r_pl, 4),
                              "p": round(spearman_p_pl, 6)},
        "pleiades_table": pl_decile.to_dict(orient="records"),
    },
    "parameters": {
        "pole_lat": POLE_LAT, "pole_lon": POLE_LON,
        "threshold_km": THRESHOLD_KM,
        "off_corridor_min_km": 200,
        "mc_trials": MC_TRIALS,
        "random_seed": 42,
    },
}

with open(OUT_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("  -> saved results.json")


# ── RESULTS.md ────────────────────────────────────────────────────────────────

def _md_res(res):
    """Format a compare_groups result for markdown."""
    if res.get("skip"):
        return f"skipped ({res.get('reason', 'N/A')}, on={res['n_on']}, off={res['n_off']})"
    return (f"median ratio={res['median_ratio']:.3f}, "
            f"p={res['mann_whitney_p']:.4g}, "
            f"CI=[{res['bootstrap_ci_95'][0]:.3f}, {res['bootstrap_ci_95'][1]:.3f}] "
            f"(on={res['n_on']}, off={res['n_off']})")


md = []
md.append("# Study 1: Publication & Excavation Density Bias Test\n")
md.append(f"**Verdict: {verdict}**\n")
md.append(f"{explanation}\n")

md.append("## Phase 1: Radiocarbon Intensity\n")
md.append("### p3k14c\n")
md.append(f"- {_md_res(phase1_p3k)}\n")
if phase1_p3k.get("skip"):
    md.append("- p3k14c geographic coverage does not reach the great circle corridor.\n")

md.append("### XRONOS (primary)\n")
md.append(f"- {_md_res(phase1_xronos)}\n")
if phase1_xronos_regional:
    md.append("\n| Country | On (n) | Off (n) | Median Ratio | p-value |")
    md.append("|---------|--------|---------|-------------|---------|")
    for country, res in sorted(phase1_xronos_regional.items()):
        md.append(f"| {country} | {res['n_on']} | {res['n_off']} | "
                  f"{res['median_ratio']:.3f} | {res['mann_whitney_p']:.4g} |")
    md.append("")

md.append("\n## Phase 2: Pleiades Attention Proxy\n")
md.append(f"- {_md_res(phase2_connections)}\n")

md.append("\n## Phase 3: Monument vs Settlement Decomposition\n")
md.append("### Pleiades\n")
for cls, res in phase3_pleiades.items():
    md.append(f"- **{cls}**: {_md_res(res)}")
md.append("\n### XRONOS\n")
for cls, res in phase3_xronos.items():
    md.append(f"- **{cls}**: {_md_res(res)}")

md.append("\n\n## Phase 4: Distance-Decile Analysis\n")
md.append(f"### XRONOS\n")
md.append(f"- Spearman rho = {spearman_r_xr:.4f}, p = {spearman_p_xr:.4g}\n")
md.append("| Decile | Mean dist (km) | Mean dates | Median dates | N sites |")
md.append("|--------|---------------|-----------|-------------|---------|")
for _, row in xr_decile.iterrows():
    md.append(f"| {int(row['decile'])} | {row['mean_dist_km']:.0f} | "
              f"{row['mean_dates']:.2f} | {row['median_dates']:.1f} | "
              f"{int(row['n_sites'])} |")

md.append(f"\n### Pleiades\n")
md.append(f"- Spearman rho = {spearman_r_pl:.4f}, p = {spearman_p_pl:.4g}\n")
md.append("| Decile | Mean dist (km) | Mean conn. | Median conn. | N sites |")
md.append("|--------|---------------|-----------|-------------|---------|")
for _, row in pl_decile.iterrows():
    md.append(f"| {int(row['decile'])} | {row['mean_dist_km']:.0f} | "
              f"{row['mean_dates']:.2f} | {row['median_dates']:.1f} | "
              f"{int(row['n_sites'])} |")

md.append(f"\n## Parameters\n")
md.append(f"- Pole: ({POLE_LAT}, {POLE_LON})")
md.append(f"- On-corridor threshold: {THRESHOLD_KM} km")
md.append(f"- Off-corridor minimum: 200 km")
md.append(f"- Monte Carlo trials: {MC_TRIALS}")
md.append(f"- Random seed: 42")
md.append("")

with open(OUT_DIR / "RESULTS.md", "w") as f:
    f.write("\n".join(md))
print("  -> saved RESULTS.md")

print("\nDone.")
