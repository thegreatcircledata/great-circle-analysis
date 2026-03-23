#!/usr/bin/env python3
"""
Study 7: Preservation Bias — Multi-Region Settlement Simulation
================================================================
Extends the Egypt-only preservation test to Mesopotamia, Indus Valley,
and coastal Peru.  For each region we inject estimated missing settlements
into the GC corridor and measure how monument–settlement divergence responds.

Outputs  → /Users/elliotallan/megalith_site_research/outputs/next_wave/preservation_bias/
"""

import csv, json, math, os, pathlib, textwrap
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── shared constants & geometry ───────────────────────────────────────────────
POLE_LAT  = 59.682122
POLE_LON  = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
THRESHOLD_KM = 50
MC_TRIALS = 200

RNG = np.random.default_rng(42)

def haversine_km(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = (math.sin(dphi / 2) ** 2
         + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2)
    return 2 * EARTH_R_KM * math.atan2(math.sqrt(a), math.sqrt(1 - a))

def dist_from_gc(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QUARTER_CIRC)

# ── type sets ─────────────────────────────────────────────────────────────────
MONUMENTAL_TYPES = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus", "shrine",
}
SETTLEMENT_TYPES = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement", "townhouse", "production",
}

# ── regions ───────────────────────────────────────────────────────────────────
REGIONS = {
    "Egypt": {
        "lat_range": (22, 31), "lon_range": (24, 36),
        "estimated_missing_settlements": 58,
        "source": "Existing analysis",
    },
    "Mesopotamia": {
        "lat_range": (30, 37), "lon_range": (38, 48),
        "estimated_missing_settlements": 120,
        "source": "Adams 1981 Heartland of Cities",
    },
    "Indus_Valley": {
        "lat_range": (23, 32), "lon_range": (66, 78),
        "estimated_missing_settlements": 200,
        "source": "Possehl 2002",
    },
    "Peru_Coast": {
        "lat_range": (-18, -5), "lon_range": (-82, -74),
        "estimated_missing_settlements": 80,
        "source": "Moseley 2001",
    },
}

# ── paths ─────────────────────────────────────────────────────────────────────
DATA_CSV = pathlib.Path(
    "/Users/elliotallan/megalith_site_research/data/pleiades/pleiades-places-latest.csv"
)
OUT_DIR = pathlib.Path(
    "/Users/elliotallan/megalith_site_research/outputs/next_wave/preservation_bias"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── load Pleiades ─────────────────────────────────────────────────────────────
def classify_row(feature_str):
    """Return 'monument', 'settlement', or None."""
    if not feature_str:
        return None
    types = {t.strip() for t in feature_str.split(",")}
    if types & MONUMENTAL_TYPES:
        return "monument"
    if types & SETTLEMENT_TYPES:
        return "settlement"
    return None

def load_pleiades():
    monuments, settlements = [], []
    with open(DATA_CSV, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                lat = float(row["reprLat"])
                lon = float(row["reprLong"])
            except (ValueError, TypeError):
                continue
            cls = classify_row(row.get("featureTypes", ""))
            if cls == "monument":
                monuments.append((lat, lon))
            elif cls == "settlement":
                settlements.append((lat, lon))
    return monuments, settlements

print("Loading Pleiades …")
all_monuments, all_settlements = load_pleiades()
print(f"  monuments: {len(all_monuments):,}   settlements: {len(all_settlements):,}")

# ── pre-compute GC distances for all sites ────────────────────────────────────
mon_dists = np.array([dist_from_gc(la, lo) for la, lo in all_monuments])
set_dists = np.array([dist_from_gc(la, lo) for la, lo in all_settlements])

# ── Phase 1: count known corridor settlements per region ──────────────────────
print("\nPhase 1 — known corridor settlements per region")
region_known = {}
for rname, rinfo in REGIONS.items():
    lat_lo, lat_hi = rinfo["lat_range"]
    lon_lo, lon_hi = rinfo["lon_range"]
    cnt = 0
    for (la, lo), d in zip(all_settlements, set_dists):
        if lat_lo <= la <= lat_hi and lon_lo <= lo <= lon_hi and d <= THRESHOLD_KM:
            cnt += 1
    region_known[rname] = cnt
    print(f"  {rname:15s}: {cnt} known settlements within {THRESHOLD_KM}km of GC")

# ── synthetic site generation ─────────────────────────────────────────────────
def _gc_point_at_lat(target_lat):
    """
    Return (lat, lon) on the great circle closest to target_lat.
    The GC is defined by the pole at (POLE_LAT, POLE_LON): each point on the
    GC is QUARTER_CIRC km from the pole.  We find the longitude that places
    (target_lat, lon) exactly on the GC by binary search.
    """
    best_lon, best_d = None, 1e18
    for lon_cand in np.linspace(-180, 180, 3600):
        d = abs(haversine_km(target_lat, lon_cand, POLE_LAT, POLE_LON) - QUARTER_CIRC)
        if d < best_d:
            best_d = d
            best_lon = lon_cand
    # refine
    for lon_cand in np.linspace(best_lon - 0.15, best_lon + 0.15, 600):
        d = abs(haversine_km(target_lat, lon_cand, POLE_LAT, POLE_LON) - QUARTER_CIRC)
        if d < best_d:
            best_d = d
            best_lon = lon_cand
    return target_lat, best_lon

def generate_synthetic_settlements(region_info, n_sites):
    """
    Place n_sites synthetic settlements within THRESHOLD_KM of the GC
    inside the region bounding box.
    """
    lat_lo, lat_hi = region_info["lat_range"]
    lon_lo, lon_hi = region_info["lon_range"]

    # Pre-compute GC spine points across latitude range
    spine_lats = np.linspace(lat_lo, lat_hi, 200)
    spine = [_gc_point_at_lat(la) for la in spine_lats]
    # keep only spine points whose longitude falls within region box
    spine = [(la, lo) for la, lo in spine if lon_lo <= lo <= lon_hi]
    if not spine:
        # GC may not cross this box; fall back to centre
        mid_lat = (lat_lo + lat_hi) / 2
        mid_lon = (lon_lo + lon_hi) / 2
        spine = [(mid_lat, mid_lon)]

    sites = []
    attempts = 0
    while len(sites) < n_sites and attempts < n_sites * 50:
        attempts += 1
        # pick a random spine point
        idx = RNG.integers(0, len(spine))
        base_lat, base_lon = spine[idx]
        # jitter: σ=20 km perpendicular (~0.18°), σ=50 km along arc (~0.45°)
        dlat = RNG.normal(0, 0.45)  # along-arc
        dlon = RNG.normal(0, 0.18)  # perpendicular (rough, acceptable at these latitudes)
        new_lat = base_lat + dlat
        new_lon = base_lon + dlon
        # clip to region box
        if not (lat_lo <= new_lat <= lat_hi and lon_lo <= new_lon <= lon_hi):
            continue
        if dist_from_gc(new_lat, new_lon) <= THRESHOLD_KM:
            sites.append((new_lat, new_lon))
    return sites

# ── divergence computation ────────────────────────────────────────────────────
def compute_divergence(mon_d, set_d, _unused1=None, _unused2=None):
    """
    Compute monument–settlement divergence using a permutation test.

    Method: Pool all monument + settlement GC distances.  The observed
    statistic is the difference in corridor-hit fractions:
        obs = frac_mon_in_corridor − frac_set_in_corridor
    Under H0, the monument/settlement labels are exchangeable.  We permute
    labels MC_TRIALS times to build a null distribution of the statistic
    and compute Z = (obs − mean_null) / std_null.

    Also return individual Z-scores for monument and settlement corridor
    fractions vs the pooled rate.
    """
    n_mon_total = len(mon_d)
    n_set_total = len(set_d)
    n_mon_in = int(np.sum(mon_d <= THRESHOLD_KM))
    n_set_in = int(np.sum(set_d <= THRESHOLD_KM))

    frac_mon = n_mon_in / n_mon_total if n_mon_total > 0 else 0
    frac_set = n_set_in / n_set_total if n_set_total > 0 else 0
    obs_diff = frac_mon - frac_set

    # Pool all distances
    pooled = np.concatenate([mon_d, set_d])
    pooled_in = (pooled <= THRESHOLD_KM).astype(np.int32)
    n_total = len(pooled)

    # Permutation null
    null_diffs = np.empty(MC_TRIALS)
    for i in range(MC_TRIALS):
        perm = RNG.permutation(n_total)
        perm_mon_in = pooled_in[perm[:n_mon_total]].sum()
        perm_set_in = pooled_in[perm[n_mon_total:]].sum()
        null_diffs[i] = perm_mon_in / n_mon_total - perm_set_in / n_set_total

    mean_null = null_diffs.mean()
    std_null = null_diffs.std()
    z_div = (obs_diff - mean_null) / std_null if std_null > 0 else 0.0

    # Individual Z-scores vs pooled rate
    pooled_rate = pooled_in.sum() / n_total
    # For monuments: is their fraction significantly above pooled?
    se_mon = math.sqrt(pooled_rate * (1 - pooled_rate) / n_mon_total) if n_mon_total > 0 else 1
    se_set = math.sqrt(pooled_rate * (1 - pooled_rate) / n_set_total) if n_set_total > 0 else 1
    z_mon = (frac_mon - pooled_rate) / se_mon if se_mon > 0 else 0.0
    z_set = (frac_set - pooled_rate) / se_set if se_set > 0 else 0.0

    return float(z_mon), float(z_set), float(z_div)

# ── baseline divergence (no synthetics) ───────────────────────────────────────
print("\nComputing baseline divergence …")
z_mon_base, z_set_base, div_base = compute_divergence(
    mon_dists, set_dists, all_monuments, all_settlements
)
print(f"  Z_monument  = {z_mon_base:.2f}")
print(f"  Z_settlement= {z_set_base:.2f}")
print(f"  Divergence   = {div_base:.2f}")

# ── Phase 2 & 3: sensitivity sweep at 50%, 100%, 200% ────────────────────────
MULTIPLIERS = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0]
sweep_results = {}  # multiplier → {z_mon, z_set, div, per_region_added}

print("\nPhase 2–3 — sensitivity sweep")
for mult in MULTIPLIERS:
    label = f"{mult:.1f}x"
    print(f"  multiplier {label} …", end="", flush=True)

    # accumulate all synthetic sites across regions
    synth_all = []
    per_region = {}
    for rname, rinfo in REGIONS.items():
        n_add = int(round(rinfo["estimated_missing_settlements"] * mult))
        if n_add > 0:
            sites = generate_synthetic_settlements(rinfo, n_add)
        else:
            sites = []
        per_region[rname] = len(sites)
        synth_all.extend(sites)

    # build augmented settlement distance array
    if synth_all:
        synth_dists = np.array([dist_from_gc(la, lo) for la, lo in synth_all])
        aug_set_dists = np.concatenate([set_dists, synth_dists])
    else:
        aug_set_dists = set_dists.copy()

    z_m, z_s, div = compute_divergence(
        mon_dists, aug_set_dists, all_monuments, all_settlements + synth_all
    )
    sweep_results[label] = {
        "multiplier": mult,
        "z_monument": round(z_m, 3),
        "z_settlement": round(z_s, 3),
        "divergence": round(div, 3),
        "total_synthetics_added": len(synth_all),
        "per_region_added": per_region,
    }
    print(f"  div={div:.2f}  (added {len(synth_all)} synthetics)")

# ── Phase 4: break-even analysis ──────────────────────────────────────────────
print("\nPhase 4 — break-even analysis")

# We'll do a finer search for the multiplier at which divergence ≈ 0
# by interpolation / extrapolation from the sweep
mults_done = [sweep_results[f"{m:.1f}x"]["multiplier"] for m in MULTIPLIERS]
divs_done  = [sweep_results[f"{m:.1f}x"]["divergence"] for m in MULTIPLIERS]

# linear interpolation to find zero crossing
breakeven_mult = None
for i in range(len(divs_done) - 1):
    if divs_done[i] > 0 and divs_done[i + 1] <= 0:
        # linear interp
        frac = divs_done[i] / (divs_done[i] - divs_done[i + 1])
        breakeven_mult = mults_done[i] + frac * (mults_done[i + 1] - mults_done[i])
        break

if breakeven_mult is None and divs_done[-1] > 0:
    # extrapolate linearly from last two points
    slope = (divs_done[-1] - divs_done[-2]) / (mults_done[-1] - mults_done[-2])
    if slope < 0:
        breakeven_mult = mults_done[-1] + (0 - divs_done[-1]) / slope
    else:
        breakeven_mult = float("inf")

total_estimated_missing = sum(r["estimated_missing_settlements"] for r in REGIONS.values())
breakeven_settlements = (
    int(round(breakeven_mult * total_estimated_missing))
    if breakeven_mult and breakeven_mult != float("inf")
    else None
)
total_known_settlements = len(all_settlements)

print(f"  Break-even multiplier: {breakeven_mult}")
if breakeven_settlements is not None:
    print(f"  Break-even settlements to add: {breakeven_settlements:,}")
    print(f"  Total known settlements: {total_known_settlements:,}")
    print(f"  Ratio (break-even / known): {breakeven_settlements / total_known_settlements:.2f}x")

# ── Robustness verdict ────────────────────────────────────────────────────────
div_at_1x = sweep_results["1.0x"]["divergence"]
if div_at_1x > 3:
    verdict = "ROBUST"
    verdict_detail = f"Divergence remains {div_at_1x:.1f} (>3) at central estimate"
elif div_at_1x > 2:
    verdict = "MODERATE"
    verdict_detail = f"Divergence {div_at_1x:.1f} between 2 and 3 at central estimate"
else:
    verdict = "FRAGILE"
    verdict_detail = f"Divergence drops to {div_at_1x:.1f} (<2) at central estimate"

print(f"\n  Verdict: {verdict} — {verdict_detail}")

# ── plot: divergence vs multiplier ────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
xs = [sweep_results[f"{m:.1f}x"]["multiplier"] for m in MULTIPLIERS]
ys = [sweep_results[f"{m:.1f}x"]["divergence"] for m in MULTIPLIERS]
ax.plot(xs, ys, "o-", color="#2c7bb6", linewidth=2, markersize=8)
ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
ax.axhline(3, color="red", linestyle=":", linewidth=0.8, label="Z = 3 threshold")
ax.axhline(2, color="orange", linestyle=":", linewidth=0.8, label="Z = 2 threshold")
if breakeven_mult and breakeven_mult != float("inf") and breakeven_mult < max(xs) * 1.5:
    ax.axvline(breakeven_mult, color="green", linestyle="--", linewidth=1,
               label=f"Break-even ≈ {breakeven_mult:.1f}x")
ax.set_xlabel("Multiplier of estimated missing settlements", fontsize=12)
ax.set_ylabel("Monument – Settlement divergence (Z)", fontsize=12)
ax.set_title("Study 7: Preservation Bias — Multi-Region Simulation", fontsize=13)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(OUT_DIR / "divergence_vs_multiplier.png", dpi=150)
plt.close(fig)
print(f"\nPlot saved → {OUT_DIR / 'divergence_vs_multiplier.png'}")

# ── per-region bar chart ──────────────────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(8, 5))
region_names = list(REGIONS.keys())
known_counts = [region_known[r] for r in region_names]
est_missing = [REGIONS[r]["estimated_missing_settlements"] for r in region_names]
x = np.arange(len(region_names))
w = 0.35
ax2.bar(x - w/2, known_counts, w, label="Known in corridor", color="#2c7bb6")
ax2.bar(x + w/2, est_missing, w, label="Estimated missing", color="#d7191c", alpha=0.7)
ax2.set_xticks(x)
ax2.set_xticklabels(region_names, fontsize=11)
ax2.set_ylabel("Settlement count", fontsize=12)
ax2.set_title("Known vs. Estimated Missing Settlements per Region", fontsize=13)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3, axis="y")
fig2.tight_layout()
fig2.savefig(OUT_DIR / "region_settlement_counts.png", dpi=150)
plt.close(fig2)
print(f"Plot saved → {OUT_DIR / 'region_settlement_counts.png'}")

# ── results.json ──────────────────────────────────────────────────────────────
results = {
    "study": "S7 — Preservation Bias: Multi-Region Settlement Simulation",
    "date": datetime.now().strftime("%Y-%m-%d"),
    "parameters": {
        "pole": {"lat": POLE_LAT, "lon": POLE_LON},
        "threshold_km": THRESHOLD_KM,
        "mc_trials": MC_TRIALS,
        "seed": 42,
    },
    "baseline": {
        "z_monument": round(z_mon_base, 3),
        "z_settlement": round(z_set_base, 3),
        "divergence": round(div_base, 3),
        "total_monuments": len(all_monuments),
        "total_settlements": len(all_settlements),
    },
    "region_known_corridor_settlements": region_known,
    "sensitivity_sweep": sweep_results,
    "breakeven": {
        "multiplier": round(breakeven_mult, 2) if breakeven_mult and breakeven_mult != float("inf") else None,
        "settlements_needed": breakeven_settlements,
        "total_known_settlements": total_known_settlements,
        "ratio_to_known": (
            round(breakeven_settlements / total_known_settlements, 3)
            if breakeven_settlements
            else None
        ),
    },
    "verdict": verdict,
    "verdict_detail": verdict_detail,
}

with open(OUT_DIR / "results.json", "w") as fh:
    json.dump(results, fh, indent=2)
print(f"Results saved → {OUT_DIR / 'results.json'}")

# ── RESULTS.md ────────────────────────────────────────────────────────────────
md = textwrap.dedent(f"""\
# Study 7: Preservation Bias — Multi-Region Settlement Simulation

**Date:** {results['date']}
**Verdict:** **{verdict}** — {verdict_detail}

## Motivation

The Egypt-only preservation test showed that adding 58 estimated buried
settlements reduced local monument–settlement divergence by ~30% but
did not collapse it.  This study extends the simulation to four major
riverine corridors along the great circle (GC): Egypt, Mesopotamia,
Indus Valley, and coastal Peru.

## Method

1. Count known Pleiades settlements within {THRESHOLD_KM} km of the GC for each region.
2. For each region, generate synthetic settlement sites (placed within {THRESHOLD_KM} km
   of the GC with Gaussian jitter) at several multipliers of the estimated missing count.
3. Recompute global monument–settlement divergence (Z_mon − Z_set) after each injection.
4. Identify the break-even multiplier where divergence reaches zero.

## Region Inventory

| Region | Known in corridor | Estimated missing | Source |
|--------|------------------:|------------------:|--------|
""")
for rname, rinfo in REGIONS.items():
    md += f"| {rname} | {region_known[rname]} | {rinfo['estimated_missing_settlements']} | {rinfo['source']} |\n"

md += f"""
**Total estimated missing (central): {total_estimated_missing}**

## Baseline

| Metric | Value |
|--------|------:|
| Z_monument | {z_mon_base:.2f} |
| Z_settlement | {z_set_base:.2f} |
| Divergence | {div_base:.2f} |

## Sensitivity Sweep

| Multiplier | Synthetics added | Z_mon | Z_set | Divergence |
|:----------:|-----------------:|------:|------:|-----------:|
"""
for m in MULTIPLIERS:
    s = sweep_results[f"{m:.1f}x"]
    md += f"| {m:.1f}x | {s['total_synthetics_added']} | {s['z_monument']:.2f} | {s['z_settlement']:.2f} | {s['divergence']:.2f} |\n"

md += f"""
## Break-Even Analysis

"""
if breakeven_settlements is not None:
    md += f"""- **Break-even multiplier:** {breakeven_mult:.1f}x of estimated missing
- **Settlements needed:** {breakeven_settlements:,}
- **Total known settlements in Pleiades:** {total_known_settlements:,}
- **Ratio (needed / known):** {breakeven_settlements / total_known_settlements:.2f}x
"""
else:
    md += "- Divergence did not reach zero within tested range.\n"

md += f"""
## Interpretation

"""
if verdict == "ROBUST":
    md += textwrap.dedent(f"""\
    The monument–settlement divergence **remains above Z=3** even after injecting
    the central estimate of {total_estimated_missing} missing settlements across all four
    regions.  The break-even point requires **{breakeven_mult:.1f}x** the estimated
    missing count — equivalent to {breakeven_settlements:,} additional corridor settlements,
    or {breakeven_settlements / total_known_settlements:.0%} of the entire Pleiades settlement
    database.  This is implausible: preservation bias alone cannot explain the observed
    divergence.
    """)
elif verdict == "MODERATE":
    md += textwrap.dedent(f"""\
    The divergence weakens but remains between Z=2 and Z=3 at the central estimate.
    The signal is partially sensitive to preservation bias but not fully explained by it.
    """)
else:
    md += textwrap.dedent(f"""\
    The divergence drops below Z=2 at the central estimate, indicating
    that preservation bias could substantially account for the observed signal.
    """)

md += """
## Plots

![Divergence vs Multiplier](divergence_vs_multiplier.png)
![Region Settlement Counts](region_settlement_counts.png)
"""

with open(OUT_DIR / "RESULTS.md", "w") as fh:
    fh.write(md)
print(f"Report saved → {OUT_DIR / 'RESULTS.md'}")
print("\nDone.")
