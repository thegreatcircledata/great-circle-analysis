#!/usr/bin/env python3
"""
Study 5: Directional Granger Causality Along the Arc
=====================================================
Tests whether monument-building activity propagates directionally along the
great-circle corridor. If diffusion occurs, activity at arc position X at
time T should predict activity at X+delta at time T+lag.

Data: Uses Pleiades (primary, good arc coverage with minDate) and p3k14c
(secondary, wider corridor needed due to sparse coverage on this arc).

Phases:
  1. Arc Position Computation
  2. Autocorrelation Analysis (temporal ACF + spatial Moran's I)
  3. Granger Causality Tests (adjacent arc segments, with null model)
  4. Diffusion Rate Estimation
  5. Sensitivity (250yr and 1000yr bin widths)
"""

import sys, os, math, json, warnings
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout = os.fdopen(sys.stdout.fileno(), "w", buffering=1)
warnings.filterwarnings("ignore", category=FutureWarning)

# ── Shared Constants ──────────────────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
THRESHOLD_KM = 100  # corridor half-width for Pleiades
P3K_THRESHOLD_KM = 500  # wider for p3k14c (sparse on this arc)
GIZA_LAT, GIZA_LON = 29.9792, 31.1342
MC_TRIALS = 1000

np.random.seed(42)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = Path("/Users/elliotallan/megalith_site_research")
P3K_PATH = BASE / "data/p3k14c/p3k14c_data.csv"
PLEIADES_PATH = BASE / "data/pleiades/pleiades-places-latest.csv"
OUT_DIR = BASE / "outputs/next_wave/granger_causality"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Geometry Functions ────────────────────────────────────────────────────────
def haversine_km(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return 2 * EARTH_R_KM * math.atan2(math.sqrt(a), math.sqrt(1 - a))

def dist_from_gc(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QUARTER_CIRC)

def forward_azimuth(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dlam = math.radians(lon2 - lon1)
    x = math.sin(dlam) * math.cos(phi2)
    y = math.cos(phi1) * math.sin(phi2) - math.sin(phi1) * math.cos(phi2) * math.cos(dlam)
    return math.degrees(math.atan2(x, y)) % 360

def arc_position_deg(lat, lon):
    """Angular position along the GC, measured as azimuth difference from Giza."""
    az_site = forward_azimuth(POLE_LAT, POLE_LON, lat, lon)
    az_giza = forward_azimuth(POLE_LAT, POLE_LON, GIZA_LAT, GIZA_LON)
    diff = (az_site - az_giza) % 360
    return diff

# ── Granger Causality (manual implementation) ────────────────────────────────
def manual_granger_test(x, y, max_lag=1):
    """
    Test if x Granger-causes y (past x helps predict y beyond past y alone).
    Restricted:   y[t] = a + b*y[t-1] + e
    Unrestricted: y[t] = a + b*y[t-1] + c*x[t-1] + e
    Returns F-statistic and p-value.
    """
    n = len(y)
    if n <= max_lag + 2:
        return np.nan, 1.0

    y_dep = y[max_lag:]
    y_lag = y[max_lag - 1:-1]
    x_lag = x[max_lag - 1:-1]
    n_obs = len(y_dep)

    X_r = np.column_stack([np.ones(n_obs), y_lag])
    X_u = np.column_stack([np.ones(n_obs), y_lag, x_lag])
    try:
        beta_r = np.linalg.lstsq(X_r, y_dep, rcond=None)[0]
        rss_r = np.sum((y_dep - X_r @ beta_r) ** 2)
        beta_u = np.linalg.lstsq(X_u, y_dep, rcond=None)[0]
        rss_u = np.sum((y_dep - X_u @ beta_u) ** 2)
    except np.linalg.LinAlgError:
        return np.nan, 1.0

    df_resid = n_obs - X_u.shape[1]
    if df_resid <= 0 or rss_u <= 0 or rss_r <= rss_u:
        return 0.0, 1.0

    f_stat = ((rss_r - rss_u) / 1) / (rss_u / df_resid)
    p_value = 1.0 - stats.f.cdf(f_stat, 1, df_resid)
    return f_stat, p_value


def compute_morans_i(values):
    """Moran's I for a 1D spatial series (adjacent weights)."""
    n = len(values)
    if n < 3:
        return np.nan
    mean_v = np.mean(values)
    diffs = values - mean_v
    denom = np.sum(diffs ** 2)
    if denom == 0:
        return np.nan
    numer = np.sum(diffs[:-1] * diffs[1:])
    w_sum = n - 1
    return (n / w_sum) * (numer / denom)


def build_matrix_and_test(sites_df, time_col, time_max, time_bin_width,
                          n_arc_segments=36, mc_trials=MC_TRIALS):
    """
    Build arc-time matrix and run all Granger tests.
    time_col: column name for age (in BP)
    time_max: maximum age in BP to consider
    """
    arc_deg = 360.0 / n_arc_segments
    time_bins_start = np.arange(time_max, 0, -time_bin_width)
    n_time = len(time_bins_start)

    sites = sites_df.copy()
    sites["arc_seg"] = (sites["arc_pos"] / arc_deg).astype(int) % n_arc_segments
    sites["time_bin"] = np.nan
    for i, t_start in enumerate(time_bins_start):
        t_end = t_start - time_bin_width
        mask = (sites[time_col] <= t_start) & (sites[time_col] > t_end)
        sites.loc[mask, "time_bin"] = i
    sites = sites.dropna(subset=["time_bin"])
    sites["time_bin"] = sites["time_bin"].astype(int)

    M = np.zeros((n_arc_segments, n_time), dtype=float)
    for _, row in sites.iterrows():
        M[int(row["arc_seg"]), int(row["time_bin"])] += 1

    # Temporal ACF per arc segment (lag-1)
    temporal_acfs = []
    for seg in range(n_arc_segments):
        ts = M[seg, :]
        if np.std(ts) > 0 and len(ts) > 2:
            temporal_acfs.append(np.corrcoef(ts[:-1], ts[1:])[0, 1])
        else:
            temporal_acfs.append(np.nan)

    # Spatial Moran's I per time step
    spatial_morans = [compute_morans_i(M[:, t]) for t in range(n_time)]

    # Granger Causality: adjacent pairs
    pair_results = []
    for i in range(n_arc_segments - 1):
        ts_i, ts_j = M[i, :], M[i + 1, :]
        f_fwd, p_fwd = manual_granger_test(ts_i, ts_j)
        f_rev, p_rev = manual_granger_test(ts_j, ts_i)
        pair_results.append({
            "seg_from": i, "seg_to": i + 1,
            "arc_from_deg": round(i * arc_deg, 1),
            "arc_to_deg": round((i + 1) * arc_deg, 1),
            "F_forward": round(float(f_fwd), 4) if not np.isnan(f_fwd) else None,
            "p_forward": round(float(p_fwd), 6),
            "F_reverse": round(float(f_rev), 4) if not np.isnan(f_rev) else None,
            "p_reverse": round(float(p_rev), 6),
            "forward_sig": bool(p_fwd < 0.05),
            "reverse_sig": bool(p_rev < 0.05),
        })

    n_fwd_only = sum(1 for r in pair_results if r["forward_sig"] and not r["reverse_sig"])
    n_rev_only = sum(1 for r in pair_results if r["reverse_sig"] and not r["forward_sig"])
    n_both = sum(1 for r in pair_results if r["forward_sig"] and r["reverse_sig"])
    n_neither = sum(1 for r in pair_results if not r["forward_sig"] and not r["reverse_sig"])
    n_pairs = len(pair_results)
    n_any_sig = n_fwd_only + n_rev_only + n_both
    frac_unidirectional = (n_fwd_only + n_rev_only) / n_pairs if n_pairs > 0 else 0
    frac_any_sig = n_any_sig / n_pairs if n_pairs > 0 else 0

    if n_fwd_only + n_rev_only > 0:
        frac_dominant = max(n_fwd_only, n_rev_only) / (n_fwd_only + n_rev_only)
        dominant_dir = "increasing_arc" if n_fwd_only >= n_rev_only else "decreasing_arc"
    else:
        frac_dominant = 0
        dominant_dir = "none"

    # Null model: permute time labels
    null_frac_uni = []
    for trial in range(mc_trials):
        M_perm = M.copy()
        for seg in range(n_arc_segments):
            M_perm[seg, :] = np.random.permutation(M_perm[seg, :])
        fwd_ct = rev_ct = 0
        for ii in range(n_arc_segments - 1):
            _, pp_fwd = manual_granger_test(M_perm[ii, :], M_perm[ii + 1, :])
            _, pp_rev = manual_granger_test(M_perm[ii + 1, :], M_perm[ii, :])
            if pp_fwd < 0.05 and pp_rev >= 0.05:
                fwd_ct += 1
            elif pp_rev < 0.05 and pp_fwd >= 0.05:
                rev_ct += 1
        null_frac_uni.append((fwd_ct + rev_ct) / n_pairs if n_pairs > 0 else 0)

    p_unidirectional = float(np.mean(np.array(null_frac_uni) >= frac_unidirectional))

    return {
        "time_bin_width": time_bin_width,
        "n_arc_segments": n_arc_segments,
        "n_time_bins": n_time,
        "n_sites_used": len(sites),
        "matrix_shape": [n_arc_segments, n_time],
        "matrix_nonzero_cells": int(np.sum(M > 0)),
        "total_counts": int(np.sum(M)),
        "occupied_segments": int(np.sum(M.sum(axis=1) > 0)),
        "temporal_acf_mean": float(np.nanmean(temporal_acfs)) if any(not np.isnan(x) for x in temporal_acfs) else None,
        "temporal_acf_median": float(np.nanmedian(temporal_acfs)) if any(not np.isnan(x) for x in temporal_acfs) else None,
        "spatial_moran_mean": float(np.nanmean(spatial_morans)) if any(not np.isnan(x) for x in spatial_morans) else None,
        "spatial_moran_median": float(np.nanmedian(spatial_morans)) if any(not np.isnan(x) for x in spatial_morans) else None,
        "n_forward_only_sig": n_fwd_only,
        "n_reverse_only_sig": n_rev_only,
        "n_both_sig": n_both,
        "n_neither_sig": n_neither,
        "n_any_sig": n_any_sig,
        "frac_unidirectional": round(frac_unidirectional, 4),
        "frac_any_significant": round(frac_any_sig, 4),
        "p_unidirectional_null": round(p_unidirectional, 4),
        "dominant_direction": dominant_dir,
        "frac_dominant": round(frac_dominant, 4),
        "pair_results": pair_results,
        "matrix": M.tolist(),
        "null_frac_unidirectional_dist": [round(x, 4) for x in sorted(null_frac_uni)],
    }


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("Study 5: Directional Granger Causality Along the Arc")
print("=" * 70)

# ── Load Pleiades (primary source) ───────────────────────────────────────────
print("\n[Phase 0a] Loading Pleiades data...")
pl = pd.read_csv(PLEIADES_PATH, low_memory=False)
pl = pl.dropna(subset=["reprLat", "reprLong", "minDate"])
pl["minDate"] = pd.to_numeric(pl["minDate"], errors="coerce")
pl = pl.dropna(subset=["minDate"])
print(f"  Total Pleiades records with coordinates and dates: {len(pl):,}")

# Convert minDate (CE, negative=BCE) to BP (before 1950)
pl["age_bp"] = 1950 - pl["minDate"]  # e.g., -3000 CE -> 4950 BP
pl = pl[pl["age_bp"] > 0]  # only pre-1950

# Distance from GC
print("  Computing distances from great circle...")
pl["dist_gc_km"] = pl.apply(lambda r: dist_from_gc(r["reprLat"], r["reprLong"]), axis=1)
pl_corridor = pl[pl["dist_gc_km"] <= THRESHOLD_KM].copy()
print(f"  Pleiades within {THRESHOLD_KM} km of GC: {len(pl_corridor):,} records")

# Arc position
pl_corridor["arc_pos"] = pl_corridor.apply(
    lambda r: arc_position_deg(r["reprLat"], r["reprLong"]), axis=1
)

# Use id as unique site identifier (already unique in Pleiades)
pl_sites = pl_corridor.rename(columns={"reprLat": "lat", "reprLong": "lon"})[
    ["id", "lat", "lon", "age_bp", "arc_pos"]
].copy()
pl_sites = pl_sites[(pl_sites["age_bp"] >= 0) & (pl_sites["age_bp"] <= 12000)]
print(f"  Pleiades sites on corridor (0-12000 BP): {len(pl_sites):,}")

# ── Load p3k14c (secondary source, wider corridor) ──────────────────────────
print("\n[Phase 0b] Loading p3k14c data...")
df = pd.read_csv(P3K_PATH, low_memory=False)
df = df.dropna(subset=["Lat", "Long", "Age", "SiteID"])
df["Age"] = pd.to_numeric(df["Age"], errors="coerce")
df = df.dropna(subset=["Age"])
df = df[df["Age"] > 0]
print(f"  Total p3k14c records: {len(df):,}")

df["dist_gc_km"] = df.apply(lambda r: dist_from_gc(r["Lat"], r["Long"]), axis=1)

# Report coverage issue
n_100 = (df["dist_gc_km"] <= 100).sum()
n_500 = (df["dist_gc_km"] <= 500).sum()
print(f"  p3k14c within 100 km of GC: {n_100} (insufficient for analysis)")
print(f"  p3k14c within 500 km of GC: {n_500}")
print(f"  NOTE: p3k14c has negligible coverage along this great circle.")
print(f"  Using {P3K_THRESHOLD_KM} km corridor for supplementary p3k14c analysis.")

p3k_corridor = df[df["dist_gc_km"] <= P3K_THRESHOLD_KM].copy()
# Aggregate by SiteID
p3k_sites = p3k_corridor.groupby("SiteID").agg(
    lat=("Lat", "mean"), lon=("Long", "mean"),
    age_bp=("Age", "min"), n_dates=("Age", "count"),
).reset_index()
p3k_sites["arc_pos"] = p3k_sites.apply(
    lambda r: arc_position_deg(r["lat"], r["lon"]), axis=1
)
p3k_sites = p3k_sites[(p3k_sites["age_bp"] >= 0) & (p3k_sites["age_bp"] <= 12000)]
print(f"  p3k14c sites within {P3K_THRESHOLD_KM} km (0-12000 BP): {len(p3k_sites):,}")

# ── Phase 1: Show arc distributions ─────────────────────────────────────────
print("\n[Phase 1] Arc position distributions...")
for label, sdf in [("Pleiades", pl_sites), ("p3k14c", p3k_sites)]:
    arc_hist, _ = np.histogram(sdf["arc_pos"], bins=36, range=(0, 360))
    n_occ = np.sum(arc_hist > 0)
    print(f"  {label}: {len(sdf)} sites, {n_occ}/36 segments occupied, "
          f"min={arc_hist.min()}, max={arc_hist.max()}, mean={arc_hist.mean():.1f}")

# ── Phase 2-3: Pleiades primary analysis (500yr bins) ────────────────────────
print("\n[Phase 2-3] Pleiades Granger causality (500yr bins)...")
results_pl_500 = build_matrix_and_test(pl_sites, "age_bp", 12000, 500)
print(f"  Matrix: {results_pl_500['matrix_shape']} with {results_pl_500['matrix_nonzero_cells']} non-zero cells")
print(f"  Occupied segments: {results_pl_500['occupied_segments']}/36")
print(f"  Temporal ACF (mean): {results_pl_500['temporal_acf_mean']}")
print(f"  Spatial Moran's I (mean): {results_pl_500['spatial_moran_mean']}")
print(f"  Forward-only sig: {results_pl_500['n_forward_only_sig']}/{results_pl_500['n_arc_segments']-1}")
print(f"  Reverse-only sig: {results_pl_500['n_reverse_only_sig']}/{results_pl_500['n_arc_segments']-1}")
print(f"  Both sig: {results_pl_500['n_both_sig']}")
print(f"  Neither sig: {results_pl_500['n_neither_sig']}")
print(f"  Frac unidirectional: {results_pl_500['frac_unidirectional']:.3f}")
print(f"  Null p-value: {results_pl_500['p_unidirectional_null']:.4f}")
print(f"  Dominant direction: {results_pl_500['dominant_direction']}")

# ── p3k14c supplementary analysis ───────────────────────────────────────────
print("\n[Phase 2-3b] p3k14c Granger causality (500yr bins, 500km corridor)...")
results_p3k_500 = build_matrix_and_test(p3k_sites, "age_bp", 12000, 500)
print(f"  Matrix: {results_p3k_500['matrix_shape']} with {results_p3k_500['matrix_nonzero_cells']} non-zero cells")
print(f"  Occupied segments: {results_p3k_500['occupied_segments']}/36")
print(f"  Frac unidirectional: {results_p3k_500['frac_unidirectional']:.3f}")
print(f"  Null p-value: {results_p3k_500['p_unidirectional_null']:.4f}")

# ── Phase 4: Diffusion Rate Estimation ───────────────────────────────────────
print("\n[Phase 4] Diffusion rate estimation (Pleiades)...")
M_primary = np.array(results_pl_500["matrix"])
lag_results = {}

for lag in [1, 2, 3]:
    f_vals, p_vals = [], []
    for i in range(M_primary.shape[0] - 1):
        ts_i, ts_j = M_primary[i, :], M_primary[i + 1, :]
        if lag == 1:
            f, p = manual_granger_test(ts_i, ts_j)
        else:
            n = len(ts_j)
            if n <= lag + 2:
                continue
            y_dep = ts_j[lag:]
            y_lag_arr = ts_j[lag - 1:-1]
            x_lag_arr = ts_i[:n - lag]
            n_obs = min(len(y_dep), len(y_lag_arr), len(x_lag_arr))
            y_dep, y_lag_arr, x_lag_arr = y_dep[:n_obs], y_lag_arr[:n_obs], x_lag_arr[:n_obs]
            X_r = np.column_stack([np.ones(n_obs), y_lag_arr])
            X_u = np.column_stack([np.ones(n_obs), y_lag_arr, x_lag_arr])
            try:
                rss_r = np.sum((y_dep - X_r @ np.linalg.lstsq(X_r, y_dep, rcond=None)[0]) ** 2)
                rss_u = np.sum((y_dep - X_u @ np.linalg.lstsq(X_u, y_dep, rcond=None)[0]) ** 2)
                df_resid = n_obs - 3
                if df_resid > 0 and rss_u > 0 and rss_r > rss_u:
                    f = ((rss_r - rss_u) / 1) / (rss_u / df_resid)
                    p = 1.0 - stats.f.cdf(f, 1, df_resid)
                else:
                    f, p = 0.0, 1.0
            except Exception:
                f, p = np.nan, 1.0
        if not np.isnan(f):
            f_vals.append(f)
            p_vals.append(p)

    avg_f = np.mean(f_vals) if f_vals else 0
    n_sig = sum(1 for p in p_vals if p < 0.05)
    lag_results[lag] = {
        "avg_F": round(float(avg_f), 3),
        "n_sig": n_sig,
        "n_tested": len(p_vals),
        "frac_sig": round(n_sig / len(p_vals), 3) if p_vals else 0,
    }

best_lag = max(lag_results, key=lambda l: lag_results[l]["avg_F"]) if lag_results else None

print(f"  Lag analysis:")
for lag, lr in lag_results.items():
    print(f"    Lag {lag}: avg_F={lr['avg_F']:.3f}, sig={lr['n_sig']}/{lr['n_tested']} ({lr['frac_sig']:.1%})")

arc_segment_deg = 360.0 / 36
arc_segment_km = (arc_segment_deg / 360.0) * 2 * math.pi * EARTH_R_KM
diffusion_rate_km_yr = None
if best_lag and lag_results[best_lag]["avg_F"] > 0:
    time_per_lag = 500 * best_lag
    diffusion_rate_km_yr = arc_segment_km / time_per_lag
    print(f"  Best lag: {best_lag} ({time_per_lag} yr)")
    print(f"  Arc segment width: {arc_segment_km:.0f} km")
    print(f"  Implied diffusion rate: {diffusion_rate_km_yr:.2f} km/yr")
    print(f"  (Neolithic farming spread: ~1 km/yr)")
else:
    print("  No significant directional lag found")

# ── Phase 5: Sensitivity ────────────────────────────────────────────────────
print("\n[Phase 5] Sensitivity analysis (Pleiades)...")
print("  250yr bins...")
results_pl_250 = build_matrix_and_test(pl_sites, "age_bp", 12000, 250)
print(f"    Frac uni: {results_pl_250['frac_unidirectional']:.3f}, "
      f"p={results_pl_250['p_unidirectional_null']:.4f}")

print("  1000yr bins...")
results_pl_1000 = build_matrix_and_test(pl_sites, "age_bp", 12000, 1000)
print(f"    Frac uni: {results_pl_1000['frac_unidirectional']:.3f}, "
      f"p={results_pl_1000['p_unidirectional_null']:.4f}")

# ── Verdict ──────────────────────────────────────────────────────────────────
print("\n[Verdict]")
frac_uni = results_pl_500["frac_unidirectional"]
p_uni = results_pl_500["p_unidirectional_null"]

if frac_uni > 0.6 and p_uni < 0.01:
    verdict = "STRONG"
    explanation = (f">{frac_uni*100:.0f}% of adjacent arc-segment pairs show significant "
                   f"unidirectional Granger causality (p_null={p_uni:.4f})")
elif frac_uni > 0.3 or (frac_uni > 0.15 and p_uni < 0.05):
    verdict = "SUGGESTIVE"
    explanation = (f"{frac_uni*100:.0f}% of pairs show unidirectional Granger causality "
                   f"(p_null={p_uni:.4f})")
else:
    verdict = "NULL"
    explanation = (f"Only {frac_uni*100:.0f}% of pairs show unidirectional causality "
                   f"(p_null={p_uni:.4f}), insufficient for diffusion claim")

print(f"  Verdict: {verdict}")
print(f"  {explanation}")

# ══════════════════════════════════════════════════════════════════════════════
# OUTPUTS
# ══════════════════════════════════════════════════════════════════════════════
print("\n[Outputs] Generating plots...")

# ── Plot 1: Arc-position vs time heatmap ─────────────────────────────────────
M_plot = np.array(results_pl_500["matrix"])
fig, ax = plt.subplots(figsize=(14, 7))
im = ax.imshow(M_plot, aspect="auto", cmap="YlOrRd", interpolation="nearest", origin="lower")
ax.set_xlabel("Time Bin (BP, descending)", fontsize=12)
ax.set_ylabel("Arc Segment (degrees from Giza)", fontsize=12)
ax.set_title(f"Pleiades Site Counts: Arc Position vs Time (500yr bins, {THRESHOLD_KM}km corridor)", fontsize=13)

time_labels = [f"{t}" for t in np.arange(12000, 0, -500)]
xtick_pos = np.arange(0, len(time_labels), 4)
ax.set_xticks(xtick_pos)
ax.set_xticklabels([time_labels[i] for i in xtick_pos], rotation=45, ha="right")
ytick_pos = np.arange(0, 36, 6)
ax.set_yticks(ytick_pos)
ax.set_yticklabels([f"{i*10}\u00b0" for i in ytick_pos])

plt.colorbar(im, ax=ax, label="Site Count")
plt.tight_layout()
plt.savefig(OUT_DIR / "arc_time_heatmap.png", dpi=150)
plt.close()
print("  Saved arc_time_heatmap.png")

# ── Plot 2: Granger causality summary ────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

pairs = results_pl_500["pair_results"]
x_pos = [r["seg_from"] for r in pairs]
p_fwd = [r["p_forward"] for r in pairs]
p_rev = [r["p_reverse"] for r in pairs]

axes[0].scatter(x_pos, p_fwd, c="blue", alpha=0.7, label="Forward (i\u2192i+1)", s=30)
axes[0].scatter(x_pos, p_rev, c="red", alpha=0.7, label="Reverse (i+1\u2192i)", s=30)
axes[0].axhline(0.05, color="black", linestyle="--", linewidth=1, label="p=0.05")
axes[0].set_xlabel("Arc Segment Index")
axes[0].set_ylabel("p-value")
axes[0].set_title("Granger Causality p-values by Segment Pair")
axes[0].legend(fontsize=9)
axes[0].set_ylim(-0.02, 1.02)

null_dist = results_pl_500["null_frac_unidirectional_dist"]
axes[1].hist(null_dist, bins=30, color="gray", alpha=0.7, edgecolor="black", label="Null")
axes[1].axvline(frac_uni, color="red", linewidth=2, label=f"Observed ({frac_uni:.3f})")
axes[1].set_xlabel("Fraction Unidirectional")
axes[1].set_ylabel("Count (null trials)")
axes[1].set_title(f"Null Model Comparison (p={p_uni:.4f})")
axes[1].legend()

plt.tight_layout()
plt.savefig(OUT_DIR / "granger_summary.png", dpi=150)
plt.close()
print("  Saved granger_summary.png")

# ── Plot 3: Diffusion lag ────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
lags = sorted(lag_results.keys())
fracs = [lag_results[l]["frac_sig"] for l in lags]
ax.bar(lags, fracs, color="steelblue", edgecolor="black")
ax.set_xlabel("Lag (time bins of 500yr)")
ax.set_ylabel("Fraction of Pairs with p < 0.05")
ax.set_title("Granger Causality Significance by Lag")
ax.set_xticks(lags)
ax.set_xticklabels([f"Lag {l}\n({l*500}yr)" for l in lags])
plt.tight_layout()
plt.savefig(OUT_DIR / "diffusion_lag_plot.png", dpi=150)
plt.close()
print("  Saved diffusion_lag_plot.png")

# ── Plot 4: Sensitivity ─────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
bin_widths = [250, 500, 1000]
frac_unis = [results_pl_250["frac_unidirectional"], results_pl_500["frac_unidirectional"],
             results_pl_1000["frac_unidirectional"]]
p_vals_sens = [results_pl_250["p_unidirectional_null"], results_pl_500["p_unidirectional_null"],
               results_pl_1000["p_unidirectional_null"]]
colors = ["green" if p < 0.05 else "gray" for p in p_vals_sens]
ax.bar(range(3), frac_unis, color=colors, edgecolor="black")
ax.set_xticks(range(3))
ax.set_xticklabels([f"{w}yr\n(p={p:.3f})" for w, p in zip(bin_widths, p_vals_sens)])
ax.set_ylabel("Fraction Unidirectional Pairs")
ax.set_title("Sensitivity: Time Bin Width (Pleiades)")
ax.axhline(0.3, color="orange", linestyle="--", linewidth=1, label="Suggestive (30%)")
ax.axhline(0.6, color="red", linestyle="--", linewidth=1, label="Strong (60%)")
ax.legend()
plt.tight_layout()
plt.savefig(OUT_DIR / "sensitivity_binwidth.png", dpi=150)
plt.close()
print("  Saved sensitivity_binwidth.png")

# ── results.json ─────────────────────────────────────────────────────────────
def clean_for_json(d):
    """Strip large arrays for JSON output."""
    return {k: v for k, v in d.items() if k not in ("matrix", "null_frac_unidirectional_dist")}

output = {
    "study": "S5: Directional Granger Causality Along the Arc",
    "verdict": verdict,
    "explanation": explanation,
    "parameters": {
        "pleiades_corridor_km": THRESHOLD_KM,
        "p3k14c_corridor_km": P3K_THRESHOLD_KM,
        "n_arc_segments": 36,
        "mc_trials": MC_TRIALS,
        "pole": [POLE_LAT, POLE_LON],
        "giza_ref": [GIZA_LAT, GIZA_LON],
    },
    "data_note": (f"p3k14c has {n_100} records within 100km of the great circle "
                  f"(effectively zero coverage). Pleiades ({len(pl_sites)} sites within "
                  f"{THRESHOLD_KM}km) is used as the primary source. p3k14c is run "
                  f"supplementarily with a {P3K_THRESHOLD_KM}km corridor."),
    "pleiades_500yr": {k: v for k, v in results_pl_500.items()
                       if k not in ("matrix",)},
    "pleiades_250yr": clean_for_json(results_pl_250),
    "pleiades_1000yr": clean_for_json(results_pl_1000),
    "p3k14c_500yr": clean_for_json(results_p3k_500),
    "diffusion": {
        "lag_results": lag_results,
        "best_lag": best_lag,
        "diffusion_rate_km_per_yr": round(diffusion_rate_km_yr, 3) if diffusion_rate_km_yr else None,
        "arc_segment_km": round(arc_segment_km, 1),
    },
}

with open(OUT_DIR / "results.json", "w") as f:
    json.dump(output, f, indent=2, default=str)
print("  Saved results.json")

# ── RESULTS.md ───────────────────────────────────────────────────────────────
n_pairs = results_pl_500["n_arc_segments"] - 1

md = f"""# Study 5: Directional Granger Causality Along the Arc

## Verdict: **{verdict}**

{explanation}

## Data Note

The p3k14c radiocarbon database has only **{n_100} records within 100 km** of this
great circle (the arc passes through the Sahara, West Africa, South America, and
the Pacific where p3k14c has minimal coverage). Pleiades ({len(pl_sites):,} sites
within {THRESHOLD_KM} km) provides the primary dataset for this analysis, with
p3k14c run supplementarily using a wider {P3K_THRESHOLD_KM} km corridor
({len(p3k_sites):,} sites).

## Method

Tests whether archaeological activity propagates directionally along the
great-circle corridor using Granger causality on adjacent arc segments.

- **Corridor**: {THRESHOLD_KM} km (Pleiades) / {P3K_THRESHOLD_KM} km (p3k14c)
- **Arc**: 36 segments of 10 degrees each (reference: Giza at 0 degrees)
- **Time**: 500-year bins from 12000-0 BP (primary analysis)
- **Granger test**: F-test comparing restricted (AR(1)) vs unrestricted (AR(1) + lagged neighbor) model

## Phase 1: Arc-Time Matrix (Pleiades, 500yr)

- Matrix shape: {results_pl_500['matrix_shape'][0]} arc segments x {results_pl_500['matrix_shape'][1]} time bins
- Non-zero cells: {results_pl_500['matrix_nonzero_cells']}
- Occupied segments: {results_pl_500['occupied_segments']}/36
- Total site counts: {results_pl_500['total_counts']:,}

## Phase 2: Autocorrelation

| Metric | Mean | Median |
|--------|------|--------|
| Temporal ACF (lag-1) | {f"{results_pl_500['temporal_acf_mean']:.3f}" if results_pl_500['temporal_acf_mean'] else 'N/A'} | {f"{results_pl_500['temporal_acf_median']:.3f}" if results_pl_500['temporal_acf_median'] else 'N/A'} |
| Spatial Moran's I | {f"{results_pl_500['spatial_moran_mean']:.3f}" if results_pl_500['spatial_moran_mean'] else 'N/A'} | {f"{results_pl_500['spatial_moran_median']:.3f}" if results_pl_500['spatial_moran_median'] else 'N/A'} |

## Phase 3: Granger Causality (Pleiades, 500yr bins)

| Category | Count | Fraction |
|----------|-------|----------|
| Forward only (i->i+1) | {results_pl_500['n_forward_only_sig']} | {results_pl_500['n_forward_only_sig']/n_pairs:.1%} |
| Reverse only (i+1->i) | {results_pl_500['n_reverse_only_sig']} | {results_pl_500['n_reverse_only_sig']/n_pairs:.1%} |
| Both directions | {results_pl_500['n_both_sig']} | {results_pl_500['n_both_sig']/n_pairs:.1%} |
| Neither | {results_pl_500['n_neither_sig']} | {results_pl_500['n_neither_sig']/n_pairs:.1%} |

- **Fraction unidirectional**: {results_pl_500['frac_unidirectional']:.3f}
- **Null model p-value**: {results_pl_500['p_unidirectional_null']:.4f} ({MC_TRIALS:,} permutations)
- **Dominant direction**: {results_pl_500['dominant_direction']}

### p3k14c Supplementary ({P3K_THRESHOLD_KM}km corridor)

| Metric | Value |
|--------|-------|
| Sites used | {results_p3k_500['n_sites_used']} |
| Occupied segments | {results_p3k_500['occupied_segments']}/36 |
| Frac unidirectional | {results_p3k_500['frac_unidirectional']:.3f} |
| Null p-value | {results_p3k_500['p_unidirectional_null']:.4f} |

## Phase 4: Diffusion Rate

| Lag | Duration | Fraction Significant | Avg F |
|-----|----------|---------------------|-------|
"""

for lag in sorted(lag_results.keys()):
    lr = lag_results[lag]
    md += f"| {lag} | {lag*500} yr | {lr['frac_sig']:.1%} ({lr['n_sig']}/{lr['n_tested']}) | {lr['avg_F']:.3f} |\n"

if diffusion_rate_km_yr:
    md += f"""
- **Best lag**: {best_lag} time bins = {best_lag * 500} years
- **Arc segment width**: {arc_segment_km:.0f} km
- **Implied diffusion rate**: {diffusion_rate_km_yr:.2f} km/yr
- **Comparison**: Neolithic farming spread ~1 km/yr (Ammerman & Cavalli-Sforza)
"""
else:
    md += "\nNo significant directional diffusion detected.\n"

md += f"""
## Phase 5: Sensitivity

| Bin Width | Frac Unidirectional | Null p-value | Dominant Dir |
|-----------|--------------------:|-------------:|--------------|
| 250 yr | {results_pl_250['frac_unidirectional']:.3f} | {results_pl_250['p_unidirectional_null']:.4f} | {results_pl_250['dominant_direction']} |
| 500 yr | {results_pl_500['frac_unidirectional']:.3f} | {results_pl_500['p_unidirectional_null']:.4f} | {results_pl_500['dominant_direction']} |
| 1000 yr | {results_pl_1000['frac_unidirectional']:.3f} | {results_pl_1000['p_unidirectional_null']:.4f} | {results_pl_1000['dominant_direction']} |

## Criteria Applied

- **STRONG**: >60% adjacent pairs show significant unidirectional Granger causality, null p < 0.01
- **SUGGESTIVE**: 30-60% significant with some directional bias
- **NULL**: <30% significant or no directional bias

## Output Files

- `arc_time_heatmap.png` -- Arc position vs time heatmap
- `granger_summary.png` -- Per-pair p-values and null distribution
- `diffusion_lag_plot.png` -- Significance by lag
- `sensitivity_binwidth.png` -- Sensitivity to bin width
- `results.json` -- Full numerical results
"""

with open(OUT_DIR / "RESULTS.md", "w") as f:
    f.write(md)
print("  Saved RESULTS.md")

print("\n" + "=" * 70)
print(f"Study 5 complete. Verdict: {verdict}")
print("=" * 70)
