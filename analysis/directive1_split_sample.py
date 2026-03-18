#!/usr/bin/env python3
"""
Directive 1: Split-Sample Blinded Validation Test
===================================================
Addresses circularity concern: Alison defined the Great Circle by observing
that specific famous sites fell on it, then we test that same circle against
databases containing those sites.

This script performs:
  A) 100 random 50/50 splits of the Megalithic Portal database
     - For each split: find optimal pole on Discovery half, test on Validation half
     - Also test Alison pole on Validation half
  B) "15 Famous Sites" control: fit optimal circles to 1,000 random sets of 15 sites
     and compare their Z-scores to Alison's

OPTIMIZED: Uses numpy vectorization for all distance calculations.
"""

import xml.etree.ElementTree as ET
import glob, math, random, json, os, re, sys, time
import numpy as np
from scipy import stats as sp_stats
from scipy.optimize import minimize
from collections import Counter

# Force unbuffered output
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

# ============================================================
# CONFIGURATION
# ============================================================
POLE_LAT = 59.682122
POLE_LNG = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * np.pi / 2
THRESHOLD_KM = 50
N_MC_TRIALS = 200
N_SPLITS = 100
N_FAMOUS_TRIALS = 200
N_FAMOUS_SITES = 15
POLE_SCAN_RES = 2  # degrees

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "split_sample_validation")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# VECTORIZED HELPER FUNCTIONS
# ============================================================
def haversine_vec(lat1, lon1, lat2_arr, lon2_arr):
    """Vectorized haversine: one point vs arrays of points. Returns km."""
    lat1r = np.radians(lat1)
    lon1r = np.radians(lon1)
    lat2r = np.radians(lat2_arr)
    lon2r = np.radians(lon2_arr)
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat/2)**2 + np.cos(lat1r)*np.cos(lat2r)*np.sin(dlon/2)**2
    return EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))

def dist_from_gc_vec(site_lats, site_lons, pole_lat, pole_lon):
    """Vectorized distance from great circle for arrays of sites."""
    d = haversine_vec(pole_lat, pole_lon, site_lats, site_lons)
    return np.abs(d - QUARTER_CIRC)

def count_within_vec(site_lats, site_lons, pole_lat, pole_lon, threshold):
    """Count sites within threshold km of the great circle defined by pole."""
    dists = dist_from_gc_vec(site_lats, site_lons, pole_lat, pole_lon)
    return int(np.sum(dists <= threshold))

def mc_zscore_vec(site_lats, site_lons, pole_lat, pole_lon, threshold, n_trials):
    """Vectorized Monte Carlo Z-score."""
    n = len(site_lats)
    observed = count_within_vec(site_lats, site_lons, pole_lat, pole_lon, threshold)

    rand_counts = np.empty(n_trials)
    for t in range(n_trials):
        idx = np.random.randint(0, n, size=n)
        rlats = site_lats[idx] + np.random.normal(0, 2, n)
        rlons = site_lons[idx] + np.random.normal(0, 2, n)
        rlats = np.clip(rlats, -90, 90)
        rlons = np.clip(rlons, -180, 180)
        rand_counts[t] = count_within_vec(rlats, rlons, pole_lat, pole_lon, threshold)

    mu = np.mean(rand_counts)
    sigma = np.std(rand_counts)
    z = (observed - mu) / sigma if sigma > 0 else 0
    return float(z), int(observed), float(mu), float(sigma)

# ============================================================
# LOAD DATA
# ============================================================
print("=" * 70)
print("LOADING MEGALITHIC PORTAL DATA")
print("=" * 70)

all_sites = []
for filepath in sorted(glob.glob(os.path.join(BASE_DIR, "MegP_*.kml"))):
    try:
        tree = ET.parse(filepath)
        root = tree.getroot()
        for pm in root.iter():
            if 'Placemark' not in pm.tag: continue
            coords_el = None
            for child in pm.iter():
                if 'coordinates' in child.tag:
                    coords_el = child
            if coords_el is not None and coords_el.text:
                try:
                    parts = coords_el.text.strip().split(',')
                    lon, lat = float(parts[0]), float(parts[1])
                    if -90 <= lat <= 90 and -180 <= lon <= 180 and (lat != 0 or lon != 0):
                        all_sites.append((lat, lon))
                except (ValueError, IndexError):
                    pass
    except ET.ParseError:
        pass

# Deduplicate
seen = set()
sites_list = []
for lat, lon in all_sites:
    key = (round(lat, 3), round(lon, 3))
    if key not in seen:
        seen.add(key)
        sites_list.append((lat, lon))

# Convert to numpy arrays for vectorized operations
sites_arr = np.array(sites_list)
all_lats = sites_arr[:, 0]
all_lons = sites_arr[:, 1]
N_SITES = len(sites_list)
print(f"Loaded {N_SITES} unique sites")

# ============================================================
# BUILD POLE SCAN GRID
# ============================================================
print(f"\nBuilding pole scan grid at {POLE_SCAN_RES}° resolution...")
pole_lats = []
pole_lons = []
for plat in range(-90, 91, POLE_SCAN_RES):
    for plon in range(-180, 180, POLE_SCAN_RES):
        pole_lats.append(plat)
        pole_lons.append(plon)
pole_lats = np.array(pole_lats)
pole_lons = np.array(pole_lons)
N_POLES = len(pole_lats)
print(f"Grid has {N_POLES} candidate poles")

# ============================================================
# PART A: 100 SPLIT-SAMPLE VALIDATION RUNS
# ============================================================
print("\n" + "=" * 70)
print(f"PART A: {N_SPLITS} SPLIT-SAMPLE VALIDATION RUNS")
print("=" * 70)

split_results = []
t_start = time.time()

for split_idx in range(N_SPLITS):
    t_split = time.time()

    # Random 50/50 split
    indices = np.random.permutation(N_SITES)
    mid = N_SITES // 2
    disc_idx = indices[:mid]
    val_idx = indices[mid:]

    disc_lats = all_lats[disc_idx]
    disc_lons = all_lons[disc_idx]
    val_lats = all_lats[val_idx]
    val_lons = all_lons[val_idx]

    # --- Find optimal pole on Discovery half (vectorized scan) ---
    best_count = -1
    best_pole_idx = 0
    for pi in range(N_POLES):
        c = count_within_vec(disc_lats, disc_lons, pole_lats[pi], pole_lons[pi], THRESHOLD_KM)
        if c > best_count:
            best_count = c
            best_pole_idx = pi

    best_plat = float(pole_lats[best_pole_idx])
    best_plon = float(pole_lons[best_pole_idx])

    # --- Test Discovery-optimal pole on Validation half ---
    disc_z, disc_obs, disc_mu, disc_sigma = mc_zscore_vec(
        val_lats, val_lons, best_plat, best_plon, THRESHOLD_KM, N_MC_TRIALS)

    # --- Test Alison pole on Validation half ---
    alison_z, alison_obs, alison_mu, alison_sigma = mc_zscore_vec(
        val_lats, val_lons, POLE_LAT, POLE_LNG, THRESHOLD_KM, N_MC_TRIALS)

    elapsed = time.time() - t_split
    split_results.append({
        "split": split_idx + 1,
        "discovery_n": int(len(disc_idx)),
        "validation_n": int(len(val_idx)),
        "discovery_optimal_pole": {"lat": best_plat, "lon": best_plon},
        "discovery_optimal_count": best_count,
        "discovery_optimal_z_on_validation": round(disc_z, 2),
        "discovery_optimal_obs_on_validation": disc_obs,
        "discovery_optimal_expected": round(disc_mu, 2),
        "alison_z_on_validation": round(alison_z, 2),
        "alison_obs_on_validation": alison_obs,
        "alison_expected": round(alison_mu, 2),
        "elapsed_sec": round(elapsed, 1)
    })

    print(f"  Split {split_idx+1}/{N_SPLITS}: "
          f"Disc-opt pole=({best_plat},{best_plon}) Z_val={disc_z:.2f} | "
          f"Alison Z_val={alison_z:.2f} [{elapsed:.1f}s]")

total_a = time.time() - t_start
print(f"\nPart A complete in {total_a:.0f}s")

# Compute summary statistics
disc_zs = [r["discovery_optimal_z_on_validation"] for r in split_results]
alison_zs = [r["alison_z_on_validation"] for r in split_results]
disc_poles = [r["discovery_optimal_pole"] for r in split_results]

summary_a = {
    "n_splits": N_SPLITS,
    "threshold_km": THRESHOLD_KM,
    "mc_trials_per_test": N_MC_TRIALS,
    "pole_scan_resolution_deg": POLE_SCAN_RES,
    "n_pole_candidates": N_POLES,
    "total_sites": N_SITES,
    "discovery_optimal_z_on_validation": {
        "mean": round(float(np.mean(disc_zs)), 2),
        "std": round(float(np.std(disc_zs)), 2),
        "min": round(float(np.min(disc_zs)), 2),
        "max": round(float(np.max(disc_zs)), 2),
        "median": round(float(np.median(disc_zs)), 2),
        "pct_above_3": round(sum(1 for z in disc_zs if z > 3) / len(disc_zs) * 100, 1),
        "pct_above_5": round(sum(1 for z in disc_zs if z > 5) / len(disc_zs) * 100, 1),
    },
    "alison_z_on_validation": {
        "mean": round(float(np.mean(alison_zs)), 2),
        "std": round(float(np.std(alison_zs)), 2),
        "min": round(float(np.min(alison_zs)), 2),
        "max": round(float(np.max(alison_zs)), 2),
        "median": round(float(np.median(alison_zs)), 2),
        "pct_above_3": round(sum(1 for z in alison_zs if z > 3) / len(alison_zs) * 100, 1),
        "pct_above_5": round(sum(1 for z in alison_zs if z > 5) / len(alison_zs) * 100, 1),
    },
    "discovery_optimal_pole_locations": {
        "unique_poles": len(set((p["lat"], p["lon"]) for p in disc_poles)),
        "most_common": None
    },
    "elapsed_seconds": round(total_a, 1)
}

# Most common Discovery-optimal pole
pole_counts = Counter((p["lat"], p["lon"]) for p in disc_poles)
most_common_pole, mc_count = pole_counts.most_common(1)[0]

def haversine_km_scalar(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1; dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))

summary_a["discovery_optimal_pole_locations"]["most_common"] = {
    "pole": {"lat": most_common_pole[0], "lon": most_common_pole[1]},
    "count": mc_count,
    "distance_from_alison_km": round(haversine_km_scalar(
        most_common_pole[0], most_common_pole[1], POLE_LAT, POLE_LNG), 1)
}

# Save Part A results
part_a_output = {"summary": summary_a, "splits": split_results}
with open(os.path.join(OUT_DIR, "split_results.json"), "w") as f:
    json.dump(part_a_output, f, indent=2)
print(f"Saved: split_results.json")

# Print Part A summary
print("\n" + "-" * 50)
print("PART A SUMMARY")
print("-" * 50)
sa = summary_a["discovery_optimal_z_on_validation"]
print(f"Discovery-optimal Z on validation:")
print(f"  Mean={sa['mean']}, Std={sa['std']}, Min={sa['min']}, Max={sa['max']}")
print(f"  {sa['pct_above_3']}% above Z=3, {sa['pct_above_5']}% above Z=5")
sb = summary_a["alison_z_on_validation"]
print(f"Alison Z on validation:")
print(f"  Mean={sb['mean']}, Std={sb['std']}, Min={sb['min']}, Max={sb['max']}")
print(f"  {sb['pct_above_3']}% above Z=3, {sb['pct_above_5']}% above Z=5")

# ============================================================
# PART B: "15 FAMOUS SITES" CONTROL
# ============================================================
print("\n" + "=" * 70)
print(f"PART B: {N_FAMOUS_TRIALS} RANDOM-15-SITE CIRCLE FITS")
print("=" * 70)

def fit_optimal_pole(subset_lats, subset_lons):
    """Find the pole that minimizes sum of squared distances from each site to the great circle."""
    def objective(params):
        plat, plon = params
        dists = dist_from_gc_vec(subset_lats, subset_lons, plat, plon)
        return float(np.sum(dists ** 2))

    best_result = None
    for start_lat in range(-60, 90, 30):
        for start_lon in range(-180, 180, 60):
            try:
                res = minimize(objective, [start_lat, start_lon],
                             method='Nelder-Mead',
                             options={'maxiter': 500, 'xatol': 0.5, 'fatol': 100})
                if best_result is None or res.fun < best_result.fun:
                    best_result = res
            except:
                pass

    if best_result is not None:
        plat = max(-90, min(90, best_result.x[0]))
        plon = ((best_result.x[1] + 180) % 360) - 180
        return plat, plon
    return 0.0, 0.0

# Full-database Alison Z for reference
alison_full_z, alison_full_obs, _, _ = mc_zscore_vec(
    all_lats, all_lons, POLE_LAT, POLE_LNG, THRESHOLD_KM, N_MC_TRIALS)
print(f"Alison full-database Z = {alison_full_z:.2f} (observed={alison_full_obs})")

famous_results = []
t_start = time.time()

for trial in range(N_FAMOUS_TRIALS):
    # Pick 15 random sites
    subset_idx = np.random.choice(N_SITES, N_FAMOUS_SITES, replace=False)
    subset_lats = all_lats[subset_idx]
    subset_lons = all_lons[subset_idx]

    # Fit optimal circle to these 15
    fit_pole_lat, fit_pole_lon = fit_optimal_pole(subset_lats, subset_lons)

    # Test fitted circle against remaining sites
    mask = np.ones(N_SITES, dtype=bool)
    mask[subset_idx] = False
    rem_lats = all_lats[mask]
    rem_lons = all_lons[mask]

    fit_z, fit_obs, fit_mu, fit_sigma = mc_zscore_vec(
        rem_lats, rem_lons, fit_pole_lat, fit_pole_lon, THRESHOLD_KM, N_MC_TRIALS)

    famous_results.append({
        "trial": trial + 1,
        "fit_pole": {"lat": round(fit_pole_lat, 2), "lon": round(fit_pole_lon, 2)},
        "z_on_remaining": round(fit_z, 2),
        "observed": fit_obs,
        "expected": round(fit_mu, 2)
    })

    if (trial + 1) % 50 == 0:
        elapsed = time.time() - t_start
        rate = (trial + 1) / elapsed
        eta = (N_FAMOUS_TRIALS - trial - 1) / rate
        print(f"  Trial {trial+1}/{N_FAMOUS_TRIALS}: Z={fit_z:.2f} "
              f"[{elapsed:.0f}s elapsed, ~{eta:.0f}s remaining]")

total_b = time.time() - t_start

# Analysis
famous_zs = [r["z_on_remaining"] for r in famous_results]
alison_exceeds = sum(1 for z in famous_zs if z >= alison_full_z)
alison_pctile = (1 - alison_exceeds / len(famous_zs)) * 100

summary_b = {
    "n_trials": N_FAMOUS_TRIALS,
    "n_sites_per_trial": N_FAMOUS_SITES,
    "alison_full_z": round(alison_full_z, 2),
    "random_15_z_distribution": {
        "mean": round(float(np.mean(famous_zs)), 2),
        "std": round(float(np.std(famous_zs)), 2),
        "min": round(float(np.min(famous_zs)), 2),
        "max": round(float(np.max(famous_zs)), 2),
        "median": round(float(np.median(famous_zs)), 2),
        "pct_5": round(float(np.percentile(famous_zs, 5)), 2),
        "pct_95": round(float(np.percentile(famous_zs, 95)), 2),
    },
    "alison_percentile": round(alison_pctile, 2),
    "n_random_exceeding_alison": alison_exceeds,
    "elapsed_seconds": round(total_b, 1)
}

part_b_output = {"summary": summary_b, "trials": famous_results}
with open(os.path.join(OUT_DIR, "famous_15_control.json"), "w") as f:
    json.dump(part_b_output, f, indent=2)
print(f"\nSaved: famous_15_control.json")

# Print Part B summary
zd = summary_b["random_15_z_distribution"]
print("\n" + "-" * 50)
print("PART B SUMMARY")
print("-" * 50)
print(f"Alison full-database Z = {alison_full_z:.2f}")
print(f"Random-15-site Z distribution:")
print(f"  Mean={zd['mean']}, Std={zd['std']}, Min={zd['min']}, Max={zd['max']}")
print(f"  5th pctile={zd['pct_5']}, 95th pctile={zd['pct_95']}")
print(f"Alison percentile: {alison_pctile:.2f}%")
print(f"Random circles exceeding Alison: {alison_exceeds}/{N_FAMOUS_TRIALS}")

# ============================================================
# PASS/FAIL ASSESSMENT
# ============================================================
print("\n" + "=" * 70)
print("PASS / FAIL ASSESSMENT")
print("=" * 70)

alison_val_mean = summary_a["alison_z_on_validation"]["mean"]
alison_val_pct5 = summary_a["alison_z_on_validation"]["pct_above_5"]
famous_pctile = summary_b["alison_percentile"]

if alison_val_mean > 5 and alison_val_pct5 > 80 and famous_pctile > 95:
    verdict = "PASS"
    detail = (f"Alison Z on validation consistently > 5 (mean={alison_val_mean:.2f}, "
              f"{alison_val_pct5}% above Z=5). "
              f"Alison exceeds {famous_pctile:.1f}% of random-15-site circles.")
elif alison_val_mean < 3:
    verdict = "FAIL"
    detail = f"Alison Z on validation drops below 3 (mean={alison_val_mean:.2f})"
else:
    verdict = "PARTIAL"
    detail = (f"Validation Z is significant but reduced (mean={alison_val_mean:.2f}). "
              f"This indicates original Z was partially inflated by circularity, "
              f"but the underlying signal survives blinded validation.")

assessment = {
    "verdict": verdict,
    "detail": detail,
    "criteria": {
        "alison_validation_mean_z": alison_val_mean,
        "pct_splits_above_z5": alison_val_pct5,
        "alison_percentile_in_famous15": famous_pctile,
    }
}

print(f"\n  VERDICT: {verdict}")
print(f"  {detail}")

# ============================================================
# SAVE RESULTS.md
# ============================================================
mc_pole = summary_a['discovery_optimal_pole_locations']['most_common']
results_md = f"""# Directive 1: Split-Sample Blinded Validation — Results

## Verdict: {verdict}
{detail}

## Part A: {N_SPLITS}-Split Cross-Validation

| Metric | Discovery-Optimal | Alison |
|--------|-------------------|--------|
| Mean Z on validation | {sa['mean']} | {sb['mean']} |
| Std Z | {sa['std']} | {sb['std']} |
| Min Z | {sa['min']} | {sb['min']} |
| Max Z | {sa['max']} | {sb['max']} |
| % above Z=3 | {sa['pct_above_3']}% | {sb['pct_above_3']}% |
| % above Z=5 | {sa['pct_above_5']}% | {sb['pct_above_5']}% |

Most common Discovery-optimal pole: ({mc_pole['pole']['lat']}, {mc_pole['pole']['lon']})
Distance from Alison pole: {mc_pole['distance_from_alison_km']} km
Appeared in {mc_pole['count']}/{N_SPLITS} splits

## Part B: 15-Famous-Sites Control ({N_FAMOUS_TRIALS} trials)

Alison full-database Z = {summary_b['alison_full_z']}

Random-15-site Z distribution:
- Mean: {zd['mean']}, Std: {zd['std']}
- Range: [{zd['min']}, {zd['max']}]
- 5th-95th percentile: [{zd['pct_5']}, {zd['pct_95']}]

**Alison percentile: {alison_pctile:.2f}%** ({alison_exceeds}/{N_FAMOUS_TRIALS} random circles exceed Alison's Z)

## Configuration
- Total sites: {N_SITES}
- Threshold: {THRESHOLD_KM} km
- MC trials per test: {N_MC_TRIALS}
- Pole scan resolution: {POLE_SCAN_RES} deg
- Part A elapsed: {summary_a['elapsed_seconds']}s
- Part B elapsed: {summary_b['elapsed_seconds']}s
"""

with open(os.path.join(OUT_DIR, "RESULTS.md"), "w") as f:
    f.write(results_md)

combined = {
    "assessment": assessment,
    "part_a_summary": summary_a,
    "part_b_summary": summary_b,
}
with open(os.path.join(OUT_DIR, "combined_summary.json"), "w") as f:
    json.dump(combined, f, indent=2)

print(f"\nAll results saved to {OUT_DIR}/")
print("Done.")
