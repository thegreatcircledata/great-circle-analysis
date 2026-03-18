#!/usr/bin/env python3
"""
Directive 2: Monument-Settlement Divergence on 10,000 Random Circles
=====================================================================
Tests whether the Alison circle's monument-settlement divergence is unique
among a large population of random great circles.

Strategy: Use 50 MC trials per circle for the 10,000-circle scan, then
re-run with 200 trials on any circles that come close to Alison's divergence.
"""

import csv, math, random, json, os, sys, time
import numpy as np
from scipy import stats as sp_stats

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
N_RANDOM_CIRCLES = 10000
N_MC_FAST = 50       # fast MC for scan
N_MC_PRECISE = 200   # precise MC for close contenders
ALISON_DIVERGENCE_ANCIENT = 12.78  # monument_Z - settlement_Z (pre-2000 BCE)
ALISON_DIVERGENCE_ALL = 7.36       # all-periods

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "divergence_10k")
os.makedirs(OUT_DIR, exist_ok=True)

# Pleiades type classifications
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

# ============================================================
# VECTORIZED HELPER FUNCTIONS
# ============================================================
def dist_from_gc_vec(site_lats, site_lons, pole_lat, pole_lon):
    """Vectorized distance from great circle for numpy arrays."""
    lat1r = np.radians(pole_lat)
    lon1r = np.radians(pole_lon)
    lat2r = np.radians(site_lats)
    lon2r = np.radians(site_lons)
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat/2)**2 + np.cos(lat1r)*np.cos(lat2r)*np.sin(dlon/2)**2
    d = EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QUARTER_CIRC)

def count_within_vec(site_lats, site_lons, pole_lat, pole_lon):
    dists = dist_from_gc_vec(site_lats, site_lons, pole_lat, pole_lon)
    return int(np.sum(dists <= THRESHOLD_KM))

def mc_zscore_vec(site_lats, site_lons, pole_lat, pole_lon, n_trials):
    """Vectorized MC Z-score."""
    n = len(site_lats)
    observed = count_within_vec(site_lats, site_lons, pole_lat, pole_lon)

    rand_counts = np.empty(n_trials)
    for t in range(n_trials):
        idx = np.random.randint(0, n, size=n)
        rlats = np.clip(site_lats[idx] + np.random.normal(0, 2, n), -90, 90)
        rlons = np.clip(site_lons[idx] + np.random.normal(0, 2, n), -180, 180)
        rand_counts[t] = count_within_vec(rlats, rlons, pole_lat, pole_lon)

    mu = np.mean(rand_counts)
    sigma = np.std(rand_counts)
    z = (observed - mu) / sigma if sigma > 0 else 0
    return float(z)

def random_pole_matched(all_lats, all_lons):
    """Generate a distribution-matched random pole."""
    lat = np.random.choice(all_lats) + np.random.normal(0, 10)
    lon = np.random.choice(all_lons) + np.random.normal(0, 10)
    return max(-90.0, min(90.0, lat)), ((lon + 180) % 360) - 180

# ============================================================
# LOAD PLEIADES DATA
# ============================================================
print("=" * 70)
print("LOADING PLEIADES DATA")
print("=" * 70)

pleiades_path = os.path.join(BASE_DIR, "pleiades-places-latest.csv")
mon_all_list = []
set_all_list = []
mon_anc_list = []
set_anc_list = []

with open(pleiades_path, encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            lat = float(row["reprLat"])
            lon = float(row["reprLong"])
            if lat == 0 and lon == 0: continue
            feature_types = {t.strip() for t in row.get("featureTypes", "").split(",")}
            is_monumental = bool(feature_types & MONUMENTAL_TYPES)
            is_settlement = bool(feature_types & SETTLEMENT_TYPES)

            try:
                min_date = int(row.get("minDate", "0"))
            except (ValueError, TypeError):
                min_date = 0

            if is_monumental:
                mon_all_list.append((lat, lon))
                if min_date < -2000:
                    mon_anc_list.append((lat, lon))
            elif is_settlement:
                set_all_list.append((lat, lon))
                if min_date < -2000:
                    set_anc_list.append((lat, lon))
        except (ValueError, KeyError, TypeError):
            pass

# Convert to numpy arrays
mon_all_arr = np.array(mon_all_list)
set_all_arr = np.array(set_all_list)
mon_anc_arr = np.array(mon_anc_list)
set_anc_arr = np.array(set_anc_list)

mon_all_lats, mon_all_lons = mon_all_arr[:, 0], mon_all_arr[:, 1]
set_all_lats, set_all_lons = set_all_arr[:, 0], set_all_arr[:, 1]
mon_anc_lats, mon_anc_lons = mon_anc_arr[:, 0], mon_anc_arr[:, 1]
set_anc_lats, set_anc_lons = set_anc_arr[:, 0], set_anc_arr[:, 1]

print(f"Monumental (all): {len(mon_all_arr)}")
print(f"Settlement (all): {len(set_all_arr)}")
print(f"Monumental (ancient, pre-2000 BCE): {len(mon_anc_arr)}")
print(f"Settlement (ancient, pre-2000 BCE): {len(set_anc_arr)}")

# Combined for random pole generation
all_pleiades_lats = np.concatenate([mon_all_lats, set_all_lats])
all_pleiades_lons = np.concatenate([mon_all_lons, set_all_lons])

# ============================================================
# VERIFY ALISON DIVERGENCE
# ============================================================
print("\n--- Verifying Alison divergence values ---")

alison_mon_ancient_z = mc_zscore_vec(mon_anc_lats, mon_anc_lons, POLE_LAT, POLE_LNG, N_MC_PRECISE)
alison_set_ancient_z = mc_zscore_vec(set_anc_lats, set_anc_lons, POLE_LAT, POLE_LNG, N_MC_PRECISE)
alison_div_ancient = alison_mon_ancient_z - alison_set_ancient_z

alison_mon_all_z = mc_zscore_vec(mon_all_lats, mon_all_lons, POLE_LAT, POLE_LNG, N_MC_PRECISE)
alison_set_all_z = mc_zscore_vec(set_all_lats, set_all_lons, POLE_LAT, POLE_LNG, N_MC_PRECISE)
alison_div_all = alison_mon_all_z - alison_set_all_z

print(f"Alison ancient: monument Z={alison_mon_ancient_z:.2f}, settlement Z={alison_set_ancient_z:.2f}, "
      f"divergence={alison_div_ancient:.2f} (reference: {ALISON_DIVERGENCE_ANCIENT})")
print(f"Alison all-period: monument Z={alison_mon_all_z:.2f}, settlement Z={alison_set_all_z:.2f}, "
      f"divergence={alison_div_all:.2f} (reference: {ALISON_DIVERGENCE_ALL})")

# ============================================================
# 10,000 RANDOM CIRCLES — FAST SCAN
# ============================================================
print("\n" + "=" * 70)
print(f"SCANNING {N_RANDOM_CIRCLES} RANDOM CIRCLES (fast MC, {N_MC_FAST} trials each)")
print("=" * 70)

random_divergences_ancient = []
random_divergences_all = []
circle_details = []
t_start = time.time()

for i in range(N_RANDOM_CIRCLES):
    # Random pole
    plat, plon = random_pole_matched(all_pleiades_lats, all_pleiades_lons)

    # Ancient divergence
    mon_z = mc_zscore_vec(mon_anc_lats, mon_anc_lons, plat, plon, N_MC_FAST)
    set_z = mc_zscore_vec(set_anc_lats, set_anc_lons, plat, plon, N_MC_FAST)
    div_ancient = mon_z - set_z

    # All-periods divergence
    mon_all_z = mc_zscore_vec(mon_all_lats, mon_all_lons, plat, plon, N_MC_FAST)
    set_all_z = mc_zscore_vec(set_all_lats, set_all_lons, plat, plon, N_MC_FAST)
    div_all = mon_all_z - set_all_z

    random_divergences_ancient.append(div_ancient)
    random_divergences_all.append(div_all)

    circle_details.append({
        "pole": {"lat": round(plat, 2), "lon": round(plon, 2)},
        "ancient": {"mon_z": round(mon_z, 2), "set_z": round(set_z, 2), "divergence": round(div_ancient, 2)},
        "all_periods": {"mon_z": round(mon_all_z, 2), "set_z": round(set_all_z, 2), "divergence": round(div_all, 2)}
    })

    if (i + 1) % 500 == 0:
        elapsed = time.time() - t_start
        rate = (i + 1) / elapsed
        eta = (N_RANDOM_CIRCLES - i - 1) / rate
        max_anc = max(random_divergences_ancient)
        print(f"  Circle {i+1}/{N_RANDOM_CIRCLES}: "
              f"max_ancient_div={max_anc:.2f} "
              f"[{elapsed:.0f}s elapsed, ~{eta:.0f}s remaining]")

total_scan = time.time() - t_start
print(f"\nScan complete in {total_scan:.0f}s")

# ============================================================
# RE-RUN CLOSE CONTENDERS WITH PRECISE MC
# ============================================================
print("\n--- Re-running close contenders with precise MC ---")

# Find circles with ancient divergence > 0.7 * Alison
threshold_rerun = alison_div_ancient * 0.7
close_indices = [i for i, d in enumerate(random_divergences_ancient) if d >= threshold_rerun]
print(f"Circles within 70% of Alison's divergence: {len(close_indices)}")

for idx in close_indices:
    plat = circle_details[idx]["pole"]["lat"]
    plon = circle_details[idx]["pole"]["lon"]

    # Re-run with precise MC
    mon_z = mc_zscore_vec(mon_anc_lats, mon_anc_lons, plat, plon, N_MC_PRECISE)
    set_z = mc_zscore_vec(set_anc_lats, set_anc_lons, plat, plon, N_MC_PRECISE)
    div_ancient = mon_z - set_z

    # Update
    old_div = random_divergences_ancient[idx]
    random_divergences_ancient[idx] = div_ancient
    circle_details[idx]["ancient"] = {
        "mon_z": round(mon_z, 2), "set_z": round(set_z, 2),
        "divergence": round(div_ancient, 2), "precise": True
    }
    print(f"  Circle {idx}: {old_div:.2f} -> {div_ancient:.2f} (precise)")

# ============================================================
# COMPUTE PERCENTILES
# ============================================================
print("\n" + "=" * 70)
print("COMPUTING PERCENTILES")
print("=" * 70)

div_anc = np.array(random_divergences_ancient)
div_all_arr = np.array(random_divergences_all)

# Ancient divergence
n_exceeding_ancient = int(np.sum(div_anc >= alison_div_ancient))
pctile_ancient = (1 - n_exceeding_ancient / N_RANDOM_CIRCLES) * 100

# All-periods divergence
n_exceeding_all = int(np.sum(div_all_arr >= alison_div_all))
pctile_all = (1 - n_exceeding_all / N_RANDOM_CIRCLES) * 100

percentile_results = {
    "ancient_pre2000bce": {
        "alison_divergence": round(alison_div_ancient, 2),
        "alison_monument_z": round(alison_mon_ancient_z, 2),
        "alison_settlement_z": round(alison_set_ancient_z, 2),
        "n_random_circles": N_RANDOM_CIRCLES,
        "n_exceeding_alison": n_exceeding_ancient,
        "alison_percentile": round(pctile_ancient, 4),
        "distribution": {
            "mean": round(float(np.mean(div_anc)), 2),
            "median": round(float(np.median(div_anc)), 2),
            "std": round(float(np.std(div_anc)), 2),
            "min": round(float(np.min(div_anc)), 2),
            "max": round(float(np.max(div_anc)), 2),
            "pct_5": round(float(np.percentile(div_anc, 5)), 2),
            "pct_25": round(float(np.percentile(div_anc, 25)), 2),
            "pct_75": round(float(np.percentile(div_anc, 75)), 2),
            "pct_95": round(float(np.percentile(div_anc, 95)), 2),
            "pct_99": round(float(np.percentile(div_anc, 99)), 2),
        }
    },
    "all_periods": {
        "alison_divergence": round(alison_div_all, 2),
        "alison_monument_z": round(alison_mon_all_z, 2),
        "alison_settlement_z": round(alison_set_all_z, 2),
        "n_random_circles": N_RANDOM_CIRCLES,
        "n_exceeding_alison": n_exceeding_all,
        "alison_percentile": round(pctile_all, 4),
        "distribution": {
            "mean": round(float(np.mean(div_all_arr)), 2),
            "median": round(float(np.median(div_all_arr)), 2),
            "std": round(float(np.std(div_all_arr)), 2),
            "min": round(float(np.min(div_all_arr)), 2),
            "max": round(float(np.max(div_all_arr)), 2),
            "pct_5": round(float(np.percentile(div_all_arr, 5)), 2),
            "pct_25": round(float(np.percentile(div_all_arr, 25)), 2),
            "pct_75": round(float(np.percentile(div_all_arr, 75)), 2),
            "pct_95": round(float(np.percentile(div_all_arr, 95)), 2),
            "pct_99": round(float(np.percentile(div_all_arr, 99)), 2),
        }
    }
}

print(f"\nAncient (pre-2000 BCE):")
print(f"  Alison divergence: {alison_div_ancient:.2f}")
print(f"  Random distribution: mean={np.mean(div_anc):.2f}, std={np.std(div_anc):.2f}, max={np.max(div_anc):.2f}")
print(f"  Exceeding Alison: {n_exceeding_ancient}/{N_RANDOM_CIRCLES}")
print(f"  Alison percentile: {pctile_ancient:.4f}%")

print(f"\nAll periods:")
print(f"  Alison divergence: {alison_div_all:.2f}")
print(f"  Random distribution: mean={np.mean(div_all_arr):.2f}, std={np.std(div_all_arr):.2f}, max={np.max(div_all_arr):.2f}")
print(f"  Exceeding Alison: {n_exceeding_all}/{N_RANDOM_CIRCLES}")
print(f"  Alison percentile: {pctile_all:.4f}%")

# Confidence interval on percentile (binomial)
if n_exceeding_ancient == 0:
    ci_upper_ancient = 3 / N_RANDOM_CIRCLES * 100  # rule of 3
    ci_note_ancient = f"0/{N_RANDOM_CIRCLES} exceeded; 95% CI upper bound: {ci_upper_ancient:.4f}%"
else:
    se = math.sqrt(pctile_ancient * (100 - pctile_ancient) / N_RANDOM_CIRCLES)
    ci_note_ancient = f"SE of percentile: {se:.2f}%"

percentile_results["ancient_pre2000bce"]["ci_note"] = ci_note_ancient

# ============================================================
# SAVE OUTPUTS
# ============================================================
with open(os.path.join(OUT_DIR, "alison_percentile.json"), "w") as f:
    json.dump(percentile_results, f, indent=2)

# Save all divergences (compact — just values, not full circle details)
with open(os.path.join(OUT_DIR, "random_divergences.json"), "w") as f:
    json.dump({
        "ancient_divergences": [round(d, 2) for d in random_divergences_ancient],
        "all_period_divergences": [round(d, 2) for d in random_divergences_all],
        "n_circles": N_RANDOM_CIRCLES,
        "mc_trials_fast": N_MC_FAST,
        "mc_trials_precise": N_MC_PRECISE,
        "n_rerun_precise": len(close_indices),
        "elapsed_seconds": round(total_scan, 1)
    }, f, indent=2)

# RESULTS.md
da = percentile_results["ancient_pre2000bce"]
dp = percentile_results["all_periods"]

results_md = f"""# Directive 2: Monument-Settlement Divergence on {N_RANDOM_CIRCLES} Random Circles — Results

## Ancient Sites (pre-2000 BCE)

Alison monument Z: {da['alison_monument_z']}
Alison settlement Z: {da['alison_settlement_z']}
**Alison divergence: {da['alison_divergence']}**

Random circle divergence distribution ({N_RANDOM_CIRCLES} circles):
- Mean: {da['distribution']['mean']}, Std: {da['distribution']['std']}
- Median: {da['distribution']['median']}
- Range: [{da['distribution']['min']}, {da['distribution']['max']}]
- 95th percentile: {da['distribution']['pct_95']}
- 99th percentile: {da['distribution']['pct_99']}

**Alison percentile: {da['alison_percentile']}%** ({da['n_exceeding_alison']}/{N_RANDOM_CIRCLES} random circles exceed Alison)
{ci_note_ancient}

## All-Period Sites

Alison monument Z: {dp['alison_monument_z']}
Alison settlement Z: {dp['alison_settlement_z']}
**Alison divergence: {dp['alison_divergence']}**

Random circle divergence distribution:
- Mean: {dp['distribution']['mean']}, Std: {dp['distribution']['std']}
- Range: [{dp['distribution']['min']}, {dp['distribution']['max']}]
- 95th percentile: {dp['distribution']['pct_95']}
- 99th percentile: {dp['distribution']['pct_99']}

**Alison percentile: {dp['alison_percentile']}%** ({dp['n_exceeding_alison']}/{N_RANDOM_CIRCLES} random circles exceed Alison)

## Configuration
- MC trials (fast scan): {N_MC_FAST}
- MC trials (precise re-run): {N_MC_PRECISE}
- Re-run threshold: 70% of Alison divergence
- Circles re-run with precise MC: {len(close_indices)}
- Total scan time: {total_scan:.0f}s
"""

with open(os.path.join(OUT_DIR, "RESULTS.md"), "w") as f:
    f.write(results_md)

print(f"\nAll results saved to {OUT_DIR}/")
print("Done.")
