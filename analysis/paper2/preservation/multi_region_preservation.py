#!/usr/bin/env python3
"""
Multi-Region Preservation Simulation (Directive 3)
====================================================
Extends the Egypt preservation test to Mesopotamia, Indus Valley, and Peru.

Question: Does the preservation argument hold beyond Egypt?
Hypothesis: Adding estimated missing settlements in riverine floodplains
should INCREASE divergence D because buried settlements cluster far from
the great circle (which clips desert/mountain edges, not floodplain cores).

Method:
  1. Per-region analysis: generate estimated missing settlements with
     region-appropriate Gaussian distributions, measure how many fall
     within 50km of the great circle.
  2. Combined Monte Carlo: add ALL estimated missing settlements
     (Egypt +350, Mesopotamia +200, Indus +500, Peru +200 = +1,250)
     and recompute D across 1,000 iterations.
  3. Compare Egypt-only vs multi-region preservation results.

Uses the GLOBAL Pleiades dataset (not Egypt-only) for the divergence
computation, matching the main paper's methodology.
"""

import csv, math, random, json, os, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
THRESHOLD_KM = 50
N_MC_TRIALS = 200           # MC trials per Z-score (matches main paper default)
N_MONTE_CARLO_ITER = 1000   # outer MC for combined test

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "extended_controls")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# VECTORIZED HELPERS (same as egypt_preservation_test)
# ============================================================
def haversine_vec(lat1, lon1, lat2_arr, lon2_arr):
    lat1_r, lon1_r = np.radians(lat1), np.radians(lon1)
    lat2_r, lon2_r = np.radians(lat2_arr), np.radians(lon2_arr)
    dlat = lat2_r - lat1_r; dlon = lon2_r - lon1_r
    a = np.sin(dlat/2)**2 + np.cos(lat1_r)*np.cos(lat2_r)*np.sin(dlon/2)**2
    return EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))

def gc_dist_vec(pole_lat, pole_lon, site_lats, site_lons):
    d = haversine_vec(pole_lat, pole_lon, site_lats, site_lons)
    return np.abs(d - QUARTER_CIRC)

# ============================================================
# CLASSIFICATION (global — matching alison_pole_divergence_scan.py)
# ============================================================
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
# MC Z-SCORE ENGINE
# ============================================================
def run_mc(site_lats, site_lons, threshold=THRESHOLD_KM, n_trials=N_MC_TRIALS):
    n = len(site_lats)
    if n < 5:
        return {"n_sites": int(n), "too_few": True}
    dists = gc_dist_vec(POLE_LAT, POLE_LON, site_lats, site_lons)
    observed = int(np.sum(dists <= threshold))
    rand_counts = []
    for _ in range(n_trials):
        idx = np.random.randint(0, n, size=n)
        rl = site_lats[idx] + np.random.normal(0, 2, n)
        rn = site_lons[idx] + np.random.normal(0, 2, n)
        rl = np.clip(rl, -90, 90)
        rn = np.clip(rn, -180, 180)
        rd = gc_dist_vec(POLE_LAT, POLE_LON, rl, rn)
        rand_counts.append(int(np.sum(rd <= threshold)))
    mu = np.mean(rand_counts); sigma = np.std(rand_counts)
    z = float((observed - mu) / sigma) if sigma > 0 else 0.0
    enrich = float(observed / mu) if mu > 0 else 0.0
    p_val = float(np.sum(np.array(rand_counts) >= observed) / n_trials)
    return {
        "n_sites": int(n), "observed": int(observed),
        "expected": round(float(mu), 2), "std": round(float(sigma), 2),
        "z_score": round(z, 2), "enrichment": round(enrich, 3),
        "p_value": round(p_val, 4)
    }

def compute_divergence(mono_res, settle_res):
    if mono_res.get("too_few") or settle_res.get("too_few"):
        return "N/A"
    return round(mono_res["z_score"] - settle_res["z_score"], 2)

# ============================================================
# REGION DEFINITIONS — estimated missing settlements
# ============================================================
REGIONS = {
    "egypt": {
        "label": "Egypt (Nile floodplain)",
        "center": (30.5, 31.0),
        "std": (0.8, 0.5),
        "n_estimated": 350,
        "rationale": "70% Delta (30-31.5N), 30% Valley (24-30N). Most buried under modern cities and alluvium.",
        "custom_generator": True,  # uses rand_nile_floodplain()
    },
    "mesopotamia": {
        "label": "Mesopotamia (Tigris-Euphrates alluvial plain)",
        "center": (32.0, 45.5),
        "std": (0.8, 1.2),
        "n_estimated": 200,
        "n_sweep": [50, 100, 200, 500],
        "rationale": "Thousands of unexcavated tells under modern Iraqi cities. Circle passes ~40km from Ur. Settlements cluster along rivers; circle clips the desert edge.",
    },
    "indus": {
        "label": "Indus Valley (Mohenjo-daro region + Punjab)",
        "center": (27.5, 69.0),
        "std": (1.5, 2.0),
        "n_estimated": 500,
        "n_sweep": [100, 500, 1000],
        "rationale": "~2,600 known Indus sites, most in Punjab/Cholistan far from circle. Circle passes ~20km from Mohenjo-daro. City walls = monumental but most settlements are domestic.",
    },
    "peru": {
        "label": "Peru (Lima coastal plain)",
        "center": (-12.0, -77.0),
        "std": (1.0, 0.5),
        "n_estimated": 200,
        "n_sweep": [100, 200, 500],
        "rationale": "Thousands of huaca sites under modern Lima. Circle passes through southern Peru (~14S). Most unexcavated sites are coastal, 100-200km from circle.",
    },
}

# ============================================================
# EGYPT FLOODPLAIN GENERATOR (from egypt_preservation_test.py)
# ============================================================
def rand_nile_floodplain():
    if random.random() < 0.70:
        lat = random.uniform(30.0, 31.5)
        lon = random.uniform(29.5, 32.5)
        center_lon = 31.0
        half_width = 0.5 + (lat - 30.0) * 1.3
        if abs(lon - center_lon) > half_width:
            lon = center_lon + random.uniform(-half_width, half_width)
    else:
        lat = random.uniform(24.0, 30.0)
        if lat < 25.5:
            nile_lon = 32.8 + (lat - 24.0) * 0.05
        elif lat < 26.5:
            nile_lon = 32.6 - (lat - 25.5) * 0.8
        elif lat < 27.5:
            nile_lon = 31.8 - (lat - 26.5) * 0.5
        elif lat < 28.5:
            nile_lon = 31.3 - (lat - 27.5) * 0.2
        else:
            nile_lon = 31.1 - (lat - 28.5) * 0.1
        lon = nile_lon + random.uniform(-0.15, 0.15)
    return lat, lon

def generate_region_settlements(region_key, n):
    """Generate n random settlement positions for a region."""
    if region_key == "egypt":
        points = [rand_nile_floodplain() for _ in range(n)]
        return np.array([p[0] for p in points]), np.array([p[1] for p in points])
    else:
        r = REGIONS[region_key]
        lats = np.random.normal(r["center"][0], r["std"][0], n)
        lons = np.random.normal(r["center"][1], r["std"][1], n)
        lats = np.clip(lats, -90, 90)
        lons = np.clip(lons, -180, 180)
        return lats, lons

# ============================================================
# MAIN ANALYSIS
# ============================================================
print("=" * 70)
print("MULTI-REGION PRESERVATION SIMULATION")
print("Directive 3: Does the preservation argument hold beyond Egypt?")
print("=" * 70)

t0 = time.time()
np.random.seed(42)
random.seed(42)

# --- Load global Pleiades dataset ---
print("\n--- Loading GLOBAL Pleiades data ---")
pleiades_csv = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")

mon_list = []
set_list = []

with open(pleiades_csv, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            lat = float(row.get('reprLat', '') or '0')
            lon = float(row.get('reprLong', '') or '0')
        except ValueError:
            continue
        if lat == 0 and lon == 0:
            continue
        feature_types = {t.strip() for t in row.get("featureTypes", "").split(",")}
        is_mon = bool(feature_types & MONUMENTAL_TYPES)
        is_set = bool(feature_types & SETTLEMENT_TYPES)
        if is_mon:
            mon_list.append((lat, lon))
        elif is_set:
            set_list.append((lat, lon))

mon_lats = np.array([s[0] for s in mon_list])
mon_lons = np.array([s[1] for s in mon_list])
set_lats = np.array([s[0] for s in set_list])
set_lons = np.array([s[1] for s in set_list])

print(f"  Monumental sites: {len(mon_list)}")
print(f"  Settlement sites: {len(set_list)}")

# ============================================================
# TEST 0: GLOBAL BASELINE (no augmentation)
# ============================================================
print("\n" + "=" * 70)
print("TEST 0: GLOBAL BASELINE")
print("=" * 70)

baseline_mono = run_mc(mon_lats, mon_lons)
baseline_settle = run_mc(set_lats, set_lons)
baseline_D = compute_divergence(baseline_mono, baseline_settle)

print(f"  Monumental: Z={baseline_mono.get('z_score')}, obs={baseline_mono.get('observed')}/{baseline_mono.get('n_sites')}")
print(f"  Settlement: Z={baseline_settle.get('z_score')}, obs={baseline_settle.get('observed')}/{baseline_settle.get('n_sites')}")
print(f"  Divergence D = {baseline_D}")

# ============================================================
# PER-REGION SWEEP: distance distributions
# ============================================================
print("\n" + "=" * 70)
print("PER-REGION ANALYSIS: How many estimated settlements fall near the circle?")
print("=" * 70)

region_results = {}

for rkey, rdef in REGIONS.items():
    if rkey == "egypt":
        sweep_ns = [50, 100, 200, 350, 500]
    else:
        sweep_ns = rdef.get("n_sweep", [rdef["n_estimated"]])

    print(f"\n  --- {rdef['label']} ---")
    print(f"  Rationale: {rdef['rationale']}")

    sweep_results = []
    all_distances = []

    for n_est in sweep_ns:
        # Run 100 samples to get stable statistics
        near_counts = []
        mean_dists = []
        for _ in range(100):
            est_lats, est_lons = generate_region_settlements(rkey, n_est)
            dists = gc_dist_vec(POLE_LAT, POLE_LON, est_lats, est_lons)
            near = int(np.sum(dists <= THRESHOLD_KM))
            near_counts.append(near)
            mean_dists.append(float(np.mean(dists)))
            if n_est == rdef["n_estimated"]:
                all_distances.extend(dists.tolist())

        avg_near = np.mean(near_counts)
        avg_dist = np.mean(mean_dists)
        pct_near = avg_near / n_est * 100

        result = {
            "n_estimated": n_est,
            "avg_within_50km": round(float(avg_near), 1),
            "pct_within_50km": round(float(pct_near), 1),
            "avg_mean_distance_km": round(float(avg_dist), 0),
        }
        sweep_results.append(result)
        print(f"    +{n_est:>4d}: {avg_near:.1f} within 50km ({pct_near:.1f}%), mean dist = {avg_dist:.0f} km")

    region_results[rkey] = {
        "label": rdef["label"],
        "center": rdef["center"],
        "std": rdef["std"],
        "n_estimated": rdef["n_estimated"],
        "rationale": rdef["rationale"],
        "sweep": sweep_results,
        "distances_sample": all_distances[:1000],  # keep sample for plotting
    }

# ============================================================
# PLOT 1: Per-region distance distributions
# ============================================================
print("\n--- Generating per-region distance distribution plot ---")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("Distance from Great Circle: Estimated Missing Settlements by Region",
             fontsize=14, fontweight='bold')

colors = {'egypt': '#d4a574', 'mesopotamia': '#8B4513', 'indus': '#2E86AB', 'peru': '#A23B72'}
region_order = ['egypt', 'mesopotamia', 'indus', 'peru']

for ax, rkey in zip(axes.flat, region_order):
    dists = region_results[rkey]["distances_sample"]
    if len(dists) == 0:
        continue
    ax.hist(dists, bins=50, color=colors[rkey], alpha=0.7, edgecolor='black', linewidth=0.5)
    ax.axvline(x=THRESHOLD_KM, color='red', linestyle='--', linewidth=2, label=f'{THRESHOLD_KM}km threshold')
    ax.set_title(region_results[rkey]["label"], fontsize=12, fontweight='bold')
    ax.set_xlabel("Distance from Great Circle (km)")
    ax.set_ylabel("Count")

    n_near = sum(1 for d in dists if d <= THRESHOLD_KM)
    pct = n_near / len(dists) * 100
    ax.text(0.95, 0.95, f"{pct:.1f}% within {THRESHOLD_KM}km\n(n={len(dists)})",
            transform=ax.transAxes, ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=10)
    ax.legend(loc='upper left')

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "preservation_by_region.png"), dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: outputs/extended_controls/preservation_by_region.png")

# ============================================================
# COMBINED MULTI-REGION MONTE CARLO
# ============================================================
print("\n" + "=" * 70)
print(f"COMBINED MULTI-REGION PRESERVATION TEST")
print(f"Egypt +350, Mesopotamia +200, Indus +500, Peru +200 = +1,250 total")
print(f"{N_MONTE_CARLO_ITER} iterations, {N_MC_TRIALS} MC trials per Z-score")
print("=" * 70)

combined_D_values = []
combined_settle_z = []
combined_mono_z = []
combined_n_near = []

for i in range(N_MONTE_CARLO_ITER):
    if (i + 1) % 100 == 0:
        elapsed = time.time() - t0
        print(f"  Iteration {i+1}/{N_MONTE_CARLO_ITER}... ({elapsed:.0f}s elapsed)")

    # Generate estimated settlements for all regions
    all_est_lats = []
    all_est_lons = []

    for rkey in region_order:
        n = REGIONS[rkey]["n_estimated"]
        rl, rn = generate_region_settlements(rkey, n)
        all_est_lats.append(rl)
        all_est_lons.append(rn)

    est_lats = np.concatenate(all_est_lats)
    est_lons = np.concatenate(all_est_lons)

    # Count how many estimated settlements are near the circle
    est_dists = gc_dist_vec(POLE_LAT, POLE_LON, est_lats, est_lons)
    n_near = int(np.sum(est_dists <= THRESHOLD_KM))
    combined_n_near.append(n_near)

    # Augment settlement dataset
    aug_set_lats = np.concatenate([set_lats, est_lats])
    aug_set_lons = np.concatenate([set_lons, est_lons])

    # Compute Z-scores
    mono_res = run_mc(mon_lats, mon_lons)
    settle_res = run_mc(aug_set_lats, aug_set_lons)

    if not mono_res.get("too_few") and not settle_res.get("too_few"):
        D = mono_res['z_score'] - settle_res['z_score']
        combined_D_values.append(D)
        combined_settle_z.append(settle_res['z_score'])
        combined_mono_z.append(mono_res['z_score'])

combined_D_values = np.array(combined_D_values)
combined_settle_z = np.array(combined_settle_z)
combined_mono_z = np.array(combined_mono_z)
combined_n_near = np.array(combined_n_near)

mean_D = float(np.mean(combined_D_values))
std_D = float(np.std(combined_D_values))
median_D = float(np.median(combined_D_values))
pct_positive = float(np.sum(combined_D_values > 0) / len(combined_D_values) * 100)
pct_gt2 = float(np.sum(combined_D_values > 2) / len(combined_D_values) * 100)

print(f"\n  Results across {len(combined_D_values)} valid iterations:")
print(f"  Mean D: {mean_D:.2f} ± {std_D:.2f}")
print(f"  Median D: {median_D:.2f}")
print(f"  D > 0 in {pct_positive:.1f}% of iterations")
print(f"  D > 2 in {pct_gt2:.1f}% of iterations")
print(f"  Mean estimated settlements within 50km: {np.mean(combined_n_near):.1f} / 1250")

# ============================================================
# EGYPT-ONLY COMPARISON (run same methodology for fair comparison)
# ============================================================
print("\n" + "=" * 70)
print("EGYPT-ONLY COMPARISON (global Pleiades, +350 Egypt only)")
print("=" * 70)

egypt_only_D_values = []
egypt_only_settle_z = []
egypt_only_mono_z = []

for i in range(N_MONTE_CARLO_ITER):
    if (i + 1) % 100 == 0:
        elapsed = time.time() - t0
        print(f"  Iteration {i+1}/{N_MONTE_CARLO_ITER}... ({elapsed:.0f}s elapsed)")

    eg_lats, eg_lons = generate_region_settlements("egypt", 350)
    aug_set_lats = np.concatenate([set_lats, eg_lats])
    aug_set_lons = np.concatenate([set_lons, eg_lons])

    mono_res = run_mc(mon_lats, mon_lons)
    settle_res = run_mc(aug_set_lats, aug_set_lons)

    if not mono_res.get("too_few") and not settle_res.get("too_few"):
        D = mono_res['z_score'] - settle_res['z_score']
        egypt_only_D_values.append(D)
        egypt_only_settle_z.append(settle_res['z_score'])
        egypt_only_mono_z.append(mono_res['z_score'])

egypt_only_D_values = np.array(egypt_only_D_values)
egypt_only_settle_z = np.array(egypt_only_settle_z)
egypt_only_mono_z = np.array(egypt_only_mono_z)

eg_mean_D = float(np.mean(egypt_only_D_values))
eg_std_D = float(np.std(egypt_only_D_values))
eg_pct_gt2 = float(np.sum(egypt_only_D_values > 2) / len(egypt_only_D_values) * 100)

print(f"  Mean D: {eg_mean_D:.2f} ± {eg_std_D:.2f}")
print(f"  D > 2 in {eg_pct_gt2:.1f}% of iterations")

# ============================================================
# PLOT 2: Combined D distribution
# ============================================================
print("\n--- Generating combined MC D distribution plot ---")

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle("Divergence D Under Preservation Correction", fontsize=14, fontweight='bold')

# Panel 1: Baseline vs Egypt-only vs Multi-region
ax = axes[0]
ax.axvline(x=baseline_D, color='black', linewidth=2, linestyle='-', label=f'Baseline D={baseline_D}')
ax.hist(egypt_only_D_values, bins=40, color='#d4a574', alpha=0.6, edgecolor='black',
        linewidth=0.5, label=f'Egypt +350 (mean={eg_mean_D:.2f})', density=True)
ax.hist(combined_D_values, bins=40, color='#2E86AB', alpha=0.6, edgecolor='black',
        linewidth=0.5, label=f'All regions +1250 (mean={mean_D:.2f})', density=True)
ax.axvline(x=2, color='red', linestyle='--', linewidth=1.5, label='D=2 threshold')
ax.set_xlabel("Divergence D", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
ax.set_title("D Distribution Comparison", fontsize=12)
ax.legend(fontsize=9)

# Panel 2: Multi-region D with percentiles
ax = axes[1]
ax.hist(combined_D_values, bins=50, color='#2E86AB', alpha=0.7, edgecolor='black', linewidth=0.5)
ax.axvline(x=2, color='red', linestyle='--', linewidth=2, label='D=2 threshold')
ax.axvline(x=mean_D, color='navy', linestyle='-', linewidth=2, label=f'Mean={mean_D:.2f}')
p5 = np.percentile(combined_D_values, 5)
p95 = np.percentile(combined_D_values, 95)
ax.axvline(x=p5, color='navy', linestyle=':', linewidth=1, label=f'5th pctile={p5:.2f}')
ax.axvline(x=p95, color='navy', linestyle=':', linewidth=1, label=f'95th pctile={p95:.2f}')
ax.set_xlabel("Divergence D", fontsize=12)
ax.set_ylabel("Count", fontsize=12)
ax.set_title(f"Multi-Region +1250: D>{2} in {pct_gt2:.0f}%", fontsize=12)
ax.legend(fontsize=9)

# Panel 3: Settlement Z-score shift
ax = axes[2]
ax.hist(egypt_only_settle_z, bins=40, color='#d4a574', alpha=0.6, edgecolor='black',
        linewidth=0.5, label=f'Egypt +350 (mean Z={np.mean(egypt_only_settle_z):.2f})')
ax.hist(combined_settle_z, bins=40, color='#2E86AB', alpha=0.6, edgecolor='black',
        linewidth=0.5, label=f'All +1250 (mean Z={np.mean(combined_settle_z):.2f})')
ax.axvline(x=baseline_settle['z_score'], color='black', linewidth=2,
           label=f'Baseline settle Z={baseline_settle["z_score"]}')
ax.set_xlabel("Settlement Z-score", fontsize=12)
ax.set_ylabel("Density")
ax.set_title("Settlement Z-Score Under Augmentation", fontsize=12)
ax.legend(fontsize=9)

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "combined_preservation_mc.png"), dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: outputs/extended_controls/combined_preservation_mc.png")

# ============================================================
# SUMMARY
# ============================================================
elapsed = time.time() - t0
print(f"\n{'=' * 70}")
print(f"MULTI-REGION PRESERVATION SIMULATION — SUMMARY ({elapsed:.0f}s)")
print(f"{'=' * 70}")
print(f"  Baseline (no augmentation):         D = {baseline_D}")
print(f"  Egypt-only (+350, global Pleiades):  D = {eg_mean_D:.2f} ± {eg_std_D:.2f}  (D>2 in {eg_pct_gt2:.0f}%)")
print(f"  Multi-region (+1250):                D = {mean_D:.2f} ± {std_D:.2f}  (D>2 in {pct_gt2:.0f}%)")
print(f"")
print(f"  Per-region estimated settlements near circle:")
for rkey in region_order:
    sweep = region_results[rkey]["sweep"]
    main = [s for s in sweep if s["n_estimated"] == REGIONS[rkey]["n_estimated"]][0]
    print(f"    {REGIONS[rkey]['label']:45s} +{main['n_estimated']:>4d}: "
          f"{main['avg_within_50km']:.1f} within 50km ({main['pct_within_50km']:.1f}%), "
          f"mean dist = {main['avg_mean_distance_km']:.0f}km")

# ============================================================
# SAVE RESULTS
# ============================================================
output = {
    "meta": {
        "date": "2026-03-22",
        "description": "Multi-region preservation simulation (Directive 3). Tests whether adding estimated missing settlements in Mesopotamia, Indus Valley, and Peru affects the monument-settlement divergence.",
        "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
        "threshold_km": THRESHOLD_KM,
        "n_mc_trials": N_MC_TRIALS,
        "n_monte_carlo_iterations": N_MONTE_CARLO_ITER,
        "regions": {rkey: {"n_estimated": r["n_estimated"], "center": r["center"], "std": r["std"]}
                    for rkey, r in REGIONS.items()},
        "total_estimated_settlements": sum(r["n_estimated"] for r in REGIONS.values()),
    },
    "global_baseline": {
        "label": "Global Pleiades — no augmentation",
        "n_monumental": len(mon_list),
        "n_settlement": len(set_list),
        "monumental": baseline_mono,
        "settlement": baseline_settle,
        "divergence_D": baseline_D,
    },
    "per_region": {rkey: {
        "label": region_results[rkey]["label"],
        "rationale": region_results[rkey]["rationale"],
        "sweep": region_results[rkey]["sweep"],
    } for rkey in region_order},
    "egypt_only_augmented": {
        "label": "Global Pleiades + 350 Egypt estimated settlements",
        "n_estimated": 350,
        "n_iterations": len(egypt_only_D_values),
        "mean_D": round(eg_mean_D, 2),
        "std_D": round(eg_std_D, 2),
        "pct_D_gt_2": round(eg_pct_gt2, 1),
        "mean_settlement_z": round(float(np.mean(egypt_only_settle_z)), 2),
        "mean_monument_z": round(float(np.mean(egypt_only_mono_z)), 2),
    },
    "multi_region_augmented": {
        "label": "Global Pleiades + 1250 estimated settlements (all 4 regions)",
        "n_estimated_total": 1250,
        "n_iterations": len(combined_D_values),
        "mean_D": round(mean_D, 2),
        "std_D": round(std_D, 2),
        "median_D": round(median_D, 2),
        "pct_D_positive": round(pct_positive, 1),
        "pct_D_gt_2": round(pct_gt2, 1),
        "mean_settlement_z": round(float(np.mean(combined_settle_z)), 2),
        "mean_monument_z": round(float(np.mean(combined_mono_z)), 2),
        "mean_estimated_near_circle": round(float(np.mean(combined_n_near)), 1),
        "D_percentiles": {
            "5th": round(float(np.percentile(combined_D_values, 5)), 2),
            "25th": round(float(np.percentile(combined_D_values, 25)), 2),
            "50th": round(float(np.percentile(combined_D_values, 50)), 2),
            "75th": round(float(np.percentile(combined_D_values, 75)), 2),
            "95th": round(float(np.percentile(combined_D_values, 95)), 2),
        },
    },
    "comparison": {
        "baseline_D": baseline_D,
        "egypt_only_mean_D": round(eg_mean_D, 2),
        "multi_region_mean_D": round(mean_D, 2),
        "change_egypt_only": round(eg_mean_D - baseline_D, 2) if isinstance(baseline_D, (int, float)) else "N/A",
        "change_multi_region": round(mean_D - baseline_D, 2) if isinstance(baseline_D, (int, float)) else "N/A",
        "egypt_only_pct_gt2": round(eg_pct_gt2, 1),
        "multi_region_pct_gt2": round(pct_gt2, 1),
        "interpretation": (
            f"Adding {1250} estimated missing settlements across 4 riverine civilizations "
            f"{'INCREASES' if mean_D > baseline_D else 'DECREASES'} divergence from D={baseline_D} to D={mean_D:.2f}. "
            f"D > 2 in {pct_gt2:.0f}% of iterations. "
            f"The preservation argument {'STRENGTHENS' if mean_D >= baseline_D else 'HOLDS' if pct_gt2 > 99 else 'WEAKENS BUT SURVIVES' if pct_gt2 > 50 else 'DOES NOT SURVIVE'} "
            f"the monument-settlement divergence."
        ),
    },
    "verdict": {
        "egypt_only_survives": eg_pct_gt2 > 50,
        "multi_region_survives": pct_gt2 > 50,
        "multi_region_pct_gt2": round(pct_gt2, 1),
    },
    "elapsed_seconds": round(elapsed, 1),
}

out_path = os.path.join(OUT_DIR, "multi_region_preservation.json")
with open(out_path, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved: {out_path}")
print("\nDONE.")
