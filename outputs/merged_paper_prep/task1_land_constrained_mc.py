#!/usr/bin/env python3
"""
Task 1: Land-Constrained MC for Portal, p3k14c, XRONOS
========================================================
Adapts the Pleiades coastal bias correction (outputs/coastal_bias_correction/)
to three additional databases.

Method: Distribution-matched MC with Gaussian jitter (σ=2°), ocean rejection
using Natural Earth 1:10M land mask. Fallback to original coordinates if
max_attempts exceeded (conservative — doesn't inflate Z).

Runs at N=1,000 trials first for quick results, then N=10,000 for definitive.
"""

import csv
import json
import math
import os
import random
import time
import sys

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from shapely.prepared import prep

# ============================================================
# CONSTANTS
# ============================================================
POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "merged_paper_prep")
os.makedirs(OUT_DIR, exist_ok=True)

LAND_MASK_PATH = os.path.join(BASE_DIR, "data", "land_mask", "ne_10m_land.shp")
PORTAL_PATH = os.path.join(BASE_DIR, "archive", "great-circle-analysis", "data", "processed", "mp_deduplicated.csv")
P3K14C_PATH = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")
XRONOS_PATH = os.path.join(BASE_DIR, "data", "xronos", "xronos_sites.csv")

SEED = 42
SIGMA_DEG = 2.0
MAX_ATTEMPTS = 100
THRESHOLDS_KM = [25, 50]

# ============================================================
# LAND MASK
# ============================================================
print("Loading Natural Earth land mask...")
t0 = time.time()
_land_gdf = gpd.read_file(LAND_MASK_PATH)
_land_union = _land_gdf.union_all()
_land_prepared = prep(_land_union)
print(f"  Land mask loaded in {time.time()-t0:.1f}s")

# Verify
for label, lat, lon, expected in [
    ("Giza", 29.98, 31.13, True),
    ("Pacific Ocean", 0.0, -150.0, False),
    ("Easter Island", -27.12, -109.35, True),
]:
    result = _land_prepared.contains(Point(lon, lat))
    status = "OK" if result == expected else "MISMATCH"
    print(f"  Verify {label}: {result} [{status}]")


def is_on_land(lat, lon):
    return _land_prepared.contains(Point(lon, lat))


# ============================================================
# GEODESIC MATH
# ============================================================
def gc_distances_vectorized(lats, lons, pole_lat, pole_lon):
    lats_r = np.radians(lats)
    lons_r = np.radians(lons)
    pole_lat_r = math.radians(pole_lat)
    pole_lon_r = math.radians(pole_lon)
    px = math.cos(pole_lat_r) * math.cos(pole_lon_r)
    py = math.cos(pole_lat_r) * math.sin(pole_lon_r)
    pz = math.sin(pole_lat_r)
    sx = np.cos(lats_r) * np.cos(lons_r)
    sy = np.cos(lats_r) * np.sin(lons_r)
    sz = np.sin(lats_r)
    dots = np.clip(sx*px + sy*py + sz*pz, -1.0, 1.0)
    angular_dist = np.arccos(dots)
    return np.abs(angular_dist - math.pi/2) * EARTH_R_KM


def count_within_corridor(lats, lons, threshold_km=50):
    dists = gc_distances_vectorized(np.array(lats), np.array(lons), POLE_LAT, POLE_LON)
    return int(np.sum(dists <= threshold_km))


# ============================================================
# JITTER
# ============================================================
def jitter_to_land(lat, lon, sigma_deg=SIGMA_DEG, max_attempts=MAX_ATTEMPTS):
    for _ in range(max_attempts):
        new_lat = lat + random.gauss(0, sigma_deg)
        new_lon = lon + random.gauss(0, sigma_deg)
        new_lat = max(-90, min(90, new_lat))
        new_lon = ((new_lon + 180) % 360) - 180
        if is_on_land(new_lat, new_lon):
            return new_lat, new_lon, False
    return lat, lon, True


def jitter_unconstrained(lat, lon, sigma_deg=SIGMA_DEG):
    new_lat = lat + random.gauss(0, sigma_deg)
    new_lon = lon + random.gauss(0, sigma_deg)
    new_lat = max(-90, min(90, new_lat))
    new_lon = ((new_lon + 180) % 360) - 180
    return new_lat, new_lon


# ============================================================
# MC BASELINES
# ============================================================
def run_mc(site_lats, site_lons, n_trials, threshold_km, land_constrained=True):
    """Run distribution-matched MC baseline."""
    n = len(site_lats)
    counts = []
    total_fallbacks = 0
    total_jitters = 0

    for trial in range(n_trials):
        if trial % 200 == 0 and trial > 0:
            elapsed = time.time() - _trial_start
            rate = trial / elapsed
            eta = (n_trials - trial) / rate
            print(f"    Trial {trial}/{n_trials} ({rate:.1f} trials/s, ETA {eta:.0f}s)")

        syn_lats = []
        syn_lons = []
        for i in range(n):
            src_idx = random.randrange(n)
            if land_constrained:
                new_lat, new_lon, fell_back = jitter_to_land(
                    site_lats[src_idx], site_lons[src_idx]
                )
                if fell_back:
                    total_fallbacks += 1
                total_jitters += 1
            else:
                new_lat, new_lon = jitter_unconstrained(
                    site_lats[src_idx], site_lons[src_idx]
                )
            syn_lats.append(new_lat)
            syn_lons.append(new_lon)

        c = count_within_corridor(syn_lats, syn_lons, threshold_km)
        counts.append(c)

    fallback_rate = total_fallbacks / total_jitters if total_jitters > 0 else 0
    return counts, fallback_rate


def compute_stats(observed, mc_counts):
    mc_mean = np.mean(mc_counts)
    mc_std = np.std(mc_counts, ddof=0)
    z = (observed - mc_mean) / mc_std if mc_std > 0 else 0.0
    enrichment = observed / mc_mean if mc_mean > 0 else float('inf')
    p = np.mean([c >= observed for c in mc_counts])
    return {
        "observed": observed,
        "mc_mean": round(mc_mean, 2),
        "mc_std": round(mc_std, 2),
        "z_score": round(z, 2),
        "enrichment": round(enrichment, 3),
        "p_value": round(p, 6),
    }


# ============================================================
# DATA LOADING
# ============================================================
def load_portal():
    print(f"Loading Megalithic Portal from {PORTAL_PATH}...")
    lats, lons = [], []
    with open(PORTAL_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row['lat'])
                lon = float(row['lon'])
                if -90 <= lat <= 90 and -180 <= lon <= 180:
                    lats.append(lat)
                    lons.append(lon)
            except (ValueError, KeyError):
                continue
    print(f"  Loaded {len(lats)} Portal sites")
    return lats, lons


def load_p3k14c_sites():
    """Load p3k14c, deduplicate by SiteID."""
    print(f"Loading p3k14c from {P3K14C_PATH}...")
    seen_sites = {}
    with open(P3K14C_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row['Lat'])
                lon = float(row['Long'])
                site_id = row.get('SiteID', '')
                if not site_id:
                    site_id = f"{lat:.4f}_{lon:.4f}"
                if site_id not in seen_sites and -90 <= lat <= 90 and -180 <= lon <= 180:
                    seen_sites[site_id] = (lat, lon)
            except (ValueError, KeyError):
                continue
    lats = [v[0] for v in seen_sites.values()]
    lons = [v[1] for v in seen_sites.values()]
    print(f"  Loaded {len(lats)} unique p3k14c sites (from {P3K14C_PATH})")
    return lats, lons


def load_xronos_sites():
    print(f"Loading XRONOS from {XRONOS_PATH}...")
    lats, lons = [], []
    with open(XRONOS_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row['lat'])
                lon = float(row['lon'])
                if -90 <= lat <= 90 and -180 <= lon <= 180:
                    lats.append(lat)
                    lons.append(lon)
            except (ValueError, KeyError):
                continue
    print(f"  Loaded {len(lats)} XRONOS sites")
    return lats, lons


# ============================================================
# MAIN
# ============================================================
def analyze_database(name, lats, lons, n_trials_list):
    """Run unconstrained and land-constrained MC at multiple trial counts."""
    global _trial_start
    n_sites = len(lats)

    # Compute observed counts at each threshold
    dists = gc_distances_vectorized(np.array(lats), np.array(lons), POLE_LAT, POLE_LON)
    observed = {}
    for thresh in THRESHOLDS_KM:
        observed[thresh] = int(np.sum(dists <= thresh))
    print(f"\n{'='*60}")
    print(f"  {name}: {n_sites} sites")
    for thresh in THRESHOLDS_KM:
        print(f"  Observed within {thresh} km: {observed[thresh]}")
    print(f"{'='*60}")

    results = {
        "database": name,
        "n_sites": n_sites,
        "thresholds_km": THRESHOLDS_KM,
        "observed": observed,
        "sigma_deg": SIGMA_DEG,
        "seed": SEED,
    }

    for n_trials in n_trials_list:
        for thresh in THRESHOLDS_KM:
            # Unconstrained
            print(f"\n  [{name}] Unconstrained MC: N={n_trials}, threshold={thresh} km")
            random.seed(SEED)
            np.random.seed(SEED)
            _trial_start = time.time()
            mc_counts_unc, _ = run_mc(lats, lons, n_trials, thresh, land_constrained=False)
            t_unc = time.time() - _trial_start
            stats_unc = compute_stats(observed[thresh], mc_counts_unc)
            print(f"    Z={stats_unc['z_score']}, enrichment={stats_unc['enrichment']}x ({t_unc:.1f}s)")

            # Land-constrained
            print(f"  [{name}] Land-constrained MC: N={n_trials}, threshold={thresh} km")
            random.seed(SEED + 1)
            np.random.seed(SEED + 1)
            _trial_start = time.time()
            mc_counts_lc, fb_rate = run_mc(lats, lons, n_trials, thresh, land_constrained=True)
            t_lc = time.time() - _trial_start
            stats_lc = compute_stats(observed[thresh], mc_counts_lc)
            print(f"    Z={stats_lc['z_score']}, enrichment={stats_lc['enrichment']}x, fallback={fb_rate:.4%} ({t_lc:.1f}s)")

            key = f"threshold_{thresh}km_N{n_trials}"
            results[key] = {
                "n_trials": n_trials,
                "threshold_km": thresh,
                "unconstrained": stats_unc,
                "land_constrained": stats_lc,
                "fallback_rate": round(fb_rate, 6),
                "z_reduction_pct": round(
                    (1 - stats_lc['z_score'] / stats_unc['z_score']) * 100, 1
                ) if stats_unc['z_score'] != 0 else None,
                "time_unconstrained_s": round(t_unc, 1),
                "time_land_constrained_s": round(t_lc, 1),
            }

    return results


if __name__ == "__main__":
    all_results = {"date": "2026-03-30", "method": "land-constrained distribution-matched MC"}

    # Determine trial counts — start with 1000, add 10000 if user requests
    n_trials_list = [1000]
    if "--full" in sys.argv:
        n_trials_list = [1000, 10000]

    # Load all databases
    portal_lats, portal_lons = load_portal()
    p3k14c_lats, p3k14c_lons = load_p3k14c_sites()
    xronos_lats, xronos_lons = load_xronos_sites()

    # Run analyses
    all_results["megalithic_portal"] = analyze_database(
        "Megalithic Portal", portal_lats, portal_lons, n_trials_list
    )
    all_results["p3k14c"] = analyze_database(
        "p3k14c", p3k14c_lats, p3k14c_lons, n_trials_list
    )
    all_results["xronos"] = analyze_database(
        "XRONOS", xronos_lats, xronos_lons, n_trials_list
    )

    # Save
    out_path = os.path.join(OUT_DIR, "land_constrained_all_databases.json")
    with open(out_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {out_path}")
