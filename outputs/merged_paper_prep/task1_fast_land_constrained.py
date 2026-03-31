#!/usr/bin/env python3
"""
Task 1 (optimized): Land-Constrained MC for Portal, p3k14c, XRONOS
===================================================================
Uses a rasterized land mask (~0.1° resolution) for fast lookups instead
of shapely point-in-polygon. Falls back to shapely for edge cases.

Runs N=200 trials at 50 km for all three databases, then N=1000 for Portal only.
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

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "merged_paper_prep")
LAND_MASK_PATH = os.path.join(BASE_DIR, "data", "land_mask", "ne_10m_land.shp")
PORTAL_PATH = os.path.join(BASE_DIR, "archive", "great-circle-analysis", "data", "processed", "mp_deduplicated.csv")
P3K14C_PATH = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")
XRONOS_PATH = os.path.join(BASE_DIR, "data", "xronos", "xronos_sites.csv")

SEED = 42
SIGMA_DEG = 2.0
MAX_ATTEMPTS = 50
THRESHOLD_KM = 50

# ============================================================
# RASTERIZED LAND MASK (fast)
# ============================================================
print("Building rasterized land mask...")
t0 = time.time()
_land_gdf = gpd.read_file(LAND_MASK_PATH)
_land_union = _land_gdf.union_all()
_land_prepared = prep(_land_union)

# Build raster at 0.1° resolution (3600 x 1800 grid)
RASTER_RES = 0.1
n_lon = int(360 / RASTER_RES)
n_lat = int(180 / RASTER_RES)
print(f"  Rasterizing at {RASTER_RES}° ({n_lon}x{n_lat})...")

# Sample grid centers
_land_raster = np.zeros((n_lat, n_lon), dtype=bool)
# Do this in chunks for speed
for i_lat in range(n_lat):
    lat = -90 + (i_lat + 0.5) * RASTER_RES
    points = [Point(-180 + (j + 0.5) * RASTER_RES, lat) for j in range(n_lon)]
    for j, pt in enumerate(points):
        _land_raster[i_lat, j] = _land_prepared.contains(pt)
    if i_lat % 200 == 0:
        print(f"    Row {i_lat}/{n_lat}")

land_frac = _land_raster.sum() / _land_raster.size
print(f"  Rasterized in {time.time()-t0:.1f}s, land fraction: {land_frac:.3f}")


def is_on_land_fast(lat, lon):
    """Fast land check using raster. Returns True if on land."""
    i_lat = int((lat + 90) / RASTER_RES)
    i_lon = int((lon + 180) / RASTER_RES)
    i_lat = max(0, min(n_lat - 1, i_lat))
    i_lon = max(0, min(n_lon - 1, i_lon))
    return bool(_land_raster[i_lat, i_lon])


# Verify raster
for label, lat, lon, expected in [
    ("Giza", 29.98, 31.13, True),
    ("Pacific", 0.0, -150.0, False),
    ("Easter Is", -27.12, -109.35, True),
    ("London", 51.51, -0.12, True),
]:
    result = is_on_land_fast(lat, lon)
    status = "OK" if result == expected else "WARN"
    print(f"  Verify {label}: {result} (expected {expected}) [{status}]")

# Benchmark
t0 = time.time()
_test = [(random.uniform(-90, 90), random.uniform(-180, 180)) for _ in range(100000)]
for lat, lon in _test:
    is_on_land_fast(lat, lon)
rate = 100000 / (time.time() - t0)
print(f"  Raster check speed: {rate:.0f} points/sec")


# ============================================================
# GEODESIC MATH
# ============================================================
def gc_distances_vectorized(lats, lons):
    lats_r = np.radians(lats)
    lons_r = np.radians(lons)
    pole_lat_r = math.radians(POLE_LAT)
    pole_lon_r = math.radians(POLE_LON)
    px = math.cos(pole_lat_r) * math.cos(pole_lon_r)
    py = math.cos(pole_lat_r) * math.sin(pole_lon_r)
    pz = math.sin(pole_lat_r)
    sx = np.cos(lats_r) * np.cos(lons_r)
    sy = np.cos(lats_r) * np.sin(lons_r)
    sz = np.sin(lats_r)
    dots = np.clip(sx*px + sy*py + sz*pz, -1.0, 1.0)
    angular_dist = np.arccos(dots)
    return np.abs(angular_dist - math.pi/2) * EARTH_R_KM


# ============================================================
# MC ENGINE
# ============================================================
def run_mc(site_lats, site_lons, n_trials, land_constrained=True):
    n = len(site_lats)
    counts = []
    total_fallbacks = 0
    total_jitters = 0
    t_start = time.time()

    for trial in range(n_trials):
        if trial > 0 and trial % 50 == 0:
            elapsed = time.time() - t_start
            rate = trial / elapsed
            eta = (n_trials - trial) / rate
            print(f"    Trial {trial}/{n_trials} ({rate:.1f}/s, ETA {eta:.0f}s)")

        syn_lats = np.empty(n)
        syn_lons = np.empty(n)

        for i in range(n):
            src = random.randrange(n)
            base_lat = site_lats[src]
            base_lon = site_lons[src]

            if land_constrained:
                placed = False
                for _ in range(MAX_ATTEMPTS):
                    new_lat = base_lat + random.gauss(0, SIGMA_DEG)
                    new_lon = base_lon + random.gauss(0, SIGMA_DEG)
                    new_lat = max(-90, min(90, new_lat))
                    new_lon = ((new_lon + 180) % 360) - 180
                    if is_on_land_fast(new_lat, new_lon):
                        syn_lats[i] = new_lat
                        syn_lons[i] = new_lon
                        placed = True
                        break
                if not placed:
                    syn_lats[i] = base_lat
                    syn_lons[i] = base_lon
                    total_fallbacks += 1
                total_jitters += 1
            else:
                new_lat = base_lat + random.gauss(0, SIGMA_DEG)
                new_lon = base_lon + random.gauss(0, SIGMA_DEG)
                syn_lats[i] = max(-90, min(90, new_lat))
                syn_lons[i] = ((new_lon + 180) % 360) - 180

        dists = gc_distances_vectorized(syn_lats, syn_lons)
        c = int(np.sum(dists <= THRESHOLD_KM))
        counts.append(c)

    fb_rate = total_fallbacks / total_jitters if total_jitters > 0 else 0
    elapsed = time.time() - t_start
    return counts, fb_rate, elapsed


def compute_stats(observed, mc_counts):
    mc_mean = float(np.mean(mc_counts))
    mc_std = float(np.std(mc_counts, ddof=0))
    z = (observed - mc_mean) / mc_std if mc_std > 0 else 0.0
    enrichment = observed / mc_mean if mc_mean > 0 else float('inf')
    p = float(np.mean([c >= observed for c in mc_counts]))
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
    print(f"Loading Portal...")
    lats, lons = [], []
    with open(PORTAL_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat, lon = float(row['lat']), float(row['lon'])
                if -90 <= lat <= 90 and -180 <= lon <= 180:
                    lats.append(lat)
                    lons.append(lon)
            except (ValueError, KeyError):
                continue
    print(f"  {len(lats)} sites")
    return lats, lons


def load_p3k14c():
    print(f"Loading p3k14c...")
    seen = {}
    with open(P3K14C_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat, lon = float(row['Lat']), float(row['Long'])
                sid = row.get('SiteID', '') or f"{lat:.4f}_{lon:.4f}"
                if sid not in seen and -90 <= lat <= 90 and -180 <= lon <= 180:
                    seen[sid] = (lat, lon)
            except (ValueError, KeyError):
                continue
    lats = [v[0] for v in seen.values()]
    lons = [v[1] for v in seen.values()]
    print(f"  {len(lats)} unique sites")
    return lats, lons


def load_xronos():
    print(f"Loading XRONOS...")
    lats, lons = [], []
    with open(XRONOS_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat, lon = float(row['lat']), float(row['lon'])
                if -90 <= lat <= 90 and -180 <= lon <= 180:
                    lats.append(lat)
                    lons.append(lon)
            except (ValueError, KeyError):
                continue
    print(f"  {len(lats)} sites")
    return lats, lons


# ============================================================
# MAIN
# ============================================================
def analyze_one(name, lats, lons, n_trials):
    n = len(lats)
    dists = gc_distances_vectorized(np.array(lats), np.array(lons))
    observed = int(np.sum(dists <= THRESHOLD_KM))
    print(f"\n{'='*60}")
    print(f"  {name}: {n} sites, {observed} within {THRESHOLD_KM} km")
    print(f"{'='*60}")

    # Unconstrained
    print(f"  Unconstrained MC (N={n_trials})...")
    random.seed(SEED)
    np.random.seed(SEED)
    mc_unc, _, t_unc = run_mc(lats, lons, n_trials, land_constrained=False)
    stats_unc = compute_stats(observed, mc_unc)
    print(f"  → Z={stats_unc['z_score']}, enrich={stats_unc['enrichment']}x ({t_unc:.0f}s)")

    # Land-constrained
    print(f"  Land-constrained MC (N={n_trials})...")
    random.seed(SEED + 1)
    np.random.seed(SEED + 1)
    mc_lc, fb_rate, t_lc = run_mc(lats, lons, n_trials, land_constrained=True)
    stats_lc = compute_stats(observed, mc_lc)
    print(f"  → Z={stats_lc['z_score']}, enrich={stats_lc['enrichment']}x, fb={fb_rate:.4%} ({t_lc:.0f}s)")

    z_reduction = None
    if stats_unc['z_score'] != 0:
        z_reduction = round((1 - stats_lc['z_score'] / stats_unc['z_score']) * 100, 1)

    return {
        "database": name,
        "n_sites": n,
        "threshold_km": THRESHOLD_KM,
        "observed_within_threshold": observed,
        "n_trials": n_trials,
        "unconstrained": stats_unc,
        "land_constrained": stats_lc,
        "fallback_rate": round(fb_rate, 6),
        "z_reduction_pct": z_reduction,
        "time_unconstrained_s": round(t_unc, 1),
        "time_land_constrained_s": round(t_lc, 1),
    }


if __name__ == "__main__":
    results = {
        "date": "2026-03-30",
        "method": "Distribution-matched MC, Gaussian jitter σ=2°, rasterized land mask (0.1° NE 1:10M)",
        "threshold_km": THRESHOLD_KM,
        "sigma_deg": SIGMA_DEG,
        "max_jitter_attempts": MAX_ATTEMPTS,
    }

    # Load data
    portal_lats, portal_lons = load_portal()
    p3k14c_lats, p3k14c_lons = load_p3k14c()
    xronos_lats, xronos_lons = load_xronos()

    # Run analyses
    # Start with smaller databases (faster), then Portal
    results["xronos"] = analyze_one("XRONOS", xronos_lats, xronos_lons, n_trials=1000)
    results["p3k14c"] = analyze_one("p3k14c", p3k14c_lats, p3k14c_lons, n_trials=1000)
    results["megalithic_portal"] = analyze_one("Megalithic Portal", portal_lats, portal_lons, n_trials=200)

    # Save
    out_path = os.path.join(OUT_DIR, "land_constrained_all_databases.json")
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n\nResults saved to {out_path}")
    print("\nDone!")
