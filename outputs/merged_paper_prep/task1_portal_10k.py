#!/usr/bin/env python3
"""
Task 1 (Portal N=10,000): Land-Constrained MC for Megalithic Portal
=====================================================================
Definitive run at N=10,000 trials, matching the Pleiades coastal correction run.
Uses rasterized land mask for fast lookups.
"""

import csv
import json
import math
import os
import random
import time

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

SEED = 42
SIGMA_DEG = 2.0
MAX_ATTEMPTS = 50
THRESHOLD_KM = 50
N_TRIALS = 10000

# ============================================================
# RASTERIZED LAND MASK (fast)
# ============================================================
print("Building rasterized land mask...")
t0 = time.time()
_land_gdf = gpd.read_file(LAND_MASK_PATH)
_land_union = _land_gdf.union_all()
_land_prepared = prep(_land_union)

RASTER_RES = 0.1
n_lon = int(360 / RASTER_RES)
n_lat = int(180 / RASTER_RES)
print(f"  Rasterizing at {RASTER_RES}° ({n_lon}x{n_lat})...")

_land_raster = np.zeros((n_lat, n_lon), dtype=bool)
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
    i_lat = int((lat + 90) / RASTER_RES)
    i_lon = int((lon + 180) / RASTER_RES)
    i_lat = max(0, min(n_lat - 1, i_lat))
    i_lon = max(0, min(n_lon - 1, i_lon))
    return bool(_land_raster[i_lat, i_lon])


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
# MC ENGINE (land-constrained only)
# ============================================================
def run_mc_land_constrained(site_lats, site_lons, n_trials):
    n = len(site_lats)
    counts = []
    total_fallbacks = 0
    total_jitters = 0
    t_start = time.time()

    for trial in range(n_trials):
        if trial > 0 and trial % 200 == 0:
            elapsed = time.time() - t_start
            rate = trial / elapsed
            eta = (n_trials - trial) / rate
            print(f"    Trial {trial}/{n_trials} ({rate:.1f}/s, ETA {eta:.0f}s)", flush=True)

        syn_lats = np.empty(n)
        syn_lons = np.empty(n)

        for i in range(n):
            src = random.randrange(n)
            base_lat = site_lats[src]
            base_lon = site_lons[src]
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
        "p_value": round(p, 8),
    }


# ============================================================
# MAIN
# ============================================================
print(f"\nLoading Portal...")
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

site_lats = np.array(lats)
site_lons = np.array(lons)
dists = gc_distances_vectorized(site_lats, site_lons)
observed = int(np.sum(dists <= THRESHOLD_KM))
print(f"  Observed within {THRESHOLD_KM} km: {observed}")

print(f"\nLand-constrained MC (N={N_TRIALS})...")
random.seed(SEED + 1)
np.random.seed(SEED + 1)
mc_counts, fb_rate, elapsed = run_mc_land_constrained(lats, lons, N_TRIALS)
stats = compute_stats(observed, mc_counts)
print(f"\n  Z={stats['z_score']}, enrichment={stats['enrichment']}x, p={stats['p_value']}")
print(f"  MC mean={stats['mc_mean']}, std={stats['mc_std']}")
print(f"  Fallback rate: {fb_rate:.4%}")
print(f"  Time: {elapsed:.0f}s")

result = {
    "date": "2026-03-30",
    "database": "Megalithic Portal",
    "n_sites": len(lats),
    "threshold_km": THRESHOLD_KM,
    "observed_within_threshold": observed,
    "n_trials": N_TRIALS,
    "sigma_deg": SIGMA_DEG,
    "max_jitter_attempts": MAX_ATTEMPTS,
    "method": "Distribution-matched MC, Gaussian jitter sigma=2deg, rasterized land mask (0.1deg NE 1:10M), land-constrained only",
    "land_constrained": stats,
    "fallback_rate": round(fb_rate, 6),
    "time_s": round(elapsed, 1),
    "note": "Definitive N=10,000 run matching Pleiades coastal correction methodology"
}

out_path = os.path.join(OUT_DIR, "portal_10k_land_constrained.json")
with open(out_path, 'w') as f:
    json.dump(result, f, indent=2)
print(f"\nResults saved to {out_path}")
print("\nDone!")
