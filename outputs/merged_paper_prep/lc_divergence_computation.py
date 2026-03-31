#!/usr/bin/env python3
"""
Land-Constrained Divergence Computation for Portal, p3k14c, DARE
=================================================================
Computes monument and settlement Z-scores separately using the land-constrained
MC baseline, then reports D = Z_monument - Z_settlement for each database.

Saves to: outputs/merged_paper_prep/lc_divergence_all_databases.json
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
POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "merged_paper_prep")
LAND_MASK_PATH = os.path.join(BASE_DIR, "data", "land_mask", "ne_10m_land.shp")
PORTAL_PATH = os.path.join(BASE_DIR, "archive", "great-circle-analysis", "data", "processed", "mp_deduplicated.csv")
P3K14C_PATH = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")
DARE_PATH   = os.path.join(BASE_DIR, "data", "dare", "places2.geojson")

SEED = 42
SIGMA_DEG = 2.0
MAX_ATTEMPTS = 50
THRESHOLD_KM = 50
N_TRIALS_DIV = 1000  # trials for monument/settlement divergence

random.seed(SEED)
np.random.seed(SEED)

# ============================================================
# BUILD RASTERIZED LAND MASK
# ============================================================
print("Building rasterized land mask...")
t0 = time.time()
_land_gdf = gpd.read_file(LAND_MASK_PATH)
_land_union = _land_gdf.union_all()
_land_prepared = prep(_land_union)

RASTER_RES = 0.1
n_lon = int(360 / RASTER_RES)
n_lat = int(180 / RASTER_RES)
_land_raster = np.zeros((n_lat, n_lon), dtype=bool)
for i_lat in range(n_lat):
    lat = -90 + (i_lat + 0.5) * RASTER_RES
    points = [Point(-180 + (j + 0.5) * RASTER_RES, lat) for j in range(n_lon)]
    for j, pt in enumerate(points):
        _land_raster[i_lat, j] = _land_prepared.contains(pt)
print(f"  Rasterized in {time.time()-t0:.1f}s")


def is_on_land(lat, lon):
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
# LC MC ENGINE
# ============================================================
def run_lc_mc(site_lats, site_lons, n_trials):
    n = len(site_lats)
    counts = []
    t_start = time.time()
    for trial in range(n_trials):
        if trial > 0 and trial % 100 == 0:
            elapsed = time.time() - t_start
            rate = trial / elapsed
            eta = (n_trials - trial) / rate
            print(f"    Trial {trial}/{n_trials} ({rate:.1f}/s, ETA {eta:.0f}s)")
        syn_lats = np.empty(n)
        syn_lons = np.empty(n)
        for i in range(n):
            base_lat = site_lats[i]
            base_lon = site_lons[i]
            placed = False
            for _ in range(MAX_ATTEMPTS):
                new_lat = base_lat + random.gauss(0, SIGMA_DEG)
                new_lon = base_lon + random.gauss(0, SIGMA_DEG)
                new_lat = max(-90, min(90, new_lat))
                new_lon = ((new_lon + 180) % 360) - 180
                if is_on_land(new_lat, new_lon):
                    syn_lats[i] = new_lat
                    syn_lons[i] = new_lon
                    placed = True
                    break
            if not placed:
                syn_lats[i] = base_lat
                syn_lons[i] = base_lon
        dists = gc_distances_vectorized(syn_lats, syn_lons)
        counts.append(int(np.sum(dists <= THRESHOLD_KM)))
    return counts


def z_score(observed, mc_counts):
    mc_mean = float(np.mean(mc_counts))
    mc_std  = float(np.std(mc_counts, ddof=1))
    if mc_std == 0:
        return 0.0, mc_mean, mc_std
    return (observed - mc_mean) / mc_std, mc_mean, mc_std


def compute_divergence(mon_lats, mon_lons, set_lats, set_lons, n_trials, label):
    print(f"\n  {label} MONUMENTS ({len(mon_lats)} sites)...")
    obs_mon = int(np.sum(gc_distances_vectorized(np.array(mon_lats), np.array(mon_lons)) <= THRESHOLD_KM))
    mc_mon = run_lc_mc(np.array(mon_lats), np.array(mon_lons), n_trials)
    z_mon, mean_mon, std_mon = z_score(obs_mon, mc_mon)
    enr_mon = obs_mon / mean_mon if mean_mon > 0 else 0.0
    print(f"    obs={obs_mon}, mc_mean={mean_mon:.1f}, Z={z_mon:.2f}, enrichment={enr_mon:.3f}x")

    print(f"  {label} SETTLEMENTS ({len(set_lats)} sites)...")
    obs_set = int(np.sum(gc_distances_vectorized(np.array(set_lats), np.array(set_lons)) <= THRESHOLD_KM))
    mc_set = run_lc_mc(np.array(set_lats), np.array(set_lons), n_trials)
    z_set, mean_set, std_set = z_score(obs_set, mc_set)
    enr_set = obs_set / mean_set if mean_set > 0 else 0.0
    print(f"    obs={obs_set}, mc_mean={mean_set:.1f}, Z={z_set:.2f}, enrichment={enr_set:.3f}x")

    D = z_mon - z_set
    print(f"  {label} D = {D:.2f}")
    return {
        "monument": {
            "n_sites": len(mon_lats), "observed": obs_mon,
            "mc_mean": mean_mon, "mc_std": std_mon,
            "z_score": round(z_mon, 3), "enrichment": round(enr_mon, 3),
        },
        "settlement": {
            "n_sites": len(set_lats), "observed": obs_set,
            "mc_mean": mean_set, "mc_std": std_set,
            "z_score": round(z_set, 3), "enrichment": round(enr_set, 3),
        },
        "divergence_D": round(D, 3),
        "n_trials": n_trials,
        "correction": "land_constrained",
    }


# ============================================================
# IS-EUROPE FILTER (for Portal non-EU subset)
# ============================================================
def is_europe(lat, lon):
    """Approximate exclusion of European sites."""
    # UK and Ireland
    if 49 <= lat <= 61 and -11 <= lon <= 2:
        return True
    # Continental Europe
    if 35 <= lat <= 72 and -10 <= lon <= 40:
        return True
    # Nordic extension
    if lat >= 55 and -5 <= lon <= 30:
        return True
    return False


# ============================================================
# DATABASE CLASSIFICATIONS
# ============================================================

# --- Megalithic Portal ---
PORTAL_SETTLEMENT_TYPES = {"Ancient Village or Settlement"}
PORTAL_MONUMENT_TYPES = {
    "Burial Chamber or Dolmen", "Standing Stone (Menhir)", "Stone Circle",
    "Passage Grave", "Chambered Tomb", "Long Barrow", "Barrow Cemetery",
    "Round Barrow(s)", "Cairn", "Hillfort", "Ancient Temple", "Ancient Palace",
    "Broch or Nuraghe", "Holy Well or Sacred Spring", "Rock Art",
    "Megalithic Tomb or Monument", "Ancient Trackway",
}

# --- p3k14c ---
P3K_MON_KWORDS = ["pyramid", "temple", "tomb", "monument", "shrine",
                  "fortress", "palace", "monastic", "dolmen", "mound",
                  "burial", "tumulus", "megalith"]
P3K_SET_KWORDS = ["settlement", "village", "farm", "dwelling", "camp",
                  "house", "habitation", "residential", "domestic", "site"]

# --- DARE ---
# From type_mapping.json: monumental = temple(22), church(24), tumulus(32)
# Settlement = city(11), town(12), civitas(13), villa(14), camp(15),
#              station(16), settlement(34)
DARE_MON_TYPES = {"22", "24", "32", "35"}   # temple, church, tumulus, + 35
DARE_SET_TYPES = {"11", "12", "13", "14", "15", "16", "34"}


# ============================================================
# LOAD DATASETS
# ============================================================
results = {}

# ---- PORTAL (non-EU) ----
print("\n=== PORTAL (non-EU subset) ===")
mon_lats, mon_lons = [], []
set_lats, set_lons = [], []
with open(PORTAL_PATH) as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            lat = float(row["lat"])
            lon = float(row["lon"])
        except (ValueError, KeyError):
            continue
        if is_europe(lat, lon):
            continue
        t = row.get("type", "").strip()
        if t in PORTAL_SETTLEMENT_TYPES:
            set_lats.append(lat); set_lons.append(lon)
        elif t in PORTAL_MONUMENT_TYPES or t:  # any non-settlement type counts as monument
            if t not in PORTAL_SETTLEMENT_TYPES:
                mon_lats.append(lat); mon_lons.append(lon)

print(f"  Portal non-EU: {len(mon_lats)} monumental, {len(set_lats)} settlement")
results["portal_non_eu"] = compute_divergence(mon_lats, mon_lons, set_lats, set_lons, N_TRIALS_DIV, "Portal non-EU")

# ---- P3K14C ----
print("\n=== P3K14C ===")
mon_lats, mon_lons = [], []
set_lats, set_lons = [], []
with open(P3K14C_PATH) as f:
    reader = csv.DictReader(f)
    seen = set()
    for row in reader:
        site_id = row.get("SiteID", "")
        if site_id in seen:
            continue
        seen.add(site_id)
        try:
            lat = float(row["Lat"])   # Lat column = latitude
            lon = float(row["Long"])  # Long column = longitude
        except (ValueError, KeyError):
            continue
        name = (row.get("SiteName", "") + " " + row.get("SiteType", "")).lower()
        is_mon = any(k in name for k in P3K_MON_KWORDS)
        is_set = any(k in name for k in P3K_SET_KWORDS)
        if is_mon and not is_set:
            mon_lats.append(lat); mon_lons.append(lon)
        elif is_set and not is_mon:
            set_lats.append(lat); set_lons.append(lon)

print(f"  p3k14c: {len(mon_lats)} monumental, {len(set_lats)} settlement")
if len(mon_lats) > 10 and len(set_lats) > 10:
    results["p3k14c"] = compute_divergence(mon_lats, mon_lons, set_lats, set_lons, N_TRIALS_DIV, "p3k14c")
else:
    print("  WARNING: too few sites for reliable divergence")
    results["p3k14c"] = {"error": "insufficient classified sites"}

# ---- DARE ----
print("\n=== DARE ===")
mon_lats, mon_lons = [], []
set_lats, set_lons = [], []
with open(DARE_PATH) as f:
    gj = json.load(f)
for feat in gj["features"]:
    props = feat["properties"]
    t = str(props.get("type", ""))
    geom = feat.get("geometry", {})
    if not geom or geom.get("type") != "Point":
        continue
    coords = geom.get("coordinates", [])
    if len(coords) < 2:
        continue
    lon, lat = float(coords[0]), float(coords[1])
    if t in DARE_MON_TYPES:
        mon_lats.append(lat); mon_lons.append(lon)
    elif t in DARE_SET_TYPES:
        set_lats.append(lat); set_lons.append(lon)

print(f"  DARE: {len(mon_lats)} monumental, {len(set_lats)} settlement")
if len(mon_lats) > 5 and len(set_lats) > 10:
    results["dare"] = compute_divergence(mon_lats, mon_lons, set_lats, set_lons, N_TRIALS_DIV, "DARE")
else:
    print("  WARNING: too few sites")
    results["dare"] = {"error": "insufficient classified sites"}

# ============================================================
# SAVE
# ============================================================
output = {
    "description": "Land-constrained divergence (D = Z_monument - Z_settlement) for Portal, p3k14c, DARE",
    "method": "land_constrained_mc",
    "sigma_deg": SIGMA_DEG,
    "threshold_km": THRESHOLD_KM,
    "max_attempts": MAX_ATTEMPTS,
    "n_trials": N_TRIALS_DIV,
    "results": results
}
out_path = os.path.join(OUT_DIR, "lc_divergence_all_databases.json")
with open(out_path, "w") as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {out_path}")
print("\nSUMMARY:")
for db, res in results.items():
    if "divergence_D" in res:
        print(f"  {db}: D(LC) = {res['divergence_D']:.2f}  "
              f"(mon Z={res['monument']['z_score']:.2f}, "
              f"set Z={res['settlement']['z_score']:.2f})")
    else:
        print(f"  {db}: {res}")
