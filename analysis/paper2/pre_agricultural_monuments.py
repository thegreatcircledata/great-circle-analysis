#!/usr/bin/env python3
"""
Pre-Agricultural Monument Geometry Test
========================================
Do pre-agricultural monumental sites follow ANY great circle alignment pattern?
Or do they cluster by ecological zone rather than geometry?

Tests:
1. Great Circle Alignment (Alison GC, top-10 circles, full pole scan)
2. Ecological Zone Clustering (elevation bands, regional patterns)
3. Terrain Signature Comparison (TPI/prominence vs civilizational monuments)
4. Taş Tepeler Specific Analysis (cluster geometry, orientation)

Output: outputs/pre_agricultural_monuments/
"""

import json, math, os, sys, time
import numpy as np
import netCDF4 as nc

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
np.random.seed(42)

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "pre_agricultural_monuments")
os.makedirs(OUT_DIR, exist_ok=True)

EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * np.pi / 2
ALISON_POLE_LAT = 59.682122
ALISON_POLE_LON = -138.646087

# ============================================================
# PRE-AGRICULTURAL SITE DATABASE
# ============================================================
PRE_AG_SITES = [
    # Taş Tepeler cluster (SE Turkey, PPNA, 9600-8000 BCE)
    {"name": "Göbekli Tepe", "lat": 37.223, "lon": 38.922, "date_bce": 9500, "region": "Near East"},
    {"name": "Karahan Tepe", "lat": 37.058, "lon": 39.198, "date_bce": 9400, "region": "Near East"},
    {"name": "Harbetsuvan Tepesi", "lat": 37.35, "lon": 38.75, "date_bce": 9000, "region": "Near East"},
    {"name": "Sayburç", "lat": 37.42, "lon": 38.60, "date_bce": 9000, "region": "Near East"},
    {"name": "Boncuklu Tarla", "lat": 37.35, "lon": 41.35, "date_bce": 9200, "region": "Near East"},
    {"name": "Çakmaktepe", "lat": 37.15, "lon": 39.05, "date_bce": 9000, "region": "Near East"},
    {"name": "Sefer Tepe", "lat": 37.20, "lon": 39.10, "date_bce": 9000, "region": "Near East"},
    {"name": "Taslı Tepe", "lat": 37.10, "lon": 38.80, "date_bce": 9000, "region": "Near East"},

    # Additional PPNA / pre-agricultural Near East
    {"name": "Tell Qaramel", "lat": 36.373, "lon": 37.272, "date_bce": 10650, "region": "Near East"},
    {"name": "Jerf el-Ahmar", "lat": 36.383, "lon": 38.181, "date_bce": 9500, "region": "Near East"},
    {"name": "Mureybet", "lat": 36.043, "lon": 38.129, "date_bce": 10200, "region": "Near East"},
    {"name": "Nevali Çori", "lat": 37.518, "lon": 38.603, "date_bce": 8600, "region": "Near East"},
    {"name": "Wadi Hammeh 27", "lat": 32.45, "lon": 35.62, "date_bce": 12500, "region": "Near East"},

    # Levantine PPNA
    {"name": "Jericho Tower", "lat": 31.871, "lon": 35.444, "date_bce": 8300, "region": "Near East"},
    {"name": "Wadi Feynan 16", "lat": 30.62, "lon": 35.42, "date_bce": 9500, "region": "Near East"},

    # North American archaic mound complexes (pre-agricultural in context)
    {"name": "Watson Brake", "lat": 32.528, "lon": -91.909, "date_bce": 3500, "region": "North America"},
    {"name": "Poverty Point", "lat": 32.637, "lon": -91.408, "date_bce": 1700, "region": "North America"},
    {"name": "Horr's Island", "lat": 25.92, "lon": -81.72, "date_bce": 3400, "region": "North America"},
    {"name": "Tomoka Mounds", "lat": 29.33, "lon": -81.10, "date_bce": 3000, "region": "North America"},
    {"name": "Frenchman's Bend", "lat": 32.640, "lon": -92.057, "date_bce": 3500, "region": "North America"},
    {"name": "Hedgepeth Mounds", "lat": 32.679, "lon": -92.653, "date_bce": 3000, "region": "North America"},
    {"name": "L'Anse Amour", "lat": 51.480, "lon": -56.867, "date_bce": 5500, "region": "North America"},

    # South America
    {"name": "Caral", "lat": -10.893, "lon": -77.520, "date_bce": 3000, "region": "South America"},
    {"name": "Nanchoc Mounds", "lat": -6.959, "lon": -79.242, "date_bce": 5700, "region": "South America"},

    # Japanese Jomon (pre-agricultural monumental)
    {"name": "Sannai-Maruyama", "lat": 40.808, "lon": 140.698, "date_bce": 3900, "region": "East Asia"},
    {"name": "Ōyu stone circles", "lat": 39.933, "lon": 140.800, "date_bce": 2000, "region": "East Asia"},
    {"name": "Isedotai", "lat": 40.201, "lon": 140.348, "date_bce": 2000, "region": "East Asia"},
    {"name": "Kiusu Earthworks", "lat": 42.885, "lon": 141.716, "date_bce": 1200, "region": "East Asia"},
    {"name": "Omori Katsuyama", "lat": 40.699, "lon": 140.358, "date_bce": 1500, "region": "East Asia"},
    {"name": "Washinoki", "lat": 42.116, "lon": 140.526, "date_bce": 2000, "region": "East Asia"},

    # European Mesolithic/early Neolithic (pre-local-agriculture)
    {"name": "Lepenski Vir", "lat": 44.558, "lon": 22.025, "date_bce": 6500, "region": "Europe"},
    {"name": "Stonehenge (phase 1)", "lat": 51.179, "lon": -1.826, "date_bce": 3000, "region": "Europe"},
    {"name": "Newgrange", "lat": 53.695, "lon": -6.475, "date_bce": 3200, "region": "Europe"},
    {"name": "Carnac", "lat": 47.584, "lon": -3.078, "date_bce": 4500, "region": "Europe"},
    {"name": "Teviec", "lat": 47.556, "lon": -3.165, "date_bce": 5300, "region": "Europe"},
    {"name": "Hoedic", "lat": 47.340, "lon": -2.878, "date_bce": 6000, "region": "Europe"},
    {"name": "Skateholm", "lat": 55.379, "lon": 13.471, "date_bce": 5250, "region": "Europe"},
    {"name": "Almendres Cromlech", "lat": 38.558, "lon": -8.062, "date_bce": 6000, "region": "Europe"},
    {"name": "Shigir Idol site", "lat": 57.45, "lon": 60.05, "date_bce": 9500, "region": "Europe"},
]

N_SITES = len(PRE_AG_SITES)
site_lats = np.array([s['lat'] for s in PRE_AG_SITES])
site_lons = np.array([s['lon'] for s in PRE_AG_SITES])

print(f"Pre-agricultural sites loaded: {N_SITES}")
for region in sorted(set(s['region'] for s in PRE_AG_SITES)):
    count = sum(1 for s in PRE_AG_SITES if s['region'] == region)
    print(f"  {region}: {count}")

# ============================================================
# TOP-10 HIGH-D CIRCLE POLES (from systematic search)
# ============================================================
TOP10_POLES = [
    {"name": "S California", "lat": 36.0, "lon": -117.5, "D": 24.28},
    {"name": "Pacific NW (Alison family)", "lat": 45.5, "lon": -141.5, "D": 22.69},
    {"name": "Pacific NW 2", "lat": 43.5, "lon": -135.0, "D": 21.82},
    {"name": "Pacific NW 3", "lat": 47.5, "lon": -157.0, "D": 20.99},
    {"name": "Sumatra/SE Asia", "lat": 0.5, "lon": 102.5, "D": 19.59},
    {"name": "SW USA", "lat": 29.5, "lon": -108.0, "D": 17.96},
    {"name": "Central America", "lat": 6.0, "lon": -82.5, "D": 16.73},
    {"name": "Philippines/Taiwan", "lat": 17.5, "lon": 119.5, "D": 15.94},
    {"name": "Japan/Pacific", "lat": 34.0, "lon": 140.5, "D": 15.84},
    {"name": "NW Pacific", "lat": 43.0, "lon": 160.0, "D": 15.80},
]

# ============================================================
# VECTORIZED HELPER FUNCTIONS (from systematic_gc_search.py)
# ============================================================
def dist_from_gc_vec(site_lats, site_lons, pole_lat, pole_lon):
    """Distance from each site to great circle defined by pole (km)."""
    lat1r = np.radians(pole_lat)
    lon1r = np.radians(pole_lon)
    lat2r = np.radians(site_lats)
    lon2r = np.radians(site_lons)
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat/2)**2 + np.cos(lat1r)*np.cos(lat2r)*np.sin(dlon/2)**2
    d = EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QUARTER_CIRC)

def count_within(site_lats, site_lons, pole_lat, pole_lon, threshold_km):
    dists = dist_from_gc_vec(site_lats, site_lons, pole_lat, pole_lon)
    return int(np.sum(dists <= threshold_km))

def random_pole():
    """Uniform random point on sphere."""
    z = np.random.uniform(-1, 1)
    lon = np.random.uniform(-180, 180)
    lat = np.degrees(np.arcsin(z))
    return lat, lon

def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km."""
    r = EARTH_R_KM
    la1, lo1, la2, lo2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = la2 - la1
    dlon = lo2 - lo1
    a = np.sin(dlat/2)**2 + np.cos(la1)*np.cos(la2)*np.sin(dlon/2)**2
    return r * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))

# ============================================================
# LOAD ETOPO1 FOR ELEVATION
# ============================================================
print("\nLoading ETOPO1...")
etopo_path = os.path.join(BASE_DIR, "data", "geophysical", "etopo", "ETOPO1_Ice_g_gmt4.grd")
ds = nc.Dataset(etopo_path)
etopo_lons_arr = ds.variables['x'][:]
etopo_lats_arr = ds.variables['y'][:]
etopo_z = ds.variables['z']
etopo_lon_min, etopo_lon_step = float(etopo_lons_arr[0]), float(etopo_lons_arr[1] - etopo_lons_arr[0])
etopo_lat_min, etopo_lat_step = float(etopo_lats_arr[0]), float(etopo_lats_arr[1] - etopo_lats_arr[0])
n_lats_etopo = len(etopo_lats_arr)
n_lons_etopo = len(etopo_lons_arr)
print(f"ETOPO1 loaded: {n_lats_etopo} x {n_lons_etopo}")

def get_elevation(lat, lon):
    """Get elevation for a single point from ETOPO1."""
    lat_idx = int(round((lat - etopo_lat_min) / etopo_lat_step))
    lon_idx = int(round((lon - etopo_lon_min) / etopo_lon_step))
    lat_idx = max(0, min(lat_idx, n_lats_etopo - 1))
    lon_idx = max(0, min(lon_idx, n_lons_etopo - 1))
    return float(etopo_z[lat_idx, lon_idx])

def get_elevations_batch(lats, lons):
    """Get elevations for arrays of points."""
    lats = np.asarray(lats, dtype=float)
    lons = np.asarray(lons, dtype=float)
    lat_idx = np.round((lats - etopo_lat_min) / etopo_lat_step).astype(int)
    lon_idx = np.round((lons - etopo_lon_min) / etopo_lon_step).astype(int)
    lat_idx = np.clip(lat_idx, 0, n_lats_etopo - 1)
    lon_idx = np.clip(lon_idx, 0, n_lons_etopo - 1)
    elevs = np.empty(len(lats))
    for i in range(len(lats)):
        elevs[i] = float(etopo_z[lat_idx[i], lon_idx[i]])
    return elevs

def compute_tpi(lat, lon, radius_km=10, n_ring=16):
    """Topographic Position Index: site elevation minus mean of ring."""
    site_elev = get_elevation(lat, lon)
    ring_elevs = []
    for i in range(n_ring):
        bearing = 2 * np.pi * i / n_ring
        # Approximate offset in degrees
        dlat = (radius_km / 111.0) * np.cos(bearing)
        dlon = (radius_km / (111.0 * np.cos(np.radians(lat)))) * np.sin(bearing)
        ring_elevs.append(get_elevation(lat + dlat, lon + dlon))
    return site_elev - np.mean(ring_elevs)

def compute_prominence(lat, lon, radius_km=25, n_ring=24):
    """Topographic prominence: site elevation minus max of surrounding ring."""
    site_elev = get_elevation(lat, lon)
    max_ring = -1e9
    for i in range(n_ring):
        bearing = 2 * np.pi * i / n_ring
        dlat = (radius_km / 111.0) * np.cos(bearing)
        dlon = (radius_km / (111.0 * np.cos(np.radians(lat)))) * np.sin(bearing)
        e = get_elevation(lat + dlat, lon + dlon)
        if e > max_ring:
            max_ring = e
    return site_elev - max_ring

# ============================================================
# TEST 1: GREAT CIRCLE ALIGNMENT
# ============================================================
print("\n" + "="*60)
print("TEST 1: GREAT CIRCLE ALIGNMENT")
print("="*60)

THRESHOLD_KM = 100  # wider threshold due to small N

alignment_results = {
    "threshold_km": THRESHOLD_KM,
    "n_sites": N_SITES,
    "sites": [{"name": s["name"], "lat": s["lat"], "lon": s["lon"],
               "date_bce": s["date_bce"], "region": s["region"]} for s in PRE_AG_SITES],
}

# 1a) Alison Great Circle test
alison_dists = dist_from_gc_vec(site_lats, site_lons, ALISON_POLE_LAT, ALISON_POLE_LON)
near_alison = [(PRE_AG_SITES[i]["name"], float(alison_dists[i]))
               for i in range(N_SITES) if alison_dists[i] < THRESHOLD_KM]
print(f"\nAlison GC: {len(near_alison)} sites within {THRESHOLD_KM}km")
for name, d in near_alison:
    print(f"  {name}: {d:.1f} km from GC")

alignment_results["alison_gc"] = {
    "pole": [ALISON_POLE_LAT, ALISON_POLE_LON],
    "count": len(near_alison),
    "sites": [{"name": n, "dist_km": round(d, 1)} for n, d in near_alison],
}

# 1b) Top-10 high-D circles test
print(f"\nTop-10 high-D circles (threshold={THRESHOLD_KM}km):")
top10_results = []
for pole in TOP10_POLES:
    dists = dist_from_gc_vec(site_lats, site_lons, pole["lat"], pole["lon"])
    near = [(PRE_AG_SITES[i]["name"], float(dists[i]))
            for i in range(N_SITES) if dists[i] < THRESHOLD_KM]
    entry = {
        "name": pole["name"],
        "pole": [pole["lat"], pole["lon"]],
        "D": pole["D"],
        "count": len(near),
        "sites": [{"name": n, "dist_km": round(d, 1)} for n, d in near],
    }
    top10_results.append(entry)
    if near:
        print(f"  {pole['name']} (D={pole['D']:.1f}): {len(near)} sites — "
              f"{[n for n, _ in near]}")
    else:
        print(f"  {pole['name']} (D={pole['D']:.1f}): 0 sites")

alignment_results["top10_circles"] = top10_results

# 1c) Full pole scan (5° resolution)
print(f"\nFull pole scan (5° resolution, threshold={THRESHOLD_KM}km)...")
best_pole = None
best_count = 0
best_sites = []
scan_results = []

for lat in range(-85, 86, 5):
    for lon in range(-180, 180, 5):
        count = count_within(site_lats, site_lons, lat, lon, THRESHOLD_KM)
        if count > best_count:
            best_count = count
            best_pole = (lat, lon)
        if count >= 3:  # Record any pole catching 3+ sites
            dists = dist_from_gc_vec(site_lats, site_lons, lat, lon)
            near = [(PRE_AG_SITES[i]["name"], float(dists[i]))
                    for i in range(N_SITES) if dists[i] < THRESHOLD_KM]
            scan_results.append({
                "pole": [lat, lon],
                "count": count,
                "sites": [{"name": n, "dist_km": round(d, 1)} for n, d in near],
            })

# Get details for best pole
if best_pole:
    dists = dist_from_gc_vec(site_lats, site_lons, best_pole[0], best_pole[1])
    best_sites = [(PRE_AG_SITES[i]["name"], float(dists[i]))
                  for i in range(N_SITES) if dists[i] < THRESHOLD_KM]

print(f"Best pole: ({best_pole[0]}, {best_pole[1]}) with {best_count} sites within {THRESHOLD_KM}km")
for name, d in best_sites:
    print(f"  {name}: {d:.1f} km")

# 1d) Monte Carlo: how significant is the best count?
print(f"\nMonte Carlo: 10,000 random poles (threshold={THRESHOLD_KM}km)...")
N_MC = 10000
random_counts = np.empty(N_MC)
for i in range(N_MC):
    plat, plon = random_pole()
    random_counts[i] = count_within(site_lats, site_lons, plat, plon, THRESHOLD_KM)

p_value = float(np.mean(random_counts >= best_count))
mean_random = float(np.mean(random_counts))
std_random = float(np.std(random_counts))
z_score = (best_count - mean_random) / std_random if std_random > 0 else 0

print(f"Random pole mean: {mean_random:.2f} ± {std_random:.2f}")
print(f"Best scan count: {best_count}, Z-score: {z_score:.2f}, p-value: {p_value:.4f}")

# Also do a full scan with finer resolution near best pole
print(f"\nFine scan (1° resolution around best pole)...")
fine_best_pole = best_pole
fine_best_count = best_count
for lat in range(best_pole[0]-10, best_pole[0]+11, 1):
    for lon in range(best_pole[1]-10, best_pole[1]+11, 1):
        lat_c = max(-85, min(85, lat))
        count = count_within(site_lats, site_lons, lat_c, lon, THRESHOLD_KM)
        if count > fine_best_count:
            fine_best_count = count
            fine_best_pole = (lat_c, lon)

if fine_best_count > best_count:
    print(f"Fine scan improved: ({fine_best_pole[0]}, {fine_best_pole[1]}) with {fine_best_count} sites")
    dists = dist_from_gc_vec(site_lats, site_lons, fine_best_pole[0], fine_best_pole[1])
    fine_sites = [(PRE_AG_SITES[i]["name"], float(dists[i]))
                  for i in range(N_SITES) if dists[i] < THRESHOLD_KM]
    for name, d in fine_sites:
        print(f"  {name}: {d:.1f} km")
    best_pole = fine_best_pole
    best_count = fine_best_count
    best_sites = fine_sites
else:
    print("Fine scan did not improve.")

alignment_results["pole_scan"] = {
    "resolution_deg": 5,
    "best_pole": list(best_pole),
    "best_count": best_count,
    "best_sites": [{"name": n, "dist_km": round(d, 1)} for n, d in best_sites],
    "mc_trials": N_MC,
    "mc_mean": round(mean_random, 3),
    "mc_std": round(std_random, 3),
    "z_score": round(z_score, 3),
    "p_value": round(p_value, 4),
    "poles_with_3plus": len(scan_results),
    "scan_top_poles": sorted(scan_results, key=lambda x: -x["count"])[:10],
}

# Save alignment test
with open(os.path.join(OUT_DIR, "alignment_test.json"), 'w') as f:
    json.dump(alignment_results, f, indent=2)
print("\nSaved alignment_test.json")

# ============================================================
# TEST 2: ECOLOGICAL ZONE CLUSTERING
# ============================================================
print("\n" + "="*60)
print("TEST 2: ECOLOGICAL ZONE CLUSTERING")
print("="*60)

# Get elevations for all pre-ag sites
elevations = get_elevations_batch(site_lats, site_lons)
for i, s in enumerate(PRE_AG_SITES):
    s['elevation_m'] = float(elevations[i])

print(f"\nElevation range: {min(elevations):.0f} to {max(elevations):.0f} m")
print(f"Mean elevation: {np.mean(elevations):.0f} m")
print(f"Median elevation: {np.median(elevations):.0f} m")

# By region
ecological_results = {
    "global": {
        "n_sites": N_SITES,
        "elevation_range": [float(min(elevations)), float(max(elevations))],
        "elevation_mean": round(float(np.mean(elevations)), 1),
        "elevation_median": round(float(np.median(elevations)), 1),
        "elevation_std": round(float(np.std(elevations)), 1),
    },
    "by_region": {},
    "sites": [],
}

print("\nBy region:")
for region in sorted(set(s['region'] for s in PRE_AG_SITES)):
    r_sites = [s for s in PRE_AG_SITES if s['region'] == region]
    r_elevs = np.array([s['elevation_m'] for s in r_sites])
    print(f"  {region} (n={len(r_sites)}):")
    print(f"    Elevation: {np.mean(r_elevs):.0f} ± {np.std(r_elevs):.0f} m "
          f"(range: {min(r_elevs):.0f}-{max(r_elevs):.0f})")

    # Distance to centroid (cluster tightness)
    c_lat = np.mean([s['lat'] for s in r_sites])
    c_lon = np.mean([s['lon'] for s in r_sites])
    dists = [haversine_km(s['lat'], s['lon'], c_lat, c_lon) for s in r_sites]
    print(f"    Cluster radius: {np.mean(dists):.0f} km (max: {max(dists):.0f} km)")

    ecological_results["by_region"][region] = {
        "n_sites": len(r_sites),
        "elevation_mean": round(float(np.mean(r_elevs)), 1),
        "elevation_std": round(float(np.std(r_elevs)), 1),
        "elevation_range": [float(min(r_elevs)), float(max(r_elevs))],
        "centroid": [round(c_lat, 3), round(c_lon, 3)],
        "mean_dist_to_centroid_km": round(float(np.mean(dists)), 1),
        "max_dist_to_centroid_km": round(float(max(dists)), 1),
    }

# Per-site data
for s in PRE_AG_SITES:
    ecological_results["sites"].append({
        "name": s["name"],
        "lat": s["lat"],
        "lon": s["lon"],
        "date_bce": s["date_bce"],
        "region": s["region"],
        "elevation_m": s["elevation_m"],
    })

# Elevation band test: generate random land points and compare distributions
print("\nGenerating random land point elevations for comparison...")
n_random_land = 2000
land_lats = []
land_lons = []
attempts = 0
while len(land_lats) < n_random_land and attempts < 50000:
    attempts += 1
    lat = np.degrees(np.arcsin(np.random.uniform(-1, 1)))
    lon = np.random.uniform(-180, 180)
    elev = get_elevation(lat, lon)
    if elev > 0:  # on land
        land_lats.append(lat)
        land_lons.append(lon)

land_elevs = get_elevations_batch(np.array(land_lats), np.array(land_lons))
print(f"Got {len(land_elevs)} random land points (mean elev: {np.mean(land_elevs):.0f} m)")

# KS test: do pre-ag sites have different elevation distribution than random land?
from scipy.stats import ks_2samp, mannwhitneyu
ks_stat, ks_p = ks_2samp(elevations, land_elevs)
print(f"\nKS test (pre-ag vs random land): stat={ks_stat:.4f}, p={ks_p:.4f}")

# Mann-Whitney U test
mw_stat, mw_p = mannwhitneyu(elevations, land_elevs, alternative='two-sided')
print(f"Mann-Whitney U (pre-ag vs random land): stat={mw_stat:.0f}, p={mw_p:.4f}")

ecological_results["elevation_vs_random_land"] = {
    "n_random_land": len(land_elevs),
    "random_land_mean_elev": round(float(np.mean(land_elevs)), 1),
    "random_land_median_elev": round(float(np.median(land_elevs)), 1),
    "ks_stat": round(float(ks_stat), 4),
    "ks_p": round(float(ks_p), 4),
    "mw_stat": round(float(mw_stat), 1),
    "mw_p": round(float(mw_p), 4),
    "interpretation": "p < 0.05 means pre-ag sites have significantly different elevation distribution than random land"
}

# Wild cereal proximity (for Near East sites)
# Wild wheat/barley progenitor range roughly: 33-38°N, 35-45°E
WILD_CEREAL_CENTER = (35.5, 40.0)
print("\nNear East sites: distance to wild cereal progenitor range center:")
ne_sites = [s for s in PRE_AG_SITES if s['region'] == 'Near East']
for s in ne_sites:
    d = haversine_km(s['lat'], s['lon'], WILD_CEREAL_CENTER[0], WILD_CEREAL_CENTER[1])
    s['wild_cereal_dist_km'] = float(d)
    # Also check if within the range box
    in_range = (33 <= s['lat'] <= 38) and (35 <= s['lon'] <= 45)
    s['in_wild_cereal_range'] = in_range
    print(f"  {s['name']}: {d:.0f} km, in range: {in_range}")

ne_in_range = sum(1 for s in ne_sites if s.get('in_wild_cereal_range', False))
print(f"Near East sites in wild cereal range: {ne_in_range}/{len(ne_sites)} "
      f"({100*ne_in_range/len(ne_sites):.0f}%)")

ecological_results["wild_cereal_proximity"] = {
    "center": list(WILD_CEREAL_CENTER),
    "range_box": "33-38N, 35-45E",
    "ne_sites_in_range": ne_in_range,
    "ne_sites_total": len(ne_sites),
    "pct_in_range": round(100*ne_in_range/len(ne_sites), 1),
}

with open(os.path.join(OUT_DIR, "ecological_clustering.json"), 'w') as f:
    json.dump(ecological_results, f, indent=2)
print("\nSaved ecological_clustering.json")

# ============================================================
# TEST 3: TERRAIN SIGNATURE COMPARISON
# ============================================================
print("\n" + "="*60)
print("TEST 3: TERRAIN SIGNATURE COMPARISON")
print("="*60)

# Compute TPI and prominence for pre-ag sites
print("\nComputing TPI and prominence for pre-ag sites...")
pre_ag_tpi = []
pre_ag_prominence = []
for s in PRE_AG_SITES:
    tpi = compute_tpi(s['lat'], s['lon'], radius_km=10)
    prom = compute_prominence(s['lat'], s['lon'], radius_km=25)
    pre_ag_tpi.append(tpi)
    pre_ag_prominence.append(prom)
    s['tpi_10km'] = float(tpi)
    s['prominence_25km'] = float(prom)

pre_ag_tpi = np.array(pre_ag_tpi)
pre_ag_prominence = np.array(pre_ag_prominence)

print(f"Pre-ag TPI (10km): mean={np.mean(pre_ag_tpi):.1f}, median={np.median(pre_ag_tpi):.1f}")
print(f"Pre-ag Prominence (25km): mean={np.mean(pre_ag_prominence):.1f}, median={np.median(pre_ag_prominence):.1f}")

# Generate "civilizational" monument terrain stats using Pleiades monuments
# Load a sample of Pleiades monumental sites for comparison
print("\nLoading Pleiades monumental sites for comparison...")
import csv

MONUMENTAL_TYPES = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus", "shrine"
}

pleiades_path = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")
civ_lats = []
civ_lons = []
with open(pleiades_path, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            lat = float(row.get('reprLat', ''))
            lon = float(row.get('reprLong', ''))
        except (ValueError, TypeError):
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue
        ft = row.get('featureTypes', '').lower()
        types = set(ft.replace('"', '').replace('[', '').replace(']', '').split(','))
        types = {t.strip() for t in types}
        if types & MONUMENTAL_TYPES:
            civ_lats.append(lat)
            civ_lons.append(lon)

print(f"Pleiades monumental sites: {len(civ_lats)}")

# Sample 200 for TPI/prominence comparison (computing for all is slow)
n_civ_sample = min(200, len(civ_lats))
civ_idx = np.random.choice(len(civ_lats), n_civ_sample, replace=False)
civ_sample_lats = [civ_lats[i] for i in civ_idx]
civ_sample_lons = [civ_lons[i] for i in civ_idx]

print(f"Computing TPI/prominence for {n_civ_sample} civilizational monument sample...")
civ_tpi = []
civ_prominence = []
for i in range(n_civ_sample):
    tpi = compute_tpi(civ_sample_lats[i], civ_sample_lons[i], radius_km=10)
    prom = compute_prominence(civ_sample_lats[i], civ_sample_lons[i], radius_km=25)
    civ_tpi.append(tpi)
    civ_prominence.append(prom)
    if (i+1) % 50 == 0:
        print(f"  {i+1}/{n_civ_sample} done")

civ_tpi = np.array(civ_tpi)
civ_prominence = np.array(civ_prominence)

print(f"\nCivilizational TPI (10km): mean={np.mean(civ_tpi):.1f}, median={np.median(civ_tpi):.1f}")
print(f"Civilizational Prominence (25km): mean={np.mean(civ_prominence):.1f}, median={np.median(civ_prominence):.1f}")

# Statistical tests
tpi_stat, tpi_p = mannwhitneyu(pre_ag_tpi, civ_tpi, alternative='two-sided')
prom_stat, prom_p = mannwhitneyu(pre_ag_prominence, civ_prominence, alternative='two-sided')
print(f"\nMann-Whitney U (TPI, pre-ag vs civ): U={tpi_stat:.0f}, p={tpi_p:.4f}")
print(f"Mann-Whitney U (Prominence, pre-ag vs civ): U={prom_stat:.0f}, p={prom_p:.4f}")

# Also compute for random land points (subset)
n_rand_terrain = 200
rand_terrain_idx = np.random.choice(len(land_lats), n_rand_terrain, replace=False)
print(f"\nComputing TPI/prominence for {n_rand_terrain} random land points...")
rand_tpi = []
rand_prominence = []
for i in rand_terrain_idx:
    tpi = compute_tpi(land_lats[i], land_lons[i], radius_km=10)
    prom = compute_prominence(land_lats[i], land_lons[i], radius_km=25)
    rand_tpi.append(tpi)
    rand_prominence.append(prom)
    if (len(rand_tpi)) % 50 == 0:
        print(f"  {len(rand_tpi)}/{n_rand_terrain} done")

rand_tpi = np.array(rand_tpi)
rand_prominence = np.array(rand_prominence)
print(f"Random land TPI: mean={np.mean(rand_tpi):.1f}")

terrain_results = {
    "pre_ag": {
        "n": N_SITES,
        "tpi_10km": {
            "mean": round(float(np.mean(pre_ag_tpi)), 1),
            "median": round(float(np.median(pre_ag_tpi)), 1),
            "std": round(float(np.std(pre_ag_tpi)), 1),
            "values": [round(float(x), 1) for x in pre_ag_tpi],
        },
        "prominence_25km": {
            "mean": round(float(np.mean(pre_ag_prominence)), 1),
            "median": round(float(np.median(pre_ag_prominence)), 1),
            "std": round(float(np.std(pre_ag_prominence)), 1),
            "values": [round(float(x), 1) for x in pre_ag_prominence],
        },
    },
    "civilizational": {
        "n": n_civ_sample,
        "tpi_10km": {
            "mean": round(float(np.mean(civ_tpi)), 1),
            "median": round(float(np.median(civ_tpi)), 1),
            "std": round(float(np.std(civ_tpi)), 1),
        },
        "prominence_25km": {
            "mean": round(float(np.mean(civ_prominence)), 1),
            "median": round(float(np.median(civ_prominence)), 1),
            "std": round(float(np.std(civ_prominence)), 1),
        },
    },
    "random_land": {
        "n": n_rand_terrain,
        "tpi_10km": {
            "mean": round(float(np.mean(rand_tpi)), 1),
            "median": round(float(np.median(rand_tpi)), 1),
        },
        "prominence_25km": {
            "mean": round(float(np.mean(rand_prominence)), 1),
            "median": round(float(np.median(rand_prominence)), 1),
        },
    },
    "tests": {
        "pre_ag_vs_civ_tpi": {
            "U": round(float(tpi_stat), 1),
            "p": round(float(tpi_p), 4),
        },
        "pre_ag_vs_civ_prominence": {
            "U": round(float(prom_stat), 1),
            "p": round(float(prom_p), 4),
        },
    },
    "per_site": [
        {"name": s["name"], "tpi_10km": s["tpi_10km"], "prominence_25km": s["prominence_25km"],
         "elevation_m": s.get("elevation_m", None)}
        for s in PRE_AG_SITES
    ],
}

with open(os.path.join(OUT_DIR, "terrain_comparison.json"), 'w') as f:
    json.dump(terrain_results, f, indent=2)
print("\nSaved terrain_comparison.json")

# ============================================================
# TEST 4: TAŞ TEPELER SPECIFIC ANALYSIS
# ============================================================
print("\n" + "="*60)
print("TEST 4: TAŞ TEPELER CLUSTER ANALYSIS")
print("="*60)

# Taş Tepeler + nearby PPNA sites (everything in the upland SE Turkey region)
tas_tepeler = [s for s in PRE_AG_SITES
               if s['region'] == 'Near East' and s['lat'] > 36 and s['lon'] > 37 and s['lon'] < 42]
print(f"\nTaş Tepeler region sites: {len(tas_tepeler)}")
for s in tas_tepeler:
    print(f"  {s['name']}: ({s['lat']:.3f}, {s['lon']:.3f}), "
          f"elev={s.get('elevation_m', '?')}m, {s['date_bce']} BCE")

# Centroid
c_lat = np.mean([s['lat'] for s in tas_tepeler])
c_lon = np.mean([s['lon'] for s in tas_tepeler])
print(f"\nCentroid: ({c_lat:.3f}, {c_lon:.3f})")

# Standard deviational ellipse
coords = np.array([[s['lon'], s['lat']] for s in tas_tepeler])
cov = np.cov(coords.T)
eigenvalues, eigenvectors = np.linalg.eig(cov)
major_idx = np.argmax(eigenvalues)
minor_idx = 1 - major_idx
major_axis = eigenvectors[:, major_idx]
minor_axis = eigenvectors[:, minor_idx]

# Major axis bearing (degrees from north, clockwise)
bearing = np.degrees(np.arctan2(major_axis[0], major_axis[1]))
if bearing < 0:
    bearing += 360
print(f"Major axis bearing: {bearing:.1f}° (from north)")
print(f"Eigenvalue ratio (elongation): {eigenvalues[major_idx]/eigenvalues[minor_idx]:.2f}")
print(f"Major semi-axis (deg): {np.sqrt(eigenvalues[major_idx]):.3f}")
print(f"Minor semi-axis (deg): {np.sqrt(eigenvalues[minor_idx]):.3f}")

# Convert to approximate km
major_km = np.sqrt(eigenvalues[major_idx]) * 111.0  # rough deg-to-km
minor_km = np.sqrt(eigenvalues[minor_idx]) * 111.0
print(f"Major semi-axis: ~{major_km:.0f} km, Minor: ~{minor_km:.0f} km")

# What bearing does the Alison GC have at this location?
# The GC bearing at a point is perpendicular to the bearing TO the pole
# Bearing from centroid to Alison pole:
dlat_r = np.radians(ALISON_POLE_LAT - c_lat)
dlon_r = np.radians(ALISON_POLE_LON - c_lon)
y = np.sin(dlon_r) * np.cos(np.radians(ALISON_POLE_LAT))
x = (np.cos(np.radians(c_lat)) * np.sin(np.radians(ALISON_POLE_LAT)) -
     np.sin(np.radians(c_lat)) * np.cos(np.radians(ALISON_POLE_LAT)) * np.cos(dlon_r))
bearing_to_pole = np.degrees(np.arctan2(y, x)) % 360
gc_bearing = (bearing_to_pole + 90) % 360  # GC is perpendicular to pole direction
print(f"\nAlison GC bearing at centroid: ~{gc_bearing:.1f}°")
print(f"Cluster major axis bearing: {bearing:.1f}°")
bearing_diff = min(abs(bearing - gc_bearing), 360 - abs(bearing - gc_bearing))
if bearing_diff > 90:
    bearing_diff = 180 - bearing_diff
print(f"Angular difference (cluster vs GC): {bearing_diff:.1f}°")

# Elevation analysis of Taş Tepeler
tas_elevs = np.array([s.get('elevation_m', 0) for s in tas_tepeler])
print(f"\nTaş Tepeler elevations:")
print(f"  Range: {min(tas_elevs):.0f} - {max(tas_elevs):.0f} m")
print(f"  Mean: {np.mean(tas_elevs):.0f} m, Std: {np.std(tas_elevs):.0f} m")

# Are they in a narrow elevation band?
elev_cv = np.std(tas_elevs) / np.mean(tas_elevs) if np.mean(tas_elevs) > 0 else float('inf')
print(f"  Coefficient of variation: {elev_cv:.3f}")
print(f"  {'TIGHT' if elev_cv < 0.3 else 'SPREAD'} elevation band")

# Elevation contour bearing at the centroid
# Sample elevations along E-W and N-S transects to estimate local gradient
gradient_step = 0.1  # degrees
elev_e = get_elevation(c_lat, c_lon + gradient_step)
elev_w = get_elevation(c_lat, c_lon - gradient_step)
elev_n = get_elevation(c_lat + gradient_step, c_lon)
elev_s = get_elevation(c_lat - gradient_step, c_lon)

grad_ew = (elev_e - elev_w) / (2 * gradient_step * 111.0 * np.cos(np.radians(c_lat)))
grad_ns = (elev_n - elev_s) / (2 * gradient_step * 111.0)

# Gradient bearing (steepest ascent direction)
gradient_bearing = np.degrees(np.arctan2(grad_ew, grad_ns)) % 360
# Contour bearing is perpendicular to gradient
contour_bearing = (gradient_bearing + 90) % 360
print(f"\nLocal elevation gradient bearing: {gradient_bearing:.1f}°")
print(f"Local contour (isohypse) bearing: {contour_bearing:.1f}°")

contour_diff = min(abs(bearing - contour_bearing), 360 - abs(bearing - contour_bearing))
if contour_diff > 90:
    contour_diff = 180 - contour_diff
print(f"Angular difference (cluster vs contour): {contour_diff:.1f}°")

# Summary
print(f"\n--- Taş Tepeler Orientation Summary ---")
print(f"Cluster major axis: {bearing:.1f}°")
print(f"Alison GC bearing:  {gc_bearing:.1f}° (diff: {bearing_diff:.1f}°)")
print(f"Elevation contour:  {contour_bearing:.1f}° (diff: {contour_diff:.1f}°)")
follows_contour = contour_diff < bearing_diff
print(f"Cluster follows: {'ELEVATION CONTOUR' if follows_contour else 'GREAT CIRCLE BEARING'} more closely")

tas_results = {
    "sites": [{"name": s["name"], "lat": s["lat"], "lon": s["lon"],
               "date_bce": s["date_bce"], "elevation_m": s.get("elevation_m", None)}
              for s in tas_tepeler],
    "n_sites": len(tas_tepeler),
    "centroid": [round(c_lat, 3), round(c_lon, 3)],
    "major_axis_bearing_deg": round(bearing, 1),
    "eigenvalue_ratio": round(float(eigenvalues[major_idx]/eigenvalues[minor_idx]), 2),
    "major_semi_axis_km": round(major_km, 0),
    "minor_semi_axis_km": round(minor_km, 0),
    "alison_gc_bearing_deg": round(gc_bearing, 1),
    "contour_bearing_deg": round(contour_bearing, 1),
    "diff_cluster_vs_gc": round(bearing_diff, 1),
    "diff_cluster_vs_contour": round(contour_diff, 1),
    "follows": "contour" if follows_contour else "gc_bearing",
    "elevation_range_m": [float(min(tas_elevs)), float(max(tas_elevs))],
    "elevation_mean_m": round(float(np.mean(tas_elevs)), 0),
    "elevation_std_m": round(float(np.std(tas_elevs)), 0),
    "elevation_cv": round(elev_cv, 3),
    "elevation_band": "tight" if elev_cv < 0.3 else "spread",
}

with open(os.path.join(OUT_DIR, "tas_tepeler_analysis.json"), 'w') as f:
    json.dump(tas_results, f, indent=2)
print("\nSaved tas_tepeler_analysis.json")

# ============================================================
# PLOTS
# ============================================================
print("\n" + "="*60)
print("GENERATING PLOTS")
print("="*60)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- Plot 1: Site map color-coded by region ---
fig, ax = plt.subplots(1, 1, figsize=(14, 7))
region_colors = {
    "Near East": "#d62728",
    "North America": "#1f77b4",
    "South America": "#2ca02c",
    "East Asia": "#ff7f0e",
    "Europe": "#9467bd",
}
for region, color in region_colors.items():
    r_sites = [s for s in PRE_AG_SITES if s['region'] == region]
    ax.scatter([s['lon'] for s in r_sites], [s['lat'] for s in r_sites],
              c=color, label=f"{region} (n={len(r_sites)})", s=60, edgecolors='k',
              linewidths=0.5, zorder=5)

# Draw Alison GC
gc_bearings = np.linspace(0, 2*np.pi, 720)
pole_lat_r = np.radians(ALISON_POLE_LAT)
pole_lon_r = np.radians(ALISON_POLE_LON)
gc_lats = []
gc_lons = []
for b in gc_bearings:
    lat = np.arcsin(np.sin(pole_lat_r)*np.cos(np.pi/2) +
                    np.cos(pole_lat_r)*np.sin(np.pi/2)*np.cos(b))
    lon = pole_lon_r + np.arctan2(np.sin(b)*np.sin(np.pi/2)*np.cos(pole_lat_r),
                                   np.cos(np.pi/2) - np.sin(pole_lat_r)*np.sin(lat))
    gc_lats.append(np.degrees(lat))
    gc_lons.append(np.degrees(lon))
gc_lons = [(l + 180) % 360 - 180 for l in gc_lons]  # normalize

# Sort by longitude for cleaner plotting, split at wraps
gc_points = sorted(zip(gc_lons, gc_lats))
# Split at large gaps
segments = []
current_seg = [gc_points[0]]
for i in range(1, len(gc_points)):
    if abs(gc_points[i][0] - gc_points[i-1][0]) > 30:
        segments.append(current_seg)
        current_seg = [gc_points[i]]
    else:
        current_seg.append(gc_points[i])
segments.append(current_seg)
for seg in segments:
    ax.plot([p[0] for p in seg], [p[1] for p in seg], 'r-', alpha=0.3, linewidth=1)

ax.set_xlim(-180, 180)
ax.set_ylim(-60, 70)
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.set_title(f"Pre-Agricultural Monumental Sites (n={N_SITES}) with Alison Great Circle")
ax.legend(loc='lower left', fontsize=8)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "site_map.png"), dpi=150)
plt.close()
print("Saved site_map.png")

# --- Plot 2: Terrain comparison boxplots ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# TPI boxplot
data_tpi = [pre_ag_tpi, civ_tpi, rand_tpi]
labels_tpi = [f"Pre-Ag\n(n={N_SITES})", f"Civilizational\n(n={n_civ_sample})", f"Random Land\n(n={n_rand_terrain})"]
bp1 = ax1.boxplot(data_tpi, labels=labels_tpi, patch_artist=True)
colors = ['#d62728', '#1f77b4', '#7f7f7f']
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
ax1.axhline(0, color='k', linestyle='--', alpha=0.3)
ax1.set_ylabel("TPI (m)")
ax1.set_title(f"Topographic Position Index (10km radius)\np={tpi_p:.4f}")
ax1.grid(True, alpha=0.3, axis='y')

# Prominence boxplot
data_prom = [pre_ag_prominence, civ_prominence, rand_prominence]
bp2 = ax2.boxplot(data_prom, labels=labels_tpi, patch_artist=True)
for patch, color in zip(bp2['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
ax2.axhline(0, color='k', linestyle='--', alpha=0.3)
ax2.set_ylabel("Prominence (m)")
ax2.set_title(f"Topographic Prominence (25km radius)\np={prom_p:.4f}")
ax2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "terrain_comparison.png"), dpi=150)
plt.close()
print("Saved terrain_comparison.png")

# --- Plot 3: Elevation distribution by region ---
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
regions = sorted(set(s['region'] for s in PRE_AG_SITES))
region_elevs = []
region_labels = []
for region in regions:
    r_elevs = [s['elevation_m'] for s in PRE_AG_SITES if s['region'] == region]
    region_elevs.append(r_elevs)
    region_labels.append(f"{region}\n(n={len(r_elevs)})")

bp3 = ax.boxplot(region_elevs, labels=region_labels, patch_artist=True)
for patch, color in zip(bp3['boxes'], [region_colors.get(r, '#333') for r in regions]):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
ax.set_ylabel("Elevation (m)")
ax.set_title("Elevation Distribution by Region — Pre-Agricultural Sites")
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "elevation_by_region.png"), dpi=150)
plt.close()
print("Saved elevation_by_region.png")

# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "="*60)
print("SUMMARY")
print("="*60)

print(f"\nTEST 1 — GREAT CIRCLE ALIGNMENT:")
print(f"  Alison GC: {len(near_alison)} sites within {THRESHOLD_KM}km")
print(f"  Best from top-10: {max(r['count'] for r in top10_results)} sites")
print(f"  Best from full scan: {best_count} sites at pole ({best_pole[0]}, {best_pole[1]})")
print(f"  MC Z-score: {z_score:.2f}, p={p_value:.4f}")
gc_significant = p_value < 0.05
print(f"  Verdict: {'SIGNIFICANT' if gc_significant else 'NOT SIGNIFICANT'} alignment")

print(f"\nTEST 2 — ECOLOGICAL CLUSTERING:")
print(f"  KS test (elevation): p={ks_p:.4f}")
print(f"  Near East in wild cereal range: {ne_in_range}/{len(ne_sites)} ({100*ne_in_range/len(ne_sites):.0f}%)")

print(f"\nTEST 3 — TERRAIN SIGNATURE:")
print(f"  Pre-ag TPI: {np.mean(pre_ag_tpi):.1f}m (civ: {np.mean(civ_tpi):.1f}m)")
print(f"  Pre-ag Prominence: {np.mean(pre_ag_prominence):.1f}m (civ: {np.mean(civ_prominence):.1f}m)")
print(f"  TPI difference p={tpi_p:.4f}, Prominence difference p={prom_p:.4f}")
tpi_higher = np.mean(pre_ag_tpi) > np.mean(civ_tpi)
print(f"  Pre-ag sites are {'HIGHER' if tpi_higher else 'LOWER'} relative to surroundings")

print(f"\nTEST 4 — TAŞ TEPELER:")
print(f"  Cluster bearing: {bearing:.1f}°")
print(f"  GC bearing: {gc_bearing:.1f}° (diff: {bearing_diff:.1f}°)")
print(f"  Contour bearing: {contour_bearing:.1f}° (diff: {contour_diff:.1f}°)")
print(f"  Cluster follows: {'CONTOUR' if follows_contour else 'GC'}")
print(f"  Elevation band: {min(tas_elevs):.0f}-{max(tas_elevs):.0f}m (CV={elev_cv:.3f})")

# Determine overall verdict
if not gc_significant and ks_p < 0.05:
    verdict = "No GC alignment + ecological clustering + different TPI"
    interpretation = ("Pre-ag monuments follow ecological logic, not geometric. "
                      "Alignment patterns are civilizational.")
elif gc_significant and ks_p < 0.05:
    verdict = "Weak GC alignment + ecological clustering"
    interpretation = "Partially geographic but mostly ecological"
elif gc_significant:
    verdict = "Strong GC alignment"
    interpretation = ("The alignment pattern extends before agriculture — "
                      "would be a major finding")
else:
    verdict = "No clear pattern"
    interpretation = "Neither alignment nor ecological clustering detected"

print(f"\n{'='*60}")
print(f"OVERALL VERDICT: {verdict}")
print(f"INTERPRETATION: {interpretation}")
print(f"{'='*60}")

# Save all results to master JSON
master = {
    "test_date": "2026-03-21",
    "n_sites": N_SITES,
    "threshold_km": THRESHOLD_KM,
    "verdict": verdict,
    "interpretation": interpretation,
    "tests": {
        "gc_alignment": {
            "alison_count": len(near_alison),
            "best_scan_count": best_count,
            "best_scan_pole": list(best_pole),
            "z_score": round(z_score, 3),
            "p_value": round(p_value, 4),
            "significant": gc_significant,
        },
        "ecological_clustering": {
            "ks_p": round(float(ks_p), 4),
            "mw_p": round(float(mw_p), 4),
            "ne_in_wild_cereal_range_pct": round(100*ne_in_range/len(ne_sites), 1),
        },
        "terrain_signature": {
            "pre_ag_tpi_mean": round(float(np.mean(pre_ag_tpi)), 1),
            "civ_tpi_mean": round(float(np.mean(civ_tpi)), 1),
            "tpi_p": round(float(tpi_p), 4),
            "pre_ag_prominence_mean": round(float(np.mean(pre_ag_prominence)), 1),
            "civ_prominence_mean": round(float(np.mean(civ_prominence)), 1),
            "prominence_p": round(float(prom_p), 4),
            "pre_ag_higher": bool(tpi_higher),
        },
        "tas_tepeler": {
            "cluster_bearing": round(bearing, 1),
            "gc_bearing": round(gc_bearing, 1),
            "contour_bearing": round(contour_bearing, 1),
            "follows": "contour" if follows_contour else "gc_bearing",
            "elevation_band": "tight" if elev_cv < 0.3 else "spread",
        },
    },
}

with open(os.path.join(OUT_DIR, "master_results.json"), 'w') as f:
    json.dump(master, f, indent=2)

print("\nAll outputs saved to:", OUT_DIR)
print("Done!")
