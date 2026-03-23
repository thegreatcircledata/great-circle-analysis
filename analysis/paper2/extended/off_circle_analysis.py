#!/usr/bin/env python3
"""
Directive 2: Off-Circle Major Sites Analysis
=============================================
Why are Göbekli Tepe, Stonehenge, Angkor Wat, Teotihuacan, Great Zimbabwe,
and other famous monuments NOT on the Alison circle? Do they follow different
geographic logic?

Computes:
  1. Distance to Alison GC for 50 famous off-circle sites
  2. Terrain features: elevation, TPI, coastline distance, river distance
  3. Which of the 12 high-D clusters each off-circle site falls on
  4. Comparison of terrain signatures: off-circle vs on-circle
  5. Search for a "second alignment" through off-circle sites
"""

import json, math, os, sys, time
import numpy as np
from scipy.stats import mannwhitneyu, percentileofscore
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

# ============================================================
# CONFIGURATION
# ============================================================
POLE_LAT = 59.682122
POLE_LON = -138.646087
R = 6371.0
QUARTER_CIRC = R * math.pi / 2  # ~10,007.5 km

BASE = Path(os.path.expanduser("~/megalith_site_research"))
OUT = BASE / "outputs" / "extended_analysis"
OUT.mkdir(parents=True, exist_ok=True)

ETOPO_PATH = BASE / "data" / "geophysical" / "etopo" / "ETOPO1_Ice_g_gmt4.grd"
NE_DIR = BASE / "data" / "natural_earth"

# ============================================================
# 50 FAMOUS OFF-CIRCLE SITES (>100km from Alison GC)
# ============================================================
OFF_CIRCLE_SITES = [
    {"name": "Göbekli Tepe", "lat": 37.223, "lon": 38.922, "region": "Middle East", "age_bce": -9500},
    {"name": "Stonehenge", "lat": 51.179, "lon": -1.826, "region": "Europe", "age_bce": -3000},
    {"name": "Angkor Wat", "lat": 13.412, "lon": 103.867, "region": "SE Asia", "age_bce": 1150},
    {"name": "Teotihuacan", "lat": 19.692, "lon": -98.844, "region": "Mesoamerica", "age_bce": -200},
    {"name": "Great Zimbabwe", "lat": -20.267, "lon": 30.933, "region": "Africa", "age_bce": 1100},
    {"name": "Borobudur", "lat": -7.608, "lon": 110.204, "region": "SE Asia", "age_bce": 800},
    {"name": "Cahokia", "lat": 38.655, "lon": -90.062, "region": "North America", "age_bce": 1050},
    {"name": "Chichen Itza", "lat": 20.683, "lon": -88.569, "region": "Mesoamerica", "age_bce": 600},
    {"name": "Tikal", "lat": 17.222, "lon": -89.623, "region": "Mesoamerica", "age_bce": -600},
    {"name": "Karnak", "lat": 25.719, "lon": 32.657, "region": "Egypt", "age_bce": -2000},
    {"name": "Petra", "lat": 30.329, "lon": 35.444, "region": "Middle East", "age_bce": -300},
    {"name": "Newgrange", "lat": 53.695, "lon": -6.475, "region": "Europe", "age_bce": -3200},
    {"name": "Carnac", "lat": 47.584, "lon": -3.078, "region": "Europe", "age_bce": -4500},
    {"name": "Troy", "lat": 39.957, "lon": 26.239, "region": "Europe", "age_bce": -3000},
    {"name": "Mycenae", "lat": 37.731, "lon": 22.756, "region": "Europe", "age_bce": -1600},
    {"name": "Knossos", "lat": 35.298, "lon": 25.163, "region": "Europe", "age_bce": -2000},
    {"name": "Pompeii", "lat": 40.751, "lon": 14.485, "region": "Europe", "age_bce": -700},
    {"name": "Colosseum", "lat": 41.890, "lon": 12.492, "region": "Europe", "age_bce": 72},
    {"name": "Parthenon", "lat": 37.972, "lon": 23.727, "region": "Europe", "age_bce": -447},
    {"name": "Delphi", "lat": 38.482, "lon": 22.501, "region": "Europe", "age_bce": -800},
    {"name": "Baalbek", "lat": 34.007, "lon": 36.204, "region": "Middle East", "age_bce": -2900},
    {"name": "Persepolis", "lat": 29.935, "lon": 52.891, "region": "Middle East", "age_bce": -515},
    {"name": "Mohenjo-daro", "lat": 27.329, "lon": 68.136, "region": "South Asia", "age_bce": -2500},
    {"name": "Harappa", "lat": 30.631, "lon": 72.864, "region": "South Asia", "age_bce": -2600},
    {"name": "Sigiriya", "lat": 7.957, "lon": 80.760, "region": "South Asia", "age_bce": 477},
    {"name": "Prambanan", "lat": -7.752, "lon": 110.491, "region": "SE Asia", "age_bce": 850},
    {"name": "Sukhothai", "lat": 17.014, "lon": 99.703, "region": "SE Asia", "age_bce": 1238},
    {"name": "Bagan", "lat": 21.172, "lon": 94.868, "region": "SE Asia", "age_bce": 849},
    {"name": "Mesa Verde", "lat": 37.184, "lon": -108.488, "region": "North America", "age_bce": 550},
    {"name": "Chaco Canyon", "lat": 36.060, "lon": -107.961, "region": "North America", "age_bce": 850},
    {"name": "Poverty Point", "lat": 32.637, "lon": -91.408, "region": "North America", "age_bce": -1700},
    {"name": "Tiwanaku", "lat": -16.556, "lon": -68.673, "region": "South America", "age_bce": -200},
    {"name": "Chan Chan", "lat": -8.106, "lon": -79.074, "region": "South America", "age_bce": 850},
    {"name": "Caral", "lat": -10.893, "lon": -77.520, "region": "South America", "age_bce": -3000},
    {"name": "Monte Albán", "lat": 17.043, "lon": -96.769, "region": "Mesoamerica", "age_bce": -500},
    {"name": "Palenque", "lat": 17.484, "lon": -92.046, "region": "Mesoamerica", "age_bce": 226},
    {"name": "Copán", "lat": 14.840, "lon": -89.141, "region": "Mesoamerica", "age_bce": -1200},
    {"name": "Great Enclosure (Hattusa)", "lat": 40.020, "lon": 34.614, "region": "Middle East", "age_bce": -1600},
    {"name": "Palmyra", "lat": 34.551, "lon": 38.267, "region": "Middle East", "age_bce": -2000},
    {"name": "Jericho", "lat": 31.871, "lon": 35.444, "region": "Middle East", "age_bce": -9000},
    {"name": "Catalhöyük", "lat": 37.666, "lon": 32.828, "region": "Middle East", "age_bce": -7500},
    {"name": "Leptis Magna", "lat": 32.638, "lon": 14.289, "region": "Africa", "age_bce": -1000},
    {"name": "Aksum", "lat": 14.131, "lon": 38.720, "region": "Africa", "age_bce": -400},
    {"name": "Lalibela", "lat": 12.030, "lon": 39.044, "region": "Africa", "age_bce": 1181},
    {"name": "Abu Simbel", "lat": 22.337, "lon": 31.626, "region": "Egypt", "age_bce": -1264},
    {"name": "Valley of the Kings", "lat": 25.740, "lon": 32.602, "region": "Egypt", "age_bce": -1539},
    {"name": "Saqqara", "lat": 29.871, "lon": 31.216, "region": "Egypt", "age_bce": -2650},
    {"name": "Derinkuyu", "lat": 38.374, "lon": 34.734, "region": "Middle East", "age_bce": -800},
    {"name": "Skara Brae", "lat": 59.049, "lon": -3.341, "region": "Europe", "age_bce": -3180},
    {"name": "Avebury", "lat": 51.429, "lon": -1.854, "region": "Europe", "age_bce": -2850},
]

# 12 high-D cluster poles from systematic search
CLUSTER_POLES = [
    {"id": 1, "name": "S California", "pole_lat": 36.0, "pole_lon": -117.5, "D": 24.28},
    {"id": 2, "name": "Sumatra/SE Asia", "pole_lat": 0.5, "pole_lon": 102.5, "D": 19.59},
    {"id": 3, "name": "Pacific NW (Alison family)", "pole_lat": 45.5, "pole_lon": -141.5, "D": 22.69},
    {"id": 4, "name": "Philippines/Taiwan", "pole_lat": 17.5, "pole_lon": 119.5, "D": 15.94},
    {"id": 5, "name": "NW Pacific", "pole_lat": 43.0, "pole_lon": 160.0, "D": 15.80},
    {"id": 6, "name": "N Alberta/Canada", "pole_lat": 60.5, "pole_lon": -117.5, "D": 9.95},
    {"id": 7, "name": "Bering", "pole_lat": 59.5, "pole_lon": 178.5, "D": 14.08},
    {"id": 8, "name": "Central America", "pole_lat": 6.0, "pole_lon": -82.5, "D": 16.73},
    {"id": 9, "name": "US East Coast", "pole_lat": 34.0, "pole_lon": -76.5, "D": 12.12},
    {"id": 10, "name": "Atlantic/Brazil", "pole_lat": 3.5, "pole_lon": -59.0, "D": 12.34},
    {"id": 11, "name": "Japan/Pacific", "pole_lat": 34.0, "pole_lon": 140.5, "D": 15.84},
    {"id": 12, "name": "Micronesia", "pole_lat": 11.5, "pole_lon": 141.5, "D": 10.53},
]


# ============================================================
# LOAD ETOPO1
# ============================================================
print("=" * 70)
print("DIRECTIVE 2: OFF-CIRCLE MAJOR SITES ANALYSIS")
print("=" * 70)

print("\nLoading ETOPO1...")
import netCDF4 as nc
ds = nc.Dataset(str(ETOPO_PATH))
etopo_lons = ds.variables['x'][:]
etopo_lats = ds.variables['y'][:]
etopo_z = ds.variables['z']
lon_min, lon_step = float(etopo_lons[0]), float(etopo_lons[1] - etopo_lons[0])
lat_min, lat_step = float(etopo_lats[0]), float(etopo_lats[1] - etopo_lats[0])
n_lats, n_lons = len(etopo_lats), len(etopo_lons)
print(f"  ETOPO1: {n_lats}x{n_lons}, step={lat_step:.4f}°")


def get_elev(lat, lon):
    li = int(round((lat - lat_min) / lat_step))
    lo = int(round((lon - lon_min) / lon_step))
    li = max(0, min(li, n_lats - 1))
    lo = max(0, min(lo, n_lons - 1))
    return float(etopo_z[li, lo])


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(min(1.0, a)))


def dist_from_gc(lat, lon, pole_lat, pole_lon):
    d = haversine_km(lat, lon, pole_lat, pole_lon)
    return abs(d - QUARTER_CIRC)


def compute_tpi(lat, lon, radius_km=10):
    """Topographic Position Index: elevation minus mean within radius."""
    site_elev = get_elev(lat, lon)
    if site_elev <= -500:  # deep ocean
        return np.nan
    dlat = radius_km / 111.0
    dlon = radius_km / (111.0 * max(math.cos(math.radians(lat)), 0.01))
    n_ring = 32
    surr = []
    for k in range(n_ring):
        angle = 2 * math.pi * k / n_ring
        slat = lat + dlat * math.sin(angle)
        slon = lon + dlon * math.cos(angle)
        surr.append(get_elev(slat, slon))
    return site_elev - np.mean(surr)


def compute_prominence(lat, lon, radius_km=25):
    """Topographic prominence: elevation minus mean surrounding elevation."""
    center_elev = get_elev(lat, lon)
    if center_elev <= 0:
        return np.nan
    n_samples = 64
    dlat = radius_km / 111.0
    dlon = radius_km / (111.0 * max(math.cos(math.radians(lat)), 0.01))
    surr_elevs = []
    for k in range(n_samples):
        angle = 2 * math.pi * k / n_samples
        slat = lat + dlat * math.sin(angle)
        slon = lon + dlon * math.cos(angle)
        e = get_elev(slat, slon)
        if e > 0:
            surr_elevs.append(e)
    if len(surr_elevs) < 10:
        return np.nan
    return center_elev - np.mean(surr_elevs)


# ============================================================
# LOAD COASTLINE AND RIVER DATA
# ============================================================
print("\nLoading coastline and river data...")
import geopandas as gpd
from shapely.geometry import MultiLineString, LineString

def load_line_features(shp_path, sample_interval_deg=0.5):
    gdf = gpd.read_file(shp_path)
    points = []
    for geom in gdf.geometry:
        if geom is None:
            continue
        lines = []
        if geom.geom_type == 'MultiLineString':
            lines = list(geom.geoms)
        elif geom.geom_type == 'LineString':
            lines = [geom]
        for line in lines:
            for c in line.coords:
                points.append((c[1], c[0]))
            length = line.length
            if length > sample_interval_deg:
                n_samples = int(length / sample_interval_deg)
                for i in range(1, n_samples):
                    pt = line.interpolate(i / n_samples, normalized=True)
                    points.append((pt.y, pt.x))
    return np.array(points)


def min_dist_to_features(qlat, qlon, feat_pts, max_deg=10.0):
    lat_mask = np.abs(feat_pts[:, 0] - qlat) < max_deg
    cands = feat_pts[lat_mask]
    if len(cands) == 0:
        return float('inf')
    lon_diff = np.abs(cands[:, 1] - qlon)
    lon_diff = np.minimum(lon_diff, 360 - lon_diff)
    lon_mask = lon_diff < max_deg
    cands = cands[lon_mask]
    if len(cands) == 0:
        return float('inf')
    lat1, lon1 = np.radians(qlat), np.radians(qlon)
    lat2, lon2 = np.radians(cands[:, 0]), np.radians(cands[:, 1])
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    dists = R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return float(np.min(dists))


coast_shp = NE_DIR / "ne_50m_coastline.shp"
river_shp = NE_DIR / "ne_50m_rivers_lake_centerlines.shp"

if coast_shp.exists() and river_shp.exists():
    coast_pts = load_line_features(str(coast_shp), 0.2)
    river_pts = load_line_features(str(river_shp), 0.3)
    print(f"  Coastline: {len(coast_pts)} pts, Rivers: {len(river_pts)} pts")
    HAS_GEO = True
else:
    print("  WARNING: Natural Earth shapefiles not found — skipping coast/river distances")
    HAS_GEO = False


# ============================================================
# STEP 1: Compute features for each off-circle site
# ============================================================
print("\n" + "=" * 70)
print("STEP 1: Computing features for 50 off-circle sites")
print("=" * 70)

off_circle_results = []
for i, site in enumerate(OFF_CIRCLE_SITES):
    lat, lon = site["lat"], site["lon"]

    # Distance to Alison GC
    alison_dist = dist_from_gc(lat, lon, POLE_LAT, POLE_LON)

    # Elevation
    elev = get_elev(lat, lon)

    # TPI (10km radius)
    tpi = compute_tpi(lat, lon, radius_km=10)

    # Prominence (25km radius)
    prom = compute_prominence(lat, lon, radius_km=25)

    # Coast/river distances
    if HAS_GEO:
        d_coast = min_dist_to_features(lat, lon, coast_pts)
        d_river = min_dist_to_features(lat, lon, river_pts)
    else:
        d_coast, d_river = np.nan, np.nan

    # Check against all 12 cluster circles
    cluster_dists = {}
    nearest_cluster = None
    nearest_cluster_dist = float('inf')
    for cp in CLUSTER_POLES:
        cd = dist_from_gc(lat, lon, cp["pole_lat"], cp["pole_lon"])
        cluster_dists[cp["name"]] = round(cd, 1)
        if cd < nearest_cluster_dist:
            nearest_cluster_dist = cd
            nearest_cluster = cp["name"]

    # How many clusters is this site within 100km of?
    on_clusters = [cp["name"] for cp in CLUSTER_POLES
                   if dist_from_gc(lat, lon, cp["pole_lat"], cp["pole_lon"]) < 100]

    result = {
        "name": site["name"],
        "lat": lat,
        "lon": lon,
        "region": site["region"],
        "age_bce": site["age_bce"],
        "alison_gc_dist_km": round(alison_dist, 1),
        "elevation_m": round(elev, 1),
        "tpi_10km_m": round(float(tpi), 1) if not np.isnan(tpi) else None,
        "prominence_25km_m": round(float(prom), 1) if not np.isnan(prom) else None,
        "coastline_dist_km": round(d_coast, 1) if not np.isnan(d_coast) else None,
        "river_dist_km": round(d_river, 1) if not np.isnan(d_river) else None,
        "nearest_cluster": nearest_cluster,
        "nearest_cluster_dist_km": round(nearest_cluster_dist, 1),
        "on_clusters_100km": on_clusters,
        "cluster_distances": cluster_dists,
    }
    off_circle_results.append(result)

    if (i + 1) % 10 == 0:
        print(f"  {i+1}/50 done")

print(f"  All 50 sites processed")

# Verify all are truly off-circle (>100km from Alison)
on_circle_count = sum(1 for r in off_circle_results if r["alison_gc_dist_km"] < 100)
if on_circle_count > 0:
    print(f"\n  WARNING: {on_circle_count} sites are actually within 100km of Alison GC!")
    for r in off_circle_results:
        if r["alison_gc_dist_km"] < 100:
            print(f"    {r['name']}: {r['alison_gc_dist_km']} km")


# ============================================================
# STEP 2: On-circle terrain reference (from GC land points)
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Computing on-circle terrain reference")
print("=" * 70)

def trace_great_circle(pole_lat, pole_lon, step_deg=1.0):
    pole = np.array([
        np.cos(np.radians(pole_lat)) * np.cos(np.radians(pole_lon)),
        np.cos(np.radians(pole_lat)) * np.sin(np.radians(pole_lon)),
        np.sin(np.radians(pole_lat))
    ])
    if abs(pole[2]) < 0.9:
        u = np.cross(pole, [0, 0, 1])
    else:
        u = np.cross(pole, [1, 0, 0])
    u = u / np.linalg.norm(u)
    v = np.cross(pole, u)
    v = v / np.linalg.norm(v)
    n_pts = int(360.0 / step_deg)
    thetas = np.linspace(0, 2*np.pi, n_pts, endpoint=False)
    pts = np.outer(np.cos(thetas), u) + np.outer(np.sin(thetas), v)
    lats = np.degrees(np.arcsin(np.clip(pts[:, 2], -1, 1)))
    lons = np.degrees(np.arctan2(pts[:, 1], pts[:, 0]))
    return lats, lons

gc_lats, gc_lons = trace_great_circle(POLE_LAT, POLE_LON, step_deg=3.0)
gc_elevs = np.array([get_elev(la, lo) for la, lo in zip(gc_lats, gc_lons)])
gc_land = gc_elevs > 0

print(f"Traced Alison GC: {len(gc_lats)} pts, {gc_land.sum()} on land")

# Compute TPI and prominence for on-circle land points (sample up to 60 for speed)
gc_land_lats = gc_lats[gc_land]
gc_land_lons = gc_lons[gc_land]
if len(gc_land_lats) > 60:
    idx = np.random.RandomState(42).choice(len(gc_land_lats), 60, replace=False)
    gc_sample_lats = gc_land_lats[idx]
    gc_sample_lons = gc_land_lons[idx]
else:
    gc_sample_lats = gc_land_lats
    gc_sample_lons = gc_land_lons

print(f"Computing TPI/prominence for {len(gc_sample_lats)} on-circle land points...")
gc_tpi = []
gc_prom = []
gc_elev_on = []
gc_coast = []
gc_river = []
for la, lo in zip(gc_sample_lats, gc_sample_lons):
    gc_tpi.append(compute_tpi(la, lo, radius_km=10))
    gc_prom.append(compute_prominence(la, lo, radius_km=25))
    gc_elev_on.append(get_elev(la, lo))
    if HAS_GEO:
        gc_coast.append(min_dist_to_features(la, lo, coast_pts))
        gc_river.append(min_dist_to_features(la, lo, river_pts))

gc_tpi = np.array(gc_tpi)
gc_prom = np.array(gc_prom)
gc_elev_on = np.array(gc_elev_on)
gc_tpi_valid = gc_tpi[~np.isnan(gc_tpi)]
gc_prom_valid = gc_prom[~np.isnan(gc_prom)]

print(f"  On-circle TPI: μ={np.mean(gc_tpi_valid):.1f}m, median={np.median(gc_tpi_valid):.1f}m (n={len(gc_tpi_valid)})")
print(f"  On-circle Prominence: μ={np.mean(gc_prom_valid):.1f}m, median={np.median(gc_prom_valid):.1f}m (n={len(gc_prom_valid)})")
print(f"  On-circle Elevation: μ={np.mean(gc_elev_on):.0f}m, median={np.median(gc_elev_on):.0f}m")


# ============================================================
# STEP 3: Compare terrain signatures
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Terrain comparison — off-circle vs on-circle")
print("=" * 70)

off_tpi = np.array([r["tpi_10km_m"] for r in off_circle_results if r["tpi_10km_m"] is not None], dtype=float)
off_prom = np.array([r["prominence_25km_m"] for r in off_circle_results if r["prominence_25km_m"] is not None], dtype=float)
off_elev = np.array([r["elevation_m"] for r in off_circle_results], dtype=float)
off_coast = np.array([r["coastline_dist_km"] for r in off_circle_results if r["coastline_dist_km"] is not None], dtype=float)
off_river = np.array([r["river_dist_km"] for r in off_circle_results if r["river_dist_km"] is not None], dtype=float)

print(f"\nOFF-CIRCLE sites (n={len(off_circle_results)}):")
print(f"  Elevation: μ={np.mean(off_elev):.0f}m, median={np.median(off_elev):.0f}m")
print(f"  TPI: μ={np.mean(off_tpi):.1f}m, median={np.median(off_tpi):.1f}m (n={len(off_tpi)})")
print(f"  Prominence: μ={np.mean(off_prom):.1f}m, median={np.median(off_prom):.1f}m (n={len(off_prom)})")
if HAS_GEO and len(off_coast) > 0:
    print(f"  Coast dist: μ={np.mean(off_coast):.0f}km, median={np.median(off_coast):.0f}km")
    print(f"  River dist: μ={np.mean(off_river):.0f}km, median={np.median(off_river):.0f}km")

print(f"\nON-CIRCLE land points (n={len(gc_tpi_valid)}):")
print(f"  Elevation: μ={np.mean(gc_elev_on):.0f}m, median={np.median(gc_elev_on):.0f}m")
print(f"  TPI: μ={np.mean(gc_tpi_valid):.1f}m, median={np.median(gc_tpi_valid):.1f}m")
print(f"  Prominence: μ={np.mean(gc_prom_valid):.1f}m, median={np.median(gc_prom_valid):.1f}m")
if HAS_GEO and len(gc_coast) > 0:
    gc_coast_arr = np.array(gc_coast)
    gc_river_arr = np.array(gc_river)
    print(f"  Coast dist: μ={np.mean(gc_coast_arr):.0f}km, median={np.median(gc_coast_arr):.0f}km")
    print(f"  River dist: μ={np.mean(gc_river_arr):.0f}km, median={np.median(gc_river_arr):.0f}km")

# Statistical tests
print("\nStatistical tests (off-circle vs on-circle):")
if len(off_tpi) > 5 and len(gc_tpi_valid) > 5:
    stat_tpi, p_tpi = mannwhitneyu(off_tpi, gc_tpi_valid, alternative='two-sided')
    print(f"  TPI: Mann-Whitney U={stat_tpi:.0f}, p={p_tpi:.4f}")
else:
    p_tpi = np.nan

if len(off_prom) > 5 and len(gc_prom_valid) > 5:
    stat_prom, p_prom = mannwhitneyu(off_prom, gc_prom_valid, alternative='two-sided')
    print(f"  Prominence: Mann-Whitney U={stat_prom:.0f}, p={p_prom:.4f}")
else:
    p_prom = np.nan


# ============================================================
# STEP 4: Cluster membership analysis
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Which high-D clusters capture off-circle sites?")
print("=" * 70)

cluster_capture_50 = {cp["name"]: [] for cp in CLUSTER_POLES}
cluster_capture_100 = {cp["name"]: [] for cp in CLUSTER_POLES}
cluster_capture_200 = {cp["name"]: [] for cp in CLUSTER_POLES}

for r in off_circle_results:
    for cp in CLUSTER_POLES:
        cd = r["cluster_distances"][cp["name"]]
        if cd < 50:
            cluster_capture_50[cp["name"]].append(r["name"])
        if cd < 100:
            cluster_capture_100[cp["name"]].append(r["name"])
        if cd < 200:
            cluster_capture_200[cp["name"]].append(r["name"])

print("\nCluster capture at 100km threshold:")
print(f"{'Cluster':<30s} {'Sites captured':>15s} {'Names'}")
print("-" * 90)
for cp in CLUSTER_POLES:
    captured = cluster_capture_100[cp["name"]]
    if captured:
        print(f"{cp['name']:<30s} {len(captured):>15d}   {', '.join(captured[:5])}")
        if len(captured) > 5:
            print(f"{'':30s} {'':>15s}   ...and {len(captured)-5} more")

print("\nCluster capture at 200km threshold:")
for cp in CLUSTER_POLES:
    captured = cluster_capture_200[cp["name"]]
    if captured:
        print(f"  {cp['name']}: {len(captured)} sites — {', '.join(captured[:8])}")


# ============================================================
# STEP 5: Göbekli Tepe special analysis
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: Göbekli Tepe — special focus")
print("=" * 70)

gt = next(r for r in off_circle_results if r["name"] == "Göbekli Tepe")
print(f"  Alison GC distance: {gt['alison_gc_dist_km']:.1f} km")
print(f"  Elevation: {gt['elevation_m']:.0f} m")
print(f"  TPI: {gt['tpi_10km_m']} m")
print(f"  Prominence: {gt['prominence_25km_m']} m")
print(f"  Nearest cluster: {gt['nearest_cluster']} ({gt['nearest_cluster_dist_km']:.1f} km)")
print(f"  On clusters (<100km): {gt['on_clusters_100km']}")
print("\n  Distances to all 12 clusters:")
for name, d in sorted(gt["cluster_distances"].items(), key=lambda x: x[1]):
    marker = " *** ON CIRCLE" if d < 100 else ""
    print(f"    {name:<30s}: {d:>8.1f} km{marker}")


# ============================================================
# STEP 6: Search for second alignment through off-circle sites
# ============================================================
print("\n" + "=" * 70)
print("STEP 6: Searching for 'second alignment' through off-circle sites")
print("=" * 70)

# Brute-force: for each pair of off-circle sites, compute the GC pole,
# then count how many other off-circle sites fall within 100km of that GC.
# This tests whether any subset of famous sites defines a GC.

print("Testing all pairs of off-circle sites for GC alignment...")
n_sites = len(OFF_CIRCLE_SITES)
best_pole = None
best_count = 0
best_sites_on = []

for i in range(n_sites):
    for j in range(i+1, n_sites):
        lat1, lon1 = OFF_CIRCLE_SITES[i]["lat"], OFF_CIRCLE_SITES[i]["lon"]
        lat2, lon2 = OFF_CIRCLE_SITES[j]["lat"], OFF_CIRCLE_SITES[j]["lon"]

        # Compute GC pole from two points
        lat1r, lon1r = math.radians(lat1), math.radians(lon1)
        lat2r, lon2r = math.radians(lat2), math.radians(lon2)

        p1 = np.array([math.cos(lat1r)*math.cos(lon1r),
                        math.cos(lat1r)*math.sin(lon1r),
                        math.sin(lat1r)])
        p2 = np.array([math.cos(lat2r)*math.cos(lon2r),
                        math.cos(lat2r)*math.sin(lon2r),
                        math.sin(lat2r)])

        pole = np.cross(p1, p2)
        norm = np.linalg.norm(pole)
        if norm < 1e-10:
            continue
        pole = pole / norm

        pole_lat = math.degrees(math.asin(max(-1, min(1, pole[2]))))
        pole_lon = math.degrees(math.atan2(pole[1], pole[0]))

        # Count how many sites are within 100km of this GC
        count = 0
        sites_on = []
        for k, site in enumerate(OFF_CIRCLE_SITES):
            d = dist_from_gc(site["lat"], site["lon"], pole_lat, pole_lon)
            if d < 100:
                count += 1
                sites_on.append(site["name"])

        if count > best_count:
            best_count = count
            best_pole = (pole_lat, pole_lon)
            best_sites_on = sites_on
            best_pair = (OFF_CIRCLE_SITES[i]["name"], OFF_CIRCLE_SITES[j]["name"])

print(f"\nBest alignment found:")
print(f"  Pole: ({best_pole[0]:.2f}, {best_pole[1]:.2f})")
print(f"  Defining pair: {best_pair[0]} — {best_pair[1]}")
print(f"  Sites within 100km: {best_count}/50")
print(f"  Sites: {', '.join(best_sites_on)}")

# Monte Carlo: is this count significant?
print(f"\nMonte Carlo: how many sites fall on a random GC through 2 random sites?")
rng = np.random.RandomState(42)
random_counts = []
n_mc = 5000
for trial in range(n_mc):
    # Pick random pole
    z = rng.uniform(-1, 1)
    theta = rng.uniform(0, 2*np.pi)
    plat = math.degrees(math.asin(z))
    plon = math.degrees(theta) - 180

    count = 0
    for site in OFF_CIRCLE_SITES:
        d = dist_from_gc(site["lat"], site["lon"], plat, plon)
        if d < 100:
            count += 1
    random_counts.append(count)

random_counts = np.array(random_counts)
pctile = percentileofscore(random_counts, best_count)
print(f"  Random GC mean sites within 100km: {np.mean(random_counts):.1f} ± {np.std(random_counts):.1f}")
print(f"  Best alignment ({best_count}) is at percentile {pctile:.1f}%")
print(f"  p-value (≥{best_count}): {np.sum(random_counts >= best_count) / n_mc:.4f}")


# ============================================================
# STEP 7: Generate outputs
# ============================================================
print("\n" + "=" * 70)
print("STEP 7: Generating outputs")
print("=" * 70)

# ── Table: which clusters each site falls on ──
print("\n--- CLUSTER MEMBERSHIP TABLE ---")
print(f"{'Site':<25s} {'Alison km':>10s} {'Nearest Cluster':<30s} {'Dist km':>8s} {'On clusters (<100km)'}")
print("-" * 110)
for r in sorted(off_circle_results, key=lambda x: x["alison_gc_dist_km"]):
    clusters_str = ", ".join(r["on_clusters_100km"]) if r["on_clusters_100km"] else "—"
    print(f"{r['name']:<25s} {r['alison_gc_dist_km']:>10.1f} {r['nearest_cluster']:<30s} {r['nearest_cluster_dist_km']:>8.1f} {clusters_str}")


# ── Plot: terrain comparison ──
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# TPI comparison
ax = axes[0]
data_tpi = [off_tpi[~np.isnan(off_tpi)], gc_tpi_valid]
labels_tpi = [f'Off-circle\n(n={len(data_tpi[0])})', f'On-circle\n(n={len(data_tpi[1])})']
bp = ax.boxplot(data_tpi, labels=labels_tpi, patch_artist=True)
bp['boxes'][0].set_facecolor('#e74c3c')
bp['boxes'][1].set_facecolor('#3498db')
ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax.set_ylabel('TPI (m)')
ax.set_title(f'TPI (10km radius)\np={p_tpi:.4f}' if not np.isnan(p_tpi) else 'TPI (10km radius)')

# Prominence comparison
ax = axes[1]
data_prom = [off_prom[~np.isnan(off_prom)], gc_prom_valid]
labels_prom = [f'Off-circle\n(n={len(data_prom[0])})', f'On-circle\n(n={len(data_prom[1])})']
bp = ax.boxplot(data_prom, labels=labels_prom, patch_artist=True)
bp['boxes'][0].set_facecolor('#e74c3c')
bp['boxes'][1].set_facecolor('#3498db')
ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax.set_ylabel('Prominence (m)')
ax.set_title(f'Prominence (25km radius)\np={p_prom:.4f}' if not np.isnan(p_prom) else 'Prominence (25km radius)')

# Elevation comparison
ax = axes[2]
off_elev_land = off_elev[off_elev > 0]
gc_elev_land = gc_elev_on[gc_elev_on > 0]
data_elev = [off_elev_land, gc_elev_land]
labels_elev = [f'Off-circle\n(n={len(off_elev_land)})', f'On-circle\n(n={len(gc_elev_land)})']
bp = ax.boxplot(data_elev, labels=labels_elev, patch_artist=True)
bp['boxes'][0].set_facecolor('#e74c3c')
bp['boxes'][1].set_facecolor('#3498db')
ax.set_ylabel('Elevation (m)')
ax.set_title('Elevation')

fig.suptitle('Directive 2: Off-Circle vs On-Circle Terrain Comparison', fontsize=13, y=1.02)
fig.tight_layout()
fig.savefig(OUT / 'off_circle_terrain_comparison.png', dpi=150)
plt.close()
print(f"  Saved: {OUT / 'off_circle_terrain_comparison.png'}")


# ── Save JSON ──
output_json = {
    "meta": {
        "name": "Directive 2: Off-Circle Major Sites Analysis",
        "date": "2026-03-21",
        "n_off_circle_sites": len(off_circle_results),
        "alison_pole": {"lat": POLE_LAT, "lon": POLE_LON},
        "n_cluster_poles": len(CLUSTER_POLES),
    },
    "off_circle_sites": off_circle_results,
    "terrain_comparison": {
        "off_circle": {
            "elevation_mean_m": round(float(np.mean(off_elev)), 1),
            "elevation_median_m": round(float(np.median(off_elev)), 1),
            "tpi_mean_m": round(float(np.mean(off_tpi)), 1),
            "tpi_median_m": round(float(np.median(off_tpi)), 1),
            "prominence_mean_m": round(float(np.mean(off_prom)), 1) if len(off_prom) > 0 else None,
            "prominence_median_m": round(float(np.median(off_prom)), 1) if len(off_prom) > 0 else None,
        },
        "on_circle": {
            "elevation_mean_m": round(float(np.mean(gc_elev_on)), 1),
            "elevation_median_m": round(float(np.median(gc_elev_on)), 1),
            "tpi_mean_m": round(float(np.mean(gc_tpi_valid)), 1),
            "tpi_median_m": round(float(np.median(gc_tpi_valid)), 1),
            "prominence_mean_m": round(float(np.mean(gc_prom_valid)), 1) if len(gc_prom_valid) > 0 else None,
            "prominence_median_m": round(float(np.median(gc_prom_valid)), 1) if len(gc_prom_valid) > 0 else None,
        },
        "mann_whitney_tpi_p": round(float(p_tpi), 6) if not np.isnan(p_tpi) else None,
        "mann_whitney_prominence_p": round(float(p_prom), 6) if not np.isnan(p_prom) else None,
    },
    "cluster_capture": {
        "threshold_50km": {cp["name"]: cluster_capture_50[cp["name"]] for cp in CLUSTER_POLES if cluster_capture_50[cp["name"]]},
        "threshold_100km": {cp["name"]: cluster_capture_100[cp["name"]] for cp in CLUSTER_POLES if cluster_capture_100[cp["name"]]},
        "threshold_200km": {cp["name"]: cluster_capture_200[cp["name"]] for cp in CLUSTER_POLES if cluster_capture_200[cp["name"]]},
    },
    "gobekli_tepe_focus": {
        "alison_dist_km": gt["alison_gc_dist_km"],
        "nearest_cluster": gt["nearest_cluster"],
        "nearest_cluster_dist_km": gt["nearest_cluster_dist_km"],
        "on_any_cluster_100km": len(gt["on_clusters_100km"]) > 0,
        "clusters_within_100km": gt["on_clusters_100km"],
    },
    "second_alignment_search": {
        "best_pole": {"lat": round(best_pole[0], 2), "lon": round(best_pole[1], 2)},
        "best_pair": list(best_pair),
        "best_count_within_100km": best_count,
        "best_sites": best_sites_on,
        "mc_mean_random_count": round(float(np.mean(random_counts)), 2),
        "mc_std_random_count": round(float(np.std(random_counts)), 2),
        "percentile": round(pctile, 1),
        "p_value": round(float(np.sum(random_counts >= best_count) / n_mc), 4),
    },
}

out_file = OUT / "off_circle_analysis.json"
with open(out_file, "w") as f:
    json.dump(output_json, f, indent=2)
print(f"  Saved: {out_file}")


# ── Final summary ──
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"\n1. OFF-CIRCLE DISTANCES TO ALISON GC:")
dists = [r["alison_gc_dist_km"] for r in off_circle_results]
print(f"   Range: {min(dists):.0f} — {max(dists):.0f} km")
print(f"   Mean: {np.mean(dists):.0f} km, Median: {np.median(dists):.0f} km")

print(f"\n2. TERRAIN SIGNATURES:")
print(f"   Off-circle monuments have {'HIGHER' if np.mean(off_tpi) > np.mean(gc_tpi_valid) else 'LOWER'} TPI "
      f"({np.mean(off_tpi):.0f}m vs {np.mean(gc_tpi_valid):.0f}m)")
print(f"   Off-circle monuments have {'HIGHER' if np.mean(off_elev) > np.mean(gc_elev_on) else 'LOWER'} elevation "
      f"({np.mean(off_elev):.0f}m vs {np.mean(gc_elev_on):.0f}m)")

n_any_cluster = sum(1 for r in off_circle_results if len(r["on_clusters_100km"]) > 0)
print(f"\n3. CLUSTER MEMBERSHIP:")
print(f"   {n_any_cluster}/50 off-circle sites fall on at least one high-D cluster circle (100km)")

print(f"\n4. GÖBEKLI TEPE:")
print(f"   {gt['alison_gc_dist_km']:.0f} km from Alison GC")
print(f"   Nearest cluster: {gt['nearest_cluster']} ({gt['nearest_cluster_dist_km']:.0f} km)")
print(f"   On any cluster (<100km): {'YES' if gt['on_clusters_100km'] else 'NO'}")

print(f"\n5. SECOND ALIGNMENT:")
print(f"   Best GC through off-circle sites captures {best_count}/50 within 100km")
print(f"   Percentile vs random: {pctile:.1f}%")
p_val_2nd = np.sum(random_counts >= best_count) / n_mc
print(f"   {'SIGNIFICANT' if p_val_2nd < 0.05 else 'NOT significant'} (p={p_val_2nd:.4f})")

print(f"\nAll outputs saved to: {OUT}")
print("Done.")
