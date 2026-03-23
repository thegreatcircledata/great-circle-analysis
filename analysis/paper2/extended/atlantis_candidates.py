#!/usr/bin/env python3
"""
Directive 5: Alison Atlantis Candidate Locations
=================================================
Tests whether Alison's three Atlantis candidate locations are:
1. On the Great Circle (distance to GC)
2. Above water during the LGM (ETOPO1 bathymetry, sea level -120m)
3. Near any radiocarbon dates in the p3k14c database

Also finds where the GC crosses the Bay of Bengal and South China Sea.
"""

import json, math, os, sys
import numpy as np
import netCDF4 as nc
import csv

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

# ============================================================
# CONFIGURATION
# ============================================================
POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2  # ~10,007.5 km
LGM_SEA_LEVEL = -120  # meters (conservative estimate)
SEARCH_RADIUS_KM = 100  # for radiocarbon date search

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "extended_analysis")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))

def dist_from_gc(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QUARTER_CIRC)

def trace_great_circle(pole_lat, pole_lon, step_deg=0.1):
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
    n_points = int(360.0 / step_deg)
    thetas = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
    points = np.outer(np.cos(thetas), u) + np.outer(np.sin(thetas), v)
    lats = np.degrees(np.arcsin(np.clip(points[:, 2], -1, 1)))
    lons = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    return lats, lons

# ============================================================
# LOAD ETOPO1
# ============================================================
print("=" * 70)
print("LOADING ETOPO1 BATHYMETRY")
print("=" * 70)

etopo_path = os.path.join(BASE_DIR, "data", "geophysical", "etopo", "ETOPO1_Ice_g_gmt4.grd")
ds = nc.Dataset(etopo_path)
etopo_lons = ds.variables['x'][:]
etopo_lats = ds.variables['y'][:]
etopo_z = ds.variables['z']
etopo_lon_min, etopo_lon_step = float(etopo_lons[0]), float(etopo_lons[1] - etopo_lons[0])
etopo_lat_min, etopo_lat_step = float(etopo_lats[0]), float(etopo_lats[1] - etopo_lats[0])
print(f"ETOPO1 loaded: {len(etopo_lats)} x {len(etopo_lons)}")

def get_depth(lat, lon):
    lat_idx = int(round((lat - etopo_lat_min) / etopo_lat_step))
    lon_idx = int(round((lon - etopo_lon_min) / etopo_lon_step))
    lat_idx = max(0, min(lat_idx, len(etopo_lats) - 1))
    lon_idx = max(0, min(lon_idx, len(etopo_lons) - 1))
    return float(etopo_z[lat_idx, lon_idx])

def get_depth_profile(center_lat, center_lon, radius_deg=2.0, step_deg=0.05):
    """Get bathymetric profile in a grid around a point."""
    lats = np.arange(center_lat - radius_deg, center_lat + radius_deg, step_deg)
    lons = np.arange(center_lon - radius_deg, center_lon + radius_deg, step_deg)
    profile = {'center_lat': center_lat, 'center_lon': center_lon,
               'min_depth': float('inf'), 'max_depth': float('-inf'),
               'shallowest_point': None}
    for lat in lats:
        for lon in lons:
            d = get_depth(lat, lon)
            if d < profile['min_depth']:
                profile['min_depth'] = d
            if d > profile['max_depth']:
                profile['max_depth'] = d
                profile['shallowest_point'] = {'lat': float(lat), 'lon': float(lon), 'depth_m': d}
    return profile

# ============================================================
# FIND GC CROSSINGS: Bay of Bengal & South China Sea
# ============================================================
print(f"\n{'='*70}")
print("FINDING GREAT CIRCLE CROSSINGS")
print(f"{'='*70}")

gc_lats, gc_lons = trace_great_circle(POLE_LAT, POLE_LON, step_deg=0.05)

# Bay of Bengal bounding box: ~5-22°N, 78-95°E
bob_mask = (gc_lats >= 5) & (gc_lats <= 22) & (gc_lons >= 78) & (gc_lons <= 95)
bob_indices = np.where(bob_mask)[0]

# South China Sea bounding box: ~0-25°N, 105-121°E
scs_mask = (gc_lats >= 0) & (gc_lats <= 25) & (gc_lons >= 105) & (gc_lons <= 121)
scs_indices = np.where(scs_mask)[0]

bob_crossing = None
scs_crossing = None

if len(bob_indices) > 0:
    mid = bob_indices[len(bob_indices)//2]
    bob_crossing = {'lat': float(gc_lats[mid]), 'lon': float(gc_lons[mid])}
    print(f"Bay of Bengal crossing: {bob_crossing['lat']:.3f}°N, {bob_crossing['lon']:.3f}°E")
    # Also get entry/exit
    print(f"  Entry: {gc_lats[bob_indices[0]]:.3f}°N, {gc_lons[bob_indices[0]]:.3f}°E")
    print(f"  Exit:  {gc_lats[bob_indices[-1]]:.3f}°N, {gc_lons[bob_indices[-1]]:.3f}°E")
else:
    print("WARNING: Great Circle does NOT cross Bay of Bengal")

if len(scs_indices) > 0:
    mid = scs_indices[len(scs_indices)//2]
    scs_crossing = {'lat': float(gc_lats[mid]), 'lon': float(gc_lons[mid])}
    print(f"South China Sea crossing: {scs_crossing['lat']:.3f}°N, {scs_crossing['lon']:.3f}°E")
    print(f"  Entry: {gc_lats[scs_indices[0]]:.3f}°N, {gc_lons[scs_indices[0]]:.3f}°E")
    print(f"  Exit:  {gc_lats[scs_indices[-1]]:.3f}°N, {gc_lons[scs_indices[-1]]:.3f}°E")
else:
    print("WARNING: Great Circle does NOT cross South China Sea")

# ============================================================
# LOAD RADIOCARBON DATABASE
# ============================================================
print(f"\n{'='*70}")
print("LOADING P3K14C RADIOCARBON DATABASE")
print(f"{'='*70}")

p3k_path = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")
c14_sites = []
with open(p3k_path, 'r', encoding='utf-8', errors='replace') as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            lat = float(row['Lat'])
            lon = float(row['Long'])
            age = int(row['Age']) if row['Age'] else None
            c14_sites.append({'lat': lat, 'lon': lon, 'age': age,
                            'site': row.get('SiteName', ''),
                            'country': row.get('Country', '')})
        except (ValueError, KeyError):
            continue

print(f"Loaded {len(c14_sites)} radiocarbon dates")

def count_c14_nearby(lat, lon, radius_km=100):
    """Count radiocarbon dates within radius_km of a point."""
    nearby = []
    for s in c14_sites:
        d = haversine_km(lat, lon, s['lat'], s['lon'])
        if d <= radius_km:
            nearby.append({**s, 'distance_km': round(d, 1)})
    # Sort by distance
    nearby.sort(key=lambda x: x['distance_km'])
    return nearby

# ============================================================
# DEFINE ALL CANDIDATE LOCATIONS
# ============================================================
print(f"\n{'='*70}")
print("EVALUATING ATLANTIS CANDIDATE LOCATIONS")
print(f"{'='*70}")

candidates = [
    {"name": "Cape Verde (south of)", "lat": 14.0, "lon": -24.0,
     "source": "Alison — approximate location south of Cape Verde islands"},
    {"name": "Mid-Atlantic (Piri Reis)", "lat": 4.317, "lon": -41.5,
     "source": "Alison — specific coordinates 4°19'N, 41°30'W from Piri Reis map analysis"},
]

# Add BoB crossing if found
if bob_crossing:
    candidates.append({
        "name": "Bay of Bengal crossing",
        "lat": bob_crossing['lat'], "lon": bob_crossing['lon'],
        "source": "Computed — where Alison GC crosses Bay of Bengal"
    })

# Add SCS crossing if found
if scs_crossing:
    candidates.append({
        "name": "South China Sea crossing",
        "lat": scs_crossing['lat'], "lon": scs_crossing['lon'],
        "source": "Computed — where Alison GC crosses South China Sea"
    })

# ============================================================
# EVALUATE EACH CANDIDATE
# ============================================================
results = []

for c in candidates:
    print(f"\n--- {c['name']} ---")
    lat, lon = c['lat'], c['lon']

    # 1. Distance to Great Circle
    gc_dist = dist_from_gc(lat, lon)
    on_circle = gc_dist < 50  # within 50km threshold
    print(f"  Coordinates: {lat:.3f}°{'N' if lat >= 0 else 'S'}, "
          f"{abs(lon):.3f}°{'W' if lon < 0 else 'E'}")
    print(f"  Distance to GC: {gc_dist:.1f} km {'✓ ON CIRCLE' if on_circle else '✗ OFF CIRCLE'}")

    # 2. Bathymetric depth
    depth = get_depth(lat, lon)
    above_lgm = depth > LGM_SEA_LEVEL  # was above water when sea level was 120m lower
    print(f"  Depth: {depth:.0f} m {'(above water)' if depth >= 0 else ''}")
    print(f"  Above water at LGM (-{abs(LGM_SEA_LEVEL)}m): {'YES' if above_lgm else 'NO'}")

    # 3. Local bathymetric profile (find shallowest point nearby)
    print(f"  Scanning local bathymetry (±2° grid)...")
    local_profile = get_depth_profile(lat, lon, radius_deg=2.0, step_deg=0.1)
    shallowest = local_profile['shallowest_point']
    print(f"  Shallowest nearby: {shallowest['depth_m']:.0f}m at "
          f"{shallowest['lat']:.2f}°, {shallowest['lon']:.2f}°")

    # 4. Radiocarbon dates nearby
    nearby_c14 = count_c14_nearby(lat, lon, SEARCH_RADIUS_KM)
    print(f"  C14 dates within {SEARCH_RADIUS_KM}km: {len(nearby_c14)}")
    if nearby_c14:
        # Show oldest few
        dated = [s for s in nearby_c14 if s['age']]
        if dated:
            dated.sort(key=lambda x: -x['age'])
            for s in dated[:3]:
                print(f"    {s['site'] or 'Unknown'} ({s['country']}): "
                      f"{s['age']} BP, {s['distance_km']} km away")

    # 5. Verdict
    if on_circle and above_lgm:
        verdict = "VIABLE — on circle and above water during LGM"
    elif on_circle and not above_lgm and depth > -500:
        verdict = "MARGINAL — on circle but submerged (shallow enough for shelf)"
    elif on_circle and not above_lgm:
        verdict = f"REJECTED — on circle but {depth:.0f}m deep (abyss)"
    elif not on_circle and above_lgm:
        verdict = f"REJECTED — {gc_dist:.0f}km off circle"
    else:
        verdict = f"REJECTED — {gc_dist:.0f}km off circle AND {depth:.0f}m deep"
    print(f"  VERDICT: {verdict}")

    result = {
        'name': c['name'],
        'lat': lat,
        'lon': lon,
        'source': c['source'],
        'gc_distance_km': round(gc_dist, 1),
        'on_circle_50km': on_circle,
        'depth_m': round(depth, 1),
        'above_lgm': above_lgm,
        'lgm_sea_level_used_m': LGM_SEA_LEVEL,
        'shallowest_nearby_m': round(shallowest['depth_m'], 1),
        'shallowest_nearby_coords': {
            'lat': round(shallowest['lat'], 3),
            'lon': round(shallowest['lon'], 3)
        },
        'c14_dates_within_100km': len(nearby_c14),
        'oldest_c14_dates': [
            {'site': s.get('site', ''), 'country': s.get('country', ''),
             'age_bp': s['age'], 'distance_km': s['distance_km']}
            for s in sorted(
                [x for x in nearby_c14 if x['age']],
                key=lambda x: -x['age']
            )[:5]
        ],
        'verdict': verdict
    }
    results.append(result)

# ============================================================
# DETAILED BATHYMETRIC TRANSECT: Mid-Atlantic (Piri Reis)
# ============================================================
print(f"\n{'='*70}")
print("DETAILED ANALYSIS: MID-ATLANTIC (PIRI REIS) SITE")
print(f"{'='*70}")

piri_lat, piri_lon = 4.317, -41.5
print(f"Piri Reis claimed island: {piri_lat}°N, {abs(piri_lon)}°W")
print(f"\nE-W transect at {piri_lat}°N from {piri_lon-3}°W to {piri_lon+3}°W:")
print(f"{'Longitude':>10} {'Depth (m)':>10}")
print("-" * 22)
ew_transect = []
for dlon in np.arange(-3, 3.1, 0.25):
    lon_t = piri_lon + dlon
    d = get_depth(piri_lat, lon_t)
    ew_transect.append({'lon': round(float(lon_t), 2), 'depth_m': round(d, 1)})
    marker = " <<<" if abs(dlon) < 0.01 else ""
    print(f"  {lon_t:>8.2f}° {d:>8.0f}m{marker}")

# Check for any seamounts or ridges (depth > -2000m in the area)
print(f"\nSearching for seamounts/ridges within ±5° of Piri Reis location...")
seamounts = []
for lat_s in np.arange(piri_lat - 5, piri_lat + 5, 0.25):
    for lon_s in np.arange(piri_lon - 5, piri_lon + 5, 0.25):
        d = get_depth(lat_s, lon_s)
        if d > -1000 and d < 0:  # shallow submarine feature
            seamounts.append({'lat': round(float(lat_s), 2),
                            'lon': round(float(lon_s), 2),
                            'depth_m': round(d, 1)})

if seamounts:
    # Find shallowest
    seamounts.sort(key=lambda x: -x['depth_m'])
    print(f"Found {len(seamounts)} shallow points (above -1000m):")
    for s in seamounts[:10]:
        gc_d = dist_from_gc(s['lat'], s['lon'])
        print(f"  {s['lat']:.2f}°N, {abs(s['lon']):.2f}°W: {s['depth_m']:.0f}m "
              f"(GC dist: {gc_d:.0f}km)")
else:
    print("No seamounts found — area is deep abyssal plain")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

for r in results:
    status = "✓" if "VIABLE" in r['verdict'] else "✗"
    print(f"{status} {r['name']}")
    print(f"    GC: {r['gc_distance_km']}km | Depth: {r['depth_m']}m | "
          f"LGM: {'dry' if r['above_lgm'] else 'submerged'} | "
          f"C14: {r['c14_dates_within_100km']} dates")
    print(f"    → {r['verdict']}")

# ============================================================
# SAVE OUTPUT
# ============================================================
output = {
    'analysis': 'Directive 5 — Alison Atlantis Candidate Locations',
    'parameters': {
        'pole_lat': POLE_LAT,
        'pole_lon': POLE_LON,
        'lgm_sea_level_m': LGM_SEA_LEVEL,
        'gc_threshold_km': 50,
        'c14_search_radius_km': SEARCH_RADIUS_KM,
    },
    'gc_crossings': {
        'bay_of_bengal': bob_crossing,
        'south_china_sea': scs_crossing,
        'bay_of_bengal_all_points': [
            {'lat': round(float(gc_lats[i]), 3), 'lon': round(float(gc_lons[i]), 3)}
            for i in bob_indices
        ] if len(bob_indices) > 0 else [],
        'south_china_sea_all_points': [
            {'lat': round(float(gc_lats[i]), 3), 'lon': round(float(gc_lons[i]), 3)}
            for i in scs_indices
        ] if len(scs_indices) > 0 else [],
    },
    'candidates': results,
    'piri_reis_transect': ew_transect,
    'nearby_seamounts': seamounts[:20] if seamounts else [],
}

out_path = os.path.join(OUT_DIR, "atlantis_candidates.json")
with open(out_path, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved: {out_path}")

ds.close()
