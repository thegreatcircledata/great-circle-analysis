#!/usr/bin/env python3
"""
Paper 3 Robustness Battery — Task 5: River-Perpendicular Control
=================================================================
Identify all points where the Great Circle crosses a major river at a steep angle (>60°).
For each crossing, compute monument density within 25 km.
Compare to the Memphis crossing. Is Memphis unique?
"""

import json, os, sys, time
import numpy as np
import fiona
import shapely.geometry as sg

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE = os.path.expanduser("~/megalith_site_research")
OUT = os.path.join(BASE, "paper3/results")
FIG = os.path.join(BASE, "paper3/figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

POLE_LAT = 59.682122
POLE_LON = -138.646087
R = 6371.0
QC = R * np.pi / 2

def gc_distance_to_circle(lats, lons):
    la1 = np.radians(POLE_LAT)
    lo1 = np.radians(POLE_LON)
    la2 = np.radians(lats)
    lo2 = np.radians(lons)
    dl = la2 - la1
    dn = lo2 - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(la2)*np.sin(dn/2)**2
    d = R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QC)

def gc_bearing_at(lat, lon):
    """Compute GC bearing at a point (perpendicular to radial from pole)."""
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    pole_lat_r, pole_lon_r = np.radians(POLE_LAT), np.radians(POLE_LON)
    dlon = pole_lon_r - lon_r
    x = np.sin(dlon) * np.cos(pole_lat_r)
    y = np.cos(lat_r) * np.sin(pole_lat_r) - np.sin(lat_r) * np.cos(pole_lat_r) * np.cos(dlon)
    bearing_to_pole = np.degrees(np.arctan2(x, y)) % 360
    return (bearing_to_pole + 90) % 360  # GC is perpendicular to radial

# ============================================================
# LOAD RIVERS
# ============================================================
print("Loading rivers (Natural Earth 50m)...")
river_path = os.path.join(BASE, "data/natural_earth/ne_50m_rivers_lake_centerlines.shp")

rivers = []
with fiona.open(river_path) as src:
    for feat in src:
        if feat['geometry'] is None:
            continue
        props = dict(feat['properties'])
        name = props.get('name', 'Unknown')
        geom = sg.shape(feat['geometry'])

        if geom.geom_type == 'LineString':
            coords = list(geom.coords)
            rivers.append({'name': name, 'coords': coords})
        elif geom.geom_type == 'MultiLineString':
            for line in geom.geoms:
                coords = list(line.coords)
                rivers.append({'name': name, 'coords': coords})

print(f"Loaded {len(rivers)} river segments")

# ============================================================
# FIND RIVER-GC CROSSINGS
# ============================================================
print("\n" + "=" * 70)
print("FINDING RIVER-GC CROSSINGS (within 5 km of GC)")
print("=" * 70)

crossings = []

for river in rivers:
    coords = river['coords']
    name = river['name']
    if not name or name == 'None':
        name = 'unnamed'

    for i in range(len(coords) - 1):
        lon1, lat1 = coords[i]
        lon2, lat2 = coords[i + 1]

        # Check if segment is near the GC
        mid_lat = (lat1 + lat2) / 2
        mid_lon = (lon1 + lon2) / 2
        d = gc_distance_to_circle(np.array([mid_lat]), np.array([mid_lon]))[0]

        if d > 10:  # skip segments far from GC
            continue

        # Compute river segment bearing
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        cos_lat = np.cos(np.radians(mid_lat))
        river_bearing = np.degrees(np.arctan2(dlon * cos_lat, dlat)) % 360

        # Compute GC bearing at midpoint
        gc_bear = gc_bearing_at(mid_lat, mid_lon)

        # Intersection angle (acute)
        angle = abs(river_bearing - gc_bear) % 180
        if angle > 90:
            angle = 180 - angle

        crossings.append({
            'name': name,
            'lat': mid_lat,
            'lon': mid_lon,
            'dist_from_gc_km': float(d),
            'river_bearing': float(river_bearing),
            'gc_bearing': float(gc_bear),
            'intersection_angle': float(angle)
        })

print(f"Found {len(crossings)} river segments within 10 km of GC")

# Filter to steep crossings (>60°)
steep = [c for c in crossings if c['intersection_angle'] > 60]
print(f"Steep crossings (>60°): {len(steep)}")

# Deduplicate by grouping nearby crossings of the same river
def deduplicate(crossings_list, min_dist_km=50):
    """Group crossings within min_dist_km of each other on the same river."""
    if not crossings_list:
        return []

    sorted_c = sorted(crossings_list, key=lambda c: (c['name'], c['lat']))
    groups = []
    current_group = [sorted_c[0]]

    for c in sorted_c[1:]:
        last = current_group[-1]
        # Same river and within min_dist_km?
        dlat = (c['lat'] - last['lat']) * 111.32
        dlon = (c['lon'] - last['lon']) * 111.32 * np.cos(np.radians(c['lat']))
        d = np.sqrt(dlat**2 + dlon**2)

        if c['name'] == last['name'] and d < min_dist_km:
            current_group.append(c)
        else:
            groups.append(current_group)
            current_group = [c]
    groups.append(current_group)

    # Take the crossing with minimum GC distance from each group
    deduped = []
    for group in groups:
        best = min(group, key=lambda c: c['dist_from_gc_km'])
        best['n_segments'] = len(group)
        deduped.append(best)

    return deduped

steep_deduped = deduplicate(steep)
print(f"Unique steep river crossings: {len(steep_deduped)}")

# ============================================================
# LOAD PLEIADES MONUMENTS FOR DENSITY COMPUTATION
# ============================================================
print("\nLoading Pleiades monuments...")
import csv

MONUMENTAL = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus",
    "shrine", "cemetery"
}

mon_lats, mon_lons = [], []
with open(os.path.join(BASE, "data/pleiades/pleiades-places-latest.csv")) as f:
    for row in csv.DictReader(f):
        try:
            lat = float(row['reprLat'])
            lon = float(row['reprLong'])
        except (ValueError, TypeError):
            continue
        if abs(lat) < 0.001 and abs(lon) < 0.001:
            continue
        fts = set(ft.strip() for ft in row['featureTypes'].split(','))
        if fts & MONUMENTAL:
            mon_lats.append(lat)
            mon_lons.append(lon)

mon_lats = np.array(mon_lats)
mon_lons = np.array(mon_lons)
print(f"Monuments loaded: {len(mon_lats)}")

# ============================================================
# COMPUTE MONUMENT DENSITY AT EACH CROSSING
# ============================================================
print("\n" + "=" * 70)
print("MONUMENT DENSITY AT STEEP RIVER CROSSINGS (25 km radius)")
print("=" * 70)

RADIUS = 25.0  # km

for crossing in steep_deduped:
    clat, clon = crossing['lat'], crossing['lon']
    dlat = (mon_lats - clat) * 111.32
    cos_lat = np.cos(np.radians(clat))
    dlon = (mon_lons - clon) * 111.32 * cos_lat
    d = np.sqrt(dlat**2 + dlon**2)
    n_in_radius = int(np.sum(d <= RADIUS))
    area_km2 = np.pi * RADIUS**2
    density = n_in_radius / area_km2
    crossing['monuments_within_25km'] = n_in_radius
    crossing['monument_density_per_km2'] = float(density)

# Sort by monument count
steep_deduped.sort(key=lambda c: c['monuments_within_25km'], reverse=True)

# Memphis reference
memphis_lat, memphis_lon = 29.87, 31.22
memphis_gc_dist = gc_distance_to_circle(np.array([memphis_lat]), np.array([memphis_lon]))[0]
dlat = (mon_lats - memphis_lat) * 111.32
cos_lat = np.cos(np.radians(memphis_lat))
dlon = (mon_lons - memphis_lon) * 111.32 * cos_lat
d = np.sqrt(dlat**2 + dlon**2)
memphis_n = int(np.sum(d <= RADIUS))
memphis_density = memphis_n / (np.pi * RADIUS**2)

print(f"\nMemphis reference: {memphis_n} monuments within {RADIUS}km "
      f"(density: {memphis_density:.4f}/km²)")
print(f"Memphis GC distance: {memphis_gc_dist:.1f} km")

print(f"\n{'Name':<30} {'Lat':>7} {'Lon':>8} {'Angle':>6} {'Mon':>5} {'Density':>10}")
print("-" * 70)
for c in steep_deduped[:20]:
    print(f"{c['name'][:30]:<30} {c['lat']:>7.2f} {c['lon']:>8.2f} "
          f"{c['intersection_angle']:>5.1f}° {c['monuments_within_25km']:>5d} "
          f"{c['monument_density_per_km2']:>.4f}")

# Compare to Memphis
print(f"\n{'Memphis (reference)':<30} {memphis_lat:>7.2f} {memphis_lon:>8.2f} "
      f"{'81.3':>6}° {memphis_n:>5d} {memphis_density:>.4f}")

# Is Memphis unique?
higher = [c for c in steep_deduped if c['monuments_within_25km'] > memphis_n]
print(f"\nCrossings with MORE monuments than Memphis: {len(higher)}")
if higher:
    for h in higher:
        print(f"  {h['name']} ({h['lat']:.2f}, {h['lon']:.2f}): "
              f"{h['monuments_within_25km']} monuments")

comparable = [c for c in steep_deduped
              if c['monuments_within_25km'] >= memphis_n * 0.5
              and c['name'] != 'Nile']
print(f"Crossings with ≥50% of Memphis density: {len(comparable)}")

# ============================================================
# SAVE
# ============================================================
results = {
    'analysis': 'River-Perpendicular Control',
    'radius_km': RADIUS,
    'total_crossings_near_gc': len(crossings),
    'steep_crossings_gt60': len(steep),
    'unique_steep_crossings': len(steep_deduped),
    'memphis_reference': {
        'lat': memphis_lat,
        'lon': memphis_lon,
        'gc_dist_km': float(memphis_gc_dist),
        'intersection_angle': 81.3,
        'monuments_within_25km': memphis_n,
        'density_per_km2': float(memphis_density)
    },
    'crossings_with_more_than_memphis': len(higher),
    'top_crossings': steep_deduped[:20],
    'memphis_is_unique': len(higher) == 0,
    'conclusion': (
        'Memphis is the densest steep river-GC crossing'
        if len(higher) == 0
        else f'{len(higher)} crossings exceed Memphis'
    )
}

with open(os.path.join(OUT, 'river_perpendicular_control.json'), 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved: {os.path.join(OUT, 'river_perpendicular_control.json')}")
print("Done.")
