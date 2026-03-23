#!/usr/bin/env python3
"""
Directive T1-4: Coastline Intersection Analysis
================================================
Identify every point where the Great Circle crosses a coastline.
At each crossing, compute intersection angle, distance to nearest river mouth,
and continental vs island classification.
"""

import json, math, os, sys
import numpy as np

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "coastline_intersections")
os.makedirs(OUT_DIR, exist_ok=True)

POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
N_RANDOM = 200

def haversine_km(lat1, lon1, lat2, lon2):
    R = EARTH_R_KM
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(min(1.0, a)))

def great_circle_points(pole_lat, pole_lon, n_points=3600):
    pole_lat_r = math.radians(pole_lat)
    pole_lon_r = math.radians(pole_lon)
    bearings = np.linspace(0, 2 * math.pi, n_points, endpoint=False)
    lats = np.empty(n_points)
    lons = np.empty(n_points)
    for i, b in enumerate(bearings):
        lat = math.asin(
            math.sin(pole_lat_r) * math.cos(math.pi/2) +
            math.cos(pole_lat_r) * math.sin(math.pi/2) * math.cos(b)
        )
        lon = pole_lon_r + math.atan2(
            math.sin(b) * math.sin(math.pi/2) * math.cos(pole_lat_r),
            math.cos(math.pi/2) - math.sin(pole_lat_r) * math.sin(lat)
        )
        lats[i] = math.degrees(lat)
        lons[i] = math.degrees(lon)
    lons = ((lons + 180) % 360) - 180
    return lats, lons

def bearing_at_point(lat1, lon1, lat2, lon2):
    """Compute bearing from point 1 to point 2 in degrees."""
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(dlon)
    return math.degrees(math.atan2(x, y)) % 360

def random_pole():
    z = np.random.uniform(0, 1)
    theta = np.random.uniform(0, 2 * math.pi)
    return math.degrees(math.asin(z)), math.degrees(theta) - 180

# ============================================================
# LOAD COASTLINE AND RIVERS
# ============================================================
print("Loading coastline and rivers...")
import geopandas as gpd
import shapely.geometry as sg

coast_shp = os.path.join(BASE_DIR, "data", "natural_earth", "ne_50m_coastline.shp")
river_shp = os.path.join(BASE_DIR, "data", "natural_earth", "ne_50m_rivers_lake_centerlines.shp")

coastlines = gpd.read_file(coast_shp)
rivers = gpd.read_file(river_shp)

# Extract river mouths (endpoints of rivers near coast, approximate as endpoints)
river_mouths = []
for _, row in rivers.iterrows():
    geom = row.geometry
    if geom is None:
        continue
    if geom.geom_type == 'LineString':
        coords = list(geom.coords)
        if coords:
            river_mouths.append(coords[-1])  # downstream end
            river_mouths.append(coords[0])
    elif geom.geom_type == 'MultiLineString':
        for line in geom.geoms:
            coords = list(line.coords)
            if coords:
                river_mouths.append(coords[-1])
                river_mouths.append(coords[0])

river_mouth_arr = np.array(river_mouths)  # (lon, lat) pairs
print(f"Loaded {len(coastlines)} coastline features, {len(river_mouths)} river endpoints")

# ============================================================
# FIND COASTLINE CROSSINGS
# ============================================================
def find_crossings(gc_lats, gc_lons, coastline_gdf, fast_mode=False):
    """Find all points where the GC crosses coastlines using spatial index."""
    crossings = []
    n = len(gc_lats)

    # Build GC as segments (avoiding antimeridian)
    segments = []
    seg_indices = []
    for i in range(n):
        j = (i + 1) % n
        if abs(gc_lons[j] - gc_lons[i]) > 180:
            continue
        segments.append(sg.LineString([(gc_lons[i], gc_lats[i]), (gc_lons[j], gc_lats[j])]))
        seg_indices.append(i)

    # Build GC as a multi-linestring and find all coastline intersections at once
    coast_geoms = [g for g in coastline_gdf.geometry if g is not None]
    coast_lengths = [g.length for g in coast_geoms]

    # Use STRtree spatial index for fast intersection
    from shapely.strtree import STRtree
    coast_tree = STRtree(coast_geoms)

    for seg_idx, seg in zip(seg_indices, segments):
        i = seg_idx
        j = (i + 1) % n

        # Query spatial index for potential intersecting coastlines
        candidates = coast_tree.query(seg)
        for cidx in candidates:
            coast_geom = coast_geoms[cidx]
            try:
                if not seg.intersects(coast_geom):
                    continue
                intersection = seg.intersection(coast_geom)
                points = []
                if intersection.geom_type == 'Point':
                    points = [intersection]
                elif intersection.geom_type == 'MultiPoint':
                    points = list(intersection.geoms)

                for pt in points:
                    cross_lon, cross_lat = pt.x, pt.y
                    gc_bear = bearing_at_point(gc_lats[i], gc_lons[i], gc_lats[j], gc_lons[j])

                    if fast_mode:
                        # Skip detailed angle/river analysis for random circles
                        crossings.append({
                            'lat': float(cross_lat),
                            'lon': float(cross_lon),
                            'intersection_angle': 45.0,  # placeholder
                            'nearest_river_mouth_km': None,
                            'is_continental': bool(coast_lengths[cidx] > 10),
                        })
                        continue

                    # Coastline tangent
                    d_along = coast_geom.project(pt)
                    d1 = max(0, d_along - 0.01)
                    d2 = min(coast_geom.length, d_along + 0.01)
                    p1 = coast_geom.interpolate(d1)
                    p2 = coast_geom.interpolate(d2)
                    coast_bear = bearing_at_point(p1.y, p1.x, p2.y, p2.x)

                    angle_diff = abs(gc_bear - coast_bear) % 180
                    if angle_diff > 90:
                        angle_diff = 180 - angle_diff

                    nearest_river_km = None
                    if len(river_mouth_arr) > 0:
                        dists = np.sqrt(
                            ((river_mouth_arr[:, 0] - cross_lon) * math.cos(math.radians(cross_lat)))**2 +
                            (river_mouth_arr[:, 1] - cross_lat)**2
                        ) * 111.0
                        nearest_river_km = float(np.min(dists))

                    crossings.append({
                        'lat': float(cross_lat),
                        'lon': float(cross_lon),
                        'gc_bearing': float(gc_bear),
                        'coast_bearing': float(coast_bear),
                        'intersection_angle': float(angle_diff),
                        'nearest_river_mouth_km': nearest_river_km,
                        'is_continental': bool(coast_lengths[cidx] > 10),
                        'arc_position': float(i * 360.0 / n)
                    })
            except Exception:
                continue

    return crossings

print("\nFinding Alison circle coastline crossings...")
alison_lats, alison_lons = great_circle_points(POLE_LAT, POLE_LON)
alison_crossings = find_crossings(alison_lats, alison_lons, coastlines)
print(f"Found {len(alison_crossings)} coastline crossings")

# ============================================================
# RANDOM CIRCLES FOR COMPARISON
# ============================================================
print(f"\nComputing crossings for {N_RANDOM} random circles...")
random_crossing_counts = []
random_angles = []
random_river_dists = []

for i in range(N_RANDOM):
    if (i + 1) % 50 == 0:
        print(f"  Random circle {i+1}/{N_RANDOM}")
    plat, plon = random_pole()
    rlats, rlons = great_circle_points(plat, plon)
    rc = find_crossings(rlats, rlons, coastlines, fast_mode=True)
    random_crossing_counts.append(len(rc))

# ============================================================
# ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("RESULTS")
print("=" * 70)

alison_angles = [c['intersection_angle'] for c in alison_crossings]
alison_river_dists = [c['nearest_river_mouth_km'] for c in alison_crossings if c['nearest_river_mouth_km'] is not None]
continental_fraction = sum(1 for c in alison_crossings if c['is_continental']) / max(1, len(alison_crossings))

results = {
    'n_crossings': len(alison_crossings),
    'random_crossing_mean': float(np.mean(random_crossing_counts)),
    'random_crossing_std': float(np.std(random_crossing_counts)),
    'mean_intersection_angle': float(np.mean(alison_angles)) if alison_angles else None,
    'random_mean_angle': float(np.mean(random_angles)) if random_angles else None,
    'mean_river_mouth_distance_km': float(np.mean(alison_river_dists)) if alison_river_dists else None,
    'random_mean_river_dist': float(np.mean(random_river_dists)) if random_river_dists else None,
    'continental_fraction': continental_fraction,
    'steep_fraction_alison': sum(1 for a in alison_angles if a > 60) / max(1, len(alison_angles)),
    'shallow_fraction_alison': sum(1 for a in alison_angles if a < 30) / max(1, len(alison_angles)),
    'steep_fraction_random': sum(1 for a in random_angles if a > 60) / max(1, len(random_angles)) if random_angles else None,
    'shallow_fraction_random': sum(1 for a in random_angles if a < 30) / max(1, len(random_angles)) if random_angles else None,
}

print(f"\n1. CROSSING COUNT")
print(f"   Alison: {results['n_crossings']}")
print(f"   Random: {results['random_crossing_mean']:.1f} ± {results['random_crossing_std']:.1f}")

print(f"\n2. INTERSECTION ANGLES")
print(f"   Alison mean angle: {results['mean_intersection_angle']:.1f}°" if results['mean_intersection_angle'] else "   Alison mean angle: N/A")
print(f"   Random mean angle: {results['random_mean_angle']:.1f}°" if results['random_mean_angle'] else "   Random mean angle: N/A (fast mode)")
print(f"   Alison steep (>60°): {results['steep_fraction_alison']*100:.1f}%")
print(f"   Alison shallow (<30°): {results['shallow_fraction_alison']*100:.1f}%")
if results['steep_fraction_random'] is not None:
    print(f"   Random steep: {results['steep_fraction_random']*100:.1f}%")
    print(f"   Random shallow: {results['shallow_fraction_random']*100:.1f}%")

print(f"\n3. RIVER MOUTH PROXIMITY")
if results['mean_river_mouth_distance_km'] is not None:
    print(f"   Alison mean distance to river mouth: {results['mean_river_mouth_distance_km']:.1f} km")
if results['random_mean_river_dist'] is not None:
    print(f"   Random mean: {results['random_mean_river_dist']:.1f} km")

print(f"\n4. CONTINENTAL VS ISLAND")
print(f"   Alison continental fraction: {results['continental_fraction']*100:.1f}%")

# ============================================================
# SAVE OUTPUTS
# ============================================================
with open(os.path.join(OUT_DIR, "crossings.json"), 'w') as f:
    json.dump({'alison_crossings': alison_crossings, 'summary': results}, f, indent=2)

# ============================================================
# PLOTS
# ============================================================
print("\nGenerating plots...")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Crossing map
fig, ax = plt.subplots(figsize=(16, 8))
# Plot coastlines
for _, row in coastlines.iterrows():
    geom = row.geometry
    if geom is None:
        continue
    if geom.geom_type == 'LineString':
        x, y = geom.xy
        ax.plot(x, y, color='gray', linewidth=0.3)
    elif geom.geom_type == 'MultiLineString':
        for line in geom.geoms:
            x, y = line.xy
            ax.plot(x, y, color='gray', linewidth=0.3)

# Plot GC
ax.plot(alison_lons, alison_lats, 'b-', linewidth=0.5, alpha=0.5, label='Great Circle')

# Plot crossings
for c in alison_crossings:
    color = 'red' if c['is_continental'] else 'orange'
    ax.plot(c['lon'], c['lat'], 'o', color=color, markersize=5)

ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title(f'Alison Great Circle — Coastline Crossings (n={len(alison_crossings)})')
ax.set_xlim(-180, 180)
ax.set_ylim(-60, 80)
from matplotlib.patches import Patch
legend_elements = [
    plt.Line2D([0], [0], color='blue', alpha=0.5, label='Great Circle'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label='Continental crossing'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=8, label='Island crossing')
]
ax.legend(handles=legend_elements)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "crossing_map.png"), dpi=150)
plt.close()
print("  Saved crossing_map.png")

# Angle histogram
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax = axes[0]
ax.hist(alison_angles, bins=18, range=(0, 90), alpha=0.7, color='red', label='Alison', density=True)
if random_angles:
    ax.hist(random_angles, bins=18, range=(0, 90), alpha=0.5, color='gray', label='Random', density=True)
ax.set_xlabel('Intersection Angle (degrees)')
ax.set_ylabel('Density')
ax.set_title('Coastline Intersection Angles')
ax.axvline(x=30, color='green', linestyle='--', alpha=0.5, label='Shallow threshold')
ax.axvline(x=60, color='orange', linestyle='--', alpha=0.5, label='Steep threshold')
ax.legend()

ax = axes[1]
ax.hist(alison_river_dists, bins=30, alpha=0.7, color='red', label='Alison crossings', density=True)
if random_river_dists:
    ax.hist(random_river_dists, bins=30, alpha=0.5, color='gray', label='Random crossings', density=True)
ax.set_xlabel('Distance to Nearest River Mouth (km)')
ax.set_ylabel('Density')
ax.set_title('Crossing Distance to River Mouths')
ax.legend()

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "angle_histogram.png"), dpi=150)
plt.close()
print("  Saved angle_histogram.png")

print("\nDone! Outputs saved to", OUT_DIR)
