#!/usr/bin/env python3
"""
Directive T1-3: Topographic/Elevation Profile
==============================================
Extract elevation along the Alison Great Circle at 0.1° intervals (~3,600 points).
Compare against 1,000 random great circles.

Tests:
1. Mean elevation comparison
2. Elevation variance comparison
3. Walkability index (consecutive elevation changes > 500m)
4. River crossing count (Strahler order >= 5 proxy via Natural Earth rivers)
5. Maximum elevation along path
6. Elevation at monument clusters vs between clusters
"""

import json, math, os, sys, struct, time
import numpy as np
import netCDF4 as nc

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "elevation_profile")
os.makedirs(OUT_DIR, exist_ok=True)

POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
ARC_STEP_DEG = 0.1
N_RANDOM = 1000

# ============================================================
# LOAD ETOPO1
# ============================================================
print("Loading ETOPO1...")
etopo_path = os.path.join(BASE_DIR, "data", "geophysical", "etopo", "ETOPO1_Ice_g_gmt4.grd")
ds = nc.Dataset(etopo_path)
etopo_lons = ds.variables['x'][:]
etopo_lats = ds.variables['y'][:]
etopo_z = ds.variables['z']
etopo_lon_min, etopo_lon_step = float(etopo_lons[0]), float(etopo_lons[1] - etopo_lons[0])
etopo_lat_min, etopo_lat_step = float(etopo_lats[0]), float(etopo_lats[1] - etopo_lats[0])
print(f"ETOPO1 loaded: {len(etopo_lats)} x {len(etopo_lons)}")

def get_elevations_batch(lats, lons):
    lat_idx = np.round((lats - etopo_lat_min) / etopo_lat_step).astype(int)
    lon_idx = np.round((lons - etopo_lon_min) / etopo_lon_step).astype(int)
    lat_idx = np.clip(lat_idx, 0, len(etopo_lats) - 1)
    lon_idx = np.clip(lon_idx, 0, len(etopo_lons) - 1)
    elevs = np.empty(len(lats))
    for i in range(len(lats)):
        elevs[i] = float(etopo_z[lat_idx[i], lon_idx[i]])
    return elevs

# ============================================================
# GREAT CIRCLE POINT GENERATION
# ============================================================
def great_circle_points(pole_lat, pole_lon, n_points=3600):
    """Generate n_points evenly spaced along the great circle defined by pole."""
    pole_lat_r = math.radians(pole_lat)
    pole_lon_r = math.radians(pole_lon)
    bearings = np.linspace(0, 2 * math.pi, n_points, endpoint=False)

    # Point on GC = rotate pole by 90° along each bearing
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

    # Normalize longitudes to [-180, 180]
    lons = ((lons + 180) % 360) - 180
    return lats, lons

def random_pole():
    """Generate a uniform random pole on the sphere (upper hemisphere for uniqueness)."""
    z = np.random.uniform(0, 1)  # upper hemisphere
    theta = np.random.uniform(0, 2 * math.pi)
    lat = math.degrees(math.asin(z))
    lon = math.degrees(theta) - 180
    return lat, lon

# ============================================================
# COMPUTE ELEVATION STATS FOR A CIRCLE
# ============================================================
def compute_circle_stats(pole_lat, pole_lon, n_points=3600):
    lats, lons = great_circle_points(pole_lat, pole_lon, n_points)
    elevs = get_elevations_batch(lats, lons)

    mean_elev = float(np.mean(elevs))
    var_elev = float(np.var(elevs))
    max_elev = float(np.max(elevs))

    # Walkability: count consecutive changes > 500m
    diffs = np.abs(np.diff(elevs))
    walkability_count = int(np.sum(diffs > 500))

    return {
        'mean_elevation': mean_elev,
        'elevation_variance': var_elev,
        'max_elevation': max_elev,
        'walkability_violations': walkability_count,
        'elevations': elevs,
        'lats': lats,
        'lons': lons
    }

# ============================================================
# LOAD MONUMENTS FOR CLUSTER ANALYSIS
# ============================================================
import csv

def haversine_km(lat1, lon1, lat2, lon2):
    R = EARTH_R_KM
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(min(1.0, a)))

def dist_from_gc(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    d = haversine_km(lat, lon, pole_lat, pole_lon)
    return abs(d - QUARTER_CIRC)

MONUMENTAL_TYPES = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus", "shrine"
}

def load_monuments():
    """Load monumental sites from Pleiades within 100km of GC."""
    path = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")
    monuments = []
    with open(path, encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row['reprLat'])
                lon = float(row['reprLong'])
            except (ValueError, KeyError, TypeError):
                continue
            if lat == 0 and lon == 0:
                continue
            ftypes = {t.strip() for t in row.get('featureTypes', '').split(',')}
            if ftypes & MONUMENTAL_TYPES:
                d = dist_from_gc(lat, lon)
                if d < 100:
                    monuments.append((lat, lon, d))
    return monuments

# ============================================================
# RIVER CROSSING ANALYSIS (Natural Earth rivers)
# ============================================================
def load_rivers_and_count_crossings(gc_lats, gc_lons):
    """Count major river crossings using Natural Earth rivers shapefile."""
    try:
        import shapely.geometry as sg
        import json as json_

        rivers_shp = os.path.join(BASE_DIR, "data", "natural_earth", "ne_50m_rivers_lake_centerlines.shp")

        # Read shapefile using basic approach
        try:
            import geopandas as gpd
            rivers = gpd.read_file(rivers_shp)
            # Filter for major rivers (scalerank <= 3 is a proxy for Strahler >= 5)
            major_rivers = rivers[rivers['scalerank'] <= 3]
            river_geoms = list(major_rivers.geometry)
        except Exception:
            print("Could not load rivers shapefile, skipping river analysis")
            return None

        # Build GC as a linestring
        gc_points = list(zip(gc_lons, gc_lats))

        # Split into segments that don't cross the antimeridian
        segments = []
        current_seg = [gc_points[0]]
        for i in range(1, len(gc_points)):
            if abs(gc_points[i][0] - gc_points[i-1][0]) > 180:
                if len(current_seg) > 1:
                    segments.append(current_seg)
                current_seg = [gc_points[i]]
            else:
                current_seg.append(gc_points[i])
        if len(current_seg) > 1:
            segments.append(current_seg)

        crossing_count = 0
        for seg in segments:
            gc_line = sg.LineString(seg)
            for river in river_geoms:
                try:
                    if gc_line.intersects(river):
                        intersection = gc_line.intersection(river)
                        if intersection.is_empty:
                            continue
                        if intersection.geom_type == 'Point':
                            crossing_count += 1
                        elif intersection.geom_type == 'MultiPoint':
                            crossing_count += len(list(intersection.geoms))
                        else:
                            crossing_count += 1
                except Exception:
                    continue

        return crossing_count
    except Exception as e:
        print(f"River analysis error: {e}")
        return None

# ============================================================
# MAIN EXECUTION
# ============================================================
print("\n" + "=" * 70)
print("COMPUTING ALISON CIRCLE ELEVATION PROFILE")
print("=" * 70)

alison = compute_circle_stats(POLE_LAT, POLE_LON)
print(f"Alison circle: mean={alison['mean_elevation']:.1f}m, "
      f"var={alison['elevation_variance']:.0f}, "
      f"max={alison['max_elevation']:.0f}m, "
      f"walkability_violations={alison['walkability_violations']}")

# River crossings for Alison circle
print("\nComputing river crossings for Alison circle...")
alison_rivers = load_rivers_and_count_crossings(alison['lats'], alison['lons'])
print(f"Alison circle river crossings: {alison_rivers}")

# Monument cluster elevation analysis
print("\nLoading monuments for cluster analysis...")
monuments = load_monuments()
print(f"Found {len(monuments)} monuments within 100km of GC")

# Find GC points near monument clusters (within 200km of a monument)
gc_bearings = np.linspace(0, 360, 3600, endpoint=False)
near_monument = np.zeros(3600, dtype=bool)
for m_lat, m_lon, m_dist in monuments:
    for i in range(3600):
        d = haversine_km(alison['lats'][i], alison['lons'][i], m_lat, m_lon)
        if d < 200:
            near_monument[i] = True

elev_at_monuments = alison['elevations'][near_monument]
elev_between = alison['elevations'][~near_monument]
print(f"GC points near monuments: {near_monument.sum()}, between: {(~near_monument).sum()}")
print(f"Mean elevation at monument clusters: {np.mean(elev_at_monuments):.1f}m")
print(f"Mean elevation between clusters: {np.mean(elev_between):.1f}m")

# ============================================================
# RANDOM CIRCLES
# ============================================================
print(f"\nComputing {N_RANDOM} random circles...")
random_stats = {
    'mean_elevations': [],
    'elevation_variances': [],
    'max_elevations': [],
    'walkability_violations': [],
    'river_crossings': []
}

for i in range(N_RANDOM):
    if (i + 1) % 100 == 0:
        print(f"  Random circle {i+1}/{N_RANDOM}")
    plat, plon = random_pole()
    stats = compute_circle_stats(plat, plon)
    random_stats['mean_elevations'].append(stats['mean_elevation'])
    random_stats['elevation_variances'].append(stats['elevation_variance'])
    random_stats['max_elevations'].append(stats['max_elevation'])
    random_stats['walkability_violations'].append(stats['walkability_violations'])

# River crossings for a subset (slow)
print("Computing river crossings for 100 random circles...")
for i in range(100):
    if (i + 1) % 20 == 0:
        print(f"  River crossing {i+1}/100")
    plat, plon = random_pole()
    rlats, rlons = great_circle_points(plat, plon)
    rc = load_rivers_and_count_crossings(rlats, rlons)
    if rc is not None:
        random_stats['river_crossings'].append(rc)

# ============================================================
# STATISTICAL TESTS
# ============================================================
print("\n" + "=" * 70)
print("RESULTS")
print("=" * 70)

def percentile_rank(value, distribution):
    return 100.0 * sum(1 for x in distribution if x <= value) / len(distribution)

results = {}

# Test 1: Mean elevation
rm = random_stats['mean_elevations']
pct = percentile_rank(alison['mean_elevation'], rm)
results['mean_elevation'] = {
    'alison': alison['mean_elevation'],
    'random_mean': float(np.mean(rm)),
    'random_std': float(np.std(rm)),
    'percentile': pct,
    'interpretation': 'lower than random' if pct < 50 else 'higher than random'
}
print(f"\n1. MEAN ELEVATION")
print(f"   Alison: {alison['mean_elevation']:.1f}m")
print(f"   Random: {np.mean(rm):.1f} ± {np.std(rm):.1f}m")
print(f"   Percentile: {pct:.1f}% — {'traces lowlands' if pct < 30 else 'no clear lowland preference'}")

# Test 2: Elevation variance
rv = random_stats['elevation_variances']
pct = percentile_rank(alison['elevation_variance'], rv)
results['elevation_variance'] = {
    'alison': alison['elevation_variance'],
    'random_mean': float(np.mean(rv)),
    'random_std': float(np.std(rv)),
    'percentile': pct,
    'interpretation': 'smoother terrain' if pct < 50 else 'rougher terrain'
}
print(f"\n2. ELEVATION VARIANCE")
print(f"   Alison: {alison['elevation_variance']:.0f}")
print(f"   Random: {np.mean(rv):.0f} ± {np.std(rv):.0f}")
print(f"   Percentile: {pct:.1f}% — {'smoother' if pct < 30 else 'not notably smooth'}")

# Test 3: Walkability
rw = random_stats['walkability_violations']
pct = percentile_rank(alison['walkability_violations'], rw)
results['walkability'] = {
    'alison_violations': alison['walkability_violations'],
    'random_mean': float(np.mean(rw)),
    'random_std': float(np.std(rw)),
    'percentile': pct,
    'interpretation': 'more walkable' if pct < 50 else 'less walkable'
}
print(f"\n3. WALKABILITY (>500m jumps)")
print(f"   Alison: {alison['walkability_violations']} violations")
print(f"   Random: {np.mean(rw):.1f} ± {np.std(rw):.1f}")
print(f"   Percentile: {pct:.1f}%")

# Test 4: River crossings
if random_stats['river_crossings']:
    rr = random_stats['river_crossings']
    pct = percentile_rank(alison_rivers, rr)
    results['river_crossings'] = {
        'alison': alison_rivers,
        'random_mean': float(np.mean(rr)),
        'random_std': float(np.std(rr)),
        'percentile': pct,
        'n_random_sampled': len(rr)
    }
    print(f"\n4. RIVER CROSSINGS (major rivers)")
    print(f"   Alison: {alison_rivers}")
    print(f"   Random: {np.mean(rr):.1f} ± {np.std(rr):.1f}")
    print(f"   Percentile: {pct:.1f}%")

# Test 5: Max elevation
rmax = random_stats['max_elevations']
pct = percentile_rank(alison['max_elevation'], rmax)
results['max_elevation'] = {
    'alison': alison['max_elevation'],
    'random_mean': float(np.mean(rmax)),
    'random_std': float(np.std(rmax)),
    'percentile': pct
}
print(f"\n5. MAX ELEVATION")
print(f"   Alison: {alison['max_elevation']:.0f}m")
print(f"   Random: {np.mean(rmax):.0f} ± {np.std(rmax):.0f}m")
print(f"   Percentile: {pct:.1f}%")

# Test 6: Elevation at monument clusters vs between
results['monument_cluster_elevation'] = {
    'near_monument_mean': float(np.mean(elev_at_monuments)),
    'between_mean': float(np.mean(elev_between)),
    'near_monument_n': int(near_monument.sum()),
    'between_n': int((~near_monument).sum()),
    'difference': float(np.mean(elev_at_monuments) - np.mean(elev_between))
}
print(f"\n6. ELEVATION AT MONUMENT CLUSTERS vs BETWEEN")
print(f"   Near monuments: {np.mean(elev_at_monuments):.1f}m (n={near_monument.sum()})")
print(f"   Between: {np.mean(elev_between):.1f}m (n={(~near_monument).sum()})")
print(f"   Difference: {np.mean(elev_at_monuments) - np.mean(elev_between):.1f}m")

# ============================================================
# SAVE OUTPUTS
# ============================================================
# Save elevation sweep
sweep_data = {
    'arc_positions_deg': gc_bearings.tolist(),
    'lats': alison['lats'].tolist(),
    'lons': alison['lons'].tolist(),
    'elevations_m': alison['elevations'].tolist()
}
with open(os.path.join(OUT_DIR, "elevation_sweep.json"), 'w') as f:
    json.dump(sweep_data, f, indent=2)

# Save walkability comparison
walkability_data = {
    'alison': results.get('walkability', {}),
    'mean_elevation': results.get('mean_elevation', {}),
    'elevation_variance': results.get('elevation_variance', {}),
    'max_elevation': results.get('max_elevation', {}),
    'river_crossings': results.get('river_crossings', {}),
    'monument_cluster_elevation': results.get('monument_cluster_elevation', {}),
    'random_distributions': {
        'mean_elevations': random_stats['mean_elevations'],
        'walkability_violations': random_stats['walkability_violations'],
        'max_elevations': random_stats['max_elevations']
    }
}
with open(os.path.join(OUT_DIR, "walkability_comparison.json"), 'w') as f:
    json.dump(walkability_data, f, indent=2)

# ============================================================
# PLOTS
# ============================================================
print("\nGenerating plots...")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Plot 1: Full 360° elevation profile with monument locations
fig, ax = plt.subplots(figsize=(16, 6))
ax.fill_between(gc_bearings, 0, alison['elevations'], alpha=0.3, color='brown', label='Elevation')
ax.plot(gc_bearings, alison['elevations'], color='brown', linewidth=0.5)
ax.axhline(y=0, color='blue', linewidth=0.5, alpha=0.5)

# Mark monument cluster regions
for i in range(3600):
    if near_monument[i]:
        ax.axvline(x=gc_bearings[i], color='red', alpha=0.02, linewidth=1)

# Add some landmark labels
ax.set_xlabel('Arc Position (degrees)')
ax.set_ylabel('Elevation (m)')
ax.set_title('Alison Great Circle — Elevation Profile (360°)')
ax.set_xlim(0, 360)
ax.legend(['Elevation', 'Sea level'], loc='upper right')

# Add red patch for legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='brown', alpha=0.3, label='Elevation'),
    plt.Line2D([0], [0], color='blue', alpha=0.5, label='Sea level'),
    Patch(facecolor='red', alpha=0.3, label='Monument cluster regions')
]
ax.legend(handles=legend_elements, loc='upper right')

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "elevation_plot.png"), dpi=150)
plt.close()
print("  Saved elevation_plot.png")

# Plot 2: Histograms vs random
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Mean elevation histogram
ax = axes[0, 0]
ax.hist(random_stats['mean_elevations'], bins=50, alpha=0.7, color='gray', label='Random circles')
ax.axvline(x=alison['mean_elevation'], color='red', linewidth=2, label=f'Alison ({alison["mean_elevation"]:.0f}m)')
ax.set_xlabel('Mean Elevation (m)')
ax.set_ylabel('Count')
ax.set_title('Mean Elevation Distribution')
ax.legend()

# Variance histogram
ax = axes[0, 1]
ax.hist(random_stats['elevation_variances'], bins=50, alpha=0.7, color='gray', label='Random circles')
ax.axvline(x=alison['elevation_variance'], color='red', linewidth=2, label=f'Alison ({alison["elevation_variance"]:.0f})')
ax.set_xlabel('Elevation Variance')
ax.set_ylabel('Count')
ax.set_title('Elevation Variance Distribution')
ax.legend()

# Walkability histogram
ax = axes[1, 0]
ax.hist(random_stats['walkability_violations'], bins=50, alpha=0.7, color='gray', label='Random circles')
ax.axvline(x=alison['walkability_violations'], color='red', linewidth=2, label=f'Alison ({alison["walkability_violations"]})')
ax.set_xlabel('>500m Elevation Jumps')
ax.set_ylabel('Count')
ax.set_title('Walkability Violations')
ax.legend()

# Max elevation histogram
ax = axes[1, 1]
ax.hist(random_stats['max_elevations'], bins=50, alpha=0.7, color='gray', label='Random circles')
ax.axvline(x=alison['max_elevation'], color='red', linewidth=2, label=f'Alison ({alison["max_elevation"]:.0f}m)')
ax.set_xlabel('Max Elevation (m)')
ax.set_ylabel('Count')
ax.set_title('Maximum Elevation Distribution')
ax.legend()

plt.suptitle('Alison Great Circle vs 1000 Random Circles — Elevation Statistics', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "histogram_vs_random.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved histogram_vs_random.png")

print("\nDone! All outputs saved to", OUT_DIR)
