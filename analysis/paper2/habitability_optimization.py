#!/usr/bin/env python3
"""
Test 2: Habitability Optimization
==================================
Does the great circle that maximizes traversal through habitable terrain
match the Alison circle?

Uses a proxy habitability score (no GAEZ/MODIS available):
  - Land fraction along the circle
  - Distance-to-coast bonus (coastal = more habitable)
  - Latitude penalty (extreme latitudes less habitable)
  - Elevation suitability (avoid extremes)
  - River crossings (water access)

For each candidate pole on a 2-degree grid:
  1. Trace the great circle (360 points at 1-degree steps)
  2. Sample habitability at each point
  3. Compute composite score

Score = mean_habitability * land_fraction * (1 + log(river_crossings + 1))

Output:
  - outputs/resolution_tests/habitability_optimization.json
  - outputs/resolution_tests/habitability_heatmap.png
"""

import json
import math
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Constants ──────────────────────────────────────────────────────────────
EARTH_R_KM = 6371.0
QUARTER_CIRC = math.pi * EARTH_R_KM / 2  # ~10,008 km
ALISON_POLE = (59.682, -138.646)
DEG2RAD = math.pi / 180.0

# ── Data paths ─────────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAND_GRID = os.path.join(BASE, 'data', 'land_grid_05deg.npy')
RIVER_SHP = os.path.join(BASE, 'data', 'natural_earth',
                         'ne_50m_rivers_lake_centerlines')
COAST_SHP = os.path.join(BASE, 'data', 'natural_earth', 'ne_50m_coastline')
ETOPO_NC  = os.path.join(BASE, 'ml_features', 'etopo2022.nc')
OUT_DIR   = os.path.join(BASE, 'outputs', 'resolution_tests')
os.makedirs(OUT_DIR, exist_ok=True)


# ── Load land mask (0.5-degree, 360x720) ──────────────────────────────────
print("Loading land mask...")
land_grid = np.load(LAND_GRID)  # shape (360, 720), uint8, 0/1
# Convention: row 0 = 89.75°N, row 359 = -89.75°N  (top-to-bottom)
# col 0 = -179.75°E, col 719 = 179.75°E


def latlon_to_grid(lat, lon):
    """Convert lat/lon to (row, col) in the 0.5-degree grid."""
    row = int((89.75 - lat) / 0.5)
    col = int((lon + 179.75) / 0.5)
    row = max(0, min(359, row))
    col = max(0, min(719, col))
    return row, col


def is_land(lat, lon):
    r, c = latlon_to_grid(lat, lon)
    return bool(land_grid[r, c])


# ── Build distance-to-coast raster ────────────────────────────────────────
print("Computing distance-to-coast raster (from land grid edges)...")
from scipy.ndimage import distance_transform_edt

# Coastal cells: land cells adjacent to ocean
land_bool = land_grid.astype(bool)
# Distance transform on ocean (distance of each ocean cell to nearest land)
# We want: for each LAND cell, how far to the nearest coast (ocean edge)?
# A land cell on the coast has distance 0.
ocean_mask = ~land_bool
# For land cells: distance to nearest ocean cell (in grid units, ~0.5°)
coast_dist_grid = distance_transform_edt(land_bool).astype(np.float32)
# Convert to rough km (0.5° ≈ 55.5 km at equator, less at poles)
# We'll use grid-unit distance as proxy — no need for exact km
# Normalize to 0-1 where 0 = coast, 1 = deep interior
max_coast_dist = coast_dist_grid.max()
if max_coast_dist > 0:
    coast_proximity = 1.0 - (coast_dist_grid / max_coast_dist)
else:
    coast_proximity = np.ones_like(coast_dist_grid)
# Only meaningful on land
coast_proximity[~land_bool] = 0.0


# ── Build latitude suitability raster ─────────────────────────────────────
print("Building latitude suitability raster...")
lats_1d = np.linspace(89.75, -89.75, 360)
# Habitability peaks in temperate/tropical zones, drops at poles & equatorial extremes
# Simple model: Gaussian centered around ±30° (agricultural sweet spots)
lat_suit = np.zeros(360, dtype=np.float32)
for i, lat in enumerate(lats_1d):
    abs_lat = abs(lat)
    # Peak suitability 15-45°, drops toward poles and deep tropics
    if abs_lat < 60:
        # Cosine-weighted with tropical/temperate bonus
        lat_suit[i] = np.cos(np.radians(abs_lat)) * (1.0 - 0.3 * np.exp(-((abs_lat - 30) / 20) ** 2))
    else:
        lat_suit[i] = max(0, np.cos(np.radians(abs_lat)) * 0.5)
lat_suit = lat_suit / lat_suit.max()  # normalize to 0-1
lat_suit_grid = np.tile(lat_suit[:, np.newaxis], (1, 720))


# ── Load elevation data for suitability ───────────────────────────────────
print("Loading elevation data...")
try:
    import netCDF4 as nc
    ds = nc.Dataset(ETOPO_NC, 'r')
    # ETOPO2022 variables
    elev_var = None
    for vname in ['z', 'elevation', 'Band1']:
        if vname in ds.variables:
            elev_var = vname
            break
    if elev_var is None:
        raise KeyError(f"No elevation variable found. Variables: {list(ds.variables.keys())}")

    elev_lats = ds.variables['lat'][:]
    elev_lons = ds.variables['lon'][:]
    # Subsample to ~0.5° for speed (ETOPO2022 is 60-arcsec = ~0.017°)
    step = max(1, len(elev_lats) // 360)
    elev_data = ds.variables[elev_var][::step, ::step]
    elev_lats_sub = elev_lats[::step]
    elev_lons_sub = elev_lons[::step]
    ds.close()

    # Build elevation suitability: penalize very high (>3000m) and below sea level
    elev_suit_raw = np.clip(elev_data, 0, 5000).astype(np.float32)
    # Suitability: 1.0 at 0-500m, linearly decreasing to 0 at 4000m+
    elev_suit_full = np.where(elev_suit_raw <= 500, 1.0,
                              np.where(elev_suit_raw >= 4000, 0.0,
                                       1.0 - (elev_suit_raw - 500) / 3500))
    HAS_ELEVATION = True
    print(f"  Elevation loaded: {elev_data.shape}, subsampled step={step}")
except Exception as e:
    print(f"  Warning: Could not load elevation: {e}")
    HAS_ELEVATION = False


def sample_elevation_suitability(lat, lon):
    """Sample elevation suitability at a point."""
    if not HAS_ELEVATION:
        return 0.7  # neutral default
    # Find nearest indices
    lat_idx = np.argmin(np.abs(elev_lats_sub - lat))
    lon_idx = np.argmin(np.abs(elev_lons_sub - lon))
    return float(elev_suit_full[lat_idx, lon_idx])


# ── Load river segments for crossing counts ───────────────────────────────
print("Loading river data...")
import shapefile
river_reader = shapefile.Reader(RIVER_SHP)
river_segments = []
for shape in river_reader.shapes():
    pts = shape.points
    for i in range(len(pts) - 1):
        # (lon1, lat1, lon2, lat2)
        river_segments.append((pts[i][0], pts[i][1], pts[i+1][0], pts[i+1][1]))
river_segments = np.array(river_segments, dtype=np.float64)
print(f"  Loaded {len(river_segments)} river segments")


# ── Build composite habitability raster at 0.5° ──────────────────────────
print("Building composite habitability raster...")
# Combine: land * (0.4*coast_proximity + 0.3*lat_suitability + 0.3*elev_suitability)
hab_raster = np.zeros((360, 720), dtype=np.float32)

if HAS_ELEVATION:
    # Resample elevation suitability to 0.5° grid
    elev_suit_05 = np.zeros((360, 720), dtype=np.float32)
    for i in range(360):
        lat = 89.75 - i * 0.5
        lat_idx = np.argmin(np.abs(elev_lats_sub - lat))
        for j in range(720):
            lon = -179.75 + j * 0.5
            lon_idx = np.argmin(np.abs(elev_lons_sub - lon))
            elev_suit_05[i, j] = elev_suit_full[lat_idx, lon_idx]
    print("  Elevation resampled to 0.5° grid")
else:
    elev_suit_05 = np.full((360, 720), 0.7, dtype=np.float32)

hab_raster = land_grid.astype(np.float32) * (
    0.35 * coast_proximity +
    0.35 * lat_suit_grid +
    0.30 * elev_suit_05
)
print(f"  Habitability raster: min={hab_raster.min():.3f}, max={hab_raster.max():.3f}, "
      f"mean(land)={hab_raster[land_bool].mean():.3f}")


# ── Great circle tracing functions ────────────────────────────────────────
def great_circle_points(pole_lat, pole_lon, n_points=360):
    """
    Generate n_points along the great circle defined by the given pole.
    The great circle is the 'equator' relative to this pole.
    """
    # Convert pole to Cartesian
    p_lat = pole_lat * DEG2RAD
    p_lon = pole_lon * DEG2RAD

    # Pole unit vector
    px = math.cos(p_lat) * math.cos(p_lon)
    py = math.cos(p_lat) * math.sin(p_lon)
    pz = math.sin(p_lat)
    pole = np.array([px, py, pz])

    # Find two orthogonal vectors in the plane perpendicular to pole
    # Pick a reference not parallel to pole
    if abs(pz) < 0.9:
        ref = np.array([0, 0, 1.0])
    else:
        ref = np.array([1.0, 0, 0])

    u = np.cross(pole, ref)
    u = u / np.linalg.norm(u)
    v = np.cross(pole, u)
    v = v / np.linalg.norm(v)

    # Parameterize: point(t) = cos(t)*u + sin(t)*v
    angles = np.linspace(0, 2 * math.pi, n_points, endpoint=False)
    points = np.outer(np.cos(angles), u) + np.outer(np.sin(angles), v)

    # Convert back to lat/lon
    lats = np.degrees(np.arcsin(np.clip(points[:, 2], -1, 1)))
    lons = np.degrees(np.arctan2(points[:, 1], points[:, 0]))

    return lats, lons


def sample_habitability_along_gc(pole_lat, pole_lon, n_points=360):
    """
    Sample the habitability raster along a great circle.
    Returns: mean_hab, land_fraction, hab_values
    """
    lats, lons = great_circle_points(pole_lat, pole_lon, n_points)

    # Convert to grid indices
    rows = np.clip(((89.75 - lats) / 0.5).astype(int), 0, 359)
    cols = np.clip(((lons + 179.75) / 0.5).astype(int), 0, 719)

    hab_values = hab_raster[rows, cols]
    land_values = land_grid[rows, cols]

    land_fraction = land_values.mean()
    # Mean habitability (over all points, including ocean=0)
    mean_hab = hab_values.mean()
    # Mean habitability over land points only
    if land_values.sum() > 0:
        mean_hab_land = hab_values[land_values > 0].mean()
    else:
        mean_hab_land = 0.0

    return mean_hab, land_fraction, mean_hab_land, lats, lons, land_values


def count_river_crossings(lats, lons, land_values, threshold_deg=1.0):
    """
    Count approximate river crossings along the great circle path.
    For each land segment of the path, check proximity to river segments.
    """
    # Only check land points
    land_mask = land_values > 0
    if land_mask.sum() < 2:
        return 0

    land_lats = lats[land_mask]
    land_lons = lons[land_mask]

    # For efficiency, use bounding-box pre-filter
    crossings = 0
    # Check each consecutive pair of GC points on land
    for i in range(len(land_lats) - 1):
        lat1, lon1 = land_lats[i], land_lons[i]
        lat2, lon2 = land_lats[i + 1], land_lons[i + 1]

        # Midpoint
        mlat = (lat1 + lat2) / 2
        mlon = (lon1 + lon2) / 2

        # Find river segments near this midpoint
        dlat = np.abs(river_segments[:, 1] - mlat) + np.abs(river_segments[:, 3] - mlat)
        dlon = np.abs(river_segments[:, 0] - mlon) + np.abs(river_segments[:, 2] - mlon)
        nearby = (dlat < threshold_deg * 4) & (dlon < threshold_deg * 4)

        if nearby.any():
            # More precise: minimum distance from midpoint to river segment midpoints
            r_mlat = (river_segments[nearby, 1] + river_segments[nearby, 3]) / 2
            r_mlon = (river_segments[nearby, 0] + river_segments[nearby, 2]) / 2
            dist = np.sqrt((r_mlat - mlat) ** 2 + (r_mlon - mlon) ** 2)
            if dist.min() < threshold_deg:
                crossings += 1

    return crossings


# ── Main pole grid search ─────────────────────────────────────────────────
def run_optimization():
    print("\n" + "=" * 70)
    print("HABITABILITY OPTIMIZATION — 2° Pole Grid Search")
    print("=" * 70)

    # Grid: all poles on 2° grid (full sphere, but only unique hemispheres)
    # A pole at (lat, lon) gives the same GC as (-lat, lon+180), so restrict
    # to lat >= 0 to avoid duplicates
    lat_range = np.arange(0, 90.01, 2.0)
    lon_range = np.arange(-180, 180, 2.0)
    n_poles = len(lat_range) * len(lon_range)
    print(f"Testing {n_poles} candidate poles ({len(lat_range)} lats x {len(lon_range)} lons)")

    results = np.zeros((len(lat_range), len(lon_range)), dtype=np.float32)
    all_scores = []

    t0 = time.time()
    best_score = -1
    best_pole = None

    for i, plat in enumerate(lat_range):
        for j, plon in enumerate(lon_range):
            mean_hab, land_frac, mean_hab_land, lats, lons, land_vals = \
                sample_habitability_along_gc(plat, plon, n_points=360)

            # River crossings (expensive — use for top candidates only in pass 2)
            # For pass 1, use land_fraction * mean_hab as score
            score = mean_hab * land_frac
            results[i, j] = score
            all_scores.append({
                'pole_lat': float(plat),
                'pole_lon': float(plon),
                'mean_habitability': float(mean_hab),
                'land_fraction': float(land_frac),
                'mean_hab_land_only': float(mean_hab_land),
                'score_pass1': float(score)
            })

            if score > best_score:
                best_score = score
                best_pole = (plat, plon)

        if (i + 1) % 5 == 0:
            elapsed = time.time() - t0
            pct = (i + 1) / len(lat_range) * 100
            print(f"  {pct:.0f}% ({i+1}/{len(lat_range)} lat bands) "
                  f"— {elapsed:.1f}s — best so far: ({best_pole[0]:.0f}°, {best_pole[1]:.0f}°) "
                  f"score={best_score:.4f}")

    elapsed = time.time() - t0
    print(f"\nPass 1 complete in {elapsed:.1f}s")
    print(f"Best pole (pass 1): ({best_pole[0]:.1f}°, {best_pole[1]:.1f}°) score={best_score:.4f}")

    # ── Pass 2: River crossings for top 100 poles ─────────────────────────
    print("\nPass 2: Computing river crossings for top 100 poles...")
    all_scores.sort(key=lambda x: x['score_pass1'], reverse=True)
    top100 = all_scores[:100]

    for k, entry in enumerate(top100):
        plat, plon = entry['pole_lat'], entry['pole_lon']
        mean_hab, land_frac, mean_hab_land, lats, lons, land_vals = \
            sample_habitability_along_gc(plat, plon, n_points=360)
        river_x = count_river_crossings(lats, lons, land_vals, threshold_deg=1.0)
        # Final score
        final_score = mean_hab * land_frac * (1 + math.log(river_x + 1))
        entry['river_crossings'] = river_x
        entry['final_score'] = float(final_score)
        if (k + 1) % 20 == 0:
            print(f"  {k+1}/100 done")

    top100.sort(key=lambda x: x['final_score'], reverse=True)

    # ── Compute Alison pole metrics ───────────────────────────────────────
    print("\nComputing Alison pole metrics...")
    a_lat, a_lon = ALISON_POLE
    mean_hab, land_frac, mean_hab_land, lats, lons, land_vals = \
        sample_habitability_along_gc(a_lat, a_lon, n_points=360)
    river_x = count_river_crossings(lats, lons, land_vals, threshold_deg=1.0)
    alison_final = mean_hab * land_frac * (1 + math.log(river_x + 1))

    alison_metrics = {
        'pole_lat': float(a_lat),
        'pole_lon': float(a_lon),
        'mean_habitability': float(mean_hab),
        'land_fraction': float(land_frac),
        'mean_hab_land_only': float(mean_hab_land),
        'river_crossings': river_x,
        'final_score': float(alison_final),
        'score_pass1': float(mean_hab * land_frac)
    }

    # ── Angular distance ──────────────────────────────────────────────────
    opt = top100[0]
    dlat = math.radians(opt['pole_lat'] - a_lat)
    dlon = math.radians(opt['pole_lon'] - a_lon)
    a_r = math.sin(dlat / 2) ** 2 + \
          math.cos(math.radians(a_lat)) * math.cos(math.radians(opt['pole_lat'])) * \
          math.sin(dlon / 2) ** 2
    angular_dist = math.degrees(2 * math.asin(math.sqrt(min(1.0, a_r))))

    # Also compute angular distance between Alison pole and optimal pass1 pole
    # (which may differ from top final-score pole)
    # Find Alison rank in pass1
    all_scores_p1 = sorted(all_scores, key=lambda x: x['score_pass1'], reverse=True)
    alison_rank_p1 = None
    for idx, s in enumerate(all_scores_p1):
        if abs(s['pole_lat'] - a_lat) < 2 and abs(s['pole_lon'] - a_lon) < 2:
            alison_rank_p1 = idx + 1
            break

    # ── Rank of Alison pole among all poles ───────────────────────────────
    # Find pass1 score for Alison pole position
    alison_p1_score = mean_hab * land_frac
    n_better = sum(1 for s in all_scores if s['score_pass1'] > alison_p1_score)
    alison_percentile = (1 - n_better / len(all_scores)) * 100

    # ── Verdict ───────────────────────────────────────────────────────────
    if angular_dist < 10:
        verdict = "HABITABLE — Habitability geometry explains the collinearity"
    elif angular_dist < 25:
        verdict = "PARTIAL — Habitability contributes but doesn't fully explain"
    else:
        verdict = "NOT EXPLAINED — Habitability doesn't explain it, something else is needed"

    # ── Results ───────────────────────────────────────────────────────────
    output = {
        'test': 'Habitability Optimization',
        'method': 'Proxy habitability (land fraction, coast proximity, latitude suitability, elevation)',
        'grid_resolution_deg': 2,
        'n_poles_tested': n_poles,
        'gc_points_per_circle': 360,
        'scoring': 'mean_hab * land_fraction * (1 + log(river_crossings + 1))',
        'optimal_pole': {
            'lat': opt['pole_lat'],
            'lon': opt['pole_lon'],
            'final_score': opt['final_score'],
            'mean_habitability': opt['mean_habitability'],
            'land_fraction': opt['land_fraction'],
            'river_crossings': opt['river_crossings']
        },
        'alison_pole': alison_metrics,
        'alison_percentile_pass1': round(alison_percentile, 2),
        'alison_rank_approx': alison_rank_p1,
        'angular_distance_deg': round(angular_dist, 2),
        'angular_distance_km': round(angular_dist * 111.195, 1),
        'verdict': verdict,
        'top_20_poles': top100[:20]
    }

    out_json = os.path.join(OUT_DIR, 'habitability_optimization.json')
    with open(out_json, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_json}")

    # ── Print summary ─────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"Optimal pole:  ({opt['pole_lat']:.1f}°N, {opt['pole_lon']:.1f}°E)")
    print(f"  Score:         {opt['final_score']:.4f}")
    print(f"  Land fraction: {opt['land_fraction']:.3f}")
    print(f"  Mean hab:      {opt['mean_habitability']:.4f}")
    print(f"  Rivers:        {opt['river_crossings']}")
    print()
    print(f"Alison pole:   ({a_lat}°N, {a_lon}°E)")
    print(f"  Score:         {alison_final:.4f}")
    print(f"  Land fraction: {alison_metrics['land_fraction']:.3f}")
    print(f"  Mean hab:      {alison_metrics['mean_habitability']:.4f}")
    print(f"  Rivers:        {alison_metrics['river_crossings']}")
    print(f"  Percentile:    {alison_percentile:.1f}%")
    print()
    print(f"Angular distance: {angular_dist:.2f}°  ({angular_dist * 111.195:.0f} km)")
    print(f"VERDICT: {verdict}")

    # ── Heatmap ───────────────────────────────────────────────────────────
    print("\nGenerating heatmap...")
    fig, ax = plt.subplots(1, 1, figsize=(14, 7))

    # Plot pass1 scores as heatmap
    extent = [-180, 180, 0, 90]
    im = ax.imshow(results, origin='upper', extent=extent,
                   aspect='auto', cmap='YlOrRd', interpolation='bilinear')
    plt.colorbar(im, ax=ax, label='Habitability Score (pass 1: mean_hab × land_frac)')

    # Mark optimal pole
    ax.plot(opt['pole_lon'], opt['pole_lat'], 'b*', markersize=18,
            markeredgecolor='white', markeredgewidth=1.5,
            label=f"Optimal ({opt['pole_lat']:.0f}°, {opt['pole_lon']:.0f}°)")

    # Mark Alison pole
    ax.plot(a_lon, a_lat, 'g^', markersize=14,
            markeredgecolor='white', markeredgewidth=1.5,
            label=f"Alison ({a_lat:.1f}°, {a_lon:.1f}°)")

    # Mark top 5
    for k in range(min(5, len(top100))):
        if k == 0:
            continue
        ax.plot(top100[k]['pole_lon'], top100[k]['pole_lat'], 'bo',
                markersize=8, markeredgecolor='white', alpha=0.7)

    ax.set_xlabel('Pole Longitude (°)')
    ax.set_ylabel('Pole Latitude (°N)')
    ax.set_title(f'Habitability Optimization: Score by Pole Position\n'
                 f'Angular distance optimal→Alison: {angular_dist:.1f}° — {verdict}')
    ax.legend(loc='lower left', fontsize=9)

    out_png = os.path.join(OUT_DIR, 'habitability_heatmap.png')
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"Heatmap saved to {out_png}")

    return output


if __name__ == '__main__':
    run_optimization()
