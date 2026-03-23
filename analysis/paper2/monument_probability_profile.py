#!/usr/bin/env python3
"""
Directive 3: Monument Probability Profile Along the Circle
==========================================================
At every 1° of arc along the Great Circle, extract geographic features
and use the trained XGBoost model to predict P(monument). This creates
a 360-point profile showing where the ML model predicts monuments are
most likely vs settlements.

Outputs:
  - outputs/extended_analysis/monument_probability_profile.json
  - outputs/extended_analysis/monument_probability_profile.png
"""

import os, sys, json, math, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

warnings.filterwarnings('ignore')
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

# ─── Paths ───
BASE = Path(os.path.expanduser("~/megalith_site_research"))
DATA = BASE / "data"
ML_DATA = BASE / "ml_features"
OUT = BASE / "outputs" / "extended_analysis"
OUT.mkdir(parents=True, exist_ok=True)

# ─── Great Circle Parameters ───
GC_POLE_LAT = 59.682122
GC_POLE_LON = -138.646087
EARTH_R = 6371.0
QUARTER_CIRC = EARTH_R * math.pi / 2

# Memphis (Giza) reference
MEMPHIS_LAT = 29.9792
MEMPHIS_LON = 31.1342


# ════════════════════════════════════════════════════════════════
# GEOMETRY UTILITIES
# ════════════════════════════════════════════════════════════════

def haversine(lat1, lon1, lat2, lon2):
    """Vectorized haversine distance in km."""
    R = EARTH_R
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    return R * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def gc_distance(lats, lons):
    """Distance from each point to the great circle defined by its pole."""
    pole_lat_r = np.radians(GC_POLE_LAT)
    pole_lon_r = np.radians(GC_POLE_LON)
    lats_r = np.radians(lats)
    lons_r = np.radians(lons)
    cos_ang = (np.sin(pole_lat_r) * np.sin(lats_r) +
               np.cos(pole_lat_r) * np.cos(lats_r) * np.cos(lons_r - pole_lon_r))
    cos_ang = np.clip(cos_ang, -1, 1)
    ang = np.arccos(cos_ang)
    gc_ang = np.abs(ang - np.pi/2)
    return gc_ang * EARTH_R


def arc_position(lat, lon):
    """Compute arc position (0-360°) along the great circle."""
    lat_r = math.radians(lat)
    lon_r = math.radians(lon)
    pole_lat_r = math.radians(GC_POLE_LAT)
    pole_lon_r = math.radians(GC_POLE_LON)

    x = math.cos(lat_r) * math.cos(lon_r)
    y = math.cos(lat_r) * math.sin(lon_r)
    z = math.sin(lat_r)

    px = math.cos(pole_lat_r) * math.cos(pole_lon_r)
    py = math.cos(pole_lat_r) * math.sin(pole_lon_r)
    pz = math.sin(pole_lat_r)

    dot = x * px + y * py + z * pz
    proj_x = x - dot * px
    proj_y = y - dot * py
    proj_z = z - dot * pz

    norm = math.sqrt(proj_x**2 + proj_y**2 + proj_z**2)
    if norm < 1e-10:
        return 0.0
    proj_x /= norm; proj_y /= norm; proj_z /= norm

    ref_lon = GC_POLE_LON + 90
    ref_lat_r = math.radians(0)
    ref_lon_r = math.radians(ref_lon)
    rx = math.cos(ref_lat_r) * math.cos(ref_lon_r)
    ry = math.cos(ref_lat_r) * math.sin(ref_lon_r)
    rz = math.sin(ref_lat_r)

    rdot = rx * px + ry * py + rz * pz
    rproj_x = rx - rdot * px
    rproj_y = ry - rdot * py
    rproj_z = rz - rdot * pz
    rnorm = math.sqrt(rproj_x**2 + rproj_y**2 + rproj_z**2)
    if rnorm < 1e-10:
        return 0.0
    rproj_x /= rnorm; rproj_y /= rnorm; rproj_z /= rnorm

    bx = py * rproj_z - pz * rproj_y
    by = pz * rproj_x - px * rproj_z
    bz = px * rproj_y - py * rproj_x

    cos_a = proj_x * rproj_x + proj_y * rproj_y + proj_z * rproj_z
    sin_a = proj_x * bx + proj_y * by + proj_z * bz

    return math.degrees(math.atan2(sin_a, cos_a)) % 360


# ════════════════════════════════════════════════════════════════
# FEATURE EXTRACTION FOR ARC POINTS
# ════════════════════════════════════════════════════════════════

def load_circle_coordinates():
    """Load pre-computed circle coordinates (360 points)."""
    cc_path = BASE / "github-repo" / "data" / "circle_coordinates.json"
    with open(cc_path) as f:
        data = json.load(f)
    return data['points']  # list of {bearing, lat, lon}


def extract_elevation_batch(lats, lons):
    """Extract elevation from ETOPO 2022."""
    import netCDF4
    etopo_path = ML_DATA / "etopo2022.nc"
    if not etopo_path.exists():
        print("  ETOPO not found, returning NaN")
        return np.full(len(lats), np.nan)

    ds = netCDF4.Dataset(str(etopo_path))
    elev_var = ds.variables['z']
    lat_var = ds.variables['lat'][:]
    lon_var = ds.variables['lon'][:]

    elevations = np.zeros(len(lats))
    for i in range(len(lats)):
        lat_idx = np.argmin(np.abs(lat_var - lats[i]))
        lon_idx = np.argmin(np.abs(lon_var - lons[i]))
        elevations[i] = float(elev_var[lat_idx, lon_idx])
    ds.close()
    return elevations


def extract_shapefile_distances(lats, lons, shapefile_path, label="feature"):
    """Compute minimum distance from each point to nearest shapefile geometry."""
    try:
        import geopandas as gpd
        from shapely.geometry import Point
        from shapely.ops import nearest_points, unary_union

        gdf = gpd.read_file(shapefile_path)
        merged = unary_union(gdf.geometry)

        distances = np.zeros(len(lats))
        for i in range(len(lats)):
            pt = Point(lons[i], lats[i])
            nearest = nearest_points(pt, merged)[1]
            distances[i] = haversine(
                np.array([lats[i]]), np.array([lons[i]]),
                np.array([nearest.y]), np.array([nearest.x])
            )[0]
        return distances
    except Exception as e:
        print(f"  Warning: could not compute {label}: {e}")
        return np.full(len(lats), np.nan)


def extract_raster_values(lats, lons, raster_path, label="raster"):
    """Extract values from a GeoTIFF at given coordinates."""
    try:
        import rasterio
        with rasterio.open(raster_path) as src:
            values = np.zeros(len(lats))
            for i in range(len(lats)):
                try:
                    row, col = src.index(lons[i], lats[i])
                    if 0 <= row < src.height and 0 <= col < src.width:
                        val = src.read(1, window=rasterio.windows.Window(col, row, 1, 1))[0, 0]
                        values[i] = float(val) if val != src.nodata else np.nan
                    else:
                        values[i] = np.nan
                except:
                    values[i] = np.nan
            return values
    except Exception as e:
        print(f"  Warning: could not extract {label}: {e}")
        return np.full(len(lats), np.nan)


def compute_arc_features(circle_points, site_lats, site_lons, site_labels):
    """Compute geographic features for 360 arc points (matching ML model features)."""
    lats = np.array([p['lat'] for p in circle_points])
    lons = np.array([p['lon'] for p in circle_points])
    n = len(lats)

    print(f"Computing features for {n} arc points...")
    features = {}

    # Feature 1: dist_to_gc_km — should be ~0 since points are ON the circle
    print("  dist_to_gc_km (should be ~0)")
    features['dist_to_gc_km'] = gc_distance(lats, lons)

    # Feature 2: elevation
    print("  elevation_m")
    features['elevation_m'] = extract_elevation_batch(lats, lons)

    # Feature 3: distance to coastline
    print("  dist_to_coast_km")
    coast_file = ML_DATA / "ne_10m_coastline.shp"
    if not coast_file.exists():
        coast_file = DATA / "natural_earth" / "ne_50m_coastline.shp"
    if coast_file.exists():
        features['dist_to_coast_km'] = extract_shapefile_distances(lats, lons, coast_file, "coastline")
    else:
        features['dist_to_coast_km'] = np.full(n, np.nan)

    # Feature 4: distance to rivers
    print("  dist_to_river_km")
    river_file = ML_DATA / "ne_10m_rivers_lake_centerlines.shp"
    if not river_file.exists():
        river_file = DATA / "natural_earth" / "ne_50m_rivers_lake_centerlines.shp"
    if river_file.exists():
        features['dist_to_river_km'] = extract_shapefile_distances(lats, lons, river_file, "rivers")
    else:
        features['dist_to_river_km'] = np.full(n, np.nan)

    # Feature 5: cloud cover
    print("  cloud_cover")
    cloud_file = ML_DATA / "cloud_mean_annual.tif"
    if cloud_file.exists():
        features['cloud_cover'] = extract_raster_values(lats, lons, cloud_file, "cloud")
    else:
        features['cloud_cover'] = np.full(n, np.nan)

    # Feature 6: seismic hazard
    print("  seismic_pga")
    seismic_file = ML_DATA / "gem_seismic_hazard.tif"
    if seismic_file.exists():
        features['seismic_pga'] = extract_raster_values(lats, lons, seismic_file, "seismic")
    else:
        features['seismic_pga'] = np.full(n, np.nan)

    # Feature 7: absolute latitude
    print("  abs_latitude")
    features['abs_latitude'] = np.abs(lats)

    # Feature 8: longitude
    print("  longitude")
    features['longitude'] = lons.copy()

    # Feature 9: site density (count of known sites within 100km)
    print("  site_density (100km radius)")
    site_lats_arr = np.array(site_lats)
    site_lons_arr = np.array(site_lons)
    density = np.zeros(n, dtype=int)
    for i in range(n):
        dists = haversine(
            np.full(len(site_lats_arr), lats[i]),
            np.full(len(site_lons_arr), lons[i]),
            site_lats_arr, site_lons_arr
        )
        density[i] = np.sum(dists <= 100)
        if (i + 1) % 60 == 0:
            print(f"    {i+1}/360 done")
    features['site_density'] = density

    # Feature 10: cross_type_dist_km (distance to nearest opposite-type site)
    # For arc points, compute distance to nearest monument AND nearest settlement
    # Use the average as a proxy (the model was trained with this feature)
    print("  cross_type_dist_km")
    mon_mask = np.array(site_labels) == 1
    set_mask = np.array(site_labels) == 0
    mon_lats = site_lats_arr[mon_mask]
    mon_lons = site_lons_arr[mon_mask]
    set_lats = site_lats_arr[set_mask]
    set_lons = site_lons_arr[set_mask]

    cross_dist = np.zeros(n)
    for i in range(n):
        # For a neutral point on the circle, use distance to nearest monument
        # (since we're predicting P(monument), the model expects settlement→monument dist
        # for settlements and monument→settlement dist for monuments; we use settlement perspective)
        d_to_mon = haversine(
            np.full(len(mon_lats), lats[i]),
            np.full(len(mon_lons), lons[i]),
            mon_lats, mon_lons
        )
        cross_dist[i] = np.min(d_to_mon) if len(d_to_mon) > 0 else 500.0
    features['cross_type_dist_km'] = cross_dist

    return pd.DataFrame(features)


# ════════════════════════════════════════════════════════════════
# MODEL TRAINING (replicate the ML feature importance pipeline)
# ════════════════════════════════════════════════════════════════

def load_pleiades():
    """Load Pleiades and classify into monuments (1) vs settlements (0)."""
    df = pd.read_csv(DATA / "pleiades" / "pleiades-places-latest.csv")
    df = df.dropna(subset=['reprLat', 'reprLong'])

    monument_types = {'temple', 'temple-2', 'sanctuary', 'pyramid', 'monument',
                      'amphitheatre', 'aqueduct', 'tumulus', 'tomb', 'theatre',
                      'bath', 'church', 'church-2', 'architecturalcomplex'}
    settlement_types = {'settlement', 'villa', 'settlement-modern', 'fort', 'fort-2',
                        'station', 'port', 'farm', 'city', 'village', 'town'}

    def classify(ft_str):
        if pd.isna(ft_str):
            return None
        types = {t.strip() for t in ft_str.split(',')}
        is_mon = bool(types & monument_types)
        is_set = bool(types & settlement_types)
        if is_mon and not is_set:
            return 1
        elif is_set and not is_mon:
            return 0
        return None

    df['label'] = df['featureTypes'].apply(classify)
    df = df.dropna(subset=['label'])
    df['label'] = df['label'].astype(int)
    print(f"Pleiades: {len(df)} sites (monuments={df['label'].sum()}, settlements={(df['label']==0).sum()})")
    return df


def train_model():
    """Train the XGBoost model using pre-computed feature matrix."""
    import xgboost as xgb

    feature_csv = BASE / "outputs" / "ml_feature_importance" / "feature_matrix.csv"
    print(f"Loading feature matrix from {feature_csv}...")
    df = pd.read_csv(feature_csv)

    # The feature columns used in the model
    feature_cols = ['dist_to_gc_km', 'elevation_m', 'dist_to_coast_km', 'dist_to_river_km',
                    'cloud_cover', 'seismic_pga', 'abs_latitude', 'longitude',
                    'site_density', 'cross_type_dist_km']

    # Filter to rows that have all features
    available_cols = [c for c in feature_cols if c in df.columns]
    print(f"  Available features: {available_cols}")

    X = df[available_cols].copy()
    y = df['label'].values

    # Impute missing values with median (same as original script)
    for col in X.columns:
        n_miss = X[col].isna().sum()
        if n_miss > 0:
            X[col].fillna(X[col].median(), inplace=True)
            print(f"  Imputed {col}: {n_miss} values")

    print(f"  Training XGBoost on {len(X)} sites with {len(available_cols)} features...")
    model = xgb.XGBClassifier(
        n_estimators=500, max_depth=6, learning_rate=0.05,
        random_state=42, eval_metric='logloss', verbosity=0
    )
    model.fit(X.values, y)
    print("  Model trained.")

    return model, available_cols, df


# ════════════════════════════════════════════════════════════════
# PLOTTING
# ════════════════════════════════════════════════════════════════

def plot_profile(arc_data, monument_arcs, settlement_arcs, memphis_arc, output_path):
    """Create the monument probability profile plot."""
    bearings = [d['bearing'] for d in arc_data]
    probs = [d['p_monument'] for d in arc_data]

    fig, ax = plt.subplots(figsize=(18, 7))

    # Main probability curve
    ax.fill_between(bearings, probs, alpha=0.3, color='crimson', label='P(monument)')
    ax.plot(bearings, probs, 'crimson', linewidth=1.5)

    # Actual monument positions (as rug ticks along top)
    if monument_arcs:
        for arc in monument_arcs:
            ax.axvline(x=arc, ymin=0.92, ymax=1.0, color='red', alpha=0.15, linewidth=0.5)
        # Add a few as visible markers
        ax.scatter(monument_arcs, [1.02] * len(monument_arcs),
                   marker='v', s=3, color='red', alpha=0.3, zorder=5,
                   clip_on=False, label=f'Monuments on circle (n={len(monument_arcs)})')

    # Settlement positions
    if settlement_arcs:
        for arc in settlement_arcs:
            ax.axvline(x=arc, ymin=0.0, ymax=0.08, color='blue', alpha=0.08, linewidth=0.5)
        ax.scatter(settlement_arcs, [-0.02] * len(settlement_arcs),
                   marker='^', s=3, color='blue', alpha=0.2, zorder=5,
                   clip_on=False, label=f'Settlements on circle (n={len(settlement_arcs)})')

    # Memphis marker
    if memphis_arc is not None:
        ax.axvline(x=memphis_arc, color='gold', linewidth=2.5, linestyle='--',
                   zorder=10, label=f'Memphis/Giza ({memphis_arc:.1f}°)')

    # Peak marker
    peak_idx = np.argmax(probs)
    peak_bearing = bearings[peak_idx]
    peak_prob = probs[peak_idx]
    ax.annotate(f'Peak: {peak_bearing}° (P={peak_prob:.3f})',
                xy=(peak_bearing, peak_prob),
                xytext=(peak_bearing + 15, peak_prob + 0.03),
                arrowprops=dict(arrowstyle='->', color='black'),
                fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.8))

    # Annotate notable regions
    ax.set_xlabel('Arc Position (degrees)', fontsize=13)
    ax.set_ylabel('P(monument)', fontsize=13)
    ax.set_title('Monument Probability Profile Along the Great Circle\n'
                 'XGBoost model predicting P(monument) from geographic features at each 1° arc position',
                 fontsize=14)
    ax.set_xlim(0, 360)
    ax.set_ylim(-0.05, max(probs) + 0.1)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.2)

    # Add geographic annotations along x-axis
    # Identify which regions correspond to which part of the world
    # (from circle_coordinates.json: bearing 0 is ~30.3°N, 41.4°E = Iraq)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved plot: {output_path}")


# ════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("DIRECTIVE 3: Monument Probability Profile Along the Circle")
    print("=" * 70)

    # Step 1: Load circle coordinates
    print("\n[1] Loading circle coordinates...")
    circle_points = load_circle_coordinates()
    print(f"  {len(circle_points)} arc points loaded")

    # Step 2: Train the model using pre-computed feature matrix
    print("\n[2] Training XGBoost model...")
    model, feature_cols, site_df = train_model()

    # Step 3: Extract features for each arc point
    print("\n[3] Extracting geographic features at each arc point...")
    site_lats = site_df['reprLat'].values
    site_lons = site_df['reprLong'].values
    site_labels = site_df['label'].values

    arc_features = compute_arc_features(circle_points, site_lats, site_lons, site_labels)

    # Ensure column order matches the model
    arc_features = arc_features[feature_cols]

    # Impute any NaN with training medians
    train_matrix = pd.read_csv(BASE / "outputs" / "ml_feature_importance" / "feature_matrix.csv")
    for col in feature_cols:
        n_miss = arc_features[col].isna().sum()
        if n_miss > 0:
            median_val = train_matrix[col].median()
            arc_features[col].fillna(median_val, inplace=True)
            print(f"  Imputed {col} for arc points: {n_miss} values (median={median_val:.2f})")

    # Step 4: Predict P(monument) at each arc position
    print("\n[4] Predicting P(monument) at each arc position...")
    probs = model.predict_proba(arc_features.values)[:, 1]
    print(f"  P(monument) range: [{probs.min():.4f}, {probs.max():.4f}]")
    print(f"  Mean P(monument): {probs.mean():.4f}")

    # Step 5: Get actual monument/settlement positions along the circle
    print("\n[5] Computing arc positions for actual sites near the circle...")
    # Filter to sites within 50km of the circle
    site_gc_dists = gc_distance(site_lats, site_lons)
    near_mask = site_gc_dists <= 50  # within 50km corridor

    monument_arcs = []
    settlement_arcs = []
    near_indices = np.where(near_mask)[0]
    print(f"  {len(near_indices)} sites within 50km of circle")

    for idx in near_indices:
        arc = arc_position(site_lats[idx], site_lons[idx])
        if site_labels[idx] == 1:
            monument_arcs.append(arc)
        else:
            settlement_arcs.append(arc)

    print(f"  Monuments on circle: {len(monument_arcs)}")
    print(f"  Settlements on circle: {len(settlement_arcs)}")

    # Step 6: Memphis arc position
    memphis_arc = arc_position(MEMPHIS_LAT, MEMPHIS_LON)
    memphis_gc_dist = gc_distance(np.array([MEMPHIS_LAT]), np.array([MEMPHIS_LON]))[0]
    print(f"\n  Memphis arc position: {memphis_arc:.2f}°")
    print(f"  Memphis distance from circle: {memphis_gc_dist:.1f} km")

    # Step 7: Find the peak
    peak_idx = np.argmax(probs)
    peak_bearing = circle_points[peak_idx]['bearing']
    peak_lat = circle_points[peak_idx]['lat']
    peak_lon = circle_points[peak_idx]['lon']
    peak_prob = probs[peak_idx]

    print(f"\n[6] P(monument) PEAK:")
    print(f"  Bearing: {peak_bearing}°")
    print(f"  Location: ({peak_lat:.4f}, {peak_lon:.4f})")
    print(f"  P(monument): {peak_prob:.4f}")

    # Distance from peak to Memphis
    peak_to_memphis_km = haversine(
        np.array([peak_lat]), np.array([peak_lon]),
        np.array([MEMPHIS_LAT]), np.array([MEMPHIS_LON])
    )[0]
    memphis_bearing_delta = abs(peak_bearing - memphis_arc)
    if memphis_bearing_delta > 180:
        memphis_bearing_delta = 360 - memphis_bearing_delta

    print(f"  Distance from peak to Memphis: {peak_to_memphis_km:.1f} km")
    print(f"  Bearing difference from Memphis: {memphis_bearing_delta:.1f}°")
    memphis_is_peak = memphis_bearing_delta <= 5
    print(f"  Memphis corresponds to peak? {'YES' if memphis_is_peak else 'NO'}")

    # P(monument) at Memphis's arc position
    memphis_arc_idx = int(round(memphis_arc)) % 360
    p_at_memphis = probs[memphis_arc_idx]
    print(f"  P(monument) at Memphis arc position ({memphis_arc_idx}°): {p_at_memphis:.4f}")

    # Rank of Memphis position among all 360 points
    rank = int(np.sum(probs >= p_at_memphis))
    print(f"  Memphis rank: #{rank} out of 360")

    # Step 8: Build output data
    arc_data = []
    for i, pt in enumerate(circle_points):
        arc_data.append({
            'bearing': pt['bearing'],
            'lat': pt['lat'],
            'lon': pt['lon'],
            'p_monument': float(probs[i]),
            'features': {col: float(arc_features.iloc[i][col]) for col in feature_cols}
        })

    # Top 10 peaks
    top_indices = np.argsort(probs)[::-1][:10]
    top_peaks = []
    for idx in top_indices:
        top_peaks.append({
            'bearing': circle_points[idx]['bearing'],
            'lat': circle_points[idx]['lat'],
            'lon': circle_points[idx]['lon'],
            'p_monument': float(probs[idx])
        })

    results = {
        'metadata': {
            'description': 'Monument probability profile along the Great Circle',
            'model': 'XGBoost (n_estimators=500, max_depth=6, lr=0.05)',
            'features': feature_cols,
            'n_arc_points': 360,
            'corridor_width_km': 50,
            'monuments_on_circle': len(monument_arcs),
            'settlements_on_circle': len(settlement_arcs)
        },
        'statistics': {
            'p_monument_mean': float(probs.mean()),
            'p_monument_std': float(probs.std()),
            'p_monument_min': float(probs.min()),
            'p_monument_max': float(probs.max()),
            'p_monument_median': float(np.median(probs))
        },
        'peak': {
            'bearing': peak_bearing,
            'lat': peak_lat,
            'lon': peak_lon,
            'p_monument': float(peak_prob)
        },
        'memphis': {
            'arc_position': float(memphis_arc),
            'p_monument_at_memphis': float(p_at_memphis),
            'rank_out_of_360': rank,
            'distance_from_peak_km': float(peak_to_memphis_km),
            'bearing_diff_from_peak': float(memphis_bearing_delta),
            'corresponds_to_peak': memphis_is_peak,
            'gc_distance_km': float(memphis_gc_dist)
        },
        'top_10_peaks': top_peaks,
        'profile': arc_data
    }

    # Save JSON
    json_path = OUT / "monument_probability_profile.json"
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved JSON: {json_path}")

    # Step 9: Plot
    print("\n[7] Generating plot...")
    plot_profile(arc_data, monument_arcs, settlement_arcs, memphis_arc,
                 OUT / "monument_probability_profile.png")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Peak P(monument): {peak_prob:.4f} at bearing {peak_bearing}° ({peak_lat:.2f}°N, {peak_lon:.2f}°E)")
    print(f"Memphis P(monument): {p_at_memphis:.4f} at bearing {memphis_arc:.1f}° (rank #{rank}/360)")
    if memphis_is_peak:
        print(">>> Memphis IS at the geographic prediction peak")
    else:
        print(f">>> Memphis is NOT at the peak (off by {memphis_bearing_delta:.1f}°, {peak_to_memphis_km:.0f} km)")
    print("=" * 70)


if __name__ == '__main__':
    main()
