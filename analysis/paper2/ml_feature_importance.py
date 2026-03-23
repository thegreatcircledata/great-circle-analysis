#!/usr/bin/env python3
"""
ML Feature Importance Test: Does "Distance to Great Circle" predict monument placement?

Tests whether distance to the Great Circle has predictive power for distinguishing
monuments from settlements after controlling for geographic, geological, climatic,
and hydrological variables.
"""

import os
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from collections import Counter

warnings.filterwarnings('ignore')

# ─── Paths ───
BASE = Path(os.path.expanduser("~/megalith_site_research"))
DATA = BASE / "data"
ML_DATA = BASE / "ml_features"
OUT = BASE / "outputs" / "ml_feature_importance"
OUT.mkdir(parents=True, exist_ok=True)

# ─── Great Circle Parameters ───
GC_POLE_LAT = 59.682122
GC_POLE_LON = -138.646087


# ════════════════════════════════════════════════════════════════
# STEP 1: Load and classify Pleiades sites
# ════════════════════════════════════════════════════════════════

def load_pleiades():
    """Load Pleiades and classify into monuments (1) vs settlements (0)."""
    df = pd.read_csv(DATA / "pleiades" / "pleiades-places-latest.csv")

    # Need valid coordinates
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

    print(f"Pleiades sites: {len(df)} (monuments={df['label'].sum()}, settlements={(df['label']==0).sum()})")
    return df[['reprLat', 'reprLong', 'label', 'featureTypes', 'minDate', 'maxDate', 'title']].copy()


# ════════════════════════════════════════════════════════════════
# STEP 2: Compute Features
# ════════════════════════════════════════════════════════════════

def haversine(lat1, lon1, lat2, lon2):
    """Vectorized haversine distance in km."""
    R = 6371.0
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    return R * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def gc_distance(lats, lons, pole_lat=GC_POLE_LAT, pole_lon=GC_POLE_LON):
    """Distance from each point to the great circle defined by its pole."""
    pole_lat_r = np.radians(pole_lat)
    pole_lon_r = np.radians(pole_lon)
    lats_r = np.radians(lats)
    lons_r = np.radians(lons)

    # Angular distance from pole
    cos_ang = (np.sin(pole_lat_r) * np.sin(lats_r) +
               np.cos(pole_lat_r) * np.cos(lats_r) * np.cos(lons_r - pole_lon_r))
    cos_ang = np.clip(cos_ang, -1, 1)
    ang = np.arccos(cos_ang)

    # Distance to great circle = |angular distance to pole - 90°|
    gc_ang = np.abs(ang - np.pi/2)
    return gc_ang * 6371.0


def compute_shapefile_distances(lats, lons, shapefile_path, label="feature"):
    """Compute minimum distance from each site to nearest shapefile geometry."""
    try:
        import geopandas as gpd
        from shapely.geometry import Point
        from shapely.ops import nearest_points

        gdf = gpd.read_file(shapefile_path)
        # Merge all geometries into a single multi-geometry for faster querying
        from shapely.ops import unary_union
        merged = unary_union(gdf.geometry)

        distances = np.zeros(len(lats))
        for i in range(len(lats)):
            pt = Point(lons.iloc[i] if hasattr(lons, 'iloc') else lons[i],
                       lats.iloc[i] if hasattr(lats, 'iloc') else lats[i])
            nearest = nearest_points(pt, merged)[1]
            distances[i] = haversine(
                lats.iloc[i] if hasattr(lats, 'iloc') else lats[i],
                lons.iloc[i] if hasattr(lons, 'iloc') else lons[i],
                nearest.y, nearest.x
            )
            if (i+1) % 2000 == 0:
                print(f"  {label}: {i+1}/{len(lats)} sites processed")

        return distances
    except Exception as e:
        print(f"  Warning: could not compute {label} distances: {e}")
        return np.full(len(lats), np.nan)


def extract_raster_values(lats, lons, raster_path, label="raster"):
    """Extract values from a GeoTIFF at given coordinates."""
    try:
        import rasterio
        with rasterio.open(raster_path) as src:
            values = np.zeros(len(lats))
            for i in range(len(lats)):
                lat = lats.iloc[i] if hasattr(lats, 'iloc') else lats[i]
                lon = lons.iloc[i] if hasattr(lons, 'iloc') else lons[i]
                try:
                    row, col = src.index(lon, lat)
                    if 0 <= row < src.height and 0 <= col < src.width:
                        val = src.read(1, window=rasterio.windows.Window(col, row, 1, 1))[0, 0]
                        if val == src.nodata:
                            values[i] = np.nan
                        else:
                            values[i] = float(val)
                    else:
                        values[i] = np.nan
                except:
                    values[i] = np.nan

            valid = np.sum(~np.isnan(values))
            print(f"  {label}: {valid}/{len(values)} valid extractions")
            return values
    except Exception as e:
        print(f"  Warning: could not extract {label}: {e}")
        return np.full(len(lats), np.nan)


def compute_elevation_etopo(lats, lons):
    """Extract elevation from ETOPO 2022 NetCDF."""
    import netCDF4

    etopo_path = ML_DATA / "etopo2022.nc"
    if not etopo_path.exists():
        print("  ETOPO 2022 not found, trying SRTM cache")
        return compute_elevation_srtm_fallback(lats, lons)

    ds = netCDF4.Dataset(str(etopo_path))
    elev_var = ds.variables['z']
    lat_var = ds.variables['lat'][:]
    lon_var = ds.variables['lon'][:]

    elevations = np.zeros(len(lats))
    for i in range(len(lats)):
        lat = lats.iloc[i] if hasattr(lats, 'iloc') else lats[i]
        lon = lons.iloc[i] if hasattr(lons, 'iloc') else lons[i]

        lat_idx = np.argmin(np.abs(lat_var - lat))
        lon_idx = np.argmin(np.abs(lon_var - lon))
        elevations[i] = float(elev_var[lat_idx, lon_idx])

    ds.close()
    valid = np.sum(~np.isnan(elevations))
    print(f"  Elevation (ETOPO 2022): {valid}/{len(elevations)} valid")
    return elevations


def compute_elevation_srtm_fallback(lats, lons):
    """Fallback: extract elevation from cached SRTM .hgt files."""
    import struct

    srtm_dir = DATA / "srtm_cache"
    elevations = np.full(len(lats), np.nan)

    if not srtm_dir.exists():
        print("  No SRTM cache found, skipping elevation")
        return elevations

    available = {f.stem: f for f in srtm_dir.glob("*.hgt")}

    for i in range(len(lats)):
        lat = lats.iloc[i] if hasattr(lats, 'iloc') else lats[i]
        lon = lons.iloc[i] if hasattr(lons, 'iloc') else lons[i]

        lat_floor = int(np.floor(lat))
        lon_floor = int(np.floor(lon))
        ns = 'N' if lat_floor >= 0 else 'S'
        ew = 'E' if lon_floor >= 0 else 'W'
        tile_name = f"{ns}{abs(lat_floor):02d}{ew}{abs(lon_floor):03d}"

        if tile_name in available:
            try:
                hgt_file = available[tile_name]
                file_size = hgt_file.stat().st_size
                if file_size == 2884802:
                    samples = 1201
                elif file_size == 25934402:
                    samples = 3601
                else:
                    continue

                frac_lat = lat - lat_floor
                frac_lon = lon - lon_floor
                row = int((1 - frac_lat) * (samples - 1))
                col = int(frac_lon * (samples - 1))

                with open(hgt_file, 'rb') as f:
                    f.seek((row * samples + col) * 2)
                    buf = f.read(2)
                    if len(buf) == 2:
                        elev = struct.unpack('>h', buf)[0]
                        if elev != -32768:
                            elevations[i] = elev
            except:
                pass

    valid = np.sum(~np.isnan(elevations))
    print(f"  Elevation (SRTM fallback): {valid}/{len(elevations)} valid")
    return elevations


def compute_site_density(lats, lons, radius_km=100):
    """Count number of sites within radius_km of each site."""
    n = len(lats)
    density = np.zeros(n, dtype=int)

    lats_arr = np.array(lats)
    lons_arr = np.array(lons)

    # Batch compute for efficiency
    for i in range(n):
        dists = haversine(lats_arr[i], lons_arr[i], lats_arr, lons_arr)
        density[i] = np.sum(dists <= radius_km) - 1  # exclude self
        if (i+1) % 3000 == 0:
            print(f"  Site density: {i+1}/{n}")

    return density


def compute_cross_type_distance(lats, lons, labels):
    """For each monument: distance to nearest settlement; vice versa."""
    lats_arr = np.array(lats)
    lons_arr = np.array(lons)
    labels_arr = np.array(labels)

    mon_mask = labels_arr == 1
    set_mask = labels_arr == 0

    mon_lats, mon_lons = lats_arr[mon_mask], lons_arr[mon_mask]
    set_lats, set_lons = lats_arr[set_mask], lons_arr[set_mask]

    distances = np.zeros(len(lats_arr))

    # For monuments: find nearest settlement
    for i in np.where(mon_mask)[0]:
        dists = haversine(lats_arr[i], lons_arr[i], set_lats, set_lons)
        distances[i] = np.min(dists) if len(dists) > 0 else np.nan

    # For settlements: find nearest monument
    for i in np.where(set_mask)[0]:
        dists = haversine(lats_arr[i], lons_arr[i], mon_lats, mon_lons)
        distances[i] = np.min(dists) if len(dists) > 0 else np.nan

    print(f"  Cross-type distance: done")
    return distances


def compute_all_features(df):
    """Compute all 12 features for the site dataframe."""
    lats = df['reprLat']
    lons = df['reprLong']
    labels = df['label']
    n = len(df)

    print(f"\nComputing features for {n} sites...")
    features = {}

    # Feature 1: Distance to Great Circle
    print("Feature 1: Distance to Great Circle")
    features['dist_to_gc_km'] = gc_distance(lats.values, lons.values)

    # Feature 2: Elevation
    print("Feature 2: Elevation")
    features['elevation_m'] = compute_elevation_etopo(lats, lons)

    # Feature 3: Distance to coastline
    print("Feature 3: Distance to coastline")
    coast_10m = ML_DATA / "ne_10m_coastline.shp"
    coast_50m = DATA / "natural_earth" / "ne_50m_coastline.shp"
    coast_file = coast_10m if coast_10m.exists() else coast_50m
    if coast_file.exists():
        features['dist_to_coast_km'] = compute_shapefile_distances(lats, lons, coast_file, "coastline")
    else:
        print("  No coastline shapefile found, using NaN")
        features['dist_to_coast_km'] = np.full(n, np.nan)

    # Feature 4: Distance to rivers
    print("Feature 4: Distance to rivers")
    river_10m = ML_DATA / "ne_10m_rivers_lake_centerlines.shp"
    river_50m = DATA / "natural_earth" / "ne_50m_rivers_lake_centerlines.shp"
    river_file = river_10m if river_10m.exists() else river_50m
    if river_file.exists():
        features['dist_to_river_km'] = compute_shapefile_distances(lats, lons, river_file, "rivers")
    else:
        print("  No rivers shapefile found, using NaN")
        features['dist_to_river_km'] = np.full(n, np.nan)

    # Feature 5: Soil productivity — use as proxy if available, else skip
    print("Feature 5: Soil productivity")
    # We don't have GAEZ or SoilGrids locally — leave as NaN for now
    features['soil_productivity'] = np.full(n, np.nan)
    print("  Soil data not available locally, will impute")

    # Feature 6: Water table depth
    print("Feature 6: Water table depth")
    features['water_table_m'] = np.full(n, np.nan)
    print("  Water table data not available locally, will impute")

    # Feature 7: Cloud cover
    print("Feature 7: Cloud cover")
    cloud_file = ML_DATA / "cloud_mean_annual.tif"
    if cloud_file.exists():
        features['cloud_cover'] = extract_raster_values(lats, lons, cloud_file, "cloud cover")
    else:
        print("  Cloud cover data not available, will impute")
        features['cloud_cover'] = np.full(n, np.nan)

    # Feature 8: Seismic hazard
    print("Feature 8: Seismic hazard")
    seismic_file = ML_DATA / "gem_seismic_hazard.tif"
    if seismic_file.exists():
        features['seismic_pga'] = extract_raster_values(lats, lons, seismic_file, "seismic PGA")
    else:
        print("  Seismic data not available, will impute")
        features['seismic_pga'] = np.full(n, np.nan)

    # Feature 9: Absolute latitude
    print("Feature 9: Absolute latitude")
    features['abs_latitude'] = np.abs(lats.values)

    # Feature 10: Longitude
    print("Feature 10: Longitude")
    features['longitude'] = lons.values.copy()

    # Feature 11: Site density
    print("Feature 11: Site density (100km radius)")
    features['site_density'] = compute_site_density(lats, lons, radius_km=100)

    # Feature 12: Cross-type distance
    print("Feature 12: Cross-type distance")
    features['cross_type_dist_km'] = compute_cross_type_distance(lats, lons, labels)

    return pd.DataFrame(features, index=df.index)


# ════════════════════════════════════════════════════════════════
# STEP 3: Imputation & Preparation
# ════════════════════════════════════════════════════════════════

def prepare_feature_matrix(feature_df):
    """Handle missing values and prepare final feature matrix."""
    X = feature_df.copy()

    imputation_log = {}
    for col in X.columns:
        n_missing = X[col].isna().sum()
        if n_missing > 0:
            pct = 100 * n_missing / len(X)
            if pct > 90:
                # Too much missing — drop this feature
                imputation_log[col] = f"DROPPED ({pct:.0f}% missing)"
                X.drop(columns=[col], inplace=True)
            else:
                median_val = X[col].median()
                X[col].fillna(median_val, inplace=True)
                imputation_log[col] = f"median imputed ({n_missing} values, {pct:.1f}%)"

    print("\nImputation summary:")
    for col, action in imputation_log.items():
        print(f"  {col}: {action}")

    return X, imputation_log


# ════════════════════════════════════════════════════════════════
# STEP 4: Models & Analyses
# ════════════════════════════════════════════════════════════════

def run_models(X, y, feature_names, tag="global"):
    """Train RF, XGBoost, Logistic Regression with cross-validation."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    import xgboost as xgb

    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    results = {}

    # Random Forest
    print(f"\n[{tag}] Training Random Forest...")
    rf = RandomForestClassifier(n_estimators=500, max_depth=10, random_state=42, n_jobs=-1)
    rf_scores = cross_val_score(rf, X, y, cv=cv, scoring='roc_auc')
    rf.fit(X, y)
    results['random_forest'] = {
        'cv_auc': float(np.mean(rf_scores)),
        'cv_std': float(np.std(rf_scores)),
        'feature_importances': dict(zip(feature_names, rf.feature_importances_.tolist()))
    }
    print(f"  RF AUC: {np.mean(rf_scores):.4f} ± {np.std(rf_scores):.4f}")

    # XGBoost
    print(f"[{tag}] Training XGBoost...")
    xgb_model = xgb.XGBClassifier(
        n_estimators=500, max_depth=6, learning_rate=0.05,
        random_state=42, eval_metric='logloss', verbosity=0
    )
    xgb_scores = cross_val_score(xgb_model, X, y, cv=cv, scoring='roc_auc')
    xgb_model.fit(X, y)
    results['xgboost'] = {
        'cv_auc': float(np.mean(xgb_scores)),
        'cv_std': float(np.std(xgb_scores)),
        'feature_importances': dict(zip(feature_names, xgb_model.feature_importances_.tolist()))
    }
    print(f"  XGB AUC: {np.mean(xgb_scores):.4f} ± {np.std(xgb_scores):.4f}")

    # Logistic Regression
    print(f"[{tag}] Training Logistic Regression...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    lr = LogisticRegression(max_iter=1000, random_state=42)
    lr_scores = cross_val_score(lr, X_scaled, y, cv=cv, scoring='roc_auc')
    lr.fit(X_scaled, y)
    results['logistic_regression'] = {
        'cv_auc': float(np.mean(lr_scores)),
        'cv_std': float(np.std(lr_scores)),
        'coefficients': dict(zip(feature_names, lr.coef_[0].tolist()))
    }
    print(f"  LR AUC: {np.mean(lr_scores):.4f} ± {np.std(lr_scores):.4f}")

    return results, rf, xgb_model, lr, scaler


def run_shap_analysis(model, X, feature_names, output_path, title="SHAP Feature Importance"):
    """Run SHAP analysis on XGBoost model."""
    try:
        import shap

        # Use subset for speed if large
        n_samples = min(2000, len(X))
        X_sample = X[:n_samples] if isinstance(X, np.ndarray) else X.iloc[:n_samples]

        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X_sample)

        # Mean absolute SHAP values
        mean_abs_shap = np.mean(np.abs(shap_values), axis=0)
        shap_dict = dict(zip(feature_names, mean_abs_shap.tolist()))

        # Rank GC distance
        ranked = sorted(shap_dict.items(), key=lambda x: x[1], reverse=True)
        gc_rank = next(i+1 for i, (name, _) in enumerate(ranked) if name == 'dist_to_gc_km')

        print(f"  SHAP ranking (GC distance = #{gc_rank} of {len(feature_names)}):")
        for i, (name, val) in enumerate(ranked):
            marker = " <<<" if name == 'dist_to_gc_km' else ""
            print(f"    #{i+1}: {name} = {val:.4f}{marker}")

        # Plot
        fig, ax = plt.subplots(figsize=(10, 8))
        shap.summary_plot(shap_values, X_sample, feature_names=feature_names, show=False)
        plt.title(title, fontsize=14)
        plt.tight_layout()
        plt.savefig(output_path, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_path}")

        return shap_dict, gc_rank, shap_values

    except Exception as e:
        print(f"  SHAP failed ({e}), falling back to permutation importance")
        from sklearn.inspection import permutation_importance
        from sklearn.model_selection import train_test_split

        X_train, X_test, y_train, y_test = train_test_split(
            X, model.predict(X), test_size=0.3, random_state=42
        )
        perm = permutation_importance(model, X, model.predict(X), n_repeats=30, random_state=42)
        shap_dict = dict(zip(feature_names, perm.importances_mean.tolist()))
        ranked = sorted(shap_dict.items(), key=lambda x: x[1], reverse=True)
        gc_rank = next(i+1 for i, (name, _) in enumerate(ranked) if name == 'dist_to_gc_km')
        return shap_dict, gc_rank, None


def ablation_test(X, y, feature_names, gc_col='dist_to_gc_km'):
    """Compare AUC with and without GC distance."""
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    import xgboost as xgb

    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

    # With GC distance
    model_with = xgb.XGBClassifier(n_estimators=500, max_depth=6, learning_rate=0.05,
                                    random_state=42, eval_metric='logloss', verbosity=0)
    auc_with = cross_val_score(model_with, X, y, cv=cv, scoring='roc_auc')

    # Without GC distance
    gc_idx = feature_names.index(gc_col)
    X_without = np.delete(X if isinstance(X, np.ndarray) else X.values, gc_idx, axis=1)
    model_without = xgb.XGBClassifier(n_estimators=500, max_depth=6, learning_rate=0.05,
                                       random_state=42, eval_metric='logloss', verbosity=0)
    auc_without = cross_val_score(model_without, X_without, y, cv=cv, scoring='roc_auc')

    delta = float(np.mean(auc_with) - np.mean(auc_without))

    result = {
        'auc_with_gc': float(np.mean(auc_with)),
        'auc_with_gc_std': float(np.std(auc_with)),
        'auc_without_gc': float(np.mean(auc_without)),
        'auc_without_gc_std': float(np.std(auc_without)),
        'delta_auc': delta,
        'interpretation': (
            'SIGNIFICANT: GC adds information beyond geography' if delta > 0.02
            else 'MODEST: GC adds some information' if delta > 0.005
            else 'MINIMAL: GC is largely a geographic proxy'
        )
    }

    print(f"\n  Ablation Test:")
    print(f"    AUC with GC:    {result['auc_with_gc']:.4f} ± {result['auc_with_gc_std']:.4f}")
    print(f"    AUC without GC: {result['auc_without_gc']:.4f} ± {result['auc_without_gc_std']:.4f}")
    print(f"    Delta AUC:      {result['delta_auc']:.4f}")
    print(f"    → {result['interpretation']}")

    return result


def partial_dependence_plot(model, X, feature_names, gc_col='dist_to_gc_km', output_path=None):
    """Generate partial dependence plot for GC distance."""
    from sklearn.inspection import partial_dependence

    gc_idx = feature_names.index(gc_col)

    pdp = partial_dependence(model, X, features=[gc_idx], kind='average', grid_resolution=100)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(pdp['grid_values'][0], pdp['average'][0], 'b-', linewidth=2)
    ax.axvline(x=25, color='r', linestyle='--', alpha=0.7, label='25 km band')
    ax.axvline(x=50, color='orange', linestyle='--', alpha=0.7, label='50 km band')
    ax.set_xlabel('Distance to Great Circle (km)', fontsize=12)
    ax.set_ylabel('Partial Dependence (monument probability)', fontsize=12)
    ax.set_title('Partial Dependence: Monument Probability vs GC Distance', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved PDP: {output_path}")

    return {
        'gc_distances': pdp['grid_values'][0].tolist(),
        'pdp_values': pdp['average'][0].tolist()
    }


def feature_interaction_test(model, X, feature_names, output_path=None):
    """Test interaction between GC distance and elevation/river distance."""
    try:
        import shap

        n_samples = min(1000, len(X))
        X_sample = X[:n_samples] if isinstance(X, np.ndarray) else X.iloc[:n_samples]

        explainer = shap.TreeExplainer(model)
        shap_interaction = explainer.shap_interaction_values(X_sample)

        gc_idx = feature_names.index('dist_to_gc_km')

        interactions = {}
        for j, fname in enumerate(feature_names):
            if fname != 'dist_to_gc_km':
                interaction_strength = float(np.mean(np.abs(shap_interaction[:, gc_idx, j])))
                interactions[fname] = interaction_strength

        # Sort by strength
        sorted_interactions = sorted(interactions.items(), key=lambda x: x[1], reverse=True)

        print("\n  GC Distance Interactions:")
        for name, val in sorted_interactions:
            print(f"    GC × {name}: {val:.4f}")

        # Plot top interactions
        fig, ax = plt.subplots(figsize=(10, 6))
        names = [x[0] for x in sorted_interactions]
        vals = [x[1] for x in sorted_interactions]
        ax.barh(names[::-1], vals[::-1], color='steelblue')
        ax.set_xlabel('Mean |SHAP interaction value|', fontsize=12)
        ax.set_title('Feature Interactions with Great Circle Distance', fontsize=14)
        plt.tight_layout()
        plt.savefig(output_path, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  Saved interaction plot: {output_path}")

        return dict(sorted_interactions)

    except Exception as e:
        print(f"  Interaction analysis failed: {e}")
        return {}


# ════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════

def main():
    # Load data
    df = load_pleiades()

    # Compute features
    feature_df = compute_all_features(df)

    # Prepare feature matrix
    X, imputation_log = prepare_feature_matrix(feature_df)
    y = df.loc[X.index, 'label'].values
    feature_names = list(X.columns)

    # Save feature matrix
    full_df = df.loc[X.index].copy()
    for col in X.columns:
        full_df[col] = X[col].values
    full_df.to_csv(OUT / "feature_matrix.csv", index=False)
    print(f"\nSaved feature matrix: {OUT / 'feature_matrix.csv'}")

    X_arr = X.values

    # ─── Analysis 1: Full Global Model ───
    print("\n" + "="*60)
    print("ANALYSIS 1: Full Global Model")
    print("="*60)
    results, rf, xgb_model, lr, scaler = run_models(X_arr, y, feature_names, tag="global")

    # SHAP analysis
    print("\nSHAP Analysis (Global)...")
    shap_dict, gc_rank_global, shap_values = run_shap_analysis(
        xgb_model, X_arr, feature_names,
        OUT / "shap_summary.png",
        "SHAP Feature Importance — Global Model"
    )
    results['shap'] = {
        'mean_abs_shap': shap_dict,
        'gc_rank': gc_rank_global,
        'n_features': len(feature_names)
    }

    # ─── Analysis 2: Egypt-Levant-Iran Corridor ───
    print("\n" + "="*60)
    print("ANALYSIS 2: Egypt-Levant-Iran Corridor (20-38°N, 25-60°E)")
    print("="*60)
    corridor_mask = (
        (df.loc[X.index, 'reprLat'] >= 20) & (df.loc[X.index, 'reprLat'] <= 38) &
        (df.loc[X.index, 'reprLong'] >= 25) & (df.loc[X.index, 'reprLong'] <= 60)
    )
    X_corridor = X_arr[corridor_mask]
    y_corridor = y[corridor_mask]
    n_corridor = len(y_corridor)
    print(f"  Corridor sites: {n_corridor} (monuments={y_corridor.sum()}, settlements={(y_corridor==0).sum()})")

    if n_corridor >= 100:
        corridor_results, _, xgb_corridor, _, _ = run_models(
            X_corridor, y_corridor, feature_names, tag="corridor"
        )

        shap_corridor, gc_rank_corridor, _ = run_shap_analysis(
            xgb_corridor, X_corridor, feature_names,
            OUT / "shap_summary_corridor.png",
            "SHAP Feature Importance — Egypt-Levant-Iran Corridor"
        )
        results['corridor'] = {
            **corridor_results,
            'shap': {'mean_abs_shap': shap_corridor, 'gc_rank': gc_rank_corridor},
            'n_sites': n_corridor
        }
    else:
        print("  Too few corridor sites for separate analysis")
        gc_rank_corridor = None

    # ─── Analysis 3: Ablation Test ───
    print("\n" + "="*60)
    print("ANALYSIS 3: Ablation Test")
    print("="*60)
    ablation = ablation_test(X_arr, y, feature_names)
    results['ablation'] = ablation

    # Corridor ablation too
    if n_corridor >= 100:
        print("\n  Corridor Ablation:")
        ablation_corridor = ablation_test(X_corridor, y_corridor, feature_names)
        results['ablation_corridor'] = ablation_corridor

    # ─── Analysis 4: Partial Dependence Plot ───
    print("\n" + "="*60)
    print("ANALYSIS 4: Partial Dependence Plot")
    print("="*60)
    xgb_model.fit(X_arr, y)  # refit on full data
    pdp_data = partial_dependence_plot(
        xgb_model, X_arr, feature_names,
        output_path=OUT / "partial_dependence.png"
    )
    results['partial_dependence'] = pdp_data

    # ─── Analysis 5: Feature Interaction ───
    print("\n" + "="*60)
    print("ANALYSIS 5: Feature Interaction Analysis")
    print("="*60)
    interactions = feature_interaction_test(
        xgb_model, X_arr, feature_names,
        output_path=OUT / "feature_interactions.png"
    )
    results['interactions'] = interactions

    # ─── Save Results ───
    print("\n" + "="*60)
    print("SAVING RESULTS")
    print("="*60)

    # Model results JSON
    with open(OUT / "model_results.json", 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {OUT / 'model_results.json'}")

    # Logistic coefficients
    lr_coefs = dict(zip(feature_names, lr.coef_[0].tolist()))
    with open(OUT / "logistic_coefficients.json", 'w') as f:
        json.dump(lr_coefs, f, indent=2)

    # Ablation test
    with open(OUT / "ablation_test.json", 'w') as f:
        json.dump({
            'global': ablation,
            'corridor': results.get('ablation_corridor', {})
        }, f, indent=2)

    # Feature importance table
    importance_table = []
    for fname in feature_names:
        row = {
            'feature': fname,
            'rf_importance': results['random_forest']['feature_importances'].get(fname, 0),
            'xgb_importance': results['xgboost']['feature_importances'].get(fname, 0),
            'lr_coefficient': lr_coefs.get(fname, 0),
            'shap_mean_abs': shap_dict.get(fname, 0)
        }
        importance_table.append(row)

    importance_table.sort(key=lambda x: x['shap_mean_abs'], reverse=True)
    with open(OUT / "feature_importance_table.json", 'w') as f:
        json.dump(importance_table, f, indent=2)

    # ─── Generate RESULTS.md ───
    generate_results_md(results, feature_names, importance_table, imputation_log,
                        gc_rank_global, gc_rank_corridor, ablation, len(y),
                        int(y.sum()), int((y==0).sum()), n_corridor)

    print(f"\nAll outputs saved to: {OUT}")
    print("DONE.")


def generate_results_md(results, feature_names, importance_table, imputation_log,
                        gc_rank_global, gc_rank_corridor, ablation, n_total,
                        n_monuments, n_settlements, n_corridor):
    """Generate the RESULTS.md summary."""

    # Determine verdict
    delta = ablation['delta_auc']
    if gc_rank_global <= 3 and delta > 0.02:
        verdict = "STRONG: The Great Circle captures something geography doesn't"
        verdict_detail = "Distance to the Great Circle ranks in the top 3 features and removing it significantly degrades model performance."
    elif gc_rank_global <= 6 and delta > 0.005:
        verdict = "MODEST: The Great Circle adds some information beyond geography"
        verdict_detail = "Distance to the Great Circle has moderate predictive power that is not fully explained by other geographic features."
    elif gc_rank_corridor and gc_rank_corridor <= 3 and gc_rank_global > 6:
        verdict = "REGIONAL: The Great Circle matters in the Egypt-Levant-Iran corridor but not globally"
        verdict_detail = "The Great Circle's predictive power is concentrated in the corridor region, consistent with the Memphis-specific finding."
    else:
        verdict = "WEAK: The Great Circle is largely a geographic proxy"
        verdict_detail = "After controlling for measurable geographic features, the Great Circle adds minimal predictive power."

    md = f"""# ML Feature Importance Test — Results

**Date:** 2026-03-21
**Test:** Does "Distance to Great Circle" predict monument placement after geographic controls?

---

## Verdict

**{verdict}**

{verdict_detail}

---

## Data Summary

| Metric | Value |
|--------|-------|
| Total sites | {n_total} |
| Monuments | {n_monuments} |
| Settlements | {n_settlements} |
| Features | {len(feature_names)} |
| Corridor sites | {n_corridor} |

Features used: {', '.join(feature_names)}

---

## Model Performance (10-fold CV AUC)

| Model | Global AUC | Corridor AUC |
|-------|-----------|-------------|
| Random Forest | {results['random_forest']['cv_auc']:.4f} ± {results['random_forest']['cv_std']:.4f} | {results.get('corridor', {}).get('random_forest', {}).get('cv_auc', 'N/A')} |
| XGBoost | {results['xgboost']['cv_auc']:.4f} ± {results['xgboost']['cv_std']:.4f} | {results.get('corridor', {}).get('xgboost', {}).get('cv_auc', 'N/A')} |
| Logistic Regression | {results['logistic_regression']['cv_auc']:.4f} ± {results['logistic_regression']['cv_std']:.4f} | {results.get('corridor', {}).get('logistic_regression', {}).get('cv_auc', 'N/A')} |

---

## SHAP Feature Importance Ranking

### Global Model

| Rank | Feature | Mean |SHAP| |
|------|---------|-------------|
"""

    ranked = sorted(importance_table, key=lambda x: x['shap_mean_abs'], reverse=True)
    for i, row in enumerate(ranked):
        marker = " **<<<**" if row['feature'] == 'dist_to_gc_km' else ""
        md += f"| {i+1} | {row['feature']} | {row['shap_mean_abs']:.4f}{marker} |\n"

    md += f"""
**GC Distance Rank: #{gc_rank_global} of {len(feature_names)}**

"""

    if gc_rank_corridor:
        md += f"### Corridor Model\n\n**GC Distance Rank: #{gc_rank_corridor} of {len(feature_names)}**\n\n"

    md += f"""---

## Ablation Test

| Metric | Global | Corridor |
|--------|--------|----------|
| AUC with GC | {ablation['auc_with_gc']:.4f} | {results.get('ablation_corridor', {}).get('auc_with_gc', 'N/A')} |
| AUC without GC | {ablation['auc_without_gc']:.4f} | {results.get('ablation_corridor', {}).get('auc_without_gc', 'N/A')} |
| Delta AUC | {ablation['delta_auc']:.4f} | {results.get('ablation_corridor', {}).get('delta_auc', 'N/A')} |

**Interpretation:** {ablation['interpretation']}

---

## Interpretation Table

| Scenario | Threshold | This Result |
|----------|-----------|-------------|
| GC ranks #1-3, ablation > 0.02 | Circle captures something geography doesn't | {'YES' if gc_rank_global <= 3 and delta > 0.02 else 'NO'} |
| GC ranks #4-6, ablation 0.005-0.02 | Circle adds modest information | {'YES' if 4 <= gc_rank_global <= 6 and 0.005 < delta <= 0.02 else 'NO'} |
| GC ranks #7+, ablation < 0.005 | Circle is a geographic proxy | {'YES' if gc_rank_global >= 7 and delta < 0.005 else 'NO'} |
| GC high SHAP in corridor only | Circle matters regionally | {'YES' if gc_rank_corridor and gc_rank_corridor <= 3 and gc_rank_global > 6 else 'NO'} |

---

## Imputation Notes

"""
    for col, action in imputation_log.items():
        md += f"- **{col}**: {action}\n"

    md += f"""
---

## Output Files

- `feature_matrix.csv` — all sites with all features
- `model_results.json` — full model results
- `shap_summary.png` — SHAP beeswarm plot (global)
- `shap_summary_corridor.png` — SHAP beeswarm plot (corridor)
- `ablation_test.json` — AUC with/without GC distance
- `partial_dependence.png` — monument probability vs GC distance
- `feature_interactions.png` — GC distance interactions with other features
- `logistic_coefficients.json` — interpretable model coefficients
- `feature_importance_table.json` — all importance metrics per feature
"""

    with open(OUT / "RESULTS.md", 'w') as f:
        f.write(md)
    print(f"Saved: {OUT / 'RESULTS.md'}")


if __name__ == '__main__':
    main()
