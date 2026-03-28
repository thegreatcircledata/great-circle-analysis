#!/usr/bin/env python3
"""
Paper 3 — Inhomogeneous Poisson Point Process Model
=====================================================
Fits monument and settlement intensity as log-linear functions of
geographic covariates, then tests whether adding GC-distance significantly
improves the monument model but NOT the settlement model.

Methodology follows Bilotti et al. (2024) PLOS ONE:
- Inhomogeneous Poisson process via maximum likelihood
- AIC and likelihood ratio test (LRT) for model comparison
- rhohat-style diagnostic: intensity vs GC distance

Implementation in Python (scipy.optimize) rather than R spatstat,
but the statistical framework is identical.

The key insight: in a Poisson point process, the log-likelihood is
  LL = Σ log λ(sᵢ) - ∫ λ(s) ds
where λ(s) = exp(β₀ + β₁x₁(s) + β₂x₂(s) + ...) is the intensity function.

For a global dataset with covariate values at each site, we approximate
the integral by quadrature over a regular grid of "dummy points."
"""

import csv, json, math, os, sys, time
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2
import netCDF4 as nc

np.random.seed(42)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE = os.path.expanduser("~/megalith_site_research")
OUT = os.path.join(BASE, "paper3/results")
FIG = os.path.join(BASE, "paper3/figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

POLE_LAT = 59.682122
POLE_LON = -138.646087
R_EARTH = 6371.0
QC = R_EARTH * np.pi / 2

MONUMENTAL = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus",
    "shrine", "cemetery"
}
SETTLEMENT = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement", "townhouse", "production"
}

# ============================================================
# GEOGRAPHIC FUNCTIONS
# ============================================================
def gc_distance_to_circle(lats, lons):
    """Distance from points to Alison's Great Circle in km."""
    la1 = np.radians(POLE_LAT)
    lo1 = np.radians(POLE_LON)
    la2 = np.radians(lats)
    lo2 = np.radians(lons)
    dl = la2 - la1
    dn = lo2 - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(la2)*np.sin(dn/2)**2
    d = R_EARTH * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QC)

# ============================================================
# LOAD ETOPO1
# ============================================================
print("Loading ETOPO1...")
t0 = time.time()
etopo_path = os.path.join(BASE, "data/geophysical/etopo/ETOPO1_Ice_g_gmt4.grd")
ds = nc.Dataset(etopo_path)
etopo_x = ds.variables['x'][:]
etopo_y = ds.variables['y'][:]
etopo_z = ds.variables['z']
e_lon_min = float(etopo_x[0])
e_lon_step = float(etopo_x[1] - etopo_x[0])
e_lat_min = float(etopo_y[0])
e_lat_step = float(etopo_y[1] - etopo_y[0])
print(f"  ETOPO1 loaded in {time.time()-t0:.1f}s, shape {etopo_z.shape}")

def get_elevation(lats, lons):
    """Extract elevation for arrays of lat/lon from ETOPO1."""
    n = len(lats)
    elevs = np.zeros(n)
    for i in range(n):
        li = int(round((lats[i] - e_lat_min) / e_lat_step))
        lo = int(round((lons[i] - e_lon_min) / e_lon_step))
        li = max(0, min(li, etopo_z.shape[0]-1))
        lo = max(0, min(lo, etopo_z.shape[1]-1))
        elevs[i] = float(etopo_z[li, lo])
    return elevs

# ============================================================
# LOAD COASTLINE AND RIVER SHAPEFILES
# ============================================================
print("Loading coastline shapefile...")
import shapely.geometry as sg
import fiona

coast_path = os.path.join(BASE, "data/natural_earth/ne_50m_coastline.shp")
coast_points = []
with fiona.open(coast_path) as src:
    for feat in src:
        geom = sg.shape(feat['geometry'])
        if geom.geom_type == 'LineString':
            coords = list(geom.coords)
            coast_points.extend(coords[::10])  # subsample for speed
        elif geom.geom_type == 'MultiLineString':
            for line in geom.geoms:
                coords = list(line.coords)
                coast_points.extend(coords[::10])
coast_arr = np.array(coast_points)
print(f"  Coast points: {len(coast_arr)}")

river_path = os.path.join(BASE, "data/natural_earth/ne_50m_rivers_lake_centerlines.shp")
river_points = []
with fiona.open(river_path) as src:
    for feat in src:
        if feat['geometry'] is None:
            continue
        geom = sg.shape(feat['geometry'])
        if geom.geom_type == 'LineString':
            coords = list(geom.coords)
            river_points.extend(coords[::5])
        elif geom.geom_type == 'MultiLineString':
            for line in geom.geoms:
                coords = list(line.coords)
                river_points.extend(coords[::5])
river_arr = np.array(river_points)
print(f"  River points: {len(river_arr)}")

def approx_dist_to_features(lats, lons, feature_arr):
    """Approximate min distance to nearest feature point (degrees → km)."""
    n = len(lats)
    dists = np.zeros(n)
    for i in range(n):
        dlat = feature_arr[:, 1] - lats[i]
        dlon = feature_arr[:, 0] - lons[i]
        # Approximate km using cos(lat) correction
        cos_lat = np.cos(np.radians(lats[i]))
        d_km = np.sqrt((dlat * 111.32)**2 + (dlon * 111.32 * cos_lat)**2)
        dists[i] = np.min(d_km)
    return dists

# ============================================================
# LOAD PLEIADES
# ============================================================
print("\nLoading Pleiades data...")
monuments = []
settlements = []

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
        is_m = bool(fts & MONUMENTAL)
        is_s = bool(fts & SETTLEMENT)
        if is_m and not is_s:
            monuments.append((lat, lon))
        elif is_s and not is_m:
            settlements.append((lat, lon))

print(f"Monuments: {len(monuments)}, Settlements: {len(settlements)}")

# ============================================================
# COMPUTE COVARIATES FOR ALL SITES
# ============================================================
print("\nComputing covariates...")

def compute_covariates(sites, label):
    """Compute all covariates for a list of (lat, lon) tuples."""
    lats = np.array([s[0] for s in sites])
    lons = np.array([s[1] for s in sites])
    n = len(sites)

    print(f"  {label}: elevation...")
    elev = get_elevation(lats, lons)

    print(f"  {label}: coast distance...")
    coast_d = approx_dist_to_features(lats, lons, coast_arr)

    print(f"  {label}: river distance...")
    river_d = approx_dist_to_features(lats, lons, river_arr)

    print(f"  {label}: GC distance...")
    gc_d = gc_distance_to_circle(lats, lons)

    abs_lat = np.abs(lats)

    return {
        'lats': lats, 'lons': lons,
        'elevation': elev,
        'coast_dist': coast_d,
        'river_dist': river_d,
        'gc_dist': gc_d,
        'abs_lat': abs_lat
    }

t0 = time.time()
mon_cov = compute_covariates(monuments, "Monuments")
set_cov = compute_covariates(settlements, "Settlements")
print(f"Covariates computed in {time.time()-t0:.1f}s")

# ============================================================
# GENERATE QUADRATURE GRID (dummy points)
# ============================================================
print("\nGenerating quadrature grid...")
# 2° resolution grid over land areas where Pleiades has sites
# Pleiades is mostly Mediterranean/Middle East: lat 20-60, lon -15 to 90
grid_lats = np.arange(20, 60, 2.0)
grid_lons = np.arange(-15, 90, 2.0)
grid_lat_mesh, grid_lon_mesh = np.meshgrid(grid_lats, grid_lons, indexing='ij')
grid_lat_flat = grid_lat_mesh.ravel()
grid_lon_flat = grid_lon_mesh.ravel()

# Filter to land only (elevation > 0)
grid_elev = get_elevation(grid_lat_flat, grid_lon_flat)
land_mask = grid_elev > 0
grid_lat_land = grid_lat_flat[land_mask]
grid_lon_land = grid_lon_flat[land_mask]
n_quad = len(grid_lat_land)
print(f"Quadrature points: {n_quad} (land cells at 2° resolution)")

# Compute covariates for grid
print("Computing grid covariates...")
grid_elev_land = grid_elev[land_mask]
grid_coast_d = approx_dist_to_features(grid_lat_land, grid_lon_land, coast_arr)
grid_river_d = approx_dist_to_features(grid_lat_land, grid_lon_land, river_arr)
grid_gc_d = gc_distance_to_circle(grid_lat_land, grid_lon_land)
grid_abs_lat = np.abs(grid_lat_land)

# Cell area (km²) for each grid cell
cell_area = (2 * 111.32) * (2 * 111.32 * np.cos(np.radians(grid_lat_land)))

# ============================================================
# STANDARDIZE COVARIATES
# ============================================================
def standardize(arr, mean=None, std=None):
    if mean is None:
        mean = np.mean(arr)
    if std is None:
        std = np.std(arr)
    return (arr - mean) / (std + 1e-10), mean, std

# Compute stats from combined monument+settlement data
all_elev = np.concatenate([mon_cov['elevation'], set_cov['elevation']])
all_coast = np.concatenate([mon_cov['coast_dist'], set_cov['coast_dist']])
all_river = np.concatenate([mon_cov['river_dist'], set_cov['river_dist']])
all_gc = np.concatenate([mon_cov['gc_dist'], set_cov['gc_dist']])
all_ablat = np.concatenate([mon_cov['abs_lat'], set_cov['abs_lat']])

stats = {}
_, stats['elev_m'], stats['elev_s'] = standardize(all_elev)
_, stats['coast_m'], stats['coast_s'] = standardize(all_coast)
_, stats['river_m'], stats['river_s'] = standardize(all_river)
_, stats['gc_m'], stats['gc_s'] = standardize(all_gc)
_, stats['ablat_m'], stats['ablat_s'] = standardize(all_ablat)

def build_design_matrix(cov_dict, include_gc=True):
    """Build standardized design matrix from covariate dict."""
    X = np.column_stack([
        standardize(cov_dict['elevation'], stats['elev_m'], stats['elev_s'])[0],
        standardize(cov_dict['coast_dist'], stats['coast_m'], stats['coast_s'])[0],
        standardize(cov_dict['river_dist'], stats['river_m'], stats['river_s'])[0],
        standardize(cov_dict['abs_lat'], stats['ablat_m'], stats['ablat_s'])[0],
    ])
    if include_gc:
        X = np.column_stack([
            X,
            standardize(cov_dict['gc_dist'], stats['gc_m'], stats['gc_s'])[0]
        ])
    return X

def build_grid_matrix(include_gc=True):
    """Build standardized design matrix for quadrature grid."""
    X = np.column_stack([
        standardize(grid_elev_land, stats['elev_m'], stats['elev_s'])[0],
        standardize(grid_coast_d, stats['coast_m'], stats['coast_s'])[0],
        standardize(grid_river_d, stats['river_m'], stats['river_s'])[0],
        standardize(grid_abs_lat, stats['ablat_m'], stats['ablat_s'])[0],
    ])
    if include_gc:
        X = np.column_stack([
            X,
            standardize(grid_gc_d, stats['gc_m'], stats['gc_s'])[0]
        ])
    return X

# ============================================================
# POISSON POINT PROCESS LOG-LIKELIHOOD
# ============================================================
def poisson_ppm_negll(beta, X_sites, X_grid, grid_areas):
    """
    Negative log-likelihood for inhomogeneous Poisson process.

    LL = Σᵢ (β₀ + X_sites[i] @ β[1:]) - Σⱼ exp(β₀ + X_grid[j] @ β[1:]) * area[j]

    Parameters:
    - beta: coefficients [intercept, β₁, β₂, ...]
    - X_sites: design matrix at data points (n_sites × p)
    - X_grid: design matrix at quadrature points (n_quad × p)
    - grid_areas: area weights for quadrature (km²)
    """
    n_sites = X_sites.shape[0]

    # Linear predictor at data points
    eta_sites = beta[0] + X_sites @ beta[1:]

    # Linear predictor at quadrature points
    eta_grid = beta[0] + X_grid @ beta[1:]

    # Clamp to prevent overflow
    eta_grid = np.minimum(eta_grid, 20)

    # Log-likelihood
    ll = np.sum(eta_sites) - np.sum(np.exp(eta_grid) * grid_areas)

    return -ll  # return negative for minimization

def fit_ppm(sites_cov, include_gc, label):
    """Fit a Poisson PPM and return results."""
    X_sites = build_design_matrix(sites_cov, include_gc)
    X_grid = build_grid_matrix(include_gc)
    p = X_sites.shape[1] + 1  # +1 for intercept

    beta0 = np.zeros(p)
    # Initialize intercept to log(n_sites / total_area)
    total_area = np.sum(cell_area)
    n_sites = X_sites.shape[0]
    beta0[0] = np.log(n_sites / total_area)

    result = minimize(
        poisson_ppm_negll,
        beta0,
        args=(X_sites, X_grid, cell_area),
        method='L-BFGS-B',
        options={'maxiter': 500, 'ftol': 1e-10}
    )

    if not result.success:
        print(f"  WARNING: optimization did not converge for {label}: {result.message}")

    beta = result.x
    negll = result.fun
    ll = -negll
    aic = 2 * p - 2 * ll

    cov_names = ['elevation', 'coast_dist', 'river_dist', 'abs_lat']
    if include_gc:
        cov_names.append('gc_dist')

    coefs = {'intercept': float(beta[0])}
    for i, name in enumerate(cov_names):
        coefs[name] = float(beta[i+1])

    return {
        'label': label,
        'include_gc': include_gc,
        'n_sites': n_sites,
        'n_params': p,
        'log_likelihood': float(ll),
        'AIC': float(aic),
        'coefficients': coefs,
        'converged': result.success
    }

# ============================================================
# FIT MODELS
# ============================================================
print("\n" + "=" * 70)
print("FITTING POINT PROCESS MODELS")
print("=" * 70)

print("\n--- Monument Models ---")
t0 = time.time()
mon_no_gc = fit_ppm(mon_cov, include_gc=False, label="Monuments without GC")
print(f"  Model A (no GC): LL={mon_no_gc['log_likelihood']:.1f}, AIC={mon_no_gc['AIC']:.1f}, "
      f"converged={mon_no_gc['converged']}")
print(f"  Coefficients: {mon_no_gc['coefficients']}")

mon_with_gc = fit_ppm(mon_cov, include_gc=True, label="Monuments with GC")
print(f"  Model B (+GC):   LL={mon_with_gc['log_likelihood']:.1f}, AIC={mon_with_gc['AIC']:.1f}, "
      f"converged={mon_with_gc['converged']}")
print(f"  Coefficients: {mon_with_gc['coefficients']}")
print(f"  Time: {time.time()-t0:.1f}s")

print("\n--- Settlement Models ---")
t0 = time.time()
set_no_gc = fit_ppm(set_cov, include_gc=False, label="Settlements without GC")
print(f"  Model A (no GC): LL={set_no_gc['log_likelihood']:.1f}, AIC={set_no_gc['AIC']:.1f}, "
      f"converged={set_no_gc['converged']}")
print(f"  Coefficients: {set_no_gc['coefficients']}")

set_with_gc = fit_ppm(set_cov, include_gc=True, label="Settlements with GC")
print(f"  Model B (+GC):   LL={set_with_gc['log_likelihood']:.1f}, AIC={set_with_gc['AIC']:.1f}, "
      f"converged={set_with_gc['converged']}")
print(f"  Coefficients: {set_with_gc['coefficients']}")
print(f"  Time: {time.time()-t0:.1f}s")

# ============================================================
# LIKELIHOOD RATIO TESTS
# ============================================================
print("\n" + "=" * 70)
print("LIKELIHOOD RATIO TESTS")
print("=" * 70)

def lrt(model_null, model_alt, label):
    """Likelihood ratio test: does adding GC distance improve fit?"""
    ll_null = model_null['log_likelihood']
    ll_alt = model_alt['log_likelihood']
    df = model_alt['n_params'] - model_null['n_params']
    lr_stat = 2 * (ll_alt - ll_null)
    p_value = chi2.sf(lr_stat, df)
    delta_aic = model_null['AIC'] - model_alt['AIC']  # positive = alt is better

    return {
        'label': label,
        'LR_statistic': float(lr_stat),
        'df': df,
        'p_value': float(p_value),
        'delta_AIC': float(delta_aic),
        'significant_005': bool(p_value < 0.05),
        'significant_001': bool(p_value < 0.01)
    }

mon_lrt = lrt(mon_no_gc, mon_with_gc, "Monuments: does GC distance improve fit?")
set_lrt = lrt(set_no_gc, set_with_gc, "Settlements: does GC distance improve fit?")

print(f"\nMonuments:")
print(f"  LR statistic: {mon_lrt['LR_statistic']:.2f}")
print(f"  p-value: {mon_lrt['p_value']:.6f}")
print(f"  ΔAIC: {mon_lrt['delta_AIC']:.1f} ({'GC model better' if mon_lrt['delta_AIC'] > 0 else 'Base model better'})")
print(f"  Significant at α=0.05: {mon_lrt['significant_005']}")
print(f"  GC distance coefficient: {mon_with_gc['coefficients']['gc_dist']:.4f}")

print(f"\nSettlements:")
print(f"  LR statistic: {set_lrt['LR_statistic']:.2f}")
print(f"  p-value: {set_lrt['p_value']:.6f}")
print(f"  ΔAIC: {set_lrt['delta_AIC']:.1f} ({'GC model better' if set_lrt['delta_AIC'] > 0 else 'Base model better'})")
print(f"  Significant at α=0.05: {set_lrt['significant_005']}")
print(f"  GC distance coefficient: {set_with_gc['coefficients']['gc_dist']:.4f}")

# ============================================================
# KEY DIAGNOSTIC: DIVERGENCE IN GC COEFFICIENT
# ============================================================
print("\n" + "=" * 70)
print("DIVERGENCE DIAGNOSTIC")
print("=" * 70)

gc_coef_mon = mon_with_gc['coefficients']['gc_dist']
gc_coef_set = set_with_gc['coefficients']['gc_dist']

print(f"\nGC distance coefficient:")
print(f"  Monuments:   β_gc = {gc_coef_mon:.4f}  ({'NEGATIVE = attracted to circle' if gc_coef_mon < 0 else 'POSITIVE = repelled from circle'})")
print(f"  Settlements: β_gc = {gc_coef_set:.4f}  ({'NEGATIVE = attracted to circle' if gc_coef_set < 0 else 'POSITIVE = repelled from circle'})")
print(f"  Difference:  Δβ  = {gc_coef_mon - gc_coef_set:.4f}")

if gc_coef_mon < 0 and (gc_coef_set >= 0 or gc_coef_set > gc_coef_mon):
    print("\n→ DIVERGENCE CONFIRMED: Monuments are attracted to the circle")
    print("  (negative coefficient = higher intensity at smaller GC distance)")
    if gc_coef_set >= 0:
        print("  while settlements are repelled or indifferent.")
    else:
        print("  more strongly than settlements.")
else:
    print("\n→ No clear divergence in GC coefficient direction.")

# ============================================================
# RHOHAT-STYLE DIAGNOSTIC
# ============================================================
print("\n" + "=" * 70)
print("RHOHAT DIAGNOSTIC: Intensity vs GC Distance")
print("=" * 70)

# Bin sites by GC distance and compute density per distance band
gc_bins = np.arange(0, 2001, 50)  # 0 to 2000 km in 50 km bins
bin_centers = (gc_bins[:-1] + gc_bins[1:]) / 2

mon_gc_dists = mon_cov['gc_dist']
set_gc_dists = set_cov['gc_dist']

# Count sites per bin
mon_counts, _ = np.histogram(mon_gc_dists, bins=gc_bins)
set_counts, _ = np.histogram(set_gc_dists, bins=gc_bins)

# Compute expected count per bin from quadrature grid
grid_gc_dists = grid_gc_d
grid_in_bins = np.zeros(len(gc_bins) - 1)
grid_area_in_bins = np.zeros(len(gc_bins) - 1)
for i in range(len(gc_bins) - 1):
    mask = (grid_gc_dists >= gc_bins[i]) & (grid_gc_dists < gc_bins[i+1])
    grid_area_in_bins[i] = np.sum(cell_area[mask])

# Density = count / area (sites per km²)
mon_density = mon_counts / (grid_area_in_bins + 1e-10)
set_density = set_counts / (grid_area_in_bins + 1e-10)

# Relative density (normalized to mean)
mon_mean_density = len(monuments) / np.sum(cell_area)
set_mean_density = len(settlements) / np.sum(cell_area)
mon_relative = mon_density / (mon_mean_density + 1e-15)
set_relative = set_density / (set_mean_density + 1e-15)

print("\nRelative density by GC distance band:")
print(f"{'Distance (km)':<15} {'Mon rel.':>10} {'Set rel.':>10} {'Mon/Set':>10}")
print("-" * 50)
for i in range(min(20, len(bin_centers))):
    ratio = mon_relative[i] / (set_relative[i] + 1e-10)
    marker = " ***" if mon_relative[i] > 2.0 and ratio > 1.5 else ""
    print(f"{bin_centers[i]:>10.0f}     {mon_relative[i]:>10.2f} {set_relative[i]:>10.2f} {ratio:>10.2f}{marker}")

# ============================================================
# SAVE
# ============================================================
results = {
    'models': {
        'monuments_no_gc': mon_no_gc,
        'monuments_with_gc': mon_with_gc,
        'settlements_no_gc': set_no_gc,
        'settlements_with_gc': set_with_gc
    },
    'likelihood_ratio_tests': {
        'monuments': mon_lrt,
        'settlements': set_lrt
    },
    'divergence': {
        'gc_coef_monuments': float(gc_coef_mon),
        'gc_coef_settlements': float(gc_coef_set),
        'difference': float(gc_coef_mon - gc_coef_set),
        'monument_attracted': bool(gc_coef_mon < 0),
        'settlement_attracted': bool(gc_coef_set < 0)
    },
    'rhohat': {
        'bin_centers_km': [float(x) for x in bin_centers],
        'mon_relative_density': [float(x) for x in mon_relative],
        'set_relative_density': [float(x) for x in set_relative]
    },
    'standardization': {k: float(v) for k, v in stats.items()},
    'quadrature': {'n_points': n_quad, 'resolution_deg': 2.0}
}

with open(os.path.join(OUT, '08_point_process_model.json'), 'w') as f:
    json.dump(results, f, indent=2)

# ============================================================
# PLOT
# ============================================================
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Left: Rhohat — relative density vs GC distance
    ax = axes[0]
    ax.plot(bin_centers, mon_relative, 'r-', linewidth=2, label='Monuments')
    ax.plot(bin_centers, set_relative, 'b-', linewidth=2, label='Settlements')
    ax.axhline(y=1.0, color='gray', linestyle=':', linewidth=1)
    ax.fill_between(bin_centers, 0, mon_relative, alpha=0.1, color='red')
    ax.set_xlabel('Distance from Great Circle (km)')
    ax.set_ylabel('Relative Density (obs / expected)')
    ax.set_title('Intensity vs GC Distance (rhohat diagnostic)')
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 1500)

    # Middle: Zoomed rhohat (0-200 km)
    ax2 = axes[1]
    close_mask = bin_centers <= 200
    ax2.plot(bin_centers[close_mask], mon_relative[close_mask], 'ro-', linewidth=2,
             markersize=4, label='Monuments')
    ax2.plot(bin_centers[close_mask], set_relative[close_mask], 'bs-', linewidth=2,
             markersize=4, label='Settlements')
    ax2.axhline(y=1.0, color='gray', linestyle=':', linewidth=1)
    ax2.fill_between(bin_centers[close_mask], 1.0, mon_relative[close_mask],
                     where=mon_relative[close_mask] > 1.0, alpha=0.2, color='red')
    ax2.set_xlabel('Distance from Great Circle (km)')
    ax2.set_ylabel('Relative Density')
    ax2.set_title('Zoomed: 0-200 km from GC')
    ax2.legend()
    ax2.grid(alpha=0.3)

    # Right: Coefficient comparison
    ax3 = axes[2]
    cov_names_base = ['elevation', 'coast_dist', 'river_dist', 'abs_lat']
    mon_coefs = [mon_with_gc['coefficients'][c] for c in cov_names_base + ['gc_dist']]
    set_coefs = [set_with_gc['coefficients'][c] for c in cov_names_base + ['gc_dist']]
    x = np.arange(5)
    w = 0.35
    ax3.bar(x - w/2, mon_coefs, w, label='Monuments', color='#c0392b', alpha=0.8)
    ax3.bar(x + w/2, set_coefs, w, label='Settlements', color='#2980b9', alpha=0.8)
    ax3.set_xticks(x)
    ax3.set_xticklabels(['Elevation', 'Coast\nDist', 'River\nDist', '|Latitude|', 'GC\nDist'],
                        fontsize=9)
    ax3.set_ylabel('Coefficient (β)')
    ax3.set_title('PPM Coefficients: Monuments vs Settlements')
    ax3.legend()
    ax3.axhline(y=0, color='black', linewidth=0.5)
    ax3.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG, '08_point_process_model.png'), dpi=150)
    print(f"\nPlot saved: {os.path.join(FIG, '08_point_process_model.png')}")
except ImportError:
    print("matplotlib not available")

print(f"Results saved: {os.path.join(OUT, '08_point_process_model.json')}")
print("Done.")
