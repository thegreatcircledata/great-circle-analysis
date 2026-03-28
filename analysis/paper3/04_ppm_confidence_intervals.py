#!/usr/bin/env python3
"""
Paper 3 PPM — Task 1: β Confidence Intervals + Task 3: Quadrature Sensitivity
===============================================================================
Extract standard errors via Hessian inversion, compute 95% CIs,
test whether monument and settlement β_gc CIs overlap.
Also rerun at 0.5° quadrature resolution to test sensitivity.
"""

import csv, json, math, os, sys, time
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2, norm
import netCDF4 as nc
import fiona
import shapely.geometry as sg

np.random.seed(42)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE = os.path.expanduser("~/megalith_site_research")
OUT = os.path.join(BASE, "paper3/results")
FIG = os.path.join(BASE, "paper3/figures")

POLE_LAT, POLE_LON = 59.682122, -138.646087
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
# LOAD DATA (same as 03_point_process_model.py)
# ============================================================
def gc_distance_to_circle(lats, lons):
    la1, lo1 = np.radians(POLE_LAT), np.radians(POLE_LON)
    la2, lo2 = np.radians(lats), np.radians(lons)
    dl, dn = la2 - la1, lo2 - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(la2)*np.sin(dn/2)**2
    d = R_EARTH * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QC)

print("Loading data...")
etopo_path = os.path.join(BASE, "data/geophysical/etopo/ETOPO1_Ice_g_gmt4.grd")
ds = nc.Dataset(etopo_path)
etopo_x, etopo_y, etopo_z = ds.variables['x'][:], ds.variables['y'][:], ds.variables['z']
e_lon_min, e_lon_step = float(etopo_x[0]), float(etopo_x[1] - etopo_x[0])
e_lat_min, e_lat_step = float(etopo_y[0]), float(etopo_y[1] - etopo_y[0])

def get_elevation(lats, lons):
    n = len(lats)
    elevs = np.zeros(n)
    for i in range(n):
        li = int(round((lats[i] - e_lat_min) / e_lat_step))
        lo = int(round((lons[i] - e_lon_min) / e_lon_step))
        li = max(0, min(li, etopo_z.shape[0]-1))
        lo = max(0, min(lo, etopo_z.shape[1]-1))
        elevs[i] = float(etopo_z[li, lo])
    return elevs

coast_path = os.path.join(BASE, "data/natural_earth/ne_50m_coastline.shp")
coast_points = []
with fiona.open(coast_path) as src:
    for feat in src:
        geom = sg.shape(feat['geometry'])
        if geom.geom_type == 'LineString':
            coast_points.extend(list(geom.coords)[::10])
        elif geom.geom_type == 'MultiLineString':
            for line in geom.geoms:
                coast_points.extend(list(line.coords)[::10])
coast_arr = np.array(coast_points)

river_path = os.path.join(BASE, "data/natural_earth/ne_50m_rivers_lake_centerlines.shp")
river_points = []
with fiona.open(river_path) as src:
    for feat in src:
        if feat['geometry'] is None:
            continue
        geom = sg.shape(feat['geometry'])
        if geom.geom_type == 'LineString':
            river_points.extend(list(geom.coords)[::5])
        elif geom.geom_type == 'MultiLineString':
            for line in geom.geoms:
                river_points.extend(list(line.coords)[::5])
river_arr = np.array(river_points)

def approx_dist_to_features(lats, lons, feature_arr):
    n = len(lats)
    dists = np.zeros(n)
    for i in range(n):
        dlat = feature_arr[:, 1] - lats[i]
        dlon = feature_arr[:, 0] - lons[i]
        cos_lat = np.cos(np.radians(lats[i]))
        d_km = np.sqrt((dlat * 111.32)**2 + (dlon * 111.32 * cos_lat)**2)
        dists[i] = np.min(d_km)
    return dists

# Load Pleiades
monuments, settlements = [], []
with open(os.path.join(BASE, "data/pleiades/pleiades-places-latest.csv")) as f:
    for row in csv.DictReader(f):
        try:
            lat, lon = float(row['reprLat']), float(row['reprLong'])
        except (ValueError, TypeError):
            continue
        if abs(lat) < 0.001 and abs(lon) < 0.001:
            continue
        fts = set(ft.strip() for ft in row['featureTypes'].split(','))
        is_m, is_s = bool(fts & MONUMENTAL), bool(fts & SETTLEMENT)
        if is_m and not is_s:
            monuments.append((lat, lon))
        elif is_s and not is_m:
            settlements.append((lat, lon))

print(f"Monuments: {len(monuments)}, Settlements: {len(settlements)}")

def compute_covariates(sites):
    lats = np.array([s[0] for s in sites])
    lons = np.array([s[1] for s in sites])
    return {
        'lats': lats, 'lons': lons,
        'elevation': get_elevation(lats, lons),
        'coast_dist': approx_dist_to_features(lats, lons, coast_arr),
        'river_dist': approx_dist_to_features(lats, lons, river_arr),
        'gc_dist': gc_distance_to_circle(lats, lons),
        'abs_lat': np.abs(lats)
    }

print("Computing site covariates...")
mon_cov = compute_covariates(monuments)
set_cov = compute_covariates(settlements)

# Standardization stats
all_vals = {}
for key in ['elevation', 'coast_dist', 'river_dist', 'gc_dist', 'abs_lat']:
    combined = np.concatenate([mon_cov[key], set_cov[key]])
    all_vals[key] = (np.mean(combined), np.std(combined))

def standardize(arr, key):
    m, s = all_vals[key]
    return (arr - m) / (s + 1e-10)

def build_X(cov, include_gc=True):
    cols = [standardize(cov['elevation'], 'elevation'),
            standardize(cov['coast_dist'], 'coast_dist'),
            standardize(cov['river_dist'], 'river_dist'),
            standardize(cov['abs_lat'], 'abs_lat')]
    if include_gc:
        cols.append(standardize(cov['gc_dist'], 'gc_dist'))
    return np.column_stack(cols)

# ============================================================
# PPM FITTING WITH HESSIAN FOR STANDARD ERRORS
# ============================================================
def poisson_negll(beta, X_sites, X_grid, grid_areas):
    eta_sites = beta[0] + X_sites @ beta[1:]
    eta_grid = beta[0] + X_grid @ beta[1:]
    eta_grid = np.minimum(eta_grid, 20)
    ll = np.sum(eta_sites) - np.sum(np.exp(eta_grid) * grid_areas)
    return -ll

def fit_ppm_with_se(sites_cov, include_gc, grid_lats, grid_lons, grid_areas, label):
    X_sites = build_X(sites_cov, include_gc)

    # Grid covariates
    g_elev = get_elevation(grid_lats, grid_lons)
    g_coast = approx_dist_to_features(grid_lats, grid_lons, coast_arr)
    g_river = approx_dist_to_features(grid_lats, grid_lons, river_arr)
    g_gc = gc_distance_to_circle(grid_lats, grid_lons)
    g_ablat = np.abs(grid_lats)

    cols = [standardize(g_elev, 'elevation'),
            standardize(g_coast, 'coast_dist'),
            standardize(g_river, 'river_dist'),
            standardize(g_ablat, 'abs_lat')]
    if include_gc:
        cols.append(standardize(g_gc, 'gc_dist'))
    X_grid = np.column_stack(cols)

    p = X_sites.shape[1] + 1
    n_sites = X_sites.shape[0]
    beta0 = np.zeros(p)
    beta0[0] = np.log(n_sites / np.sum(grid_areas))

    result = minimize(poisson_negll, beta0, args=(X_sites, X_grid, grid_areas),
                      method='L-BFGS-B', options={'maxiter': 500, 'ftol': 1e-12})

    beta = result.x
    ll = -result.fun
    aic = 2 * p - 2 * ll

    # Hessian via finite differences for standard errors
    from scipy.optimize import approx_fprime
    eps = 1e-5
    hess = np.zeros((p, p))
    for i in range(p):
        def grad_i(b):
            g = approx_fprime(b, poisson_negll, eps, X_sites, X_grid, grid_areas)
            return g[i]
        hess[i, :] = approx_fprime(beta, grad_i, eps)

    try:
        cov_matrix = np.linalg.inv(hess)  # inverse of Hessian = variance-covariance
        se = np.sqrt(np.abs(np.diag(cov_matrix)))
    except np.linalg.LinAlgError:
        se = np.full(p, np.nan)

    cov_names = ['intercept', 'elevation', 'coast_dist', 'river_dist', 'abs_lat']
    if include_gc:
        cov_names.append('gc_dist')

    coefs = {}
    for i, name in enumerate(cov_names):
        z_val = beta[i] / (se[i] + 1e-15)
        p_val = 2 * norm.sf(abs(z_val))
        coefs[name] = {
            'estimate': float(beta[i]),
            'se': float(se[i]),
            'z': float(z_val),
            'p_value': float(p_val),
            'ci_lower': float(beta[i] - 1.96 * se[i]),
            'ci_upper': float(beta[i] + 1.96 * se[i])
        }

    return {
        'label': label, 'include_gc': include_gc,
        'n_sites': n_sites, 'n_params': p,
        'log_likelihood': float(ll), 'AIC': float(aic),
        'coefficients': coefs, 'converged': result.success
    }

# ============================================================
# TASK 1: CIs at 2° resolution (original)
# ============================================================
print("\n" + "=" * 70)
print("TASK 1: β CONFIDENCE INTERVALS (2° grid)")
print("=" * 70)

def make_grid(step):
    g_lats = np.arange(20, 60, step)
    g_lons = np.arange(-15, 90, step)
    gl, gn = np.meshgrid(g_lats, g_lons, indexing='ij')
    gl_f, gn_f = gl.ravel(), gn.ravel()
    elev = get_elevation(gl_f, gn_f)
    land = elev > 0
    lat_land, lon_land = gl_f[land], gn_f[land]
    areas = (step * 111.32) * (step * 111.32 * np.cos(np.radians(lat_land)))
    return lat_land, lon_land, areas

print("Building 2° grid...")
g2_lat, g2_lon, g2_area = make_grid(2.0)
print(f"  Grid points: {len(g2_lat)}")

print("\nFitting models with SEs (2° grid)...")
t0 = time.time()
mon_2deg = fit_ppm_with_se(mon_cov, True, g2_lat, g2_lon, g2_area, "Mon 2°")
set_2deg = fit_ppm_with_se(set_cov, True, g2_lat, g2_lon, g2_area, "Set 2°")
mon_2deg_noGC = fit_ppm_with_se(mon_cov, False, g2_lat, g2_lon, g2_area, "Mon 2° noGC")
set_2deg_noGC = fit_ppm_with_se(set_cov, False, g2_lat, g2_lon, g2_area, "Set 2° noGC")
print(f"  Fitted in {time.time()-t0:.1f}s")

# Print table
print("\n--- Monument Model (with GC) ---")
print(f"{'Covariate':<15} {'β':>8} {'SE':>8} {'z':>8} {'p':>10} {'95% CI':>22}")
print("-" * 75)
for name, c in mon_2deg['coefficients'].items():
    sig = "***" if c['p_value'] < 0.001 else "**" if c['p_value'] < 0.01 else "*" if c['p_value'] < 0.05 else ""
    print(f"{name:<15} {c['estimate']:>8.4f} {c['se']:>8.4f} {c['z']:>8.2f} {c['p_value']:>10.2e} [{c['ci_lower']:>8.4f}, {c['ci_upper']:>8.4f}] {sig}")

print("\n--- Settlement Model (with GC) ---")
print(f"{'Covariate':<15} {'β':>8} {'SE':>8} {'z':>8} {'p':>10} {'95% CI':>22}")
print("-" * 75)
for name, c in set_2deg['coefficients'].items():
    sig = "***" if c['p_value'] < 0.001 else "**" if c['p_value'] < 0.01 else "*" if c['p_value'] < 0.05 else ""
    print(f"{name:<15} {c['estimate']:>8.4f} {c['se']:>8.4f} {c['z']:>8.2f} {c['p_value']:>10.2e} [{c['ci_lower']:>8.4f}, {c['ci_upper']:>8.4f}] {sig}")

# CI overlap test for gc_dist
mon_gc = mon_2deg['coefficients']['gc_dist']
set_gc = set_2deg['coefficients']['gc_dist']
overlap = mon_gc['ci_upper'] >= set_gc['ci_lower'] and set_gc['ci_upper'] >= mon_gc['ci_lower']
print(f"\n--- GC Distance Coefficient Comparison ---")
print(f"  Monuments:   β_gc = {mon_gc['estimate']:.4f} [{mon_gc['ci_lower']:.4f}, {mon_gc['ci_upper']:.4f}]")
print(f"  Settlements: β_gc = {set_gc['estimate']:.4f} [{set_gc['ci_lower']:.4f}, {set_gc['ci_upper']:.4f}]")
print(f"  CIs overlap: {overlap}")

# Formal z-test for difference
diff = mon_gc['estimate'] - set_gc['estimate']
se_diff = np.sqrt(mon_gc['se']**2 + set_gc['se']**2)
z_diff = diff / se_diff
p_diff = 2 * norm.sf(abs(z_diff))
print(f"  Difference:  Δβ = {diff:.4f} ± {se_diff:.4f}")
print(f"  z-test:      z = {z_diff:.2f}, p = {p_diff:.2e}")
print(f"  Formally significant (p < 0.05): {p_diff < 0.05}")

# LRT
mon_lrt_stat = 2 * (mon_2deg['log_likelihood'] - mon_2deg_noGC['log_likelihood'])
set_lrt_stat = 2 * (set_2deg['log_likelihood'] - set_2deg_noGC['log_likelihood'])
mon_lrt_p = chi2.sf(mon_lrt_stat, 1)
set_lrt_p = chi2.sf(set_lrt_stat, 1)

print(f"\n--- LRT: Does GC distance improve fit? ---")
print(f"  Monuments:   LR = {mon_lrt_stat:.1f}, p = {mon_lrt_p:.2e}, ΔAIC = {mon_2deg_noGC['AIC'] - mon_2deg['AIC']:.1f}")
print(f"  Settlements: LR = {set_lrt_stat:.1f}, p = {set_lrt_p:.2e}, ΔAIC = {set_2deg_noGC['AIC'] - set_2deg['AIC']:.1f}")

# ============================================================
# TASK 3: 0.5° QUADRATURE GRID
# ============================================================
print("\n" + "=" * 70)
print("TASK 3: FINER QUADRATURE GRID (0.5°)")
print("=" * 70)

print("Building 0.5° grid...")
g05_lat, g05_lon, g05_area = make_grid(0.5)
print(f"  Grid points: {len(g05_lat)}")

print("Fitting models with SEs (0.5° grid)...")
t0 = time.time()
mon_05deg = fit_ppm_with_se(mon_cov, True, g05_lat, g05_lon, g05_area, "Mon 0.5°")
set_05deg = fit_ppm_with_se(set_cov, True, g05_lat, g05_lon, g05_area, "Set 0.5°")
print(f"  Fitted in {time.time()-t0:.1f}s")

mon_gc_05 = mon_05deg['coefficients']['gc_dist']
set_gc_05 = set_05deg['coefficients']['gc_dist']

print(f"\n--- 0.5° Grid Results ---")
print(f"  Mon β_gc = {mon_gc_05['estimate']:.4f} [{mon_gc_05['ci_lower']:.4f}, {mon_gc_05['ci_upper']:.4f}]")
print(f"  Set β_gc = {set_gc_05['estimate']:.4f} [{set_gc_05['ci_lower']:.4f}, {set_gc_05['ci_upper']:.4f}]")

print(f"\n--- Sensitivity Comparison ---")
print(f"  {'Metric':<25} {'2° grid':>15} {'0.5° grid':>15}")
print(f"  {'-'*55}")
print(f"  {'Mon β_gc':<25} {mon_gc['estimate']:>15.4f} {mon_gc_05['estimate']:>15.4f}")
print(f"  {'Set β_gc':<25} {set_gc['estimate']:>15.4f} {set_gc_05['estimate']:>15.4f}")
print(f"  {'Mon SE':<25} {mon_gc['se']:>15.4f} {mon_gc_05['se']:>15.4f}")
print(f"  {'Set SE':<25} {set_gc['se']:>15.4f} {set_gc_05['se']:>15.4f}")
print(f"  {'Δβ (mon-set)':<25} {diff:>15.4f} {mon_gc_05['estimate']-set_gc_05['estimate']:>15.4f}")

# Rhohat at 0.5° resolution
gc_bins = np.arange(0, 2001, 50)
bin_centers = (gc_bins[:-1] + gc_bins[1:]) / 2

mon_gc_dists = mon_cov['gc_dist']
set_gc_dists = set_cov['gc_dist']
grid_gc_dists = gc_distance_to_circle(g05_lat, g05_lon)

mon_counts, _ = np.histogram(mon_gc_dists, bins=gc_bins)
set_counts, _ = np.histogram(set_gc_dists, bins=gc_bins)

grid_area_in_bins = np.zeros(len(gc_bins) - 1)
for i in range(len(gc_bins) - 1):
    mask = (grid_gc_dists >= gc_bins[i]) & (grid_gc_dists < gc_bins[i+1])
    grid_area_in_bins[i] = np.sum(g05_area[mask])

mon_density = mon_counts / (grid_area_in_bins + 1e-10)
set_density = set_counts / (grid_area_in_bins + 1e-10)
mon_mean_d = len(monuments) / np.sum(g05_area)
set_mean_d = len(settlements) / np.sum(g05_area)
mon_rel = mon_density / (mon_mean_d + 1e-15)
set_rel = set_density / (set_mean_d + 1e-15)

print(f"\n--- Rhohat at 0.5° (first 10 bins) ---")
print(f"{'Dist km':<10} {'Mon rel':>10} {'Set rel':>10} {'Mon/Set':>10}")
for i in range(min(10, len(bin_centers))):
    ratio = mon_rel[i] / (set_rel[i] + 1e-10)
    print(f"{bin_centers[i]:>8.0f}   {mon_rel[i]:>10.2f} {set_rel[i]:>10.2f} {ratio:>10.2f}")

# ============================================================
# SAVE
# ============================================================
ci_results = {
    'task1_confidence_intervals': {
        'grid_resolution': '2°',
        'monuments': mon_2deg,
        'settlements': set_2deg,
        'gc_dist_comparison': {
            'mon_beta': mon_gc['estimate'], 'mon_se': mon_gc['se'],
            'mon_ci': [mon_gc['ci_lower'], mon_gc['ci_upper']],
            'set_beta': set_gc['estimate'], 'set_se': set_gc['se'],
            'set_ci': [set_gc['ci_lower'], set_gc['ci_upper']],
            'cis_overlap': bool(overlap),
            'difference': float(diff), 'se_diff': float(se_diff),
            'z_test': float(z_diff), 'p_value_diff': float(p_diff),
            'formally_significant': bool(p_diff < 0.05)
        },
        'lrt_monuments': {'LR': float(mon_lrt_stat), 'p': float(mon_lrt_p)},
        'lrt_settlements': {'LR': float(set_lrt_stat), 'p': float(set_lrt_p)}
    },
    'task3_quadrature_sensitivity': {
        'fine_grid_resolution': '0.5°',
        'fine_grid_points': len(g05_lat),
        'coarse_grid_points': len(g2_lat),
        'monuments_05': {'beta_gc': mon_gc_05['estimate'], 'se': mon_gc_05['se'],
                         'ci': [mon_gc_05['ci_lower'], mon_gc_05['ci_upper']]},
        'settlements_05': {'beta_gc': set_gc_05['estimate'], 'se': set_gc_05['se'],
                           'ci': [set_gc_05['ci_lower'], set_gc_05['ci_upper']]},
        'rhohat_05': {
            'bin_centers_km': [float(x) for x in bin_centers],
            'mon_relative': [float(x) for x in mon_rel],
            'set_relative': [float(x) for x in set_rel]
        }
    }
}

with open(os.path.join(OUT, 'ppm_coefficients_with_ci.json'), 'w') as f:
    json.dump(ci_results, f, indent=2)

# CSV table
with open(os.path.join(OUT, 'ppm_coefficients_with_ci.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['model', 'covariate', 'beta', 'se', 'z', 'p_value', 'ci_lower', 'ci_upper'])
    for model_name, model in [('monuments_2deg', mon_2deg), ('settlements_2deg', set_2deg),
                                ('monuments_05deg', mon_05deg), ('settlements_05deg', set_05deg)]:
        for cov_name, c in model['coefficients'].items():
            w.writerow([model_name, cov_name, f"{c['estimate']:.6f}", f"{c['se']:.6f}",
                        f"{c['z']:.4f}", f"{c['p_value']:.2e}",
                        f"{c['ci_lower']:.6f}", f"{c['ci_upper']:.6f}"])

with open(os.path.join(OUT, 'ppm_quadrature_sensitivity.json'), 'w') as f:
    json.dump(ci_results['task3_quadrature_sensitivity'], f, indent=2)

# ============================================================
# PLOT
# ============================================================
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Left: CI comparison
    ax = axes[0]
    labels = ['Monuments\n(2°)', 'Settlements\n(2°)', 'Monuments\n(0.5°)', 'Settlements\n(0.5°)']
    betas = [mon_gc['estimate'], set_gc['estimate'], mon_gc_05['estimate'], set_gc_05['estimate']]
    ci_lo = [mon_gc['ci_lower'], set_gc['ci_lower'], mon_gc_05['ci_lower'], set_gc_05['ci_lower']]
    ci_hi = [mon_gc['ci_upper'], set_gc['ci_upper'], mon_gc_05['ci_upper'], set_gc_05['ci_upper']]
    colors = ['#c0392b', '#2980b9', '#e74c3c', '#3498db']

    x_pos = [0, 1, 2.5, 3.5]
    for i in range(4):
        ax.errorbar(x_pos[i], betas[i], yerr=[[betas[i]-ci_lo[i]], [ci_hi[i]-betas[i]]],
                    fmt='o', color=colors[i], markersize=8, capsize=6, capthick=2, linewidth=2)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=9)
    ax.axhline(y=0, color='gray', linestyle=':', linewidth=1)
    ax.set_ylabel('β_gc (GC distance coefficient)')
    ax.set_title('GC Distance Coefficient with 95% CIs')
    ax.grid(axis='y', alpha=0.3)

    # Middle: Rhohat at 0.5° resolution
    ax2 = axes[1]
    ax2.plot(bin_centers, mon_rel, 'r-', linewidth=2, label='Monuments')
    ax2.plot(bin_centers, set_rel, 'b-', linewidth=2, label='Settlements')
    ax2.axhline(y=1.0, color='gray', linestyle=':', linewidth=1)
    ax2.fill_between(bin_centers, 1.0, mon_rel, where=mon_rel > 1, alpha=0.15, color='red')
    ax2.set_xlabel('Distance from Great Circle (km)')
    ax2.set_ylabel('Relative Density')
    ax2.set_title('Rhohat at 0.5° Grid Resolution')
    ax2.legend()
    ax2.grid(alpha=0.3)
    ax2.set_xlim(0, 1500)

    # Right: Zoomed rhohat 0-300 km
    ax3 = axes[2]
    close = bin_centers <= 300
    ax3.plot(bin_centers[close], mon_rel[close], 'ro-', linewidth=2, markersize=5, label='Monuments')
    ax3.plot(bin_centers[close], set_rel[close], 'bs-', linewidth=2, markersize=5, label='Settlements')
    ax3.axhline(y=1.0, color='gray', linestyle=':', linewidth=1)
    ax3.fill_between(bin_centers[close], 1.0, mon_rel[close],
                     where=mon_rel[close] > 1, alpha=0.2, color='red')
    ax3.set_xlabel('Distance from Great Circle (km)')
    ax3.set_ylabel('Relative Density')
    ax3.set_title('Zoomed Rhohat: 0-300 km')
    ax3.legend()
    ax3.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG, '09_ppm_ci_and_sensitivity.png'), dpi=150)
    print(f"\nPlot saved: {FIG}/09_ppm_ci_and_sensitivity.png")
except ImportError:
    print("matplotlib not available")

print(f"\nResults saved to {OUT}/")
print("Done.")
