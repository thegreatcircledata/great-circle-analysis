#!/usr/bin/env python3
"""
Paper 3 Robustness Battery — Task 3: PPM Residual Spatial Autocorrelation
=========================================================================
After fitting the inhomogeneous Poisson PPM, compute residual Moran's I
for both monument and settlement models. If significant, fit a Thomas
cluster process as a robustness check.
"""

import csv, json, math, os, sys, time
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2, norm
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
R = 6371.0
QC = R * np.pi / 2

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

# ============================================================
# LOAD DATA (same as PPM code)
# ============================================================
print("Loading ETOPO1...")
etopo_path = os.path.join(BASE, "data/geophysical/etopo/ETOPO1_Ice_g_gmt4.grd")
ds = nc.Dataset(etopo_path)
etopo_x = ds.variables['x'][:]
etopo_y = ds.variables['y'][:]
etopo_z = ds.variables['z']
e_lon_min = float(etopo_x[0])
e_lon_step = float(etopo_x[1] - etopo_x[0])
e_lat_min = float(etopo_y[0])
e_lat_step = float(etopo_y[1] - etopo_y[0])

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

print("Loading coastline and rivers...")
import shapely.geometry as sg
import fiona

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
        if feat['geometry'] is None: continue
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

print("Loading Pleiades...")
monuments, settlements = [], []
with open(os.path.join(BASE, "data/pleiades/pleiades-places-latest.csv")) as f:
    for row in csv.DictReader(f):
        try:
            lat, lon = float(row['reprLat']), float(row['reprLong'])
        except (ValueError, TypeError):
            continue
        if abs(lat) < 0.001 and abs(lon) < 0.001: continue
        fts = set(ft.strip() for ft in row['featureTypes'].split(','))
        if (fts & MONUMENTAL) and not (fts & SETTLEMENT):
            monuments.append((lat, lon))
        elif (fts & SETTLEMENT) and not (fts & MONUMENTAL):
            settlements.append((lat, lon))

print(f"Monuments: {len(monuments)}, Settlements: {len(settlements)}")

# ============================================================
# FIT PPM (reuse methodology)
# ============================================================
print("\nSetting up PPM...")

# Grid at 2° resolution
grid_lats = np.arange(20, 60, 2.0)
grid_lons = np.arange(-15, 90, 2.0)
glm, glo = np.meshgrid(grid_lats, grid_lons, indexing='ij')
gf_lat, gf_lon = glm.ravel(), glo.ravel()
g_elev = get_elevation(gf_lat, gf_lon)
land = g_elev > 0
gl_lat, gl_lon = gf_lat[land], gf_lon[land]
gl_elev = g_elev[land]
gl_coast = approx_dist_to_features(gl_lat, gl_lon, coast_arr)
gl_river = approx_dist_to_features(gl_lat, gl_lon, river_arr)
gl_gc = gc_distance_to_circle(gl_lat, gl_lon)
gl_ablat = np.abs(gl_lat)
cell_area = (2 * 111.32) * (2 * 111.32 * np.cos(np.radians(gl_lat)))

def compute_site_covs(sites):
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

print("Computing monument covariates...")
mon_cov = compute_site_covs(monuments)
print("Computing settlement covariates...")
set_cov = compute_site_covs(settlements)

# Standardization stats from combined data
all_data = {}
for key in ['elevation', 'coast_dist', 'river_dist', 'gc_dist', 'abs_lat']:
    combined = np.concatenate([mon_cov[key], set_cov[key]])
    all_data[key] = (np.mean(combined), np.std(combined))

def std(arr, key):
    m, s = all_data[key]
    return (arr - m) / (s + 1e-10)

def build_X(cov):
    return np.column_stack([
        std(cov['elevation'], 'elevation'),
        std(cov['coast_dist'], 'coast_dist'),
        std(cov['river_dist'], 'river_dist'),
        std(cov['abs_lat'], 'abs_lat'),
        std(cov['gc_dist'], 'gc_dist')
    ])

def build_grid_X():
    return np.column_stack([
        std(gl_elev, 'elevation'),
        std(gl_coast, 'coast_dist'),
        std(gl_river, 'river_dist'),
        std(gl_ablat, 'abs_lat'),
        std(gl_gc, 'gc_dist')
    ])

def ppm_negll(beta, X_s, X_g, areas):
    eta_s = beta[0] + X_s @ beta[1:]
    eta_g = np.minimum(beta[0] + X_g @ beta[1:], 20)
    return -(np.sum(eta_s) - np.sum(np.exp(eta_g) * areas))

def fit(cov, label):
    X_s = build_X(cov)
    X_g = build_grid_X()
    p = X_s.shape[1] + 1
    b0 = np.zeros(p)
    b0[0] = np.log(len(cov['lats']) / np.sum(cell_area))
    res = minimize(ppm_negll, b0, args=(X_s, X_g, cell_area),
                   method='L-BFGS-B', options={'maxiter': 500, 'ftol': 1e-10})
    print(f"  {label}: converged={res.success}, LL={-res.fun:.1f}")
    return res.x, X_s, X_g

print("\nFitting PPMs...")
beta_mon, X_mon, X_grid = fit(mon_cov, "Monuments")
beta_set, X_set, _ = fit(set_cov, "Settlements")

# ============================================================
# COMPUTE PEARSON RESIDUALS ON GRID
# ============================================================
print("\n" + "=" * 70)
print("COMPUTING PEARSON RESIDUALS")
print("=" * 70)

def compute_grid_residuals(beta, sites_cov, X_grid, label):
    """
    Compute Pearson residuals on the quadrature grid.
    For each grid cell: observed count - expected count, normalized.
    """
    lats = sites_cov['lats']
    lons = sites_cov['lons']
    n_grid = len(gl_lat)

    # Expected intensity at each grid cell
    eta = beta[0] + X_grid @ beta[1:]
    eta = np.minimum(eta, 20)
    lambda_g = np.exp(eta)
    expected = lambda_g * cell_area  # expected count per cell

    # Count observed sites in each grid cell
    observed = np.zeros(n_grid)
    for k in range(len(lats)):
        # Find nearest grid cell
        dlat = np.abs(gl_lat - lats[k])
        dlon = np.abs(gl_lon - lons[k])
        idx = np.argmin(dlat + dlon)
        observed[idx] += 1

    # Pearson residuals: (O - E) / sqrt(E)
    pearson = (observed - expected) / np.sqrt(expected + 1e-10)

    print(f"\n{label}:")
    print(f"  Grid cells: {n_grid}")
    print(f"  Total observed: {int(np.sum(observed))}")
    print(f"  Total expected: {np.sum(expected):.1f}")
    print(f"  Pearson residual range: [{np.min(pearson):.2f}, {np.max(pearson):.2f}]")
    print(f"  Mean |residual|: {np.mean(np.abs(pearson)):.3f}")

    return observed, expected, pearson

obs_mon, exp_mon, res_mon = compute_grid_residuals(beta_mon, mon_cov, X_grid, "Monuments")
obs_set, exp_set, res_set = compute_grid_residuals(beta_set, set_cov, X_grid, "Settlements")

# ============================================================
# MORAN'S I ON RESIDUALS
# ============================================================
print("\n" + "=" * 70)
print("MORAN'S I — RESIDUAL SPATIAL AUTOCORRELATION")
print("=" * 70)

def morans_I(residuals, lats, lons, max_dist_km=500):
    """
    Compute Moran's I for spatial autocorrelation of residuals.
    Uses inverse-distance weighting with a cutoff.
    """
    n = len(residuals)
    z = residuals - np.mean(residuals)

    # Build weight matrix (inverse distance, truncated)
    # For efficiency, use vectorized distance computation
    print(f"  Computing Moran's I (n={n}, max_dist={max_dist_km}km)...")

    W_sum = 0.0
    numerator = 0.0
    denominator = np.sum(z**2)

    # For large n, subsample pairs
    if n > 2000:
        print(f"  Subsampling pairs (n={n} too large for full matrix)...")
        n_pairs = 500000
        rng = np.random.RandomState(42)
        idx_i = rng.randint(0, n, n_pairs)
        idx_j = rng.randint(0, n, n_pairs)
        # Remove self-pairs
        mask = idx_i != idx_j
        idx_i = idx_i[mask]
        idx_j = idx_j[mask]

        dlat = (lats[idx_i] - lats[idx_j]) * 111.32
        cos_lat = np.cos(np.radians((lats[idx_i] + lats[idx_j]) / 2))
        dlon = (lons[idx_i] - lons[idx_j]) * 111.32 * cos_lat
        d = np.sqrt(dlat**2 + dlon**2)

        within = d <= max_dist_km
        w = np.zeros(len(d))
        w[within] = 1.0 / (d[within] + 1e-3)  # inverse distance

        numerator = np.sum(w * z[idx_i] * z[idx_j])
        W_sum = np.sum(w)
    else:
        for i in range(n):
            dlat = (lats - lats[i]) * 111.32
            cos_lat = np.cos(np.radians(lats[i]))
            dlon = (lons - lons[i]) * 111.32 * cos_lat
            d = np.sqrt(dlat**2 + dlon**2)
            d[i] = np.inf  # exclude self
            within = d <= max_dist_km
            w = np.zeros(n)
            w[within] = 1.0 / (d[within] + 1e-3)
            numerator += np.sum(w * z[i] * z)
            W_sum += np.sum(w)

    if W_sum == 0 or denominator == 0:
        return 0.0, 0.0, 1.0

    I = (n / W_sum) * (numerator / denominator)

    # Expected I under null
    E_I = -1.0 / (n - 1)

    # Variance under randomization (simplified)
    V_I = 1.0 / (n - 1)  # approximate
    Z_I = (I - E_I) / np.sqrt(V_I)
    p_value = 2 * (1 - norm.cdf(abs(Z_I)))

    return I, Z_I, p_value

I_mon, Z_mon, p_mon = morans_I(res_mon, gl_lat, gl_lon)
print(f"\nMonuments: Moran's I = {I_mon:.4f}, Z = {Z_mon:.2f}, p = {p_mon:.6f}")
print(f"  {'SIGNIFICANT' if p_mon < 0.05 else 'Not significant'} residual autocorrelation")

I_set, Z_set, p_set = morans_I(res_set, gl_lat, gl_lon)
print(f"\nSettlements: Moran's I = {I_set:.4f}, Z = {Z_set:.2f}, p = {p_set:.6f}")
print(f"  {'SIGNIFICANT' if p_set < 0.05 else 'Not significant'} residual autocorrelation")

# ============================================================
# THOMAS CLUSTER PROCESS (if autocorrelation significant)
# ============================================================
print("\n" + "=" * 70)
print("THOMAS CLUSTER PROCESS ROBUSTNESS CHECK")
print("=" * 70)

def thomas_cluster_enrichment(sites_cov, label, n_mc=1000):
    """
    Test whether the GC enrichment signal persists when we account
    for spatial clustering by using a Thomas process null model.

    Instead of comparing to CSR, we fit a Thomas process (Poisson
    cluster process) to the observed pattern and test whether
    GC distance is still a significant predictor after accounting
    for the clustering.

    Approach: compute the proportion of sites within 50km of GC,
    then simulate from a Thomas process with matched intensity
    and cluster parameters, and compare.
    """
    lats = sites_cov['lats']
    lons = sites_cov['lons']
    gc_dists = sites_cov['gc_dist']
    n = len(lats)

    obs_frac = np.sum(gc_dists <= 50) / n

    # Estimate Thomas process parameters from nearest-neighbor distances
    nn_dists = np.zeros(n)
    for i in range(n):
        dlat = (lats - lats[i]) * 111.32
        cos_lat = np.cos(np.radians(lats[i]))
        dlon = (lons - lons[i]) * 111.32 * cos_lat
        d = np.sqrt(dlat**2 + dlon**2)
        d[i] = np.inf
        nn_dists[i] = np.min(d)

    mean_nn = np.mean(nn_dists)
    sigma_c = mean_nn / 2  # cluster scale ≈ half mean NN distance
    # Estimate mean cluster size from variance/mean ratio of cell counts
    cell_counts = obs_mon if 'Mon' in label else obs_set
    mu_c = max(1, np.var(cell_counts) / (np.mean(cell_counts) + 1e-10))

    print(f"\n{label}:")
    print(f"  Mean NN distance: {mean_nn:.1f} km")
    print(f"  Estimated cluster scale (σ): {sigma_c:.1f} km")
    print(f"  Estimated mean cluster size (μ): {mu_c:.1f}")
    print(f"  Observed fraction within 50km of GC: {obs_frac:.4f}")

    # Simulate Thomas process n_mc times
    # For each simulation: generate parent points from fitted intensity,
    # then scatter offspring around each parent
    rng = np.random.RandomState(42)
    sim_fracs = []

    # Use the grid-based intensity for parent placement
    lambda_parent = np.exp(beta_mon[0] + X_grid @ beta_mon[1:]) if 'Mon' in label else \
                    np.exp(beta_set[0] + X_grid @ beta_set[1:])
    lambda_parent = np.minimum(lambda_parent, np.exp(20))
    parent_expected = lambda_parent * cell_area / mu_c  # parents per cell

    for trial in range(n_mc):
        # Generate parents (Poisson per cell)
        n_parents_per_cell = rng.poisson(parent_expected)
        total_parents = np.sum(n_parents_per_cell)

        if total_parents == 0:
            sim_fracs.append(0.0)
            continue

        # Place parents at grid cell centers
        parent_lats = np.repeat(gl_lat, n_parents_per_cell)
        parent_lons = np.repeat(gl_lon, n_parents_per_cell)

        # Generate offspring around each parent
        n_offspring_per_parent = rng.poisson(mu_c, len(parent_lats))
        total_offspring = np.sum(n_offspring_per_parent)

        if total_offspring == 0:
            sim_fracs.append(0.0)
            continue

        off_parent_idx = np.repeat(np.arange(len(parent_lats)), n_offspring_per_parent)
        # Gaussian displacement in km
        dx = rng.normal(0, sigma_c, total_offspring)
        dy = rng.normal(0, sigma_c, total_offspring)

        off_lats = parent_lats[off_parent_idx] + dy / 111.32
        off_lons = parent_lons[off_parent_idx] + dx / (111.32 * np.cos(np.radians(parent_lats[off_parent_idx])) + 1e-10)

        # GC distance for offspring
        off_gc = gc_distance_to_circle(off_lats, off_lons)
        sim_frac = np.sum(off_gc <= 50) / len(off_gc) if len(off_gc) > 0 else 0
        sim_fracs.append(sim_frac)

    sim_fracs = np.array(sim_fracs)
    p_thomas = np.mean(sim_fracs >= obs_frac)

    print(f"  Thomas process simulations: {n_mc}")
    print(f"  Simulated mean fraction in corridor: {np.mean(sim_fracs):.4f}")
    print(f"  Observed fraction: {obs_frac:.4f}")
    print(f"  p-value (Thomas null): {p_thomas:.4f}")
    print(f"  Enrichment {'PERSISTS' if p_thomas < 0.05 else 'does not persist'} under Thomas null")

    return {
        'label': label,
        'mean_nn_dist_km': float(mean_nn),
        'cluster_scale_km': float(sigma_c),
        'cluster_size': float(mu_c),
        'observed_frac_in_corridor': float(obs_frac),
        'thomas_mean_frac': float(np.mean(sim_fracs)),
        'thomas_p_value': float(p_thomas),
        'enrichment_persists': bool(p_thomas < 0.05),
        'n_simulations': n_mc
    }

thomas_mon = thomas_cluster_enrichment(mon_cov, "Monuments (Thomas)", n_mc=500)
thomas_set = thomas_cluster_enrichment(set_cov, "Settlements (Thomas)", n_mc=500)

# ============================================================
# SAVE
# ============================================================
results = {
    'analysis': 'PPM Residual Spatial Autocorrelation',
    'morans_i': {
        'monuments': {
            'I': float(I_mon), 'Z': float(Z_mon), 'p_value': float(p_mon),
            'significant': bool(p_mon < 0.05)
        },
        'settlements': {
            'I': float(I_set), 'Z': float(Z_set), 'p_value': float(p_set),
            'significant': bool(p_set < 0.05)
        }
    },
    'thomas_process': {
        'monuments': thomas_mon,
        'settlements': thomas_set
    },
    'conclusion': (
        'Monument enrichment persists under Thomas cluster null'
        if thomas_mon.get('enrichment_persists', False)
        else 'Monument enrichment does NOT persist under Thomas cluster null'
    )
}

with open(os.path.join(OUT, 'ppm_residual_autocorrelation.json'), 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved: {os.path.join(OUT, 'ppm_residual_autocorrelation.json')}")
print("Done.")
