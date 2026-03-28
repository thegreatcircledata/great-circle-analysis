#!/usr/bin/env python3
"""
Paper 3 — Task 2: Portal PPM Replication
==========================================
Rerun PPM pipeline on Megalithic Portal (61,913 sites) instead of Pleiades.
Same covariates, same model comparison, CIs on β_gc.
Flag geographic concentration (65% UK/Ireland).
"""

import csv, json, os, sys, time
import numpy as np
from scipy.optimize import minimize, approx_fprime
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
R_EARTH, QC = 6371.0, 6371.0 * np.pi / 2

MONUMENT_TYPES = {
    'Burial Chamber or Dolmen', 'Standing Stone (Menhir)', 'Round Barrow(s)',
    'Stone Circle', 'Chambered Tomb', 'Long Barrow', 'Passage Grave',
    'Cairn', 'Chambered Cairn', 'Ring Cairn', 'Round Cairn',
    'Stone Row  Alignment', 'Stone Fort or Dun', 'Cursus',
    'Ancient Temple', 'Pyramid  Mastaba', 'Artificial Mound',
    'Rock Cut Tomb', 'Barrow Cemetery', 'Tumulus', 'Henge',
    'Promontory Fort  Cliff Castle', 'Ancient Cross',
    'Standing Stones', 'Broch or Nuraghe', 'Dolmen',
    'Sculptured Stone', 'Misc Earthwork'
}
SETTLEMENT_TYPES = {
    'Ancient Village or Settlement', 'Hillfort', 'Castro or Chafurdão',
    'Cave or Rock Shelter', 'Ancient Mine, Quarry or other Industry',
    'Crannog or Lake Dwelling'
}

# ============================================================
# LOAD INFRASTRUCTURE (same as other PPM scripts)
# ============================================================
def gc_dist_vec(lats, lons):
    la1, lo1 = np.radians(POLE_LAT), np.radians(POLE_LON)
    la2, lo2 = np.radians(lats), np.radians(lons)
    dl, dn = la2 - la1, lo2 - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(la2)*np.sin(dn/2)**2
    d = R_EARTH * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QC)

print("Loading ETOPO1...")
ds = nc.Dataset(os.path.join(BASE, "data/geophysical/etopo/ETOPO1_Ice_g_gmt4.grd"))
etopo_x, etopo_y, etopo_z = ds.variables['x'][:], ds.variables['y'][:], ds.variables['z']
e_lon_min, e_lon_step = float(etopo_x[0]), float(etopo_x[1] - etopo_x[0])
e_lat_min, e_lat_step = float(etopo_y[0]), float(etopo_y[1] - etopo_y[0])

def get_elevation(lats, lons):
    n = len(lats)
    elevs = np.zeros(n)
    for i in range(n):
        li = max(0, min(int(round((lats[i] - e_lat_min) / e_lat_step)), etopo_z.shape[0]-1))
        lo = max(0, min(int(round((lons[i] - e_lon_min) / e_lon_step)), etopo_z.shape[1]-1))
        elevs[i] = float(etopo_z[li, lo])
    return elevs

print("Loading shapefiles...")
coast_points = []
with fiona.open(os.path.join(BASE, "data/natural_earth/ne_50m_coastline.shp")) as src:
    for feat in src:
        geom = sg.shape(feat['geometry'])
        if geom.geom_type == 'LineString':
            coast_points.extend(list(geom.coords)[::10])
        elif geom.geom_type == 'MultiLineString':
            for line in geom.geoms:
                coast_points.extend(list(line.coords)[::10])
coast_arr = np.array(coast_points)

river_points = []
with fiona.open(os.path.join(BASE, "data/natural_earth/ne_50m_rivers_lake_centerlines.shp")) as src:
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

# ============================================================
# LOAD MEGALITHIC PORTAL
# ============================================================
print("\nLoading Megalithic Portal...")
mp_path = os.path.join(BASE, "archive/great-circle-analysis/data/processed/mp_deduplicated.csv")
mp_mon, mp_set = [], []
from collections import Counter
region_counter = Counter()

with open(mp_path) as f:
    for row in csv.DictReader(f):
        try:
            lat, lon = float(row['lat']), float(row['lon'])
        except (ValueError, TypeError):
            continue
        st = row['type']

        # Track geographic concentration
        if 50 < lat < 60 and -10 < lon < 2:
            region_counter['UK/Ireland'] += 1
        elif 45 < lat < 52 and -5 < lon < 10:
            region_counter['France/Belgium'] += 1
        elif 35 < lat < 45 and -10 < lon < 5:
            region_counter['Iberia'] += 1
        elif 35 < lat < 48 and 5 < lon < 20:
            region_counter['Italy/Alps'] += 1
        elif 50 < lat < 70 and 5 < lon < 30:
            region_counter['Scandinavia'] += 1
        else:
            region_counter['Other'] += 1

        if st in MONUMENT_TYPES:
            mp_mon.append((lat, lon))
        elif st in SETTLEMENT_TYPES:
            mp_set.append((lat, lon))

total = sum(region_counter.values())
print(f"Portal sites: {len(mp_mon)} monuments, {len(mp_set)} settlements")
print(f"\nGeographic concentration:")
for region, count in region_counter.most_common():
    print(f"  {region:<20}: {count:>6} ({100*count/total:.1f}%)")

# ============================================================
# COMPUTE COVARIATES
# ============================================================
print("\nComputing covariates...")

def compute_cov(sites):
    lats = np.array([s[0] for s in sites])
    lons = np.array([s[1] for s in sites])
    return {
        'lats': lats, 'lons': lons,
        'elevation': get_elevation(lats, lons),
        'coast_dist': approx_dist_to_features(lats, lons, coast_arr),
        'river_dist': approx_dist_to_features(lats, lons, river_arr),
        'gc_dist': gc_dist_vec(lats, lons),
        'abs_lat': np.abs(lats)
    }

t0 = time.time()
mon_c = compute_cov(mp_mon)
set_c = compute_cov(mp_set)
print(f"  Covariates in {time.time()-t0:.1f}s")

# Standardize using combined stats
all_vals = {}
for key in ['elevation', 'coast_dist', 'river_dist', 'gc_dist', 'abs_lat']:
    combined = np.concatenate([mon_c[key], set_c[key]])
    all_vals[key] = (np.mean(combined), np.std(combined))

def std(arr, key):
    m, s = all_vals[key]
    return (arr - m) / (s + 1e-10)

def build_X(cov, gc=True):
    cols = [std(cov['elevation'], 'elevation'), std(cov['coast_dist'], 'coast_dist'),
            std(cov['river_dist'], 'river_dist'), std(cov['abs_lat'], 'abs_lat')]
    if gc:
        cols.append(std(cov['gc_dist'], 'gc_dist'))
    return np.column_stack(cols)

# ============================================================
# BUILD GRID — Portal is Europe-heavy so use Europe-focused grid
# ============================================================
print("Building quadrature grid (1° resolution, European focus)...")
g_lats = np.arange(30, 72, 1.0)
g_lons = np.arange(-15, 45, 1.0)
gl, gn = np.meshgrid(g_lats, g_lons, indexing='ij')
gl_f, gn_f = gl.ravel(), gn.ravel()
g_elev = get_elevation(gl_f, gn_f)
land = g_elev > 0
g_lat_l, g_lon_l = gl_f[land], gn_f[land]
g_areas = (1 * 111.32) * (1 * 111.32 * np.cos(np.radians(g_lat_l)))
print(f"  Grid points: {len(g_lat_l)}")

# Grid covariates
g_el = get_elevation(g_lat_l, g_lon_l)
g_co = approx_dist_to_features(g_lat_l, g_lon_l, coast_arr)
g_ri = approx_dist_to_features(g_lat_l, g_lon_l, river_arr)
g_gc = gc_dist_vec(g_lat_l, g_lon_l)
g_ab = np.abs(g_lat_l)

def build_grid_X(gc=True):
    cols = [std(g_el, 'elevation'), std(g_co, 'coast_dist'),
            std(g_ri, 'river_dist'), std(g_ab, 'abs_lat')]
    if gc:
        cols.append(std(g_gc, 'gc_dist'))
    return np.column_stack(cols)

# ============================================================
# FIT PPM WITH SEs
# ============================================================
def poisson_negll(beta, X_s, X_g, areas):
    eta_s = beta[0] + X_s @ beta[1:]
    eta_g = beta[0] + X_g @ beta[1:]
    eta_g = np.minimum(eta_g, 20)
    return -(np.sum(eta_s) - np.sum(np.exp(eta_g) * areas))

def fit_with_se(cov, gc, label):
    X_s = build_X(cov, gc)
    X_g = build_grid_X(gc)
    p = X_s.shape[1] + 1
    b0 = np.zeros(p)
    b0[0] = np.log(len(cov['lats']) / np.sum(g_areas))

    res = minimize(poisson_negll, b0, args=(X_s, X_g, g_areas),
                   method='L-BFGS-B', options={'maxiter': 500, 'ftol': 1e-12})
    beta = res.x
    ll = -res.fun
    aic = 2*p - 2*ll

    # Hessian
    eps = 1e-5
    hess = np.zeros((p, p))
    for i in range(p):
        def gi(b):
            return approx_fprime(b, poisson_negll, eps, X_s, X_g, g_areas)[i]
        hess[i, :] = approx_fprime(beta, gi, eps)
    try:
        se = np.sqrt(np.abs(np.diag(np.linalg.inv(hess))))
    except:
        se = np.full(p, np.nan)

    names = ['intercept', 'elevation', 'coast_dist', 'river_dist', 'abs_lat']
    if gc: names.append('gc_dist')

    coefs = {}
    for i, n in enumerate(names):
        zv = beta[i] / (se[i] + 1e-15)
        pv = 2 * norm.sf(abs(zv))
        coefs[n] = {
            'estimate': float(beta[i]), 'se': float(se[i]),
            'z': float(zv), 'p_value': float(pv),
            'ci_lower': float(beta[i] - 1.96*se[i]),
            'ci_upper': float(beta[i] + 1.96*se[i])
        }

    return {'label': label, 'n_sites': len(cov['lats']), 'n_params': p,
            'log_likelihood': float(ll), 'AIC': float(aic),
            'coefficients': coefs, 'converged': res.success}

print("\n" + "=" * 70)
print("FITTING PORTAL PPM MODELS")
print("=" * 70)

t0 = time.time()
pm_gc = fit_with_se(mon_c, True, "Portal Mon +GC")
pm_no = fit_with_se(mon_c, False, "Portal Mon -GC")
ps_gc = fit_with_se(set_c, True, "Portal Set +GC")
ps_no = fit_with_se(set_c, False, "Portal Set -GC")
print(f"Fitted in {time.time()-t0:.1f}s")

# Results
mg = pm_gc['coefficients']['gc_dist']
sg_ = ps_gc['coefficients']['gc_dist']

print(f"\n--- Portal Monument β_gc ---")
print(f"  β = {mg['estimate']:.4f} [{mg['ci_lower']:.4f}, {mg['ci_upper']:.4f}], p = {mg['p_value']:.2e}")

print(f"\n--- Portal Settlement β_gc ---")
print(f"  β = {sg_['estimate']:.4f} [{sg_['ci_lower']:.4f}, {sg_['ci_upper']:.4f}], p = {sg_['p_value']:.2e}")

# LRT
m_lr = 2 * (pm_gc['log_likelihood'] - pm_no['log_likelihood'])
s_lr = 2 * (ps_gc['log_likelihood'] - ps_no['log_likelihood'])
print(f"\n--- LRT ---")
print(f"  Monuments: LR={m_lr:.1f}, p={chi2.sf(m_lr,1):.2e}")
print(f"  Settlements: LR={s_lr:.1f}, p={chi2.sf(s_lr,1):.2e}")

# Divergence
diff = mg['estimate'] - sg_['estimate']
se_d = np.sqrt(mg['se']**2 + sg_['se']**2)
z_d = diff / se_d
p_d = 2 * norm.sf(abs(z_d))
overlap = mg['ci_upper'] >= sg_['ci_lower'] and sg_['ci_upper'] >= mg['ci_lower']
print(f"\n--- Divergence ---")
print(f"  Δβ = {diff:.4f}, z = {z_d:.2f}, p = {p_d:.2e}")
print(f"  CIs overlap: {overlap}")
print(f"  Direction: Mon {'more attracted' if mg['estimate'] < sg_['estimate'] else 'less attracted'}")

# Rhohat
gc_bins = np.arange(0, 2001, 50)
bin_c = (gc_bins[:-1] + gc_bins[1:]) / 2
m_gc_d = mon_c['gc_dist']
s_gc_d = set_c['gc_dist']
g_gc_d = gc_dist_vec(g_lat_l, g_lon_l)
m_cnt, _ = np.histogram(m_gc_d, bins=gc_bins)
s_cnt, _ = np.histogram(s_gc_d, bins=gc_bins)
g_area_bins = np.array([np.sum(g_areas[(g_gc_d >= gc_bins[i]) & (g_gc_d < gc_bins[i+1])]) for i in range(len(gc_bins)-1)])
m_dens = m_cnt / (g_area_bins + 1e-10)
s_dens = s_cnt / (g_area_bins + 1e-10)
m_mean_d = len(mp_mon) / np.sum(g_areas)
s_mean_d = len(mp_set) / np.sum(g_areas)
m_rel = m_dens / (m_mean_d + 1e-15)
s_rel = s_dens / (s_mean_d + 1e-15)

# ============================================================
# SAVE
# ============================================================
results = {
    'dataset': 'Megalithic Portal',
    'n_monuments': len(mp_mon), 'n_settlements': len(mp_set),
    'geographic_concentration': dict(region_counter.most_common()),
    'warning': 'Dataset is ~65% UK/Ireland — European bias may inflate/deflate β_gc',
    'monuments_with_gc': pm_gc,
    'settlements_with_gc': ps_gc,
    'lrt_monuments': {'LR': float(m_lr), 'p': float(chi2.sf(m_lr, 1))},
    'lrt_settlements': {'LR': float(s_lr), 'p': float(chi2.sf(s_lr, 1))},
    'divergence': {
        'delta_beta': float(diff), 'se': float(se_d),
        'z': float(z_d), 'p': float(p_d),
        'cis_overlap': bool(overlap)
    },
    'rhohat': {
        'bin_centers': [float(x) for x in bin_c],
        'mon_relative': [float(x) for x in m_rel],
        'set_relative': [float(x) for x in s_rel]
    }
}

with open(os.path.join(OUT, 'ppm_portal_replication.json'), 'w') as f:
    json.dump(results, f, indent=2)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    close = bin_c <= 500
    ax.plot(bin_c[close], m_rel[close], 'r-', linewidth=2, label='Monuments')
    ax.plot(bin_c[close], s_rel[close], 'b-', linewidth=2, label='Settlements')
    ax.axhline(y=1.0, color='gray', linestyle=':', linewidth=1)
    ax.set_xlabel('Distance from Great Circle (km)')
    ax.set_ylabel('Relative Density')
    ax.set_title('Portal Rhohat: 0-500 km')
    ax.legend()
    ax.grid(alpha=0.3)

    ax2 = axes[1]
    names_plot = ['elevation', 'coast_dist', 'river_dist', 'abs_lat', 'gc_dist']
    mc = [pm_gc['coefficients'][n]['estimate'] for n in names_plot]
    sc = [ps_gc['coefficients'][n]['estimate'] for n in names_plot]
    x = np.arange(5)
    w = 0.35
    ax2.bar(x-w/2, mc, w, label='Monuments', color='#c0392b', alpha=0.8)
    ax2.bar(x+w/2, sc, w, label='Settlements', color='#2980b9', alpha=0.8)
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Elev', 'Coast', 'River', '|Lat|', 'GC Dist'], fontsize=9)
    ax2.set_ylabel('β coefficient')
    ax2.set_title('Portal PPM Coefficients')
    ax2.legend()
    ax2.axhline(y=0, color='black', linewidth=0.5)
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG, '11_portal_ppm.png'), dpi=150)
    print(f"\nPlot: {FIG}/11_portal_ppm.png")
except ImportError:
    print("matplotlib not available")

print(f"Results: {OUT}/ppm_portal_replication.json")
print("Done.")
