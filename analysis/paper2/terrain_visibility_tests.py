#!/usr/bin/env python3
"""
Terrain Boundary, Visibility, and Navigational Waypoint Tests
=============================================================
Tests A-E from the terrain-visibility directive.

A: Terrain Transition Zone — does the GC trace terrain gradient boundaries?
B: Viewshed/Visibility — do GC locations have higher topographic prominence?
C: Navigational Waypoint — are GC coast crossings at prominent headlands?
D: Nile Valley Constriction — does the GC cross at the valley's narrowest?
E: Monument Elevation Differential — monuments on high ground, settlements low?
"""

import json, math, os, sys, time
import numpy as np
import netCDF4 as nc
import pandas as pd
from scipy.ndimage import sobel
from scipy.stats import percentileofscore, mannwhitneyu
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

# ============================================================
# CONFIGURATION
# ============================================================
POLE_LAT = 59.682122
POLE_LON = -138.646087
R = 6371.0

BASE = Path(os.path.expanduser("~/megalith_site_research"))
OUT = BASE / "outputs" / "terrain_visibility_tests"
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# LOAD ETOPO1
# ============================================================
print("=" * 70)
print("LOADING ETOPO1")
print("=" * 70)

etopo_path = BASE / "data" / "geophysical" / "etopo" / "ETOPO1_Ice_g_gmt4.grd"
ds = nc.Dataset(str(etopo_path))
etopo_lons = ds.variables['x'][:]
etopo_lats = ds.variables['y'][:]
etopo_z = ds.variables['z']

lon_min, lon_step = float(etopo_lons[0]), float(etopo_lons[1] - etopo_lons[0])
lat_min, lat_step = float(etopo_lats[0]), float(etopo_lats[1] - etopo_lats[0])
n_lats = len(etopo_lats)
n_lons = len(etopo_lons)
print(f"ETOPO1: {n_lats} x {n_lons}, step={lat_step:.4f}°")


def get_elev(lat, lon):
    """Single point elevation lookup."""
    li = int(round((lat - lat_min) / lat_step))
    lo = int(round((lon - lon_min) / lon_step))
    li = max(0, min(li, n_lats - 1))
    lo = max(0, min(lo, n_lons - 1))
    return float(etopo_z[li, lo])


def get_elev_batch(lats, lons):
    """Vectorized elevation lookup."""
    lat_idx = np.round((np.asarray(lats) - lat_min) / lat_step).astype(int)
    lon_idx = np.round((np.asarray(lons) - lon_min) / lon_step).astype(int)
    lat_idx = np.clip(lat_idx, 0, n_lats - 1)
    lon_idx = np.clip(lon_idx, 0, n_lons - 1)
    elevs = np.empty(len(lat_idx))
    for i in range(len(lat_idx)):
        elevs[i] = float(etopo_z[lat_idx[i], lon_idx[i]])
    return elevs


def sample_elevation_grid(lat, lon, radius_km=50, n=20):
    """Sample an n×n elevation grid centered on (lat, lon) within radius_km."""
    dlat = radius_km / 111.0
    dlon = radius_km / (111.0 * max(np.cos(np.radians(lat)), 0.01))
    lat_arr = np.linspace(lat - dlat, lat + dlat, n)
    lon_arr = np.linspace(lon - dlon, lon + dlon, n)
    grid = np.zeros((n, n))
    for i, la in enumerate(lat_arr):
        for j, lo in enumerate(lon_arr):
            grid[i, j] = get_elev(la, lo)
    return grid


def haversine_km(lat1, lon1, lat2, lon2):
    """Great circle distance in km."""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    return R * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def gc_distance_scalar(lat, lon):
    """Distance from point to GC in km."""
    pole_lat_r = np.radians(POLE_LAT)
    pole_lon_r = np.radians(POLE_LON)
    lat_r = np.radians(lat)
    lon_r = np.radians(lon)
    cos_ang = (np.sin(pole_lat_r)*np.sin(lat_r) +
               np.cos(pole_lat_r)*np.cos(lat_r)*np.cos(lon_r - pole_lon_r))
    ang = np.arccos(np.clip(cos_ang, -1, 1))
    return abs(ang - np.pi/2) * R


def trace_great_circle(pole_lat, pole_lon, step_deg=1.0):
    """Trace GC at given step. Returns (lats, lons)."""
    pole = np.array([
        np.cos(np.radians(pole_lat)) * np.cos(np.radians(pole_lon)),
        np.cos(np.radians(pole_lat)) * np.sin(np.radians(pole_lon)),
        np.sin(np.radians(pole_lat))
    ])
    if abs(pole[2]) < 0.9:
        u = np.cross(pole, [0, 0, 1])
    else:
        u = np.cross(pole, [1, 0, 0])
    u = u / np.linalg.norm(u)
    v = np.cross(pole, u)
    v = v / np.linalg.norm(v)

    n_pts = int(360.0 / step_deg)
    thetas = np.linspace(0, 2*np.pi, n_pts, endpoint=False)
    pts = np.outer(np.cos(thetas), u) + np.outer(np.sin(thetas), v)
    lats = np.degrees(np.arcsin(np.clip(pts[:, 2], -1, 1)))
    lons = np.degrees(np.arctan2(pts[:, 1], pts[:, 0]))
    return lats, lons


def random_pole():
    """Uniform random pole on sphere."""
    z = np.random.uniform(-1, 1)
    theta = np.random.uniform(0, 2*np.pi)
    lat = np.degrees(np.arcsin(z))
    lon = np.degrees(theta) - 180
    return lat, lon


def is_land(lat, lon):
    """Quick check if point is land (elevation > 0)."""
    return get_elev(lat, lon) > 0


# ============================================================
# TRACE THE GREAT CIRCLE
# ============================================================
print("\nTracing Great Circle at 1° intervals...")
gc_lats, gc_lons = trace_great_circle(POLE_LAT, POLE_LON, step_deg=1.0)
print(f"  {len(gc_lats)} points")

# Identify land points
gc_elevs = get_elev_batch(gc_lats, gc_lons)
gc_land_mask = gc_elevs > 0
gc_land_lats = gc_lats[gc_land_mask]
gc_land_lons = gc_lons[gc_land_mask]
print(f"  {gc_land_mask.sum()} land points, {(~gc_land_mask).sum()} ocean points")


# ================================================================
# TEST A: TERRAIN TRANSITION ZONE
# ================================================================
print(f"\n{'='*70}")
print("TEST A: TERRAIN TRANSITION ZONE")
print(f"{'='*70}")

def compute_gradient(lat, lon, radius_km=50, n=20):
    """Compute mean terrain gradient magnitude at a point."""
    grid = sample_elevation_grid(lat, lon, radius_km, n)
    # Set ocean cells to 0 for gradient purposes
    grid = np.maximum(grid, 0)
    sx = sobel(grid, axis=0)
    sy = sobel(grid, axis=1)
    mag = np.sqrt(sx**2 + sy**2)
    return float(mag.mean())

# Compute gradient along the GC (land points only to make fair comparison)
print("Computing terrain gradient along GC land points...")
t0 = time.time()
gc_gradients = []
for i, (la, lo) in enumerate(zip(gc_land_lats, gc_land_lons)):
    gc_gradients.append(compute_gradient(la, lo))
    if (i+1) % 20 == 0:
        print(f"  {i+1}/{len(gc_land_lats)} ({time.time()-t0:.0f}s)")
gc_mean_gradient = np.mean(gc_gradients)
print(f"GC mean gradient: {gc_mean_gradient:.2f}")

# Random circle comparison (reduce to 200 for speed, still statistically adequate)
N_RANDOM = 200
print(f"\nComputing gradients for {N_RANDOM} random circles...")
random_gradients = []
for trial in range(N_RANDOM):
    plat, plon = random_pole()
    rlats, rlons = trace_great_circle(plat, plon, step_deg=3.0)  # coarser for speed
    relevs = get_elev_batch(rlats, rlons)
    land_mask = relevs > 0
    if land_mask.sum() < 10:
        continue
    rl, ro = rlats[land_mask], rlons[land_mask]
    # Sample up to 30 land points for speed
    if len(rl) > 30:
        idx = np.random.choice(len(rl), 30, replace=False)
        rl, ro = rl[idx], ro[idx]
    grads = [compute_gradient(la, lo) for la, lo in zip(rl, ro)]
    random_gradients.append(np.mean(grads))
    if (trial+1) % 50 == 0:
        print(f"  {trial+1}/{N_RANDOM} ({time.time()-t0:.0f}s)")

gc_gradient_pct = percentileofscore(random_gradients, gc_mean_gradient)
print(f"GC gradient percentile: {gc_gradient_pct:.1f}%")

# Gradient at monument cluster locations
monument_clusters = {
    'Memphis/Giza': (29.88, 31.13),
    'Nazca': (-14.73, -75.13),
    'Easter Island': (-27.12, -109.37),
    'Persepolis': (29.93, 52.89),
    'Angkor': (13.41, 103.87),
    'Mohenjo-daro': (27.32, 68.14),
}
print("\nGradients at monument clusters:")
cluster_gradients = {}
for name, (la, lo) in monument_clusters.items():
    g = compute_gradient(la, lo)
    cluster_gradients[name] = g
    print(f"  {name}: {g:.2f}")

# Save Test A results
test_a = {
    'gc_mean_gradient': gc_mean_gradient,
    'gc_gradient_percentile': gc_gradient_pct,
    'n_random_circles': len(random_gradients),
    'random_mean': float(np.mean(random_gradients)),
    'random_std': float(np.std(random_gradients)),
    'cluster_gradients': cluster_gradients,
    'gc_gradient_values': [float(g) for g in gc_gradients],
}
with open(OUT / 'terrain_transition.json', 'w') as f:
    json.dump(test_a, f, indent=2)

# Plot: gradient along the circle
fig, ax = plt.subplots(figsize=(14, 4))
ax.plot(range(len(gc_gradients)), gc_gradients, 'k-', alpha=0.3, linewidth=0.5)
# Running mean
window = 10
if len(gc_gradients) >= window:
    rm = np.convolve(gc_gradients, np.ones(window)/window, mode='valid')
    ax.plot(range(window//2, window//2+len(rm)), rm, 'b-', linewidth=1.5, label='10-pt running mean')
ax.axhline(gc_mean_gradient, color='red', linestyle='--', label=f'GC mean={gc_mean_gradient:.1f}')
ax.axhline(np.mean(random_gradients), color='gray', linestyle=':', label=f'Random mean={np.mean(random_gradients):.1f}')
ax.set_xlabel('Position along GC (land points)')
ax.set_ylabel('Terrain gradient magnitude')
ax.set_title(f'Test A: Terrain Transition Along Great Circle (pct={gc_gradient_pct:.0f}%)')
ax.legend()
fig.tight_layout()
fig.savefig(OUT / 'terrain_profile.png', dpi=150)
plt.close()
print(f"\nTest A complete. Percentile: {gc_gradient_pct:.1f}%")


# ================================================================
# TEST B: VIEWSHED / VISIBILITY (TOPOGRAPHIC PROMINENCE)
# ================================================================
print(f"\n{'='*70}")
print("TEST B: VIEWSHED / VISIBILITY (TOPOGRAPHIC PROMINENCE)")
print(f"{'='*70}")

def compute_prominence(lat, lon, radius_km=25):
    """Topographic prominence: elevation minus mean surrounding elevation."""
    center_elev = get_elev(lat, lon)
    if center_elev <= 0:
        return np.nan  # ocean
    # Sample surrounding in a circle
    n_samples = 64
    angles = np.linspace(0, 2*np.pi, n_samples, endpoint=False)
    dlat_per_km = 1.0 / 111.0
    dlon_per_km = 1.0 / (111.0 * max(np.cos(np.radians(lat)), 0.01))
    surr_elevs = []
    for a in angles:
        for r_frac in [0.33, 0.67, 1.0]:
            r = radius_km * r_frac
            slat = lat + r * np.sin(a) * dlat_per_km
            slon = lon + r * np.cos(a) * dlon_per_km
            e = get_elev(slat, slon)
            if e > -100:  # allow near-coast
                surr_elevs.append(max(e, 0))
    if len(surr_elevs) < 10:
        return np.nan
    return center_elev - np.mean(surr_elevs)

# Prominence along GC land points
print("Computing prominence along GC land points...")
t0 = time.time()
gc_prom = []
for i, (la, lo) in enumerate(zip(gc_land_lats, gc_land_lons)):
    p = compute_prominence(la, lo)
    gc_prom.append(p)
    if (i+1) % 20 == 0:
        print(f"  {i+1}/{len(gc_land_lats)} ({time.time()-t0:.0f}s)")
gc_prom = np.array(gc_prom)
gc_prom_valid = gc_prom[~np.isnan(gc_prom)]
gc_mean_prom = float(np.mean(gc_prom_valid))
print(f"GC mean prominence: {gc_mean_prom:.1f}m (n={len(gc_prom_valid)})")

# Control: offset 500km perpendicular
print("\nComputing prominence for 500km offset control points...")

def offset_perpendicular(lat, lon, dist_km, pole_lat, pole_lon):
    """Move point ~dist_km perpendicular to the GC (toward/away from pole)."""
    # Direction toward pole
    dlat = pole_lat - lat
    dlon = pole_lon - lon
    # Normalize
    d = haversine_km(lat, lon, pole_lat, pole_lon)
    if d < 1:
        return lat, lon
    # Move dist_km toward the pole
    frac = dist_km / d
    new_lat = lat + dlat * frac
    new_lon = lon + dlon * frac
    # Clamp
    new_lat = max(-89, min(89, new_lat))
    return new_lat, new_lon

control_prom = []
for la, lo in zip(gc_land_lats, gc_land_lons):
    ola, olo = offset_perpendicular(la, lo, 500, POLE_LAT, POLE_LON)
    if is_land(ola, olo):
        p = compute_prominence(ola, olo)
        control_prom.append(p)
    else:
        control_prom.append(np.nan)
control_prom = np.array(control_prom)
control_valid = control_prom[~np.isnan(control_prom)]
ctrl_mean_prom = float(np.mean(control_valid)) if len(control_valid) > 0 else 0
print(f"Control mean prominence: {ctrl_mean_prom:.1f}m (n={len(control_valid)})")

# Mann-Whitney test
if len(gc_prom_valid) > 5 and len(control_valid) > 5:
    stat, p_val = mannwhitneyu(gc_prom_valid, control_valid, alternative='greater')
    print(f"Mann-Whitney U: stat={stat:.0f}, p={p_val:.4f}")
else:
    stat, p_val = np.nan, np.nan

# Random circle comparison (prominence)
print(f"\nComputing prominence for {N_RANDOM} random circles...")
random_prom_means = []
for trial in range(N_RANDOM):
    plat, plon = random_pole()
    rlats, rlons = trace_great_circle(plat, plon, step_deg=3.0)
    relevs = get_elev_batch(rlats, rlons)
    land_mask = relevs > 0
    if land_mask.sum() < 10:
        continue
    rl, ro = rlats[land_mask], rlons[land_mask]
    if len(rl) > 25:
        idx = np.random.choice(len(rl), 25, replace=False)
        rl, ro = rl[idx], ro[idx]
    proms = [compute_prominence(la, lo) for la, lo in zip(rl, ro)]
    proms = [p for p in proms if not np.isnan(p)]
    if proms:
        random_prom_means.append(np.mean(proms))
    if (trial+1) % 50 == 0:
        print(f"  {trial+1}/{N_RANDOM} ({time.time()-t0:.0f}s)")

gc_prom_pct = percentileofscore(random_prom_means, gc_mean_prom)
print(f"GC prominence percentile: {gc_prom_pct:.1f}%")

# Save Test B
test_b = {
    'gc_mean_prominence_m': gc_mean_prom,
    'gc_n_land_points': int(len(gc_prom_valid)),
    'control_mean_prominence_m': ctrl_mean_prom,
    'control_n_points': int(len(control_valid)),
    'mann_whitney_stat': float(stat) if not np.isnan(stat) else None,
    'mann_whitney_p': float(p_val) if not np.isnan(p_val) else None,
    'gc_prominence_percentile_vs_random': gc_prom_pct,
    'n_random_circles': len(random_prom_means),
}
with open(OUT / 'visibility_analysis.json', 'w') as f:
    json.dump(test_b, f, indent=2)

# Box plot
fig, ax = plt.subplots(figsize=(8, 5))
data = [gc_prom_valid, control_valid]
labels = [f'GC (μ={gc_mean_prom:.0f}m)', f'Control 500km offset (μ={ctrl_mean_prom:.0f}m)']
bp = ax.boxplot(data, labels=labels, patch_artist=True)
bp['boxes'][0].set_facecolor('#4488cc')
bp['boxes'][1].set_facecolor('#cc8844')
ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax.set_ylabel('Topographic Prominence (m)')
ax.set_title(f'Test B: Topographic Prominence — GC vs Control\n(pct vs random={gc_prom_pct:.0f}%, MW p={p_val:.4f})')
fig.tight_layout()
fig.savefig(OUT / 'prominence_comparison.png', dpi=150)
plt.close()


# ================================================================
# TEST B EXTENSION: MONUMENT VS SETTLEMENT PROMINENCE
# ================================================================
print(f"\n{'='*70}")
print("TEST B (ext): MONUMENT VS SETTLEMENT PROMINENCE (Pleiades)")
print(f"{'='*70}")

# Load Pleiades
pleiades = pd.read_csv(BASE / "data" / "pleiades" / "pleiades-places-latest.csv")
pleiades = pleiades.dropna(subset=['reprLat', 'reprLong'])

monument_types = {'temple', 'temple-2', 'sanctuary', 'pyramid', 'monument',
                  'amphitheatre', 'aqueduct', 'tumulus', 'tomb', 'theatre',
                  'bath', 'church', 'church-2', 'architecturalcomplex'}
settlement_types = {'settlement', 'villa', 'settlement-modern', 'fort', 'fort-2',
                    'station', 'port', 'farm', 'city', 'village', 'town'}

def classify(ft_str):
    if pd.isna(ft_str):
        return None
    types = {t.strip() for t in ft_str.split(',')}
    is_m = bool(types & monument_types)
    is_s = bool(types & settlement_types)
    if is_m and not is_s:
        return 'monument'
    elif is_s and not is_m:
        return 'settlement'
    return None

pleiades['site_class'] = pleiades['featureTypes'].apply(classify)
pleiades = pleiades.dropna(subset=['site_class'])

# GC distance for all Pleiades
pole_lat_r = np.radians(POLE_LAT)
pole_lon_r = np.radians(POLE_LON)
lats_r = np.radians(pleiades['reprLat'].values)
lons_r = np.radians(pleiades['reprLong'].values)
cos_ang = (np.sin(pole_lat_r)*np.sin(lats_r) +
           np.cos(pole_lat_r)*np.cos(lats_r)*np.cos(lons_r - pole_lon_r))
ang = np.arccos(np.clip(cos_ang, -1, 1))
pleiades['gc_dist_km'] = np.abs(ang - np.pi/2) * R

# Filter to within 100km of GC
near_gc = pleiades[pleiades['gc_dist_km'] < 100].copy()
print(f"Pleiades within 100km of GC: {len(near_gc)} ({(near_gc['site_class']=='monument').sum()} monuments, {(near_gc['site_class']=='settlement').sum()} settlements)")

# Compute prominence for near-GC sites (cap at 500 for speed)
mon_near = near_gc[near_gc['site_class'] == 'monument']
set_near = near_gc[near_gc['site_class'] == 'settlement']

# Sample if too many
if len(mon_near) > 250:
    mon_near = mon_near.sample(250, random_state=42)
if len(set_near) > 250:
    set_near = set_near.sample(250, random_state=42)

print(f"Computing prominence for {len(mon_near)} monuments and {len(set_near)} settlements...")
mon_proms = []
for _, row in mon_near.iterrows():
    p = compute_prominence(row['reprLat'], row['reprLong'], radius_km=10)
    if not np.isnan(p):
        mon_proms.append(p)

set_proms = []
for _, row in set_near.iterrows():
    p = compute_prominence(row['reprLat'], row['reprLong'], radius_km=10)
    if not np.isnan(p):
        set_proms.append(p)

print(f"Monuments: μ={np.mean(mon_proms):.1f}m (n={len(mon_proms)})")
print(f"Settlements: μ={np.mean(set_proms):.1f}m (n={len(set_proms)})")

if len(mon_proms) > 5 and len(set_proms) > 5:
    ms_stat, ms_p = mannwhitneyu(mon_proms, set_proms, alternative='greater')
    print(f"Mann-Whitney (monuments > settlements): stat={ms_stat:.0f}, p={ms_p:.4f}")
else:
    ms_stat, ms_p = np.nan, np.nan

test_b['monument_mean_prominence_m'] = float(np.mean(mon_proms)) if mon_proms else None
test_b['settlement_mean_prominence_m'] = float(np.mean(set_proms)) if set_proms else None
test_b['monument_n'] = len(mon_proms)
test_b['settlement_n'] = len(set_proms)
test_b['mon_vs_set_mann_whitney_stat'] = float(ms_stat) if not np.isnan(ms_stat) else None
test_b['mon_vs_set_mann_whitney_p'] = float(ms_p) if not np.isnan(ms_p) else None
with open(OUT / 'visibility_analysis.json', 'w') as f:
    json.dump(test_b, f, indent=2)

# Box plot monuments vs settlements
fig, ax = plt.subplots(figsize=(7, 5))
data = [mon_proms, set_proms]
labels = [f'Monuments (μ={np.mean(mon_proms):.0f}m, n={len(mon_proms)})',
          f'Settlements (μ={np.mean(set_proms):.0f}m, n={len(set_proms)})']
bp = ax.boxplot(data, labels=labels, patch_artist=True)
bp['boxes'][0].set_facecolor('#cc4444')
bp['boxes'][1].set_facecolor('#4488cc')
ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax.set_ylabel('Topographic Prominence (m)')
p_str = f'{ms_p:.2e}' if ms_p < 0.01 else f'{ms_p:.4f}'
ax.set_title(f'Test B: Monument vs Settlement Prominence (within 100km of GC)\nMann-Whitney p={p_str}')
fig.tight_layout()
fig.savefig(OUT / 'monument_vs_settlement_prominence.png', dpi=150)
plt.close()


# ================================================================
# TEST C: NAVIGATIONAL WAYPOINT ANALYSIS
# ================================================================
print(f"\n{'='*70}")
print("TEST C: NAVIGATIONAL WAYPOINT ANALYSIS")
print(f"{'='*70}")

# Load coastline crossings
with open(BASE / "outputs" / "coastline_intersections" / "crossings.json") as f:
    crossings_data = json.load(f)
crossings = crossings_data['alison_crossings']
print(f"Loaded {len(crossings)} coastline crossings")

# Method 1: Coastal Prominence at Crossing Points
print("\nMethod 1: Prominence at GC coast crossings vs random coastal points...")
crossing_proms = []
for c in crossings:
    p = compute_prominence(c['lat'], c['lon'], radius_km=50)
    crossing_proms.append({'lat': c['lat'], 'lon': c['lon'],
                           'prominence': float(p) if not np.isnan(p) else None,
                           'is_continental': bool(c.get('is_continental', False))})

crossing_prom_vals = [c['prominence'] for c in crossing_proms if c['prominence'] is not None]
print(f"GC crossing prominence: μ={np.mean(crossing_prom_vals):.1f}m (n={len(crossing_prom_vals)})")

# Generate random coastal points for comparison
# Use a simple approach: find coastal points (land adjacent to ocean) at similar latitudes
print("Generating random coastal reference points...")
random_coast_proms = []
crossing_lats_used = [c['lat'] for c in crossings]
for lat_ref in crossing_lats_used:
    # Scan along the latitude for coast points
    for lon_try in np.random.uniform(-180, 180, 50):
        e = get_elev(lat_ref, lon_try)
        if 0 < e < 200:  # likely coastal
            # Check neighbor is ocean
            e_east = get_elev(lat_ref, lon_try + 0.1)
            e_west = get_elev(lat_ref, lon_try - 0.1)
            if e_east <= 0 or e_west <= 0:
                p = compute_prominence(lat_ref, lon_try, radius_km=50)
                if not np.isnan(p):
                    random_coast_proms.append(p)
                    if len(random_coast_proms) >= len(crossing_prom_vals) * 10:
                        break

print(f"Random coastal prominence: μ={np.mean(random_coast_proms):.1f}m (n={len(random_coast_proms)})")

if len(crossing_prom_vals) > 5 and len(random_coast_proms) > 5:
    coast_stat, coast_p = mannwhitneyu(crossing_prom_vals, random_coast_proms, alternative='greater')
    print(f"Mann-Whitney (GC crossings > random coast): stat={coast_stat:.0f}, p={coast_p:.4f}")
else:
    coast_stat, coast_p = np.nan, np.nan

# Method 2: Inter-Visibility Chain
print("\nMethod 2: Inter-visibility chain analysis...")
# Filter to continental crossings only and sort by position
continental_crossings = [c for c in crossings if c.get('is_continental', False)]
continental_crossings.sort(key=lambda x: x.get('arc_position', 0))

intervis_results = []
for i in range(len(continental_crossings) - 1):
    c1 = continental_crossings[i]
    c2 = continental_crossings[i+1]
    dist = haversine_km(c1['lat'], c1['lon'], c2['lat'], c2['lon'])

    # Find max elevation along the segment between crossings
    n_sample = max(10, int(dist / 50))  # sample every ~50km
    seg_lats = np.linspace(c1['lat'], c2['lat'], n_sample)
    seg_lons = np.linspace(c1['lon'], c2['lon'], n_sample)
    seg_elevs = get_elev_batch(seg_lats, seg_lons)
    max_elev = float(np.max(np.maximum(seg_elevs, 0)))

    # Theoretical line-of-sight distance from max elevation
    los_km = 3.57 * np.sqrt(max(max_elev, 0)) if max_elev > 0 else 0
    can_see = los_km > dist / 2  # from midpoint high point can see both ends?

    intervis_results.append({
        'from': [c1['lat'], c1['lon']],
        'to': [c2['lat'], c2['lon']],
        'distance_km': dist,
        'max_elevation_m': max_elev,
        'theoretical_los_km': los_km,
        'potentially_intervisible': bool(can_see),
    })

n_visible = sum(1 for r in intervis_results if r['potentially_intervisible'])
print(f"Inter-visibility: {n_visible}/{len(intervis_results)} segments potentially intervisible")

# Method 3: Easter Island visibility
print("\nMethod 3: Easter Island visibility check...")
easter_elev = get_elev(-27.12, -109.37)
easter_vis_km = 3.57 * np.sqrt(max(easter_elev, 0))
print(f"Easter Island max elevation: {easter_elev:.0f}m → visible from {easter_vis_km:.0f}km")

# Save Test C
test_c = {
    'crossing_mean_prominence_m': float(np.mean(crossing_prom_vals)),
    'crossing_n': len(crossing_prom_vals),
    'random_coast_mean_prominence_m': float(np.mean(random_coast_proms)) if random_coast_proms else None,
    'random_coast_n': len(random_coast_proms),
    'coast_mann_whitney_stat': float(coast_stat) if not np.isnan(coast_stat) else None,
    'coast_mann_whitney_p': float(coast_p) if not np.isnan(coast_p) else None,
    'intervisibility_results': intervis_results,
    'n_intervisible_segments': n_visible,
    'n_total_segments': len(intervis_results),
    'easter_island_elevation_m': float(easter_elev),
    'easter_island_visibility_km': float(easter_vis_km),
    'crossing_details': crossing_proms,
}
with open(OUT / 'navigational_analysis.json', 'w') as f:
    json.dump(test_c, f, indent=2)

# Plot crossing prominence
fig, ax = plt.subplots(figsize=(8, 5))
if crossing_prom_vals and random_coast_proms:
    bp = ax.boxplot([crossing_prom_vals, random_coast_proms],
                    labels=[f'GC crossings (n={len(crossing_prom_vals)})',
                            f'Random coast (n={len(random_coast_proms)})'],
                    patch_artist=True)
    bp['boxes'][0].set_facecolor('#44aa44')
    bp['boxes'][1].set_facecolor('#aaaaaa')
    ax.set_ylabel('Topographic Prominence (m)')
    p_str = f'{coast_p:.2e}' if coast_p < 0.01 else f'{coast_p:.4f}'
    ax.set_title(f'Test C: Coastal Crossing Prominence\nMann-Whitney p={p_str}')
fig.tight_layout()
fig.savefig(OUT / 'crossing_prominence.png', dpi=150)
plt.close()

# Save intervisibility
with open(OUT / 'intervisibility_chain.json', 'w') as f:
    json.dump(intervis_results, f, indent=2)


# ================================================================
# TEST D: NILE VALLEY CONSTRICTION AT MEMPHIS
# ================================================================
print(f"\n{'='*70}")
print("TEST D: NILE VALLEY CONSTRICTION AT MEMPHIS")
print(f"{'='*70}")

# Compute valley width at different latitudes
valley_widths = []
elev_threshold = 50  # meters — escarpment edge
for lat in np.arange(29.0, 30.5, 0.05):
    # East-west profile at this latitude
    lons_profile = np.arange(30.5, 31.8, 0.005)
    elevs = []
    for lo in lons_profile:
        elevs.append(get_elev(lat, lo))
    elevs = np.array(elevs)

    # Find western escarpment (first point > threshold from west)
    west_idx = None
    east_idx = None
    for j in range(len(elevs)):
        if elevs[j] < elev_threshold:
            if west_idx is None:
                west_idx = j
            east_idx = j

    if west_idx is not None and east_idx is not None and east_idx > west_idx:
        width_km = haversine_km(lat, lons_profile[west_idx], lat, lons_profile[east_idx])
        valley_widths.append({
            'lat': float(lat),
            'width_km': float(width_km),
            'west_lon': float(lons_profile[west_idx]),
            'east_lon': float(lons_profile[east_idx]),
        })
    else:
        valley_widths.append({
            'lat': float(lat),
            'width_km': float('nan'),
            'west_lon': None,
            'east_lon': None,
        })

# Where does the GC cross the Nile?
# Find GC latitude at Nile longitude (~31.2°E)
# From traced GC, find points near lon=31.2
gc_nile_mask = (gc_lons > 30.5) & (gc_lons < 32.0) & (gc_lats > 28) & (gc_lats < 31)
if gc_nile_mask.any():
    gc_nile_lat = float(gc_lats[gc_nile_mask].mean())
    gc_nile_lon = float(gc_lons[gc_nile_mask].mean())
else:
    gc_nile_lat = 29.88  # approximate from directive
    gc_nile_lon = 31.13

print(f"GC crosses Nile region at approximately {gc_nile_lat:.2f}°N, {gc_nile_lon:.2f}°E")

# Find minimum valley width
valid_widths = [v for v in valley_widths if not np.isnan(v['width_km'])]
if valid_widths:
    min_width = min(valid_widths, key=lambda x: x['width_km'])
    print(f"Minimum valley width: {min_width['width_km']:.1f}km at lat={min_width['lat']:.2f}°N")

    # Width at GC crossing
    gc_width = None
    for v in valley_widths:
        if abs(v['lat'] - gc_nile_lat) < 0.03:
            gc_width = v
            break
    if gc_width:
        print(f"Valley width at GC crossing ({gc_nile_lat:.2f}°N): {gc_width['width_km']:.1f}km")

    # Memphis location
    memphis_lat = 29.85
    memphis_width = None
    for v in valley_widths:
        if abs(v['lat'] - memphis_lat) < 0.03:
            memphis_width = v
            break
    if memphis_width:
        print(f"Valley width at Memphis ({memphis_lat:.2f}°N): {memphis_width['width_km']:.1f}km")

    print(f"\nGC crossing lat ({gc_nile_lat:.2f}°N) vs minimum width lat ({min_width['lat']:.2f}°N): "
          f"offset = {abs(gc_nile_lat - min_width['lat']):.2f}° ({abs(gc_nile_lat - min_width['lat'])*111:.0f}km)")
    print(f"Memphis lat ({memphis_lat:.2f}°N) vs minimum width lat ({min_width['lat']:.2f}°N): "
          f"offset = {abs(memphis_lat - min_width['lat']):.2f}° ({abs(memphis_lat - min_width['lat'])*111:.0f}km)")

# Save Test D
test_d = {
    'gc_nile_crossing_lat': gc_nile_lat,
    'gc_nile_crossing_lon': gc_nile_lon,
    'memphis_lat': memphis_lat,
    'valley_widths': valley_widths,
    'min_width_lat': min_width['lat'] if valid_widths else None,
    'min_width_km': min_width['width_km'] if valid_widths else None,
    'gc_crossing_width_km': gc_width['width_km'] if gc_width else None,
    'memphis_width_km': memphis_width['width_km'] if memphis_width else None,
}
with open(OUT / 'nile_valley_width.json', 'w') as f:
    json.dump(test_d, f, indent=2)

# Plot: valley width vs latitude
fig, ax = plt.subplots(figsize=(8, 6))
valid_lats = [v['lat'] for v in valid_widths]
valid_ws = [v['width_km'] for v in valid_widths]
ax.plot(valid_lats, valid_ws, 'b-o', markersize=3, label='Valley width')
if min_width:
    ax.axvline(min_width['lat'], color='green', linestyle='--', alpha=0.7,
               label=f'Min width ({min_width["lat"]:.2f}°N, {min_width["width_km"]:.1f}km)')
ax.axvline(gc_nile_lat, color='red', linestyle='--', alpha=0.7,
           label=f'GC crossing ({gc_nile_lat:.2f}°N)')
ax.axvline(memphis_lat, color='orange', linestyle=':', alpha=0.7,
           label=f'Memphis ({memphis_lat:.2f}°N)')
ax.set_xlabel('Latitude (°N)')
ax.set_ylabel('Valley Width (km)')
ax.set_title('Test D: Nile Valley Width Profile')
ax.legend(fontsize=8)
ax.invert_xaxis()  # North at left
fig.tight_layout()
fig.savefig(OUT / 'nile_valley_profile.png', dpi=150)
plt.close()
print("Test D complete.")


# ================================================================
# TEST E: MONUMENT ELEVATION DIFFERENTIAL (TPI)
# ================================================================
print(f"\n{'='*70}")
print("TEST E: MONUMENT ELEVATION DIFFERENTIAL (TPI)")
print(f"{'='*70}")

def compute_tpi(lat, lon, radius_km=10):
    """Topographic Position Index: elevation minus mean within radius."""
    site_elev = get_elev(lat, lon)
    if site_elev <= -100:
        return np.nan
    n_samples = 48
    angles = np.linspace(0, 2*np.pi, n_samples, endpoint=False)
    dlat_per_km = 1.0 / 111.0
    dlon_per_km = 1.0 / (111.0 * max(np.cos(np.radians(lat)), 0.01))
    surr = []
    for a in angles:
        for r_frac in [0.5, 1.0]:
            r = radius_km * r_frac
            slat = lat + r * np.sin(a) * dlat_per_km
            slon = lon + r * np.cos(a) * dlon_per_km
            e = get_elev(slat, slon)
            surr.append(e)
    return site_elev - np.mean(surr)

# Use already-filtered near-GC Pleiades (within 100km)
print(f"Computing TPI for near-GC monuments and settlements...")
# Re-filter from full near_gc set (before sampling)
near_gc_full = pleiades[pleiades['gc_dist_km'] < 100].copy()
mon_all = near_gc_full[near_gc_full['site_class'] == 'monument']
set_all = near_gc_full[near_gc_full['site_class'] == 'settlement']
if len(mon_all) > 300:
    mon_all = mon_all.sample(300, random_state=42)
if len(set_all) > 300:
    set_all = set_all.sample(300, random_state=42)

print(f"  Monuments: {len(mon_all)}, Settlements: {len(set_all)}")

mon_tpi = []
for _, row in mon_all.iterrows():
    t = compute_tpi(row['reprLat'], row['reprLong'])
    if not np.isnan(t):
        mon_tpi.append(t)

set_tpi = []
for _, row in set_all.iterrows():
    t = compute_tpi(row['reprLat'], row['reprLong'])
    if not np.isnan(t):
        set_tpi.append(t)

print(f"Monument TPI: μ={np.mean(mon_tpi):.1f}m, median={np.median(mon_tpi):.1f}m (n={len(mon_tpi)})")
print(f"Settlement TPI: μ={np.mean(set_tpi):.1f}m, median={np.median(set_tpi):.1f}m (n={len(set_tpi)})")

if len(mon_tpi) > 5 and len(set_tpi) > 5:
    tpi_stat, tpi_p = mannwhitneyu(mon_tpi, set_tpi, alternative='greater')
    print(f"Mann-Whitney (monuments > settlements): stat={tpi_stat:.0f}, p={tpi_p:.4f}")
else:
    tpi_stat, tpi_p = np.nan, np.nan

# Save Test E
test_e = {
    'monument_mean_tpi_m': float(np.mean(mon_tpi)),
    'monument_median_tpi_m': float(np.median(mon_tpi)),
    'monument_n': len(mon_tpi),
    'settlement_mean_tpi_m': float(np.mean(set_tpi)),
    'settlement_median_tpi_m': float(np.median(set_tpi)),
    'settlement_n': len(set_tpi),
    'mann_whitney_stat': float(tpi_stat) if not np.isnan(tpi_stat) else None,
    'mann_whitney_p': float(tpi_p) if not np.isnan(tpi_p) else None,
}
with open(OUT / 'elevation_differential.json', 'w') as f:
    json.dump(test_e, f, indent=2)

# Histogram/box plot
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Box plot
ax = axes[0]
bp = ax.boxplot([mon_tpi, set_tpi],
                labels=[f'Monuments (n={len(mon_tpi)})', f'Settlements (n={len(set_tpi)})'],
                patch_artist=True)
bp['boxes'][0].set_facecolor('#cc4444')
bp['boxes'][1].set_facecolor('#4488cc')
ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax.set_ylabel('TPI (m)')
p_str = f'{tpi_p:.2e}' if tpi_p < 0.01 else f'{tpi_p:.4f}'
ax.set_title(f'TPI: Monuments vs Settlements\nMW p={p_str}')

# Histogram
ax = axes[1]
bins = np.linspace(min(min(mon_tpi), min(set_tpi)), max(max(mon_tpi), max(set_tpi)), 40)
ax.hist(mon_tpi, bins=bins, alpha=0.6, color='#cc4444', label=f'Monuments (μ={np.mean(mon_tpi):.0f}m)')
ax.hist(set_tpi, bins=bins, alpha=0.6, color='#4488cc', label=f'Settlements (μ={np.mean(set_tpi):.0f}m)')
ax.axvline(0, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('TPI (m)')
ax.set_ylabel('Count')
ax.legend()
ax.set_title('TPI Distribution')

fig.suptitle('Test E: Monument Elevation Differential', fontsize=13, y=1.02)
fig.tight_layout()
fig.savefig(OUT / 'tpi_comparison.png', dpi=150)
plt.close()


# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'='*70}")
print("SUMMARY OF ALL TESTS")
print(f"{'='*70}")

summary = {
    'test_a_terrain_transition': {
        'gc_mean_gradient': gc_mean_gradient,
        'gc_percentile': gc_gradient_pct,
        'verdict': 'SUPPORTED' if gc_gradient_pct > 80 else ('MARGINAL' if gc_gradient_pct > 60 else 'NOT SUPPORTED'),
    },
    'test_b_visibility': {
        'gc_mean_prominence_m': gc_mean_prom,
        'gc_percentile': gc_prom_pct,
        'gc_vs_control_p': float(p_val) if not np.isnan(p_val) else None,
        'monument_prominence_m': float(np.mean(mon_proms)) if mon_proms else None,
        'settlement_prominence_m': float(np.mean(set_proms)) if set_proms else None,
        'mon_vs_set_p': float(ms_p) if not np.isnan(ms_p) else None,
    },
    'test_c_navigation': {
        'crossing_prominence_m': float(np.mean(crossing_prom_vals)),
        'random_coast_prominence_m': float(np.mean(random_coast_proms)) if random_coast_proms else None,
        'crossing_vs_coast_p': float(coast_p) if not np.isnan(coast_p) else None,
        'intervisible_segments': f'{n_visible}/{len(intervis_results)}',
    },
    'test_d_nile_constriction': {
        'gc_crossing_lat': gc_nile_lat,
        'min_width_lat': min_width['lat'] if valid_widths else None,
        'memphis_lat': memphis_lat,
        'offset_gc_to_min_km': abs(gc_nile_lat - min_width['lat'])*111 if valid_widths else None,
    },
    'test_e_elevation_differential': {
        'monument_mean_tpi_m': float(np.mean(mon_tpi)),
        'settlement_mean_tpi_m': float(np.mean(set_tpi)),
        'mann_whitney_p': float(tpi_p) if not np.isnan(tpi_p) else None,
        'verdict': 'SUPPORTED' if (not np.isnan(tpi_p) and tpi_p < 0.05) else 'NOT SUPPORTED',
    },
}

with open(OUT / 'summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print(json.dumps(summary, indent=2))

# Interpretation table
print(f"\n{'='*70}")
print("INTERPRETATION")
print(f"{'='*70}")
print(f"\nTest A — Terrain Transition Zone:")
print(f"  GC gradient percentile: {gc_gradient_pct:.0f}%")
print(f"  {'→ Circle traces high-gradient terrain boundaries' if gc_gradient_pct > 70 else '→ Circle does NOT preferentially trace terrain transitions'}")

print(f"\nTest B — Visibility:")
print(f"  GC prominence percentile: {gc_prom_pct:.0f}%")
print(f"  GC vs control offset p-value: {p_val:.4f}" if not np.isnan(p_val) else "  GC vs control: insufficient data")
if mon_proms and set_proms:
    print(f"  Monument prominence: {np.mean(mon_proms):.0f}m vs Settlement: {np.mean(set_proms):.0f}m (p={ms_p:.4f})")

print(f"\nTest C — Navigation:")
print(f"  Crossing prominence vs random coast p-value: {coast_p:.4f}" if not np.isnan(coast_p) else "  Insufficient data")
print(f"  Inter-visible segments: {n_visible}/{len(intervis_results)}")

print(f"\nTest D — Nile Constriction:")
print(f"  GC crossing at {gc_nile_lat:.2f}°N")
if valid_widths:
    print(f"  Minimum width at {min_width['lat']:.2f}°N")
    print(f"  Offset: {abs(gc_nile_lat - min_width['lat'])*111:.0f}km")
    print(f"  {'→ GC crosses near the valley constriction' if abs(gc_nile_lat - min_width['lat']) < 0.3 else '→ GC does NOT cross at the constriction'}")

print(f"\nTest E — Monument TPI:")
if not np.isnan(tpi_p):
    print(f"  Monument TPI: {np.mean(mon_tpi):.0f}m vs Settlement TPI: {np.mean(set_tpi):.0f}m")
    print(f"  p-value: {tpi_p:.4f}")
    print(f"  {'→ Monuments are on systematically HIGHER ground' if tpi_p < 0.05 else '→ No significant elevation difference'}")

print(f"\nAll outputs saved to: {OUT}")
print("Done.")
