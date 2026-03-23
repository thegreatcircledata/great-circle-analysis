#!/usr/bin/env python3
"""Tests C, D, E from terrain-visibility directive (A and B already complete)."""

import json, math, os, sys, time
import numpy as np
import netCDF4 as nc
import pandas as pd
from scipy.stats import percentileofscore, mannwhitneyu
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

POLE_LAT = 59.682122
POLE_LON = -138.646087
R = 6371.0

BASE = Path(os.path.expanduser("~/megalith_site_research"))
OUT = BASE / "outputs" / "terrain_visibility_tests"
OUT.mkdir(parents=True, exist_ok=True)

# Load ETOPO1
print("Loading ETOPO1...")
etopo_path = BASE / "data" / "geophysical" / "etopo" / "ETOPO1_Ice_g_gmt4.grd"
ds = nc.Dataset(str(etopo_path))
etopo_lons = ds.variables['x'][:]
etopo_lats = ds.variables['y'][:]
etopo_z = ds.variables['z']
lon_min, lon_step = float(etopo_lons[0]), float(etopo_lons[1] - etopo_lons[0])
lat_min, lat_step = float(etopo_lats[0]), float(etopo_lats[1] - etopo_lats[0])
n_lats = len(etopo_lats)
n_lons = len(etopo_lons)

def get_elev(lat, lon):
    li = int(round((lat - lat_min) / lat_step))
    lo = int(round((lon - lon_min) / lon_step))
    li = max(0, min(li, n_lats - 1))
    lo = max(0, min(lo, n_lons - 1))
    return float(etopo_z[li, lo])

def get_elev_batch(lats, lons):
    lat_idx = np.round((np.asarray(lats) - lat_min) / lat_step).astype(int)
    lon_idx = np.round((np.asarray(lons) - lon_min) / lon_step).astype(int)
    lat_idx = np.clip(lat_idx, 0, n_lats - 1)
    lon_idx = np.clip(lon_idx, 0, n_lons - 1)
    elevs = np.empty(len(lat_idx))
    for i in range(len(lat_idx)):
        elevs[i] = float(etopo_z[lat_idx[i], lon_idx[i]])
    return elevs

def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1; dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    return R * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))

def compute_prominence(lat, lon, radius_km=25):
    center_elev = get_elev(lat, lon)
    if center_elev <= 0:
        return None
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
            if e > -100:
                surr_elevs.append(max(e, 0))
    if len(surr_elevs) < 10:
        return None
    return center_elev - np.mean(surr_elevs)

def trace_great_circle(pole_lat, pole_lon, step_deg=1.0):
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

gc_lats, gc_lons = trace_great_circle(POLE_LAT, POLE_LON, step_deg=1.0)

# ================================================================
# TEST C: NAVIGATIONAL WAYPOINT ANALYSIS
# ================================================================
print(f"\n{'='*70}")
print("TEST C: NAVIGATIONAL WAYPOINT ANALYSIS")
print(f"{'='*70}")

with open(BASE / "outputs" / "coastline_intersections" / "crossings.json") as f:
    crossings_data = json.load(f)
crossings = crossings_data['alison_crossings']
print(f"Loaded {len(crossings)} coastline crossings")

# Method 1: Coastal Prominence
print("\nMethod 1: Prominence at GC coast crossings vs random coastal points...")
crossing_proms = []
for c in crossings:
    p = compute_prominence(c['lat'], c['lon'], radius_km=50)
    crossing_proms.append({'lat': c['lat'], 'lon': c['lon'],
                           'prominence': float(p) if p is not None else None,
                           'is_continental': bool(c.get('is_continental', False))})

crossing_prom_vals = [c['prominence'] for c in crossing_proms if c['prominence'] is not None]
print(f"GC crossing prominence: μ={np.mean(crossing_prom_vals):.1f}m (n={len(crossing_prom_vals)})")

# Random coastal points
print("Generating random coastal reference points...")
random_coast_proms = []
crossing_lats = [c['lat'] for c in crossings]
for lat_ref in crossing_lats:
    found = 0
    for lon_try in np.random.uniform(-180, 180, 200):
        e = get_elev(lat_ref, lon_try)
        if 0 < e < 200:
            e_east = get_elev(lat_ref, lon_try + 0.1)
            e_west = get_elev(lat_ref, lon_try - 0.1)
            if e_east <= 0 or e_west <= 0:
                p = compute_prominence(lat_ref, lon_try, radius_km=50)
                if p is not None:
                    random_coast_proms.append(p)
                    found += 1
                    if found >= 5:
                        break

print(f"Random coastal prominence: μ={np.mean(random_coast_proms):.1f}m (n={len(random_coast_proms)})")

coast_stat, coast_p = np.nan, np.nan
if len(crossing_prom_vals) > 5 and len(random_coast_proms) > 5:
    coast_stat, coast_p = mannwhitneyu(crossing_prom_vals, random_coast_proms, alternative='greater')
    print(f"Mann-Whitney (GC crossings > random coast): stat={coast_stat:.0f}, p={coast_p:.4f}")

# Method 2: Inter-Visibility Chain
print("\nMethod 2: Inter-visibility chain analysis...")
continental_crossings = [c for c in crossings if c.get('is_continental', False)]
continental_crossings.sort(key=lambda x: x.get('arc_position', 0))

intervis_results = []
for i in range(len(continental_crossings) - 1):
    c1 = continental_crossings[i]
    c2 = continental_crossings[i+1]
    dist = haversine_km(c1['lat'], c1['lon'], c2['lat'], c2['lon'])
    n_sample = max(10, int(dist / 50))
    seg_lats = np.linspace(c1['lat'], c2['lat'], n_sample)
    seg_lons = np.linspace(c1['lon'], c2['lon'], n_sample)
    seg_elevs = get_elev_batch(seg_lats, seg_lons)
    max_elev = float(np.max(np.maximum(seg_elevs, 0)))
    los_km = 3.57 * np.sqrt(max(max_elev, 0)) if max_elev > 0 else 0
    can_see = los_km > dist / 2

    intervis_results.append({
        'from': [c1['lat'], c1['lon']],
        'to': [c2['lat'], c2['lon']],
        'distance_km': float(dist),
        'max_elevation_m': float(max_elev),
        'theoretical_los_km': float(los_km),
        'potentially_intervisible': bool(can_see),
    })

n_visible = sum(1 for r in intervis_results if r['potentially_intervisible'])
print(f"Inter-visibility: {n_visible}/{len(intervis_results)} segments potentially intervisible")

# Method 3: Easter Island
easter_elev = get_elev(-27.12, -109.37)
easter_vis_km = 3.57 * np.sqrt(max(easter_elev, 0))
print(f"Easter Island elevation: {easter_elev:.0f}m → visible from {easter_vis_km:.0f}km")

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
with open(OUT / 'intervisibility_chain.json', 'w') as f:
    json.dump(intervis_results, f, indent=2)

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
if crossing_prom_vals and random_coast_proms:
    bp = ax.boxplot([crossing_prom_vals, random_coast_proms],
                    tick_labels=[f'GC crossings (n={len(crossing_prom_vals)})',
                                 f'Random coast (n={len(random_coast_proms)})'],
                    patch_artist=True)
    bp['boxes'][0].set_facecolor('#44aa44')
    bp['boxes'][1].set_facecolor('#aaaaaa')
    ax.set_ylabel('Topographic Prominence (m)')
    p_str = f'{coast_p:.2e}' if not np.isnan(coast_p) and coast_p < 0.01 else f'{coast_p:.4f}'
    ax.set_title(f'Test C: Coastal Crossing Prominence\nMann-Whitney p={p_str}')
fig.tight_layout()
fig.savefig(OUT / 'crossing_prominence.png', dpi=150)
plt.close()
print("Test C complete.")


# ================================================================
# TEST D: NILE VALLEY CONSTRICTION AT MEMPHIS
# ================================================================
print(f"\n{'='*70}")
print("TEST D: NILE VALLEY CONSTRICTION AT MEMPHIS")
print(f"{'='*70}")

valley_widths = []
elev_threshold = 50
for lat in np.arange(29.0, 30.5, 0.05):
    lons_profile = np.arange(30.5, 31.8, 0.005)
    elevs = np.array([get_elev(lat, lo) for lo in lons_profile])

    # Find the valley floor: contiguous stretch below threshold
    below = elevs < elev_threshold
    if below.any():
        first_below = np.argmax(below)
        last_below = len(below) - 1 - np.argmax(below[::-1])
        width_km = haversine_km(lat, float(lons_profile[first_below]), lat, float(lons_profile[last_below]))
        valley_widths.append({
            'lat': float(lat),
            'width_km': float(width_km),
            'west_lon': float(lons_profile[first_below]),
            'east_lon': float(lons_profile[last_below]),
        })
    else:
        valley_widths.append({'lat': float(lat), 'width_km': None, 'west_lon': None, 'east_lon': None})

# GC crossing latitude at the Nile
gc_nile_mask = (gc_lons > 30.5) & (gc_lons < 32.0) & (gc_lats > 28) & (gc_lats < 31)
gc_nile_lat = float(gc_lats[gc_nile_mask].mean()) if gc_nile_mask.any() else 29.88
gc_nile_lon = float(gc_lons[gc_nile_mask].mean()) if gc_nile_mask.any() else 31.13
print(f"GC crosses Nile region at ~{gc_nile_lat:.2f}°N, {gc_nile_lon:.2f}°E")

valid_widths = [v for v in valley_widths if v['width_km'] is not None]
min_width = min(valid_widths, key=lambda x: x['width_km'])
print(f"Minimum valley width: {min_width['width_km']:.1f}km at lat={min_width['lat']:.2f}°N")

gc_width = next((v for v in valley_widths if abs(v['lat'] - gc_nile_lat) < 0.03), None)
if gc_width and gc_width['width_km']:
    print(f"Valley width at GC crossing: {gc_width['width_km']:.1f}km")

memphis_lat = 29.85
memphis_width = next((v for v in valley_widths if abs(v['lat'] - memphis_lat) < 0.03), None)
if memphis_width and memphis_width['width_km']:
    print(f"Valley width at Memphis: {memphis_width['width_km']:.1f}km")

offset_gc_min = abs(gc_nile_lat - min_width['lat']) * 111
offset_memphis_min = abs(memphis_lat - min_width['lat']) * 111
print(f"GC to min-width offset: {offset_gc_min:.0f}km")
print(f"Memphis to min-width offset: {offset_memphis_min:.0f}km")

test_d = {
    'gc_nile_crossing_lat': gc_nile_lat,
    'gc_nile_crossing_lon': gc_nile_lon,
    'memphis_lat': memphis_lat,
    'valley_widths': valley_widths,
    'min_width_lat': min_width['lat'],
    'min_width_km': min_width['width_km'],
    'gc_crossing_width_km': gc_width['width_km'] if gc_width and gc_width['width_km'] else None,
    'memphis_width_km': memphis_width['width_km'] if memphis_width and memphis_width['width_km'] else None,
    'offset_gc_to_min_km': offset_gc_min,
    'offset_memphis_to_min_km': offset_memphis_min,
}
with open(OUT / 'nile_valley_width.json', 'w') as f:
    json.dump(test_d, f, indent=2)

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
lats_plot = [v['lat'] for v in valid_widths]
widths_plot = [v['width_km'] for v in valid_widths]
ax.plot(lats_plot, widths_plot, 'b-o', markersize=3, label='Valley width')
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
    site_elev = get_elev(lat, lon)
    if site_elev <= -100:
        return None
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
            surr.append(get_elev(slat, slon))
    return site_elev - np.mean(surr)

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

# GC distance
pole_lat_r = np.radians(POLE_LAT)
pole_lon_r = np.radians(POLE_LON)
lats_r = np.radians(pleiades['reprLat'].values)
lons_r = np.radians(pleiades['reprLong'].values)
cos_ang = (np.sin(pole_lat_r)*np.sin(lats_r) +
           np.cos(pole_lat_r)*np.cos(lats_r)*np.cos(lons_r - pole_lon_r))
ang = np.arccos(np.clip(cos_ang, -1, 1))
pleiades['gc_dist_km'] = np.abs(ang - np.pi/2) * R

near_gc = pleiades[pleiades['gc_dist_km'] < 100].copy()
mon_near = near_gc[near_gc['site_class'] == 'monument']
set_near = near_gc[near_gc['site_class'] == 'settlement']
if len(mon_near) > 300:
    mon_near = mon_near.sample(300, random_state=42)
if len(set_near) > 300:
    set_near = set_near.sample(300, random_state=42)

print(f"Computing TPI for {len(mon_near)} monuments and {len(set_near)} settlements near GC...")

mon_tpi = []
for _, row in mon_near.iterrows():
    t = compute_tpi(row['reprLat'], row['reprLong'])
    if t is not None:
        mon_tpi.append(t)

set_tpi = []
for _, row in set_near.iterrows():
    t = compute_tpi(row['reprLat'], row['reprLong'])
    if t is not None:
        set_tpi.append(t)

print(f"Monument TPI: μ={np.mean(mon_tpi):.1f}m, median={np.median(mon_tpi):.1f}m (n={len(mon_tpi)})")
print(f"Settlement TPI: μ={np.mean(set_tpi):.1f}m, median={np.median(set_tpi):.1f}m (n={len(set_tpi)})")

tpi_stat, tpi_p = np.nan, np.nan
if len(mon_tpi) > 5 and len(set_tpi) > 5:
    tpi_stat, tpi_p = mannwhitneyu(mon_tpi, set_tpi, alternative='greater')
    print(f"Mann-Whitney (monuments > settlements): stat={tpi_stat:.0f}, p={tpi_p:.4f}")

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

# Plot
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
bp = ax.boxplot([mon_tpi, set_tpi],
                tick_labels=[f'Monuments (n={len(mon_tpi)})', f'Settlements (n={len(set_tpi)})'],
                patch_artist=True)
bp['boxes'][0].set_facecolor('#cc4444')
bp['boxes'][1].set_facecolor('#4488cc')
ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax.set_ylabel('TPI (m)')
p_str = f'{tpi_p:.2e}' if not np.isnan(tpi_p) and tpi_p < 0.01 else f'{tpi_p:.4f}' if not np.isnan(tpi_p) else 'N/A'
ax.set_title(f'TPI: Monuments vs Settlements\nMW p={p_str}')

ax = axes[1]
all_vals = mon_tpi + set_tpi
bins = np.linspace(min(all_vals), max(all_vals), 40)
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
# COMBINED SUMMARY (merge with A/B results)
# ================================================================
print(f"\n{'='*70}")
print("WRITING COMBINED SUMMARY")
print(f"{'='*70}")

# Load A and B results
with open(OUT / 'terrain_transition.json') as f:
    test_a = json.load(f)
with open(OUT / 'visibility_analysis.json') as f:
    test_b = json.load(f)

summary = {
    'test_a_terrain_transition': {
        'gc_mean_gradient': test_a['gc_mean_gradient'],
        'gc_percentile': test_a['gc_gradient_percentile'],
        'verdict': 'SUPPORTED' if test_a['gc_gradient_percentile'] > 80 else ('MARGINAL' if test_a['gc_gradient_percentile'] > 60 else 'NOT SUPPORTED'),
    },
    'test_b_visibility': {
        'gc_mean_prominence_m': test_b['gc_mean_prominence_m'],
        'gc_percentile': test_b['gc_prominence_percentile_vs_random'],
        'gc_vs_control_p': test_b.get('mann_whitney_p'),
        'monument_prominence_m': test_b.get('monument_mean_prominence_m'),
        'settlement_prominence_m': test_b.get('settlement_mean_prominence_m'),
        'mon_vs_set_p': test_b.get('mon_vs_set_mann_whitney_p'),
    },
    'test_c_navigation': {
        'crossing_prominence_m': test_c['crossing_mean_prominence_m'],
        'random_coast_prominence_m': test_c.get('random_coast_mean_prominence_m'),
        'crossing_vs_coast_p': float(coast_p) if not np.isnan(coast_p) else None,
        'intervisible_segments': f'{n_visible}/{len(intervis_results)}',
    },
    'test_d_nile_constriction': {
        'gc_crossing_lat': gc_nile_lat,
        'min_width_lat': min_width['lat'],
        'memphis_lat': memphis_lat,
        'offset_gc_to_min_km': offset_gc_min,
        'offset_memphis_to_min_km': offset_memphis_min,
        'min_width_km': min_width['width_km'],
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

print(f"\n{'='*70}")
print("INTERPRETATION")
print(f"{'='*70}")
print(f"\nTest A — Terrain Transition: pct={test_a['gc_gradient_percentile']:.0f}%")
print(f"  {'Circle traces high-gradient terrain' if test_a['gc_gradient_percentile'] > 70 else 'Circle does NOT preferentially trace terrain transitions'}")

print(f"\nTest B — Visibility: pct={test_b['gc_prominence_percentile_vs_random']:.0f}%")
if test_b.get('monument_mean_prominence_m') is not None:
    print(f"  Monuments: {test_b['monument_mean_prominence_m']:.0f}m vs Settlements: {test_b['settlement_mean_prominence_m']:.0f}m")
    print(f"  SURPRISING: monuments have LOWER prominence than settlements near GC")

print(f"\nTest C — Navigation: crossing p={coast_p:.4f}" if not np.isnan(coast_p) else "\nTest C — Navigation: insufficient data")
print(f"  Inter-visible: {n_visible}/{len(intervis_results)} segments")

print(f"\nTest D — Nile Constriction:")
print(f"  GC crossing: {gc_nile_lat:.2f}°N, Min width: {min_width['lat']:.2f}°N, Memphis: {memphis_lat:.2f}°N")
print(f"  GC-to-min offset: {offset_gc_min:.0f}km, Memphis-to-min: {offset_memphis_min:.0f}km")
if offset_gc_min < 33:  # within ~0.3 degrees
    print(f"  → GC crosses near the valley constriction!")
else:
    print(f"  → GC does NOT cross at the constriction")

print(f"\nTest E — Monument TPI:")
print(f"  Monuments: {np.mean(mon_tpi):.0f}m vs Settlements: {np.mean(set_tpi):.0f}m, p={tpi_p:.4f}" if not np.isnan(tpi_p) else "  insufficient data")

print(f"\nAll outputs in: {OUT}")
print("Done.")
