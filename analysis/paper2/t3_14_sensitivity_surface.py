#!/usr/bin/env python3
"""
Directive T3-14: Sensitivity Surface Map
=========================================
Shift the Alison pole by 1° increments in all directions (±30° lat/lon).
At each position, compute monument-settlement divergence D.
Produce a 2D heatmap of D values.
"""

import json, math, os, sys, csv
import numpy as np

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "sensitivity_surface")
os.makedirs(OUT_DIR, exist_ok=True)

POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
THRESHOLD_KM = 50
N_TRIALS = 100  # Reduced for speed (61x61 grid = 3,721 evaluations)
GRID_STEP = 1   # degrees
GRID_RANGE = 30 # ±30 degrees

MONUMENTAL_TYPES = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus", "shrine"
}
SETTLEMENT_TYPES = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement", "townhouse", "production"
}

# ============================================================
# VECTORIZED FUNCTIONS
# ============================================================
def haversine_km_vec(lat1, lon1, lat2, lon2):
    """Vectorized haversine distance in km."""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    return EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))

def dist_from_gc_vec(lats, lons, pole_lat, pole_lon):
    d = haversine_km_vec(lats, lons, pole_lat, pole_lon)
    return np.abs(d - QUARTER_CIRC)

def count_within(dists, threshold):
    return int(np.sum(dists < threshold))

# ============================================================
# LOAD DATA
# ============================================================
print("Loading Pleiades data...")
path = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")
mon_lats, mon_lons = [], []
set_lats, set_lons = [], []

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
            mon_lats.append(lat)
            mon_lons.append(lon)
        elif ftypes & SETTLEMENT_TYPES:
            set_lats.append(lat)
            set_lons.append(lon)

mon_lats = np.array(mon_lats)
mon_lons = np.array(mon_lons)
set_lats = np.array(set_lats)
set_lons = np.array(set_lons)
all_lats = np.concatenate([mon_lats, set_lats])
all_lons = np.concatenate([mon_lons, set_lons])

print(f"Loaded {len(mon_lats)} monuments, {len(set_lats)} settlements")

# ============================================================
# FAST D COMPUTATION
# ============================================================
def compute_D_fast(pole_lat, pole_lon, n_trials=N_TRIALS):
    """Fast divergence computation with vectorized operations."""
    mon_dists = dist_from_gc_vec(mon_lats, mon_lons, pole_lat, pole_lon)
    set_dists = dist_from_gc_vec(set_lats, set_lons, pole_lat, pole_lon)

    obs_mon = count_within(mon_dists, THRESHOLD_KM)
    obs_set = count_within(set_dists, THRESHOLD_KM)

    mc_mon = np.empty(n_trials)
    mc_set = np.empty(n_trials)

    for t in range(n_trials):
        jitter_lat = np.random.normal(0, 2, len(all_lats))
        jitter_lon = np.random.normal(0, 2, len(all_lons))
        j_lats = all_lats + jitter_lat
        j_lons = all_lons + jitter_lon

        idx_m = np.random.choice(len(j_lats), len(mon_lats), replace=True)
        idx_s = np.random.choice(len(j_lats), len(set_lats), replace=True)

        jm_dists = dist_from_gc_vec(j_lats[idx_m], j_lons[idx_m], pole_lat, pole_lon)
        js_dists = dist_from_gc_vec(j_lats[idx_s], j_lons[idx_s], pole_lat, pole_lon)

        mc_mon[t] = count_within(jm_dists, THRESHOLD_KM)
        mc_set[t] = count_within(js_dists, THRESHOLD_KM)

    z_mon = (obs_mon - mc_mon.mean()) / max(mc_mon.std(), 1)
    z_set = (obs_set - mc_set.mean()) / max(mc_set.std(), 1)

    return float(z_mon - z_set)

# ============================================================
# GRID COMPUTATION
# ============================================================
lat_offsets = np.arange(-GRID_RANGE, GRID_RANGE + 1, GRID_STEP)
lon_offsets = np.arange(-GRID_RANGE, GRID_RANGE + 1, GRID_STEP)
n_lat = len(lat_offsets)
n_lon = len(lon_offsets)
total = n_lat * n_lon

print(f"\nComputing D on {n_lat}x{n_lon} = {total} grid points...")
print(f"Pole center: ({POLE_LAT:.2f}, {POLE_LON:.2f})")
print(f"Grid: lat [{POLE_LAT-GRID_RANGE:.0f}, {POLE_LAT+GRID_RANGE:.0f}], "
      f"lon [{POLE_LON-GRID_RANGE:.0f}, {POLE_LON+GRID_RANGE:.0f}]")

D_grid = np.zeros((n_lat, n_lon))
count = 0

for i, dlat in enumerate(lat_offsets):
    plat = POLE_LAT + dlat
    if plat > 90:
        plat = 180 - plat  # wrap over pole
    if plat < -90:
        plat = -180 - plat

    for j, dlon in enumerate(lon_offsets):
        plon = POLE_LON + dlon
        plon = ((plon + 180) % 360) - 180

        D = compute_D_fast(plat, plon)
        D_grid[i, j] = D
        count += 1

        if count % 100 == 0:
            print(f"  {count}/{total} ({100*count/total:.0f}%) — "
                  f"current: dlat={dlat:+.0f}, dlon={dlon:+.0f}, D={D:.2f}")

# ============================================================
# RESULTS
# ============================================================
center_i = n_lat // 2
center_j = n_lon // 2
alison_D = D_grid[center_i, center_j]
max_D = float(np.max(D_grid))
max_pos = np.unravel_index(np.argmax(D_grid), D_grid.shape)
max_dlat = lat_offsets[max_pos[0]]
max_dlon = lon_offsets[max_pos[1]]

print(f"\n{'='*70}")
print("SENSITIVITY SURFACE RESULTS")
print(f"{'='*70}")
print(f"Alison pole D: {alison_D:.2f}")
print(f"Maximum D: {max_D:.2f} at offset ({max_dlat:+.0f}°, {max_dlon:+.0f}°)")
print(f"Maximum D pole: ({POLE_LAT+max_dlat:.2f}, {POLE_LON+max_dlon:.2f})")

# Gradient analysis
grad_lat = np.gradient(D_grid, axis=0)
grad_lon = np.gradient(D_grid, axis=1)
grad_mag = np.sqrt(grad_lat**2 + grad_lon**2)

# Sensitivity at Alison pole
grad_at_alison = grad_mag[center_i, center_j]
grad_lat_at = grad_lat[center_i, center_j]
grad_lon_at = grad_lon[center_i, center_j]
steepest_dir = math.degrees(math.atan2(grad_lon_at, grad_lat_at))

print(f"\nGradient at Alison pole:")
print(f"  Magnitude: {grad_at_alison:.4f} D/degree")
print(f"  Lat sensitivity: {abs(grad_lat_at):.4f}")
print(f"  Lon sensitivity: {abs(grad_lon_at):.4f}")
print(f"  Steepest descent direction: {steepest_dir:.0f}°")

# Check for multiple peaks
# Find local maxima
from scipy import ndimage
local_max = ndimage.maximum_filter(D_grid, size=5) == D_grid
local_max_D = D_grid[local_max]
peaks = [(float(lat_offsets[i]), float(lon_offsets[j]), float(D_grid[i, j]))
         for i in range(n_lat) for j in range(n_lon)
         if local_max[i, j] and D_grid[i, j] > alison_D * 0.5]
peaks.sort(key=lambda x: -x[2])

print(f"\nSignificant peaks (D > {alison_D*0.5:.1f}):")
for dlat, dlon, d in peaks[:10]:
    print(f"  Offset ({dlat:+.0f}°, {dlon:+.0f}°): D = {d:.2f}")

# ============================================================
# SAVE
# ============================================================
sensitivity_data = {
    'pole_center': {'lat': POLE_LAT, 'lon': POLE_LON},
    'grid_range_deg': GRID_RANGE,
    'grid_step_deg': GRID_STEP,
    'lat_offsets': lat_offsets.tolist(),
    'lon_offsets': lon_offsets.tolist(),
    'D_grid': D_grid.tolist(),
    'alison_D': alison_D,
    'max_D': max_D,
    'max_offset': {'dlat': float(max_dlat), 'dlon': float(max_dlon)},
    'n_trials_per_point': N_TRIALS,
    'threshold_km': THRESHOLD_KM
}
with open(os.path.join(OUT_DIR, "sensitivity_data.json"), 'w') as f:
    json.dump(sensitivity_data, f, indent=2)

gradient_data = {
    'alison_gradient_magnitude': float(grad_at_alison),
    'lat_sensitivity': float(abs(grad_lat_at)),
    'lon_sensitivity': float(abs(grad_lon_at)),
    'steepest_descent_direction_deg': float(steepest_dir),
    'peaks': [{'dlat': p[0], 'dlon': p[1], 'D': p[2]} for p in peaks[:20]]
}
with open(os.path.join(OUT_DIR, "gradient_analysis.json"), 'w') as f:
    json.dump(gradient_data, f, indent=2)

# ============================================================
# HEATMAP PLOT
# ============================================================
print("\nGenerating heatmap...")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(14, 12))
im = ax.imshow(D_grid, origin='lower', aspect='auto',
               extent=[lon_offsets[0], lon_offsets[-1], lat_offsets[0], lat_offsets[-1]],
               cmap='RdYlBu_r', interpolation='bilinear')

# Mark Alison pole (center)
ax.plot(0, 0, 'k*', markersize=15, label=f'Alison pole (D={alison_D:.2f})')

# Mark maximum
ax.plot(max_dlon, max_dlat, 'w^', markersize=12, label=f'Maximum D={max_D:.2f}')

# Contours
contour_levels = np.linspace(np.min(D_grid), np.max(D_grid), 10)
cs = ax.contour(lon_offsets, lat_offsets, D_grid, levels=contour_levels,
                colors='black', linewidths=0.5, alpha=0.5)
ax.clabel(cs, fontsize=7, fmt='%.1f')

plt.colorbar(im, ax=ax, label='Monument-Settlement Divergence D', shrink=0.8)
ax.set_xlabel('Longitude Offset from Alison Pole (degrees)')
ax.set_ylabel('Latitude Offset from Alison Pole (degrees)')
ax.set_title('Sensitivity Surface: Divergence D vs Pole Position\n'
             f'Center: Alison pole ({POLE_LAT:.1f}°N, {POLE_LON:.1f}°W)')
ax.legend(loc='upper right')
ax.grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "heatmap.png"), dpi=150)
plt.close()
print("  Saved heatmap.png")

print(f"\nDone! Outputs saved to {OUT_DIR}")
