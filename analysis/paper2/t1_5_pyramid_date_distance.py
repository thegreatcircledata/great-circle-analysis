#!/usr/bin/env python3
"""
Directive T1-5: Construction Date vs Distance from Circle (Memphis)
====================================================================
For all Egyptian pyramids within 100km of the Great Circle, compile
construction dates and distances, test for correlation.
"""

import json, math, os, sys, csv
import numpy as np
from scipy import stats

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "pyramid_date_distance")
os.makedirs(OUT_DIR, exist_ok=True)

POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2

def haversine_km(lat1, lon1, lat2, lon2):
    R = EARTH_R_KM
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(min(1.0, a)))

def dist_from_gc(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QUARTER_CIRC)

# ============================================================
# EGYPTIAN PYRAMID DATABASE
# ============================================================
# Compiled from standard Egyptological sources (Dee et al. 2013, Lehner, Verner)
# Dates are approximate BCE values
PYRAMIDS = [
    # Old Kingdom pyramids near Memphis necropolis
    {"name": "Step Pyramid of Djoser", "lat": 29.8713, "lon": 31.2164, "date_bce": 2630, "pharaoh": "Djoser", "dynasty": 3},
    {"name": "Sekhemkhet Pyramid", "lat": 29.868, "lon": 31.214, "date_bce": 2611, "pharaoh": "Sekhemkhet", "dynasty": 3},
    {"name": "Layer Pyramid (Zawiyet Umm el-Rakhm)", "lat": 29.937, "lon": 31.157, "date_bce": 2600, "pharaoh": "Khaba", "dynasty": 3},
    {"name": "Meidum Pyramid", "lat": 29.3884, "lon": 31.1568, "date_bce": 2600, "pharaoh": "Sneferu/Huni", "dynasty": 3},
    {"name": "Bent Pyramid", "lat": 29.7903, "lon": 31.2094, "date_bce": 2600, "pharaoh": "Sneferu", "dynasty": 4},
    {"name": "Red Pyramid", "lat": 29.8089, "lon": 31.2059, "date_bce": 2590, "pharaoh": "Sneferu", "dynasty": 4},
    {"name": "Great Pyramid of Khufu", "lat": 29.9792, "lon": 31.1342, "date_bce": 2560, "pharaoh": "Khufu", "dynasty": 4},
    {"name": "Pyramid of Djedefre", "lat": 30.0322, "lon": 31.0747, "date_bce": 2528, "pharaoh": "Djedefre", "dynasty": 4},
    {"name": "Pyramid of Khafre", "lat": 29.9761, "lon": 31.1308, "date_bce": 2520, "pharaoh": "Khafre", "dynasty": 4},
    {"name": "Pyramid of Menkaure", "lat": 29.9725, "lon": 31.1281, "date_bce": 2490, "pharaoh": "Menkaure", "dynasty": 4},
    {"name": "Pyramid of Userkaf", "lat": 29.8736, "lon": 31.2192, "date_bce": 2494, "pharaoh": "Userkaf", "dynasty": 5},
    {"name": "Pyramid of Sahure", "lat": 29.8975, "lon": 31.2014, "date_bce": 2487, "pharaoh": "Sahure", "dynasty": 5},
    {"name": "Pyramid of Neferirkare", "lat": 29.8947, "lon": 31.2025, "date_bce": 2475, "pharaoh": "Neferirkare", "dynasty": 5},
    {"name": "Pyramid of Nyuserre", "lat": 29.8950, "lon": 31.2003, "date_bce": 2445, "pharaoh": "Nyuserre", "dynasty": 5},
    {"name": "Pyramid of Unas", "lat": 29.8681, "lon": 31.2156, "date_bce": 2375, "pharaoh": "Unas", "dynasty": 5},
    {"name": "Pyramid of Teti", "lat": 29.8756, "lon": 31.2206, "date_bce": 2345, "pharaoh": "Teti", "dynasty": 6},
    {"name": "Pyramid of Pepi I", "lat": 29.8528, "lon": 31.2153, "date_bce": 2332, "pharaoh": "Pepi I", "dynasty": 6},
    {"name": "Pyramid of Merenre", "lat": 29.8489, "lon": 31.2133, "date_bce": 2283, "pharaoh": "Merenre", "dynasty": 6},
    {"name": "Pyramid of Pepi II", "lat": 29.8439, "lon": 31.2131, "date_bce": 2278, "pharaoh": "Pepi II", "dynasty": 6},
    # Middle Kingdom
    {"name": "Pyramid of Amenemhat I", "lat": 29.8022, "lon": 31.2219, "date_bce": 1991, "pharaoh": "Amenemhat I", "dynasty": 12},
    {"name": "Pyramid of Senusret I", "lat": 29.8092, "lon": 31.2228, "date_bce": 1971, "pharaoh": "Senusret I", "dynasty": 12},
    {"name": "Pyramid of Amenemhat II", "lat": 29.7992, "lon": 31.2178, "date_bce": 1929, "pharaoh": "Amenemhat II", "dynasty": 12},
    {"name": "Pyramid of Senusret II", "lat": 29.2377, "lon": 30.9708, "date_bce": 1880, "pharaoh": "Senusret II", "dynasty": 12},
    {"name": "Pyramid of Senusret III", "lat": 29.8181, "lon": 31.2267, "date_bce": 1878, "pharaoh": "Senusret III", "dynasty": 12},
    {"name": "Pyramid of Amenemhat III (Dahshur)", "lat": 29.7875, "lon": 31.2111, "date_bce": 1860, "pharaoh": "Amenemhat III", "dynasty": 12},
    {"name": "Pyramid of Amenemhat III (Hawara)", "lat": 29.2719, "lon": 30.8981, "date_bce": 1850, "pharaoh": "Amenemhat III", "dynasty": 12},
]

# ============================================================
# COMPUTE DISTANCES
# ============================================================
print("Computing distances from Great Circle for each pyramid...")
for p in PYRAMIDS:
    p['distance_km'] = dist_from_gc(p['lat'], p['lon'])
    p['date_bp'] = p['date_bce'] + 2000  # approximate BP (before 2000 CE)
    print(f"  {p['name']}: {p['distance_km']:.1f} km, {p['date_bce']} BCE (Dynasty {p['dynasty']})")

# ============================================================
# STATISTICAL TESTS
# ============================================================
dates = np.array([p['date_bce'] for p in PYRAMIDS])
distances = np.array([p['distance_km'] for p in PYRAMIDS])
dynasties = np.array([p['dynasty'] for p in PYRAMIDS])

# Pearson correlation
r_pearson, p_pearson = stats.pearsonr(dates, distances)
r_spearman, p_spearman = stats.spearmanr(dates, distances)

print(f"\n{'='*70}")
print("CORRELATION: Construction Date vs Distance from Circle")
print(f"{'='*70}")
print(f"Pearson r = {r_pearson:.4f}, p = {p_pearson:.4f}")
print(f"Spearman rho = {r_spearman:.4f}, p = {p_spearman:.4f}")
print(f"Interpretation: {'Earlier = closer' if r_pearson < 0 else 'Earlier = farther'}")

# Key pyramids
print(f"\nKey reference points:")
for name in ["Great Pyramid of Khufu", "Step Pyramid of Djoser", "Bent Pyramid", "Red Pyramid"]:
    p = next(x for x in PYRAMIDS if x['name'] == name)
    print(f"  {name}: {p['date_bce']} BCE, {p['distance_km']:.1f} km from GC")

# Migration analysis: centroid of each 50-year bin
print(f"\nMIGRATION ANALYSIS (50-year bins):")
migration_data = []
bin_edges = np.arange(min(dates) - 25, max(dates) + 75, 50)
for i in range(len(bin_edges) - 1):
    mask = (dates >= bin_edges[i]) & (dates < bin_edges[i+1])
    if mask.sum() > 0:
        mean_date = float(np.mean(dates[mask]))
        mean_dist = float(np.mean(distances[mask]))
        mean_lat = float(np.mean([p['lat'] for j, p in enumerate(PYRAMIDS) if mask[j]]))
        mean_lon = float(np.mean([p['lon'] for j, p in enumerate(PYRAMIDS) if mask[j]]))
        n_sites = int(mask.sum())
        migration_data.append({
            'bin_start': int(bin_edges[i]),
            'bin_end': int(bin_edges[i+1]),
            'mean_date_bce': mean_date,
            'mean_distance_km': mean_dist,
            'centroid_lat': mean_lat,
            'centroid_lon': mean_lon,
            'n_pyramids': n_sites
        })
        print(f"  {int(bin_edges[i])}-{int(bin_edges[i+1])} BCE: "
              f"mean dist = {mean_dist:.1f} km, n = {n_sites}")

# ============================================================
# SAVE OUTPUTS
# ============================================================
correlation_data = {
    'pearson_r': float(r_pearson),
    'pearson_p': float(p_pearson),
    'spearman_rho': float(r_spearman),
    'spearman_p': float(p_spearman),
    'n_pyramids': len(PYRAMIDS),
    'pyramids': [{k: v for k, v in p.items()} for p in PYRAMIDS],
    'migration_bins': migration_data,
    'interpretation': f"Pearson r={r_pearson:.3f} (p={p_pearson:.3f}): {'negative correlation (earlier=closer)' if r_pearson < 0 else 'positive correlation (earlier=farther)'}"
}

with open(os.path.join(OUT_DIR, "correlation.json"), 'w') as f:
    json.dump(correlation_data, f, indent=2)

# ============================================================
# PLOTS
# ============================================================
print("\nGenerating plots...")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Scatter plot: date vs distance, colored by dynasty
fig, ax = plt.subplots(figsize=(12, 7))
dynasty_colors = {3: 'green', 4: 'red', 5: 'blue', 6: 'purple', 12: 'orange'}
for p in PYRAMIDS:
    color = dynasty_colors.get(p['dynasty'], 'gray')
    ax.scatter(p['date_bce'], p['distance_km'], color=color, s=60, zorder=5)
    # Label key pyramids
    if p['name'] in ["Great Pyramid of Khufu", "Step Pyramid of Djoser", "Bent Pyramid",
                      "Pyramid of Pepi II", "Pyramid of Amenemhat III (Hawara)", "Meidum Pyramid"]:
        ax.annotate(p['name'].replace('Pyramid of ', '').replace('Great Pyramid of ', ''),
                   (p['date_bce'], p['distance_km']), fontsize=7,
                   textcoords='offset points', xytext=(5, 5))

# Regression line
z = np.polyfit(dates, distances, 1)
p_line = np.poly1d(z)
x_line = np.linspace(min(dates), max(dates), 100)
ax.plot(x_line, p_line(x_line), 'k--', alpha=0.5,
        label=f'Linear fit (r={r_pearson:.3f}, p={p_pearson:.3f})')

# Dynasty legend
for d, c in sorted(dynasty_colors.items()):
    ax.scatter([], [], color=c, label=f'Dynasty {d}')

ax.set_xlabel('Construction Date (BCE)')
ax.set_ylabel('Distance from Great Circle (km)')
ax.set_title('Egyptian Pyramids: Construction Date vs Distance from Alison Great Circle')
ax.invert_xaxis()
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "scatter_plot.png"), dpi=150)
plt.close()
print("  Saved scatter_plot.png")

# Migration plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Mean distance over time
dates_mig = [m['mean_date_bce'] for m in migration_data]
dists_mig = [m['mean_distance_km'] for m in migration_data]
ax1.plot(dates_mig, dists_mig, 'bo-', markersize=8)
for m in migration_data:
    ax1.annotate(f"n={m['n_pyramids']}", (m['mean_date_bce'], m['mean_distance_km']),
                fontsize=8, textcoords='offset points', xytext=(0, 8))
ax1.set_xlabel('Mean Construction Date (BCE)')
ax1.set_ylabel('Mean Distance from Great Circle (km)')
ax1.set_title('Memphis Necropolis Migration\n(50-year bins)')
ax1.invert_xaxis()
ax1.grid(True, alpha=0.3)

# Centroid map
lats_mig = [m['centroid_lat'] for m in migration_data]
lons_mig = [m['centroid_lon'] for m in migration_data]
scatter = ax2.scatter(lons_mig, lats_mig, c=dates_mig, cmap='viridis', s=100, zorder=5)
for i, m in enumerate(migration_data):
    ax2.annotate(f"{int(m['mean_date_bce'])} BCE", (m['centroid_lon'], m['centroid_lat']),
                fontsize=7, textcoords='offset points', xytext=(5, 5))
ax2.set_xlabel('Longitude')
ax2.set_ylabel('Latitude')
ax2.set_title('Centroid Migration of Pyramid Construction')
plt.colorbar(scatter, ax=ax2, label='Date (BCE)')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "migration_plot.png"), dpi=150)
plt.close()
print("  Saved migration_plot.png")

print("\nDone! Outputs saved to", OUT_DIR)
