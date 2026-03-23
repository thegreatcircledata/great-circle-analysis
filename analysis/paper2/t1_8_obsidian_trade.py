#!/usr/bin/env python3
"""
Directive T1-8: Obsidian Trade Network Alignment
=================================================
Test whether Neolithic Near Eastern obsidian trade routes align with the Great Circle.
Uses published obsidian sourcing data (Cauvin 1998, Renfrew, Frahm 2012).
"""

import json, math, os, sys
import numpy as np

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "obsidian_trade")
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

def dist_from_gc(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    d = haversine_km(lat, lon, pole_lat, pole_lon)
    return abs(d - QUARTER_CIRC)

def bearing(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(dlon)
    return math.degrees(math.atan2(x, y)) % 360

def gc_bearing_at(lat, lon):
    """Approximate GC tangent bearing at a point near the circle."""
    pole_lat_r = math.radians(POLE_LAT)
    pole_lon_r = math.radians(POLE_LON)
    lat_r = math.radians(lat)
    lon_r = math.radians(lon)
    # Bearing from point to pole
    b = bearing(lat, lon, POLE_LAT, POLE_LON)
    # GC tangent is perpendicular to the bearing-to-pole direction
    return (b + 90) % 360

# ============================================================
# OBSIDIAN SOURCE AND DESTINATION DATABASE
# ============================================================
# Compiled from Frahm 2012, Cauvin 1998, Renfrew et al.

OBSIDIAN_SOURCES = {
    'Göllü Dağ': {'lat': 38.27, 'lon': 34.07, 'region': 'Central Anatolia'},
    'Nenezi Dağ': {'lat': 38.45, 'lon': 34.30, 'region': 'Central Anatolia'},
    'Acıgöl': {'lat': 38.53, 'lon': 34.52, 'region': 'Central Anatolia'},
    'Bingöl A': {'lat': 39.10, 'lon': 40.50, 'region': 'Eastern Anatolia'},
    'Bingöl B': {'lat': 39.15, 'lon': 40.55, 'region': 'Eastern Anatolia'},
    'Nemrut Dağ': {'lat': 38.65, 'lon': 42.23, 'region': 'Eastern Anatolia'},
    'Süphan Dağ': {'lat': 38.93, 'lon': 42.83, 'region': 'Eastern Anatolia'},
    'Arteni': {'lat': 40.52, 'lon': 44.17, 'region': 'Armenia'},
    'Gutansar': {'lat': 40.38, 'lon': 44.58, 'region': 'Armenia'},
    'Melos (Aegean)': {'lat': 36.68, 'lon': 24.42, 'region': 'Aegean'},
}

# Destination sites with known obsidian imports and approximate primary sources
OBSIDIAN_DESTINATIONS = [
    # Site, lat, lon, primary source(s), approximate date BCE, notes
    {'name': 'Çatalhöyük', 'lat': 37.67, 'lon': 32.83, 'sources': ['Göllü Dağ', 'Nenezi Dağ'], 'date_bce': 7000},
    {'name': 'Aşıklı Höyük', 'lat': 38.34, 'lon': 34.22, 'sources': ['Göllü Dağ'], 'date_bce': 8200},
    {'name': 'Abu Hureyra', 'lat': 35.87, 'lon': 38.40, 'sources': ['Bingöl A', 'Bingöl B', 'Nemrut Dağ'], 'date_bce': 9500},
    {'name': 'Jericho (Tell es-Sultan)', 'lat': 31.87, 'lon': 35.44, 'sources': ['Göllü Dağ', 'Nemrut Dağ'], 'date_bce': 8500},
    {'name': 'Byblos', 'lat': 34.12, 'lon': 35.65, 'sources': ['Göllü Dağ', 'Nemrut Dağ'], 'date_bce': 7000},
    {'name': 'Ain Ghazal', 'lat': 31.98, 'lon': 35.93, 'sources': ['Göllü Dağ', 'Nemrut Dağ'], 'date_bce': 7500},
    {'name': 'Jarmo', 'lat': 35.55, 'lon': 44.95, 'sources': ['Nemrut Dağ'], 'date_bce': 7000},
    {'name': 'Çayönü', 'lat': 38.22, 'lon': 39.73, 'sources': ['Bingöl A', 'Bingöl B'], 'date_bce': 8600},
    {'name': 'Göbekli Tepe', 'lat': 37.22, 'lon': 38.92, 'sources': ['Bingöl A', 'Nemrut Dağ'], 'date_bce': 9500},
    {'name': 'Nevali Çori', 'lat': 37.52, 'lon': 38.61, 'sources': ['Bingöl A'], 'date_bce': 8500},
    {'name': 'Tell Halula', 'lat': 36.43, 'lon': 38.17, 'sources': ['Bingöl A', 'Göllü Dağ'], 'date_bce': 7500},
    {'name': 'Mersin-Yumuktepe', 'lat': 36.78, 'lon': 34.60, 'sources': ['Göllü Dağ', 'Acıgöl'], 'date_bce': 6000},
    {'name': 'Tell Aswad', 'lat': 33.60, 'lon': 36.60, 'sources': ['Göllü Dağ', 'Nemrut Dağ'], 'date_bce': 8500},
    {'name': 'Bouqras', 'lat': 35.27, 'lon': 40.37, 'sources': ['Nemrut Dağ', 'Bingöl A'], 'date_bce': 6400},
    {'name': 'Tell Ramad', 'lat': 33.50, 'lon': 36.00, 'sources': ['Göllü Dağ'], 'date_bce': 6000},
    {'name': 'Franchthi Cave', 'lat': 37.42, 'lon': 23.13, 'sources': ['Melos (Aegean)'], 'date_bce': 9000},
    {'name': 'Ali Kosh', 'lat': 32.40, 'lon': 47.55, 'sources': ['Nemrut Dağ'], 'date_bce': 7500},
    {'name': 'Hajji Firuz', 'lat': 37.00, 'lon': 45.50, 'sources': ['Nemrut Dağ', 'Süphan Dağ'], 'date_bce': 5500},
    {'name': 'Domuztepe', 'lat': 37.25, 'lon': 36.87, 'sources': ['Göllü Dağ', 'Bingöl A'], 'date_bce': 5500},
    {'name': 'Knossos', 'lat': 35.30, 'lon': 25.16, 'sources': ['Melos (Aegean)'], 'date_bce': 7000},
]

# ============================================================
# ANALYSIS
# ============================================================
print("Computing trade route alignments with Great Circle...\n")

trade_routes = []
for dest in OBSIDIAN_DESTINATIONS:
    dest_dist = dist_from_gc(dest['lat'], dest['lon'])
    for src_name in dest['sources']:
        if src_name not in OBSIDIAN_SOURCES:
            continue
        src = OBSIDIAN_SOURCES[src_name]
        src_dist = dist_from_gc(src['lat'], src['lon'])

        # Midpoint of trade route
        mid_lat = (dest['lat'] + src['lat']) / 2
        mid_lon = (dest['lon'] + src['lon']) / 2
        mid_dist = dist_from_gc(mid_lat, mid_lon)

        # Trade route bearing
        trade_bearing = bearing(src['lat'], src['lon'], dest['lat'], dest['lon'])

        # GC bearing at midpoint
        gc_bear = gc_bearing_at(mid_lat, mid_lon)

        # Angle between trade route and GC
        angle = abs(trade_bearing - gc_bear) % 180
        if angle > 90:
            angle = 180 - angle

        # Distance of trade route
        trade_dist_km = haversine_km(src['lat'], src['lon'], dest['lat'], dest['lon'])

        # Does the route cross the GC?
        # Approximate: check if source and dest are on opposite sides
        src_side = haversine_km(src['lat'], src['lon'], POLE_LAT, POLE_LON) - QUARTER_CIRC
        dest_side = haversine_km(dest['lat'], dest['lon'], POLE_LAT, POLE_LON) - QUARTER_CIRC
        crosses_gc = (src_side * dest_side) < 0

        route = {
            'source': src_name,
            'destination': dest['name'],
            'source_lat': src['lat'],
            'source_lon': src['lon'],
            'dest_lat': dest['lat'],
            'dest_lon': dest['lon'],
            'source_dist_gc': float(src_dist),
            'dest_dist_gc': float(dest_dist),
            'midpoint_dist_gc': float(mid_dist),
            'trade_bearing': float(trade_bearing),
            'gc_bearing_at_mid': float(gc_bear),
            'angle_with_gc': float(angle),
            'trade_distance_km': float(trade_dist_km),
            'crosses_gc': bool(crosses_gc),
            'date_bce': dest['date_bce']
        }
        trade_routes.append(route)

        print(f"  {src_name} → {dest['name']}: {trade_dist_km:.0f}km, "
              f"angle={angle:.0f}°, mid_dist_gc={mid_dist:.0f}km, crosses={crosses_gc}")

# ============================================================
# RANDOM CIRCLE COMPARISON
# ============================================================
print(f"\nComparing crossing rates to random circles...")
N_RANDOM = 200

def random_pole():
    z = np.random.uniform(0, 1)
    theta = np.random.uniform(0, 2 * math.pi)
    return math.degrees(math.asin(z)), math.degrees(theta) - 180

alison_crossing_rate = sum(1 for r in trade_routes if r['crosses_gc']) / len(trade_routes)
random_crossing_rates = []

for i in range(N_RANDOM):
    plat, plon = random_pole()
    qc = EARTH_R_KM * math.pi / 2
    crosses = 0
    for r in trade_routes:
        src_d = haversine_km(r['source_lat'], r['source_lon'], plat, plon) - qc
        dst_d = haversine_km(r['dest_lat'], r['dest_lon'], plat, plon) - qc
        if src_d * dst_d < 0:
            crosses += 1
    random_crossing_rates.append(crosses / len(trade_routes))

# ============================================================
# RESULTS
# ============================================================
print(f"\n{'='*70}")
print("OBSIDIAN TRADE NETWORK ANALYSIS RESULTS")
print(f"{'='*70}")

n_crossing = sum(1 for r in trade_routes if r['crosses_gc'])
n_total = len(trade_routes)
mean_angle = float(np.mean([r['angle_with_gc'] for r in trade_routes]))
mean_mid_dist = float(np.mean([r['midpoint_dist_gc'] for r in trade_routes]))

print(f"\n1. GC CROSSING RATE")
print(f"   {n_crossing}/{n_total} trade routes cross the Alison GC ({alison_crossing_rate*100:.1f}%)")
print(f"   Random mean: {np.mean(random_crossing_rates)*100:.1f}% ± {np.std(random_crossing_rates)*100:.1f}%")
pct = sum(1 for r in random_crossing_rates if r >= alison_crossing_rate) / len(random_crossing_rates) * 100
print(f"   Percentile: {100-pct:.1f}th")

print(f"\n2. TRADE ROUTE ANGLE WITH GC")
print(f"   Mean angle: {mean_angle:.1f}°")
print(f"   Parallel (<30°): {sum(1 for r in trade_routes if r['angle_with_gc'] < 30)}/{n_total}")
print(f"   Perpendicular (>60°): {sum(1 for r in trade_routes if r['angle_with_gc'] > 60)}/{n_total}")

print(f"\n3. MIDPOINT PROXIMITY TO GC")
print(f"   Mean midpoint distance: {mean_mid_dist:.0f} km")

# Source vs destination proximity
src_dists = [r['source_dist_gc'] for r in trade_routes]
dst_dists = [r['dest_dist_gc'] for r in trade_routes]
print(f"\n4. SOURCE vs DESTINATION PROXIMITY")
print(f"   Mean source distance from GC: {np.mean(src_dists):.0f} km")
print(f"   Mean destination distance from GC: {np.mean(dst_dists):.0f} km")
near_gc_dests = [r for r in trade_routes if r['dest_dist_gc'] < 200]
far_gc_dests = [r for r in trade_routes if r['dest_dist_gc'] >= 200]
if near_gc_dests:
    near_sources = set(r['source'] for r in near_gc_dests)
    far_sources = set(r['source'] for r in far_gc_dests) if far_gc_dests else set()
    print(f"   Near-GC sites import from: {near_sources}")
    print(f"   Far-GC sites import from: {far_sources}")
    overlap = near_sources & far_sources
    unique_near = near_sources - far_sources
    print(f"   Shared sources: {overlap}")
    print(f"   Unique to near-GC: {unique_near}")

# ============================================================
# SAVE
# ============================================================
results = {
    'n_sources': len(OBSIDIAN_SOURCES),
    'n_destinations': len(OBSIDIAN_DESTINATIONS),
    'n_trade_routes': n_total,
    'crossing_rate': alison_crossing_rate,
    'random_crossing_mean': float(np.mean(random_crossing_rates)),
    'random_crossing_std': float(np.std(random_crossing_rates)),
    'mean_angle_with_gc': mean_angle,
    'mean_midpoint_dist': mean_mid_dist,
    'sources': {k: v for k, v in OBSIDIAN_SOURCES.items()},
}

with open(os.path.join(OUT_DIR, "trade_routes_vs_gc.json"), 'w') as f:
    json.dump({'trade_routes': trade_routes}, f, indent=2)

with open(os.path.join(OUT_DIR, "crossing_analysis.json"), 'w') as f:
    json.dump(results, f, indent=2)

# Plot
print("\nGenerating plot...")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(14, 8))

# Plot GC (just the Near East segment)
def great_circle_points_local(pole_lat, pole_lon, n_points=3600):
    pole_lat_r = math.radians(pole_lat)
    pole_lon_r = math.radians(pole_lon)
    bearings = np.linspace(0, 2 * math.pi, n_points, endpoint=False)
    lats = np.empty(n_points)
    lons = np.empty(n_points)
    for i, b in enumerate(bearings):
        lat = math.asin(math.sin(pole_lat_r)*math.cos(math.pi/2)+math.cos(pole_lat_r)*math.sin(math.pi/2)*math.cos(b))
        lon = pole_lon_r + math.atan2(math.sin(b)*math.sin(math.pi/2)*math.cos(pole_lat_r), math.cos(math.pi/2)-math.sin(pole_lat_r)*math.sin(lat))
        lats[i] = math.degrees(lat)
        lons[i] = math.degrees(lon)
    lons = ((lons + 180) % 360) - 180
    return lats, lons

gc_lats, gc_lons = great_circle_points_local(POLE_LAT, POLE_LON)
mask = (gc_lons > 15) & (gc_lons < 55) & (gc_lats > 25) & (gc_lats < 45)
ax.plot(gc_lons[mask], gc_lats[mask], 'b-', linewidth=2, alpha=0.5, label='Great Circle')

# Plot sources
for name, src in OBSIDIAN_SOURCES.items():
    if 15 < src['lon'] < 55 and 25 < src['lat'] < 45:
        ax.plot(src['lon'], src['lat'], '^', color='red', markersize=12, zorder=10)
        ax.annotate(name, (src['lon'], src['lat']), fontsize=7, textcoords='offset points', xytext=(5, 5))

# Plot trade routes
for r in trade_routes:
    if r['dest_lon'] < 15 or r['dest_lon'] > 55:
        continue
    color = 'green' if r['crosses_gc'] else 'orange'
    ax.plot([r['source_lon'], r['dest_lon']], [r['source_lat'], r['dest_lat']],
            '-', color=color, alpha=0.5, linewidth=1)
    ax.plot(r['dest_lon'], r['dest_lat'], 'o', color=color, markersize=5)

ax.set_xlim(15, 55)
ax.set_ylim(28, 42)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Obsidian Trade Routes vs Alison Great Circle')
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='blue', alpha=0.5, linewidth=2, label='Great Circle'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='red', markersize=10, label='Obsidian source'),
    Line2D([0], [0], color='green', alpha=0.5, label='Route crosses GC'),
    Line2D([0], [0], color='orange', alpha=0.5, label='Route does not cross GC'),
]
ax.legend(handles=legend_elements)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "obsidian_trade_map.png"), dpi=150)
plt.close()
print("  Saved obsidian_trade_map.png")

print(f"\nDone! Outputs saved to {OUT_DIR}")
