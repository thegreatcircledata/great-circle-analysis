#!/usr/bin/env python3
"""
Memphis Necropolis Deep Dive — Ahramat Branch × Great Circle Geometry
======================================================================
WHY do 38 Old Kingdom pyramids cluster within 10 km of Alison's Great Circle?

Analysis:
1. Compute the Great Circle's bearing at Memphis (~30°N, 31°E)
2. Establish the Ahramat Branch azimuth (Lisht → Giza, ~SSE-NNW)
3. Compute the angle between them
4. Map every known Old Kingdom pyramid, compute distance to GC
5. Counterfactual: shift the Ahramat Branch east/west and recount
6. Test: does pyramid building stop along the circle after 4.2 ka event?

Key references:
- Ghoneim et al. 2024, Comms Earth & Env: Ahramat Branch 64 km, Lisht to Giza
- Sheisha et al. 2022, PNAS: Khufu Branch, water levels during construction
- Lehner GPMP: causeway orientations
"""

import math, json, os, sys
import numpy as np

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE = os.path.expanduser("~/megalith_site_research")
OUT = os.path.join(BASE, "outputs/memphis_deep_dive")
os.makedirs(OUT, exist_ok=True)

# ============================================================
# GREAT CIRCLE PARAMETERS
# ============================================================
POLE_LAT = 59.682122
POLE_LON = -138.646087
R = 6371.0
QC = R * np.pi / 2

# ============================================================
# KNOWN PYRAMID COORDINATES
# ============================================================
# Major Old Kingdom and Middle Kingdom pyramids along the Memphis necropolis
# Sources: Lehner "The Complete Pyramids" (1997), Wikipedia (verified against
# satellite imagery), Mark & Jarus compilations
# Format: (name, lat, lon, dynasty, approx_date_BCE)

PYRAMIDS = [
    # ABU RAWASH
    ("Djedefre (Abu Rawash)", 30.0322, 31.0750, "4th", -2560),

    # GIZA
    ("Great Pyramid (Khufu)", 29.9792, 31.1342, "4th", -2560),
    ("Khafre", 29.9761, 31.1306, "4th", -2530),
    ("Menkaure", 29.9725, 31.1278, "4th", -2510),
    ("Queen Hetepheres I (G1a)", 29.9806, 31.1319, "4th", -2560),
    ("Queen Meritites I (G1b)", 29.9800, 31.1331, "4th", -2560),
    ("Queen Henutsen (G1c)", 29.9794, 31.1342, "4th", -2560),
    ("Queen Khamerernebty II (G3a)", 29.9711, 31.1272, "4th", -2510),

    # ZAWYET EL-ARYAN
    ("Unfinished Northern (Zawyet)", 29.9353, 31.1611, "3rd/4th", -2570),
    ("Layer Pyramid (Zawyet)", 29.9336, 31.1594, "3rd", -2630),

    # ABUSIR
    ("Sahure (Abusir)", 29.8975, 31.2028, "5th", -2480),
    ("Neferirkare (Abusir)", 29.8950, 31.2025, "5th", -2470),
    ("Niuserre (Abusir)", 29.8964, 31.2008, "5th", -2440),
    ("Neferefre (Abusir)", 29.8936, 31.2036, "5th", -2460),
    ("Raneferef (Abusir)", 29.8936, 31.2036, "5th", -2450),

    # SAQQARA
    ("Step Pyramid (Djoser)", 29.8714, 31.2164, "3rd", -2650),
    ("Sekhemkhet (Buried Pyramid)", 29.8683, 31.2142, "3rd", -2640),
    ("Userkaf", 29.8731, 31.2183, "5th", -2490),
    ("Teti", 29.8756, 31.2217, "6th", -2340),
    ("Pepi I", 29.8567, 31.2172, "6th", -2320),
    ("Merenre I", 29.8544, 31.2164, "6th", -2280),
    ("Pepi II", 29.8453, 31.2147, "6th", -2250),
    ("Ibi", 29.8464, 31.2153, "8th", -2170),
    ("Unas", 29.8681, 31.2181, "5th", -2360),
    ("Djedkare Isesi", 29.8508, 31.2156, "5th", -2380),

    # DAHSHUR
    ("Red Pyramid (Sneferu)", 29.8092, 31.2061, "4th", -2590),
    ("Bent Pyramid (Sneferu)", 29.7903, 31.2094, "4th", -2600),
    ("Black Pyramid (Amenemhat III)", 29.8014, 31.2264, "12th", -1850),
    ("White Pyramid (Amenemhat II)", 29.8053, 31.2208, "12th", -1900),
    ("Pyramid of Senusret III", 29.8019, 31.2222, "12th", -1870),
    ("Pyramid of Khendjer", 29.8019, 31.2347, "13th", -1760),

    # MAZGHUNA
    ("South Mazghuna", 29.7722, 31.2264, "12th/13th", -1800),
    ("North Mazghuna", 29.7753, 31.2256, "12th/13th", -1800),

    # LISHT
    ("Amenemhat I (Lisht)", 29.5756, 31.2250, "12th", -1970),
    ("Senusret I (Lisht)", 29.5647, 31.2228, "12th", -1940),

    # MEIDUM (further south, sometimes included)
    ("Meidum Pyramid (Sneferu/Huni)", 29.3883, 31.1569, "3rd/4th", -2610),
]

# ============================================================
# AHRAMAT BRANCH APPROXIMATE CENTERLINE
# ============================================================
# Ghoneim et al. 2024: runs from Lisht (south) to Giza (north),
# 2.5-10.25 km west of modern Nile, ~64 km long
# The branch follows the foot of the Western Desert Plateau escarpment
# Approximate waypoints from the paper's Figure 2:
AHRAMAT_WAYPOINTS = [
    (29.57, 31.20),   # Southern end: Lisht area
    (29.65, 31.19),   # Between Lisht and Dahshur
    (29.72, 31.19),   # Mazghuna area
    (29.79, 31.18),   # Dahshur area
    (29.85, 31.18),   # South Saqqara
    (29.87, 31.18),   # Saqqara
    (29.90, 31.17),   # Abusir
    (29.93, 31.15),   # Zawyet el-Aryan
    (29.97, 31.12),   # Giza area (branch widens into inlet)
    (30.03, 31.07),   # Northern end: Abu Rawash area
]

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def haversine(lat1, lon1, lat2, lon2):
    la1, lo1, la2, lo2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dl, dn = la2 - la1, lo2 - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(la2)*np.sin(dn/2)**2
    return R * 2 * np.arcsin(np.sqrt(min(1.0, a)))

def gc_dist_point(lat, lon):
    """Distance from a single point to Alison's Great Circle."""
    d = haversine(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QC)

def gc_bearing_at(lat, lon):
    """
    Compute the bearing (azimuth) of the Great Circle at a given point.
    The GC is the equator of the pole (POLE_LAT, POLE_LON).
    The bearing at any point is perpendicular to the direction toward the pole.
    """
    lat_r = np.radians(lat)
    lon_r = np.radians(lon)
    pole_lat_r = np.radians(POLE_LAT)
    pole_lon_r = np.radians(POLE_LON)

    # Bearing from point to pole
    dlon = pole_lon_r - lon_r
    x = np.sin(dlon) * np.cos(pole_lat_r)
    y = np.cos(lat_r) * np.sin(pole_lat_r) - np.sin(lat_r) * np.cos(pole_lat_r) * np.cos(dlon)
    bearing_to_pole = np.degrees(np.arctan2(x, y)) % 360

    # GC bearing is perpendicular to bearing-to-pole
    # The circle goes in two directions; take the one going roughly east
    gc_bearing_1 = (bearing_to_pole + 90) % 360
    gc_bearing_2 = (bearing_to_pole - 90) % 360

    return gc_bearing_1, gc_bearing_2, bearing_to_pole

def line_bearing(lat1, lon1, lat2, lon2):
    """Initial bearing from point 1 to point 2."""
    la1, lo1, la2, lo2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlon = lo2 - lo1
    x = np.sin(dlon) * np.cos(la2)
    y = np.cos(la1) * np.sin(la2) - np.sin(la1) * np.cos(la2) * np.cos(dlon)
    return np.degrees(np.arctan2(x, y)) % 360

def angle_between(b1, b2):
    """Smallest angle between two bearings."""
    diff = abs(b1 - b2) % 360
    return min(diff, 360 - diff)

# ============================================================
# ANALYSIS 1: GREAT CIRCLE BEARING AT MEMPHIS
# ============================================================
print("=" * 70)
print("ANALYSIS 1: GREAT CIRCLE BEARING AT MEMPHIS")
print("=" * 70)

# Memphis is approximately at Saqqara (29.87°N, 31.22°E)
memphis_lat, memphis_lon = 29.87, 31.22

gc_b1, gc_b2, bearing_to_pole = gc_bearing_at(memphis_lat, memphis_lon)
print(f"\nMemphis location: {memphis_lat}°N, {memphis_lon}°E")
print(f"Distance to Great Circle: {gc_dist_point(memphis_lat, memphis_lon):.1f} km")
print(f"Bearing to pole: {bearing_to_pole:.1f}°")
print(f"Great Circle bearing at Memphis: {gc_b1:.1f}° / {gc_b2:.1f}°")

# Also compute at Giza and Lisht
for name, lat, lon in [("Giza", 29.98, 31.13), ("Lisht", 29.57, 31.22),
                         ("Dahshur", 29.80, 31.21), ("Abu Rawash", 30.03, 31.07)]:
    b1, b2, btp = gc_bearing_at(lat, lon)
    d = gc_dist_point(lat, lon)
    print(f"  {name:12s}: {lat:.2f}°N, {lon:.2f}°E  →  GC bearing {b1:.1f}°/{b2:.1f}°, dist {d:.1f} km")

# ============================================================
# ANALYSIS 2: AHRAMAT BRANCH AZIMUTH
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 2: AHRAMAT BRANCH AZIMUTH")
print("=" * 70)

# Overall bearing: Lisht to Giza
ahramat_bearing = line_bearing(
    AHRAMAT_WAYPOINTS[0][0], AHRAMAT_WAYPOINTS[0][1],
    AHRAMAT_WAYPOINTS[-1][0], AHRAMAT_WAYPOINTS[-1][1]
)
print(f"\nAhramat Branch overall bearing (Lisht → Abu Rawash): {ahramat_bearing:.1f}°")
print(f"Reciprocal (Abu Rawash → Lisht): {(ahramat_bearing + 180) % 360:.1f}°")

# Segment-by-segment bearings
print("\nSegment bearings:")
for i in range(len(AHRAMAT_WAYPOINTS) - 1):
    lat1, lon1 = AHRAMAT_WAYPOINTS[i]
    lat2, lon2 = AHRAMAT_WAYPOINTS[i+1]
    seg_bearing = line_bearing(lat1, lon1, lat2, lon2)
    seg_dist = haversine(lat1, lon1, lat2, lon2)
    print(f"  Seg {i}: ({lat1:.2f},{lon1:.2f}) → ({lat2:.2f},{lon2:.2f}): "
          f"bearing {seg_bearing:.1f}°, dist {seg_dist:.1f} km")

# ============================================================
# ANALYSIS 3: ANGLE BETWEEN GC AND AHRAMAT BRANCH
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 3: ANGLE BETWEEN GREAT CIRCLE AND AHRAMAT BRANCH")
print("=" * 70)

gc_bearing_at_memphis = gc_b1  # the eastward-going bearing
angle = angle_between(gc_bearing_at_memphis, ahramat_bearing)
angle2 = angle_between(gc_b2, ahramat_bearing)
best_angle = min(angle, angle2)

print(f"\nGreat Circle bearing at Memphis: {gc_b1:.1f}° and {gc_b2:.1f}°")
print(f"Ahramat Branch bearing: {ahramat_bearing:.1f}°")
print(f"Angle between GC and Ahramat: {best_angle:.1f}°")

if best_angle < 15:
    print("→ NEARLY PARALLEL: The Great Circle runs almost parallel to the Ahramat Branch")
elif best_angle < 30:
    print("→ OBLIQUE: The Great Circle crosses the Ahramat Branch at a shallow angle")
elif best_angle < 60:
    print("→ MODERATE ANGLE: The Great Circle crosses at a significant angle")
elif best_angle < 80:
    print("→ NEAR-PERPENDICULAR: The Great Circle is nearly perpendicular to the branch")
else:
    print("→ PERPENDICULAR: The Great Circle crosses the Ahramat Branch at ~90°")

# Causeway direction: causeways run perpendicular to the branch (east-west)
# toward the river. If GC is parallel to branch, causeways are perpendicular to GC.
causeway_bearing = (ahramat_bearing + 90) % 360
causeway_angle_to_gc = angle_between(causeway_bearing, gc_bearing_at_memphis)
print(f"\nCauseway bearing (perpendicular to branch): ~{causeway_bearing:.1f}°")
print(f"Angle between causeways and GC: {causeway_angle_to_gc:.1f}°")

# ============================================================
# ANALYSIS 4: EVERY PYRAMID'S DISTANCE TO THE GREAT CIRCLE
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 4: PYRAMID DISTANCES TO GREAT CIRCLE")
print("=" * 70)

pyramid_results = []
print(f"\n{'Pyramid':<45} {'Dist km':>8} {'Dynasty':>8} {'Date BCE':>10}")
print("-" * 75)

for name, lat, lon, dynasty, date in sorted(PYRAMIDS, key=lambda x: gc_dist_point(x[1], x[2])):
    dist = gc_dist_point(lat, lon)
    pyramid_results.append({
        'name': name, 'lat': lat, 'lon': lon, 'dynasty': dynasty,
        'date_bce': date, 'dist_to_gc_km': round(dist, 2)
    })
    marker = " ***" if dist < 10 else ""
    print(f"{name:<45} {dist:8.2f} {dynasty:>8} {abs(date):>10}{marker}")

within_10 = sum(1 for p in pyramid_results if p['dist_to_gc_km'] < 10)
within_25 = sum(1 for p in pyramid_results if p['dist_to_gc_km'] < 25)
within_50 = sum(1 for p in pyramid_results if p['dist_to_gc_km'] < 50)
print(f"\nPyramids within 10 km: {within_10}/{len(PYRAMIDS)}")
print(f"Pyramids within 25 km: {within_25}/{len(PYRAMIDS)}")
print(f"Pyramids within 50 km: {within_50}/{len(PYRAMIDS)}")

# ============================================================
# ANALYSIS 5: COUNTERFACTUAL — SHIFT THE NILE
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 5: COUNTERFACTUAL — WHAT IF THE NILE WERE ELSEWHERE?")
print("=" * 70)

# If the Ahramat Branch ran 10 km further east or west,
# the pyramids would have been built at a different longitude.
# How would this affect their distance to the Great Circle?

print("\nShifting all pyramid longitudes by ΔE and recomputing distances:")
print(f"{'Shift (km E)':>12} {'Within 10km':>12} {'Within 25km':>12} {'Within 50km':>12}")
print("-" * 52)

counterfactual_results = []
for shift_km in np.arange(-30, 31, 5):
    # 1 km ≈ 0.0114° longitude at 30°N
    lon_shift = shift_km * 0.0114
    count_10 = 0
    count_25 = 0
    count_50 = 0
    for name, lat, lon, dynasty, date in PYRAMIDS:
        d = gc_dist_point(lat, lon + lon_shift)
        if d < 10: count_10 += 1
        if d < 25: count_25 += 1
        if d < 50: count_50 += 1

    counterfactual_results.append({
        'shift_km': float(shift_km),
        'within_10': count_10,
        'within_25': count_25,
        'within_50': count_50
    })
    marker = " ← ACTUAL" if shift_km == 0 else ""
    print(f"{shift_km:>12.0f} {count_10:>12} {count_25:>12} {count_50:>12}{marker}")

# ============================================================
# ANALYSIS 6: TEMPORAL — CONSTRUCTION BEFORE & AFTER 4.2 ka
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 6: PYRAMID CONSTRUCTION vs 4.2 ka EVENT")
print("=" * 70)

pre_42ka = [p for p in PYRAMIDS if p[4] < -2200]  # before 2200 BCE
post_42ka = [p for p in PYRAMIDS if p[4] >= -2200]

pre_dists = [gc_dist_point(p[1], p[2]) for p in pre_42ka]
post_dists = [gc_dist_point(p[1], p[2]) for p in post_42ka]

print(f"\nPre-4.2 ka (before 2200 BCE): {len(pre_42ka)} pyramids")
print(f"  Mean distance to GC: {np.mean(pre_dists):.1f} km")
print(f"  Within 10 km: {sum(1 for d in pre_dists if d < 10)}")
print(f"  Within 25 km: {sum(1 for d in pre_dists if d < 25)}")

print(f"\nPost-4.2 ka (2200 BCE and later): {len(post_42ka)} pyramids")
print(f"  Mean distance to GC: {np.mean(post_dists):.1f} km")
print(f"  Within 10 km: {sum(1 for d in post_dists if d < 10)}")
print(f"  Within 25 km: {sum(1 for d in post_dists if d < 25)}")

print("\nPost-4.2 ka pyramids:")
for name, lat, lon, dynasty, date in sorted(post_42ka, key=lambda x: x[4]):
    d = gc_dist_point(lat, lon)
    print(f"  {abs(date)} BCE  {dynasty:>5}  {d:6.1f} km  {name}")

# ============================================================
# ANALYSIS 7: WHICH WAY DOES THE CIRCLE CROSS EGYPT?
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 7: GREAT CIRCLE PATH THROUGH EGYPT")
print("=" * 70)

# Sample the GC at points across Egypt
print("\nGreat Circle trajectory through Egypt (sampled every 0.5° longitude):")
print(f"{'Lon':>8} {'Lat on GC':>10} {'Location':>30}")
print("-" * 50)

for lon_test in np.arange(25, 37, 0.5):
    # Find lat where GC passes at this longitude
    # The GC is the set of points exactly QC from the pole
    # Solve: haversine(lat, lon_test, pole_lat, pole_lon) = QC
    # Binary search
    best_lat = None
    best_err = 999
    for lat_try in np.arange(20, 40, 0.01):
        d = haversine(lat_try, lon_test, POLE_LAT, POLE_LON)
        err = abs(d - QC)
        if err < best_err:
            best_err = err
            best_lat = lat_try

    if best_lat and best_err < 1:
        label = ""
        if 31.0 <= lon_test <= 31.5:
            label = "← Memphis necropolis"
        elif 32.5 <= lon_test <= 33.0:
            label = "← Nile Delta"
        elif 29.5 <= lon_test <= 30.0:
            label = "← Western Desert"
        print(f"{lon_test:>8.1f} {best_lat:>10.2f} {label:>30}")

# ============================================================
# SAVE & PLOT
# ============================================================
output = {
    'gc_bearing_memphis': {'b1': float(gc_b1), 'b2': float(gc_b2)},
    'ahramat_bearing': float(ahramat_bearing),
    'angle_gc_ahramat': float(best_angle),
    'pyramids': pyramid_results,
    'within_10km': within_10,
    'within_25km': within_25,
    'within_50km': within_50,
    'counterfactual': counterfactual_results,
    'pre_42ka': {
        'count': len(pre_42ka),
        'mean_dist': float(np.mean(pre_dists)),
        'within_10': sum(1 for d in pre_dists if d < 10),
    },
    'post_42ka': {
        'count': len(post_42ka),
        'mean_dist': float(np.mean(post_dists)),
        'within_10': sum(1 for d in post_dists if d < 10),
    },
}

with open(os.path.join(OUT, 'results.json'), 'w') as f:
    json.dump(output, f, indent=2)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # Top-left: Map of pyramids and GC path
    ax = axes[0, 0]
    pyr_lats = [p[1] for p in PYRAMIDS]
    pyr_lons = [p[2] for p in PYRAMIDS]
    pyr_dists = [gc_dist_point(p[1], p[2]) for p in PYRAMIDS]

    # Color by distance to GC
    sc = ax.scatter(pyr_lons, pyr_lats, c=pyr_dists, cmap='RdYlGn_r',
                    s=50, edgecolors='black', linewidths=0.5, vmin=0, vmax=50, zorder=5)
    plt.colorbar(sc, ax=ax, label='Distance to GC (km)')

    # Ahramat Branch
    ah_lats = [w[0] for w in AHRAMAT_WAYPOINTS]
    ah_lons = [w[1] for w in AHRAMAT_WAYPOINTS]
    ax.plot(ah_lons, ah_lats, 'b-', linewidth=3, alpha=0.5, label='Ahramat Branch (approx)')

    # GC path through Egypt
    gc_path_lons = []
    gc_path_lats = []
    for lon_test in np.arange(30.5, 32.0, 0.02):
        best_lat = None
        best_err = 999
        for lat_try in np.arange(28, 32, 0.005):
            d = haversine(lat_try, lon_test, POLE_LAT, POLE_LON)
            err = abs(d - QC)
            if err < best_err:
                best_err = err
                best_lat = lat_try
        if best_lat and best_err < 0.5:
            gc_path_lons.append(lon_test)
            gc_path_lats.append(best_lat)

    ax.plot(gc_path_lons, gc_path_lats, 'r-', linewidth=2, label='Great Circle', zorder=4)

    ax.set_xlabel('Longitude (°E)')
    ax.set_ylabel('Latitude (°N)')
    ax.set_title('Memphis Necropolis: Pyramids, Ahramat Branch & Great Circle')
    ax.legend(fontsize=8)
    ax.set_xlim(30.9, 31.5)
    ax.set_ylim(29.3, 30.15)
    ax.set_aspect('equal')
    ax.grid(alpha=0.3)

    # Top-right: Pyramid distance histogram
    ax2 = axes[0, 1]
    ax2.hist(pyr_dists, bins=20, color='#c0392b', alpha=0.7, edgecolor='black')
    ax2.axvline(x=10, color='green', linestyle='--', label='10 km')
    ax2.axvline(x=25, color='orange', linestyle='--', label='25 km')
    ax2.axvline(x=50, color='red', linestyle='--', label='50 km')
    ax2.set_xlabel('Distance to Great Circle (km)')
    ax2.set_ylabel('Number of Pyramids')
    ax2.set_title(f'Pyramid Distance Distribution (n={len(PYRAMIDS)})')
    ax2.legend()
    ax2.grid(alpha=0.3)

    # Bottom-left: Counterfactual
    ax3 = axes[1, 0]
    shifts = [c['shift_km'] for c in counterfactual_results]
    w10 = [c['within_10'] for c in counterfactual_results]
    w25 = [c['within_25'] for c in counterfactual_results]
    w50 = [c['within_50'] for c in counterfactual_results]
    ax3.plot(shifts, w10, 'g-o', linewidth=2, markersize=4, label='Within 10 km')
    ax3.plot(shifts, w25, 'orange', linewidth=2, marker='s', markersize=4, label='Within 25 km')
    ax3.plot(shifts, w50, 'r-^', linewidth=2, markersize=4, label='Within 50 km')
    ax3.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5, label='Actual position')
    ax3.set_xlabel('Nile Shift (km east)')
    ax3.set_ylabel('Pyramids within threshold')
    ax3.set_title('Counterfactual: If the Nile Were Shifted East/West')
    ax3.legend(fontsize=8)
    ax3.grid(alpha=0.3)

    # Bottom-right: Timeline
    ax4 = axes[1, 1]
    dates = [abs(p[4]) for p in PYRAMIDS]
    dists_timeline = [gc_dist_point(p[1], p[2]) for p in PYRAMIDS]
    dynasties = [p[3] for p in PYRAMIDS]

    colors_dyn = {'3rd': '#e74c3c', '3rd/4th': '#e67e22', '4th': '#f39c12',
                  '5th': '#2ecc71', '6th': '#3498db', '8th': '#9b59b6',
                  '12th': '#1abc9c', '12th/13th': '#16a085', '13th': '#2c3e50'}

    for d_name in sorted(set(dynasties)):
        mask = [i for i, dy in enumerate(dynasties) if dy == d_name]
        ax4.scatter([dates[i] for i in mask], [dists_timeline[i] for i in mask],
                   c=colors_dyn.get(d_name, 'gray'), s=50, edgecolors='black',
                   linewidths=0.5, label=f'Dynasty {d_name}', zorder=5)

    ax4.axhline(y=10, color='green', linestyle='--', alpha=0.5)
    ax4.axvline(x=2200, color='purple', linestyle='--', linewidth=2, alpha=0.7, label='4.2 ka event')
    ax4.set_xlabel('Date (BCE)')
    ax4.set_ylabel('Distance to Great Circle (km)')
    ax4.set_title('Pyramid Distance to GC Over Time')
    ax4.legend(fontsize=7, ncol=2)
    ax4.grid(alpha=0.3)
    ax4.invert_xaxis()

    plt.tight_layout()
    plt.savefig(os.path.join(OUT, 'memphis_deep_dive.png'), dpi=150)
    print(f"\nPlot saved: {os.path.join(OUT, 'memphis_deep_dive.png')}")
except ImportError:
    print("matplotlib not available")

print(f"Results saved: {os.path.join(OUT, 'results.json')}")
print("Done.")
