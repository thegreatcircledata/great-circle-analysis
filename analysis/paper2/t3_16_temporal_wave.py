#!/usr/bin/env python3
"""
Directive T3-16: Temporal Wave Analysis
========================================
Test whether monument construction along the circle shows a directional wave.
Did construction start at one end and propagate along the line?
"""

import json, math, os, sys
import numpy as np
from scipy import stats

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "temporal_wave")
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

def arc_position(lat, lon):
    """Compute arc position (0-360°) along the great circle."""
    lat_r = math.radians(lat)
    lon_r = math.radians(lon)
    pole_lat_r = math.radians(POLE_LAT)
    pole_lon_r = math.radians(POLE_LON)

    # Convert to Cartesian
    x = math.cos(lat_r) * math.cos(lon_r)
    y = math.cos(lat_r) * math.sin(lon_r)
    z = math.sin(lat_r)

    # Pole in Cartesian
    px = math.cos(pole_lat_r) * math.cos(pole_lon_r)
    py = math.cos(pole_lat_r) * math.sin(pole_lon_r)
    pz = math.sin(pole_lat_r)

    # Project point onto plane perpendicular to pole
    dot = x * px + y * py + z * pz
    proj_x = x - dot * px
    proj_y = y - dot * py
    proj_z = z - dot * pz

    # Normalize
    norm = math.sqrt(proj_x**2 + proj_y**2 + proj_z**2)
    if norm < 1e-10:
        return 0.0
    proj_x /= norm
    proj_y /= norm
    proj_z /= norm

    # Need a reference direction in the GC plane
    # Use the point where GC crosses the prime meridian-ish
    ref_lat = 0  # equatorial
    ref_lon = POLE_LON + 90  # 90° from pole longitude
    ref_lat_r = math.radians(ref_lat)
    ref_lon_r = math.radians(ref_lon)
    rx = math.cos(ref_lat_r) * math.cos(ref_lon_r)
    ry = math.cos(ref_lat_r) * math.sin(ref_lon_r)
    rz = math.sin(ref_lat_r)

    # Project reference onto GC plane
    rdot = rx * px + ry * py + rz * pz
    rproj_x = rx - rdot * px
    rproj_y = ry - rdot * py
    rproj_z = rz - rdot * pz
    rnorm = math.sqrt(rproj_x**2 + rproj_y**2 + rproj_z**2)
    if rnorm < 1e-10:
        return 0.0
    rproj_x /= rnorm
    rproj_y /= rnorm
    rproj_z /= rnorm

    # Second basis vector (cross product of pole and ref)
    bx = py * rproj_z - pz * rproj_y
    by = pz * rproj_x - px * rproj_z
    bz = px * rproj_y - py * rproj_x

    # Angle
    cos_a = proj_x * rproj_x + proj_y * rproj_y + proj_z * rproj_z
    sin_a = proj_x * bx + proj_y * by + proj_z * bz

    angle = math.degrees(math.atan2(sin_a, cos_a)) % 360
    return angle

# ============================================================
# MONUMENT DATABASE WITH CONSTRUCTION DATES
# ============================================================
# Compiled from published sources: Dee et al. 2013, standard references
# Dates are approximate BCE values

MONUMENTS = [
    # Egypt
    {"name": "Step Pyramid of Djoser", "lat": 29.8713, "lon": 31.2164, "date_bce": 2630},
    {"name": "Great Pyramid of Khufu", "lat": 29.9792, "lon": 31.1342, "date_bce": 2560},
    {"name": "Pyramid of Khafre", "lat": 29.9761, "lon": 31.1308, "date_bce": 2520},
    {"name": "Temple of Karnak (earliest)", "lat": 25.7188, "lon": 32.6573, "date_bce": 2055},
    {"name": "Temple of Luxor", "lat": 25.6995, "lon": 32.6392, "date_bce": 1400},
    {"name": "Abu Simbel", "lat": 22.3369, "lon": 31.6256, "date_bce": 1264},

    # Levant / Anatolia
    {"name": "Göbekli Tepe", "lat": 37.2233, "lon": 38.9225, "date_bce": 9500},
    {"name": "Jericho Tower", "lat": 31.871, "lon": 35.444, "date_bce": 8000},
    {"name": "Ain Ghazal statues", "lat": 31.98, "lon": 35.93, "date_bce": 7500},
    {"name": "Çatalhöyük shrine", "lat": 37.67, "lon": 32.83, "date_bce": 7000},
    {"name": "Byblos temple", "lat": 34.12, "lon": 35.65, "date_bce": 3000},

    # Mesopotamia
    {"name": "Eridu temple", "lat": 30.82, "lon": 45.99, "date_bce": 5400},
    {"name": "Uruk (White Temple)", "lat": 31.32, "lon": 45.64, "date_bce": 3500},
    {"name": "Ur Ziggurat", "lat": 30.96, "lon": 46.10, "date_bce": 2100},

    # Iran
    {"name": "Tall-e Bakun", "lat": 30.04, "lon": 52.39, "date_bce": 4000},
    {"name": "Chogha Zanbil Ziggurat", "lat": 32.01, "lon": 48.52, "date_bce": 1250},
    {"name": "Persepolis", "lat": 29.93, "lon": 52.89, "date_bce": 515},
    {"name": "Naqsh-e Rostam", "lat": 29.99, "lon": 52.87, "date_bce": 500},

    # Indus Valley
    {"name": "Mehrgarh", "lat": 29.37, "lon": 67.62, "date_bce": 7000},
    {"name": "Mohenjo-daro", "lat": 27.33, "lon": 68.14, "date_bce": 2500},
    {"name": "Harappa", "lat": 30.63, "lon": 72.87, "date_bce": 2600},

    # South/Southeast Asia
    {"name": "Sanchi Stupa", "lat": 23.48, "lon": 77.74, "date_bce": 250},
    {"name": "Borobudur", "lat": -7.61, "lon": 110.20, "date_bce": -800},  # 800 CE
    {"name": "Angkor Wat", "lat": 13.41, "lon": 103.87, "date_bce": -1113},  # 1113 CE
    {"name": "Prambanan", "lat": -7.75, "lon": 110.49, "date_bce": -850},  # 850 CE

    # Pacific
    {"name": "Easter Island Ahu (earliest)", "lat": -27.12, "lon": -109.35, "date_bce": -800},  # 800 CE
    {"name": "Nan Madol", "lat": 6.84, "lon": 158.33, "date_bce": -1100},  # 1100 CE

    # Americas
    {"name": "Caral/Norte Chico", "lat": -10.89, "lon": -77.52, "date_bce": 3000},
    {"name": "Nazca Lines", "lat": -14.74, "lon": -75.13, "date_bce": 200},
    {"name": "Machu Picchu", "lat": -13.16, "lon": -72.55, "date_bce": -1450},  # 1450 CE
    {"name": "Serpent Mound", "lat": 39.025, "lon": -83.431, "date_bce": 300},
    {"name": "Poverty Point", "lat": 32.63, "lon": -91.41, "date_bce": 1700},

    # Europe (near circle)
    {"name": "Newgrange", "lat": 53.69, "lon": -6.48, "date_bce": 3200},
    {"name": "Stonehenge", "lat": 51.18, "lon": -1.83, "date_bce": 3000},
    {"name": "Carnac stones", "lat": 47.59, "lon": -3.07, "date_bce": 4500},
]

# ============================================================
# COMPUTE ARC POSITIONS AND DISTANCES
# ============================================================
print("Computing arc positions and distances from Great Circle...")
for m in MONUMENTS:
    m['dist_from_gc'] = dist_from_gc(m['lat'], m['lon'])
    m['arc_position'] = arc_position(m['lat'], m['lon'])
    # Convert negative BCE to actual year for consistent plotting
    m['year_bce'] = m['date_bce']  # negative means CE
    print(f"  {m['name']}: arc={m['arc_position']:.1f}°, dist={m['dist_from_gc']:.0f}km, "
          f"date={'%d BCE' % m['date_bce'] if m['date_bce'] > 0 else '%d CE' % -m['date_bce']}")

# Filter to sites within 500km of GC (generous threshold for global analysis)
near_gc = [m for m in MONUMENTS if m['dist_from_gc'] < 500]
print(f"\nSites within 500km of GC: {len(near_gc)}/{len(MONUMENTS)}")

# ============================================================
# STATISTICAL TESTS
# ============================================================
arc_pos = np.array([m['arc_position'] for m in near_gc])
dates = np.array([m['date_bce'] for m in near_gc])

# Test 1: Pearson/Spearman correlation
r_pearson, p_pearson = stats.pearsonr(arc_pos, dates)
r_spearman, p_spearman = stats.spearmanr(arc_pos, dates)

print(f"\n{'='*70}")
print("TEMPORAL WAVE ANALYSIS RESULTS")
print(f"{'='*70}")
print(f"\n1. ARC POSITION vs CONSTRUCTION DATE CORRELATION")
print(f"   Pearson r = {r_pearson:.4f}, p = {p_pearson:.4f}")
print(f"   Spearman rho = {r_spearman:.4f}, p = {p_spearman:.4f}")
if p_pearson < 0.05:
    direction = "construction progresses in one direction" if r_pearson > 0 else "construction progresses in reverse direction"
    print(f"   Significant! {direction}")
else:
    print(f"   Not significant — no clear directional wave")

# Test 2: Wavefront analysis (500-year bins)
print(f"\n2. WAVEFRONT ANALYSIS (500-year bins)")
bins = np.arange(-2000, 10000, 500)  # CE to 10000 BCE
wavefront = []
for i in range(len(bins) - 1):
    mask = (dates >= bins[i]) & (dates < bins[i+1])
    if mask.sum() > 0:
        mean_arc = float(np.mean(arc_pos[mask]))
        mean_date = float(np.mean(dates[mask]))
        n_sites = int(mask.sum())
        names_in_bin = [m['name'] for j, m in enumerate(near_gc) if mask[j]]
        wavefront.append({
            'bin_start_bce': int(bins[i]),
            'bin_end_bce': int(bins[i+1]),
            'mean_arc_position': mean_arc,
            'mean_date_bce': mean_date,
            'n_sites': n_sites,
            'sites': names_in_bin
        })
        bin_label = f"{bins[i]}-{bins[i+1]} BCE" if bins[i] > 0 else f"{-bins[i+1]}-{-bins[i]} CE"
        print(f"   {bin_label}: mean arc = {mean_arc:.0f}°, n = {n_sites}")

# Test 3: Speed estimate
if p_pearson < 0.05 and len(near_gc) > 5:
    slope = np.polyfit(dates, arc_pos, 1)[0]  # degrees per year
    arc_km = slope * (QUARTER_CIRC * 4 / 360)  # convert degrees to km
    speed_km_per_century = abs(arc_km) * 100
    print(f"\n3. PROPAGATION SPEED ESTIMATE")
    print(f"   {speed_km_per_century:.0f} km/century")
else:
    speed_km_per_century = None
    print(f"\n3. PROPAGATION SPEED: not computed (no significant wave)")

# Test 4: Compare to off-circle monuments
off_gc = [m for m in MONUMENTS if m['dist_from_gc'] >= 500]
if len(off_gc) > 5:
    off_arc = np.array([m['arc_position'] for m in off_gc])
    off_dates = np.array([m['date_bce'] for m in off_gc])
    r_off, p_off = stats.pearsonr(off_arc, off_dates) if len(off_gc) > 2 else (0, 1)
    print(f"\n4. OFF-CIRCLE COMPARISON (>500km from GC)")
    print(f"   n = {len(off_gc)}")
    print(f"   Pearson r = {r_off:.4f}, p = {p_off:.4f}")
else:
    r_off, p_off = None, None
    print(f"\n4. OFF-CIRCLE: insufficient data ({len(off_gc)} sites)")

# ============================================================
# SAVE
# ============================================================
# Monument dates CSV
import csv
with open(os.path.join(OUT_DIR, "monument_dates.csv"), 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['name', 'lat', 'lon', 'date_bce', 'arc_position', 'dist_from_gc'])
    writer.writeheader()
    for m in MONUMENTS:
        writer.writerow({
            'name': m['name'], 'lat': m['lat'], 'lon': m['lon'],
            'date_bce': m['date_bce'], 'arc_position': m['arc_position'],
            'dist_from_gc': m['dist_from_gc']
        })

wave_results = {
    'pearson_r': float(r_pearson),
    'pearson_p': float(p_pearson),
    'spearman_rho': float(r_spearman),
    'spearman_p': float(p_spearman),
    'n_on_circle': len(near_gc),
    'n_off_circle': len(off_gc),
    'propagation_speed_km_per_century': speed_km_per_century,
    'wavefront_bins': wavefront,
    'off_circle_pearson_r': float(r_off) if r_off is not None else None,
    'off_circle_pearson_p': float(p_off) if p_off is not None else None,
    'interpretation': (
        f"Significant directional wave detected (r={r_pearson:.3f}, p={p_pearson:.4f}). "
        f"Speed ≈ {speed_km_per_century:.0f} km/century."
        if p_pearson < 0.05 else
        f"No significant directional wave (r={r_pearson:.3f}, p={p_pearson:.4f}). "
        "Construction appears to occur independently at multiple locations."
    )
}
with open(os.path.join(OUT_DIR, "wave_analysis.json"), 'w') as f:
    json.dump(wave_results, f, indent=2)

# ============================================================
# PLOTS
# ============================================================
print("\nGenerating plots...")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# Plot 1: Arc position vs construction date
colors_region = []
for m in near_gc:
    if m['lon'] < 35:
        colors_region.append('blue')  # Egypt/Europe
    elif m['lon'] < 55:
        colors_region.append('green')  # Levant/Mesopotamia
    elif m['lon'] < 90:
        colors_region.append('orange')  # Iran/India
    elif m['lon'] < 130 or m['lon'] < -60:
        colors_region.append('red')  # SE Asia/Americas
    else:
        colors_region.append('purple')  # Pacific

ax1.scatter(arc_pos, dates, c=colors_region, s=80, zorder=5)
for i, m in enumerate(near_gc):
    ax1.annotate(m['name'].split('(')[0].strip()[:20], (arc_pos[i], dates[i]),
                fontsize=6, textcoords='offset points', xytext=(3, 3), rotation=20)

# Regression line if significant
if p_pearson < 0.1:
    z = np.polyfit(arc_pos, dates, 1)
    p_line = np.poly1d(z)
    x_line = np.linspace(0, 360, 100)
    ax1.plot(x_line, p_line(x_line), 'k--', alpha=0.5,
             label=f'Linear fit (r={r_pearson:.3f}, p={p_pearson:.3f})')
    ax1.legend()

ax1.set_xlabel('Arc Position Along Great Circle (degrees)')
ax1.set_ylabel('Construction Date (BCE, negative=CE)')
ax1.set_title('Monument Construction Date vs Arc Position')
ax1.grid(True, alpha=0.3)

# Plot 2: Wavefront over time
if wavefront:
    wf_dates = [w['mean_date_bce'] for w in wavefront]
    wf_arcs = [w['mean_arc_position'] for w in wavefront]
    wf_sizes = [w['n_sites'] * 20 for w in wavefront]
    ax2.scatter(wf_dates, wf_arcs, s=wf_sizes, c='steelblue', alpha=0.7, zorder=5)
    ax2.plot(wf_dates, wf_arcs, 'b-', alpha=0.3)
    for w in wavefront:
        ax2.annotate(f"n={w['n_sites']}", (w['mean_date_bce'], w['mean_arc_position']),
                    fontsize=7, textcoords='offset points', xytext=(5, 5))

ax2.set_xlabel('Mean Construction Date (BCE)')
ax2.set_ylabel('Mean Arc Position (degrees)')
ax2.set_title('Wavefront: Mean Arc Position Over Time\n(500-year bins, size ∝ site count)')
ax2.invert_xaxis()
ax2.grid(True, alpha=0.3)

plt.suptitle('Temporal Wave Analysis — Alison Great Circle', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "arc_position_vs_date.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved arc_position_vs_date.png")

print(f"\nDone! Outputs saved to {OUT_DIR}")
