#!/usr/bin/env python3
"""
Analysis 3: Easter Island Ahu Distribution vs Great Circle Bearing
====================================================================
Directive 09, Script 04

Tests whether ahu (platform) positions around Easter Island's coastline
cluster near the points where the Great Circle intersects the island.

Methods:
  A. Compute GC bearing at Easter Island center
  B. Find GC-coastline intersection azimuths
  C. Compute bearing from island center to each ahu
  D. Test: does ahu density peak near GC intersection bearings?
  E. Rayleigh test on angular offsets from GC crossing points

Output:
  - ahu_orientation_test.json
  - ahu_density_by_bearing.json
  - easter_island_circle_overlay.png
"""

import json
import math
import os
import sys

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ─── Import data ─────────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from easter_island_ahu_data import AHU_DATA

# ─── Constants ────────────────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
R_EARTH = 6371.0
QUARTER_CIRC = R_EARTH * math.pi / 2

# Easter Island center (approximate)
EI_CENTER_LAT = -27.1127
EI_CENTER_LON = -109.3497

N_MC = 10000


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlat = la2 - la1
    dlon = lo2 - lo1
    a = math.sin(dlat/2)**2 + math.cos(la1)*math.cos(la2)*math.sin(dlon/2)**2
    return R_EARTH * 2 * math.asin(math.sqrt(min(1.0, a)))


def initial_bearing(lat1, lon1, lat2, lon2):
    """Forward azimuth from point 1 to point 2, in degrees [0, 360)."""
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlon = lo2 - lo1
    x = math.sin(dlon) * math.cos(la2)
    y = math.cos(la1)*math.sin(la2) - math.sin(la1)*math.cos(la2)*math.cos(dlon)
    return math.degrees(math.atan2(x, y)) % 360


def angular_diff(a, b):
    """Minimum angular difference between two bearings (0-180)."""
    d = abs(a - b) % 360
    return min(d, 360 - d)


# ═══════════════════════════════════════════════════════════════════════
# A. Great Circle bearing at Easter Island
# ═══════════════════════════════════════════════════════════════════════

print("=" * 70)
print("EASTER ISLAND AHU DISTRIBUTION vs GREAT CIRCLE")
print("=" * 70)

# Bearing computation (from Nazca script pattern)
plat = math.radians(POLE_LAT)
plon = math.radians(POLE_LON)
qlat = math.radians(EI_CENTER_LAT)
qlon = math.radians(EI_CENTER_LON)
dlon = plon - qlon
x = math.sin(dlon) * math.cos(plat)
y = (math.cos(qlat) * math.sin(plat) -
     math.sin(qlat) * math.cos(plat) * math.cos(dlon))
bearing_to_pole = math.degrees(math.atan2(x, y)) % 360

gc_bearing_1 = (bearing_to_pole + 90) % 360
gc_bearing_2 = (bearing_to_pole - 90) % 360

# Distance to circle
dist_to_pole = haversine_km(EI_CENTER_LAT, EI_CENTER_LON, POLE_LAT, POLE_LON)
dist_to_circle = abs(dist_to_pole - QUARTER_CIRC)

print(f"\nEaster Island center: ({EI_CENTER_LAT}, {EI_CENTER_LON})")
print(f"Distance to Great Circle: {dist_to_circle:.1f} km")
print(f"GC bearing at Easter Island: {gc_bearing_1:.2f}° / {gc_bearing_2:.2f}°")

# The two GC crossing azimuths from island center
# (the circle enters from one side and exits the other)
gc_crossing_1 = gc_bearing_1
gc_crossing_2 = gc_bearing_2
print(f"\nGC crosses the island at bearings: {gc_crossing_1:.1f}° and {gc_crossing_2:.1f}°")


# ═══════════════════════════════════════════════════════════════════════
# B. Compute bearing from island center to each ahu
# ═══════════════════════════════════════════════════════════════════════

print(f"\n--- Ahu bearings from island center ---")
print(f"{'Ahu':25s}  {'Bearing':>8s}  {'Dist(km)':>9s}  {'ΔGC1':>6s}  {'ΔGC2':>6s}")
print(f"{'-'*25}  {'-'*8}  {'-'*9}  {'-'*6}  {'-'*6}")

ahu_bearings = []
ahu_dists = []
ahu_names = []

for ahu in AHU_DATA:
    b = initial_bearing(EI_CENTER_LAT, EI_CENTER_LON, ahu["lat"], ahu["lon"])
    d = haversine_km(EI_CENTER_LAT, EI_CENTER_LON, ahu["lat"], ahu["lon"])
    diff1 = angular_diff(b, gc_crossing_1)
    diff2 = angular_diff(b, gc_crossing_2)

    ahu_bearings.append(b)
    ahu_dists.append(d)
    ahu_names.append(ahu["name"])

    print(f"  {ahu['name']:23s}  {b:>8.1f}°  {d:>9.2f}  {diff1:>6.1f}  {diff2:>6.1f}")

ahu_bearings = np.array(ahu_bearings)
ahu_dists = np.array(ahu_dists)
n_ahu = len(ahu_bearings)


# ═══════════════════════════════════════════════════════════════════════
# C. Density by bearing sector
# ═══════════════════════════════════════════════════════════════════════

print(f"\n--- Ahu density by 30° bearing sector ---")

n_sectors = 12
sector_width = 360 / n_sectors
sector_counts = np.zeros(n_sectors, dtype=int)
sector_labels = []

for i in range(n_sectors):
    lo = i * sector_width
    hi = lo + sector_width
    sector_labels.append(f"{lo:.0f}-{hi:.0f}°")
    for b in ahu_bearings:
        if lo <= b < hi:
            sector_counts[i] += 1

# Which sectors contain the GC crossings?
gc_sector_1 = int(gc_crossing_1 // sector_width) % n_sectors
gc_sector_2 = int(gc_crossing_2 // sector_width) % n_sectors

expected = n_ahu / n_sectors

print(f"\n  {'Sector':12s}  {'Count':>6s}  {'Ratio':>6s}  {'Note':>10s}")
print(f"  {'-'*12}  {'-'*6}  {'-'*6}  {'-'*10}")
for i in range(n_sectors):
    ratio = sector_counts[i] / expected if expected > 0 else 0
    note = ""
    if i == gc_sector_1:
        note = "← GC entry"
    elif i == gc_sector_2:
        note = "← GC exit"
    print(f"  {sector_labels[i]:12s}  {sector_counts[i]:>6d}  {ratio:>6.2f}  {note:>10s}")

# Ahu count in GC sectors vs others
gc_sector_count = sector_counts[gc_sector_1] + sector_counts[gc_sector_2]
other_sector_count = n_ahu - gc_sector_count
gc_sector_expected = 2 * expected
other_expected = n_ahu - gc_sector_expected

print(f"\n  GC crossing sectors ({sector_labels[gc_sector_1]}, {sector_labels[gc_sector_2]}): "
      f"{gc_sector_count} ahu (expected {gc_sector_expected:.1f})")
print(f"  Other sectors: {other_sector_count} ahu (expected {other_expected:.1f})")
gc_enrichment = gc_sector_count / gc_sector_expected if gc_sector_expected > 0 else 0
print(f"  GC sector enrichment: {gc_enrichment:.2f}x")


# ═══════════════════════════════════════════════════════════════════════
# D. Circular statistics — Rayleigh and V-test
# ═══════════════════════════════════════════════════════════════════════

print(f"\n--- Circular statistics ---")

# Angular offsets from nearest GC crossing
offsets = []
for b in ahu_bearings:
    d1 = angular_diff(b, gc_crossing_1)
    d2 = angular_diff(b, gc_crossing_2)
    offsets.append(min(d1, d2))
offsets = np.array(offsets)

print(f"\n  Angular offsets from nearest GC crossing point:")
print(f"    Mean offset: {np.mean(offsets):.1f}°")
print(f"    Median offset: {np.median(offsets):.1f}°")
print(f"    Std offset: {np.std(offsets):.1f}°")
print(f"    If uniform around circle, expected mean offset: 45.0°")

# Rayleigh test: do ahu cluster in any direction?
# Use doubled angles for bimodal test (cluster at both GC crossings)
doubled_rad = 2 * np.radians(ahu_bearings)
C = np.mean(np.cos(doubled_rad))
S = np.mean(np.sin(doubled_rad))
R_bar = np.sqrt(C**2 + S**2)
mean_dir_doubled = np.degrees(np.arctan2(S, C)) % 360
mean_dir = (mean_dir_doubled / 2) % 180

n = len(doubled_rad)
Z = R_bar**2 * n
rayleigh_p = np.exp(-Z) * (1 + (2*Z - Z**2)/(4*n))
rayleigh_p = max(rayleigh_p, 0)

print(f"\n  Rayleigh test (bimodal / doubled angles):")
print(f"    R̄ = {R_bar:.4f}")
print(f"    Mean direction (axial): {mean_dir:.1f}°")
print(f"    Z = {Z:.2f}")
print(f"    p-value ≈ {rayleigh_p:.6f}")
print(f"    {'→ Significant clustering' if rayleigh_p < 0.05 else '→ No significant clustering'}")

# V-test directed at GC crossing bearings
gc_doubled = 2 * math.radians(gc_crossing_1)
V = R_bar * math.cos(np.arctan2(S, C) - gc_doubled)
u = V * math.sqrt(2 * n)
v_p = 1 - stats.norm.cdf(u)

print(f"\n  V-test (directed at GC bearing {gc_crossing_1:.1f}°):")
print(f"    V = {V:.4f}")
print(f"    u = {u:.2f}")
print(f"    p-value = {v_p:.6f}")
print(f"    {'→ Ahu cluster toward GC crossing' if v_p < 0.05 else '→ No directed clustering toward GC'}")


# ═══════════════════════════════════════════════════════════════════════
# E. Monte Carlo: random bearing test
# ═══════════════════════════════════════════════════════════════════════

print(f"\n--- Monte Carlo: {N_MC} random test bearings ---")

np.random.seed(42)
mc_mean_offsets = np.zeros(N_MC)

for i in range(N_MC):
    # Random pair of antipodal crossing bearings
    random_b1 = np.random.uniform(0, 180)
    random_b2 = (random_b1 + 180) % 360

    offsets_mc = []
    for b in ahu_bearings:
        d1 = angular_diff(b, random_b1)
        d2 = angular_diff(b, random_b2)
        offsets_mc.append(min(d1, d2))
    mc_mean_offsets[i] = np.mean(offsets_mc)

actual_mean_offset = np.mean(offsets)
mc_p = np.mean(mc_mean_offsets <= actual_mean_offset)

print(f"  Actual mean offset from GC crossing: {actual_mean_offset:.1f}°")
print(f"  MC mean offset (random crossings): {np.mean(mc_mean_offsets):.1f}° ± {np.std(mc_mean_offsets):.1f}°")
print(f"  p-value (actual ≤ MC): {mc_p:.4f}")
print(f"  {'→ Ahu are significantly closer to GC crossings' if mc_p < 0.05 else '→ No significant clustering at GC crossings'}")


# ═══════════════════════════════════════════════════════════════════════
# Generate map
# ═══════════════════════════════════════════════════════════════════════

print(f"\nGenerating Easter Island map...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# --- Left: Map of Easter Island with ahu and GC line ---
ahu_lons = [a["lon"] for a in AHU_DATA]
ahu_lats = [a["lat"] for a in AHU_DATA]

ax1.scatter(ahu_lons, ahu_lats, s=60, c="red", edgecolors="black",
            linewidth=0.5, zorder=5, label=f"Ahu (n={n_ahu})")

# Draw GC line through the island
# Extend the GC bearing as a line segment ~20km in each direction
for bearing, linestyle in [(gc_bearing_1, '-'), (gc_bearing_2, '--')]:
    b_rad = math.radians(bearing)
    for dist_km in [0.5, 1, 2, 4, 8, 12, 16, 20]:
        # Destination point at distance along bearing
        lat1 = math.radians(EI_CENTER_LAT)
        lon1 = math.radians(EI_CENTER_LON)
        d_r = dist_km / R_EARTH
        lat2 = math.asin(math.sin(lat1)*math.cos(d_r) +
                          math.cos(lat1)*math.sin(d_r)*math.cos(b_rad))
        lon2 = lon1 + math.atan2(math.sin(b_rad)*math.sin(d_r)*math.cos(lat1),
                                   math.cos(d_r) - math.sin(lat1)*math.sin(lat2))

gc_line_lats = []
gc_line_lons = []
for sign in [-1, 1]:
    for dist_km in np.linspace(0, 20, 50):
        bearing_use = gc_bearing_1 if sign > 0 else gc_bearing_2
        b_rad = math.radians(bearing_use)
        lat1 = math.radians(EI_CENTER_LAT)
        lon1 = math.radians(EI_CENTER_LON)
        d_r = dist_km / R_EARTH
        lat2 = math.asin(math.sin(lat1)*math.cos(d_r) +
                          math.cos(lat1)*math.sin(d_r)*math.cos(b_rad))
        lon2 = lon1 + math.atan2(math.sin(b_rad)*math.sin(d_r)*math.cos(lat1),
                                   math.cos(d_r) - math.sin(lat1)*math.sin(lat2))
        gc_line_lats.append(math.degrees(lat2))
        gc_line_lons.append(math.degrees(lon2))

ax1.plot(gc_line_lons, gc_line_lats, 'b-', linewidth=2, alpha=0.6,
         label=f'Great Circle ({gc_bearing_1:.0f}°/{gc_bearing_2:.0f}°)', zorder=3)

# Mark center
ax1.scatter([EI_CENTER_LON], [EI_CENTER_LAT], s=40, c="blue", marker="+",
            linewidth=2, zorder=6)

# Labels for major ahu
for a in AHU_DATA:
    if a["n_moai"] >= 5 or a["name"] in ["Ahu Tongariki", "Ahu Akivi", "Ahu Nau Nau",
                                            "Ahu Vinapu", "Ahu Te Pito Kura"]:
        ax1.annotate(a["name"].replace("Ahu ", ""), (a["lon"], a["lat"]),
                     xytext=(3, 3), textcoords="offset points", fontsize=6, zorder=7)

ax1.set_xlabel("Longitude")
ax1.set_ylabel("Latitude")
ax1.set_title("Easter Island: Ahu and Great Circle Path")
ax1.legend(fontsize=8, loc="lower left")
ax1.grid(True, alpha=0.3)
ax1.set_aspect("equal")

# --- Right: Polar plot of ahu density by bearing ---
theta = np.radians(np.arange(0, 360, sector_width) + sector_width/2)
width = np.radians(sector_width)

ax2 = fig.add_subplot(122, projection='polar')
ax2.set_theta_zero_location('N')
ax2.set_theta_direction(-1)

bars = ax2.bar(theta, sector_counts, width=width, alpha=0.6,
               edgecolor='black', linewidth=0.5)

# Color GC sectors differently
for i, bar in enumerate(bars):
    if i == gc_sector_1 or i == gc_sector_2:
        bar.set_facecolor('red')
        bar.set_alpha(0.8)
    else:
        bar.set_facecolor('steelblue')

# Mark GC crossing bearings
ax2.axvline(math.radians(gc_crossing_1), color='red', linewidth=2, linestyle='--',
            label=f'GC crossing ({gc_crossing_1:.0f}°)')
ax2.axvline(math.radians(gc_crossing_2), color='red', linewidth=2, linestyle='--',
            label=f'GC crossing ({gc_crossing_2:.0f}°)')

# Expected uniform level
ax2.axhline(expected, color='gray', linewidth=1, linestyle=':', alpha=0.7)

ax2.set_title("Ahu Density by Bearing from Island Center\n(red sectors = GC crossings)", pad=20)
ax2.legend(fontsize=7, loc='upper right', bbox_to_anchor=(1.3, 1.1))

plt.tight_layout()
fig.savefig("easter_island_circle_overlay.png", dpi=150)
print("Saved easter_island_circle_overlay.png")


# ═══════════════════════════════════════════════════════════════════════
# Save JSON results
# ═══════════════════════════════════════════════════════════════════════

orientation_test = {
    "analysis": "Easter Island ahu distribution vs Great Circle",
    "island_center": {"lat": EI_CENTER_LAT, "lon": EI_CENTER_LON},
    "distance_to_circle_km": round(dist_to_circle, 1),
    "gc_bearing": {
        "bearing_1": round(gc_bearing_1, 2),
        "bearing_2": round(gc_bearing_2, 2),
    },
    "n_ahu": n_ahu,
    "rayleigh_test": {
        "R_bar": round(float(R_bar), 4),
        "mean_direction_axial": round(float(mean_dir), 2),
        "Z": round(float(Z), 2),
        "p_value": round(float(rayleigh_p), 6),
        "significant": bool(rayleigh_p < 0.05),
    },
    "v_test": {
        "test_direction": round(gc_crossing_1, 2),
        "V": round(float(V), 4),
        "u": round(float(u), 2),
        "p_value": round(float(v_p), 6),
        "significant": bool(v_p < 0.05),
    },
    "monte_carlo": {
        "n_trials": N_MC,
        "actual_mean_offset_deg": round(float(actual_mean_offset), 2),
        "mc_mean_offset_deg": round(float(np.mean(mc_mean_offsets)), 2),
        "mc_std_offset_deg": round(float(np.std(mc_mean_offsets)), 2),
        "p_value": round(float(mc_p), 4),
        "significant": bool(mc_p < 0.05),
    },
    "gc_sector_enrichment": round(gc_enrichment, 3),
}

with open("ahu_orientation_test.json", "w") as f:
    json.dump(orientation_test, f, indent=2)

print("Saved ahu_orientation_test.json")

# Density by bearing
density_data = {
    "sector_width_deg": sector_width,
    "sectors": [
        {
            "label": sector_labels[i],
            "center_deg": round(i * sector_width + sector_width/2, 1),
            "count": int(sector_counts[i]),
            "is_gc_crossing_sector": bool(i == gc_sector_1 or i == gc_sector_2),
        }
        for i in range(n_sectors)
    ],
    "expected_per_sector": round(expected, 2),
}

with open("ahu_density_by_bearing.json", "w") as f:
    json.dump(density_data, f, indent=2)

print("Saved ahu_density_by_bearing.json")
