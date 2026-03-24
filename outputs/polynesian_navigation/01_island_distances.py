#!/usr/bin/env python3
"""
Analysis 1a: Pacific Island Distances to the Great Circle
==========================================================
Directive 09, Script 01

Compute the distance from every major Polynesian island group to the
Great Circle. Identify which islands fall within each distance band
(50, 100, 200, 500 km). Produce a CSV and a Pacific map.

Output:
  - island_distances.csv
  - pacific_circle_map.png
"""

import csv
import json
import math
import sys
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ─── Import data module ──────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from polynesian_voyaging_data import VOYAGING_ROUTES, COLONIZATION_CHRONOLOGY

# ─── Constants ────────────────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
R_EARTH = 6371.0
QUARTER_CIRC = R_EARTH * math.pi / 2
THRESHOLDS_KM = [50, 100, 200, 500]

# ─── Island groups (comprehensive list) ──────────────────────────────
ISLAND_GROUPS = [
    {"name": "Easter Island (Rapa Nui)", "lat": -27.11, "lon": -109.35},
    {"name": "Pitcairn Islands", "lat": -25.07, "lon": -130.10},
    {"name": "Henderson Island", "lat": -24.37, "lon": -128.32},
    {"name": "Mangareva / Gambier", "lat": -23.13, "lon": -134.97},
    {"name": "Marquesas Islands", "lat": -9.43, "lon": -140.07},
    {"name": "Tuamotu Archipelago", "lat": -17.00, "lon": -144.00},
    {"name": "Society Islands", "lat": -17.32, "lon": -149.75},
    {"name": "Austral Islands", "lat": -22.44, "lon": -151.35},
    {"name": "Cook Islands", "lat": -19.00, "lon": -159.50},
    {"name": "Samoa", "lat": -13.76, "lon": -172.10},
    {"name": "Tonga", "lat": -20.42, "lon": -175.20},
    {"name": "Fiji", "lat": -17.71, "lon": 178.07},
    {"name": "Hawaiian Islands", "lat": 20.79, "lon": -156.33},
    {"name": "New Zealand (Aotearoa)", "lat": -41.27, "lon": 174.78},
    {"name": "Niue", "lat": -19.05, "lon": -169.87},
    {"name": "Tokelau", "lat": -9.17, "lon": -171.82},
    {"name": "Rapa Iti", "lat": -27.62, "lon": -144.34},
    {"name": "Line Islands (Kiribati)", "lat": 1.87, "lon": -157.47},
    {"name": "Phoenix Islands", "lat": -3.85, "lon": -171.73},
    {"name": "Wallis & Futuna", "lat": -13.28, "lon": -176.18},
]


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlat = la2 - la1
    dlon = lo2 - lo1
    a = math.sin(dlat/2)**2 + math.cos(la1)*math.cos(la2)*math.sin(dlon/2)**2
    return R_EARTH * 2 * math.asin(math.sqrt(min(1.0, a)))


def dist_from_great_circle(lat, lon):
    dist_to_pole = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(dist_to_pole - QUARTER_CIRC)


def bearing_at_point(point_lat, point_lon):
    plat = math.radians(POLE_LAT)
    plon = math.radians(POLE_LON)
    qlat = math.radians(point_lat)
    qlon = math.radians(point_lon)
    dlon = plon - qlon
    x = math.sin(dlon) * math.cos(plat)
    y = (math.cos(qlat) * math.sin(plat) -
         math.sin(qlat) * math.cos(plat) * math.cos(dlon))
    bearing_to_pole = math.degrees(math.atan2(x, y)) % 360
    return (bearing_to_pole + 90) % 360, (bearing_to_pole - 90) % 360


# ═══════════════════════════════════════════════════════════════════════
# Compute distances
# ═══════════════════════════════════════════════════════════════════════

print("=" * 70)
print("PACIFIC ISLAND DISTANCES TO THE GREAT CIRCLE")
print("=" * 70)
print(f"\nPole: ({POLE_LAT}, {POLE_LON})")

results = []
for ig in ISLAND_GROUPS:
    d = dist_from_great_circle(ig["lat"], ig["lon"])
    b1, b2 = bearing_at_point(ig["lat"], ig["lon"])
    bands = [t for t in THRESHOLDS_KM if d <= t]
    results.append({
        "name": ig["name"],
        "lat": ig["lat"],
        "lon": ig["lon"],
        "distance_km": round(d, 1),
        "bearing_east": round(b1, 2),
        "bearing_west": round(b2, 2),
        "within_bands_km": bands,
    })

results.sort(key=lambda x: x["distance_km"])

print(f"\n{'Island Group':40s}  {'Dist (km)':>10s}  {'Bands':>20s}")
print(f"{'-'*40}  {'-'*10}  {'-'*20}")
for r in results:
    bands_str = ", ".join(str(b) for b in r["within_bands_km"]) if r["within_bands_km"] else "—"
    print(f"  {r['name']:38s}  {r['distance_km']:>10.1f}  {bands_str:>20s}")

# Band summary
print(f"\n--- Islands within each distance band ---")
for t in THRESHOLDS_KM:
    within = [r["name"] for r in results if r["distance_km"] <= t]
    print(f"  ≤{t:4d} km: {len(within)} islands — {', '.join(within) if within else 'none'}")


# ═══════════════════════════════════════════════════════════════════════
# Save CSV
# ═══════════════════════════════════════════════════════════════════════

csv_path = "island_distances.csv"
with open(csv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "name", "lat", "lon", "distance_km", "bearing_east", "bearing_west",
    ])
    writer.writeheader()
    for r in results:
        writer.writerow({k: r[k] for k in writer.fieldnames})

print(f"\nSaved {csv_path}")


# ═══════════════════════════════════════════════════════════════════════
# Pacific Circle Map
# ═══════════════════════════════════════════════════════════════════════

print("\nGenerating Pacific circle map...")

fig, ax = plt.subplots(1, 1, figsize=(16, 10))

# Draw the Great Circle across the Pacific
# Generate circle points by rotating around the pole
circle_lats = []
circle_lons = []
for theta_deg in range(360):
    theta = math.radians(theta_deg)
    # Point on great circle at angular position theta
    # Circle is the locus of points at 90° from the pole
    plat = math.radians(POLE_LAT)
    plon = math.radians(POLE_LON)

    # Unit vector of pole
    px = math.cos(plat) * math.cos(plon)
    py = math.cos(plat) * math.sin(plon)
    pz = math.sin(plat)

    # Two orthogonal vectors in the plane perpendicular to pole
    # v1 = pole × (0,0,1) normalized, v2 = pole × v1
    if abs(pz) < 0.999:
        v1x = -py
        v1y = px
        v1z = 0
    else:
        v1x = 1
        v1y = 0
        v1z = 0
    norm = math.sqrt(v1x**2 + v1y**2 + v1z**2)
    v1x /= norm; v1y /= norm; v1z /= norm

    v2x = py*v1z - pz*v1y
    v2y = pz*v1x - px*v1z
    v2z = px*v1y - py*v1x

    # Point on circle
    qx = v1x * math.cos(theta) + v2x * math.sin(theta)
    qy = v1y * math.cos(theta) + v2y * math.sin(theta)
    qz = v1z * math.cos(theta) + v2z * math.sin(theta)

    lat_q = math.degrees(math.asin(max(-1, min(1, qz))))
    lon_q = math.degrees(math.atan2(qy, qx))

    circle_lats.append(lat_q)
    circle_lons.append(lon_q)

# Filter to Pacific region and handle dateline wrapping
# Shift longitudes to center on Pacific (0-360 system centered at 180)
def to_pacific_lon(lon):
    """Convert -180..180 to a Pacific-centered view (roughly 100E to 60W)."""
    if lon < -30:
        return lon + 360
    return lon

pac_circle_x = [to_pacific_lon(lo) for lo in circle_lons]
pac_circle_y = circle_lats

# Sort by x to avoid jumps (approximately)
pairs = sorted(zip(pac_circle_x, pac_circle_y))
# Break into segments where x-gap > 30°
segments_x = [[]]
segments_y = [[]]
for i, (x, y) in enumerate(pairs):
    if i > 0 and abs(x - pairs[i-1][0]) > 30:
        segments_x.append([])
        segments_y.append([])
    segments_x[-1].append(x)
    segments_y[-1].append(y)

for sx, sy in zip(segments_x, segments_y):
    ax.plot(sx, sy, 'r-', linewidth=2, alpha=0.7, zorder=2)

# Plot island groups
for r in results:
    px = to_pacific_lon(r["lon"])
    py = r["lat"]
    if r["distance_km"] <= 200:
        color = "red"
        size = 120
    elif r["distance_km"] <= 500:
        color = "orange"
        size = 80
    else:
        color = "dodgerblue"
        size = 60

    ax.scatter(px, py, s=size, c=color, edgecolors="black", linewidth=0.5, zorder=5)
    # Label offset
    offset_x = 1.5
    offset_y = 1.5
    ax.annotate(r["name"], (px, py), xytext=(px + offset_x, py + offset_y),
                fontsize=6.5, ha="left", va="bottom",
                arrowprops=dict(arrowstyle="-", lw=0.3, color="gray"),
                zorder=6)

# Distance bands around the circle (200 km and 500 km)
# (simplified as visual reference)

ax.set_xlim(100, 310)
ax.set_ylim(-55, 35)
ax.set_xlabel("Longitude (Pacific-centered)", fontsize=11)
ax.set_ylabel("Latitude", fontsize=11)
ax.set_title("The Great Circle Across the Pacific\nwith Polynesian Island Groups", fontsize=14)

# Custom x-tick labels
xticks = [100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300]
xlabels = []
for x in xticks:
    if x <= 180:
        xlabels.append(f"{x}°E")
    else:
        xlabels.append(f"{360-x}°W")
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='red', linewidth=2, label='Great Circle'),
    plt.scatter([], [], s=120, c='red', edgecolors='black', label='≤200 km from circle'),
    plt.scatter([], [], s=80, c='orange', edgecolors='black', label='200–500 km'),
    plt.scatter([], [], s=60, c='dodgerblue', edgecolors='black', label='>500 km'),
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=9)

ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig("pacific_circle_map.png", dpi=150)
print("Saved pacific_circle_map.png")

# ═══════════════════════════════════════════════════════════════════════
# Save JSON summary
# ═══════════════════════════════════════════════════════════════════════

summary = {
    "pole": {"lat": POLE_LAT, "lon": POLE_LON},
    "n_island_groups": len(results),
    "band_counts": {str(t): len([r for r in results if r["distance_km"] <= t])
                    for t in THRESHOLDS_KM},
    "closest_islands": [
        {"name": r["name"], "distance_km": r["distance_km"]}
        for r in results[:5]
    ],
    "results": results,
}

with open("island_distances.json", "w") as f:
    json.dump(summary, f, indent=2)

print("Saved island_distances.json")
