#!/usr/bin/env python3
"""
Analysis 1: Great Circle bearing at Nazca
==========================================
Compute the bearing (azimuth) of Alison's Great Circle as it passes through
the Nazca Pampa region.

The Great Circle is defined by its pole at (59.682122, -138.646087).
At any point on the circle, the bearing is perpendicular to the direction
toward the pole.
"""

import json
import math
import numpy as np

# ─── Constants ──────────────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
R_EARTH = 6371.0  # km

# 5 reference points across the Nazca Pampa (~50km extent)
NAZCA_POINTS = {
    "Nazca_Pampa_center":     (-14.735, -75.130),
    "Nazca_Pampa_NW":         (-14.60,  -75.25),
    "Nazca_Pampa_SE":         (-14.85,  -75.00),
    "Nazca_Pampa_NE":         (-14.60,  -75.00),
    "Nazca_Pampa_SW":         (-14.85,  -75.25),
    "Cahuachi_temple":        (-14.817, -75.120),
    "Hummingbird_geoglyph":   (-14.692, -75.149),
}

OUTPUT_DIR = "."


def bearing_of_great_circle_at_point(pole_lat, pole_lon, point_lat, point_lon):
    """
    Compute the bearing of a great circle (defined by its pole) at a given point.
    Returns two bearings (the circle goes both ways), both in [0, 360).
    """
    plat = math.radians(pole_lat)
    plon = math.radians(pole_lon)
    qlat = math.radians(point_lat)
    qlon = math.radians(point_lon)

    # Direction from point to pole
    dlon = plon - qlon
    x = math.sin(dlon) * math.cos(plat)
    y = (math.cos(qlat) * math.sin(plat) -
         math.sin(qlat) * math.cos(plat) * math.cos(dlon))
    bearing_to_pole = math.degrees(math.atan2(x, y)) % 360

    # Circle bearing is perpendicular to direction to pole
    bearing1 = (bearing_to_pole + 90) % 360
    bearing2 = (bearing_to_pole - 90) % 360
    return bearing1, bearing2


def dist_from_great_circle(lat, lon, pole_lat, pole_lon):
    """Distance (km) from a point to the great circle defined by the pole."""
    lat1, lon1 = math.radians(lat), math.radians(lon)
    lat2, lon2 = math.radians(pole_lat), math.radians(pole_lon)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    dist_to_pole = R_EARTH * 2 * math.asin(math.sqrt(min(1.0, a)))
    quarter_circ = R_EARTH * math.pi / 2
    return abs(dist_to_pole - quarter_circ)


# ═══════════════════════════════════════════════════════════════════
# Compute bearing at each reference point
# ═══════════════════════════════════════════════════════════════════

print("=" * 65)
print("GREAT CIRCLE BEARING AT NAZCA")
print("=" * 65)
print(f"\nPole: ({POLE_LAT}, {POLE_LON})")
print()

results = {}
for name, (lat, lon) in NAZCA_POINTS.items():
    b1, b2 = bearing_of_great_circle_at_point(POLE_LAT, POLE_LON, lat, lon)
    dist = dist_from_great_circle(lat, lon, POLE_LAT, POLE_LON)

    # Normalize: pick the bearing that points roughly E (0-180)
    eastward = b1 if 0 <= b1 <= 180 else b2
    westward = b2 if 0 <= b1 <= 180 else b1

    results[name] = {
        "lat": lat,
        "lon": lon,
        "bearing_east": round(eastward, 3),
        "bearing_west": round(westward, 3),
        "distance_to_circle_km": round(dist, 2),
    }

    print(f"  {name:30s}  bearing = {eastward:7.2f}° / {westward:7.2f}°  "
          f"dist_to_circle = {dist:.1f} km")

# Summary
bearings_east = [r["bearing_east"] for r in results.values()]
mean_bearing = np.mean(bearings_east)
std_bearing = np.std(bearings_east)
min_dist = min(r["distance_to_circle_km"] for r in results.values())

print(f"\n  Mean bearing (eastward):  {mean_bearing:.2f}° ± {std_bearing:.2f}°")
print(f"  Mean bearing (westward):  {(mean_bearing + 180) % 360:.2f}°")
print(f"  Closest approach to circle: {min_dist:.1f} km")

# The Great Circle bearing at Nazca for the orientation test
# Use the center point value
center = results["Nazca_Pampa_center"]
print(f"\n  *** Test bearing (at Pampa center): {center['bearing_east']:.2f}° / "
      f"{center['bearing_west']:.2f}° ***")

# ═══════════════════════════════════════════════════════════════════
# Cardinal direction analysis
# ═══════════════════════════════════════════════════════════════════

print("\n" + "=" * 65)
print("CARDINAL DIRECTION ANALYSIS")
print("=" * 65)

# What direction is the bearing?
b = center["bearing_east"]
if 67.5 <= b < 112.5:
    direction = "roughly E-W"
elif 22.5 <= b < 67.5:
    direction = "roughly NE-SW"
elif 112.5 <= b < 157.5:
    direction = "roughly SE-NW"
elif 157.5 <= b < 202.5:
    direction = "roughly N-S"
else:
    direction = f"bearing {b:.1f}°"

print(f"\n  The Great Circle runs {direction} through Nazca.")
print(f"  Bearing: {b:.2f}° (east) / {center['bearing_west']:.2f}° (west)")
print(f"  This is {b:.1f}° east of north, or {90-b:.1f}° north of east.")

# Known Nazca orientation references
print("\n  Comparison with known Nazca orientations:")
print(f"    Aveni's bimodal peaks: ~10° and ~80° (from north)")
print(f"    Circle bearing:        {b:.1f}°")
diff_10 = min(abs(b - 10), abs(b - 190))
diff_80 = min(abs(b - 80), abs(b - 260))
print(f"    Offset from ~10° peak: {diff_10:.1f}°")
print(f"    Offset from ~80° peak: {diff_80:.1f}°")

# ═══════════════════════════════════════════════════════════════════
# Save results
# ═══════════════════════════════════════════════════════════════════

output = {
    "pole": {"lat": POLE_LAT, "lon": POLE_LON},
    "reference_points": results,
    "summary": {
        "mean_bearing_east": round(mean_bearing, 3),
        "mean_bearing_west": round((mean_bearing + 180) % 360, 3),
        "bearing_std": round(std_bearing, 3),
        "closest_approach_km": round(min_dist, 2),
        "cardinal_direction": direction,
    },
    "test_bearing": {
        "location": "Nazca_Pampa_center",
        "bearing_east": center["bearing_east"],
        "bearing_west": center["bearing_west"],
    },
    "aveni_comparison": {
        "peak_1_deg": 10,
        "peak_2_deg": 80,
        "offset_from_peak_1": round(diff_10, 2),
        "offset_from_peak_2": round(diff_80, 2),
        "note": "Aveni peaks are approximate; exact values from books not digitally available",
    }
}

with open("circle_bearing_at_nazca.json", "w") as f:
    json.dump(output, f, indent=2)

print(f"\nSaved to circle_bearing_at_nazca.json")
