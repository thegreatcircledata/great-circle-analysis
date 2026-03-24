#!/usr/bin/env python3
"""
Analysis 2: Polynesian Navigation Route Geometry
==================================================
Directive 09, Script 03

For each documented/reconstructed voyaging route:
  1. Fit a great circle arc to the route waypoints
  2. Compute RMS deviation from the best-fit great circle
  3. Compare: is the route closer to a great circle than a rhumb line?
  4. Compute angular separation between each route's best-fit GC and THE Great Circle
  5. Star compass bearing analysis at Easter Island

Output:
  - route_geometries.csv
  - route_circle_comparison.json
  - star_compass_bearing.json
"""

import csv
import json
import math
import os
import sys

import numpy as np

# ─── Import data ─────────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from polynesian_voyaging_data import VOYAGING_ROUTES

# ─── Constants ────────────────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
R_EARTH = 6371.0
QUARTER_CIRC = R_EARTH * math.pi / 2


def to_cartesian(lat, lon):
    """Convert lat/lon to unit sphere cartesian."""
    la, lo = math.radians(lat), math.radians(lon)
    return np.array([math.cos(la)*math.cos(lo),
                     math.cos(la)*math.sin(lo),
                     math.sin(la)])


def to_latlon(xyz):
    """Convert unit cartesian to lat/lon."""
    x, y, z = xyz
    lat = math.degrees(math.asin(max(-1, min(1, z))))
    lon = math.degrees(math.atan2(y, x))
    return lat, lon


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlat = la2 - la1
    dlon = lo2 - lo1
    a = math.sin(dlat/2)**2 + math.cos(la1)*math.cos(la2)*math.sin(dlon/2)**2
    return R_EARTH * 2 * math.asin(math.sqrt(min(1.0, a)))


def dist_from_gc_pole(lat, lon, pole_lat, pole_lon):
    """Distance (km) from a point to a great circle defined by its pole."""
    d = haversine_km(lat, lon, pole_lat, pole_lon)
    return abs(d - QUARTER_CIRC)


def fit_great_circle_to_points(waypoints):
    """
    Fit the best great circle (pole) to a set of waypoints.
    Uses SVD: the pole is the singular vector with smallest singular value
    of the matrix of cartesian waypoint positions.
    """
    cart = np.array([to_cartesian(lat, lon) for lat, lon in waypoints])
    # SVD
    U, S, Vt = np.linalg.svd(cart, full_matrices=True)
    # Pole is the row of Vt corresponding to smallest singular value
    pole_xyz = Vt[-1]
    pole_lat, pole_lon = to_latlon(pole_xyz)
    return pole_lat, pole_lon


def rms_deviation_from_gc(waypoints, pole_lat, pole_lon):
    """RMS deviation (km) of waypoints from a great circle."""
    deviations = []
    for lat, lon in waypoints:
        d = dist_from_gc_pole(lat, lon, pole_lat, pole_lon)
        deviations.append(d ** 2)
    return math.sqrt(np.mean(deviations))


def rhumb_line_deviation(waypoints):
    """
    Compute the RMS deviation of waypoints from a rhumb line
    (constant bearing) connecting first and last point.
    """
    lat1, lon1 = waypoints[0]
    lat2, lon2 = waypoints[-1]

    # Rhumb bearing
    dlon = math.radians(lon2 - lon1)
    dphi = math.log(math.tan(math.pi/4 + math.radians(lat2)/2) /
                     math.tan(math.pi/4 + math.radians(lat1)/2))
    bearing = math.atan2(dlon, dphi)

    # For each intermediate point, compute distance from rhumb line
    # Approximate: project point onto the rhumb line and compute perpendicular distance
    deviations = []
    total_dist = haversine_km(lat1, lon1, lat2, lon2)

    for lat, lon in waypoints[1:-1]:
        # Distance from point to first waypoint along rhumb bearing
        d_start = haversine_km(lat1, lon1, lat, lon)
        d_end = haversine_km(lat, lon, lat2, lon2)

        # Fraction along route
        if total_dist > 0:
            frac = d_start / (d_start + d_end)
        else:
            frac = 0.5

        # Expected position on rhumb line at this fraction
        exp_lat = lat1 + frac * (lat2 - lat1)
        exp_lon = lon1 + frac * (lon2 - lon1)

        dev = haversine_km(lat, lon, exp_lat, exp_lon)
        deviations.append(dev ** 2)

    if deviations:
        return math.sqrt(np.mean(deviations))
    return 0.0


def angular_separation_of_poles(pole1_lat, pole1_lon, pole2_lat, pole2_lon):
    """
    Angular separation (degrees) between two great circles defined by their poles.
    This is the angle between the pole vectors (or 180 - that if >90).
    """
    p1 = to_cartesian(pole1_lat, pole1_lon)
    p2 = to_cartesian(pole2_lat, pole2_lon)
    cos_angle = np.dot(p1, p2)
    cos_angle = max(-1, min(1, cos_angle))
    angle = math.degrees(math.acos(abs(cos_angle)))  # abs because poles can be antipodal
    return angle


# ═══════════════════════════════════════════════════════════════════════
# Analyze each route
# ═══════════════════════════════════════════════════════════════════════

print("=" * 75)
print("POLYNESIAN NAVIGATION ROUTE GEOMETRY")
print("=" * 75)
print(f"\nGreat Circle pole: ({POLE_LAT}, {POLE_LON})")
print(f"Routes to analyze: {len(VOYAGING_ROUTES)}")

route_results = []

for route in VOYAGING_ROUTES:
    name = route["name"]
    wps = route["waypoints"]
    n_wp = len(wps)

    if n_wp < 2:
        continue

    # Total route distance
    total_dist = sum(haversine_km(wps[i][0], wps[i][1], wps[i+1][0], wps[i+1][1])
                     for i in range(n_wp - 1))

    # Fit great circle to waypoints
    fit_pole_lat, fit_pole_lon = fit_great_circle_to_points(wps)

    # RMS deviation from best-fit GC
    rms_gc = rms_deviation_from_gc(wps, fit_pole_lat, fit_pole_lon)

    # RMS deviation from rhumb line
    rms_rhumb = rhumb_line_deviation(wps)

    # Angular separation between route's GC and THE Great Circle
    ang_sep = angular_separation_of_poles(fit_pole_lat, fit_pole_lon, POLE_LAT, POLE_LON)

    # Is the route closer to a great circle than a rhumb line?
    gc_preferred = rms_gc < rms_rhumb if n_wp > 2 else None

    result = {
        "name": name,
        "n_waypoints": n_wp,
        "route_type": route["route_type"],
        "total_distance_km": round(total_dist, 0),
        "fit_pole_lat": round(fit_pole_lat, 4),
        "fit_pole_lon": round(fit_pole_lon, 4),
        "rms_from_gc_km": round(rms_gc, 1),
        "rms_from_rhumb_km": round(rms_rhumb, 1),
        "gc_preferred": gc_preferred,
        "angular_sep_from_great_circle_deg": round(ang_sep, 2),
        "within_5deg": ang_sep <= 5.0,
        "within_10deg": ang_sep <= 10.0,
    }
    route_results.append(result)

    print(f"\n  {name}")
    print(f"    Waypoints: {n_wp}, Total distance: {total_dist:.0f} km")
    print(f"    Best-fit GC pole: ({fit_pole_lat:.2f}, {fit_pole_lon:.2f})")
    print(f"    RMS from GC: {rms_gc:.1f} km, from rhumb: {rms_rhumb:.1f} km"
          + (f" → {'GC' if gc_preferred else 'rhumb'} preferred" if gc_preferred is not None else ""))
    print(f"    Angular sep from Great Circle: {ang_sep:.2f}°"
          + (" ← WITHIN 5°!" if ang_sep <= 5 else " ← within 10°" if ang_sep <= 10 else ""))


# ═══════════════════════════════════════════════════════════════════════
# Star compass bearing at key locations
# ═══════════════════════════════════════════════════════════════════════

print(f"\n{'='*75}")
print("STAR COMPASS BEARING ANALYSIS")
print(f"{'='*75}")

# Compute Great Circle bearing at key Pacific locations
key_locations = {
    "Easter Island": (-27.11, -109.35),
    "Tahiti": (-17.53, -149.57),
    "Raiatea": (-16.73, -151.44),
    "Nuku Hiva (Marquesas)": (-8.92, -140.07),
    "Hawaii (Big Island)": (19.72, -155.08),
    "Rarotonga": (-21.24, -159.77),
    "Tongatapu": (-21.13, -175.20),
}

# Star compass sectors (Hawaiian tradition, approximate)
STAR_SECTORS = [
    (0, 11.25, "N (Hoku Pa'a / Polaris)"),
    (11.25, 33.75, "NNE (Na Leo)"),
    (33.75, 56.25, "NE (Manu)"),
    (56.25, 67.5, "ENE (Noio)"),
    (67.5, 78.75, "E-NE (Hokule'a rising / Arcturus)"),
    (78.75, 90, "E (La / Sun equinox rise)"),
    (90, 101.25, "E (Hikianalia / Spica rising)"),
    (101.25, 123.75, "ESE-SE"),
    (123.75, 146.25, "SE (Ka Makau Nui / Scorpius)"),
    (146.25, 168.75, "SSE (Newe / S. Cross)"),
    (168.75, 191.25, "S"),
    (191.25, 213.75, "SSW"),
    (213.75, 236.25, "SW"),
    (236.25, 258.75, "WSW"),
    (258.75, 281.25, "W"),
    (281.25, 303.75, "WNW-NW"),
    (303.75, 326.25, "NW"),
    (326.25, 348.75, "NNW"),
    (348.75, 360, "N (Hoku Pa'a / Polaris)"),
]

star_compass_results = {}

for loc_name, (lat, lon) in key_locations.items():
    plat = math.radians(POLE_LAT)
    plon = math.radians(POLE_LON)
    qlat = math.radians(lat)
    qlon = math.radians(lon)
    dlon = plon - qlon
    x = math.sin(dlon) * math.cos(plat)
    y = (math.cos(qlat) * math.sin(plat) -
         math.sin(qlat) * math.cos(plat) * math.cos(dlon))
    bearing_to_pole = math.degrees(math.atan2(x, y)) % 360
    b1 = (bearing_to_pole + 90) % 360
    b2 = (bearing_to_pole - 90) % 360

    d = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    dist_to_circle = abs(d - QUARTER_CIRC)

    # Find star compass sector
    sector1 = "unknown"
    sector2 = "unknown"
    for lo, hi, name in STAR_SECTORS:
        if lo <= b1 < hi:
            sector1 = name
        if lo <= b2 < hi:
            sector2 = name

    star_compass_results[loc_name] = {
        "lat": lat,
        "lon": lon,
        "distance_to_circle_km": round(dist_to_circle, 1),
        "bearing_1": round(b1, 2),
        "bearing_2": round(b2, 2),
        "star_sector_1": sector1,
        "star_sector_2": sector2,
    }

    print(f"\n  {loc_name}:")
    print(f"    Distance to circle: {dist_to_circle:.1f} km")
    print(f"    GC bearing: {b1:.1f}° / {b2:.1f}°")
    print(f"    Star compass: {sector1} / {sector2}")


# ═══════════════════════════════════════════════════════════════════════
# Save outputs
# ═══════════════════════════════════════════════════════════════════════

# Route geometries CSV
with open("route_geometries.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "name", "n_waypoints", "route_type", "total_distance_km",
        "fit_pole_lat", "fit_pole_lon", "rms_from_gc_km", "rms_from_rhumb_km",
        "gc_preferred", "angular_sep_from_great_circle_deg",
    ])
    writer.writeheader()
    for r in route_results:
        writer.writerow({k: r[k] for k in writer.fieldnames})

print(f"\nSaved route_geometries.csv")

# Route-circle comparison JSON
comparison = {
    "analysis": "Voyaging route great circle comparison",
    "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
    "n_routes": len(route_results),
    "routes": route_results,
    "routes_within_5deg": [r["name"] for r in route_results if r["within_5deg"]],
    "routes_within_10deg": [r["name"] for r in route_results if r["within_10deg"]],
    "summary": {
        "mean_angular_sep": round(np.mean([r["angular_sep_from_great_circle_deg"] for r in route_results]), 2),
        "min_angular_sep": round(min(r["angular_sep_from_great_circle_deg"] for r in route_results), 2),
        "min_angular_sep_route": min(route_results, key=lambda r: r["angular_sep_from_great_circle_deg"])["name"],
        "routes_gc_preferred_count": sum(1 for r in route_results if r["gc_preferred"]),
        "routes_rhumb_preferred_count": sum(1 for r in route_results if r["gc_preferred"] is False),
    },
}

with open("route_circle_comparison.json", "w") as f:
    json.dump(comparison, f, indent=2)

print("Saved route_circle_comparison.json")

# Star compass JSON
with open("star_compass_bearing.json", "w") as f:
    json.dump(star_compass_results, f, indent=2)

print("Saved star_compass_bearing.json")
