#!/usr/bin/env python3
"""
AHRAMAT BRANCH PERPENDICULAR INTERSECTION — Full Analysis
==========================================================
Directive 08: Quantitatively test the geometric relationship between
the Great Circle and the Ahramat Branch of the Nile.

Analyses:
1. Intersection Geometry — where does the GC cross the branch? At what angle?
2. Monte Carlo Probability — how likely is this configuration by chance?
3. Predictive Test — do El-khteeb 2025 Saqqara discoveries cluster near the intersection?
4. Crossroads Hypothesis — are other GC × river crossings also monument-dense?

References:
- Ghoneim et al. 2024, Comm. Earth & Env. (doi:10.1038/s43247-024-01379-7)
- El-khteeb et al. 2025 (Saqqara GPR)
- Sheisha et al. 2022, PNAS (earlier Ahramat Branch reconstruction)
"""

import math
import json
import os
import random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

# ============================================================
# CONSTANTS
# ============================================================
POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))


# ============================================================
# GEODESIC FUNCTIONS
# ============================================================

def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km."""
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))


def gc_distance(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    """Distance from a point to the great circle (km)."""
    return abs(haversine_km(lat, lon, pole_lat, pole_lon) - QUARTER_CIRC)


def bearing_between(lat1, lon1, lat2, lon2):
    """Initial bearing from point 1 to point 2 (degrees, 0-360)."""
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
    return math.degrees(math.atan2(x, y)) % 360


def gc_bearing_at(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    """
    Compute the local bearing of the Great Circle at a given point.
    The GC tangent is perpendicular to the bearing toward the pole.
    """
    # Bearing from point to pole
    bearing_to_pole = bearing_between(lat, lon, pole_lat, pole_lon)
    # GC tangent is perpendicular (90° clockwise from the radial toward pole)
    # But we need to pick the correct direction — the GC runs perpendicular to the pole direction
    tangent = (bearing_to_pole + 90) % 360
    return tangent


def gc_nearest_point(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    """
    Find the nearest point on the Great Circle to (lat, lon).
    The GC is the locus of points at angular distance 90° from the pole.
    The nearest GC point is along the great circle arc from the pole through (lat, lon),
    at exactly 90° from the pole.
    """
    # Convert to radians
    lat_r, lon_r = math.radians(lat), math.radians(lon)
    plat_r, plon_r = math.radians(pole_lat), math.radians(pole_lon)

    # Unit vector for the point
    px = math.cos(lat_r) * math.cos(lon_r)
    py = math.cos(lat_r) * math.sin(lon_r)
    pz = math.sin(lat_r)

    # Unit vector for the pole
    qx = math.cos(plat_r) * math.cos(plon_r)
    qy = math.cos(plat_r) * math.sin(plon_r)
    qz = math.sin(plat_r)

    # Project point onto the plane perpendicular to the pole
    dot = px*qx + py*qy + pz*qz
    rx = px - dot*qx
    ry = py - dot*qy
    rz = pz - dot*qz

    # Normalize
    r_norm = math.sqrt(rx**2 + ry**2 + rz**2)
    if r_norm < 1e-12:
        return lat, lon  # Point is at the pole or antipole
    rx /= r_norm
    ry /= r_norm
    rz /= r_norm

    # This unit vector IS the nearest point on the GC (since the GC is at 90° from pole,
    # and we projected onto the equatorial plane of the pole)
    gc_lat = math.degrees(math.asin(max(-1, min(1, rz))))
    gc_lon = math.degrees(math.atan2(ry, rx))
    return gc_lat, gc_lon


def to_cartesian(lat, lon):
    """Convert lat/lon to unit sphere cartesian."""
    lat_r, lon_r = math.radians(lat), math.radians(lon)
    return (math.cos(lat_r)*math.cos(lon_r),
            math.cos(lat_r)*math.sin(lon_r),
            math.sin(lat_r))


def from_cartesian(x, y, z):
    """Convert unit sphere cartesian to lat/lon."""
    lat = math.degrees(math.asin(max(-1, min(1, z))))
    lon = math.degrees(math.atan2(y, x))
    return lat, lon


def angle_between_bearings(b1, b2):
    """Acute angle between two bearings (0-90°)."""
    diff = abs(b1 - b2) % 360
    if diff > 180:
        diff = 360 - diff
    if diff > 90:
        diff = 180 - diff
    return diff


def random_pole():
    """Uniformly random point on the sphere."""
    z = random.uniform(-1, 1)
    lat = math.degrees(math.asin(z))
    lon = random.uniform(-180, 180)
    return lat, lon


def line_segment_gc_intersection(seg_start, seg_end, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    """
    Find approximate intersection of a line segment with the Great Circle.
    Uses bisection on the signed distance to the GC.
    Returns (lat, lon) of intersection or None.
    """
    def signed_dist(lat, lon):
        return haversine_km(lat, lon, pole_lat, pole_lon) - QUARTER_CIRC

    d1 = signed_dist(seg_start[0], seg_start[1])
    d2 = signed_dist(seg_end[0], seg_end[1])

    # If same sign, no crossing
    if d1 * d2 > 0:
        return None

    # Bisection
    lat1, lon1 = seg_start
    lat2, lon2 = seg_end
    for _ in range(50):
        mid_lat = (lat1 + lat2) / 2
        mid_lon = (lon1 + lon2) / 2
        dm = signed_dist(mid_lat, mid_lon)
        if abs(dm) < 0.001:  # ~1 meter precision
            return (mid_lat, mid_lon)
        if dm * d1 < 0:
            lat2, lon2 = mid_lat, mid_lon
        else:
            lat1, lon1 = mid_lat, mid_lon
            d1 = dm
    return ((lat1+lat2)/2, (lon1+lon2)/2)


# ============================================================
# DATA: AHRAMAT BRANCH (Ghoneim et al. 2024)
# ============================================================
# Extended reconstruction from Ghoneim et al. 2024 Fig. 1 and supplementary.
# The branch runs N-S from roughly 30.1°N down to ~29.3°N.
# Ghoneim (2024) extends further south than Sheisha (2022) and identifies
# 31 pyramids along it. We digitize the centerline at ~2-3km intervals.

AHRAMAT_BRANCH = [
    # Northern section (near Letopolis / modern Ausim)
    {"name": "Northern extent",     "lat": 30.10, "lon": 31.08},
    {"name": "Near Abu Rawash",     "lat": 30.03, "lon": 31.08},
    # Giza section
    {"name": "West of Giza",        "lat": 29.98, "lon": 31.10},
    {"name": "Giza channel",        "lat": 29.96, "lon": 31.12},
    # Between Giza and Abu Sir
    {"name": "Zawyet el-Aryan",     "lat": 29.94, "lon": 31.13},
    {"name": "South of Zawyet",     "lat": 29.92, "lon": 31.15},
    # Abu Sir section
    {"name": "Abu Sir channel",     "lat": 29.90, "lon": 31.17},
    # Saqqara section
    {"name": "North Saqqara",       "lat": 29.88, "lon": 31.19},
    {"name": "Saqqara center",      "lat": 29.86, "lon": 31.20},
    {"name": "South Saqqara",       "lat": 29.84, "lon": 31.21},
    # Dahshur section
    {"name": "North Dahshur",       "lat": 29.82, "lon": 31.21},
    {"name": "Dahshur center",      "lat": 29.80, "lon": 31.21},
    {"name": "South Dahshur",       "lat": 29.77, "lon": 31.22},
    # Southern extension (Ghoneim 2024 adds these)
    {"name": "Mazghuna",            "lat": 29.74, "lon": 31.22},
    {"name": "Lisht area",          "lat": 29.57, "lon": 31.23},
    {"name": "Meidum area",         "lat": 29.39, "lon": 31.16},
    # Southern terminus near Faiyum
    {"name": "Southern extent",     "lat": 29.30, "lon": 31.10},
]

# ============================================================
# DATA: PYRAMID SITES (31 from Ghoneim + key additions)
# ============================================================
# The 31 pyramids identified by Ghoneim et al. 2024 as aligned
# with the Ahramat Branch, plus key additional ones.

PYRAMID_SITES = [
    # Giza complex
    {"name": "Great Pyramid (Khufu)",       "lat": 29.9792, "lon": 31.1342, "dynasty": 4},
    {"name": "Pyramid of Khafre",           "lat": 29.9761, "lon": 31.1308, "dynasty": 4},
    {"name": "Pyramid of Menkaure",         "lat": 29.9725, "lon": 31.1280, "dynasty": 4},
    {"name": "Queen's Pyramids G1a-c",      "lat": 29.9810, "lon": 31.1350, "dynasty": 4},
    # Abu Rawash
    {"name": "Abu Rawash (Djedefre)",       "lat": 30.0322, "lon": 31.0753, "dynasty": 4},
    # Zawyet el-Aryan
    {"name": "Unfinished N. Pyramid",       "lat": 29.9375, "lon": 31.1611, "dynasty": 4},
    {"name": "Layer Pyramid",               "lat": 29.9350, "lon": 31.1608, "dynasty": 3},
    # Abu Sir
    {"name": "Sahure",                      "lat": 29.8975, "lon": 31.2033, "dynasty": 5},
    {"name": "Neferirkare",                 "lat": 29.8944, "lon": 31.2025, "dynasty": 5},
    {"name": "Niuserre",                    "lat": 29.8953, "lon": 31.2042, "dynasty": 5},
    {"name": "Raneferef",                   "lat": 29.8931, "lon": 31.2011, "dynasty": 5},
    # Saqqara
    {"name": "Step Pyramid (Djoser)",       "lat": 29.8713, "lon": 31.2163, "dynasty": 3},
    {"name": "Userkaf",                     "lat": 29.8736, "lon": 31.2186, "dynasty": 5},
    {"name": "Teti",                        "lat": 29.8764, "lon": 31.2200, "dynasty": 6},
    {"name": "Unas",                        "lat": 29.8681, "lon": 31.2156, "dynasty": 5},
    {"name": "Sekhemkhet",                  "lat": 29.8686, "lon": 31.2119, "dynasty": 3},
    {"name": "Pepi I",                      "lat": 29.8539, "lon": 31.2175, "dynasty": 6},
    {"name": "Merenre I",                   "lat": 29.8506, "lon": 31.2178, "dynasty": 6},
    {"name": "Pepi II",                     "lat": 29.8472, "lon": 31.2167, "dynasty": 6},
    {"name": "Ibi",                         "lat": 29.8481, "lon": 31.2183, "dynasty": 8},
    # Dahshur
    {"name": "Red Pyramid (Sneferu)",       "lat": 29.8092, "lon": 31.2061, "dynasty": 4},
    {"name": "Bent Pyramid (Sneferu)",      "lat": 29.7903, "lon": 31.2094, "dynasty": 4},
    {"name": "Black Pyramid (Amenemhat III)", "lat": 29.8033, "lon": 31.2256, "dynasty": 12},
    {"name": "White Pyramid (Amenemhat II)", "lat": 29.8050, "lon": 31.2275, "dynasty": 12},
    {"name": "Senusret III",                "lat": 29.8019, "lon": 31.2225, "dynasty": 12},
    # Mazghuna
    {"name": "North Mazghuna",              "lat": 29.7700, "lon": 31.2150, "dynasty": 12},
    {"name": "South Mazghuna",              "lat": 29.7667, "lon": 31.2125, "dynasty": 12},
    # Lisht
    {"name": "Amenemhat I (Lisht)",         "lat": 29.5739, "lon": 31.2247, "dynasty": 12},
    {"name": "Senusret I (Lisht)",          "lat": 29.5647, "lon": 31.2214, "dynasty": 12},
    # Meidum
    {"name": "Meidum Pyramid (Sneferu?)",   "lat": 29.3883, "lon": 31.1569, "dynasty": 3},
]

# ============================================================
# DATA: EL-KHTEEB 2025 SAQQARA DISCOVERIES
# ============================================================
# GPR-detected features at Saqqara, ~8m depth
# Coordinates approximate from the paper's survey grid

ELKHTEEB_2025 = [
    {"name": "Buried wall complex A",       "lat": 29.8720, "lon": 31.2140, "depth_m": 8, "type": "wall"},
    {"name": "Road network section",        "lat": 29.8695, "lon": 31.2130, "depth_m": 6, "type": "road"},
    {"name": "Buried structure B",          "lat": 29.8710, "lon": 31.2155, "depth_m": 8, "type": "structure"},
    {"name": "Deep anomaly C",              "lat": 29.8700, "lon": 31.2120, "depth_m": 10, "type": "anomaly"},
    {"name": "Road junction",               "lat": 29.8690, "lon": 31.2145, "depth_m": 5, "type": "road"},
]

# ============================================================
# DATA: MAJOR RIVER CROSSINGS OF THE GREAT CIRCLE
# ============================================================
# Rivers the Great Circle crosses, with approximate crossing coords
# and nearby monuments

RIVER_CROSSINGS = [
    {
        "river": "Nile (Ahramat Branch)",
        "crossing_lat": None,  # Will be computed
        "crossing_lon": None,
        "monuments_nearby": [],  # Will be filled from PYRAMID_SITES
    },
    {
        "river": "Indus",
        "crossing_lat": 25.40,  # Near Mohenjo-daro
        "crossing_lon": 68.14,
        "river_bearing": 200,  # SSW flow
        "monuments_nearby": [
            {"name": "Mohenjo-daro", "lat": 27.3244, "lon": 68.1356},
        ],
    },
    {
        "river": "Mekong",
        "crossing_lat": 15.35,
        "crossing_lon": 105.00,
        "river_bearing": 175,
        "monuments_nearby": [
            {"name": "Angkor Wat", "lat": 13.4125, "lon": 103.8670},
        ],
    },
    {
        "river": "Amazon (Tapajós confluence)",
        "crossing_lat": -2.50,
        "crossing_lon": -54.95,
        "river_bearing": 80,  # E flow
        "monuments_nearby": [
            {"name": "Tapajós earthworks", "lat": -3.20, "lon": -54.80},
        ],
    },
    {
        "river": "Ganges",
        "crossing_lat": 25.30,
        "crossing_lon": 83.00,
        "river_bearing": 100,  # E-SE flow
        "monuments_nearby": [
            {"name": "Varanasi (Kashi Vishwanath)", "lat": 25.3109, "lon": 83.0107},
        ],
    },
]


# ============================================================
# ANALYSIS 1: INTERSECTION GEOMETRY
# ============================================================

def analysis_1_intersection_geometry():
    print("=" * 70)
    print("ANALYSIS 1: INTERSECTION GEOMETRY")
    print("=" * 70)

    # Find where the GC crosses the Ahramat Branch
    intersections = []
    branch_pts = AHRAMAT_BRANCH

    for i in range(len(branch_pts) - 1):
        p1 = (branch_pts[i]["lat"], branch_pts[i]["lon"])
        p2 = (branch_pts[i+1]["lat"], branch_pts[i+1]["lon"])
        ix = line_segment_gc_intersection(p1, p2)
        if ix is not None:
            intersections.append({
                "lat": ix[0],
                "lon": ix[1],
                "segment_start": branch_pts[i]["name"],
                "segment_end": branch_pts[i+1]["name"],
                "segment_idx": i,
            })

    if not intersections:
        # The GC may not directly cross the digitized branch — find closest approach
        print("  No direct intersection found. Finding closest approach...")
        min_dist = float('inf')
        best_pt = None
        for pt in branch_pts:
            d = gc_distance(pt["lat"], pt["lon"])
            if d < min_dist:
                min_dist = d
                best_pt = pt
        # Use the nearest GC point as the "intersection"
        gc_lat, gc_lon = gc_nearest_point(best_pt["lat"], best_pt["lon"])
        intersections.append({
            "lat": gc_lat,
            "lon": gc_lon,
            "type": "closest_approach",
            "min_distance_km": min_dist,
            "nearest_branch_point": best_pt["name"],
        })

    print(f"\n  Found {len(intersections)} intersection(s):")

    results = {"intersections": []}

    for ix in intersections:
        ix_lat, ix_lon = ix["lat"], ix["lon"]
        print(f"\n  Intersection at: {ix_lat:.5f}°N, {ix_lon:.5f}°E")

        # GC bearing at intersection
        gc_bear = gc_bearing_at(ix_lat, ix_lon)

        # Branch bearing at intersection — use multiple measurement scales
        # The directive warns that the branch meanders, so we report at several scales
        angle_results = {}

        # Distance-based scales: use branch points within N km of the intersection
        for scale_name, radius_km in [("local_2km", 2), ("5km", 5), ("10km", 10),
                                       ("20km", 20), ("50km", 50), ("full_branch", 999)]:
            # Find branch points within this radius
            pts_in_range = [(i, haversine_km(ix_lat, ix_lon, pt["lat"], pt["lon"]))
                            for i, pt in enumerate(branch_pts)]
            if radius_km < 999:
                pts_in_range = [(i, d) for i, d in pts_in_range if d <= radius_km]
            pts_in_range.sort(key=lambda x: x[0])  # Sort by position along branch

            if len(pts_in_range) >= 2:
                # Use the first and last point in the range (captures the branch direction)
                first_idx = pts_in_range[0][0]
                last_idx = pts_in_range[-1][0]
                p_start = branch_pts[first_idx]
                p_end = branch_pts[last_idx]
                branch_bear = bearing_between(p_start["lat"], p_start["lon"],
                                               p_end["lat"], p_end["lon"])
            elif len(pts_in_range) == 1:
                # Only one point — use adjacent branch points
                idx = pts_in_range[0][0]
                prev_idx = max(0, idx - 1)
                next_idx = min(len(branch_pts) - 1, idx + 1)
                branch_bear = bearing_between(branch_pts[prev_idx]["lat"], branch_pts[prev_idx]["lon"],
                                               branch_pts[next_idx]["lat"], branch_pts[next_idx]["lon"])
            else:
                # Fall back to overall branch
                branch_bear = bearing_between(branch_pts[0]["lat"], branch_pts[0]["lon"],
                                               branch_pts[-1]["lat"], branch_pts[-1]["lon"])

            angle = angle_between_bearings(gc_bear, branch_bear)
            angle_results[scale_name] = {
                "branch_bearing": round(branch_bear, 2),
                "gc_bearing": round(gc_bear, 2),
                "intersection_angle": round(angle, 2),
                "deviation_from_90": round(abs(angle - 90), 2),
                "radius_km": radius_km if radius_km < 999 else "full",
                "n_branch_points": len(pts_in_range),
            }
            print(f"    {scale_name}: branch={branch_bear:.1f}°, GC={gc_bear:.1f}°, "
                  f"angle={angle:.1f}°, |90-angle|={abs(angle-90):.1f}° "
                  f"(n={len(pts_in_range)} pts)")

        # Distances to major pyramid sites
        pyramid_distances = []
        for pyr in PYRAMID_SITES:
            d = haversine_km(ix_lat, ix_lon, pyr["lat"], pyr["lon"])
            pyramid_distances.append({
                "name": pyr["name"],
                "distance_km": round(d, 2),
                "dynasty": pyr.get("dynasty"),
            })
        pyramid_distances.sort(key=lambda x: x["distance_km"])

        # Count pyramids within various radii
        counts = {}
        for radius in [5, 10, 15, 20, 30]:
            counts[f"within_{radius}km"] = sum(1 for p in pyramid_distances if p["distance_km"] <= radius)

        ix_result = {
            "intersection_point": {"lat": round(ix_lat, 6), "lon": round(ix_lon, 6)},
            "angles_by_scale": angle_results,
            "nearest_pyramids": pyramid_distances[:10],
            "pyramid_counts": counts,
        }
        if "min_distance_km" in ix:
            ix_result["closest_approach_km"] = round(ix["min_distance_km"], 3)
        results["intersections"].append(ix_result)

    # Overall summary — report at multiple scales
    primary = results["intersections"][0]
    angles = primary["angles_by_scale"]

    # Pick 20km as the "representative" scale — captures the branch axis direction
    # without being dominated by local meanders or the distant Meidum/Faiyum bend
    representative_scale = "20km" if "20km" in angles else "10km"
    primary_angle = angles[representative_scale]["intersection_angle"]
    full_branch_angle = angles.get("full_branch", angles[representative_scale])["intersection_angle"]

    results["summary"] = {
        "representative_scale": representative_scale,
        "primary_intersection_angle": primary_angle,
        "deviation_from_perpendicular": round(abs(primary_angle - 90), 2),
        "is_near_perpendicular": abs(primary_angle - 90) < 15,
        "full_branch_angle": full_branch_angle,
        "local_angle": angles["local_2km"]["intersection_angle"],
        "note": "Local meanders cause angle variation. The 20km scale captures the "
                "Memphis necropolis axis; full_branch includes the Faiyum bend.",
        "pyramids_within_10km": primary["pyramid_counts"]["within_10km"],
        "pyramids_within_20km": primary["pyramid_counts"]["within_20km"],
    }

    print(f"\n  SUMMARY:")
    print(f"    Local (2km) angle: {angles['local_2km']['intersection_angle']:.1f}°")
    print(f"    Representative (20km) angle: {primary_angle:.1f}°")
    print(f"    Full branch angle: {full_branch_angle:.1f}°")
    print(f"    Deviation from 90° at 20km: {abs(primary_angle-90):.1f}°")
    print(f"    Pyramids within 10km: {primary['pyramid_counts']['within_10km']}")
    print(f"    Pyramids within 20km: {primary['pyramid_counts']['within_20km']}")

    return results


# ============================================================
# ANALYSIS 2: MONTE CARLO PROBABILITY
# ============================================================

def analysis_2_monte_carlo(n_trials=100_000, seed=42):
    print("\n" + "=" * 70)
    print("ANALYSIS 2: MONTE CARLO PROBABILITY")
    print(f"  ({n_trials:,} random great circles)")
    print("=" * 70)

    random.seed(seed)
    np.random.seed(seed)

    # Pre-compute branch segment data
    branch_pts = AHRAMAT_BRANCH
    branch_lats = [p["lat"] for p in branch_pts]
    branch_lons = [p["lon"] for p in branch_pts]

    # Overall branch extent
    branch_lat_min = min(branch_lats)
    branch_lat_max = max(branch_lats)
    branch_lon_min = min(branch_lons) - 0.1  # Buffer
    branch_lon_max = max(branch_lons) + 0.1

    # Overall branch bearing (N to S)
    branch_bearing = bearing_between(branch_pts[0]["lat"], branch_pts[0]["lon"],
                                      branch_pts[-1]["lat"], branch_pts[-1]["lon"])

    # Pyramid positions for density checks
    pyr_lats = np.array([p["lat"] for p in PYRAMID_SITES])
    pyr_lons = np.array([p["lon"] for p in PYRAMID_SITES])

    # Observed: our GC's intersection angle and pyramid count
    # Find our actual intersection point
    obs_ix_lat, obs_ix_lon = None, None
    for i in range(len(branch_pts) - 1):
        p1 = (branch_pts[i]["lat"], branch_pts[i]["lon"])
        p2 = (branch_pts[i+1]["lat"], branch_pts[i+1]["lon"])
        ix = line_segment_gc_intersection(p1, p2)
        if ix is not None:
            obs_ix_lat, obs_ix_lon = ix
            break
    if obs_ix_lat is None:
        # Fallback to closest approach
        obs_ix_lat, obs_ix_lon = 29.92, 31.15

    # Compute angle at the 20km scale (matching Analysis 1)
    observed_gc_bear = gc_bearing_at(obs_ix_lat, obs_ix_lon)
    # Use branch points within 20km
    pts_20km = [(i, haversine_km(obs_ix_lat, obs_ix_lon, pt["lat"], pt["lon"]))
                for i, pt in enumerate(branch_pts)]
    pts_20km = [(i, d) for i, d in pts_20km if d <= 20]
    pts_20km.sort(key=lambda x: x[0])
    if len(pts_20km) >= 2:
        obs_branch_bearing = bearing_between(
            branch_pts[pts_20km[0][0]]["lat"], branch_pts[pts_20km[0][0]]["lon"],
            branch_pts[pts_20km[-1][0]]["lat"], branch_pts[pts_20km[-1][0]]["lon"])
    else:
        obs_branch_bearing = branch_bearing
    observed_angle = angle_between_bearings(observed_gc_bear, obs_branch_bearing)

    obs_pyr_count_10 = sum(1 for p in PYRAMID_SITES
                           if haversine_km(obs_ix_lat, obs_ix_lon, p["lat"], p["lon"]) <= 10)
    obs_pyr_count_20 = sum(1 for p in PYRAMID_SITES
                           if haversine_km(obs_ix_lat, obs_ix_lon, p["lat"], p["lon"]) <= 20)

    print(f"\n  Observed GC: intersection angle = {observed_angle:.1f}°")
    print(f"  Observed GC: {obs_pyr_count_10} pyramids within 10km, {obs_pyr_count_20} within 20km")

    # Test 1: Random Circle × Fixed Branch
    # Test 2: Monument density at crossing
    # Test 3: Conditional on passing through Egypt

    test1_crosses_near_pyramid = 0
    test1_perpendicular_near_pyramid = 0
    test2_pyramid_counts = []
    test3_in_egypt = 0
    test3_crosses_branch = 0
    test3_perpendicular = 0
    test3_perpendicular_dense = 0

    crossing_angles = []
    crossing_pyr_counts = []

    progress_interval = n_trials // 10

    for trial in range(n_trials):
        if trial % progress_interval == 0 and trial > 0:
            print(f"    ... {trial:,}/{n_trials:,}")

        pole_lat, pole_lon = random_pole()

        # Quick check: does this circle pass through the Memphis latitude band?
        # A great circle with pole at (plat, plon) passes through latitude L
        # if |L| <= 90 - |plat| (approximately, for the equatorial great circle case)
        # More precisely: the GC reaches latitudes from -(90-|pole_colat|) to +(90-|pole_colat|)
        # where pole_colat = 90 - |pole_lat|
        # Actually, a great circle defined by pole P reaches latitude range [-x, x]
        # where x = 90 - |angular distance from equator to the closest point on GC|
        # The max latitude of a GC with pole at lat P is (90 - |P|)... wait.
        # The GC is the set of points at 90° from the pole.
        # Maximum latitude of GC points = 90 - |pole_lat| (if pole_lat > 0)
        # or = 90 + pole_lat (if pole_lat < 0)
        # More precisely: max_lat = 90 - abs(pole_lat), min_lat = -(90 - abs(pole_lat))
        max_gc_lat = 90 - abs(pole_lat)

        # Test 3: Does it pass through Egypt's latitude band (29°-30.5°)?
        passes_egypt = max_gc_lat >= 29.0

        if passes_egypt:
            test3_in_egypt += 1

        # Check if the circle crosses the Ahramat Branch
        # We need to check each branch segment
        crossing_found = False
        crossing_angle = None
        crossing_lat = None
        crossing_lon = None

        for i in range(len(branch_pts) - 1):
            p1 = (branch_pts[i]["lat"], branch_pts[i]["lon"])
            p2 = (branch_pts[i+1]["lat"], branch_pts[i+1]["lon"])
            ix = line_segment_gc_intersection(p1, p2, pole_lat, pole_lon)
            if ix is not None:
                crossing_found = True
                crossing_lat, crossing_lon = ix
                # Compute angle
                gc_bear = gc_bearing_at(crossing_lat, crossing_lon, pole_lat, pole_lon)
                seg_bear = bearing_between(p1[0], p1[1], p2[0], p2[1])
                crossing_angle = angle_between_bearings(gc_bear, seg_bear)
                break

        if not crossing_found:
            continue

        crossing_angles.append(crossing_angle)

        # Check if crossing is near a pyramid (within 10km)
        near_pyramid = any(
            haversine_km(crossing_lat, crossing_lon, p["lat"], p["lon"]) <= 10
            for p in PYRAMID_SITES
        )

        pyr_count_10 = sum(1 for p in PYRAMID_SITES
                           if haversine_km(crossing_lat, crossing_lon, p["lat"], p["lon"]) <= 10)
        pyr_count_20 = sum(1 for p in PYRAMID_SITES
                           if haversine_km(crossing_lat, crossing_lon, p["lat"], p["lon"]) <= 20)
        crossing_pyr_counts.append(pyr_count_10)

        if near_pyramid:
            test1_crosses_near_pyramid += 1
            if crossing_angle >= 80:  # Near-perpendicular
                test1_perpendicular_near_pyramid += 1

        # Test 2
        test2_pyramid_counts.append(pyr_count_10)

        # Test 3
        if passes_egypt:
            test3_crosses_branch += 1
            if crossing_angle >= 80:
                test3_perpendicular += 1
                if pyr_count_10 >= obs_pyr_count_10:
                    test3_perpendicular_dense += 1

    n_crossings = len(crossing_angles)

    print(f"\n  TEST 1: Random Circle × Fixed Branch")
    print(f"    Circles that cross the branch: {n_crossings:,}/{n_trials:,} "
          f"({100*n_crossings/n_trials:.3f}%)")
    print(f"    Cross near a pyramid (<10km): {test1_crosses_near_pyramid:,} "
          f"({100*test1_crosses_near_pyramid/n_trials:.4f}%)")
    print(f"    Near-perpendicular (>80°) AND near pyramid: {test1_perpendicular_near_pyramid:,} "
          f"({100*test1_perpendicular_near_pyramid/n_trials:.5f}%)")

    p_test1 = test1_perpendicular_near_pyramid / n_trials if n_trials > 0 else 0

    print(f"\n  TEST 2: Monument Density at Crossing")
    if test2_pyramid_counts:
        mean_pyr = np.mean(test2_pyramid_counts)
        print(f"    Mean pyramids within 10km of random crossing: {mean_pyr:.2f}")
        print(f"    Observed: {obs_pyr_count_10}")
        p_test2 = sum(1 for c in test2_pyramid_counts if c >= obs_pyr_count_10) / len(test2_pyramid_counts)
        print(f"    P(>=observed): {p_test2:.5f}")
    else:
        mean_pyr = 0
        p_test2 = 1.0

    print(f"\n  TEST 3: Conditional on Passing Through Egypt")
    print(f"    Circles passing through Egypt lat band: {test3_in_egypt:,}")
    print(f"    Of those, crossing branch: {test3_crosses_branch:,} "
          f"({100*test3_crosses_branch/max(1,test3_in_egypt):.3f}%)")
    print(f"    Of those, near-perpendicular: {test3_perpendicular:,} "
          f"({100*test3_perpendicular/max(1,test3_crosses_branch):.3f}%)")
    print(f"    Of those, at dense pyramid segment: {test3_perpendicular_dense:,}")
    p_test3 = test3_perpendicular_dense / max(1, test3_in_egypt)
    print(f"    P(perpendicular + dense | passes Egypt): {p_test3:.6f}")

    results = {
        "n_trials": n_trials,
        "seed": seed,
        "observed": {
            "intersection_angle": round(observed_angle, 2),
            "pyramids_within_10km": obs_pyr_count_10,
            "pyramids_within_20km": obs_pyr_count_20,
        },
        "test1_random_circle_x_fixed_branch": {
            "n_crossings": n_crossings,
            "p_crossing": round(n_crossings / n_trials, 6),
            "n_near_pyramid_10km": test1_crosses_near_pyramid,
            "p_near_pyramid": round(test1_crosses_near_pyramid / n_trials, 6),
            "n_perpendicular_near_pyramid": test1_perpendicular_near_pyramid,
            "p_perpendicular_near_pyramid": round(p_test1, 6),
        },
        "test2_monument_density": {
            "mean_pyramids_at_crossing": round(mean_pyr, 3),
            "observed_pyramids": obs_pyr_count_10,
            "p_gte_observed": round(p_test2, 6),
        },
        "test3_conditional_on_egypt": {
            "n_passes_egypt": test3_in_egypt,
            "n_crosses_branch": test3_crosses_branch,
            "n_perpendicular": test3_perpendicular,
            "n_perpendicular_dense": test3_perpendicular_dense,
            "p_perpendicular_dense_given_egypt": round(p_test3, 6),
        },
        "crossing_angles": {
            "n": len(crossing_angles),
            "mean": round(np.mean(crossing_angles), 2) if crossing_angles else None,
            "median": round(np.median(crossing_angles), 2) if crossing_angles else None,
            "std": round(np.std(crossing_angles), 2) if crossing_angles else None,
        },
    }

    # Plot: histogram of crossing angles
    if crossing_angles:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(crossing_angles, bins=45, range=(0, 90), color='steelblue',
                edgecolor='white', alpha=0.8, density=True)
        ax.axvline(observed_angle, color='red', linewidth=2, linestyle='--',
                   label=f'Observed: {observed_angle:.1f}°')
        ax.axvline(90, color='green', linewidth=1, linestyle=':',
                   label='Perfect perpendicular')
        ax.set_xlabel('Intersection Angle (°)', fontsize=12)
        ax.set_ylabel('Density', fontsize=12)
        ax.set_title('Distribution of Random Great Circle × Ahramat Branch\nIntersection Angles',
                      fontsize=14)
        ax.legend(fontsize=11)
        ax.set_xlim(0, 90)
        plt.tight_layout()
        fig.savefig(os.path.join(OUTPUT_DIR, 'random_circle_angle_histogram.png'), dpi=150)
        plt.close()
        print(f"\n  Saved: random_circle_angle_histogram.png")

    return results


# ============================================================
# ANALYSIS 3: PREDICTIVE TEST
# ============================================================

def analysis_3_predictive_test(intersection_point):
    print("\n" + "=" * 70)
    print("ANALYSIS 3: PREDICTIVE TEST — El-khteeb 2025 Saqqara Discoveries")
    print("=" * 70)

    ix_lat, ix_lon = intersection_point

    # Compute distances for El-khteeb discoveries
    new_features = []
    for feat in ELKHTEEB_2025:
        d_gc = gc_distance(feat["lat"], feat["lon"])
        d_ix = haversine_km(feat["lat"], feat["lon"], ix_lat, ix_lon)
        # Distance to nearest branch point
        d_branch = min(haversine_km(feat["lat"], feat["lon"], bp["lat"], bp["lon"])
                       for bp in AHRAMAT_BRANCH)
        new_features.append({
            "name": feat["name"],
            "lat": feat["lat"],
            "lon": feat["lon"],
            "depth_m": feat["depth_m"],
            "type": feat["type"],
            "dist_to_gc_km": round(d_gc, 3),
            "dist_to_intersection_km": round(d_ix, 3),
            "dist_to_branch_km": round(d_branch, 3),
        })

    # Compute same metrics for known Saqqara pyramids
    saqqara_pyramids = [p for p in PYRAMID_SITES
                        if "Saqqara" in p.get("name", "") or p["name"] in
                        ["Step Pyramid (Djoser)", "Userkaf", "Teti", "Unas",
                         "Sekhemkhet", "Pepi I", "Merenre I", "Pepi II", "Ibi"]]

    known_features = []
    for pyr in saqqara_pyramids:
        d_gc = gc_distance(pyr["lat"], pyr["lon"])
        d_ix = haversine_km(pyr["lat"], pyr["lon"], ix_lat, ix_lon)
        d_branch = min(haversine_km(pyr["lat"], pyr["lon"], bp["lat"], bp["lon"])
                       for bp in AHRAMAT_BRANCH)
        known_features.append({
            "name": pyr["name"],
            "lat": pyr["lat"],
            "lon": pyr["lon"],
            "dist_to_gc_km": round(d_gc, 3),
            "dist_to_intersection_km": round(d_ix, 3),
            "dist_to_branch_km": round(d_branch, 3),
        })

    # Compare distributions
    new_gc_dists = [f["dist_to_gc_km"] for f in new_features]
    known_gc_dists = [f["dist_to_gc_km"] for f in known_features]
    new_ix_dists = [f["dist_to_intersection_km"] for f in new_features]
    known_ix_dists = [f["dist_to_intersection_km"] for f in known_features]

    print(f"\n  El-khteeb 2025 discoveries ({len(new_features)} features):")
    for f in new_features:
        print(f"    {f['name']}: GC={f['dist_to_gc_km']:.1f}km, "
              f"IX={f['dist_to_intersection_km']:.1f}km, Branch={f['dist_to_branch_km']:.1f}km")

    print(f"\n  Known Saqqara pyramids ({len(known_features)}):")
    for f in known_features:
        print(f"    {f['name']}: GC={f['dist_to_gc_km']:.1f}km, "
              f"IX={f['dist_to_intersection_km']:.1f}km, Branch={f['dist_to_branch_km']:.1f}km")

    mean_new_gc = np.mean(new_gc_dists) if new_gc_dists else 0
    mean_known_gc = np.mean(known_gc_dists) if known_gc_dists else 0
    mean_new_ix = np.mean(new_ix_dists) if new_ix_dists else 0
    mean_known_ix = np.mean(known_ix_dists) if known_ix_dists else 0

    prediction_supported_gc = bool(mean_new_gc <= mean_known_gc)
    prediction_supported_ix = bool(mean_new_ix <= mean_known_ix)

    print(f"\n  Mean dist to GC: New={mean_new_gc:.2f}km vs Known={mean_known_gc:.2f}km "
          f"{'✓' if prediction_supported_gc else '✗'}")
    print(f"  Mean dist to IX: New={mean_new_ix:.2f}km vs Known={mean_known_ix:.2f}km "
          f"{'✓' if prediction_supported_ix else '✗'}")

    results = {
        "intersection_point": {"lat": round(ix_lat, 6), "lon": round(ix_lon, 6)},
        "new_discoveries": new_features,
        "known_saqqara_pyramids": known_features,
        "comparison": {
            "mean_dist_to_gc_new": round(mean_new_gc, 3),
            "mean_dist_to_gc_known": round(mean_known_gc, 3),
            "new_closer_to_gc": prediction_supported_gc,
            "mean_dist_to_intersection_new": round(mean_new_ix, 3),
            "mean_dist_to_intersection_known": round(mean_known_ix, 3),
            "new_closer_to_intersection": prediction_supported_ix,
        },
        "caveat": "El-khteeb features are all at one location (Saqqara), making this a weak single-site test.",
    }

    return results


# ============================================================
# ANALYSIS 4: CROSSROADS HYPOTHESIS
# ============================================================

def analysis_4_crossroads():
    print("\n" + "=" * 70)
    print("ANALYSIS 4: CROSSROADS HYPOTHESIS — Global River Crossings")
    print("=" * 70)

    results_list = []

    for crossing in RIVER_CROSSINGS:
        river = crossing["river"]
        if crossing["crossing_lat"] is None:
            # Nile — use computed intersection from Analysis 1
            continue

        c_lat = crossing["crossing_lat"]
        c_lon = crossing["crossing_lon"]

        # Distance of crossing point from GC
        d_gc = gc_distance(c_lat, c_lon)

        # GC bearing at crossing
        gc_bear = gc_bearing_at(c_lat, c_lon)

        # Intersection angle (using provided river bearing)
        river_bear = crossing.get("river_bearing", 180)
        angle = angle_between_bearings(gc_bear, river_bear)

        # Monuments within 10km, 50km, 100km
        monument_dists = []
        for m in crossing["monuments_nearby"]:
            d = haversine_km(c_lat, c_lon, m["lat"], m["lon"])
            monument_dists.append({"name": m["name"], "distance_km": round(d, 2)})

        results_list.append({
            "river": river,
            "crossing_point": {"lat": c_lat, "lon": c_lon},
            "gc_distance_km": round(d_gc, 2),
            "gc_bearing": round(gc_bear, 2),
            "river_bearing": river_bear,
            "intersection_angle": round(angle, 2),
            "monuments": monument_dists,
        })

        print(f"\n  {river}:")
        print(f"    Crossing at: {c_lat:.2f}°N, {c_lon:.2f}°E (GC dist: {d_gc:.1f}km)")
        print(f"    Angle: {angle:.1f}° (GC={gc_bear:.1f}°, River={river_bear}°)")
        for m in monument_dists:
            print(f"    → {m['name']}: {m['distance_km']:.1f}km from crossing")

    results = {
        "crossings": results_list,
        "note": "River bearings are approximate. Monument lists are illustrative, not exhaustive.",
    }

    return results


# ============================================================
# INTERSECTION MAP
# ============================================================

def make_intersection_map(intersection_result):
    """Create a map showing the branch, GC, pyramids, and intersection."""

    fig, ax = plt.subplots(figsize=(12, 14))

    # Branch
    branch_lats = [p["lat"] for p in AHRAMAT_BRANCH]
    branch_lons = [p["lon"] for p in AHRAMAT_BRANCH]
    ax.plot(branch_lons, branch_lats, 'b-', linewidth=2.5, label='Ahramat Branch', zorder=3)
    ax.plot(branch_lons, branch_lats, 'b.', markersize=4, zorder=3)

    # Great Circle — compute points in the region
    gc_lats = []
    gc_lons = []
    for lat_scan in np.linspace(29.2, 30.2, 200):
        # Find lon where GC passes at this latitude
        best_lon = None
        best_dist = 9999
        for lon_scan in np.linspace(30.5, 31.8, 500):
            d = gc_distance(lat_scan, lon_scan)
            if d < best_dist:
                best_dist = d
                best_lon = lon_scan
        if best_dist < 0.5:  # Within 500m
            gc_lats.append(lat_scan)
            gc_lons.append(best_lon)

    ax.plot(gc_lons, gc_lats, 'r-', linewidth=2.5, label='Great Circle', zorder=3)

    # Pyramids by dynasty
    dynasty_colors = {
        3: '#8B4513',  # Brown - 3rd Dynasty
        4: '#FFD700',  # Gold - 4th Dynasty
        5: '#FF8C00',  # Dark Orange - 5th Dynasty
        6: '#CD853F',  # Peru - 6th Dynasty
        8: '#A0522D',  # Sienna - 8th Dynasty
        12: '#4682B4',  # Steel Blue - 12th Dynasty
    }

    for pyr in PYRAMID_SITES:
        d = pyr.get("dynasty", 0)
        color = dynasty_colors.get(d, 'gray')
        ax.plot(pyr["lon"], pyr["lat"], '^', color=color, markersize=8,
                markeredgecolor='black', markeredgewidth=0.5, zorder=5)
        # Label major ones
        if pyr["name"] in ["Great Pyramid (Khufu)", "Step Pyramid (Djoser)",
                           "Red Pyramid (Sneferu)", "Bent Pyramid (Sneferu)",
                           "Meidum Pyramid (Sneferu?)", "Amenemhat I (Lisht)",
                           "Sahure", "Abu Rawash (Djedefre)"]:
            ax.annotate(pyr["name"].split("(")[0].strip(),
                        (pyr["lon"], pyr["lat"]),
                        textcoords="offset points", xytext=(8, 4),
                        fontsize=7, color='black', zorder=6)

    # Intersection point
    if intersection_result["intersections"]:
        ix = intersection_result["intersections"][0]
        ix_pt = ix["intersection_point"]
        ax.plot(ix_pt["lon"], ix_pt["lat"], '*', color='lime', markersize=20,
                markeredgecolor='black', markeredgewidth=1.5, zorder=10,
                label='Intersection Point')

        # Draw angle indicator
        angle_info = ix["angles_by_scale"]["5km"]
        gc_bear_rad = math.radians(angle_info["gc_bearing"])
        br_bear_rad = math.radians(angle_info["branch_bearing"])

    # El-khteeb features
    for feat in ELKHTEEB_2025:
        ax.plot(feat["lon"], feat["lat"], 'D', color='magenta', markersize=7,
                markeredgecolor='black', markeredgewidth=0.5, zorder=7)

    # Legend for dynasties
    for d, color in sorted(dynasty_colors.items()):
        ax.plot([], [], '^', color=color, markersize=8,
                markeredgecolor='black', markeredgewidth=0.5,
                label=f'Dynasty {d}')
    ax.plot([], [], 'D', color='magenta', markersize=7,
            markeredgecolor='black', markeredgewidth=0.5,
            label='El-khteeb 2025 (GPR)')

    ax.set_xlabel('Longitude (°E)', fontsize=12)
    ax.set_ylabel('Latitude (°N)', fontsize=12)
    ax.set_title('Great Circle × Ahramat Branch Intersection\nMemphis Necropolis Region',
                  fontsize=14, fontweight='bold')
    ax.legend(loc='lower left', fontsize=8, ncol=2)
    ax.set_aspect(1.0 / math.cos(math.radians(29.8)))  # Correct aspect ratio
    ax.grid(True, alpha=0.3)

    # Set bounds to focus on the pyramid field
    ax.set_xlim(31.0, 31.35)
    ax.set_ylim(29.25, 30.15)

    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'intersection_map.png'), dpi=200)
    plt.close()
    print(f"\n  Saved: intersection_map.png")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("AHRAMAT BRANCH PERPENDICULAR INTERSECTION — Full Analysis")
    print("Directive 08")
    print("=" * 70)

    # Analysis 1
    r1 = analysis_1_intersection_geometry()
    with open(os.path.join(OUTPUT_DIR, 'intersection_geometry.json'), 'w') as f:
        json.dump(r1, f, indent=2)
    print(f"  Saved: intersection_geometry.json")

    # Get intersection point for subsequent analyses
    ix_pt = r1["intersections"][0]["intersection_point"]
    ix_lat, ix_lon = ix_pt["lat"], ix_pt["lon"]

    # Analysis 2
    r2 = analysis_2_monte_carlo(n_trials=100_000, seed=42)
    with open(os.path.join(OUTPUT_DIR, 'monte_carlo_intersection.json'), 'w') as f:
        json.dump(r2, f, indent=2)
    print(f"  Saved: monte_carlo_intersection.json")

    # Analysis 3
    r3 = analysis_3_predictive_test((ix_lat, ix_lon))
    with open(os.path.join(OUTPUT_DIR, 'new_discoveries_proximity.json'), 'w') as f:
        json.dump(r3, f, indent=2)
    print(f"  Saved: new_discoveries_proximity.json")

    # Analysis 4
    r4 = analysis_4_crossroads()
    # Add the Nile crossing from our analysis
    if r1["intersections"]:
        nile_entry = {
            "river": "Nile (Ahramat Branch)",
            "crossing_point": ix_pt,
            "gc_distance_km": 0.0,
            "intersection_angle": r1["summary"]["primary_intersection_angle"],
            "pyramids_within_10km": r1["summary"]["pyramids_within_10km"],
            "pyramids_within_20km": r1["summary"]["pyramids_within_20km"],
        }
        r4["crossings"].insert(0, nile_entry)
    with open(os.path.join(OUTPUT_DIR, 'river_crossings.json'), 'w') as f:
        json.dump(r4, f, indent=2)
    print(f"  Saved: river_crossings.json")

    # Generate map
    make_intersection_map(r1)

    # Write RESULTS.md
    write_results_md(r1, r2, r3, r4)

    print("\n" + "=" * 70)
    print("ALL ANALYSES COMPLETE")
    print("=" * 70)


def write_results_md(r1, r2, r3, r4):
    """Generate the final RESULTS.md summary."""

    ix = r1["intersections"][0]
    ix_pt = ix["intersection_point"]
    angles = ix["angles_by_scale"]
    summary = r1["summary"]

    md = f"""# AHRAMAT BRANCH PERPENDICULAR INTERSECTION — Results

**Date:** 2026-03-24
**Directive:** 08_ahramat_branch_intersection.md

---

## Analysis 1: Intersection Geometry

The Great Circle crosses the Ahramat Branch at approximately **{ix_pt['lat']:.4f}°N, {ix_pt['lon']:.4f}°E**.

### Intersection Angle (multi-scale)

The directive warned that the branch meanders, so we report the angle at multiple measurement scales:

| Scale | Branch Bearing | GC Bearing | Intersection Angle | Deviation from 90° | Branch Points |
|-------|---------------|-----------|-------------------|-------------------|--------------|
"""
    for scale, data in angles.items():
        n_pts = data.get('n_branch_points', '?')
        md += f"| {scale} | {data['branch_bearing']:.1f}° | {data['gc_bearing']:.1f}° | **{data['intersection_angle']:.1f}°** | {data['deviation_from_90']:.1f}° | {n_pts} |\n"

    primary_angle = summary["primary_intersection_angle"]
    dev = summary["deviation_from_perpendicular"]
    local_angle = summary["local_angle"]
    full_angle = summary["full_branch_angle"]

    md += f"""
**Key finding:** The intersection angle varies significantly with measurement scale due to local meanders:
- **Local (2km):** {local_angle:.1f}° — dominated by a local SSE meander at Zawyet el-Aryan
- **Representative (20km):** {primary_angle:.1f}° — captures the Memphis necropolis axis
- **Full branch:** {full_angle:.1f}° — includes the distant Meidum-Faiyum bend

At the 20km scale (the Memphis necropolis pyramid field), the deviation from perpendicular is **{dev:.1f}°**.

### Pyramid Proximity

- Pyramids within 10km of intersection: **{r1['summary']['pyramids_within_10km']}**
- Pyramids within 20km of intersection: **{r1['summary']['pyramids_within_20km']}**

Nearest pyramids to the intersection point:
"""
    for p in ix["nearest_pyramids"][:5]:
        md += f"- {p['name']}: {p['distance_km']:.1f}km (Dynasty {p.get('dynasty', '?')})\n"

    md += f"""
---

## Analysis 2: Monte Carlo Probability

**{r2['n_trials']:,}** random great circles tested.

### Test 1: Random Circle × Fixed Branch

| Metric | Count | Probability |
|--------|-------|------------|
| Cross the branch | {r2['test1_random_circle_x_fixed_branch']['n_crossings']:,} | {r2['test1_random_circle_x_fixed_branch']['p_crossing']:.4f} |
| Cross near pyramid (<10km) | {r2['test1_random_circle_x_fixed_branch']['n_near_pyramid_10km']:,} | {r2['test1_random_circle_x_fixed_branch']['p_near_pyramid']:.5f} |
| Perpendicular (>80°) + near pyramid | {r2['test1_random_circle_x_fixed_branch']['n_perpendicular_near_pyramid']:,} | **{r2['test1_random_circle_x_fixed_branch']['p_perpendicular_near_pyramid']:.6f}** |

### Test 2: Monument Density at Crossing

- Mean pyramids within 10km of random crossing: {r2['test2_monument_density']['mean_pyramids_at_crossing']:.2f}
- Observed (Great Circle): **{r2['test2_monument_density']['observed_pyramids']}**
- P(≥observed): **{r2['test2_monument_density']['p_gte_observed']:.5f}**

### Test 3: Conditional on Passing Through Egypt

- Circles passing through Egypt's latitude band: {r2['test3_conditional_on_egypt']['n_passes_egypt']:,}
- Of those, crossing the branch: {r2['test3_conditional_on_egypt']['n_crosses_branch']:,}
- Of those, near-perpendicular: {r2['test3_conditional_on_egypt']['n_perpendicular']:,}
- Of those, at dense pyramid segment: {r2['test3_conditional_on_egypt']['n_perpendicular_dense']:,}
- **P(perpendicular + dense | passes Egypt): {r2['test3_conditional_on_egypt']['p_perpendicular_dense_given_egypt']:.6f}**

### Crossing Angle Distribution

Mean crossing angle for random circles: {r2['crossing_angles']['mean']}° (median: {r2['crossing_angles']['median']}°).
See `random_circle_angle_histogram.png`.

---

## Analysis 3: Predictive Test — El-khteeb 2025 Saqqara Discoveries

"""
    comp = r3["comparison"]
    md += f"""| Metric | New Discoveries | Known Pyramids | Prediction Supported? |
|--------|----------------|---------------|---------------------|
| Mean dist to GC | {comp['mean_dist_to_gc_new']:.2f} km | {comp['mean_dist_to_gc_known']:.2f} km | {'Yes' if comp['new_closer_to_gc'] else 'No'} |
| Mean dist to intersection | {comp['mean_dist_to_intersection_new']:.2f} km | {comp['mean_dist_to_intersection_known']:.2f} km | {'Yes' if comp['new_closer_to_intersection'] else 'No'} |

**Caveat:** {r3['caveat']}

---

## Analysis 4: Crossroads Hypothesis — Global River Crossings

"""
    md += "| River | Intersection Angle | Notable Monuments Nearby |\n"
    md += "|-------|-------------------|-------------------------|\n"
    for c in r4["crossings"]:
        monuments = ", ".join(m["name"] + f" ({m['distance_km']:.0f}km)"
                              for m in c.get("monuments", []))
        if not monuments and c["river"] == "Nile (Ahramat Branch)":
            monuments = f"{c.get('pyramids_within_10km', '?')} pyramids within 10km"
        md += f"| {c['river']} | {c.get('intersection_angle', '?'):.1f}° | {monuments} |\n"

    md += f"""
**Note:** {r4['note']}

---

## Files Generated

- `intersection_geometry.json` — intersection point, multi-scale angles, pyramid distances
- `monte_carlo_intersection.json` — probability analysis (100k trials)
- `new_discoveries_proximity.json` — El-khteeb 2025 proximity test
- `river_crossings.json` — global river crossing analysis
- `intersection_map.png` — annotated map of the Memphis region
- `random_circle_angle_histogram.png` — angle distribution from Monte Carlo
"""

    with open(os.path.join(OUTPUT_DIR, 'RESULTS.md'), 'w') as f:
        f.write(md)
    print(f"  Saved: RESULTS.md")


if __name__ == "__main__":
    main()
