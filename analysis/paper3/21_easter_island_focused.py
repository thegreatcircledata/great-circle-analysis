#!/usr/bin/env python3
"""
Easter Island Focused Analysis
===============================
Easter Island sits ~6-10 km from the Alison Great Circle and is the terminus
of the longest open-ocean voyage in human history.

Task 1: Position test — how anomalous is the circle passing this close?
Task 2: Ahu distribution — do ahu cluster on the circle-side of the island?
Task 3: Route context — Tahiti-Easter Island approach direction
Task 4: Joint probability — Egypt + Easter Island on the same circle

CULTURAL NOTE: Uses only published, publicly available coordinate data.
Results describe geometric properties, not claims about intent.
"""

import json, math, os, sys, time
import numpy as np
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================
POLE_LAT = 59.682122
POLE_LON = -138.646087
R = 6371.0  # km

# Key locations
EASTER_ISLAND = (-27.1127, -109.3497)  # centroid
MEMPHIS = (29.87, 31.22)  # Memphis/Saqqara

OUT_DIR = Path(__file__).resolve().parent.parent.parent / "paper4" / "easter_island"
os.makedirs(OUT_DIR, exist_ok=True)

# Ahu data (from easter_island_ahu_data.py)
AHU_DATA = [
    {"name": "Ahu Tongariki", "lat": -27.1222, "lon": -109.2728, "n_moai": 15, "platform_m": 220},
    {"name": "Ahu Akivi", "lat": -27.1150, "lon": -109.3950, "n_moai": 7, "platform_m": 70},
    {"name": "Ahu Nau Nau", "lat": -27.0736, "lon": -109.3221, "n_moai": 7, "platform_m": 65},
    {"name": "Ahu Tahai", "lat": -27.1395, "lon": -109.4280, "n_moai": 1, "platform_m": 8},
    {"name": "Ahu Ko Te Riku", "lat": -27.1390, "lon": -109.4290, "n_moai": 1, "platform_m": 10},
    {"name": "Ahu Vai Uri", "lat": -27.1398, "lon": -109.4275, "n_moai": 5, "platform_m": 25},
    {"name": "Ahu Hanga Te'e", "lat": -27.1453, "lon": -109.3428, "n_moai": 8, "platform_m": 86},
    {"name": "Ahu Vinapu", "lat": -27.1764, "lon": -109.4064, "n_moai": 6, "platform_m": 40},
    {"name": "Ahu Te Pito Kura", "lat": -27.0883, "lon": -109.2996, "n_moai": 1, "platform_m": 45},
    {"name": "Ahu Huri a Urenga", "lat": -27.1525, "lon": -109.3949, "n_moai": 1, "platform_m": 12},
    {"name": "Ahu Ature Huki", "lat": -27.0801, "lon": -109.3196, "n_moai": 1, "platform_m": 12},
    {"name": "Ahu Akahanga", "lat": -27.1488, "lon": -109.3367, "n_moai": 13, "platform_m": 60},
    {"name": "Ahu One Makihi", "lat": -27.1411, "lon": -109.2913, "n_moai": 1, "platform_m": 15},
    {"name": "Ahu Hanga Poukura", "lat": -27.1635, "lon": -109.3793, "n_moai": 5, "platform_m": 30},
    {"name": "Ahu Runga Va'e", "lat": -27.1667, "lon": -109.3667, "n_moai": 2, "platform_m": 15},
    {"name": "Ahu Vai Teka", "lat": -27.1155, "lon": -109.3945, "n_moai": 1, "platform_m": 16},
    {"name": "Ahu Vinapu II", "lat": -27.1760, "lon": -109.4058, "n_moai": 5, "platform_m": 35},
    {"name": "Ahu Hanga Te Tenga", "lat": -27.1650, "lon": -109.3600, "n_moai": 2, "platform_m": 20},
]

# ============================================================
# GEOMETRY HELPERS
# ============================================================
def dist_to_gc_km(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    """Distance from point to great circle defined by pole (km)."""
    la, lo = math.radians(lat), math.radians(lon)
    pla, plo = math.radians(pole_lat), math.radians(pole_lon)
    sx = math.cos(la) * math.cos(lo)
    sy = math.cos(la) * math.sin(lo)
    sz = math.sin(la)
    px = math.cos(pla) * math.cos(plo)
    py = math.cos(pla) * math.sin(plo)
    pz = math.sin(pla)
    dot = max(-1.0, min(1.0, sx*px + sy*py + sz*pz))
    ang = math.acos(dot)
    return abs(ang - math.pi/2) * R


def dist_to_gc_vec(lats, lons, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    """Vectorized distance from points to great circle (km)."""
    lats_r = np.radians(lats)
    lons_r = np.radians(lons)
    pla, plo = math.radians(pole_lat), math.radians(pole_lon)
    sx = np.cos(lats_r) * np.cos(lons_r)
    sy = np.cos(lats_r) * np.sin(lons_r)
    sz = np.sin(lats_r)
    px = math.cos(pla) * math.cos(plo)
    py = math.cos(pla) * math.sin(plo)
    pz = math.sin(pla)
    dot = np.clip(sx*px + sy*py + sz*pz, -1.0, 1.0)
    ang = np.arccos(dot)
    return np.abs(ang - math.pi/2) * R


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance between two points in km."""
    la1, lo1, la2, lo2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = la2 - la1
    dlon = lo2 - lo1
    a = math.sin(dlat/2)**2 + math.cos(la1)*math.cos(la2)*math.sin(dlon/2)**2
    return 2 * R * math.asin(math.sqrt(min(a, 1.0)))


def random_pole_uniform(rng):
    """Generate uniformly random pole on sphere. Returns (lat_rad, lon_rad)."""
    z = rng.uniform(-1, 1)
    theta = rng.uniform(0, 2 * math.pi)
    return math.asin(z), theta - math.pi


def gc_bearing_at_point(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    """Local bearing (azimuth) of the great circle at a point. Returns degrees [0,360)."""
    la, lo = math.radians(lat), math.radians(lon)
    pla, plo = math.radians(pole_lat), math.radians(pole_lon)
    p = np.array([math.cos(la)*math.cos(lo), math.cos(la)*math.sin(lo), math.sin(la)])
    pole = np.array([math.cos(pla)*math.cos(plo), math.cos(pla)*math.sin(plo), math.sin(pla)])
    tangent = np.cross(pole, p)
    norm = np.linalg.norm(tangent)
    if norm < 1e-10:
        return 0.0
    tangent /= norm
    north = np.array([-math.sin(la)*math.cos(lo), -math.sin(la)*math.sin(lo), math.cos(la)])
    east = np.array([-math.sin(lo), math.cos(lo), 0])
    bearing = math.degrees(math.atan2(np.dot(tangent, east), np.dot(tangent, north)))
    return bearing % 360


def signed_dist_to_gc(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    """Signed distance from point to GC (positive = toward pole, negative = away).
    Used to determine which side of the circle a point is on."""
    la, lo = math.radians(lat), math.radians(lon)
    pla, plo = math.radians(pole_lat), math.radians(pole_lon)
    sx = math.cos(la) * math.cos(lo)
    sy = math.cos(la) * math.sin(lo)
    sz = math.sin(la)
    px = math.cos(pla) * math.cos(plo)
    py = math.cos(pla) * math.sin(plo)
    pz = math.sin(pla)
    dot = max(-1.0, min(1.0, sx*px + sy*py + sz*pz))
    ang = math.acos(dot)
    return (ang - math.pi/2) * R  # positive = pole side


# ============================================================
# TASK 1: EASTER ISLAND POSITION TEST
# ============================================================
def run_task1():
    """Test how anomalous it is for a great circle to pass within ~10km of Easter Island."""
    print("="*70)
    print("TASK 1: EASTER ISLAND POSITION TEST")
    print("="*70)

    ei_lat, ei_lon = EASTER_ISLAND
    ei_dist = dist_to_gc_km(ei_lat, ei_lon)
    print(f"Easter Island centroid: {ei_lat}°, {ei_lon}°")
    print(f"Distance to Alison Great Circle: {ei_dist:.1f} km")

    # Test 1a: Among 10,000 random great circles, how many pass within this distance?
    N_MC = 10000
    threshold = ei_dist  # use actual distance as threshold
    thresholds = [6, 10, 25, 50]

    print(f"\nMonte Carlo: {N_MC} uniformly random great circles")
    rng = np.random.RandomState(42)

    mc_dists = []
    for i in range(N_MC):
        plat_r, plon_r = random_pole_uniform(rng)
        # Compute distance from Easter Island to this random circle
        pla, plo = plat_r, plon_r  # already in radians
        la, lo = math.radians(ei_lat), math.radians(ei_lon)
        sx = math.cos(la) * math.cos(lo)
        sy = math.cos(la) * math.sin(lo)
        sz = math.sin(la)
        px = math.cos(pla) * math.cos(plo)
        py = math.cos(pla) * math.sin(plo)
        pz = math.sin(pla)
        dot = max(-1.0, min(1.0, sx*px + sy*py + sz*pz))
        ang = math.acos(dot)
        d = abs(ang - math.pi/2) * R
        mc_dists.append(d)

    mc_dists = np.array(mc_dists)

    results_1a = {}
    for t in thresholds:
        n_within = int(np.sum(mc_dists <= t))
        p_val = n_within / N_MC
        pctile = float(np.mean(mc_dists >= ei_dist) * 100) if t == thresholds[0] else None
        print(f"  Random circles within {t}km of Easter Island: {n_within}/{N_MC} ({100*p_val:.2f}%)")
        results_1a[f"{t}km"] = {
            "n_within": n_within,
            "p_value": round(p_val, 6),
            "fraction": round(p_val, 6),
        }

    # Percentile of Alison circle
    pctile = float(np.mean(mc_dists >= ei_dist) * 100)
    print(f"\n  Alison circle distance: {ei_dist:.1f} km")
    print(f"  Percentile (closer = higher): {pctile:.1f}%")
    print(f"  Median random distance: {np.median(mc_dists):.0f} km")
    print(f"  Mean random distance: {np.mean(mc_dists):.0f} km")

    # Test 1b: Among random circles passing through Egypt, how many also pass near Easter Island?
    print(f"\n--- Joint Test: Circles through Egypt AND near Easter Island ---")
    memphis_lat, memphis_lon = MEMPHIS
    N_MC_JOINT = 100000  # more trials for joint probability

    n_through_egypt = 0
    n_joint = 0
    egypt_threshold = 100  # km
    ei_threshold = 50  # km

    rng2 = np.random.RandomState(123)
    for i in range(N_MC_JOINT):
        plat_r, plon_r = random_pole_uniform(rng2)
        pla, plo = plat_r, plon_r

        # Distance to Memphis
        la_m, lo_m = math.radians(memphis_lat), math.radians(memphis_lon)
        sx = math.cos(la_m) * math.cos(lo_m)
        sy = math.cos(la_m) * math.sin(lo_m)
        sz = math.sin(la_m)
        px = math.cos(pla) * math.cos(plo)
        py = math.cos(pla) * math.sin(plo)
        pz = math.sin(pla)
        dot_m = max(-1.0, min(1.0, sx*px + sy*py + sz*pz))
        d_egypt = abs(math.acos(dot_m) - math.pi/2) * R

        if d_egypt > egypt_threshold:
            continue
        n_through_egypt += 1

        # Distance to Easter Island
        la_e, lo_e = math.radians(ei_lat), math.radians(ei_lon)
        sx = math.cos(la_e) * math.cos(lo_e)
        sy = math.cos(la_e) * math.sin(lo_e)
        sz = math.sin(la_e)
        dot_e = max(-1.0, min(1.0, sx*px + sy*py + sz*pz))
        d_ei = abs(math.acos(dot_e) - math.pi/2) * R

        if d_ei <= ei_threshold:
            n_joint += 1

        if (i + 1) % 20000 == 0:
            print(f"  Trial {i+1}/{N_MC_JOINT} (Egypt hits: {n_through_egypt}, joint: {n_joint})")

    p_egypt = n_through_egypt / N_MC_JOINT
    p_joint = n_joint / N_MC_JOINT
    p_ei_given_egypt = n_joint / n_through_egypt if n_through_egypt > 0 else 0

    print(f"\n  Circles within {egypt_threshold}km of Memphis: {n_through_egypt}/{N_MC_JOINT} ({100*p_egypt:.2f}%)")
    print(f"  Of those, also within {ei_threshold}km of Easter Island: {n_joint}/{n_through_egypt}")
    print(f"  Joint probability: {p_joint:.6f} ({100*p_joint:.4f}%)")
    print(f"  P(Easter Island | Egypt): {p_ei_given_egypt:.6f}")

    # Alison circle performance
    alison_egypt = dist_to_gc_km(memphis_lat, memphis_lon)
    alison_ei = dist_to_gc_km(ei_lat, ei_lon)
    print(f"\n  Alison circle: Memphis={alison_egypt:.1f}km, Easter Island={alison_ei:.1f}km")

    result = {
        "description": "Easter Island position anomaly test",
        "easter_island_centroid": {"lat": ei_lat, "lon": ei_lon},
        "alison_distance_km": round(ei_dist, 2),
        "random_circle_test": {
            "n_trials": N_MC,
            "thresholds": results_1a,
            "percentile_closer": round(pctile, 1),
            "median_random_dist_km": round(float(np.median(mc_dists)), 1),
            "mean_random_dist_km": round(float(np.mean(mc_dists)), 1),
        },
        "joint_egypt_test": {
            "n_trials": N_MC_JOINT,
            "egypt_threshold_km": egypt_threshold,
            "ei_threshold_km": ei_threshold,
            "n_through_egypt": n_through_egypt,
            "n_joint": n_joint,
            "p_egypt": round(p_egypt, 6),
            "p_joint": round(p_joint, 6),
            "p_ei_given_egypt": round(p_ei_given_egypt, 6),
            "alison_memphis_km": round(alison_egypt, 1),
            "alison_easter_island_km": round(alison_ei, 1),
        },
    }

    out_path = OUT_DIR / "easter_island_position_test.json"
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {out_path}")

    return result


# ============================================================
# TASK 2: AHU DISTRIBUTION ANALYSIS
# ============================================================
def run_task2():
    """Test whether ahu cluster on the circle-side of the island."""
    print("\n" + "="*70)
    print("TASK 2: AHU DISTRIBUTION ANALYSIS")
    print("="*70)

    n_ahu = len(AHU_DATA)
    ei_lat, ei_lon = EASTER_ISLAND

    # Compute distances and signed distances for each ahu
    ahu_results = []
    for ahu in AHU_DATA:
        d = dist_to_gc_km(ahu["lat"], ahu["lon"])
        sd = signed_dist_to_gc(ahu["lat"], ahu["lon"])
        ahu_results.append({
            "name": ahu["name"],
            "lat": ahu["lat"],
            "lon": ahu["lon"],
            "n_moai": ahu["n_moai"],
            "platform_m": ahu["platform_m"],
            "dist_gc_km": round(d, 2),
            "signed_dist_km": round(sd, 2),  # positive = pole side (N of circle)
        })

    ahu_results.sort(key=lambda x: x["dist_gc_km"])

    print(f"Total ahu: {n_ahu}")
    print(f"Island centroid distance: {dist_to_gc_km(ei_lat, ei_lon):.1f} km")
    print(f"Island diameter: ~24 km (N-S) x ~12 km (E-W)")
    print(f"\nGC bearing at Easter Island: {gc_bearing_at_point(ei_lat, ei_lon):.1f}°")

    # Tight threshold counts
    tight_thresholds = [1, 5, 10, 15, 20]
    print(f"\nAhu within tight thresholds:")
    threshold_counts = {}
    for t in tight_thresholds:
        n = sum(1 for a in ahu_results if a["dist_gc_km"] <= t)
        print(f"  {t}km: {n}/{n_ahu} ({100*n/n_ahu:.0f}%)")
        threshold_counts[f"{t}km"] = {"count": n, "fraction": round(n/n_ahu, 3)}

    # Signed distance analysis: which side of the circle?
    signed_dists = [a["signed_dist_km"] for a in ahu_results]
    n_pole_side = sum(1 for d in signed_dists if d > 0)  # north of circle
    n_anti_side = sum(1 for d in signed_dists if d <= 0)  # south of circle

    print(f"\nSided analysis (GC crosses island NNE-SSW):")
    print(f"  Pole side (N of circle): {n_pole_side}")
    print(f"  Anti-pole side (S of circle): {n_anti_side}")
    print(f"  Mean signed distance: {np.mean(signed_dists):.1f} km")

    # Ahu distances
    dists = np.array([a["dist_gc_km"] for a in ahu_results])
    print(f"\nAhu distance statistics:")
    print(f"  Min: {dists.min():.1f} km ({ahu_results[0]['name']})")
    print(f"  Max: {dists.max():.1f} km ({ahu_results[-1]['name']})")
    print(f"  Mean: {dists.mean():.1f} km")
    print(f"  Median: {np.median(dists):.1f} km")

    # Monte Carlo: random great circles through Easter Island
    # A "circle through Easter Island" = pole perpendicular to island
    # We test: is the ahu distribution more concentrated near the Alison circle
    # than near a random circle passing through the island?
    print(f"\n--- Monte Carlo: Random circles through Easter Island ---")
    N_MC = 10000
    rng = np.random.RandomState(42)

    ahu_lats = np.array([a["lat"] for a in AHU_DATA])
    ahu_lons = np.array([a["lon"] for a in AHU_DATA])

    observed_mean_dist = float(dists.mean())
    observed_median_dist = float(np.median(dists))

    mc_mean_dists = []
    mc_within_10 = []

    for trial in range(N_MC):
        # Generate random circle that passes near Easter Island
        # Method: random pole constrained so circle passes within ~12km of island center
        # This means pole must be ~90° ± 0.1° from Easter Island
        # Simpler: generate random bearing through island, compute implied pole

        # Random azimuth through the island
        az = rng.uniform(0, 360)
        az_rad = math.radians(az)

        # Compute pole from point + bearing
        # Pole is 90° away from point, in direction perpendicular to bearing
        ei_la = math.radians(ei_lat)
        ei_lo = math.radians(ei_lon)

        # Perpendicular direction (bearing + 90°)
        perp_az = az_rad + math.pi/2

        # Pole at 90° distance in perpendicular direction
        pole_lat_r = math.asin(
            math.sin(ei_la) * math.cos(math.pi/2) +
            math.cos(ei_la) * math.sin(math.pi/2) * math.cos(perp_az)
        )
        pole_lon_r = ei_lo + math.atan2(
            math.sin(perp_az) * math.sin(math.pi/2) * math.cos(ei_la),
            math.cos(math.pi/2) - math.sin(ei_la) * math.sin(pole_lat_r)
        )

        pole_lat_d = math.degrees(pole_lat_r)
        pole_lon_d = math.degrees(pole_lon_r)

        # Compute ahu distances to this random circle
        trial_dists = dist_to_gc_vec(ahu_lats, ahu_lons, pole_lat_d, pole_lon_d)
        mc_mean_dists.append(float(np.mean(trial_dists)))
        mc_within_10.append(int(np.sum(trial_dists <= 10)))

    mc_mean_arr = np.array(mc_mean_dists)
    mc_w10_arr = np.array(mc_within_10)

    # Compare Alison circle to random circles through island
    z_mean = (observed_mean_dist - np.mean(mc_mean_arr)) / np.std(mc_mean_arr) if np.std(mc_mean_arr) > 0 else 0
    p_mean = float(np.mean(mc_mean_arr <= observed_mean_dist))

    observed_w10 = int(np.sum(dists <= 10))
    z_w10 = (observed_w10 - np.mean(mc_w10_arr)) / np.std(mc_w10_arr) if np.std(mc_w10_arr) > 0 else 0
    p_w10 = float(np.mean(mc_w10_arr >= observed_w10))

    print(f"  Alison mean ahu distance: {observed_mean_dist:.1f} km")
    print(f"  Random mean: {np.mean(mc_mean_arr):.1f} ± {np.std(mc_mean_arr):.1f} km")
    print(f"  Z-score (lower = closer): {z_mean:.2f}, p = {p_mean:.4f}")
    print(f"\n  Alison ahu within 10km: {observed_w10}")
    print(f"  Random mean within 10km: {np.mean(mc_w10_arr):.1f} ± {np.std(mc_w10_arr):.1f}")
    print(f"  Z-score (higher = more): {z_w10:.2f}, p = {p_w10:.4f}")

    # Bearing analysis: bearing from island center to each ahu
    gc_bearing = gc_bearing_at_point(ei_lat, ei_lon)
    print(f"\n--- Bearing Analysis ---")
    print(f"GC bearing at island: {gc_bearing:.1f}°")

    ahu_bearings = []
    for ahu in AHU_DATA:
        # Bearing from island center to ahu
        la1, lo1 = math.radians(ei_lat), math.radians(ei_lon)
        la2, lo2 = math.radians(ahu["lat"]), math.radians(ahu["lon"])
        dlon = lo2 - lo1
        x = math.sin(dlon) * math.cos(la2)
        y = math.cos(la1) * math.sin(la2) - math.sin(la1) * math.cos(la2) * math.cos(dlon)
        bearing = math.degrees(math.atan2(x, y)) % 360
        ahu_bearings.append(bearing)

    ahu_bearings = np.array(ahu_bearings)

    # Rayleigh test for directional preference
    theta = np.radians(ahu_bearings)
    C = np.mean(np.cos(theta))
    S = np.mean(np.sin(theta))
    R_bar = math.sqrt(C**2 + S**2)
    mean_dir = math.degrees(math.atan2(S, C)) % 360

    # Rayleigh test statistic
    n = len(theta)
    Z = n * R_bar**2
    # Rayleigh p-value approximation
    p_rayleigh = math.exp(-Z) * (1 + (2*Z - Z**2) / (4*n) - (24*Z - 132*Z**2 + 76*Z**3 - 9*Z**4) / (288*n**2))
    p_rayleigh = max(0, min(1, p_rayleigh))

    print(f"  Mean bearing from center to ahu: {mean_dir:.1f}°")
    print(f"  Mean resultant length R: {R_bar:.3f}")
    print(f"  Rayleigh Z: {Z:.2f}, p = {p_rayleigh:.4f}")
    print(f"  {'SIGNIFICANT' if p_rayleigh < 0.05 else 'Not significant'} directional preference")

    # Check if mean ahu direction aligns with GC bearing
    bearing_diff = abs(mean_dir - gc_bearing)
    if bearing_diff > 180:
        bearing_diff = 360 - bearing_diff
    print(f"  Difference from GC bearing: {bearing_diff:.1f}°")

    # Which ahu are closest to the GC line (within 10km)?
    print(f"\nClosest ahu to great circle:")
    for a in ahu_results[:5]:
        print(f"  {a['name']}: {a['dist_gc_km']:.1f} km ({a['n_moai']} moai, {a['platform_m']}m platform)")

    result = {
        "description": "Easter Island ahu distribution relative to the Alison Great Circle",
        "cultural_note": "Uses only published coordinate data. Results describe geometric properties.",
        "n_ahu": n_ahu,
        "island_centroid_dist_km": round(dist_to_gc_km(ei_lat, ei_lon), 2),
        "gc_bearing_at_island_deg": round(gc_bearing, 1),
        "ahu_distances": ahu_results,
        "distance_statistics": {
            "min_km": round(float(dists.min()), 2),
            "max_km": round(float(dists.max()), 2),
            "mean_km": round(float(dists.mean()), 2),
            "median_km": round(float(np.median(dists)), 2),
            "std_km": round(float(dists.std()), 2),
        },
        "threshold_counts": threshold_counts,
        "sided_analysis": {
            "n_pole_side": n_pole_side,
            "n_anti_pole_side": n_anti_side,
            "mean_signed_dist_km": round(float(np.mean(signed_dists)), 2),
        },
        "random_circle_comparison": {
            "description": "10,000 random great circles passing through Easter Island",
            "n_trials": N_MC,
            "mean_distance": {
                "observed": round(observed_mean_dist, 1),
                "random_mean": round(float(np.mean(mc_mean_arr)), 1),
                "random_std": round(float(np.std(mc_mean_arr)), 1),
                "z_score": round(z_mean, 2),
                "p_value": round(p_mean, 4),
            },
            "within_10km": {
                "observed": observed_w10,
                "random_mean": round(float(np.mean(mc_w10_arr)), 1),
                "random_std": round(float(np.std(mc_w10_arr)), 1),
                "z_score": round(z_w10, 2),
                "p_value": round(p_w10, 4),
            },
        },
        "bearing_analysis": {
            "gc_bearing_deg": round(gc_bearing, 1),
            "mean_ahu_bearing_deg": round(mean_dir, 1),
            "bearing_difference_deg": round(bearing_diff, 1),
            "rayleigh_R": round(R_bar, 4),
            "rayleigh_Z": round(Z, 2),
            "rayleigh_p": round(p_rayleigh, 4),
            "significant": bool(p_rayleigh < 0.05),
        },
    }

    out_path = OUT_DIR / "easter_island_ahu_analysis.json"
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {out_path}")

    return result


# ============================================================
# TASK 3: POLYNESIAN EXPANSION ROUTE CONTEXT
# ============================================================
def run_task3():
    """Analyze Tahiti-Easter Island route relative to great circle."""
    print("\n" + "="*70)
    print("TASK 3: POLYNESIAN EXPANSION ROUTE CONTEXT")
    print("="*70)

    ei_lat, ei_lon = EASTER_ISLAND
    gc_bearing = gc_bearing_at_point(ei_lat, ei_lon)

    # Tahiti-Easter Island route (from polynesian_navigation_analysis.py)
    route = [
        {"name": "Mangareva (Gambier)", "lat": -23.10, "lon": -134.97},
        {"name": "Henderson Island", "lat": -24.37, "lon": -128.33},
        {"name": "Pitcairn", "lat": -25.07, "lon": -130.10},
        {"name": "Easter Island (Anakena)", "lat": -27.07, "lon": -109.32},
    ]

    print(f"Tahiti-Easter Island route ({len(route)} waypoints)")
    print(f"GC bearing at Easter Island: {gc_bearing:.1f}°")

    # Compute distance of each waypoint to GC
    route_results = []
    for wp in route:
        d = dist_to_gc_km(wp["lat"], wp["lon"])
        print(f"  {wp['name']}: {d:.0f} km from circle")
        route_results.append({**wp, "dist_gc_km": round(d, 1)})

    # Interpolate route and find closest approach
    print(f"\n--- Route interpolation (50km steps) ---")
    interp_points = []
    for i in range(len(route) - 1):
        lat1, lon1 = route[i]["lat"], route[i]["lon"]
        lat2, lon2 = route[i+1]["lat"], route[i+1]["lon"]
        seg_dist = haversine_km(lat1, lon1, lat2, lon2)
        n_steps = max(2, int(seg_dist / 50))

        for j in range(n_steps + 1):
            frac = j / n_steps
            # Simple linear interpolation (good enough for this scale)
            lat = lat1 + frac * (lat2 - lat1)
            lon = lon1 + frac * (lon2 - lon1)
            d = dist_to_gc_km(lat, lon)
            interp_points.append({"lat": round(lat, 3), "lon": round(lon, 3),
                                  "dist_gc_km": round(d, 1)})

    closest = min(interp_points, key=lambda x: x["dist_gc_km"])
    print(f"  Closest approach: {closest['dist_gc_km']} km at ({closest['lat']}°, {closest['lon']}°)")

    # Approach bearing: bearing from Mangareva to Easter Island
    la1 = math.radians(route[0]["lat"])
    lo1 = math.radians(route[0]["lon"])
    la2 = math.radians(route[-1]["lat"])
    lo2 = math.radians(route[-1]["lon"])
    dlon = lo2 - lo1
    x = math.sin(dlon) * math.cos(la2)
    y = math.cos(la1) * math.sin(la2) - math.sin(la1) * math.cos(la2) * math.cos(dlon)
    approach_bearing = math.degrees(math.atan2(x, y)) % 360

    bearing_diff = abs(approach_bearing - gc_bearing)
    if bearing_diff > 180:
        bearing_diff = 360 - bearing_diff

    print(f"\n  Approach bearing (Mangareva → Easter Island): {approach_bearing:.1f}°")
    print(f"  GC bearing at Easter Island: {gc_bearing:.1f}°")
    print(f"  Difference: {bearing_diff:.1f}°")

    # Does the circle predict the approach direction?
    # The circle runs roughly ENE-WSW through Easter Island
    # Polynesian approach was from the WNW (Mangareva direction)
    approach_matches = bearing_diff < 30 or (180 - bearing_diff) < 30
    print(f"  Approach {'ROUGHLY ALIGNS' if approach_matches else 'does NOT align'} with GC bearing (±30°)")

    result = {
        "description": "Tahiti-Easter Island route analysis relative to great circle",
        "route_waypoints": route_results,
        "route_interpolation": {
            "n_points": len(interp_points),
            "closest_approach_km": closest["dist_gc_km"],
            "closest_point": closest,
        },
        "bearings": {
            "gc_bearing_at_ei_deg": round(gc_bearing, 1),
            "approach_bearing_deg": round(approach_bearing, 1),
            "difference_deg": round(bearing_diff, 1),
            "approach_aligns_with_gc": bool(approach_matches),
        },
        "interpretation": f"The Polynesian approach to Easter Island ({approach_bearing:.0f}°) "
                         f"{'roughly aligns' if approach_matches else 'does not align'} with the "
                         f"great circle bearing ({gc_bearing:.0f}°). Difference: {bearing_diff:.0f}°.",
    }

    out_path = OUT_DIR / "easter_island_route_context.json"
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {out_path}")

    return result


# ============================================================
# TASK 4: JOINT PROBABILITY — EGYPT + EASTER ISLAND
# ============================================================
def run_task4():
    """Among 100,000 random great circles, how many pass near both Memphis AND Easter Island?"""
    print("\n" + "="*70)
    print("TASK 4: JOINT PROBABILITY — EGYPT + EASTER ISLAND")
    print("="*70)

    ei_lat, ei_lon = EASTER_ISLAND
    mem_lat, mem_lon = MEMPHIS

    # Pre-compute target unit vectors
    ei_la, ei_lo = math.radians(ei_lat), math.radians(ei_lon)
    ei_x = math.cos(ei_la) * math.cos(ei_lo)
    ei_y = math.cos(ei_la) * math.sin(ei_lo)
    ei_z = math.sin(ei_la)

    mem_la, mem_lo = math.radians(mem_lat), math.radians(mem_lon)
    mem_x = math.cos(mem_la) * math.cos(mem_lo)
    mem_y = math.cos(mem_la) * math.sin(mem_lo)
    mem_z = math.sin(mem_la)

    # Test multiple threshold combinations
    test_configs = [
        {"egypt_km": 50, "ei_km": 10},
        {"egypt_km": 50, "ei_km": 50},
        {"egypt_km": 100, "ei_km": 50},
        {"egypt_km": 100, "ei_km": 100},
    ]

    N_MC = 100000
    print(f"Monte Carlo: {N_MC} uniformly random great circles\n")
    rng = np.random.RandomState(42)

    # Pre-generate all random poles
    z_vals = rng.uniform(-1, 1, N_MC)
    theta_vals = rng.uniform(0, 2 * math.pi, N_MC)
    pole_lats_r = np.arcsin(z_vals)
    pole_lons_r = theta_vals - math.pi

    # Vectorized distance computation
    px = np.cos(pole_lats_r) * np.cos(pole_lons_r)
    py = np.cos(pole_lats_r) * np.sin(pole_lons_r)
    pz = np.sin(pole_lats_r)

    # Distances to Egypt
    dot_egypt = np.clip(mem_x*px + mem_y*py + mem_z*pz, -1, 1)
    d_egypt = np.abs(np.arccos(dot_egypt) - math.pi/2) * R

    # Distances to Easter Island
    dot_ei = np.clip(ei_x*px + ei_y*py + ei_z*pz, -1, 1)
    d_ei = np.abs(np.arccos(dot_ei) - math.pi/2) * R

    # Alison circle performance
    alison_egypt = dist_to_gc_km(mem_lat, mem_lon)
    alison_ei = dist_to_gc_km(ei_lat, ei_lon)
    print(f"Alison circle: Memphis = {alison_egypt:.1f} km, Easter Island = {alison_ei:.1f} km\n")

    results = {}
    for cfg in test_configs:
        e_thresh = cfg["egypt_km"]
        ei_thresh = cfg["ei_km"]

        n_egypt = int(np.sum(d_egypt <= e_thresh))
        n_ei = int(np.sum(d_ei <= ei_thresh))
        n_joint = int(np.sum((d_egypt <= e_thresh) & (d_ei <= ei_thresh)))

        p_egypt = n_egypt / N_MC
        p_ei = n_ei / N_MC
        p_joint = n_joint / N_MC
        p_independent = p_egypt * p_ei  # expected if independent

        # Ratio: is joint more or less likely than independent expectation?
        ratio = p_joint / p_independent if p_independent > 0 else float('inf')

        label = f"Egypt<{e_thresh}km_EI<{ei_thresh}km"
        print(f"--- {label} ---")
        print(f"  P(Egypt): {n_egypt}/{N_MC} = {p_egypt:.4f}")
        print(f"  P(Easter Island): {n_ei}/{N_MC} = {p_ei:.4f}")
        print(f"  P(joint): {n_joint}/{N_MC} = {p_joint:.6f}")
        print(f"  P(independent): {p_independent:.6f}")
        print(f"  Ratio joint/independent: {ratio:.2f}")
        alison_passes = (alison_egypt <= e_thresh) and (alison_ei <= ei_thresh)
        print(f"  Alison passes both: {'YES' if alison_passes else 'NO'}")
        print()

        results[label] = {
            "egypt_threshold_km": e_thresh,
            "ei_threshold_km": ei_thresh,
            "n_egypt": n_egypt,
            "n_ei": n_ei,
            "n_joint": n_joint,
            "n_trials": N_MC,
            "p_egypt": round(p_egypt, 6),
            "p_ei": round(p_ei, 6),
            "p_joint": round(p_joint, 6),
            "p_independent": round(p_independent, 6),
            "ratio_joint_to_independent": round(ratio, 3),
            "alison_passes_both": bool(alison_passes),
        }

    # Key result for paper: tightest threshold where Alison passes
    out = {
        "description": "Joint probability test: random great circles passing near both Memphis (Egypt) AND Easter Island",
        "memphis": {"lat": mem_lat, "lon": mem_lon},
        "easter_island": {"lat": ei_lat, "lon": ei_lon},
        "alison_performance": {
            "memphis_km": round(alison_egypt, 1),
            "easter_island_km": round(alison_ei, 1),
        },
        "n_mc_trials": N_MC,
        "tests": results,
        "interpretation": "The Egypt-Easter Island joint constraint is geometrically non-trivial: "
                         "the two points are separated by ~180° of longitude and ~57° of latitude. "
                         "Any great circle through Egypt can only reach Easter Island if its bearing "
                         "at Memphis is approximately correct. Most circles through Egypt go elsewhere.",
    }

    out_path = OUT_DIR / "egypt_easter_island_joint_probability.json"
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"Saved to {out_path}")

    return out


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("Easter Island Focused Analysis")
    print("="*70)
    t0 = time.time()

    r1 = run_task1()
    r2 = run_task2()
    r3 = run_task3()
    r4 = run_task4()

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"COMPLETE — Total time: {elapsed:.0f}s")
    print(f"{'='*70}")
