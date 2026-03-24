#!/usr/bin/env python3
"""
Analysis 1b: Marae Enrichment Along the Great Circle (Monte Carlo)
====================================================================
Directive 09, Script 02

For each distance threshold (50, 100, 200, 500 km): count sacred sites
within the band. Compare against 10,000 random great circles crossing
the Pacific to compute enrichment significance.

Also tests: on islands near the circle, are sacred sites denser than
on islands farther away?

Output:
  - marae_enrichment.json
"""

import json
import math
import os
import sys
import random

import numpy as np

# ─── Import data module ──────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from easter_island_ahu_data import AHU_DATA
from polynesian_sacred_sites import POLYNESIAN_SACRED_SITES

# ─── Constants ────────────────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
R_EARTH = 6371.0
QUARTER_CIRC = R_EARTH * math.pi / 2
N_MC = 10000
THRESHOLDS_KM = [50, 100, 200, 500]

# ─── All sacred sites ────────────────────────────────────────────────
ALL_SITES = []
for ahu in AHU_DATA:
    ALL_SITES.append({"name": ahu["name"], "lat": ahu["lat"], "lon": ahu["lon"],
                       "island_group": "Easter Island"})
for site in POLYNESIAN_SACRED_SITES:
    ALL_SITES.append({"name": site["name"], "lat": site["lat"], "lon": site["lon"],
                       "island_group": site["island_group"]})

N_SITES = len(ALL_SITES)
SITE_LATS = np.array([s["lat"] for s in ALL_SITES])
SITE_LONS = np.array([s["lon"] for s in ALL_SITES])


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlat = la2 - la1
    dlon = lo2 - lo1
    a = math.sin(dlat/2)**2 + math.cos(la1)*math.cos(la2)*math.sin(dlon/2)**2
    return R_EARTH * 2 * math.asin(math.sqrt(min(1.0, a)))


def dist_from_gc(lat, lon, pole_lat, pole_lon):
    """Distance (km) from a point to a great circle defined by its pole."""
    dist_to_pole = haversine_km(lat, lon, pole_lat, pole_lon)
    return abs(dist_to_pole - QUARTER_CIRC)


def count_within_band(pole_lat, pole_lon, threshold_km):
    """Count sites within threshold_km of the great circle defined by pole."""
    count = 0
    for i in range(N_SITES):
        d = dist_from_gc(SITE_LATS[i], SITE_LONS[i], pole_lat, pole_lon)
        if d <= threshold_km:
            count += 1
    return count


def is_pacific_crossing(pole_lat, pole_lon):
    """
    Check if a great circle crosses the Pacific (rough filter).
    A circle crosses the Pacific if it has points with longitude between
    -180 and -100 (eastern Pacific) and latitude between -50 and 30.
    We test by checking if the circle passes near Easter Island's latitude band.
    """
    # Sample 36 points on the circle and check if any are in the Pacific
    plat = math.radians(pole_lat)
    plon = math.radians(pole_lon)
    px = math.cos(plat) * math.cos(plon)
    py = math.cos(plat) * math.sin(plon)
    pz = math.sin(plat)

    if abs(pz) < 0.999:
        v1x, v1y, v1z = -py, px, 0.0
    else:
        v1x, v1y, v1z = 1.0, 0.0, 0.0
    norm = math.sqrt(v1x**2 + v1y**2 + v1z**2)
    v1x /= norm; v1y /= norm; v1z /= norm
    v2x = py*v1z - pz*v1y
    v2y = pz*v1x - px*v1z
    v2z = px*v1y - py*v1x

    pacific_count = 0
    for i in range(36):
        theta = math.radians(i * 10)
        qx = v1x * math.cos(theta) + v2x * math.sin(theta)
        qy = v1y * math.cos(theta) + v2y * math.sin(theta)
        qz = v1z * math.cos(theta) + v2z * math.sin(theta)
        lat_q = math.degrees(math.asin(max(-1, min(1, qz))))
        lon_q = math.degrees(math.atan2(qy, qx))
        # Pacific: roughly lon in (-180, -100) and lat in (-50, 30)
        if -180 <= lon_q <= -100 and -50 <= lat_q <= 30:
            pacific_count += 1

    return pacific_count >= 3  # at least 3 of 36 sample points in Pacific


# ═══════════════════════════════════════════════════════════════════════
# Actual counts for the Great Circle
# ═══════════════════════════════════════════════════════════════════════

print("=" * 70)
print("MARAE ENRICHMENT ALONG THE GREAT CIRCLE")
print("=" * 70)
print(f"\nTotal sacred sites: {N_SITES}")
print(f"Monte Carlo trials: {N_MC}")

actual_counts = {}
print(f"\n--- Great Circle site counts ---")
for t in THRESHOLDS_KM:
    c = count_within_band(POLE_LAT, POLE_LON, t)
    actual_counts[t] = c
    # List the sites
    sites_in_band = []
    for s in ALL_SITES:
        d = dist_from_gc(s["lat"], s["lon"], POLE_LAT, POLE_LON)
        if d <= t:
            sites_in_band.append(f"{s['name']} ({d:.0f} km)")
    print(f"\n  ≤{t:4d} km: {c} sites")
    for s_str in sites_in_band:
        print(f"    {s_str}")


# ═══════════════════════════════════════════════════════════════════════
# Monte Carlo: random great circles crossing the Pacific
# ═══════════════════════════════════════════════════════════════════════

print(f"\n--- Monte Carlo: generating {N_MC} random Pacific-crossing circles ---")

np.random.seed(42)
mc_counts = {t: [] for t in THRESHOLDS_KM}
n_generated = 0
n_attempts = 0

while n_generated < N_MC:
    # Random pole on the unit sphere
    z = np.random.uniform(-1, 1)
    phi = np.random.uniform(0, 2 * math.pi)
    pole_lat_r = math.degrees(math.asin(z))
    pole_lon_r = math.degrees(phi) - 180

    n_attempts += 1

    # Filter: must cross the Pacific
    if not is_pacific_crossing(pole_lat_r, pole_lon_r):
        continue

    for t in THRESHOLDS_KM:
        mc_counts[t].append(count_within_band(pole_lat_r, pole_lon_r, t))

    n_generated += 1
    if n_generated % 1000 == 0:
        print(f"  {n_generated}/{N_MC} circles generated ({n_attempts} attempts)")

print(f"  Done. {n_attempts} attempts for {N_MC} Pacific-crossing circles.")

# ═══════════════════════════════════════════════════════════════════════
# Compute enrichment and p-values
# ═══════════════════════════════════════════════════════════════════════

print(f"\n--- Enrichment results ---")
enrichment_results = {}

for t in THRESHOLDS_KM:
    mc = np.array(mc_counts[t])
    actual = actual_counts[t]
    mc_mean = np.mean(mc)
    mc_std = np.std(mc)
    mc_median = np.median(mc)

    # p-value: fraction of MC trials with count >= actual
    p_value = np.mean(mc >= actual)

    # Enrichment ratio
    enrichment = actual / mc_mean if mc_mean > 0 else float('inf')

    enrichment_results[str(t)] = {
        "threshold_km": t,
        "actual_count": int(actual),
        "mc_mean": round(float(mc_mean), 2),
        "mc_median": round(float(mc_median), 1),
        "mc_std": round(float(mc_std), 2),
        "mc_max": int(mc.max()),
        "mc_min": int(mc.min()),
        "enrichment_ratio": round(enrichment, 3),
        "p_value": round(float(p_value), 4),
        "significant_05": bool(p_value < 0.05),
        "significant_01": bool(p_value < 0.01),
    }

    sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
    print(f"\n  ≤{t} km band:")
    print(f"    Actual:     {actual} sites")
    print(f"    MC mean:    {mc_mean:.1f} ± {mc_std:.1f}")
    print(f"    MC range:   {mc.min()} – {mc.max()}")
    print(f"    Enrichment: {enrichment:.2f}x")
    print(f"    p-value:    {p_value:.4f} {sig}")


# ═══════════════════════════════════════════════════════════════════════
# Per-island-group density analysis
# ═══════════════════════════════════════════════════════════════════════

print(f"\n--- Per-island-group site density vs circle distance ---")

group_stats = {}
for s in ALL_SITES:
    g = s["island_group"]
    if g not in group_stats:
        group_stats[g] = {"count": 0, "min_dist_km": float('inf')}
    group_stats[g]["count"] += 1

# Compute distance for each group centroid
group_centroids = {
    "Easter Island": (-27.11, -109.35),
    "Society Islands": (-17.32, -149.75),
    "Hawaiian Islands": (20.79, -156.33),
    "Marquesas Islands": (-9.43, -140.07),
    "Tonga": (-20.42, -175.20),
    "Samoa": (-13.76, -172.10),
    "Cook Islands": (-19.00, -159.50),
    "New Zealand": (-41.27, 174.78),
}

print(f"\n  {'Group':30s}  {'Sites':>6s}  {'Dist (km)':>10s}")
print(f"  {'-'*30}  {'-'*6}  {'-'*10}")

density_data = []
for g, (lat, lon) in sorted(group_centroids.items()):
    d = dist_from_gc(lat, lon, POLE_LAT, POLE_LON)
    n = group_stats.get(g, {}).get("count", 0)
    print(f"  {g:30s}  {n:>6d}  {d:>10.1f}")
    density_data.append({"group": g, "n_sites": n, "distance_km": round(d, 1)})

# ═══════════════════════════════════════════════════════════════════════
# Save results
# ═══════════════════════════════════════════════════════════════════════

output = {
    "analysis": "Marae enrichment along the Great Circle",
    "n_sites": N_SITES,
    "n_mc_trials": N_MC,
    "pole": {"lat": POLE_LAT, "lon": POLE_LON},
    "enrichment_by_threshold": enrichment_results,
    "per_group_density": density_data,
    "caveat": (
        "The Polynesian sacred site database is incomplete (~41 sites). "
        "Many marae have been destroyed or are unrecorded. Easter Island "
        "(18 ahu) dominates the dataset and sits directly on the circle, "
        "which inflates enrichment at narrow bands. The 500 km band is "
        "more meaningful for Pacific-scale analysis."
    ),
}

with open("marae_enrichment.json", "w") as f:
    json.dump(output, f, indent=2)

print(f"\nSaved marae_enrichment.json")
