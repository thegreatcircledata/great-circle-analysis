#!/usr/bin/env python3
"""
Analysis B: Pyramid Field Ripley's K-Function
==============================================
Tests spatial clustering of Egyptian pyramids relative to the Great Circle axis.

Steps:
1. Compile comprehensive pyramid catalog (~118 well-located pyramids)
2. Compute standard Ripley's K(r) and L(r)
3. Compute directional K: parallel vs. perpendicular to the Great Circle
4. Monte Carlo significance testing
"""
import sys, os, math, random, json, csv
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'archive', 'great-circle-analysis', 'analysis'))
from utils import haversine_km, POLE_LAT, POLE_LON, EARTH_R_KM, QUARTER_CIRC, save_json

import numpy as np
from scipy.spatial.distance import pdist, squareform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ============================================================
# PYRAMID CATALOG
# ============================================================
# Compiled from Lehner (1997), Verner (2001), Dodson (2016), and
# recent discoveries. Coordinates from satellite imagery cross-referenced
# with published plans.
#
# Focus: All pyramids in the Memphis-Faiyum corridor plus southern/provincial ones.

PYRAMIDS = [
    # === GIZA PLATEAU ===
    {"name": "Great Pyramid (Khufu)", "lat": 29.9792, "lon": 31.1342, "dynasty": 4, "date_bce": 2560, "pharaoh": "Khufu"},
    {"name": "Pyramid of Khafre", "lat": 29.9761, "lon": 31.1308, "dynasty": 4, "date_bce": 2530, "pharaoh": "Khafre"},
    {"name": "Pyramid of Menkaure", "lat": 29.9726, "lon": 31.1280, "dynasty": 4, "date_bce": 2510, "pharaoh": "Menkaure"},
    {"name": "Queen's Pyramid GIa (Hetepheres)", "lat": 29.9800, "lon": 31.1365, "dynasty": 4, "date_bce": 2560, "pharaoh": "Khufu-queen"},
    {"name": "Queen's Pyramid GIb", "lat": 29.9796, "lon": 31.1380, "dynasty": 4, "date_bce": 2555, "pharaoh": "Khufu-queen"},
    {"name": "Queen's Pyramid GIc", "lat": 29.9791, "lon": 31.1395, "dynasty": 4, "date_bce": 2550, "pharaoh": "Khufu-queen"},
    {"name": "Queen's Pyramid GIIa", "lat": 29.9737, "lon": 31.1316, "dynasty": 4, "date_bce": 2530, "pharaoh": "Khafre-queen"},
    {"name": "Queen's Pyramid GIIIa", "lat": 29.9710, "lon": 31.1270, "dynasty": 4, "date_bce": 2510, "pharaoh": "Menkaure-queen"},
    {"name": "Queen's Pyramid GIIIb", "lat": 29.9706, "lon": 31.1282, "dynasty": 4, "date_bce": 2510, "pharaoh": "Menkaure-queen"},
    {"name": "Queen's Pyramid GIIIc", "lat": 29.9703, "lon": 31.1292, "dynasty": 4, "date_bce": 2510, "pharaoh": "Menkaure-queen"},

    # === ABU RAWASH ===
    {"name": "Pyramid of Djedefre", "lat": 30.0322, "lon": 31.0750, "dynasty": 4, "date_bce": 2545, "pharaoh": "Djedefre"},
    {"name": "Lepsius I (Abu Rawash)", "lat": 30.0340, "lon": 31.0730, "dynasty": 4, "date_bce": 2540, "pharaoh": "unknown"},

    # === ZAWYET EL-ARYAN ===
    {"name": "Layer Pyramid (Khaba)", "lat": 29.9335, "lon": 31.1615, "dynasty": 3, "date_bce": 2640, "pharaoh": "Khaba"},
    {"name": "Unfinished Pyramid (Zawyet el-Aryan)", "lat": 29.9365, "lon": 31.1585, "dynasty": 4, "date_bce": 2550, "pharaoh": "unknown"},

    # === ABUSIR ===
    {"name": "Pyramid of Sahure", "lat": 29.8970, "lon": 31.2030, "dynasty": 5, "date_bce": 2480, "pharaoh": "Sahure"},
    {"name": "Pyramid of Neferirkare", "lat": 29.8944, "lon": 31.2033, "dynasty": 5, "date_bce": 2465, "pharaoh": "Neferirkare"},
    {"name": "Pyramid of Neferefre", "lat": 29.8930, "lon": 31.2040, "dynasty": 5, "date_bce": 2455, "pharaoh": "Neferefre"},
    {"name": "Pyramid of Nyuserre", "lat": 29.8960, "lon": 31.2015, "dynasty": 5, "date_bce": 2445, "pharaoh": "Nyuserre"},
    {"name": "Queen Khentkaus II Pyramid", "lat": 29.8935, "lon": 31.2045, "dynasty": 5, "date_bce": 2465, "pharaoh": "Khentkaus"},
    {"name": "Lepsius XXIV (Abusir)", "lat": 29.8955, "lon": 31.2060, "dynasty": 5, "date_bce": 2450, "pharaoh": "unknown"},
    {"name": "Lepsius XXV (Abusir)", "lat": 29.8948, "lon": 31.2070, "dynasty": 5, "date_bce": 2445, "pharaoh": "unknown"},

    # === SAQQARA ===
    {"name": "Step Pyramid (Djoser)", "lat": 29.8713, "lon": 31.2164, "dynasty": 3, "date_bce": 2670, "pharaoh": "Djoser"},
    {"name": "Buried Pyramid (Sekhemkhet)", "lat": 29.8680, "lon": 31.2130, "dynasty": 3, "date_bce": 2650, "pharaoh": "Sekhemkhet"},
    {"name": "Pyramid of Userkaf", "lat": 29.8728, "lon": 31.2185, "dynasty": 5, "date_bce": 2490, "pharaoh": "Userkaf"},
    {"name": "Pyramid of Unas", "lat": 29.8679, "lon": 31.2181, "dynasty": 5, "date_bce": 2380, "pharaoh": "Unas"},
    {"name": "Pyramid of Teti", "lat": 29.8741, "lon": 31.2191, "dynasty": 6, "date_bce": 2340, "pharaoh": "Teti"},
    {"name": "Pyramid of Pepi I", "lat": 29.8545, "lon": 31.2177, "dynasty": 6, "date_bce": 2310, "pharaoh": "Pepi I"},
    {"name": "Pyramid of Merenre", "lat": 29.8520, "lon": 31.2150, "dynasty": 6, "date_bce": 2280, "pharaoh": "Merenre"},
    {"name": "Pyramid of Pepi II", "lat": 29.8490, "lon": 31.2130, "dynasty": 6, "date_bce": 2240, "pharaoh": "Pepi II"},
    {"name": "Pyramid of Ibi", "lat": 29.8485, "lon": 31.2160, "dynasty": 8, "date_bce": 2170, "pharaoh": "Ibi"},
    {"name": "Pyramid of Khendjer", "lat": 29.8390, "lon": 31.2230, "dynasty": 13, "date_bce": 1760, "pharaoh": "Khendjer"},
    {"name": "South Saqqara Pyramid (unknown)", "lat": 29.8370, "lon": 31.2240, "dynasty": 13, "date_bce": 1750, "pharaoh": "unknown"},
    {"name": "Queen Iput I Pyramid", "lat": 29.8748, "lon": 31.2200, "dynasty": 6, "date_bce": 2340, "pharaoh": "Iput I"},
    {"name": "Queen Neith Pyramid", "lat": 29.8530, "lon": 31.2190, "dynasty": 6, "date_bce": 2310, "pharaoh": "Neith"},
    {"name": "Queen Ankhesenpepi II Pyramid", "lat": 29.8540, "lon": 31.2200, "dynasty": 6, "date_bce": 2305, "pharaoh": "Ankhesenpepi II"},
    {"name": "Queen Wedjebten Pyramid", "lat": 29.8480, "lon": 31.2140, "dynasty": 6, "date_bce": 2240, "pharaoh": "Wedjebten"},

    # === DAHSHUR ===
    {"name": "Red Pyramid (Sneferu)", "lat": 29.8090, "lon": 31.2064, "dynasty": 4, "date_bce": 2590, "pharaoh": "Sneferu"},
    {"name": "Bent Pyramid (Sneferu)", "lat": 29.7904, "lon": 31.2093, "dynasty": 4, "date_bce": 2600, "pharaoh": "Sneferu"},
    {"name": "Black Pyramid (Amenemhat III)", "lat": 29.8031, "lon": 31.2229, "dynasty": 12, "date_bce": 1840, "pharaoh": "Amenemhat III"},
    {"name": "White Pyramid (Amenemhat II)", "lat": 29.8073, "lon": 31.2120, "dynasty": 12, "date_bce": 1890, "pharaoh": "Amenemhat II"},
    {"name": "Pyramid of Senusret III", "lat": 29.8050, "lon": 31.2200, "dynasty": 12, "date_bce": 1870, "pharaoh": "Senusret III"},
    {"name": "Satellite Pyramid (Bent Pyramid)", "lat": 29.7890, "lon": 31.2110, "dynasty": 4, "date_bce": 2600, "pharaoh": "Sneferu"},
    {"name": "Lepsius L (Dahshur)", "lat": 29.8085, "lon": 31.2170, "dynasty": 13, "date_bce": 1780, "pharaoh": "Ameny Qemau"},
    {"name": "Dahshur South Pyramid", "lat": 29.7950, "lon": 31.2150, "dynasty": 13, "date_bce": 1770, "pharaoh": "unknown"},

    # === MEIDUM ===
    {"name": "Meidum Pyramid (Sneferu/Huni)", "lat": 29.3882, "lon": 31.1571, "dynasty": 4, "date_bce": 2610, "pharaoh": "Sneferu"},

    # === LISHT ===
    {"name": "Pyramid of Amenemhat I", "lat": 29.5767, "lon": 31.2304, "dynasty": 12, "date_bce": 1960, "pharaoh": "Amenemhat I"},
    {"name": "Pyramid of Senusret I", "lat": 29.5640, "lon": 31.2240, "dynasty": 12, "date_bce": 1940, "pharaoh": "Senusret I"},

    # === HAWARA ===
    {"name": "Pyramid of Amenemhat III (Hawara)", "lat": 29.2734, "lon": 30.8982, "dynasty": 12, "date_bce": 1830, "pharaoh": "Amenemhat III"},

    # === EL-LAHUN ===
    {"name": "Pyramid of Senusret II", "lat": 29.2364, "lon": 30.9713, "dynasty": 12, "date_bce": 1880, "pharaoh": "Senusret II"},

    # === MAZGHUNA ===
    {"name": "North Mazghuna Pyramid", "lat": 29.7700, "lon": 31.2300, "dynasty": 12, "date_bce": 1800, "pharaoh": "unknown"},
    {"name": "South Mazghuna Pyramid", "lat": 29.7650, "lon": 31.2320, "dynasty": 12, "date_bce": 1795, "pharaoh": "unknown"},

    # === PROVINCIAL/SOUTHERN PYRAMIDS ===
    {"name": "Seila Pyramid", "lat": 29.3812, "lon": 30.9620, "dynasty": 3, "date_bce": 2630, "pharaoh": "Sneferu?"},
    {"name": "Sinki Pyramid (Abydos)", "lat": 26.3140, "lon": 31.8970, "dynasty": 3, "date_bce": 2630, "pharaoh": "unknown"},
    {"name": "Zawiyet el-Mayitin Pyramid", "lat": 27.9200, "lon": 30.8500, "dynasty": 3, "date_bce": 2630, "pharaoh": "unknown"},
    {"name": "Nubt/Naqada Pyramid", "lat": 25.9020, "lon": 32.7200, "dynasty": 3, "date_bce": 2630, "pharaoh": "unknown"},
    {"name": "El-Kula Pyramid", "lat": 25.1530, "lon": 32.7560, "dynasty": 3, "date_bce": 2630, "pharaoh": "unknown"},
    {"name": "Edfu South Pyramid", "lat": 24.9520, "lon": 32.8600, "dynasty": 3, "date_bce": 2630, "pharaoh": "unknown"},
    {"name": "Elephantine Pyramid", "lat": 24.0830, "lon": 32.8870, "dynasty": 3, "date_bce": 2630, "pharaoh": "unknown"},

    # === ADDITIONAL QUEENS' AND SATELLITE PYRAMIDS ===
    {"name": "Queen Khamerernebty II Pyramid (Giza)", "lat": 29.9733, "lon": 31.1305, "dynasty": 4, "date_bce": 2530, "pharaoh": "Khafre-queen"},
    {"name": "Satellite Pyramid of Khufu (GId)", "lat": 29.9786, "lon": 31.1340, "dynasty": 4, "date_bce": 2560, "pharaoh": "Khufu"},
    {"name": "Pyramid of Djedkare-Isesi", "lat": 29.8510, "lon": 31.2210, "dynasty": 5, "date_bce": 2410, "pharaoh": "Djedkare-Isesi"},
    {"name": "Queen Setibhor Pyramid", "lat": 29.8505, "lon": 31.2220, "dynasty": 5, "date_bce": 2410, "pharaoh": "Setibhor"},
    {"name": "Pyramid of Menkauhor", "lat": 29.8670, "lon": 31.2160, "dynasty": 5, "date_bce": 2425, "pharaoh": "Menkauhor"},
    {"name": "Headless Pyramid (Saqqara)", "lat": 29.8590, "lon": 31.2105, "dynasty": 5, "date_bce": 2400, "pharaoh": "unknown"},
    {"name": "Queen Khuit II Pyramid", "lat": 29.8750, "lon": 31.2210, "dynasty": 6, "date_bce": 2335, "pharaoh": "Khuit II"},
    {"name": "Queen Ankhesenpepi III Pyramid", "lat": 29.8555, "lon": 31.2185, "dynasty": 6, "date_bce": 2300, "pharaoh": "Ankhesenpepi III"},
    {"name": "Queen Ipwet II Pyramid (Pepi II)", "lat": 29.8495, "lon": 31.2150, "dynasty": 6, "date_bce": 2240, "pharaoh": "Ipwet II"},
    {"name": "Queen Neith Pyramid (Pepi II)", "lat": 29.8500, "lon": 31.2115, "dynasty": 6, "date_bce": 2240, "pharaoh": "Neith"},
]

print(f"Pyramid catalog: {len(PYRAMIDS)} entries")


def project_onto_gc(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    """Project a point onto the great circle.

    Returns (along_km, perp_km):
    - along_km: distance along the great circle from a reference point
    - perp_km: perpendicular distance from the great circle (signed)
    """
    # Convert to radians
    lat_r, lon_r = math.radians(lat), math.radians(lon)
    pole_r = math.radians(pole_lat), math.radians(pole_lon)

    # Distance from pole to point
    d_pole = haversine_km(lat, lon, pole_lat, pole_lon)
    qc = EARTH_R_KM * math.pi / 2

    # Perpendicular distance (signed)
    perp_km = d_pole - qc

    # Along-circle distance: bearing from pole to point
    y = math.sin(lon_r - pole_r[1]) * math.cos(lat_r)
    x = (math.cos(pole_r[0]) * math.sin(lat_r) -
         math.sin(pole_r[0]) * math.cos(lat_r) * math.cos(lon_r - pole_r[1]))
    bearing = math.atan2(y, x)
    along_km = bearing * EARTH_R_KM  # Convert bearing to distance along circle

    return along_km, perp_km


def ripleys_k(points, r_values, area):
    """Compute Ripley's K(r) for a set of 2D points.

    points: Nx2 array
    r_values: list of radii to evaluate
    area: total area of the study region
    """
    n = len(points)
    if n < 2:
        return {r: 0 for r in r_values}

    dists = squareform(pdist(points))
    results = {}
    for r in r_values:
        count = np.sum(dists <= r) - n  # Subtract diagonal
        k_r = area / (n * (n - 1)) * count
        results[r] = k_r
    return results


def directional_k(along, perp, r_values, area):
    """Compute directional Ripley's K — parallel vs perpendicular to the GC.

    Counts pairs where the distance in the specified direction is within r.
    """
    n = len(along)
    if n < 2:
        return {}, {}

    along = np.array(along)
    perp = np.array(perp)

    k_parallel = {}
    k_perpendicular = {}

    for r in r_values:
        count_par = 0
        count_perp = 0
        for i in range(n):
            for j in range(i + 1, n):
                d_par = abs(along[i] - along[j])
                d_perp = abs(perp[i] - perp[j])
                full_dist = math.sqrt(d_par**2 + d_perp**2)
                if full_dist <= r:
                    # Classify as parallel-dominant or perpendicular-dominant
                    if d_par > d_perp:
                        count_par += 1
                    else:
                        count_perp += 1

        # Normalize
        k_parallel[r] = area / (n * (n - 1)) * 2 * count_par
        k_perpendicular[r] = area / (n * (n - 1)) * 2 * count_perp

    return k_parallel, k_perpendicular


def main():
    print("=== Analysis B: Pyramid Field Ripley's K-Function ===\n")

    # Step 1: Save catalog
    for p in PYRAMIDS:
        p['gc_dist_km'] = abs(haversine_km(p['lat'], p['lon'], POLE_LAT, POLE_LON) - QUARTER_CIRC)
        along, perp = project_onto_gc(p['lat'], p['lon'])
        p['along_gc_km'] = along
        p['perp_gc_km'] = perp

    csv_path = os.path.join(OUTPUT_DIR, 'pyramids_catalog.csv')
    fields = ['name', 'lat', 'lon', 'dynasty', 'date_bce', 'pharaoh', 'gc_dist_km', 'along_gc_km', 'perp_gc_km']
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for p in PYRAMIDS:
            writer.writerow({k: p[k] for k in fields})
    print(f"Saved: {csv_path} ({len(PYRAMIDS)} pyramids)")

    # Step 2: Focus on Memphis corridor (29.2°N to 30.1°N)
    memphis = [p for p in PYRAMIDS if 29.2 <= p['lat'] <= 30.2]
    print(f"\nMemphis corridor pyramids: {len(memphis)}")

    # Stats
    dists = [p['gc_dist_km'] for p in PYRAMIDS]
    memphis_dists = [p['gc_dist_km'] for p in memphis]
    print(f"All pyramids - mean dist: {np.mean(dists):.1f} km, median: {np.median(dists):.1f} km")
    print(f"Memphis corridor - mean dist: {np.mean(memphis_dists):.1f} km, median: {np.median(memphis_dists):.1f} km")
    print(f"Within 10 km: {sum(1 for d in dists if d <= 10)}/{len(dists)}")
    print(f"Within 25 km: {sum(1 for d in dists if d <= 25)}/{len(dists)}")
    print(f"Within 50 km: {sum(1 for d in dists if d <= 50)}/{len(dists)}")

    # Step 3: Ripley's K
    print("\nComputing Ripley's K-function...")
    along = np.array([p['along_gc_km'] for p in memphis])
    perp = np.array([p['perp_gc_km'] for p in memphis])
    points_2d = np.column_stack([along, perp])

    # Study area: bounding box of the Memphis corridor
    along_range = np.max(along) - np.min(along)
    perp_range = np.max(perp) - np.min(perp)
    # Add 10% buffer
    area = (along_range * 1.1) * (perp_range * 1.1)
    print(f"Study area: {along_range:.1f} x {perp_range:.1f} km = {area:.0f} km²")

    r_values = [0.5, 1, 2, 5, 10, 20, 50]
    k_obs = ripleys_k(points_2d, r_values, area)

    # L(r) = sqrt(K(r)/pi) - r
    l_obs = {r: math.sqrt(max(0, k_obs[r]) / math.pi) - r for r in r_values}

    print("\nObserved Ripley's K and L:")
    for r in r_values:
        print(f"  r={r:5.1f} km: K(r)={k_obs[r]:10.2f}, L(r)={l_obs[r]:8.2f}")

    # Step 4: Directional K
    print("\nComputing directional K-function...")
    k_par, k_perp = directional_k(
        [p['along_gc_km'] for p in memphis],
        [p['perp_gc_km'] for p in memphis],
        r_values, area
    )

    print("\nDirectional K (parallel=along GC, perpendicular=across GC):")
    for r in r_values:
        ratio = k_par.get(r, 0) / k_perp.get(r, 1e-10) if k_perp.get(r, 0) > 0 else float('inf')
        print(f"  r={r:5.1f}: K_par={k_par.get(r,0):10.2f}, K_perp={k_perp.get(r,0):10.2f}, ratio={ratio:.2f}")

    # Step 5: Monte Carlo for K and directional K
    print("\nMonte Carlo (10,000 random point clouds in Memphis corridor)...")
    random.seed(42)
    np.random.seed(42)
    n_trials = 10000
    n_memphis = len(memphis)

    mc_l = {r: [] for r in r_values}
    mc_ratio = {r: [] for r in r_values}

    for trial in range(n_trials):
        if (trial + 1) % 2000 == 0:
            print(f"  Trial {trial+1}/{n_trials}...")

        # Random points in the bounding box
        rand_along = np.random.uniform(np.min(along) - along_range * 0.05,
                                        np.max(along) + along_range * 0.05, n_memphis)
        rand_perp = np.random.uniform(np.min(perp) - perp_range * 0.05,
                                       np.max(perp) + perp_range * 0.05, n_memphis)
        rand_pts = np.column_stack([rand_along, rand_perp])
        k_rand = ripleys_k(rand_pts, r_values, area)

        for r in r_values:
            l_r = math.sqrt(max(0, k_rand[r]) / math.pi) - r
            mc_l[r].append(l_r)

        # Directional K for a subset of trials (expensive)
        if trial < 2000:
            kp, kpp = directional_k(rand_along.tolist(), rand_perp.tolist(), [5, 10, 20], area)
            for r in [5, 10, 20]:
                if kpp.get(r, 0) > 0:
                    mc_ratio[r].append(kp.get(r, 0) / kpp[r])
                else:
                    mc_ratio[r].append(1.0)

    # Compute significance
    print("\n--- Ripley's K Results ---")
    k_results = {}
    for r in r_values:
        mc_arr = np.array(mc_l[r])
        z = (l_obs[r] - np.mean(mc_arr)) / np.std(mc_arr) if np.std(mc_arr) > 0 else 0
        p = np.sum(mc_arr >= l_obs[r]) / n_trials
        lo_95 = np.percentile(mc_arr, 2.5)
        hi_95 = np.percentile(mc_arr, 97.5)
        k_results[str(r)] = {
            'K_observed': round(k_obs[r], 4),
            'L_observed': round(l_obs[r], 4),
            'L_baseline_mean': round(float(np.mean(mc_arr)), 4),
            'L_baseline_std': round(float(np.std(mc_arr)), 4),
            'L_envelope_2.5': round(float(lo_95), 4),
            'L_envelope_97.5': round(float(hi_95), 4),
            'z_score': round(z, 3),
            'p_value': round(p, 4),
        }
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        print(f"  r={r:5.1f}: L={l_obs[r]:8.2f}, baseline={np.mean(mc_arr):8.2f}±{np.std(mc_arr):.2f}, "
              f"z={z:.2f}, p={p:.4f} {sig}")

    # Directional K results
    dir_results = {}
    for r in [5, 10, 20]:
        obs_ratio = k_par.get(r, 0) / k_perp.get(r, 1e-10) if k_perp.get(r, 0) > 0 else float('inf')
        mc_arr = np.array(mc_ratio.get(r, [1.0]))
        z = (obs_ratio - np.mean(mc_arr)) / np.std(mc_arr) if np.std(mc_arr) > 0 else 0
        p = np.sum(mc_arr >= obs_ratio) / len(mc_arr)
        dir_results[str(r)] = {
            'K_parallel': round(k_par.get(r, 0), 4),
            'K_perpendicular': round(k_perp.get(r, 0), 4),
            'ratio_par_perp': round(obs_ratio, 3),
            'baseline_mean_ratio': round(float(np.mean(mc_arr)), 3),
            'z_score': round(z, 3),
            'p_value': round(p, 4),
        }
        print(f"  Directional r={r}: par/perp ratio={obs_ratio:.2f}, baseline={np.mean(mc_arr):.2f}, "
              f"z={z:.2f}, p={p:.4f}")

    # Save results
    save_json(k_results, os.path.join(OUTPUT_DIR, 'ripleys_k_results.json'))
    save_json(dir_results, os.path.join(OUTPUT_DIR, 'directional_k.json'))
    print(f"\nSaved: ripleys_k_results.json, directional_k.json")

    # Step 6: Plots
    print("\nGenerating plots...")

    # L(r) plot with confidence envelope
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    r_plot = np.array(r_values)
    l_plot = np.array([l_obs[r] for r in r_values])
    lo_env = np.array([np.percentile(mc_l[r], 2.5) for r in r_values])
    hi_env = np.array([np.percentile(mc_l[r], 97.5) for r in r_values])
    mean_env = np.array([np.mean(mc_l[r]) for r in r_values])

    ax.fill_between(r_plot, lo_env, hi_env, alpha=0.2, color='gray', label='95% CI (random)')
    ax.plot(r_plot, mean_env, '--', color='gray', label='Expected (random)')
    ax.plot(r_plot, l_plot, 'ro-', markersize=8, linewidth=2, label='Observed (pyramids)')
    ax.axhline(0, color='black', linewidth=0.5, linestyle=':')
    ax.set_xlabel('Distance r (km)')
    ax.set_ylabel('L(r) = sqrt(K(r)/π) - r')
    ax.set_title("Ripley's L-Function: Memphis Corridor Pyramids\n(L > 0 indicates clustering)")
    ax.legend()
    ax.set_xscale('log')
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'ripleys_k_plot.png'), dpi=150)
    plt.close()

    # Scatter plot: along-circle vs perpendicular
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    dynasties = [p.get('dynasty', 0) for p in memphis]
    sc = ax.scatter([p['along_gc_km'] for p in memphis],
                    [p['perp_gc_km'] for p in memphis],
                    c=dynasties, cmap='RdYlBu_r', s=60, edgecolors='black', linewidth=0.5, zorder=3)
    ax.axhline(0, color='red', linewidth=1, linestyle='--', alpha=0.5, label='Great Circle axis')
    ax.set_xlabel('Distance along Great Circle (km)')
    ax.set_ylabel('Distance perpendicular to Great Circle (km)')
    ax.set_title(f'Memphis Corridor Pyramids in Great Circle Coordinates\n({len(memphis)} pyramids)')
    cbar = plt.colorbar(sc, ax=ax, label='Dynasty')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'pyramid_scatter.png'), dpi=150)
    plt.close()

    print(f"Saved: ripleys_k_plot.png, pyramid_scatter.png")
    print("\n=== Analysis B Complete ===")


if __name__ == '__main__':
    main()
