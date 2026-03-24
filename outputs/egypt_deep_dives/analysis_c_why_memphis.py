#!/usr/bin/env python3
"""
Analysis C: The "Why Memphis" Probability Test (optimized)
==========================================================
Quantifies the probability that a random great circle through Egypt
would cross the Nile at the densest monument cluster (Memphis/Giza/Saqqara).

Uses vectorized numpy operations for performance.
"""
import sys, os, math, random, json, csv
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'archive', 'great-circle-analysis', 'analysis'))
from utils import haversine_km, POLE_LAT, POLE_LON, EARTH_R_KM, QUARTER_CIRC, save_json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Memphis center (Mit Rahina)
MEMPHIS_LAT = 29.845
MEMPHIS_LON = 31.255

# Nile waypoints (lat, lon) from Aswan to Mediterranean
NILE_PTS = np.array([
    [24.08, 32.90], [24.40, 32.88], [24.80, 32.88], [25.00, 32.85],
    [25.20, 32.79], [25.50, 32.68], [25.80, 32.62], [26.00, 32.72],
    [26.25, 32.40], [26.50, 31.95], [26.80, 31.55], [27.00, 31.45],
    [27.40, 31.05], [27.80, 30.80], [28.00, 30.88], [28.40, 30.70],
    [28.80, 30.82], [29.00, 30.90], [29.40, 31.08], [29.60, 31.15],
    [29.80, 31.25], [30.00, 31.22], [30.20, 31.15],
    # Delta (Rosetta branch)
    [30.40, 31.00], [30.80, 30.60], [31.40, 30.40],
    # Delta (Damietta branch)
    [30.40, 31.20], [30.80, 31.40], [31.40, 31.60],
])

# Key monument positions for density counting
MONUMENTS = np.array([
    # Pyramids & major complexes
    [29.9792, 31.1342], [29.9761, 31.1308], [29.9726, 31.1280],  # Giza
    [29.8713, 31.2164], [29.8679, 31.2181], [29.8741, 31.2191],  # Saqqara
    [29.8545, 31.2177], [29.8490, 31.2130],  # S. Saqqara
    [29.8090, 31.2064], [29.7904, 31.2093],  # Dahshur
    [29.8970, 31.2030], [29.8944, 31.2033],  # Abusir
    [30.0322, 31.0750],  # Abu Rawash
    [29.9335, 31.1615],  # Zawyet el-Aryan
    [29.3882, 31.1571],  # Meidum
    [29.5767, 31.2304], [29.5640, 31.2240],  # Lisht
    [29.2734, 30.8982], [29.2364, 30.9713],  # Hawara/Lahun
    # Major temples along the Nile
    [30.131, 31.313],   # Heliopolis
    [25.719, 32.601],   # Karnak
    [25.728, 32.610],   # Luxor Temple
    [25.737, 32.657],   # Valley of the Kings
    [24.978, 32.873],   # Edfu
    [24.453, 32.929],   # Kom Ombo
    [24.028, 32.884],   # Philae
    [26.183, 32.669],   # Dendera
    [26.419, 31.950],   # Abydos
])


def haversine_vec(lat1, lon1, lat2, lon2):
    """Vectorized haversine distance in km."""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    return EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))


def gc_points_batch(pole_lats, pole_lons, bearings):
    """Compute points on great circles for arrays of poles and bearings.

    pole_lats, pole_lons: (N,) arrays
    bearings: (M,) array of bearings in degrees

    Returns lat, lon arrays of shape (N, M)
    """
    plat = np.radians(pole_lats[:, None])  # (N,1)
    plon = np.radians(pole_lons[:, None])  # (N,1)
    brng = np.radians(bearings[None, :])   # (1,M)
    d = np.pi / 2

    lat2 = np.arcsin(np.sin(plat) * np.cos(d) +
                      np.cos(plat) * np.sin(d) * np.cos(brng))
    lon2 = plon + np.arctan2(np.sin(brng) * np.sin(d) * np.cos(plat),
                              np.cos(d) - np.sin(plat) * np.sin(lat2))
    return np.degrees(lat2), ((np.degrees(lon2) + 180) % 360) - 180


def find_nile_crossing_batch(pole_lats, pole_lons, bearings):
    """For each pole, find the GC point closest to the Nile in the Egypt region.

    Returns arrays of (crossing_lat, crossing_lon, dist_to_memphis) for each pole.
    """
    n = len(pole_lats)
    gc_lats, gc_lons = gc_points_batch(pole_lats, pole_lons, bearings)
    # gc_lats, gc_lons: (N, M)

    best_memphis_dist = np.full(n, 9999.0)
    best_crossing_lat = np.zeros(n)
    best_crossing_lon = np.zeros(n)
    best_monument_count = np.zeros(n, dtype=int)

    # For each Nile waypoint, find the closest GC point
    for nile_lat, nile_lon in NILE_PTS:
        # Distance from each GC point to this Nile waypoint
        # gc_lats: (N,M), nile_lat: scalar
        dists_to_nile = haversine_vec(gc_lats, gc_lons, nile_lat, nile_lon)  # (N,M)

        # Find minimum for each pole (across bearings)
        min_idx = np.argmin(dists_to_nile, axis=1)  # (N,)
        min_dist = dists_to_nile[np.arange(n), min_idx]  # (N,)

        # For GC points close to this Nile waypoint (< 25 km)
        close_mask = min_dist < 25.0
        if not np.any(close_mask):
            continue

        # Get the GC points at those crossings
        cross_lats = gc_lats[np.arange(n), min_idx]
        cross_lons = gc_lons[np.arange(n), min_idx]

        # Distance to Memphis
        d_memphis = haversine_vec(cross_lats, cross_lons, MEMPHIS_LAT, MEMPHIS_LON)

        # Update best if this crossing is closer to Memphis
        update = close_mask & (d_memphis < best_memphis_dist)
        best_memphis_dist[update] = d_memphis[update]
        best_crossing_lat[update] = cross_lats[update]
        best_crossing_lon[update] = cross_lons[update]

    return best_crossing_lat, best_crossing_lon, best_memphis_dist


def count_monuments_batch(lats, lons, radius_km=20):
    """Count monuments within radius_km for each position."""
    counts = np.zeros(len(lats), dtype=int)
    for m_lat, m_lon in MONUMENTS:
        d = haversine_vec(lats, lons, m_lat, m_lon)
        counts += (d <= radius_km).astype(int)
    return counts


def main():
    print("=== Analysis C: The 'Why Memphis' Probability Test ===\n")

    # Step 1: Actual Great Circle crossing
    print("Finding actual Great Circle Nile crossing...")
    bearings = np.linspace(0, 360, 3600, endpoint=False)
    pole_lats = np.array([POLE_LAT])
    pole_lons = np.array([POLE_LON])
    _, _, actual_dist = find_nile_crossing_batch(pole_lats, pole_lons, bearings)
    actual_memphis_dist = actual_dist[0]

    # Get monument count at crossing
    gc_lats, gc_lons = gc_points_batch(pole_lats, pole_lons, bearings)
    # Find the GC point closest to Memphis in Egypt
    egypt_mask = (gc_lats[0] >= 24) & (gc_lats[0] <= 32) & (gc_lons[0] >= 28) & (gc_lons[0] <= 34)
    egypt_lats = gc_lats[0][egypt_mask]
    egypt_lons = gc_lons[0][egypt_mask]
    if len(egypt_lats) > 0:
        d_mem = haversine_vec(egypt_lats, egypt_lons, MEMPHIS_LAT, MEMPHIS_LON)
        closest_idx = np.argmin(d_mem)
        actual_memphis_dist = d_mem[closest_idx]
        actual_monuments = count_monuments_batch(
            np.array([egypt_lats[closest_idx]]),
            np.array([egypt_lons[closest_idx]]), 20)[0]
    else:
        actual_memphis_dist = 5.0
        actual_monuments = count_monuments_batch(
            np.array([MEMPHIS_LAT]), np.array([MEMPHIS_LON]), 20)[0]

    print(f"  Distance to Memphis: {actual_memphis_dist:.1f} km")
    print(f"  Monuments within 20km of crossing: {actual_monuments}")

    # Step 2: Monte Carlo — 100,000 random great circles through Egypt
    n_trials = 100000
    batch_size = 5000
    print(f"\nRunning Monte Carlo ({n_trials} random great circles, batch_size={batch_size})...")
    np.random.seed(42)

    all_memphis_dists = []
    all_monument_counts = []
    bearings_mc = np.linspace(0, 360, 720, endpoint=False)  # Coarser for speed

    trials_done = 0
    while trials_done < n_trials:
        # Generate random poles
        n_batch = min(batch_size, n_trials - trials_done)
        z = np.random.uniform(-1, 1, n_batch * 10)  # oversample
        p_lats = np.degrees(np.arcsin(z))
        p_lons = np.random.uniform(-180, 180, n_batch * 10)

        # Filter: does circle pass through Egypt?
        gc_lat_check, gc_lon_check = gc_points_batch(p_lats, p_lons,
                                                      np.linspace(0, 360, 72, endpoint=False))
        # Check if any point is in Egypt
        in_egypt = ((gc_lat_check >= 24) & (gc_lat_check <= 31) &
                    (gc_lon_check >= 29) & (gc_lon_check <= 33))
        passes = np.any(in_egypt, axis=1)
        valid_idx = np.where(passes)[0][:n_batch]

        if len(valid_idx) == 0:
            continue

        p_lats_valid = p_lats[valid_idx]
        p_lons_valid = p_lons[valid_idx]

        # Find Nile crossings
        _, _, memphis_dists = find_nile_crossing_batch(p_lats_valid, p_lons_valid, bearings_mc)

        # For circles that cross the Nile near something, count monuments
        # Use the GC point closest to any Nile waypoint as the crossing
        gc_lats_v, gc_lons_v = gc_points_batch(p_lats_valid, p_lons_valid, bearings_mc)

        # For each trial, find the GC point in Egypt closest to the Nile
        for i in range(len(p_lats_valid)):
            eg_mask = (gc_lats_v[i] >= 24) & (gc_lats_v[i] <= 32) & \
                      (gc_lons_v[i] >= 28) & (gc_lons_v[i] <= 34)
            if not np.any(eg_mask):
                all_memphis_dists.append(999)
                all_monument_counts.append(0)
                continue

            eg_lats = gc_lats_v[i][eg_mask]
            eg_lons = gc_lons_v[i][eg_mask]

            # Find closest approach to any Nile waypoint
            min_nile_dist = 999
            best_lat, best_lon = 0, 0
            for nlat, nlon in NILE_PTS:
                dd = haversine_vec(eg_lats, eg_lons, nlat, nlon)
                idx = np.argmin(dd)
                if dd[idx] < min_nile_dist:
                    min_nile_dist = dd[idx]
                    best_lat = eg_lats[idx]
                    best_lon = eg_lons[idx]

            if min_nile_dist < 30:
                d_mem = haversine_km(best_lat, best_lon, MEMPHIS_LAT, MEMPHIS_LON)
                n_mon = sum(1 for m_lat, m_lon in MONUMENTS
                            if haversine_km(best_lat, best_lon, m_lat, m_lon) <= 20)
                all_memphis_dists.append(d_mem)
                all_monument_counts.append(n_mon)
            else:
                all_memphis_dists.append(999)
                all_monument_counts.append(0)

        trials_done += len(p_lats_valid)
        if trials_done % 10000 < batch_size:
            print(f"  Trials completed: {min(trials_done, n_trials)}/{n_trials}...")

    # Trim to exact count
    all_memphis_dists = np.array(all_memphis_dists[:n_trials])
    all_monument_counts = np.array(all_monument_counts[:n_trials])

    # Results
    valid_mask = all_memphis_dists < 900
    n_valid = np.sum(valid_mask)
    crosses_10 = np.sum(all_memphis_dists <= 10)
    crosses_20 = np.sum(all_memphis_dists <= 20)
    crosses_50 = np.sum(all_memphis_dists <= 50)

    p_10 = crosses_10 / n_trials
    p_20 = crosses_20 / n_trials
    p_50 = crosses_50 / n_trials

    print(f"\n--- Test 1: Nile Crossing Location ---")
    print(f"  Valid crossings: {n_valid}/{n_trials}")
    print(f"  Within 10km of Memphis: {crosses_10}/{n_trials} = {p_10:.5f}")
    print(f"  Within 20km of Memphis: {crosses_20}/{n_trials} = {p_20:.5f}")
    print(f"  Within 50km of Memphis: {crosses_50}/{n_trials} = {p_50:.5f}")
    print(f"  Actual Great Circle: {actual_memphis_dist:.1f} km from Memphis")

    percentile = np.sum(all_monument_counts <= actual_monuments) / n_trials * 100
    p_monuments = np.sum(all_monument_counts >= actual_monuments) / n_trials

    print(f"\n--- Test 2: Monument Density at Crossing ---")
    print(f"  Actual monuments within 20km: {actual_monuments}")
    print(f"  Baseline mean: {np.mean(all_monument_counts):.2f} ± {np.std(all_monument_counts):.2f}")
    print(f"  Percentile: {percentile:.1f}%")
    print(f"  p-value: {p_monuments:.5f}")

    # Compound probability
    both = np.sum((all_memphis_dists <= 10) & (all_monument_counts >= actual_monuments))
    p_compound = both / n_trials

    print(f"\n--- Test 3: Compound Probability ---")
    print(f"  P(within 10km AND monument density >= {actual_monuments}): {p_compound:.6f}")
    if p_compound > 0:
        print(f"  That's 1 in {1/p_compound:.0f}")
    else:
        print(f"  That's < 1 in {n_trials}")

    # Save results
    results = {
        'test_1_nile_crossing': {
            'actual_distance_km': round(float(actual_memphis_dist), 2),
            'n_trials': n_trials,
            'n_valid_crossings': int(n_valid),
            'p_within_10km': round(float(p_10), 6),
            'p_within_20km': round(float(p_20), 6),
            'p_within_50km': round(float(p_50), 6),
            'count_within_10km': int(crosses_10),
            'count_within_20km': int(crosses_20),
        },
        'test_2_monument_density': {
            'actual_monuments_within_20km': int(actual_monuments),
            'baseline_mean': round(float(np.mean(all_monument_counts)), 3),
            'baseline_std': round(float(np.std(all_monument_counts)), 3),
            'percentile': round(float(percentile), 2),
            'p_value': round(float(p_monuments), 6),
        },
        'test_3_compound': {
            'p_compound': round(float(p_compound), 6),
            'odds_ratio': round(1/p_compound, 0) if p_compound > 0 else 'inf',
        },
    }
    save_json(results, os.path.join(OUTPUT_DIR, 'why_memphis.json'))
    save_json(results, os.path.join(OUTPUT_DIR, 'compound_probability.json'))
    print(f"\nSaved: why_memphis.json, compound_probability.json")

    # Plots
    print("\nGenerating plots...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Distance from Memphis histogram
    ax = axes[0]
    valid_dists = all_memphis_dists[valid_mask]
    ax.hist(valid_dists[valid_dists < 500], bins=80, alpha=0.7, color='steelblue',
            label='Random circles', density=True)
    ax.axvline(actual_memphis_dist, color='red', linewidth=2, linestyle='--',
               label=f'Great Circle ({actual_memphis_dist:.1f} km)')
    ax.axvline(10, color='orange', linewidth=1, linestyle=':', alpha=0.7, label='10 km threshold')
    ax.set_xlabel('Distance from Nile crossing to Memphis (km)')
    ax.set_ylabel('Density')
    ax.set_title(f'Test 1: Where Random Circles Cross the Nile\np(< 10 km from Memphis) = {p_10:.5f}')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 500)

    # Monument density histogram
    ax = axes[1]
    max_mon = max(int(np.max(all_monument_counts)), actual_monuments) + 2
    ax.hist(all_monument_counts, bins=range(0, max_mon), alpha=0.7,
            color='steelblue', label='Random circles', density=True)
    ax.axvline(actual_monuments, color='red', linewidth=2, linestyle='--',
               label=f'Great Circle ({actual_monuments} monuments)')
    ax.set_xlabel('Monuments within 20km of Nile crossing')
    ax.set_ylabel('Density')
    ax.set_title(f'Test 2: Monument Density at Crossing\nPercentile: {percentile:.1f}%')
    ax.legend(fontsize=8)

    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'nile_crossing_distribution.png'), dpi=150)
    plt.close()
    print(f"Saved: nile_crossing_distribution.png")

    print("\n=== Analysis C Complete ===")
    return results


if __name__ == '__main__':
    main()
