#!/usr/bin/env python3
"""
Independent Origins Formal P-Value Test
========================================
Of 100,000 uniformly random great circles, what fraction pass within 200 km
of 4 or more of the 7 independently originated civilizations?

This is the headline statistic for Paper 3.
"""

import numpy as np
import json
import os
import time

# Civilization centers (canonical archaeological capitals)
ORIGINS = {
    "Egypt":       (29.87,  31.22),   # Memphis/Saqqara
    "Mesopotamia": (31.33,  45.64),   # Ur/Uruk region
    "Indus":       (27.32,  68.14),   # Mohenjo-daro
    "China":       (34.27, 108.94),   # Erlitou/Xi'an
    "Mesoamerica": (19.69, -98.84),   # Teotihuacan
    "Andes":       (-13.16, -72.55),  # Cusco/Caral
    "West_Africa": (11.50,  -4.00),   # Jenne-jeno
}

# Alison circle pole
ALISON_POLE = (59.682122, -138.646087)

EARTH_R = 6371.0
N_TRIALS = 100_000
THRESHOLD_KM = 200

def to_xyz(lat_deg, lon_deg):
    """Convert lat/lon in degrees to unit sphere xyz."""
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    return np.array([
        np.cos(lat) * np.cos(lon),
        np.cos(lat) * np.sin(lon),
        np.sin(lat)
    ])

def gc_distance_from_pole(pole_xyz, point_xyz):
    """
    Distance from a point to the great circle defined by a pole.
    The great circle is the set of points 90° from the pole.
    Distance = |angular_distance_to_pole - 90°| * R
    """
    dot = np.clip(np.dot(pole_xyz, point_xyz), -1.0, 1.0)
    angular_dist = np.arccos(dot)  # radians from pole
    # Distance from the great circle = deviation from 90°
    return abs(angular_dist - np.pi/2) * EARTH_R

def gc_distance_from_pole_batch(pole_xyz, points_xyz):
    """Vectorized version for multiple points."""
    dots = np.clip(points_xyz @ pole_xyz, -1.0, 1.0)
    angular_dists = np.arccos(dots)
    return np.abs(angular_dists - np.pi/2) * EARTH_R

def uniform_random_poles(n, rng):
    """Generate n uniformly distributed points on the unit sphere."""
    # Use the standard method: z uniform in [-1,1], phi uniform in [0,2pi]
    z = rng.uniform(-1, 1, n)
    phi = rng.uniform(0, 2*np.pi, n)
    r = np.sqrt(1 - z**2)
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return np.column_stack([x, y, z])

def main():
    print("=" * 60)
    print("INDEPENDENT ORIGINS P-VALUE TEST")
    print(f"N_TRIALS = {N_TRIALS:,}, THRESHOLD = {THRESHOLD_KM} km")
    print("=" * 60)

    # Pre-compute origin XYZ coordinates
    origin_names = list(ORIGINS.keys())
    origin_xyz = np.array([to_xyz(lat, lon) for lat, lon in ORIGINS.values()])
    n_origins = len(origin_names)

    # Compute Alison circle distances
    alison_xyz = to_xyz(*ALISON_POLE)
    alison_dists = gc_distance_from_pole_batch(alison_xyz, origin_xyz)
    alison_within = int(np.sum(alison_dists <= THRESHOLD_KM))

    print(f"\nAlison circle distances to origins:")
    for i, name in enumerate(origin_names):
        flag = "✓" if alison_dists[i] <= THRESHOLD_KM else "✗"
        print(f"  {name:15s}: {alison_dists[i]:8.1f} km {flag}")
    print(f"  Origins within {THRESHOLD_KM} km: {alison_within}")

    # Monte Carlo: uniform random poles
    print(f"\nRunning {N_TRIALS:,} random circles...")
    t0 = time.time()

    rng = np.random.default_rng(42)

    # Process in batches for memory efficiency
    BATCH = 10_000
    counts_ge_threshold = np.zeros(n_origins + 1, dtype=int)  # histogram of n_within
    best_n_within = 0
    best_pole = None
    best_dists = None

    for batch_start in range(0, N_TRIALS, BATCH):
        batch_size = min(BATCH, N_TRIALS - batch_start)
        poles = uniform_random_poles(batch_size, rng)

        # For each pole, compute distance to each origin
        # poles: (batch, 3), origin_xyz: (n_origins, 3)
        # dots: (batch, n_origins)
        dots = np.clip(poles @ origin_xyz.T, -1.0, 1.0)
        angular_dists = np.arccos(dots)
        gc_dists = np.abs(angular_dists - np.pi/2) * EARTH_R  # (batch, n_origins)

        # Count origins within threshold for each pole
        n_within = np.sum(gc_dists <= THRESHOLD_KM, axis=1)  # (batch,)

        # Update histogram
        for k in range(n_origins + 1):
            counts_ge_threshold[k] += int(np.sum(n_within == k))

        # Track best
        batch_max = np.max(n_within)
        if batch_max > best_n_within:
            best_idx = np.argmax(n_within)
            best_n_within = int(batch_max)
            best_pole = poles[best_idx]
            best_dists = gc_dists[best_idx]

    elapsed = time.time() - t0
    print(f"Completed in {elapsed:.1f}s")

    # Compute p-values for each threshold count
    print(f"\nDistribution of n_origins within {THRESHOLD_KM} km:")
    total = N_TRIALS
    cumulative_ge = np.zeros(n_origins + 1, dtype=int)
    for k in range(n_origins, -1, -1):
        cumulative_ge[k] = counts_ge_threshold[k] + (cumulative_ge[k+1] if k < n_origins else 0)

    for k in range(n_origins + 1):
        pct = counts_ge_threshold[k] / total * 100
        cum_pct = cumulative_ge[k] / total * 100
        print(f"  n={k}: {counts_ge_threshold[k]:>7,} ({pct:6.3f}%)  P(≥{k}): {cumulative_ge[k]:>7,} ({cum_pct:6.3f}%)")

    # The key p-value: P(n_within >= alison_within)
    p_value = cumulative_ge[alison_within] / total

    print(f"\n{'='*60}")
    print(f"HEADLINE RESULT:")
    print(f"  Alison circle passes within {THRESHOLD_KM} km of {alison_within}/{n_origins} origins")
    print(f"  P(random circle hits ≥{alison_within} origins) = {p_value:.6f}")
    print(f"  = {cumulative_ge[alison_within]:,} / {total:,}")
    if p_value > 0:
        print(f"  = 1 in {1/p_value:.0f}")
    else:
        print(f"  < 1/{total:,} (zero hits in {total:,} trials)")
    print(f"{'='*60}")

    # Also test at different thresholds
    print(f"\nSensitivity to threshold distance:")
    for thresh_km in [50, 100, 150, 200, 300, 500]:
        rng2 = np.random.default_rng(42)
        n_hits = 0
        alison_n = int(np.sum(alison_dists <= thresh_km))
        if alison_n < 2:
            print(f"  {thresh_km:4d} km: Alison hits {alison_n} origins (skip)")
            continue

        for batch_start in range(0, N_TRIALS, BATCH):
            batch_size = min(BATCH, N_TRIALS - batch_start)
            poles = uniform_random_poles(batch_size, rng2)
            dots = np.clip(poles @ origin_xyz.T, -1.0, 1.0)
            angular_dists = np.arccos(dots)
            gc_dists = np.abs(angular_dists - np.pi/2) * EARTH_R
            n_within = np.sum(gc_dists <= thresh_km, axis=1)
            n_hits += int(np.sum(n_within >= alison_n))

        p = n_hits / total
        print(f"  {thresh_km:4d} km: Alison hits {alison_n}, P(≥{alison_n}) = {p:.6f} ({n_hits:,}/{total:,})")

    # Best random circle details
    if best_pole is not None:
        best_lat = np.degrees(np.arcsin(best_pole[2]))
        best_lon = np.degrees(np.arctan2(best_pole[1], best_pole[0]))
        print(f"\nBest random circle: pole ({best_lat:.2f}°N, {best_lon:.2f}°E)")
        print(f"  Hits {best_n_within} origins within {THRESHOLD_KM} km:")
        for i, name in enumerate(origin_names):
            flag = "✓" if best_dists[i] <= THRESHOLD_KM else "✗"
            print(f"    {name:15s}: {best_dists[i]:8.1f} km {flag}")

    # Save results
    results = {
        "analysis": "Independent Origins Formal P-Value",
        "n_trials": N_TRIALS,
        "threshold_km": THRESHOLD_KM,
        "civilizations": {name: {"lat": ORIGINS[name][0], "lon": ORIGINS[name][1]}
                         for name in origin_names},
        "alison_pole": {"lat": ALISON_POLE[0], "lon": ALISON_POLE[1]},
        "alison_distances_km": {name: float(alison_dists[i]) for i, name in enumerate(origin_names)},
        "alison_origins_within": alison_within,
        "alison_origins_hit": [name for i, name in enumerate(origin_names) if alison_dists[i] <= THRESHOLD_KM],
        "distribution": {str(k): int(counts_ge_threshold[k]) for k in range(n_origins + 1)},
        "cumulative_ge": {str(k): int(cumulative_ge[k]) for k in range(n_origins + 1)},
        "p_value": float(p_value),
        "p_value_description": f"P(random great circle passes within {THRESHOLD_KM}km of >={alison_within} of {n_origins} independent origins)",
        "sensitivity": {},
        "conclusion": ""
    }

    # Add sensitivity results
    for thresh_km in [50, 100, 150, 200, 300, 500]:
        alison_n = int(np.sum(alison_dists <= thresh_km))
        if alison_n >= 2:
            rng3 = np.random.default_rng(42)
            n_hits = 0
            for batch_start in range(0, N_TRIALS, BATCH):
                batch_size = min(BATCH, N_TRIALS - batch_start)
                poles = uniform_random_poles(batch_size, rng3)
                dots = np.clip(poles @ origin_xyz.T, -1.0, 1.0)
                angular_dists = np.arccos(dots)
                gc_dists = np.abs(angular_dists - np.pi/2) * EARTH_R
                n_within = np.sum(gc_dists <= thresh_km, axis=1)
                n_hits += int(np.sum(n_within >= alison_n))
            results["sensitivity"][str(thresh_km)] = {
                "alison_hits": alison_n,
                "p_value": n_hits / total,
                "n_random_hits": n_hits
            }

    # Set conclusion
    if p_value < 0.001:
        results["conclusion"] = f"Highly significant: p = {p_value:.6f}. The probability of a random great circle threading {alison_within}+ of 7 independent origins within {THRESHOLD_KM} km is extremely low."
    elif p_value < 0.05:
        results["conclusion"] = f"Significant at alpha=0.05: p = {p_value:.6f}."
    else:
        results["conclusion"] = f"Not significant: p = {p_value:.6f}. Random great circles frequently hit {alison_within}+ origins within {THRESHOLD_KM} km."

    out_path = os.path.expanduser("~/megalith_site_research/paper3/results/independent_origins_pvalue.json")
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")

if __name__ == "__main__":
    main()
