#!/usr/bin/env python3
"""
Part A: φ Ratio Multiple Testing Correction

Tests all C(17,3) = 680 triplets from Alison's 17 great-circle sites
for golden-ratio distance proportions, then applies Bonferroni correction.
"""

import json
import math
import numpy as np
from itertools import combinations
from pathlib import Path

# ---------------------------------------------------------------------------
# Alison's 17 sites
# ---------------------------------------------------------------------------
ALISON_SITES = [
    {"name": "Giza", "lat": 29.983, "lon": 31.150},
    {"name": "Siwa", "lat": 29.233, "lon": 25.517},
    {"name": "Tassili n'Ajjer", "lat": 26.533, "lon": 9.833},
    {"name": "Paratoari", "lat": -12.800, "lon": -71.417},
    {"name": "Ollantaytambo", "lat": -13.250, "lon": -72.267},
    {"name": "Machu Picchu", "lat": -13.100, "lon": -72.583},
    {"name": "Nazca", "lat": -14.700, "lon": -75.100},
    {"name": "Easter Island", "lat": -27.100, "lon": -109.333},
    {"name": "Aneityum Island", "lat": -20.167, "lon": 169.800},
    {"name": "Preah Vihear", "lat": 14.400, "lon": 104.667},
    {"name": "Sukhothai", "lat": 17.017, "lon": 99.700},
    {"name": "Pyay", "lat": 19.250, "lon": 95.083},
    {"name": "Khajuraho", "lat": 24.850, "lon": 79.933},
    {"name": "Mohenjo Daro", "lat": 27.250, "lon": 68.283},
    {"name": "Persepolis", "lat": 29.933, "lon": 52.917},
    {"name": "Ur", "lat": 30.950, "lon": 46.117},
    {"name": "Petra", "lat": 30.317, "lon": 35.467},
]

EARTH_RADIUS_KM = 6371.0
PHI = (1 + math.sqrt(5)) / 2


def great_circle_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_KM * math.asin(min(1.0, math.sqrt(a)))


def monte_carlo_phi_pvalue(tolerance_frac, n_trials=200_000, rng=None):
    """p-value: fraction of random 3-point great-circle divisions with
    any pair of gaps having ratio within tolerance_frac of φ."""
    if rng is None:
        rng = np.random.default_rng(42)
    hits = 0
    for _ in range(n_trials):
        positions = np.sort(rng.uniform(0, 360, 3))
        gaps = np.diff(np.append(positions, positions[0] + 360))
        found = False
        for i in range(3):
            if found:
                break
            for j in range(3):
                if i != j and gaps[j] > 0:
                    ratio = gaps[i] / gaps[j]
                    if abs(ratio - PHI) / PHI <= tolerance_frac:
                        found = True
                        break
        if found:
            hits += 1
    return hits / n_trials


def main():
    out_dir = Path("outputs/extended_analysis/followup")
    out_dir.mkdir(parents=True, exist_ok=True)

    n_sites = len(ALISON_SITES)
    n_triplets = math.comb(n_sites, 3)
    print(f"Testing all C({n_sites},3) = {n_triplets} triplets for φ ratio\n")

    # Precompute pairwise distances
    dist_cache = {}
    for i, s1 in enumerate(ALISON_SITES):
        for j, s2 in enumerate(ALISON_SITES):
            if i < j:
                d = great_circle_km(s1["lat"], s1["lon"], s2["lat"], s2["lon"])
                dist_cache[(i, j)] = d

    def get_dist(i, j):
        if i == j:
            return 0
        key = (min(i, j), max(i, j))
        return dist_cache[key]

    # Test all triplets
    results = []
    site_indices = list(range(n_sites))

    for i, j, k in combinations(site_indices, 3):
        d12 = get_dist(i, j)
        d23 = get_dist(j, k)
        d13 = get_dist(i, k)
        distances = sorted([d12, d23, d13])

        # Find best φ-ratio match among all 6 ordered pairs
        best_error = float("inf")
        best_pair = None
        for a in range(3):
            for b in range(3):
                if a != b and distances[b] > 0:
                    ratio = distances[a] / distances[b]
                    error = abs(ratio - PHI) / PHI
                    if error < best_error:
                        best_error = error
                        best_pair = (a, b)
                        best_ratio = ratio

        names = (
            ALISON_SITES[i]["name"],
            ALISON_SITES[j]["name"],
            ALISON_SITES[k]["name"],
        )

        results.append({
            "sites": names,
            "best_phi_error_frac": best_error,
            "best_phi_error_pct": best_error * 100,
            "best_ratio": round(best_ratio, 6),
            "distances_km": [round(d, 1) for d in sorted([d12, d23, d13])],
        })

    # Sort by best phi error
    results.sort(key=lambda x: x["best_phi_error_frac"])

    # Reference: Angkor-Giza-Nazca error from Directive 1
    # Find it in results (Angkor is not in Alison's list by that name -
    # but Preah Vihear is nearby; the original test used Angkor Wat directly)
    # The actual test used Angkor Wat (13.4125, 103.867) which isn't in this list.
    # Giza-Nazca is present. Let's find the best Giza-Nazca triplet.
    giza_nazca_triplets = [
        r for r in results
        if "Giza" in r["sites"] and "Nazca" in r["sites"]
    ]

    # The Directive 1 error was 0.2895%
    d1_error = 0.002895
    threshold = d1_error  # 0.29% relative error

    n_hits_at_threshold = sum(1 for r in results if r["best_phi_error_frac"] <= threshold)

    # Compute individual p-values for top 20 using MC
    print("Computing Monte Carlo p-values for top 20 triplets...")
    rng = np.random.default_rng(42)
    top20 = results[:20]
    for r in top20:
        r["p_value"] = monte_carlo_phi_pvalue(r["best_phi_error_frac"], n_trials=200_000, rng=rng)
        print(f"  {r['sites']}: error={r['best_phi_error_pct']:.4f}%, p={r['p_value']:.6f}")

    # Bonferroni correction
    best_p = top20[0]["p_value"]
    bonferroni_p = min(1.0, best_p * n_triplets)

    # Also compute for the original Directive 1 result
    d1_p = 0.0083  # from Directive 1
    d1_bonferroni = min(1.0, d1_p * n_triplets)

    # Find rank of best Giza-Nazca triplet
    giza_nazca_rank = None
    for idx, r in enumerate(results):
        if "Giza" in r["sites"] and "Nazca" in r["sites"]:
            giza_nazca_rank = idx + 1
            giza_nazca_best = r
            break

    print(f"\n{'=' * 70}")
    print(f"PART A: φ MULTIPLE TESTING RESULTS")
    print(f"{'=' * 70}")
    print(f"\nTotal triplets tested: {n_triplets}")
    print(f"Triplets with φ error ≤ 0.29% (Directive 1 threshold): {n_hits_at_threshold}")
    print(f"\nTop 5 φ-ratio triplets:")
    for i, r in enumerate(results[:5]):
        p_str = f"p={r.get('p_value', 'N/A')}" if 'p_value' in r else ""
        print(f"  #{i+1}: {r['sites']} — error={r['best_phi_error_pct']:.4f}%, ratio={r['best_ratio']}, {p_str}")

    print(f"\nBest Giza-Nazca triplet: rank #{giza_nazca_rank}")
    if giza_nazca_best:
        print(f"  {giza_nazca_best['sites']} — error={giza_nazca_best['best_phi_error_pct']:.4f}%")

    print(f"\nBonferroni correction:")
    print(f"  Best triplet p = {best_p:.6f}, corrected = {bonferroni_p:.4f}")
    print(f"  Directive 1 (Angkor-Giza-Nazca) p = {d1_p}, corrected = {d1_bonferroni:.4f}")

    if bonferroni_p < 0.05:
        verdict = "SURVIVES CORRECTION — best triplet is significant even after testing all 680."
    else:
        verdict = "DOES NOT SURVIVE — expected under multiple testing."

    if n_hits_at_threshold > 3:
        verdict += f" Multiple triplets ({n_hits_at_threshold}) match as well or better than Angkor-Giza-Nazca."
    elif n_hits_at_threshold == 1:
        verdict += " Angkor-Giza-Nazca is uniquely close."

    print(f"\nVERDICT: {verdict}")
    print(f"{'=' * 70}")

    # Save outputs
    # All 680 triplets
    all_results_serializable = []
    for r in results:
        entry = dict(r)
        entry["sites"] = list(entry["sites"])
        all_results_serializable.append(entry)

    with open(out_dir / "phi_all_triplets.json", "w") as f:
        json.dump({
            "description": "All C(17,3)=680 triplets from Alison's sites, sorted by φ-ratio error",
            "n_sites": n_sites,
            "n_triplets": n_triplets,
            "threshold_error_frac": threshold,
            "n_hits_at_threshold": n_hits_at_threshold,
            "giza_nazca_best_rank": giza_nazca_rank,
            "bonferroni_p_best": bonferroni_p,
            "bonferroni_p_directive1": d1_bonferroni,
            "verdict": verdict,
            "triplets": all_results_serializable,
        }, f, indent=2)

    # Top 20
    top20_serializable = []
    for r in top20:
        entry = dict(r)
        entry["sites"] = list(entry["sites"])
        top20_serializable.append(entry)

    with open(out_dir / "phi_top20.json", "w") as f:
        json.dump({
            "description": "Top 20 φ-ratio triplets from Alison's 17 sites",
            "bonferroni_p": bonferroni_p,
            "verdict": verdict,
            "triplets": top20_serializable,
        }, f, indent=2)

    print(f"\nSaved: {out_dir / 'phi_all_triplets.json'}")
    print(f"Saved: {out_dir / 'phi_top20.json'}")


if __name__ == "__main__":
    main()
