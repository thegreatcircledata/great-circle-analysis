#!/usr/bin/env python3
"""
Bandwidth sensitivity check: repeat the enrichment test at 100, 200, 300, 500 km.
Also run a simple (non-latitude-matched) MC for comparison.
"""

import json
import os
import sys

# Import from main analysis
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oral_tradition_analysis import (
    load_data, TARGET_MOTIFS, compute_enrichment, monte_carlo_enrichment,
    gc_distance, CORRIDOR_WIDTH_KM, N_MONTE_CARLO
)

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    print("Bandwidth Sensitivity Analysis")
    print("=" * 70)

    groups, motif_cols = load_data()
    motif_set = set(motif_cols)

    bandwidths = [100, 200, 300, 500]
    results = {}

    for bw in bandwidths:
        print(f"\n--- Corridor width: {bw} km ---")
        n_on = sum(1 for g in groups if g["gc_dist"] <= bw)
        print(f"  Groups on-corridor: {n_on}/{len(groups)} ({n_on/len(groups)*100:.1f}%)")

        results[bw] = {}
        for cat_name, cat_info in TARGET_MOTIFS.items():
            valid_codes = [c for c in cat_info["codes"] if c in motif_set]
            if not valid_codes:
                continue
            # Quick enrichment (no MC) for speed
            r = compute_enrichment(groups, valid_codes, bw)
            results[bw][cat_name] = {
                "enrichment": r["enrichment_ratio"],
                "on_prev": r["on_corridor_prevalence"],
                "off_prev": r["off_corridor_prevalence"],
                "on_n": r["on_corridor_total"],
            }

        # Print table
        print(f"  {'Category':<30} {'Enrich':>8} {'On-prev':>8} {'Off-prev':>9}")
        print("  " + "-" * 60)
        for cat, r in sorted(results[bw].items(), key=lambda x: -x[1]["enrichment"]):
            print(f"  {cat:<30} {r['enrichment']:>8.3f} {r['on_prev']:>8.3f} {r['off_prev']:>9.3f}")

    # Run MC for 500km bandwidth to check if wider band changes conclusions
    print(f"\n--- Monte Carlo at 500km bandwidth (1000 trials for speed) ---")
    mc_500 = {}
    for cat_name, cat_info in TARGET_MOTIFS.items():
        valid_codes = [c for c in cat_info["codes"] if c in motif_set]
        if not valid_codes:
            continue
        r = monte_carlo_enrichment(groups, valid_codes, n_trials=1000, corridor_km=500)
        mc_500[cat_name] = {
            "enrichment": r["enrichment_ratio"],
            "p_value": r["mc_p_value"],
            "z_score": r["mc_z_score"],
        }
        sig = "***" if r["mc_p_value"] < 0.001 else "**" if r["mc_p_value"] < 0.01 else "*" if r["mc_p_value"] < 0.05 else "ns"
        print(f"  {cat_name:<30} enrich={r['enrichment_ratio']:.3f} p={r['mc_p_value']:.4f} {sig}")

    # Save
    output = {"bandwidths": {str(k): v for k, v in results.items()}, "mc_500km": mc_500}
    path = os.path.join(OUTPUT_DIR, "bandwidth_sensitivity.json")
    with open(path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\n  Written: {path}")


if __name__ == "__main__":
    main()
