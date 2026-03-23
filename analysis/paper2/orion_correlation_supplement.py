#!/usr/bin/env python3
"""
Supplement to Orion Correlation Test:
1. Fix the -5000 epoch anomaly (precession pole crossing)
2. Add collinearity-matched random triangle test
3. Better characterize what's actually varying
"""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# Import from main script
from orion_correlation import (
    pyramid_triangle_km, belt_positions_at_epoch,
    procrustes_distance, procrustes_distance_with_mirror,
    OUT, FIG, PYRAMIDS, BELT_STARS
)


def collinearity_ratio(pts):
    """
    Measure how collinear a triangle is.
    Returns ratio of smallest to largest singular value of centered points.
    0 = perfectly collinear, 1 = equilateral-like.
    """
    centered = pts - pts.mean(axis=0)
    _, s, _ = np.linalg.svd(centered)
    return s[1] / s[0] if s[0] > 0 else 0


def generate_collinear_triangle(target_ratio):
    """
    Generate a random near-collinear triangle with similar collinearity
    to the target ratio, by starting with 3 points on a line and
    perturbing perpendicular to it.
    """
    # Three points along a line
    t = np.sort(np.random.uniform(-1, 1, 3))
    angle = np.random.uniform(0, 2*np.pi)
    direction = np.array([np.cos(angle), np.sin(angle)])
    perp = np.array([-np.sin(angle), np.cos(angle)])

    pts = np.outer(t, direction)

    # Add perpendicular perturbation scaled to match target collinearity
    perp_scale = target_ratio * np.std(t)
    perp_offsets = np.random.randn(3) * perp_scale
    pts += np.outer(perp_offsets, perp)

    pts -= pts.mean(axis=0)
    return pts


def collinearity_matched_significance(n_trials=10000):
    """
    Monte Carlo with near-collinear triangles matching the pyramids' aspect ratio.
    """
    pyr = pyramid_triangle_km()
    pyr_collinearity = collinearity_ratio(pyr)
    print(f"  Pyramid collinearity ratio: {pyr_collinearity:.6f}")

    # Pre-compute Belt positions — skip -5000 epoch (anomalous)
    test_epochs = [y for y in range(-15000, 2501, 500) if y != -5000]
    belt_shapes = {y: belt_positions_at_epoch(y) for y in test_epochs}

    # Get pyramid best fit
    pyr_best = min(procrustes_distance(pyr, belt_shapes[y]) for y in test_epochs)
    print(f"  Pyramid best Procrustes distance: {pyr_best:.6f}")

    # For the Belt stars, what's their collinearity?
    belt_col = collinearity_ratio(belt_shapes[-10500])
    print(f"  Belt stars collinearity ratio at -10500: {belt_col:.6f}")

    np.random.seed(42)
    count_better = 0
    random_dists = []
    random_cols = []

    for i in range(n_trials):
        tri = generate_collinear_triangle(pyr_collinearity)
        best_d = min(procrustes_distance(tri, belt_shapes[y]) for y in test_epochs)
        random_dists.append(best_d)
        random_cols.append(collinearity_ratio(tri))

        if best_d <= pyr_best:
            count_better += 1

        if (i + 1) % 2000 == 0:
            print(f"  Monte Carlo: {i+1}/{n_trials}")

    p_value = count_better / n_trials

    results = {
        "description": "Significance test using near-collinear random triangles matching pyramid aspect ratio",
        "n_trials": n_trials,
        "pyramid_collinearity_ratio": round(float(pyr_collinearity), 6),
        "belt_collinearity_ratio_10500BC": round(float(belt_col), 6),
        "pyramid_best_procrustes": round(float(pyr_best), 6),
        "random_mean_collinearity": round(float(np.mean(random_cols)), 6),
        "random_mean_best_procrustes": round(float(np.mean(random_dists)), 6),
        "random_std_best_procrustes": round(float(np.std(random_dists)), 6),
        "random_median_best_procrustes": round(float(np.median(random_dists)), 6),
        "random_min_best_procrustes": round(float(np.min(random_dists)), 6),
        "p_value": round(p_value, 4),
        "count_better_or_equal": count_better,
    }

    print(f"\n  Collinearity-matched results:")
    print(f"  Random mean best-fit: {np.mean(random_dists):.6f} ± {np.std(random_dists):.6f}")
    print(f"  Random minimum: {np.min(random_dists):.6f}")
    print(f"  p-value: {p_value:.4f} ({count_better}/{n_trials} better)")

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(random_dists, bins=60, alpha=0.7, color='steelblue', edgecolor='black', linewidth=0.5)
    ax.axvline(pyr_best, color='red', linewidth=2.5, label=f'Pyramids: {pyr_best:.4f}')
    ax.set_xlabel('Best Procrustes Distance (across all epochs)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'Collinearity-Matched Significance Test\n'
                 f'{n_trials:,} random near-collinear triangles — p = {p_value}', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "significance_collinear_matched.png", dpi=150)
    plt.close(fig)
    print(f"  Saved collinearity-matched significance plot")

    return results


def investigate_variation():
    """
    Understand what's actually changing: the triangle shape barely varies
    with precession since the stars are close together on the sky.
    Compute the actual angular separations and triangle angles at key epochs.
    """
    pyr = pyramid_triangle_km()
    pyr_sides = np.array([
        np.linalg.norm(pyr[1] - pyr[0]),
        np.linalg.norm(pyr[2] - pyr[1]),
        np.linalg.norm(pyr[2] - pyr[0]),
    ])
    pyr_ratios = pyr_sides / pyr_sides.max()

    results = {
        "pyramid_side_ratios": pyr_ratios.tolist(),
        "pyramid_collinearity": float(collinearity_ratio(pyr)),
        "epochs": {}
    }

    for year in [-15000, -12500, -10500, -7500, -5000, -2550, 0, 2026]:
        try:
            stars = belt_positions_at_epoch(year)
            sides = np.array([
                np.linalg.norm(stars[1] - stars[0]),
                np.linalg.norm(stars[2] - stars[1]),
                np.linalg.norm(stars[2] - stars[0]),
            ])
            ratios = sides / sides.max()
            col = collinearity_ratio(stars)
            d = procrustes_distance(pyr, stars)
            results["epochs"][str(year)] = {
                "side_ratios": ratios.tolist(),
                "collinearity": float(col),
                "procrustes_to_pyramids": float(d),
            }
        except Exception as e:
            results["epochs"][str(year)] = {"error": str(e)}

    return results


def main():
    print("=" * 60)
    print("ORION CORRELATION — SUPPLEMENTARY ANALYSIS")
    print("=" * 60)

    # 1. Check the -5000 epoch anomaly
    print("\n--- Investigating -5000 epoch anomaly ---")
    for y in [-5500, -5250, -5000, -4750, -4500]:
        try:
            stars = belt_positions_at_epoch(y)
            col = collinearity_ratio(stars)
            print(f"  {y}: collinearity={col:.6f}, points={stars}")
        except Exception as e:
            print(f"  {y}: ERROR — {e}")

    # 2. Triangle variation analysis
    print("\n--- Triangle variation across epochs ---")
    var = investigate_variation()
    print(f"  Pyramid side ratios: {[f'{r:.4f}' for r in var['pyramid_side_ratios']]}")
    for year, data in var["epochs"].items():
        if "error" not in data:
            print(f"  {year}: ratios={[f'{r:.4f}' for r in data['side_ratios']]}, "
                  f"col={data['collinearity']:.4f}, procrustes={data['procrustes_to_pyramids']:.6f}")

    # 3. Collinearity-matched significance test
    print("\n--- Collinearity-matched Monte Carlo ---")
    sig = collinearity_matched_significance(n_trials=10000)

    # Save
    with open(OUT / "significance_collinear_matched.json", 'w') as f:
        json.dump(sig, f, indent=2)
    print(f"\n  Saved {OUT / 'significance_collinear_matched.json'}")

    with open(OUT / "triangle_variation.json", 'w') as f:
        json.dump(var, f, indent=2)
    print(f"  Saved {OUT / 'triangle_variation.json'}")

    # 4. Summary
    print("\n" + "=" * 60)
    print("SUPPLEMENTARY FINDINGS")
    print("=" * 60)
    print(f"Pyramid collinearity ratio: {var['pyramid_collinearity']:.4f}")
    print(f"Belt star collinearity at -10500: {sig['belt_collinearity_ratio_10500BC']:.4f}")
    print(f"\nWith fully random triangles: p = 0.0000 (trivially significant)")
    print(f"With collinearity-matched triangles: p = {sig['p_value']:.4f}")

    if sig['p_value'] > 0.05:
        print("\nKEY FINDING: When comparing against triangles of similar elongation,")
        print("the pyramids' match to Orion is NOT significant. The apparent match")
        print("is primarily due to both patterns being near-collinear, not a specific")
        print("shape correspondence.")
    elif sig['p_value'] <= 0.05:
        print(f"\nKEY FINDING: The match remains significant (p={sig['p_value']:.4f})")
        print("even against collinearity-matched random triangles.")
        print("This suggests a genuine shape correspondence beyond mere collinearity.")


if __name__ == "__main__":
    main()
