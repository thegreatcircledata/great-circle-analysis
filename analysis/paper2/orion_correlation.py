#!/usr/bin/env python3
"""
Orion Correlation Formal Fit Test
=================================
Tests the Bauval/Hancock claim that the three Giza pyramids mirror
Orion's Belt stars as they appeared at ~10,500 BC.

Computes:
  1. Procrustes distance (shape similarity) across epochs
  2. Mirror vs non-mirror fit
  3. Significance via random triangle Monte Carlo
  4. Best-fit epoch identification
"""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from astropy.coordinates import SkyCoord, FK5, AltAz, EarthLocation
from astropy.time import Time
from astropy import units as u
from scipy.spatial import procrustes as scipy_procrustes

OUT = Path("outputs/orion_correlation")
FIG = OUT / "figures"

# ============================================================
# 1. DATA
# ============================================================

# Pyramid positions (WGS84)
PYRAMIDS = {
    "Khufu":   (29.9792, 31.1342),
    "Khafre":  (29.9761, 31.1308),
    "Menkaure": (29.9725, 31.1280),
}

# Orion's Belt — J2000 ICRS coordinates (SIMBAD)
BELT_STARS = {
    "Alnitak":  SkyCoord(ra="05h40m45.527s", dec="-01d56m33.26s", frame="icrs"),
    "Alnilam":  SkyCoord(ra="05h36m12.813s", dec="-01d12m06.91s", frame="icrs"),
    "Mintaka":  SkyCoord(ra="05h32m00.400s", dec="-00d17m56.74s", frame="icrs"),
}

# Giza location for alt-az computation
GIZA = EarthLocation(lat=29.9792*u.deg, lon=31.1342*u.deg, height=60*u.m)


def pyramid_triangle_km():
    """Convert pyramid lat/lon to local XY in km, centered on centroid."""
    from math import radians, cos
    coords = list(PYRAMIDS.values())
    lats = [c[0] for c in coords]
    lons = [c[1] for c in coords]

    # Local tangent plane approximation (fine for <1 km scale)
    mean_lat = np.mean(lats)
    km_per_deg_lat = 111.132
    km_per_deg_lon = 111.320 * cos(radians(mean_lat))

    pts = np.array([
        [(lat - mean_lat) * km_per_deg_lat, (lon - np.mean(lons)) * km_per_deg_lon]
        for lat, lon in coords
    ])
    pts -= pts.mean(axis=0)  # center on centroid
    return pts


def belt_positions_at_epoch(year):
    """
    Get Orion's Belt star positions at a given epoch.

    For the shape comparison we need angular separations projected
    onto a 2D plane. We use the FK5 frame precessed to the target
    epoch to get RA/Dec, then project to a local tangent plane
    centered on the middle star (Alnilam).
    """
    t = Time(year, format='jyear', scale='tt')
    target_frame = FK5(equinox=t)

    positions = {}
    for name, coord in BELT_STARS.items():
        precessed = coord.transform_to(target_frame)
        positions[name] = (precessed.ra.deg, precessed.dec.deg)

    # Project to tangent plane centered on Alnilam
    ra_c, dec_c = positions["Alnilam"]
    pts = []
    for name in ["Alnitak", "Alnilam", "Mintaka"]:
        ra, dec = positions[name]
        # Gnomonic projection (small field, fine for ~3 deg)
        dra = np.radians(ra - ra_c) * np.cos(np.radians(dec_c))
        ddec = np.radians(dec - dec_c)
        pts.append([dra, ddec])

    pts = np.array(pts)
    pts -= pts.mean(axis=0)
    return pts


def procrustes_distance(A, B):
    """
    Compute Procrustes distance between two 2D point sets.
    Removes translation, rotation, and uniform scale.
    Returns the residual (0 = identical shape, 1 = maximally different).
    """
    # Normalize each to unit Frobenius norm (removes scale)
    A = A - A.mean(axis=0)
    B = B - B.mean(axis=0)
    A = A / np.linalg.norm(A)
    B = B / np.linalg.norm(B)

    # Optimal rotation via SVD
    M = A.T @ B
    U, S, Vt = np.linalg.svd(M)
    d = np.linalg.det(U @ Vt)
    D = np.diag([1, d])  # correct for reflection — this version does NOT allow reflection
    R = U @ D @ Vt
    B_rot = B @ R.T

    dist = np.sqrt(np.sum((A - B_rot)**2))
    return dist


def procrustes_distance_with_mirror(A, B):
    """
    Compute Procrustes distance allowing reflection (mirror).
    """
    A = A - A.mean(axis=0)
    B = B - B.mean(axis=0)
    A = A / np.linalg.norm(A)
    B = B / np.linalg.norm(B)

    M = A.T @ B
    U, S, Vt = np.linalg.svd(M)
    # Allow reflection: don't correct for det
    R = U @ Vt
    B_rot = B @ R.T

    dist = np.sqrt(np.sum((A - B_rot)**2))
    return dist


# ============================================================
# 2. TEST 1 & 2: Procrustes sweep + mirror test
# ============================================================

def run_sweep():
    """Compute Procrustes distance at 500-year intervals, with and without mirror."""
    pyr = pyramid_triangle_km()

    # Note the ordering: pyramids are [Khufu, Khafre, Menkaure] = roughly N to S
    # Stars are [Alnitak, Alnilam, Mintaka] = roughly S to N on the sky
    # The "natural" correspondence (Khufu=Alnitak, etc.) is itself part of the claim

    epochs = list(range(-15000, 2501, 500))
    results = []

    for year in epochs:
        stars = belt_positions_at_epoch(year)

        d_no_mirror = procrustes_distance(pyr, stars)
        d_mirror = procrustes_distance_with_mirror(pyr, stars)

        results.append({
            "epoch": year,
            "procrustes_no_mirror": round(float(d_no_mirror), 6),
            "procrustes_mirror": round(float(d_mirror), 6),
        })

    # Also compute at specific dates of interest
    for year in [-10500, -2550, 2026]:
        stars = belt_positions_at_epoch(year)
        d_no = procrustes_distance(pyr, stars)
        d_mi = procrustes_distance_with_mirror(pyr, stars)
        # Check if already in list
        if not any(r["epoch"] == year for r in results):
            results.append({
                "epoch": year,
                "procrustes_no_mirror": round(float(d_no), 6),
                "procrustes_mirror": round(float(d_mi), 6),
            })

    results.sort(key=lambda x: x["epoch"])
    return results


def find_best_epoch(results):
    """Find the epoch with minimum Procrustes distance."""
    best_no_mirror = min(results, key=lambda x: x["procrustes_no_mirror"])
    best_mirror = min(results, key=lambda x: x["procrustes_mirror"])
    return best_no_mirror, best_mirror


# ============================================================
# 3. TEST 3: Significance (Monte Carlo)
# ============================================================

def random_triangle():
    """Generate a random triangle (3 points in 2D, centered)."""
    pts = np.random.randn(3, 2)
    pts -= pts.mean(axis=0)
    return pts


def significance_test(n_trials=10000):
    """
    Generate random triangles and find their best Procrustes match
    to Orion's Belt across epochs. Compare to the pyramids' fit.
    """
    pyr = pyramid_triangle_km()

    # Pre-compute Belt positions at each epoch
    test_epochs = list(range(-15000, 2501, 500))
    belt_shapes = {}
    for year in test_epochs:
        belt_shapes[year] = belt_positions_at_epoch(year)

    # Get pyramid's best fits
    pyr_best_no_mirror = min(
        procrustes_distance(pyr, belt_shapes[y]) for y in test_epochs
    )
    pyr_best_mirror = min(
        procrustes_distance_with_mirror(pyr, belt_shapes[y]) for y in test_epochs
    )

    # Also get fits at specific epochs
    pyr_10500_no = procrustes_distance(pyr, belt_positions_at_epoch(-10500))
    pyr_10500_mi = procrustes_distance_with_mirror(pyr, belt_positions_at_epoch(-10500))
    pyr_2550_no = procrustes_distance(pyr, belt_positions_at_epoch(-2550))
    pyr_2550_mi = procrustes_distance_with_mirror(pyr, belt_positions_at_epoch(-2550))

    # Monte Carlo
    np.random.seed(42)
    count_better_no_mirror = 0
    count_better_mirror = 0
    count_better_10500_no = 0
    count_better_10500_mi = 0
    count_better_2550_no = 0
    count_better_2550_mi = 0
    random_best_distances = []

    for i in range(n_trials):
        tri = random_triangle()

        best_no = min(procrustes_distance(tri, belt_shapes[y]) for y in test_epochs)
        best_mi = min(procrustes_distance_with_mirror(tri, belt_shapes[y]) for y in test_epochs)

        random_best_distances.append({"no_mirror": float(best_no), "mirror": float(best_mi)})

        if best_no <= pyr_best_no_mirror:
            count_better_no_mirror += 1
        if best_mi <= pyr_best_mirror:
            count_better_mirror += 1

        # At specific epochs
        d_10500 = procrustes_distance(tri, belt_shapes[-10500])
        d_10500m = procrustes_distance_with_mirror(tri, belt_shapes[-10500])
        d_2550 = procrustes_distance(tri, belt_shapes[-2500])  # closest epoch in grid
        d_2550m = procrustes_distance_with_mirror(tri, belt_shapes[-2500])

        if d_10500 <= pyr_10500_no:
            count_better_10500_no += 1
        if d_10500m <= pyr_10500_mi:
            count_better_10500_mi += 1
        if d_2550 <= pyr_2550_no:
            count_better_2550_no += 1
        if d_2550m <= pyr_2550_mi:
            count_better_2550_mi += 1

        if (i + 1) % 2000 == 0:
            print(f"  Monte Carlo: {i+1}/{n_trials}")

    return {
        "n_trials": n_trials,
        "pyramid_best_fit": {
            "no_mirror": round(float(pyr_best_no_mirror), 6),
            "mirror": round(float(pyr_best_mirror), 6),
        },
        "pyramid_at_10500BC": {
            "no_mirror": round(float(pyr_10500_no), 6),
            "mirror": round(float(pyr_10500_mi), 6),
        },
        "pyramid_at_2550BC": {
            "no_mirror": round(float(pyr_2550_no), 6),
            "mirror": round(float(pyr_2550_mi), 6),
        },
        "p_value_best_fit": {
            "no_mirror": round(count_better_no_mirror / n_trials, 4),
            "mirror": round(count_better_mirror / n_trials, 4),
        },
        "p_value_at_10500BC": {
            "no_mirror": round(count_better_10500_no / n_trials, 4),
            "mirror": round(count_better_10500_mi / n_trials, 4),
        },
        "p_value_at_2550BC": {
            "no_mirror": round(count_better_2550_no / n_trials, 4),
            "mirror": round(count_better_2550_mi / n_trials, 4),
        },
        "random_best_distances_summary": {
            "no_mirror_mean": round(float(np.mean([d["no_mirror"] for d in random_best_distances])), 6),
            "no_mirror_std": round(float(np.std([d["no_mirror"] for d in random_best_distances])), 6),
            "mirror_mean": round(float(np.mean([d["mirror"] for d in random_best_distances])), 6),
            "mirror_std": round(float(np.std([d["mirror"] for d in random_best_distances])), 6),
        }
    }


# ============================================================
# 4. PLOTTING
# ============================================================

def plot_sweep(results):
    """Plot Procrustes distance over time."""
    epochs = [r["epoch"] for r in results]
    d_no = [r["procrustes_no_mirror"] for r in results]
    d_mi = [r["procrustes_mirror"] for r in results]

    fig, ax = plt.subplots(figsize=(14, 7))
    ax.plot(epochs, d_no, 'b-o', markersize=4, label='No mirror (rotation only)', linewidth=1.5)
    ax.plot(epochs, d_mi, 'r-s', markersize=4, label='With mirror allowed', linewidth=1.5)

    # Mark claimed dates
    ax.axvline(-10500, color='green', linestyle='--', alpha=0.7, label='10,500 BC (Hancock/Bauval)')
    ax.axvline(-2550, color='orange', linestyle='--', alpha=0.7, label='2,550 BC (accepted construction)')

    # Mark best fits
    best_no = min(results, key=lambda x: x["procrustes_no_mirror"])
    best_mi = min(results, key=lambda x: x["procrustes_mirror"])
    ax.annotate(f'Best no-mirror: {best_no["epoch"]}',
                xy=(best_no["epoch"], best_no["procrustes_no_mirror"]),
                xytext=(best_no["epoch"]+1500, best_no["procrustes_no_mirror"]+0.05),
                arrowprops=dict(arrowstyle='->', color='blue'),
                fontsize=10, color='blue')
    ax.annotate(f'Best mirror: {best_mi["epoch"]}',
                xy=(best_mi["epoch"], best_mi["procrustes_mirror"]),
                xytext=(best_mi["epoch"]+1500, best_mi["procrustes_mirror"]-0.05),
                arrowprops=dict(arrowstyle='->', color='red'),
                fontsize=10, color='red')

    ax.set_xlabel('Epoch (year)', fontsize=12)
    ax.set_ylabel('Procrustes Distance (0 = perfect match)', fontsize=12)
    ax.set_title('Giza Pyramids vs. Orion\'s Belt: Shape Similarity Over Time', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-15500, 3000)

    fig.tight_layout()
    fig.savefig(FIG / "procrustes_sweep.png", dpi=150)
    plt.close(fig)
    print(f"  Saved {FIG / 'procrustes_sweep.png'}")


def plot_triangle_comparison(year=-10500):
    """Plot the pyramid triangle vs Belt triangle at a given epoch."""
    pyr = pyramid_triangle_km()
    stars = belt_positions_at_epoch(year)

    # Normalize both for visual comparison
    pyr_n = pyr - pyr.mean(axis=0)
    pyr_n = pyr_n / np.linalg.norm(pyr_n)
    stars_n = stars - stars.mean(axis=0)
    stars_n = stars_n / np.linalg.norm(stars_n)

    # Optimal rotation (no mirror)
    M = pyr_n.T @ stars_n
    U, S, Vt = np.linalg.svd(M)
    d = np.linalg.det(U @ Vt)
    R = U @ np.diag([1, d]) @ Vt
    stars_aligned = stars_n @ R.T

    # Mirror version
    M2 = pyr_n.T @ stars_n
    U2, S2, Vt2 = np.linalg.svd(M2)
    R2 = U2 @ Vt2
    stars_mirror = stars_n @ R2.T

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    pyr_labels = list(PYRAMIDS.keys())
    star_labels = list(BELT_STARS.keys())

    # Panel 1: Raw shapes
    ax = axes[0]
    for i, label in enumerate(pyr_labels):
        ax.plot(pyr_n[i, 0], pyr_n[i, 1], 'bs', markersize=12)
        ax.annotate(label, (pyr_n[i, 0]+0.02, pyr_n[i, 1]+0.02), color='blue', fontsize=9)
    for i, label in enumerate(star_labels):
        ax.plot(stars_n[i, 0], stars_n[i, 1], 'r*', markersize=15)
        ax.annotate(label, (stars_n[i, 0]+0.02, stars_n[i, 1]-0.05), color='red', fontsize=9)
    # Connect triangles
    tri_pyr = np.vstack([pyr_n, pyr_n[0]])
    tri_star = np.vstack([stars_n, stars_n[0]])
    ax.plot(tri_pyr[:, 0], tri_pyr[:, 1], 'b-', alpha=0.5)
    ax.plot(tri_star[:, 0], tri_star[:, 1], 'r-', alpha=0.5)
    ax.set_title(f'Raw normalized shapes ({year})', fontsize=11)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    # Panel 2: Best rotation (no mirror)
    ax = axes[1]
    for i, label in enumerate(pyr_labels):
        ax.plot(pyr_n[i, 0], pyr_n[i, 1], 'bs', markersize=12)
        ax.annotate(label, (pyr_n[i, 0]+0.02, pyr_n[i, 1]+0.02), color='blue', fontsize=9)
    for i, label in enumerate(star_labels):
        ax.plot(stars_aligned[i, 0], stars_aligned[i, 1], 'r*', markersize=15)
        ax.annotate(label, (stars_aligned[i, 0]+0.02, stars_aligned[i, 1]-0.05), color='red', fontsize=9)
    tri_aligned = np.vstack([stars_aligned, stars_aligned[0]])
    ax.plot(tri_pyr[:, 0], tri_pyr[:, 1], 'b-', alpha=0.5)
    ax.plot(tri_aligned[:, 0], tri_aligned[:, 1], 'r-', alpha=0.5)
    d_no = procrustes_distance(pyr, stars)
    ax.set_title(f'Best rotation (no mirror)\nProcrustes d = {d_no:.4f}', fontsize=11)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    # Panel 3: Best fit allowing mirror
    ax = axes[2]
    for i, label in enumerate(pyr_labels):
        ax.plot(pyr_n[i, 0], pyr_n[i, 1], 'bs', markersize=12)
        ax.annotate(label, (pyr_n[i, 0]+0.02, pyr_n[i, 1]+0.02), color='blue', fontsize=9)
    for i, label in enumerate(star_labels):
        ax.plot(stars_mirror[i, 0], stars_mirror[i, 1], 'r*', markersize=15)
        ax.annotate(label, (stars_mirror[i, 0]+0.02, stars_mirror[i, 1]-0.05), color='red', fontsize=9)
    tri_mirror = np.vstack([stars_mirror, stars_mirror[0]])
    ax.plot(tri_pyr[:, 0], tri_pyr[:, 1], 'b-', alpha=0.5)
    ax.plot(tri_mirror[:, 0], tri_mirror[:, 1], 'r-', alpha=0.5)
    d_mi = procrustes_distance_with_mirror(pyr, stars)
    ax.set_title(f'Best fit (mirror allowed)\nProcrustes d = {d_mi:.4f}', fontsize=11)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    fig.suptitle(f'Pyramid vs. Orion\'s Belt Triangle Comparison at {abs(year)} BC', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / f"triangle_comparison_{abs(year)}BC.png", dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved triangle comparison at {year}")


def plot_significance(sig_results):
    """Plot histogram of random triangle best-fit distances vs pyramid."""
    # Regenerate random best distances for histogram
    np.random.seed(42)
    test_epochs = list(range(-15000, 2501, 500))
    belt_shapes = {y: belt_positions_at_epoch(y) for y in test_epochs}

    rand_no = []
    rand_mi = []
    for _ in range(sig_results["n_trials"]):
        tri = random_triangle()
        best_no = min(procrustes_distance(tri, belt_shapes[y]) for y in test_epochs)
        best_mi = min(procrustes_distance_with_mirror(tri, belt_shapes[y]) for y in test_epochs)
        rand_no.append(best_no)
        rand_mi.append(best_mi)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    ax.hist(rand_no, bins=50, alpha=0.7, color='steelblue', edgecolor='black', linewidth=0.5)
    pyr_val = sig_results["pyramid_best_fit"]["no_mirror"]
    ax.axvline(pyr_val, color='red', linewidth=2, label=f'Pyramids: {pyr_val:.4f}')
    ax.set_xlabel('Best Procrustes Distance (no mirror)')
    ax.set_ylabel('Count')
    ax.set_title(f'No Mirror — p = {sig_results["p_value_best_fit"]["no_mirror"]}')
    ax.legend()

    ax = axes[1]
    ax.hist(rand_mi, bins=50, alpha=0.7, color='salmon', edgecolor='black', linewidth=0.5)
    pyr_val = sig_results["pyramid_best_fit"]["mirror"]
    ax.axvline(pyr_val, color='red', linewidth=2, label=f'Pyramids: {pyr_val:.4f}')
    ax.set_xlabel('Best Procrustes Distance (mirror allowed)')
    ax.set_ylabel('Count')
    ax.set_title(f'Mirror Allowed — p = {sig_results["p_value_best_fit"]["mirror"]}')
    ax.legend()

    fig.suptitle('Significance Test: Pyramids vs. 10,000 Random Triangles', fontsize=14)
    fig.tight_layout()
    fig.savefig(FIG / "significance_test.png", dpi=150)
    plt.close(fig)
    print(f"  Saved significance plot")


# ============================================================
# 5. MAIN
# ============================================================

def main():
    print("=" * 60)
    print("ORION CORRELATION FORMAL FIT TEST")
    print("=" * 60)

    # --- Pyramid triangle ---
    pyr = pyramid_triangle_km()
    print(f"\nPyramid triangle (km from centroid):")
    for name, pt in zip(PYRAMIDS.keys(), pyr):
        print(f"  {name}: ({pt[0]:.4f}, {pt[1]:.4f}) km")

    # Compute side lengths
    sides_pyr = [
        np.linalg.norm(pyr[1] - pyr[0]),
        np.linalg.norm(pyr[2] - pyr[1]),
        np.linalg.norm(pyr[2] - pyr[0]),
    ]
    print(f"  Side lengths: {[f'{s:.4f} km' for s in sides_pyr]}")

    # --- Belt star separations at key epochs ---
    for year in [-10500, -2550, 2026]:
        stars = belt_positions_at_epoch(year)
        seps = [
            np.linalg.norm(stars[1] - stars[0]),
            np.linalg.norm(stars[2] - stars[1]),
            np.linalg.norm(stars[2] - stars[0]),
        ]
        print(f"\nOrion's Belt at {year} (normalized angular separations):")
        for name, pt in zip(BELT_STARS.keys(), stars):
            print(f"  {name}: ({pt[0]:.6f}, {pt[1]:.6f})")
        print(f"  Side ratios: {[f'{s:.6f}' for s in seps]}")

    # --- Test 1 & 2: Procrustes sweep ---
    print("\n" + "-" * 60)
    print("TEST 1 & 2: Procrustes Sweep (500-year intervals, -15000 to 2500)")
    print("-" * 60)
    results = run_sweep()

    best_no, best_mi = find_best_epoch(results)
    print(f"\nBest fit (no mirror):    epoch {best_no['epoch']}, d = {best_no['procrustes_no_mirror']:.6f}")
    print(f"Best fit (mirror OK):    epoch {best_mi['epoch']}, d = {best_mi['procrustes_mirror']:.6f}")

    # Key epochs
    for r in results:
        if r["epoch"] in [-10500, -2550, 2026]:
            print(f"  Epoch {r['epoch']:>6}: no_mirror={r['procrustes_no_mirror']:.6f}, mirror={r['procrustes_mirror']:.6f}")

    # Save sweep
    with open(OUT / "procrustes_sweep.json", 'w') as f:
        json.dump({
            "description": "Procrustes distance between Giza pyramids and Orion's Belt at each epoch",
            "method": "Shape comparison after removing translation, rotation, and scale. Mirror column additionally removes reflection.",
            "pyramid_positions_km": {name: pt.tolist() for name, pt in zip(PYRAMIDS.keys(), pyr)},
            "epochs": results,
            "best_fit_no_mirror": best_no,
            "best_fit_mirror": best_mi,
        }, f, indent=2)
    print(f"  Saved {OUT / 'procrustes_sweep.json'}")

    # --- Mirror test ---
    print("\n" + "-" * 60)
    print("TEST 2: Mirror Inversion Analysis")
    print("-" * 60)

    mirror_results = {}
    for year in [-10500, -2550, 2026]:
        stars = belt_positions_at_epoch(year)
        d_no = procrustes_distance(pyr, stars)
        d_mi = procrustes_distance_with_mirror(pyr, stars)
        improvement = (d_no - d_mi) / d_no * 100 if d_no > 0 else 0
        mirror_results[str(year)] = {
            "no_mirror": round(float(d_no), 6),
            "mirror": round(float(d_mi), 6),
            "improvement_pct": round(float(improvement), 2),
            "mirror_needed": bool(d_mi < d_no * 0.95),  # >5% improvement means mirror helps significantly
        }
        print(f"  {year}: no_mirror={d_no:.6f}, mirror={d_mi:.6f}, improvement={improvement:.1f}%")
        print(f"         Mirror significantly helps: {d_mi < d_no * 0.95}")

    with open(OUT / "mirror_test.json", 'w') as f:
        json.dump({
            "description": "Comparison of fit with and without allowing mirror reflection",
            "interpretation": "If mirror_needed=true, the pyramids are reflected relative to Orion (weakens the claim)",
            "epochs": mirror_results,
        }, f, indent=2)
    print(f"  Saved {OUT / 'mirror_test.json'}")

    # --- Test 3: Significance ---
    print("\n" + "-" * 60)
    print("TEST 3: Monte Carlo Significance (10,000 random triangles)")
    print("-" * 60)
    sig = significance_test(n_trials=10000)

    print(f"\nPyramid best-fit Procrustes distance:")
    print(f"  No mirror: {sig['pyramid_best_fit']['no_mirror']:.6f}")
    print(f"  Mirror:    {sig['pyramid_best_fit']['mirror']:.6f}")
    print(f"\nRandom triangle mean best-fit:")
    print(f"  No mirror: {sig['random_best_distances_summary']['no_mirror_mean']:.6f} ± {sig['random_best_distances_summary']['no_mirror_std']:.6f}")
    print(f"  Mirror:    {sig['random_best_distances_summary']['mirror_mean']:.6f} ± {sig['random_best_distances_summary']['mirror_std']:.6f}")
    print(f"\np-values (fraction of random triangles with equal or better fit):")
    print(f"  Best fit anywhere: no_mirror p={sig['p_value_best_fit']['no_mirror']}, mirror p={sig['p_value_best_fit']['mirror']}")
    print(f"  At 10,500 BC:      no_mirror p={sig['p_value_at_10500BC']['no_mirror']}, mirror p={sig['p_value_at_10500BC']['mirror']}")
    print(f"  At 2,550 BC:       no_mirror p={sig['p_value_at_2550BC']['no_mirror']}, mirror p={sig['p_value_at_2550BC']['mirror']}")

    with open(OUT / "significance.json", 'w') as f:
        json.dump(sig, f, indent=2)
    print(f"  Saved {OUT / 'significance.json'}")

    # --- Plots ---
    print("\n" + "-" * 60)
    print("GENERATING FIGURES")
    print("-" * 60)
    plot_sweep(results)
    plot_triangle_comparison(-10500)
    plot_triangle_comparison(-2550)
    plot_significance(sig)

    # --- Summary ---
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Best shape match (no mirror): {best_no['epoch']} (d={best_no['procrustes_no_mirror']:.6f})")
    print(f"Best shape match (mirror OK): {best_mi['epoch']} (d={best_mi['procrustes_mirror']:.6f})")
    print(f"Match at 10,500 BC: {[r for r in results if r['epoch']==-10500][0]['procrustes_no_mirror']:.6f} (no mirror)")
    print(f"Match at 2,550 BC:  {[r for r in results if r['epoch']==-2550][0]['procrustes_no_mirror']:.6f} (no mirror)")

    if sig['p_value_best_fit']['mirror'] > 0.05:
        print("\nCONCLUSION: The pyramid-Orion match is NOT statistically significant.")
        print("Random triangles achieve equally good fits frequently.")
    elif sig['p_value_best_fit']['mirror'] <= 0.05 and sig['p_value_best_fit']['no_mirror'] > 0.05:
        print("\nCONCLUSION: The match is only significant when allowing mirror reflection,")
        print("which weakens the 'intentional design' interpretation.")
    else:
        print("\nCONCLUSION: The match IS statistically significant even without mirror.")
        print("However, the best-fit epoch should be checked against claimed dates.")

    best_at_10500 = best_no['epoch'] == -10500 or best_mi['epoch'] == -10500
    print(f"\nDoes the best match occur at 10,500 BC? {'YES' if best_at_10500 else 'NO'}")
    if not best_at_10500:
        print(f"  Best match is at {best_no['epoch']} (no mirror) / {best_mi['epoch']} (mirror)")
        print(f"  This undermines the specific 10,500 BC claim even if the overall match is good.")


if __name__ == "__main__":
    main()
