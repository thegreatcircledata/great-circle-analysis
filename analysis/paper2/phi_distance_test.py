#!/usr/bin/env python3
"""
Directive 1: Alison Golden Section (φ) Distance Test

Tests whether the great-circle distances between Giza, Nazca, and Angkor Wat
fall in golden-ratio proportion by design or by chance.

Claim (Jim Alison):
  Angkor→Giza (4,754 mi) × φ ≈ Giza→Nazca (7,692 mi)
  The three arcs divide the circumference at ~19.1%, 30.9%, 50.0%
  These percentages match first digits of Fibonacci #137-139.
"""

import json
import math
import numpy as np
from pathlib import Path

# ---------------------------------------------------------------------------
# 1. Site coordinates (WGS-84)
# ---------------------------------------------------------------------------
SITES = {
    "giza":   (29.9792, 31.1342),   # Great Pyramid
    "nazca":  (-14.7350, -75.1300),  # Nazca Lines center
    "angkor": (13.4125, 103.8670),   # Angkor Wat
}

EARTH_RADIUS_KM = 6371.0
KM_TO_MI = 0.621371

# ---------------------------------------------------------------------------
# 2. Great-circle distance (Vincenty formula for numerical stability)
# ---------------------------------------------------------------------------
def great_circle_km(lat1, lon1, lat2, lon2):
    """Haversine great-circle distance in km."""
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_KM * math.asin(math.sqrt(a))


# ---------------------------------------------------------------------------
# 3. Compute actual distances
# ---------------------------------------------------------------------------
def compute_distances():
    d_ag = great_circle_km(*SITES["angkor"], *SITES["giza"])
    d_gn = great_circle_km(*SITES["giza"], *SITES["nazca"])
    d_na = great_circle_km(*SITES["nazca"], *SITES["angkor"])
    circumference = 2 * math.pi * EARTH_RADIUS_KM
    return d_ag, d_gn, d_na, circumference


# ---------------------------------------------------------------------------
# 4. Fibonacci first-digits helper
# ---------------------------------------------------------------------------
def fibonacci_first_digits(n, num_digits=3):
    """Return first `num_digits` digits of F(n) using Binet + log trick."""
    # log10(F(n)) ≈ n*log10(φ) - log10(√5)  for large n
    phi = (1 + math.sqrt(5)) / 2
    log_fn = n * math.log10(phi) - 0.5 * math.log10(5)
    frac = log_fn - math.floor(log_fn)
    leading = 10 ** frac
    return int(leading * 10 ** (num_digits - 1))


# ---------------------------------------------------------------------------
# 5. Monte Carlo: φ-ratio test for 3 random points on a great circle
# ---------------------------------------------------------------------------
def monte_carlo_phi(n_trials, tolerance_frac, rng):
    """
    Place 3 random points on a great circle, compute the 3 arc gaps,
    check if any pair of gaps has ratio within `tolerance_frac` of φ.
    Returns fraction of trials where at least one pair matches.
    """
    phi = (1 + math.sqrt(5)) / 2
    hits = 0

    for _ in range(n_trials):
        positions = np.sort(rng.uniform(0, 360, 3))
        gaps = np.diff(np.append(positions, positions[0] + 360))
        # gaps has 3 elements summing to 360

        found = False
        for i in range(3):
            if found:
                break
            for j in range(3):
                if i != j and gaps[j] > 0:
                    ratio = gaps[i] / gaps[j]
                    if abs(ratio - phi) / phi <= tolerance_frac:
                        found = True
                        break
        if found:
            hits += 1

    return hits / n_trials


# ---------------------------------------------------------------------------
# 6. Monte Carlo: Fibonacci digit match test
# ---------------------------------------------------------------------------
def monte_carlo_fibonacci_digits(n_trials, rng):
    """
    Place 3 random points on a circle, compute arc percentages,
    check if the sorted percentages' leading 3 digits match any
    3 consecutive Fibonacci numbers' leading 3 digits.

    We pre-compute leading-3-digit sets for consecutive Fibonacci triples
    up to F(200).
    """
    # Build lookup of consecutive Fibonacci first-digit triples
    fib_triples = set()
    for n in range(10, 200):
        triple = tuple(sorted([
            fibonacci_first_digits(n, 3),
            fibonacci_first_digits(n + 1, 3),
            fibonacci_first_digits(n + 2, 3),
        ]))
        fib_triples.add(triple)

    hits = 0
    for _ in range(n_trials):
        positions = np.sort(rng.uniform(0, 360, 3))
        gaps = np.diff(np.append(positions, positions[0] + 360))
        pcts = np.sort(gaps / 360 * 100)
        # Extract leading 3 digits of each percentage
        digits = []
        for p in pcts:
            if p < 1:
                continue
            s = f"{p:.10f}".replace(".", "")
            s = s.lstrip("0")
            if len(s) >= 3:
                digits.append(int(s[:3]))
        if len(digits) == 3:
            if tuple(sorted(digits)) in fib_triples:
                hits += 1

    return hits / n_trials


# ---------------------------------------------------------------------------
# 7. Main
# ---------------------------------------------------------------------------
def main():
    phi = (1 + math.sqrt(5)) / 2
    rng = np.random.default_rng(42)
    n_trials = 200_000

    # --- Compute distances ---
    d_ag, d_gn, d_na, circumference = compute_distances()

    d_ag_mi = d_ag * KM_TO_MI
    d_gn_mi = d_gn * KM_TO_MI
    d_na_mi = d_na * KM_TO_MI

    actual_ratio = d_gn / d_ag
    phi_error_frac = abs(actual_ratio - phi) / phi

    pct_ag = d_ag / circumference * 100
    pct_gn = d_gn / circumference * 100
    pct_na = d_na / circumference * 100

    print("=" * 65)
    print("DIRECTIVE 1: Alison Golden Section (φ) Distance Test")
    print("=" * 65)
    print(f"\nSite coordinates:")
    for name, (lat, lon) in SITES.items():
        print(f"  {name:8s}: {lat:+8.4f}°, {lon:+9.4f}°")

    print(f"\nGreat-circle distances:")
    print(f"  Angkor → Giza : {d_ag:,.1f} km  ({d_ag_mi:,.1f} mi)")
    print(f"  Giza   → Nazca: {d_gn:,.1f} km  ({d_gn_mi:,.1f} mi)")
    print(f"  Nazca  → Angkor: {d_na:,.1f} km  ({d_na_mi:,.1f} mi)")
    print(f"  Sum           : {d_ag + d_gn + d_na:,.1f} km")
    print(f"  Circumference : {circumference:,.1f} km")

    print(f"\nCircumference fractions:")
    print(f"  Angkor→Giza  : {pct_ag:.2f}%")
    print(f"  Giza→Nazca   : {pct_gn:.2f}%")
    print(f"  Nazca→Angkor : {pct_na:.2f}%")

    print(f"\nGolden ratio test:")
    print(f"  Giza→Nazca / Angkor→Giza = {actual_ratio:.6f}")
    print(f"  φ                         = {phi:.6f}")
    print(f"  Relative error            = {phi_error_frac:.4%}")

    # --- Fibonacci digit check ---
    fib137 = fibonacci_first_digits(137, 3)
    fib138 = fibonacci_first_digits(138, 3)
    fib139 = fibonacci_first_digits(139, 3)
    print(f"\nFibonacci leading-digit check:")
    print(f"  Arc percentages (sorted): {sorted([pct_ag, pct_gn, pct_na])}")
    print(f"  F(137) first 3 digits: {fib137}")
    print(f"  F(138) first 3 digits: {fib138}")
    print(f"  F(139) first 3 digits: {fib139}")

    # --- Monte Carlo: φ ratio ---
    print(f"\nMonte Carlo: φ-ratio test ({n_trials:,} trials)...")
    p_phi = monte_carlo_phi(n_trials, phi_error_frac, rng)
    print(f"  Tolerance (relative): {phi_error_frac:.4%}")
    print(f"  p-value (φ match): {p_phi:.6f}")

    # --- Monte Carlo: Fibonacci digits ---
    print(f"\nMonte Carlo: Fibonacci digit test ({n_trials:,} trials)...")
    p_fib = monte_carlo_fibonacci_digits(n_trials, rng)
    print(f"  p-value (Fibonacci triple match): {p_fib:.6f}")

    # --- Verdict ---
    # Important context: Giza, Nazca, Angkor are NOT on the same great circle.
    # Alison's claim requires cherry-picking the direction of measurement.
    # Also check how far off the three sites are from a shared great circle.
    #
    # Fit the best great circle to the 3 sites and measure max deviation.
    coords_rad = []
    for name in ["giza", "nazca", "angkor"]:
        lat, lon = SITES[name]
        coords_rad.append((math.radians(lat), math.radians(lon)))

    # Convert to 3D unit vectors
    def to_xyz(lat_r, lon_r):
        return (
            math.cos(lat_r) * math.cos(lon_r),
            math.cos(lat_r) * math.sin(lon_r),
            math.sin(lat_r),
        )

    pts = [to_xyz(*c) for c in coords_rad]
    # Any 3 non-collinear points define a unique great circle (plane through origin)
    # The "pole" of the great circle through 2 points is their cross product
    # Check deviation of each point from the great circle defined by the other two
    def cross(a, b):
        return (
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        )

    def dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    def norm(a):
        return math.sqrt(dot(a, a))

    # Great circle through all 3 is trivially exact (3 points define a plane).
    # The real question is whether they lie on a meaningful great circle (e.g.,
    # Alison's "equator line"). The distances along the great circle through
    # all three give us the arcs used above. This is always possible for any 3
    # non-antipodal points, which weakens the claim.

    # The total arc is d_ag + d_gn + d_na. If this equals circumference,
    # they split the full circle. If not, they form a "triangle" and the
    # measurement direction is ambiguous.
    total_arc = d_ag + d_gn + d_na
    arc_vs_circ = total_arc / circumference
    arc_deficit_pct = abs(1 - arc_vs_circ) * 100

    print(f"\nGreat-circle geometry check:")
    print(f"  Sum of 3 arcs / circumference = {arc_vs_circ:.4f}")
    print(f"  Deficit from full circle: {arc_deficit_pct:.2f}%")
    if arc_deficit_pct > 1:
        print(f"  ⚠ The three sites do NOT lie on a single great circle.")
        print(f"    The arcs sum to {total_arc:,.0f} km vs circumference {circumference:,.0f} km.")
        print(f"    Alison's claim requires choosing measurement direction to force fit.")

    # --- Compile verdict ---
    notes = []
    if phi_error_frac < 0.01:
        notes.append(f"φ ratio holds to {phi_error_frac:.2%} — impressively close.")
    else:
        notes.append(f"φ ratio error is {phi_error_frac:.2%} — moderate fit.")

    if p_phi < 0.05:
        notes.append(f"Monte Carlo p={p_phi:.4f}: φ match is statistically unusual.")
    else:
        notes.append(f"Monte Carlo p={p_phi:.4f}: φ match is NOT statistically unusual.")

    if arc_deficit_pct > 1:
        notes.append(
            "Critical flaw: sites are not on a shared great circle. "
            "Arc directions were chosen post-hoc to produce the ratio."
        )

    if p_phi < 0.05 and arc_deficit_pct <= 1:
        verdict = "POSSIBLY MEANINGFUL — φ ratio is statistically unusual and geometry is valid."
    elif p_phi < 0.05 and arc_deficit_pct > 1:
        verdict = "LIKELY COINCIDENCE — φ ratio is tight but geometry is cherry-picked."
    else:
        verdict = "COINCIDENCE — φ ratio is within normal random range for 3 points on a circle."

    print(f"\n{'=' * 65}")
    print(f"VERDICT: {verdict}")
    for n in notes:
        print(f"  • {n}")
    print(f"{'=' * 65}")

    # --- Save results ---
    results = {
        "directive": "Directive 1: Alison Golden Section (φ) Distance Test",
        "sites": {k: {"lat": v[0], "lon": v[1]} for k, v in SITES.items()},
        "distances_km": {
            "angkor_to_giza": round(d_ag, 1),
            "giza_to_nazca": round(d_gn, 1),
            "nazca_to_angkor": round(d_na, 1),
            "sum": round(d_ag + d_gn + d_na, 1),
            "earth_circumference": round(circumference, 1),
        },
        "distances_mi": {
            "angkor_to_giza": round(d_ag_mi, 1),
            "giza_to_nazca": round(d_gn_mi, 1),
            "nazca_to_angkor": round(d_na_mi, 1),
        },
        "circumference_fractions_pct": {
            "angkor_to_giza": round(pct_ag, 4),
            "giza_to_nazca": round(pct_gn, 4),
            "nazca_to_angkor": round(pct_na, 4),
        },
        "phi_ratio_test": {
            "actual_ratio": round(actual_ratio, 6),
            "phi": round(phi, 6),
            "relative_error": round(phi_error_frac, 6),
            "relative_error_pct": f"{phi_error_frac:.4%}",
        },
        "fibonacci_digits": {
            "F137_first3": fib137,
            "F138_first3": fib138,
            "F139_first3": fib139,
            "arc_pcts_sorted": sorted([round(pct_ag, 2), round(pct_gn, 2), round(pct_na, 2)]),
        },
        "geometry_check": {
            "arc_sum_vs_circumference": round(arc_vs_circ, 4),
            "deficit_pct": round(arc_deficit_pct, 2),
            "on_shared_great_circle": arc_deficit_pct <= 1,
        },
        "monte_carlo": {
            "n_trials": n_trials,
            "p_value_phi_ratio": round(p_phi, 6),
            "p_value_fibonacci_digits": round(p_fib, 6),
            "tolerance_relative": round(phi_error_frac, 6),
        },
        "verdict": verdict,
        "notes": notes,
    }

    out_path = Path("outputs/extended_analysis/phi_distance_test.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
