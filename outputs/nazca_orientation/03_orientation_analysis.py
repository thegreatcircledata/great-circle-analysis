#!/usr/bin/env python3
"""
Analysis 2: Nazca Line Orientations vs Great Circle Bearing
============================================================
Uses the Richter et al. 2021 (Applied Sciences 11(4):1637) histogram of
2,308 straight line azimuths from the NascaGIS database.

Tests whether the GC bearing at Nazca (63.14°) is a preferred orientation
of the Nazca lines.

Methods:
  A. Bin-level comparison: where does the GC bearing fall in the distribution?
  B. Rayleigh test on doubled angles (axial test) relative to GC bearing
  C. Monte Carlo: random test bearings → percentile rank of actual GC result
  D. Kuiper test for non-uniformity of the distribution itself
"""

import json
import math
import numpy as np
from scipy import stats

# ─── Data ────────────────────────────────────────────────────────
# Richter et al. 2021 Figure 13 histogram (10° bins, 0-180° axial range).
# Each bin center and count. The 180-360° half is the symmetric mirror.
# Total = 2,308 lines.

BINS_CENTER = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115,
               125, 135, 145, 155, 165, 175]
BIN_COUNTS = [91, 126, 103, 100, 153, 127, 177, 116, 129, 76,
              104, 86, 165, 140, 172, 177, 124, 142]

N_LINES = sum(BIN_COUNTS)  # Should be 2308
BIN_WIDTH = 10.0

# Great Circle bearing at Nazca
GC_BEARING = 63.14  # degrees from north (Analysis 1 result)
GC_BEARING_AXIAL = GC_BEARING % 180  # 63.14° (axial orientation)

# Known reference azimuths at Nazca latitude (-14.7°)
SOLAR_AZIMUTHS = {
    "June_solstice_sunrise":    65.6,
    "Equinox_sunrise":          90.0,
    "December_solstice_sunrise": 114.4,
    "December_solstice_sunset": 245.6 - 180,  # axial: 65.6
    "Equinox_sunset":           270.0 - 180,   # axial: 90.0
    "June_solstice_sunset":     294.4 - 180,    # axial: 114.4
}

print("=" * 65)
print("NAZCA LINE ORIENTATIONS vs GREAT CIRCLE BEARING")
print(f"GC bearing at Nazca: {GC_BEARING:.2f}° (axial: {GC_BEARING_AXIAL:.2f}°)")
print(f"Total lines (Richter et al. 2021): {N_LINES}")
print("=" * 65)

# ═══════════════════════════════════════════════════════════════════
# STEP A: Bin-level comparison
# ═══════════════════════════════════════════════════════════════════

print("\n--- A. BIN-LEVEL COMPARISON ---")

# Which bin does the GC bearing fall in?
gc_bin_idx = int(GC_BEARING_AXIAL // BIN_WIDTH)
gc_bin_center = BINS_CENTER[gc_bin_idx]
gc_bin_count = BIN_COUNTS[gc_bin_idx]
expected_per_bin = N_LINES / 18  # uniform expectation

print(f"  GC bearing {GC_BEARING_AXIAL:.1f}° falls in bin {gc_bin_center-5}-{gc_bin_center+5}°")
print(f"  Bin count: {gc_bin_count} lines")
print(f"  Uniform expectation: {expected_per_bin:.1f} lines/bin")
print(f"  Ratio: {gc_bin_count / expected_per_bin:.2f}x expected")

# Rank of this bin among all bins
sorted_counts = sorted(BIN_COUNTS, reverse=True)
rank = sorted_counts.index(gc_bin_count) + 1
print(f"  Bin rank: {rank} of {len(BIN_COUNTS)} (1 = highest count)")

# Chi-squared test for this bin vs uniform
chi2_bin = (gc_bin_count - expected_per_bin)**2 / expected_per_bin
print(f"  Chi-squared (this bin): {chi2_bin:.2f}")

# Overall chi-squared for uniformity
chi2_total = sum((c - expected_per_bin)**2 / expected_per_bin for c in BIN_COUNTS)
chi2_p = 1 - stats.chi2.cdf(chi2_total, df=len(BIN_COUNTS)-1)
print(f"\n  Overall chi-squared for uniformity: {chi2_total:.2f} (df={len(BIN_COUNTS)-1})")
print(f"  p-value: {chi2_p:.6f}")
if chi2_p < 0.001:
    print("  → Distribution is SIGNIFICANTLY non-uniform (p < 0.001)")
elif chi2_p < 0.05:
    print("  → Distribution is significantly non-uniform (p < 0.05)")
else:
    print("  → Distribution is consistent with uniform")

# Print full distribution
print(f"\n  Full distribution (bin center: count, ratio):")
for bc, cnt in zip(BINS_CENTER, BIN_COUNTS):
    ratio = cnt / expected_per_bin
    bar = "#" * int(ratio * 20)
    marker = " ← GC" if bc == gc_bin_center else ""
    print(f"    {bc:4d}°: {cnt:4d}  ({ratio:.2f}x)  {bar}{marker}")


# ═══════════════════════════════════════════════════════════════════
# STEP B: Circular statistics — reconstruct individual angles
# ═══════════════════════════════════════════════════════════════════

print("\n--- B. CIRCULAR STATISTICS ---")

# Reconstruct individual angles from histogram (jittered within bins)
np.random.seed(42)
angles = []
for bc, cnt in zip(BINS_CENTER, BIN_COUNTS):
    # Place angles uniformly within each bin
    angles.extend(np.random.uniform(bc - BIN_WIDTH/2, bc + BIN_WIDTH/2, cnt))
angles = np.array(angles)

# For axial data, double the angles (map 0-180 to 0-360)
doubled = 2 * np.radians(angles)

# Rayleigh test for non-uniformity
C = np.mean(np.cos(doubled))
S = np.mean(np.sin(doubled))
R_bar = np.sqrt(C**2 + S**2)
mean_direction_doubled = np.degrees(np.arctan2(S, C)) % 360
mean_direction = (mean_direction_doubled / 2) % 180  # Back to axial

# Rayleigh statistic
n = len(angles)
R = n * R_bar
Z = R_bar**2 * n  # Rayleigh Z statistic
rayleigh_p = np.exp(-Z) * (1 + (2*Z - Z**2) / (4*n) - (24*Z - 132*Z**2 + 76*Z**3 - 9*Z**4) / (288*n**2))
rayleigh_p = max(rayleigh_p, 0)

print(f"  Rayleigh test (axial, doubled angles):")
print(f"    Mean resultant length R̄ = {R_bar:.4f}")
print(f"    Mean direction (axial) = {mean_direction:.1f}°")
print(f"    Rayleigh Z = {Z:.2f}")
print(f"    p-value ≈ {rayleigh_p:.6f}")
if rayleigh_p < 0.05:
    print(f"    → Significant preferred direction at {mean_direction:.1f}°")
else:
    print(f"    → No significant preferred direction")

# V-test against GC bearing
gc_rad_doubled = 2 * math.radians(GC_BEARING_AXIAL)
V = R_bar * math.cos(np.arctan2(S, C) - gc_rad_doubled)
u = V * math.sqrt(2 * n)
v_p = 1 - stats.norm.cdf(u)

print(f"\n  V-test (directed test against GC bearing {GC_BEARING_AXIAL:.1f}°):")
print(f"    V = {V:.4f}")
print(f"    u = {u:.2f}")
print(f"    p-value = {v_p:.6f}")
if v_p < 0.05:
    print(f"    → Lines significantly cluster toward GC bearing")
else:
    print(f"    → No significant clustering toward GC bearing")


# ═══════════════════════════════════════════════════════════════════
# STEP C: Monte Carlo — percentile rank of GC bearing
# ═══════════════════════════════════════════════════════════════════

print("\n--- C. MONTE CARLO: RANDOM TEST BEARINGS ---")

N_MC = 10000
np.random.seed(42)

# For each random test bearing, compute the V-test statistic
random_bearings = np.random.uniform(0, 180, N_MC)
v_stats = np.zeros(N_MC)

for i, test_b in enumerate(random_bearings):
    test_rad_doubled = 2 * math.radians(test_b)
    v_stats[i] = R_bar * math.cos(np.arctan2(S, C) - test_rad_doubled)

# GC's V-statistic
gc_v = R_bar * math.cos(np.arctan2(S, C) - gc_rad_doubled)
percentile = np.mean(v_stats <= gc_v) * 100

print(f"  Random test bearings: {N_MC}")
print(f"  GC V-statistic: {gc_v:.4f}")
print(f"  Percentile rank: {percentile:.1f}%")
print(f"  (100% = GC is perfectly aligned with the preferred direction)")

# Also test: what fraction of random bearings fall in a bin ≥ as populated?
gc_bin_count_rank = np.zeros(N_MC)
for i, test_b in enumerate(random_bearings):
    bin_idx = min(int(test_b // BIN_WIDTH), 17)
    gc_bin_count_rank[i] = BIN_COUNTS[bin_idx]

pct_higher_density = np.mean(gc_bin_count_rank >= gc_bin_count) * 100
print(f"\n  Bin density test:")
print(f"    GC bin count: {gc_bin_count}")
print(f"    % of random bearings in bins ≥ as dense: {pct_higher_density:.1f}%")


# ═══════════════════════════════════════════════════════════════════
# STEP D: Comparison with solar azimuths
# ═══════════════════════════════════════════════════════════════════

print("\n--- D. COMPARISON WITH SOLAR AZIMUTHS ---")

print(f"\n  {'Azimuth':35s}  {'Value':>7s}  {'Bin count':>9s}  {'Rank':>4s}  {'Δ from GC':>9s}")
print(f"  {'-'*35}  {'-'*7}  {'-'*9}  {'-'*4}  {'-'*9}")

all_refs = {"Great_Circle": GC_BEARING_AXIAL}
all_refs.update(SOLAR_AZIMUTHS)

for name, az in all_refs.items():
    ax = az % 180
    bin_idx = min(int(ax // BIN_WIDTH), 17)
    cnt = BIN_COUNTS[bin_idx]
    rank = sorted_counts.index(cnt) + 1
    diff_gc = abs(ax - GC_BEARING_AXIAL)
    diff_gc = min(diff_gc, 180 - diff_gc)
    marker = " ← TEST" if name == "Great_Circle" else ""
    print(f"  {name:35s}  {ax:7.1f}°  {cnt:9d}  {rank:4d}  {diff_gc:9.1f}°{marker}")


# ═══════════════════════════════════════════════════════════════════
# STEP E: Peak identification
# ═══════════════════════════════════════════════════════════════════

print("\n--- E. PEAK IDENTIFICATION ---")

# Find peaks (bins higher than both neighbors, wrapping)
peaks = []
for i in range(len(BIN_COUNTS)):
    left = BIN_COUNTS[(i-1) % len(BIN_COUNTS)]
    right = BIN_COUNTS[(i+1) % len(BIN_COUNTS)]
    if BIN_COUNTS[i] > left and BIN_COUNTS[i] > right:
        peaks.append((BINS_CENTER[i], BIN_COUNTS[i]))

print(f"  Peaks (bins higher than both neighbors):")
for center, count in sorted(peaks, key=lambda x: -x[1]):
    diff_gc = abs(center - GC_BEARING_AXIAL)
    diff_gc = min(diff_gc, 180 - diff_gc)
    print(f"    {center}°: {count} lines  (Δ from GC: {diff_gc:.1f}°)")

# Troughs
troughs = []
for i in range(len(BIN_COUNTS)):
    left = BIN_COUNTS[(i-1) % len(BIN_COUNTS)]
    right = BIN_COUNTS[(i+1) % len(BIN_COUNTS)]
    if BIN_COUNTS[i] < left and BIN_COUNTS[i] < right:
        troughs.append((BINS_CENTER[i], BIN_COUNTS[i]))

print(f"\n  Troughs:")
for center, count in sorted(troughs, key=lambda x: x[1]):
    print(f"    {center}°: {count} lines")


# ═══════════════════════════════════════════════════════════════════
# SAVE RESULTS
# ═══════════════════════════════════════════════════════════════════

results = {
    "source": "Richter et al. 2021, Applied Sciences 11(4):1637",
    "n_lines": N_LINES,
    "gc_bearing_at_nazca": GC_BEARING,
    "gc_bearing_axial": GC_BEARING_AXIAL,
    "histogram": {
        "bin_centers": BINS_CENTER,
        "bin_counts": BIN_COUNTS,
        "bin_width": BIN_WIDTH,
    },
    "gc_bin": {
        "bin_center": gc_bin_center,
        "count": gc_bin_count,
        "uniform_expected": round(expected_per_bin, 1),
        "ratio_to_expected": round(gc_bin_count / expected_per_bin, 3),
        "rank_of_18": rank,
    },
    "chi_squared": {
        "statistic": round(chi2_total, 2),
        "df": len(BIN_COUNTS) - 1,
        "p_value": round(float(chi2_p), 6),
        "is_uniform": bool(chi2_p > 0.05),
    },
    "rayleigh_test": {
        "mean_resultant_length": round(float(R_bar), 4),
        "mean_direction_axial_deg": round(float(mean_direction), 2),
        "Z_statistic": round(float(Z), 2),
        "p_value": round(float(rayleigh_p), 6),
        "significant": bool(rayleigh_p < 0.05),
    },
    "v_test": {
        "test_direction_deg": GC_BEARING_AXIAL,
        "V": round(float(V), 4),
        "u": round(float(u), 2),
        "p_value": round(float(v_p), 6),
        "significant": bool(v_p < 0.05),
    },
    "monte_carlo": {
        "n_trials": N_MC,
        "gc_v_statistic": round(gc_v, 4),
        "percentile_rank": round(percentile, 1),
        "bin_density_pct_higher": round(pct_higher_density, 1),
    },
    "peaks": [{"center_deg": p[0], "count": p[1]} for p in sorted(peaks, key=lambda x: -x[1])],
    "troughs": [{"center_deg": t[0], "count": t[1]} for t in sorted(troughs, key=lambda x: x[1])],
}

with open("orientation_analysis.json", "w") as f:
    json.dump(results, f, indent=2)

print(f"\nSaved to orientation_analysis.json")
