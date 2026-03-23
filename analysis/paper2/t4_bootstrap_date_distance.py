#!/usr/bin/env python3
"""
Test 4: Bootstrap Date-Distance Significance
=============================================
Given the spatial layout of the Memphis necropolis, how unlikely is r = -0.43 by chance?

Two permutation tests:
  A) Shuffle DATES, keep positions fixed — tests whether temporal ordering
     is non-random with respect to distance from the GC.
  B) Shuffle POSITIONS (distances), keep dates fixed — tests whether the
     spatial arrangement is non-random with respect to temporal ordering.
"""

import json, math, os, sys
import numpy as np
from scipy.stats import pearsonr

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "resolution_tests")
os.makedirs(OUT_DIR, exist_ok=True)

# Load pyramid data from prior analysis
with open(os.path.join(BASE_DIR, "outputs", "pyramid_date_distance", "correlation.json")) as f:
    data = json.load(f)

pyramids = data['pyramids']
dates = np.array([p['date_bce'] for p in pyramids])
distances = np.array([p['distance_km'] for p in pyramids])
n = len(pyramids)

observed_r, observed_p = pearsonr(dates, distances)
print(f"Observed: r = {observed_r:.4f}, parametric p = {observed_p:.4f}, n = {n}")

N_BOOTSTRAP = 10000
rng = np.random.default_rng(seed=42)

# ============================================================
# Test A: Shuffle dates, keep positions fixed
# ============================================================
print(f"\nTest A: Shuffling dates ({N_BOOTSTRAP:,} permutations)...")
count_a = 0
null_rs_a = np.empty(N_BOOTSTRAP)

for i in range(N_BOOTSTRAP):
    shuffled_dates = rng.permutation(dates)
    r, _ = pearsonr(shuffled_dates, distances)
    null_rs_a[i] = r
    if r <= observed_r:
        count_a += 1

p_shuffle_dates = count_a / N_BOOTSTRAP
print(f"  Bootstrap p (shuffle dates) = {p_shuffle_dates:.4f}")
print(f"  ({count_a} of {N_BOOTSTRAP} had r <= {observed_r:.4f})")

# ============================================================
# Test B: Shuffle positions (distances), keep dates fixed
# ============================================================
print(f"\nTest B: Shuffling positions ({N_BOOTSTRAP:,} permutations)...")
count_b = 0
null_rs_b = np.empty(N_BOOTSTRAP)

for i in range(N_BOOTSTRAP):
    shuffled_distances = rng.permutation(distances)
    r, _ = pearsonr(dates, shuffled_distances)
    null_rs_b[i] = r
    if r <= observed_r:
        count_b += 1

p_shuffle_positions = count_b / N_BOOTSTRAP
print(f"  Bootstrap p (shuffle positions) = {p_shuffle_positions:.4f}")
print(f"  ({count_b} of {N_BOOTSTRAP} had r <= {observed_r:.4f})")

# ============================================================
# Summary statistics
# ============================================================
print(f"\n{'='*60}")
print("SUMMARY")
print(f"{'='*60}")
print(f"Observed Pearson r:        {observed_r:.4f}")
print(f"Parametric p:              {observed_p:.4f}")
print(f"Bootstrap p (dates):       {p_shuffle_dates:.4f}")
print(f"Bootstrap p (positions):   {p_shuffle_positions:.4f}")
print(f"Null mean r (dates):       {null_rs_a.mean():.4f} +/- {null_rs_a.std():.4f}")
print(f"Null mean r (positions):   {null_rs_b.mean():.4f} +/- {null_rs_b.std():.4f}")

# ============================================================
# Save JSON
# ============================================================
results = {
    'observed_r': float(observed_r),
    'observed_p_parametric': float(observed_p),
    'n_pyramids': n,
    'n_permutations': N_BOOTSTRAP,
    'test_a_shuffle_dates': {
        'description': 'Keep pyramid positions fixed, shuffle dates randomly',
        'p_value': float(p_shuffle_dates),
        'count_more_extreme': int(count_a),
        'null_mean_r': float(null_rs_a.mean()),
        'null_std_r': float(null_rs_a.std()),
        'null_percentiles': {
            '2.5': float(np.percentile(null_rs_a, 2.5)),
            '5': float(np.percentile(null_rs_a, 5)),
            '50': float(np.percentile(null_rs_a, 50)),
            '95': float(np.percentile(null_rs_a, 95)),
            '97.5': float(np.percentile(null_rs_a, 97.5)),
        }
    },
    'test_b_shuffle_positions': {
        'description': 'Keep dates fixed, shuffle which pyramid is at which location',
        'p_value': float(p_shuffle_positions),
        'count_more_extreme': int(count_b),
        'null_mean_r': float(null_rs_b.mean()),
        'null_std_r': float(null_rs_b.std()),
        'null_percentiles': {
            '2.5': float(np.percentile(null_rs_b, 2.5)),
            '5': float(np.percentile(null_rs_b, 5)),
            '50': float(np.percentile(null_rs_b, 50)),
            '95': float(np.percentile(null_rs_b, 95)),
            '97.5': float(np.percentile(null_rs_b, 97.5)),
        }
    },
    'interpretation': (
        f"Shuffling dates: p={p_shuffle_dates:.4f} — "
        f"{'significant' if p_shuffle_dates < 0.05 else 'not significant'} at alpha=0.05. "
        f"Shuffling positions: p={p_shuffle_positions:.4f} — "
        f"{'significant' if p_shuffle_positions < 0.05 else 'not significant'} at alpha=0.05."
    )
}

with open(os.path.join(OUT_DIR, "bootstrap_date_distance.json"), 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved bootstrap_date_distance.json")

# ============================================================
# Histogram plot
# ============================================================
print("Generating histogram...")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Test A histogram
ax1.hist(null_rs_a, bins=80, color='steelblue', alpha=0.7, edgecolor='white', linewidth=0.3)
ax1.axvline(observed_r, color='red', linewidth=2, linestyle='--',
            label=f'Observed r = {observed_r:.3f}')
ax1.axvline(np.percentile(null_rs_a, 5), color='orange', linewidth=1, linestyle=':',
            label='5th percentile')
ax1.set_xlabel('Pearson r (null distribution)')
ax1.set_ylabel('Count')
ax1.set_title(f'Test A: Shuffle Dates (p = {p_shuffle_dates:.4f})')
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)

# Test B histogram
ax2.hist(null_rs_b, bins=80, color='seagreen', alpha=0.7, edgecolor='white', linewidth=0.3)
ax2.axvline(observed_r, color='red', linewidth=2, linestyle='--',
            label=f'Observed r = {observed_r:.3f}')
ax2.axvline(np.percentile(null_rs_b, 5), color='orange', linewidth=1, linestyle=':',
            label='5th percentile')
ax2.set_xlabel('Pearson r (null distribution)')
ax2.set_ylabel('Count')
ax2.set_title(f'Test B: Shuffle Positions (p = {p_shuffle_positions:.4f})')
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

fig.suptitle(f'Bootstrap Significance: Date-Distance Correlation (n={n}, {N_BOOTSTRAP:,} permutations)',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "bootstrap_histogram.png"), dpi=150)
plt.close()
print("Saved bootstrap_histogram.png")

print("\nDone!")
