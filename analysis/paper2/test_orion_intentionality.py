#!/usr/bin/env python3
"""
Test Claim 2: Orion Intentionality
Is the Giza-Orion shape match intentional, or an accident of near-collinear placement?

Three sub-tests:
A. Random arrangement test (are random triangles on the plateau as good?)
B. Alternative asterism test (do other star triplets match as well?)
C. Epoch sensitivity with error bars (is any epoch uniquely good?)
"""
import numpy as np
from scipy.spatial import procrustes
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from itertools import combinations
import time

OUT = Path('/Users/elliotallan/megalith_site_research/outputs/ambiguous_claims')
OUT.mkdir(parents=True, exist_ok=True)
np.random.seed(42)

print("=" * 70)
print("CLAIM 2: ORION INTENTIONALITY TEST")
print("=" * 70)

# ===================================================================
# Data: Pyramid positions (WGS84)
# ===================================================================
# Actual GPS coordinates
pyramids_wgs84 = {
    'khufu':   (31.1342, 29.9792),
    'khafre':  (31.1308, 29.9761),
    'menkaure': (31.1281, 29.9725),
}

# Convert to local planar coordinates (degrees, relative to Khafre)
ref = pyramids_wgs84['khafre']
pyramids_xy = np.array([
    [pyramids_wgs84['khufu'][0] - ref[0], pyramids_wgs84['khufu'][1] - ref[1]],
    [0.0, 0.0],
    [pyramids_wgs84['menkaure'][0] - ref[0], pyramids_wgs84['menkaure'][1] - ref[1]],
])

# Scale to common units (doesn't matter for Procrustes, but useful for display)
# 1 degree ≈ 111km at this latitude, so these are ~0.006 degree = ~660m offsets

# ===================================================================
# Orion's Belt star positions (J2000 equatorial, converted to planar)
# ===================================================================
# RA (hours), Dec (degrees)
orion_stars_j2000 = {
    'alnitak':  (5.679, -1.943),   # ζ Ori
    'alnilam':  (5.603, -1.202),   # ε Ori
    'mintaka':  (5.533, -0.299),   # δ Ori
}

# Convert to planar: RA hours → degrees (×15), relative to Alnilam
ref_star = orion_stars_j2000['alnilam']
orion_xy = np.array([
    [(orion_stars_j2000['alnitak'][0] - ref_star[0]) * 15 * np.cos(np.radians(ref_star[1])),
     orion_stars_j2000['alnitak'][1] - ref_star[1]],
    [0.0, 0.0],
    [(orion_stars_j2000['mintaka'][0] - ref_star[0]) * 15 * np.cos(np.radians(ref_star[1])),
     orion_stars_j2000['mintaka'][1] - ref_star[1]],
])

# Actual Procrustes distance
_, _, actual_d = procrustes(orion_xy, pyramids_xy)
print(f"\nActual Procrustes distance (Orion→Giza): {actual_d:.6f}")

# ===================================================================
# METHOD A: Random Arrangement Test
# ===================================================================
print("\n[A] Random Arrangement Test (10,000 random triangles on plateau)...")
t0 = time.time()

# Plateau constraints: ~2km × 1km rectangle, minimum 200m spacing
# In degrees: ~0.018° × 0.009°, min spacing ~0.002°
plateau_size = np.array([0.018, 0.009])
min_spacing_deg = 0.002

N_RANDOM = 10000
random_distances = []

for _ in range(N_RANDOM):
    while True:
        pts = np.random.uniform([0, 0], plateau_size, size=(3, 2))
        # Center relative to middle point
        pts_centered = pts - pts[1]
        dists = [np.linalg.norm(pts[i] - pts[j]) for i, j in [(0,1), (0,2), (1,2)]]
        if min(dists) > min_spacing_deg:
            break
    _, _, d = procrustes(orion_xy, pts_centered)
    random_distances.append(d)

random_distances = np.array(random_distances)
p_value_random = np.mean(random_distances <= actual_d)
print(f"   Done in {time.time()-t0:.1f}s")
print(f"   p-value (random ≤ actual): {p_value_random:.6f}")
print(f"   Actual d = {actual_d:.6f}, Random median = {np.median(random_distances):.6f}")

# ===================================================================
# But also test against COLLINEARITY-MATCHED random triangles
# ===================================================================
print("\n[A2] Collinearity-matched random triangles...")

# Measure collinearity of the pyramids
def collinearity(pts):
    """Area of triangle formed by 3 points, normalized by max side length squared."""
    v1 = pts[1] - pts[0]
    v2 = pts[2] - pts[0]
    area = abs(np.cross(v1, v2)) / 2
    max_side = max(np.linalg.norm(pts[i]-pts[j]) for i,j in [(0,1),(0,2),(1,2)])
    return area / (max_side**2) if max_side > 0 else 0

pyramid_collinearity = collinearity(pyramids_xy)
print(f"   Pyramid collinearity metric: {pyramid_collinearity:.6f}")

# Generate collinearity-matched random triangles
col_matched_distances = []
col_tolerance = 0.3  # within 30% of pyramid collinearity
attempts = 0
while len(col_matched_distances) < N_RANDOM:
    attempts += 1
    pts = np.random.uniform([0, 0], plateau_size, size=(3, 2))
    pts_centered = pts - pts[1]
    dists = [np.linalg.norm(pts[i] - pts[j]) for i, j in [(0,1), (0,2), (1,2)]]
    if min(dists) < min_spacing_deg:
        continue
    c = collinearity(pts_centered)
    if abs(c - pyramid_collinearity) / pyramid_collinearity < col_tolerance:
        _, _, d = procrustes(orion_xy, pts_centered)
        col_matched_distances.append(d)

col_matched_distances = np.array(col_matched_distances)
p_value_colmatched = np.mean(col_matched_distances <= actual_d)
print(f"   Generated {N_RANDOM} collinearity-matched triangles (from {attempts} attempts)")
print(f"   p-value (col-matched ≤ actual): {p_value_colmatched:.6f}")

# ===================================================================
# METHOD B: Alternative Asterism Test
# ===================================================================
print("\n[B] Alternative Asterism Test...")

# Bright star triplets (RA hours, Dec degrees, J2000)
# Using the 30 brightest stars and testing all triplet combinations
bright_stars = {
    'Sirius':       (6.752, -16.716),
    'Canopus':      (6.399, -52.696),
    'Arcturus':     (14.261, 19.182),
    'Vega':         (18.616, 38.784),
    'Capella':      (5.278, 45.998),
    'Rigel':        (5.242, -8.202),
    'Procyon':      (7.655, 5.225),
    'Betelgeuse':   (5.919, 7.407),
    'Altair':       (19.846, 8.868),
    'Aldebaran':    (4.599, 16.509),
    'Spica':        (13.420, -11.161),
    'Antares':      (16.490, -26.432),
    'Pollux':       (7.755, 28.026),
    'Fomalhaut':    (22.961, -29.622),
    'Deneb':        (20.690, 45.280),
    'Regulus':      (10.139, 11.967),
    'Castor':       (7.577, 31.888),
    'Bellatrix':    (5.419, 6.350),
    'Alnath':       (5.438, 28.608),
    'Alhena':       (6.629, 16.399),
}

# Named asterisms (traditional groupings)
named_asterisms = {
    'Orion Belt': ['alnitak', 'alnilam', 'mintaka'],  # already tested
    'Summer Triangle': ['Vega', 'Deneb', 'Altair'],
    'Winter Triangle': ['Sirius', 'Betelgeuse', 'Procyon'],
}

def star_to_xy(stars_dict, star_names):
    """Convert 3 star positions to normalized planar coordinates."""
    positions = []
    for name in star_names:
        if name in stars_dict:
            ra, dec = stars_dict[name]
        elif name.lower() in {k.lower(): k for k in stars_dict}:
            key = {k.lower(): k for k in stars_dict}[name.lower()]
            ra, dec = stars_dict[key]
        else:
            return None
        positions.append([ra * 15 * np.cos(np.radians(dec)), dec])
    positions = np.array(positions)
    # Center on middle star
    return positions - positions[1]

# Test named asterisms
asterism_results = {}

# Add Orion Belt
asterism_results['Orion Belt'] = actual_d

# Named ones
for name, stars in named_asterisms.items():
    if name == 'Orion Belt':
        continue
    xy = star_to_xy(bright_stars, stars)
    if xy is not None:
        _, _, d = procrustes(xy, pyramids_xy)
        asterism_results[name] = d

# Test ALL triplet combinations of the 20 brightest stars
star_names = list(bright_stars.keys())
all_triplet_distances = []
best_triplets = []

for combo in combinations(range(len(star_names)), 3):
    names = [star_names[i] for i in combo]
    xy = star_to_xy(bright_stars, names)
    if xy is not None:
        # Try all 3 possible "middle star" orderings
        for order in [[0,1,2], [1,0,2], [0,2,1]]:
            reordered = xy[order]
            reordered_centered = reordered - reordered[1]
            _, _, d = procrustes(reordered_centered, pyramids_xy)
            all_triplet_distances.append(d)
            best_triplets.append(('-'.join([names[i] for i in order]), d))

all_triplet_distances = np.array(all_triplet_distances)
best_triplets.sort(key=lambda x: x[1])

# Where does Orion rank?
orion_rank = np.sum(all_triplet_distances < actual_d) + 1
total_triplets = len(all_triplet_distances)

print(f"   Tested {total_triplets} star triplets (all combos of 20 brightest stars)")
print(f"   Orion Belt Procrustes distance: {actual_d:.6f}")
print(f"   Orion Belt rank: #{orion_rank} out of {total_triplets}")
print(f"   Top 10 matching triplets:")
for name, d in best_triplets[:10]:
    marker = " ← ORION" if 'alnitak' in name.lower() or actual_d == d else ""
    print(f"      {name}: {d:.6f}{marker}")

asterism_results['Best non-Orion'] = best_triplets[0][1] if best_triplets[0][1] != actual_d else best_triplets[1][1]

# ===================================================================
# METHOD C: Epoch Sensitivity
# ===================================================================
print("\n[C] Epoch Sensitivity (precession sweep)...")

# Proper motion of Orion's Belt stars (mas/yr from Hipparcos)
proper_motions = {
    'alnitak':  {'pm_ra': 3.19, 'pm_dec': 2.03},    # mas/yr
    'alnilam':  {'pm_ra': 1.49, 'pm_dec': -1.06},
    'mintaka':  {'pm_ra': 0.56, 'pm_dec': -0.69},
}

# Also need precession: ~50.3"/yr in RA
# For Orion's Belt, the dominant effect over 10,000 years is precession
# But the SHAPE of the triangle changes only due to proper motion
# Precession moves the whole pattern uniformly

# Compute triangle shape at different epochs
epochs = np.arange(-15000, 3001, 50)
epoch_distances = []

for epoch in epochs:
    dt = epoch - 2000  # years from J2000
    positions = {}
    for star, (ra_h, dec_d) in orion_stars_j2000.items():
        pm = proper_motions[star]
        # Proper motion in degrees
        new_ra_h = ra_h + (pm['pm_ra'] / 1000 / 3600 / 15) * dt  # mas → hours
        new_dec_d = dec_d + (pm['pm_dec'] / 1000 / 3600) * dt      # mas → degrees
        positions[star] = (new_ra_h, new_dec_d)

    # Build triangle
    ref_s = positions['alnilam']
    triangle = np.array([
        [(positions['alnitak'][0] - ref_s[0]) * 15 * np.cos(np.radians(ref_s[1])),
         positions['alnitak'][1] - ref_s[1]],
        [0.0, 0.0],
        [(positions['mintaka'][0] - ref_s[0]) * 15 * np.cos(np.radians(ref_s[1])),
         positions['mintaka'][1] - ref_s[1]],
    ])

    _, _, d = procrustes(triangle, pyramids_xy)
    epoch_distances.append(d)

epoch_distances = np.array(epoch_distances)

# Find best epoch
best_idx = np.argmin(epoch_distances)
best_epoch = epochs[best_idx]
best_d = epoch_distances[best_idx]
worst_d = np.max(epoch_distances)

# Total variation
total_variation = worst_d - best_d

print(f"   Epochs tested: {len(epochs)} ({epochs[0]} to {epochs[-1]} CE)")
print(f"   Best epoch: {best_epoch} CE (d = {best_d:.6f})")
print(f"   Worst epoch: {epochs[np.argmax(epoch_distances)]} CE (d = {worst_d:.6f})")
print(f"   Total variation: {total_variation:.6f}")
print(f"   Variation as % of mean: {total_variation/np.mean(epoch_distances)*100:.2f}%")

# At 2560 BCE (traditional construction date)
idx_2560 = np.argmin(np.abs(epochs - (-2560)))
d_2560 = epoch_distances[idx_2560]
print(f"   At 2560 BCE (construction): d = {d_2560:.6f}")

# At 10500 BCE (Hancock's claim)
idx_10500 = np.argmin(np.abs(epochs - (-10500)))
d_10500 = epoch_distances[idx_10500]
print(f"   At 10,500 BCE (Hancock): d = {d_10500:.6f}")

# Is the variation statistically meaningful?
# Compare to bootstrap uncertainty of pyramid positions (~50m = 0.00045°)
survey_uncertainty_deg = 0.00045
# Estimate Procrustes uncertainty via perturbation
n_boot = 1000
boot_distances = []
for _ in range(n_boot):
    perturbed = pyramids_xy + np.random.normal(0, survey_uncertainty_deg, pyramids_xy.shape)
    _, _, d = procrustes(orion_xy, perturbed)
    boot_distances.append(d)
boot_std = np.std(boot_distances)
print(f"   Procrustes uncertainty from survey error: ±{boot_std:.6f}")
print(f"   Total epoch variation / survey uncertainty: {total_variation/boot_std:.2f}×")

# ===================================================================
# PLOTS
# ===================================================================

# Plot 1: Random arrangement histogram
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
ax.hist(random_distances, bins=50, alpha=0.5, color='blue', label='Unconstrained random', density=True)
ax.hist(col_matched_distances, bins=50, alpha=0.5, color='orange', label='Collinearity-matched', density=True)
ax.axvline(actual_d, color='red', linewidth=2, label=f'Actual Giza-Orion (d={actual_d:.4f})')
ax.set_xlabel('Procrustes Distance')
ax.set_ylabel('Density')
ax.set_title('Random Triangle Arrangements vs Giza-Orion Match')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 2: Epoch sweep
ax = axes[1]
ax.plot(epochs, epoch_distances, 'b-', linewidth=1)
ax.axvline(-10500, color='red', linestyle='--', alpha=0.7, label='10,500 BCE (Hancock)')
ax.axvline(-2560, color='green', linestyle='--', alpha=0.7, label='2560 BCE (construction)')
ax.axhline(best_d, color='gray', linestyle=':', alpha=0.5, label=f'Best: {best_epoch} CE')
ax.fill_between(epochs, epoch_distances - boot_std, epoch_distances + boot_std, alpha=0.2, color='blue', label='±1σ survey uncertainty')
ax.set_xlabel('Epoch (CE)')
ax.set_ylabel('Procrustes Distance')
ax.set_title('Giza-Orion Match Quality Across Epochs')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUT / 'orion_random_arrangements.png', dpi=150)
plt.close()

# Plot 3: Asterism comparison
fig, ax = plt.subplots(figsize=(12, 6))
top_n = min(30, len(best_triplets))
names = [bt[0] for bt in best_triplets[:top_n]]
dists = [bt[1] for bt in best_triplets[:top_n]]
colors = ['red' if actual_d == d else 'steelblue' for d in dists]
# Mark Orion
for i, (n, d) in enumerate(best_triplets[:top_n]):
    if abs(d - actual_d) < 1e-10:
        colors[i] = 'red'
        names[i] = names[i] + ' ★'

ax.barh(range(top_n), dists, color=colors)
ax.set_yticks(range(top_n))
ax.set_yticklabels(names, fontsize=7)
ax.set_xlabel('Procrustes Distance (lower = better match)')
ax.set_title(f'Top {top_n} Star Triplet Matches to Giza Layout')
ax.invert_yaxis()
ax.grid(True, alpha=0.3, axis='x')
plt.tight_layout()
plt.savefig(OUT / 'orion_asterism_comparison.png', dpi=150)
plt.close()

# ===================================================================
# Save results
# ===================================================================
results = {
    'actual_procrustes_distance': round(float(actual_d), 6),
    'random_arrangement_test': {
        'n_random': N_RANDOM,
        'p_value_unconstrained': round(float(p_value_random), 6),
        'p_value_collinearity_matched': round(float(p_value_colmatched), 6),
        'random_median': round(float(np.median(random_distances)), 6),
        'colmatched_median': round(float(np.median(col_matched_distances)), 6),
    },
    'asterism_comparison': {
        'total_triplets_tested': total_triplets,
        'orion_rank': int(orion_rank),
        'top_10': [{'name': n, 'distance': round(float(d), 6)} for n, d in best_triplets[:10]],
    },
    'epoch_sensitivity': {
        'best_epoch': int(best_epoch),
        'best_distance': round(float(best_d), 6),
        'worst_distance': round(float(worst_d), 6),
        'total_variation': round(float(total_variation), 6),
        'variation_pct_of_mean': round(float(total_variation/np.mean(epoch_distances)*100), 2),
        'at_2560_bce': round(float(d_2560), 6),
        'at_10500_bce': round(float(d_10500), 6),
        'survey_uncertainty': round(float(boot_std), 6),
        'variation_over_uncertainty': round(float(total_variation/boot_std), 2),
    },
    'pyramid_collinearity': round(float(pyramid_collinearity), 6),
}

with open(OUT / 'orion_intentionality.json', 'w') as f:
    json.dump(results, f, indent=2)

# ===================================================================
# Verdict
# ===================================================================
print("\n" + "=" * 70)
print("VERDICT")
print("=" * 70)

# Determine verdict components
shape_match_real = p_value_colmatched < 0.05
orion_unique = orion_rank <= max(5, total_triplets * 0.01)  # top 1%
epoch_meaningful = total_variation > 3 * boot_std

if shape_match_real and orion_unique:
    shape_verdict = "SUPPORTED — the shape match is real AND Orion is among the best-matching asterisms"
elif shape_match_real and not orion_unique:
    shape_verdict = "MIXED — the shape match is real but many other star triplets match equally well"
else:
    shape_verdict = "NOT SUPPORTED — random triangles on the plateau match Orion just as well"

if epoch_meaningful:
    epoch_verdict = f"Epoch variation IS meaningful (variation = {total_variation/boot_std:.1f}× survey uncertainty). Best epoch: {best_epoch} CE."
else:
    epoch_verdict = f"Epoch variation is NOT meaningful (variation = {total_variation/boot_std:.1f}× survey uncertainty). All epochs produce essentially the same match."

# Overall
if shape_match_real and orion_unique and not epoch_meaningful:
    overall = "PARTIALLY SUPPORTED"
    summary = ("The Giza layout does match Orion's Belt significantly better than random plateau arrangements, "
              "and Orion ranks highly among tested asterisms. However, the match quality barely varies across "
              "epochs, so the claim that 10,500 BCE is a special date is NOT supported by the geometry alone.")
elif shape_match_real and not orion_unique:
    overall = "GEOMETRIC COINCIDENCE"
    summary = ("The shape match is better than random, but many other star triplets produce equal or better "
              "matches. The near-collinear geometry of both Orion's Belt and the Giza pyramids makes this "
              "correspondence expected rather than remarkable.")
elif not shape_match_real:
    overall = "NOT SUPPORTED"
    summary = "Random triangles on the plateau match Orion as well as the actual pyramid layout."
else:
    overall = "SUPPORTED"
    summary = "All three conditions met — shape match, Orion uniqueness, and epoch sensitivity."

print(f"  Shape: {shape_verdict}")
print(f"  Epoch: {epoch_verdict}")
print(f"  OVERALL: {overall}")
print(f"  {summary}")

verdict_md = f"""# Orion Intentionality — Verdict

## Claim
The three Giza pyramids were intentionally arranged to mirror Orion's Belt,
specifically as it appeared at 10,500 BCE.

## Test Results

### A. Random Arrangement Test
| Test | p-value | Interpretation |
|---|---|---|
| Unconstrained random triangles | {p_value_random:.6f} | {'Significant' if p_value_random < 0.05 else 'Not significant'} |
| Collinearity-matched triangles | {p_value_colmatched:.6f} | {'Significant' if p_value_colmatched < 0.05 else 'Not significant'} |

The collinearity-matched test is the fair comparison: it asks whether the specific shape
(not just the near-collinearity) matches Orion unusually well.

### B. Alternative Asterism Test
- Tested {total_triplets} star triplets (all combinations of 20 brightest stars)
- Orion Belt rank: **#{orion_rank}** out of {total_triplets}
- Top 5 matches:
{chr(10).join(f'  {i+1}. {n}: {d:.6f}' for i, (n, d) in enumerate(best_triplets[:5]))}

### C. Epoch Sensitivity
| Metric | Value |
|---|---|
| Best matching epoch | {best_epoch} CE (d = {best_d:.6f}) |
| At 2560 BCE (construction) | d = {d_2560:.6f} |
| At 10,500 BCE (Hancock) | d = {d_10500:.6f} |
| Total variation across all epochs | {total_variation:.6f} |
| Survey measurement uncertainty | ±{boot_std:.6f} |
| Variation / uncertainty ratio | {total_variation/boot_std:.1f}× |

{'The epoch variation is smaller than or comparable to measurement uncertainty — ALL epochs produce essentially the same match. The 10,500 BCE date claim has no geometric support.' if not epoch_meaningful else f'The epoch variation exceeds measurement uncertainty. The best-matching epoch is {best_epoch} CE.'}

## Verdict: **{overall}**

{summary}

### Key Findings
1. **Shape match quality**: Procrustes distance = {actual_d:.6f} (lower = better)
2. **vs unconstrained random**: p = {p_value_random:.6f}
3. **vs collinearity-matched random**: p = {p_value_colmatched:.6f}
4. **Orion's rank among all bright star triplets**: #{orion_rank}/{total_triplets}
5. **Epoch sensitivity**: Variation of {total_variation:.6f} across 18,000 years
   (survey uncertainty = {boot_std:.6f})
"""

with open(OUT / 'ORION_VERDICT.md', 'w') as f:
    f.write(verdict_md)

print(f"\nSaved: orion_intentionality.json, orion_random_arrangements.png, orion_asterism_comparison.png, ORION_VERDICT.md")
