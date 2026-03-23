#!/usr/bin/env python3
"""
Test Claim 3: Pillar 43 at Göbekli Tepe
Sweatman & Tsikritsis (2017) claim animal carvings encode stellar positions
corresponding to the summer solstice sky at ~10,950 BCE.

Three sub-tests:
A. Temporal sweep (is 10,950 BCE the best-matching epoch?)
B. Random arrangement test (do random arrangements match as well?)
C. Alternative identification test (is the claimed animal→constellation mapping optimal?)
"""
import numpy as np
from scipy.spatial import procrustes
from itertools import permutations
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import time
import warnings
warnings.filterwarnings('ignore')

OUT = Path('/Users/elliotallan/megalith_site_research/outputs/ambiguous_claims')
OUT.mkdir(parents=True, exist_ok=True)
np.random.seed(42)

print("=" * 70)
print("CLAIM 3: PILLAR 43 STELLAR ENCODING TEST")
print("=" * 70)

# ===================================================================
# Data: Carving positions from Sweatman & Tsikritsis (2017)
# Approximate positions extracted from their Figure 2
# (x, y in relative units on the pillar face, 0-1 range)
# ===================================================================
carvings = {
    'vulture':           {'x': 0.45, 'y': 0.65, 'claimed': 'Sagittarius'},
    'scorpion':          {'x': 0.55, 'y': 0.25, 'claimed': 'Scorpius'},
    'bird_above_scorp':  {'x': 0.50, 'y': 0.45, 'claimed': 'Libra'},
    'circle_disk':       {'x': 0.48, 'y': 0.55, 'claimed': 'Sun_solstice'},
    'fox':               {'x': 0.70, 'y': 0.50, 'claimed': None},
    'serpent':           {'x': 0.30, 'y': 0.40, 'claimed': None},
    'ibis_like_bird':    {'x': 0.35, 'y': 0.70, 'claimed': None},
    'headless_man':      {'x': 0.60, 'y': 0.35, 'claimed': None},
}

# Only carvings with claimed stellar identifications
claimed_carvings = {k: v for k, v in carvings.items() if v['claimed'] is not None}
carving_positions = np.array([[v['x'], v['y']] for v in claimed_carvings.values()])
carving_names = list(claimed_carvings.keys())
constellation_names = [v['claimed'] for v in claimed_carvings.values()]

print(f"\nCarvings with claimed identifications: {len(claimed_carvings)}")
for name, info in claimed_carvings.items():
    print(f"   {name} → {info['claimed']} at ({info['x']:.2f}, {info['y']:.2f})")

# ===================================================================
# Göbekli Tepe location
# ===================================================================
GT_LAT = 37.223  # degrees N
GT_LON = 38.922  # degrees E

# ===================================================================
# Stellar positions using astropy
# ===================================================================
try:
    from astropy.coordinates import SkyCoord, AltAz, EarthLocation, FK5, Galactic
    from astropy.time import Time
    import astropy.units as u
    HAS_ASTROPY = True
except ImportError:
    HAS_ASTROPY = False
    print("WARNING: astropy not available, using simplified precession model")

# Reference star positions (J2000 equatorial)
# Using brightest star of each claimed constellation as proxy
reference_stars = {
    'Sagittarius': {'ra': 283.816, 'dec': -26.296, 'name': 'Kaus Australis (ε Sgr)'},
    'Scorpius':    {'ra': 247.352, 'dec': -26.432, 'name': 'Antares (α Sco)'},
    'Libra':       {'ra': 229.252, 'dec': -9.383,  'name': 'Zubeneschamali (β Lib)'},
    'Sun_solstice': {'ra': 90.0,   'dec': 23.44,   'name': 'Sun at June solstice (approx)'},
}

# For the Sun position: at summer solstice, the Sun is at RA ≈ 6h = 90°, Dec ≈ +23.44°
# But this changes with precession! The ecliptic longitude is always 90° (by definition),
# but the RA/Dec changes due to precession of the equinoxes.

def get_star_positions_at_epoch(epoch_year, location_lat, location_lon):
    """
    Get projected sky positions of reference stars at summer solstice midnight
    for a given epoch, as seen from Göbekli Tepe.

    Uses precession of equatorial coordinates. For Procrustes comparison,
    we only need RELATIVE positions, so a simplified precession model
    (linear RA shift + obliquity-corrected Dec) is sufficient.
    """
    dt = epoch_year - 2000  # years from J2000
    precession_rate = 50.3 / 3600  # degrees per year in RA

    # Obliquity changes slowly: ~23.44° - 0.013°/century
    obliquity = 23.44 - 0.013 * dt / 100

    # Local sidereal time at midnight on June 21
    # LST ≈ RA of Sun + 12h (Sun is at RA~6h at solstice in J2000)
    # Sun's RA at solstice precesses with the equinox
    sun_ra_solstice = 90.0 + precession_rate * dt  # degrees
    lst = (sun_ra_solstice + 180) % 360  # midnight is 12h from Sun

    positions = {}
    for const_name, star_info in reference_stars.items():
        if const_name == 'Sun_solstice':
            # Sun at solstice: ecliptic lon=90°, convert to equatorial at epoch
            # RA = sun_ra_solstice, Dec ≈ obliquity
            ra = sun_ra_solstice
            dec = obliquity
        else:
            # Precess star: RA shifts, Dec changes slowly due to obliquity
            # More accurate: proper precession in ecliptic coords
            # But for relative positions, linear RA shift is dominant
            ra = star_info['ra'] + precession_rate * dt
            dec = star_info['dec']  # Dec change is second-order for these stars

        # Convert to hour angle
        ha = lst - ra  # degrees

        # Convert to alt/az
        ha_rad = np.radians(ha)
        dec_rad = np.radians(dec)
        lat_rad = np.radians(location_lat)

        alt = np.degrees(np.arcsin(
            np.sin(dec_rad) * np.sin(lat_rad) +
            np.cos(dec_rad) * np.cos(lat_rad) * np.cos(ha_rad)
        ))
        az = np.degrees(np.arctan2(
            -np.cos(dec_rad) * np.sin(ha_rad),
            np.sin(dec_rad) * np.cos(lat_rad) - np.cos(dec_rad) * np.sin(lat_rad) * np.cos(ha_rad)
        )) % 360

        positions[const_name] = {'az': az, 'alt': alt}

    return positions


# ===================================================================
# METHOD A: Temporal Sweep
# ===================================================================
print("\n[A] Temporal Sweep (-15000 to -5000 CE, 100-year steps)...")
t0 = time.time()

epochs = np.arange(-15000, -4900, 100)
match_history = []

for epoch in epochs:
    positions = get_star_positions_at_epoch(epoch, GT_LAT, GT_LON)

    # Get stellar positions as 2D array (using az, alt)
    star_xy = np.array([[positions[cn]['az'], positions[cn]['alt']]
                        for cn in constellation_names])

    # Procrustes match to carving positions
    try:
        _, _, d = procrustes(star_xy, carving_positions)
        match_history.append({'epoch': int(epoch), 'procrustes_d': float(d)})
    except Exception as e:
        match_history.append({'epoch': int(epoch), 'procrustes_d': float('nan'), 'error': str(e)})

match_array = np.array([(m['epoch'], m['procrustes_d']) for m in match_history if not np.isnan(m['procrustes_d'])])

if len(match_array) > 0:
    best_idx = np.argmin(match_array[:, 1])
    best_epoch = int(match_array[best_idx, 0])
    best_d = match_array[best_idx, 1]
    worst_d = np.max(match_array[:, 1])
    mean_d = np.mean(match_array[:, 1])
    total_variation = worst_d - best_d
    variation_pct = total_variation / mean_d * 100

    # Find match at claimed epoch (10,950 BCE = -10950)
    idx_claimed = np.argmin(np.abs(match_array[:, 0] - (-10950)))
    d_claimed = match_array[idx_claimed, 1]

    print(f"   Done in {time.time()-t0:.1f}s ({len(match_array)} valid epochs)")
    print(f"   Best epoch: {best_epoch} CE (d = {best_d:.6f})")
    print(f"   Worst epoch: {int(match_array[np.argmax(match_array[:,1]), 0])} CE (d = {worst_d:.6f})")
    print(f"   At claimed epoch (-10950): d = {d_claimed:.6f}")
    print(f"   Total variation: {total_variation:.6f} ({variation_pct:.1f}% of mean)")
    print(f"   Is -10950 within 10% of best? {abs(d_claimed - best_d) / total_variation < 0.1 if total_variation > 0 else 'N/A'}")
else:
    print("   ERROR: No valid epoch calculations")
    best_epoch = -10950
    best_d = 0
    d_claimed = 0
    total_variation = 0
    variation_pct = 0
    worst_d = 0
    mean_d = 0

# ===================================================================
# METHOD B: Random Arrangement Test
# ===================================================================
print("\n[B] Random Arrangement Test (10,000 random carvings)...")
t0 = time.time()

# Get stellar positions at the claimed epoch
positions_claimed = get_star_positions_at_epoch(-10950, GT_LAT, GT_LON)
star_xy_claimed = np.array([[positions_claimed[cn]['az'], positions_claimed[cn]['alt']]
                            for cn in constellation_names])

N_RANDOM = 10000
random_distances = []
n_carvings = len(claimed_carvings)

for _ in range(N_RANDOM):
    # Random positions on a pillar face (unit square)
    random_pos = np.random.uniform(0, 1, size=(n_carvings, 2))
    try:
        _, _, d = procrustes(star_xy_claimed, random_pos)
        random_distances.append(d)
    except:
        pass

random_distances = np.array(random_distances)
p_value_random = np.mean(random_distances <= d_claimed)
print(f"   Done in {time.time()-t0:.1f}s")
print(f"   Actual match at -10950: d = {d_claimed:.6f}")
print(f"   Random median: {np.median(random_distances):.6f}")
print(f"   p-value: {p_value_random:.6f}")

# ===================================================================
# Also test: random arrangements × temporal sweep
# For each random arrangement, find its BEST epoch and distance
# ===================================================================
print("\n[B2] Random arrangements with temporal optimization (1000 × epoch sweep)...")
t0 = time.time()

N_RANDOM_SWEEP = 1000
# Pre-compute stellar positions at each epoch
epoch_stellar = {}
for epoch in epochs:
    pos = get_star_positions_at_epoch(epoch, GT_LAT, GT_LON)
    epoch_stellar[epoch] = np.array([[pos[cn]['az'], pos[cn]['alt']]
                                     for cn in constellation_names])

random_best_distances = []
for _ in range(N_RANDOM_SWEEP):
    random_pos = np.random.uniform(0, 1, size=(n_carvings, 2))
    best_d_rand = float('inf')
    for epoch in epochs:
        try:
            _, _, d = procrustes(epoch_stellar[epoch], random_pos)
            if d < best_d_rand:
                best_d_rand = d
        except:
            pass
    random_best_distances.append(best_d_rand)

random_best_distances = np.array(random_best_distances)
# Compare: the actual carving's best epoch match vs random best epoch matches
actual_best_d = best_d  # from the temporal sweep
p_value_sweep = np.mean(random_best_distances <= actual_best_d)
print(f"   Done in {time.time()-t0:.1f}s")
print(f"   Actual best match (across all epochs): d = {actual_best_d:.6f}")
print(f"   Random best match: median = {np.median(random_best_distances):.6f}")
print(f"   p-value (random best ≤ actual best): {p_value_sweep:.6f}")

# ===================================================================
# METHOD C: Alternative Identification Test
# ===================================================================
print("\n[C] Alternative Identification Test (permutation test)...")

# Only 3 constellation-matched carvings (excluding Sun which is the 4th)
# Test all permutations of animal→constellation assignment
# For 4 items: 4! = 24 permutations

n_items = len(constellation_names)
all_perms = list(permutations(range(n_items)))

perm_results = []
for perm in all_perms:
    # Reorder constellation positions according to this permutation
    reordered_stars = star_xy_claimed[list(perm)]
    try:
        _, _, d = procrustes(reordered_stars, carving_positions)
        mapping = {carving_names[i]: constellation_names[perm[i]] for i in range(n_items)}
        perm_results.append({
            'permutation': list(perm),
            'mapping': mapping,
            'distance': float(d),
        })
    except:
        pass

perm_results.sort(key=lambda x: x['distance'])

# Where does the claimed assignment rank?
claimed_perm = tuple(range(n_items))  # identity = claimed mapping
claimed_rank = next(i+1 for i, p in enumerate(perm_results) if tuple(p['permutation']) == claimed_perm)

print(f"   Total permutations tested: {len(perm_results)}")
print(f"   Claimed mapping rank: #{claimed_rank} out of {len(perm_results)}")
print(f"   Best permutation:")
print(f"      Distance: {perm_results[0]['distance']:.6f}")
print(f"      Mapping: {perm_results[0]['mapping']}")
print(f"   Claimed permutation:")
print(f"      Distance: {perm_results[claimed_rank-1]['distance']:.6f}")
print(f"      Mapping: {perm_results[claimed_rank-1]['mapping']}")

# ===================================================================
# METHOD D: Expanded constellation test
# ===================================================================
print("\n[D] Expanded constellation test (testing other constellation assignments)...")

# What if the carvings map to DIFFERENT constellations entirely?
# Test 20 prominent constellations visible from GT at various epochs
alt_constellations = {
    'Sagittarius': {'ra': 283.816, 'dec': -26.296},
    'Scorpius':    {'ra': 247.352, 'dec': -26.432},
    'Libra':       {'ra': 229.252, 'dec': -9.383},
    'Virgo':       {'ra': 201.298, 'dec': -11.161},  # Spica
    'Leo':         {'ra': 152.093, 'dec': 11.967},    # Regulus
    'Cancer':      {'ra': 130.154, 'dec': 21.468},    # Al Tarf
    'Gemini':      {'ra': 116.329, 'dec': 28.026},    # Pollux
    'Taurus':      {'ra': 68.980, 'dec': 16.509},     # Aldebaran
    'Orion':       {'ra': 88.793, 'dec': 7.407},      # Betelgeuse
    'Canis_Major': {'ra': 101.287, 'dec': -16.716},   # Sirius
    'Aquila':      {'ra': 297.696, 'dec': 8.868},     # Altair
    'Lyra':        {'ra': 279.235, 'dec': 38.784},    # Vega
    'Cygnus':      {'ra': 310.358, 'dec': 45.280},    # Deneb
    'Pegasus':     {'ra': 346.190, 'dec': 15.205},    # Enif
    'Capricornus': {'ra': 305.253, 'dec': -14.782},   # Deneb Algedi
    'Aquarius':    {'ra': 331.446, 'dec': -0.320},    # Sadalsuud
    'Pisces':      {'ra': 2.097, 'dec': 29.091},      # Eta Psc
    'Aries':       {'ra': 31.793, 'dec': 23.462},     # Hamal
    'Ophiuchus':   {'ra': 263.734, 'dec': -24.167},   # Yed Prior
    'Centaurus':   {'ra': 219.902, 'dec': -60.835},   # Alpha Cen
}

# For the 3 non-Sun carvings, test all possible assignments from 20 constellations
# C(20,3) × 3! = 1140 × 6 = 6840 combinations (plus the Sun)
from itertools import combinations as combs

# Use positions at -10950
dt = -10950 - 2000
prec_rate = 50.3 / 3600

alt_results = []
const_names_list = list(alt_constellations.keys())

# For the 3 animal carvings (excluding circle/Sun)
animal_carvings = {k: v for k, v in claimed_carvings.items() if v['claimed'] != 'Sun_solstice'}
animal_positions = np.array([[v['x'], v['y']] for v in animal_carvings.values()])
animal_names = list(animal_carvings.keys())

for combo in combs(range(len(const_names_list)), 3):
    combo_names = [const_names_list[i] for i in combo]
    combo_stars = []
    for cn in combo_names:
        info = alt_constellations[cn]
        new_ra = info['ra'] + prec_rate * dt
        combo_stars.append([new_ra % 360, info['dec']])
    combo_stars = np.array(combo_stars)

    for perm in permutations(range(3)):
        reordered = combo_stars[list(perm)]
        try:
            _, _, d = procrustes(reordered, animal_positions)
            mapping = {animal_names[i]: combo_names[perm[i]] for i in range(3)}
            is_claimed = set(combo_names) == {'Sagittarius', 'Scorpius', 'Libra'}
            alt_results.append({
                'constellations': combo_names,
                'mapping': mapping,
                'distance': float(d),
                'is_claimed_set': is_claimed,
            })
        except:
            pass

alt_results.sort(key=lambda x: x['distance'])
# Find rank of the claimed set
claimed_set_results = [r for r in alt_results if r['is_claimed_set']]
claimed_best = min(claimed_set_results, key=lambda x: x['distance']) if claimed_set_results else None
claimed_set_rank = next((i+1 for i, r in enumerate(alt_results) if r['is_claimed_set']), len(alt_results))

print(f"   Tested {len(alt_results)} constellation assignments")
print(f"   Best overall: {alt_results[0]['mapping']} (d={alt_results[0]['distance']:.6f})")
if claimed_best:
    print(f"   Claimed set (Sgr/Sco/Lib) best: rank #{claimed_set_rank} (d={claimed_best['distance']:.6f})")
print(f"   Top 5:")
for i, r in enumerate(alt_results[:5]):
    flag = " ← CLAIMED" if r['is_claimed_set'] else ""
    print(f"      {i+1}. {r['mapping']} d={r['distance']:.6f}{flag}")

# ===================================================================
# PLOTS
# ===================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Temporal sweep
ax = axes[0, 0]
if len(match_array) > 0:
    ax.plot(match_array[:, 0], match_array[:, 1], 'b-', linewidth=1)
    ax.axvline(-10950, color='red', linestyle='--', alpha=0.7, label='10,950 BCE (claimed)')
    ax.axvline(best_epoch, color='green', linestyle=':', alpha=0.7, label=f'Best: {best_epoch} CE')
    ax.set_xlabel('Epoch (CE)')
    ax.set_ylabel('Procrustes Distance')
    ax.set_title('Pillar 43 — Stellar Match vs Epoch')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

# Plot 2: Random arrangement histogram
ax = axes[0, 1]
ax.hist(random_distances, bins=50, alpha=0.7, color='steelblue', density=True)
ax.axvline(d_claimed, color='red', linewidth=2, label=f'Actual (d={d_claimed:.4f})')
ax.set_xlabel('Procrustes Distance')
ax.set_ylabel('Density')
ax.set_title(f'Random Arrangements vs Actual (p={p_value_random:.4f})')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 3: Random best-epoch histogram
ax = axes[1, 0]
ax.hist(random_best_distances, bins=50, alpha=0.7, color='steelblue', density=True)
ax.axvline(actual_best_d, color='red', linewidth=2, label=f'Actual best (d={actual_best_d:.4f})')
ax.set_xlabel('Best Procrustes Distance (across all epochs)')
ax.set_ylabel('Density')
ax.set_title(f'Random Best-Epoch Match vs Actual (p={p_value_sweep:.4f})')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 4: Permutation distances
ax = axes[1, 1]
perm_dists = [p['distance'] for p in perm_results]
ax.bar(range(len(perm_dists)), sorted(perm_dists), color='steelblue', alpha=0.7)
# Highlight claimed
claimed_d_val = perm_results[claimed_rank-1]['distance']
ax.axhline(claimed_d_val, color='red', linestyle='--', label=f'Claimed mapping (rank #{claimed_rank})')
ax.set_xlabel('Permutation (sorted)')
ax.set_ylabel('Procrustes Distance')
ax.set_title(f'All {len(perm_results)} Permutations of Animal→Constellation')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUT / 'pillar43_epoch_plot.png', dpi=150)
plt.close()

# ===================================================================
# Save results
# ===================================================================
results = {
    'carving_positions': {k: {'x': v['x'], 'y': v['y'], 'claimed': v['claimed']}
                         for k, v in carvings.items()},
    'temporal_sweep': {
        'epochs_tested': len(match_array) if len(match_array) > 0 else 0,
        'best_epoch': int(best_epoch),
        'best_distance': round(float(best_d), 6),
        'at_claimed_epoch': round(float(d_claimed), 6),
        'total_variation': round(float(total_variation), 6),
        'variation_pct_of_mean': round(float(variation_pct), 2),
    },
    'random_arrangement_test': {
        'n_random': N_RANDOM,
        'p_value_fixed_epoch': round(float(p_value_random), 6),
        'random_median': round(float(np.median(random_distances)), 6),
    },
    'random_sweep_test': {
        'n_random': N_RANDOM_SWEEP,
        'p_value_best_epoch': round(float(p_value_sweep), 6),
        'random_best_median': round(float(np.median(random_best_distances)), 6),
    },
    'permutation_test': {
        'n_permutations': len(perm_results),
        'claimed_rank': claimed_rank,
        'best_permutation': perm_results[0] if perm_results else None,
        'claimed_permutation': perm_results[claimed_rank-1] if perm_results else None,
    },
    'expanded_constellation_test': {
        'total_assignments_tested': len(alt_results),
        'claimed_set_rank': claimed_set_rank,
        'top_5': alt_results[:5],
    },
}

with open(OUT / 'pillar43_temporal_sweep.json', 'w') as f:
    json.dump(results, f, indent=2)

# ===================================================================
# Verdict
# ===================================================================
print("\n" + "=" * 70)
print("VERDICT")
print("=" * 70)

# Epoch test
epoch_close = len(match_array) > 0 and abs(best_epoch - (-10950)) < 500
if epoch_close:
    epoch_verdict = f"SUPPORTED — best epoch ({best_epoch}) is close to claimed 10,950 BCE"
elif total_variation > 0 and (d_claimed - best_d) / total_variation < 0.2:
    epoch_verdict = f"PARTIALLY — 10,950 BCE is near the best epoch ({best_epoch}), within 20% of total variation"
else:
    epoch_verdict = f"NOT SUPPORTED — best epoch is {best_epoch}, not 10,950 BCE; variation is {'minimal' if variation_pct < 5 else 'significant'}"

# Random test
if p_value_random < 0.01:
    random_verdict = "SUPPORTED — actual arrangement matches significantly better than random (p < 0.01)"
elif p_value_random < 0.05:
    random_verdict = f"MARGINAL — p = {p_value_random:.4f}, weakly significant"
else:
    random_verdict = f"NOT SUPPORTED — random arrangements match just as well (p = {p_value_random:.4f})"

# Permutation test
if claimed_rank == 1:
    perm_verdict = "SUPPORTED — claimed animal→constellation mapping is the BEST permutation"
elif claimed_rank <= 3:
    perm_verdict = f"PARTIALLY — claimed mapping ranks #{claimed_rank}, near-optimal"
else:
    perm_verdict = f"NOT SUPPORTED — claimed mapping ranks #{claimed_rank}/{len(perm_results)}, not optimal"

# Expanded test
if claimed_set_rank <= max(5, len(alt_results) * 0.01):
    expanded_verdict = "SUPPORTED — claimed constellation set ranks in top 1%"
else:
    expanded_verdict = f"NOT SUPPORTED — claimed set ranks #{claimed_set_rank}/{len(alt_results)}"

# Overall
conditions_met = sum([
    epoch_close or (total_variation > 0 and (d_claimed - best_d) / total_variation < 0.2),
    p_value_random < 0.05,
    claimed_rank <= 3,
])

if conditions_met >= 3:
    overall = "SUPPORTED"
    summary = "All conditions met — epoch, arrangement significance, and mapping optimality all support the hypothesis."
elif conditions_met == 2:
    overall = "PARTIALLY SUPPORTED"
    summary = "Some conditions met but not all. The evidence is suggestive but not conclusive."
elif conditions_met == 1:
    overall = "WEAKLY SUPPORTED"
    summary = "Only one condition met. The evidence does not convincingly support the stellar encoding hypothesis."
else:
    overall = "NOT SUPPORTED"
    summary = "No conditions met. The carving arrangement is consistent with chance."

print(f"  Epoch: {epoch_verdict}")
print(f"  Random: {random_verdict}")
print(f"  Permutation: {perm_verdict}")
print(f"  Expanded: {expanded_verdict}")
print(f"  Conditions met: {conditions_met}/3")
print(f"  OVERALL: {overall}")
print(f"  {summary}")

verdict_md = f"""# Pillar 43 Stellar Encoding — Verdict

## Claim
Sweatman & Tsikritsis (2017) claim that animal carvings on Pillar 43 at Göbekli Tepe
encode stellar positions corresponding to the summer solstice sky at ~10,950 BCE.

Key identifications:
- Vulture → Sagittarius
- Scorpion → Scorpius
- Bird above scorpion → Libra
- Circle/disk → Sun at solstice

## Test Results

### A. Temporal Sweep (is 10,950 BCE the best epoch?)
| Metric | Value |
|---|---|
| Best matching epoch | {best_epoch} CE |
| Best distance | {best_d:.6f} |
| Distance at 10,950 BCE | {d_claimed:.6f} |
| Total variation | {total_variation:.6f} ({variation_pct:.1f}% of mean) |

**{epoch_verdict}**

### B. Random Arrangement Test
| Test | p-value | Interpretation |
|---|---|---|
| Fixed epoch (10,950 BCE) | {p_value_random:.6f} | {'Significant' if p_value_random < 0.05 else 'Not significant'} |
| Best epoch (temporal optimization) | {p_value_sweep:.6f} | {'Significant' if p_value_sweep < 0.05 else 'Not significant'} |

**{random_verdict}**

Note: The "best epoch" test is the fairer comparison — it asks whether the actual
carving arrangement, optimized across all epochs, matches better than random arrangements
similarly optimized.

### C. Permutation Test (is the claimed animal→constellation mapping optimal?)
- Tested all {len(perm_results)} permutations of {n_items} items
- Claimed mapping rank: **#{claimed_rank}**
- Best mapping: {perm_results[0]['mapping'] if perm_results else 'N/A'}

**{perm_verdict}**

### D. Expanded Constellation Test
- Tested {len(alt_results)} possible constellation assignments (20 constellations, all combos)
- Claimed set (Sgr/Sco/Lib) rank: **#{claimed_set_rank}**

**{expanded_verdict}**

## Verdict: **{overall}**

{summary}

### Methodological Caveats
1. **Carving positions are approximate** — extracted from published figures, not measured in situ.
   Precise positions could change results, though the random arrangement test partially controls for this.
2. **Constellation identification is subjective** — the animal→constellation mapping is an assumption
   of the hypothesis, not independently derived.
3. **The "Sun" position** — using the solstice Sun as one of the matched points builds in a
   constraint that may artificially improve the match for epochs near the solstice alignment.
4. **Multiple testing** — testing many epochs and finding the best match inflates significance
   if not corrected for. The "best epoch" random test addresses this.
"""

with open(OUT / 'PILLAR43_VERDICT.md', 'w') as f:
    f.write(verdict_md)

print(f"\nSaved: pillar43_temporal_sweep.json, pillar43_epoch_plot.png, PILLAR43_VERDICT.md")
