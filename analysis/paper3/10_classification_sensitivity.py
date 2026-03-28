#!/usr/bin/env python3
"""
Paper 3 Robustness Battery — Task 4: Classification Sensitivity
================================================================
Identify sites with ambiguous monument/settlement classification.
Randomly reassign borderline cases 1,000 times.
Report distribution of divergence D across reassignments.
"""

import csv, json, os, sys, time
import numpy as np

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE = os.path.expanduser("~/megalith_site_research")
OUT = os.path.join(BASE, "paper3/results")
FIG = os.path.join(BASE, "paper3/figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

POLE_LAT = 59.682122
POLE_LON = -138.646087
R = 6371.0
QC = R * np.pi / 2
CORRIDOR = 50.0

MONUMENTAL = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus",
    "shrine", "cemetery"
}
SETTLEMENT = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement", "townhouse", "production"
}

# Types that could reasonably be classified either way
AMBIGUOUS = {
    "fortification", "fort", "fortress", "palace", "citadel",
    "villa", "bath", "bridge", "cistern",
    "stadium", "gymnasium", "forum", "agora",
    "necropolis", "mausoleum", "heroon",
    "altar", "nymphaeum", "odeum",
    "harbor", "wharf", "dock",
    "tower", "gate", "wall",
    "road", "milestone",
    "workshop", "factory", "kiln",
    "granary", "warehouse", "market",
    "camp", "castra"
}

def gc_dist(lats, lons):
    la1 = np.radians(POLE_LAT)
    lo1 = np.radians(POLE_LON)
    la2 = np.radians(lats)
    lo2 = np.radians(lons)
    dl = la2 - la1
    dn = lo2 - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(la2)*np.sin(dn/2)**2
    d = R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QC)

# ============================================================
# LOAD AND CLASSIFY SITES
# ============================================================
print("Loading Pleiades data...")

# Track three categories: clear monument, clear settlement, ambiguous
clear_mon_lats, clear_mon_lons = [], []
clear_set_lats, clear_set_lons = [], []
ambig_lats, ambig_lons = [], []
ambig_types = []  # track what type they are
ambig_originally = []  # 'monument' or 'settlement' or 'neither'

with open(os.path.join(BASE, "data/pleiades/pleiades-places-latest.csv")) as f:
    for row in csv.DictReader(f):
        try:
            lat = float(row['reprLat'])
            lon = float(row['reprLong'])
        except (ValueError, TypeError):
            continue
        if abs(lat) < 0.001 and abs(lon) < 0.001:
            continue

        fts = set(ft.strip() for ft in row['featureTypes'].split(','))
        is_m = bool(fts & MONUMENTAL)
        is_s = bool(fts & SETTLEMENT)
        is_a = bool(fts & AMBIGUOUS)

        if is_m and is_s:
            # Dual-classified — these are inherently ambiguous
            ambig_lats.append(lat)
            ambig_lons.append(lon)
            ambig_types.append(','.join(fts & (MONUMENTAL | SETTLEMENT)))
            ambig_originally.append('both')
        elif is_a and not is_m and not is_s:
            # Has an ambiguous type but not clearly classified
            ambig_lats.append(lat)
            ambig_lons.append(lon)
            ambig_types.append(','.join(fts & AMBIGUOUS))
            ambig_originally.append('neither')
        elif is_m and is_a:
            # Monument with ambiguous features
            ambig_lats.append(lat)
            ambig_lons.append(lon)
            ambig_types.append(','.join(fts & (MONUMENTAL | AMBIGUOUS)))
            ambig_originally.append('monument')
        elif is_s and is_a:
            # Settlement with ambiguous features
            ambig_lats.append(lat)
            ambig_lons.append(lon)
            ambig_types.append(','.join(fts & (SETTLEMENT | AMBIGUOUS)))
            ambig_originally.append('settlement')
        elif is_m and not is_s:
            clear_mon_lats.append(lat)
            clear_mon_lons.append(lon)
        elif is_s and not is_m:
            clear_set_lats.append(lat)
            clear_set_lons.append(lon)

clear_mon_lats = np.array(clear_mon_lats)
clear_mon_lons = np.array(clear_mon_lons)
clear_set_lats = np.array(clear_set_lats)
clear_set_lons = np.array(clear_set_lons)
ambig_lats = np.array(ambig_lats)
ambig_lons = np.array(ambig_lons)

n_clear_mon = len(clear_mon_lats)
n_clear_set = len(clear_set_lats)
n_ambig = len(ambig_lats)
n_total = n_clear_mon + n_clear_set + n_ambig

print(f"Clear monuments: {n_clear_mon}")
print(f"Clear settlements: {n_clear_set}")
print(f"Ambiguous sites: {n_ambig}")
print(f"Ambiguous fraction: {100*n_ambig/n_total:.1f}%")

# Count ambiguous by original classification
from collections import Counter
orig_counts = Counter(ambig_originally)
print(f"\nAmbiguous breakdown:")
for k, v in orig_counts.items():
    print(f"  {k}: {v}")

# ============================================================
# COMPUTE BASELINE D
# ============================================================
print("\n" + "=" * 70)
print("BASELINE DIVERGENCE")
print("=" * 70)

# Baseline: ambiguous sites classified as per their primary type
# 'monument' → monument, 'settlement' → settlement,
# 'both' and 'neither' → excluded (conservative)

def compute_D(mon_lats, mon_lons, set_lats, set_lons):
    """Compute divergence D = mon_frac / set_frac."""
    mon_gc = gc_dist(mon_lats, mon_lons)
    set_gc = gc_dist(set_lats, set_lons)
    mon_frac = np.sum(mon_gc <= CORRIDOR) / len(mon_gc) if len(mon_gc) > 0 else 0
    set_frac = np.sum(set_gc <= CORRIDOR) / len(set_gc) if len(set_gc) > 0 else 0
    return mon_frac / set_frac if set_frac > 0 else float('inf'), mon_frac, set_frac

D_baseline, mf_base, sf_base = compute_D(clear_mon_lats, clear_mon_lons,
                                           clear_set_lats, clear_set_lons)
print(f"D (clear sites only, ambiguous excluded): {D_baseline:.4f}")
print(f"  Mon fraction: {mf_base:.4f}, Set fraction: {sf_base:.4f}")

# Include ambiguous in their original categories
all_mon_lats = list(clear_mon_lats)
all_mon_lons = list(clear_mon_lons)
all_set_lats = list(clear_set_lats)
all_set_lons = list(clear_set_lons)
for i in range(n_ambig):
    if ambig_originally[i] in ('monument', 'both'):
        all_mon_lats.append(ambig_lats[i])
        all_mon_lons.append(ambig_lons[i])
    elif ambig_originally[i] in ('settlement',):
        all_set_lats.append(ambig_lats[i])
        all_set_lons.append(ambig_lons[i])

D_with_ambig, mf_wa, sf_wa = compute_D(np.array(all_mon_lats), np.array(all_mon_lons),
                                         np.array(all_set_lats), np.array(all_set_lons))
print(f"D (ambiguous included as-classified): {D_with_ambig:.4f}")

# ============================================================
# RANDOM REASSIGNMENT (1000 trials)
# ============================================================
print("\n" + "=" * 70)
print("RANDOM REASSIGNMENT — 1000 TRIALS")
print("=" * 70)

N_TRIALS = 1000
rng = np.random.RandomState(42)

# Precompute GC distances for all sites
clear_mon_gc = gc_dist(clear_mon_lats, clear_mon_lons)
clear_set_gc = gc_dist(clear_set_lats, clear_set_lons)
ambig_gc = gc_dist(ambig_lats, ambig_lons)

n_clear_mon_in = np.sum(clear_mon_gc <= CORRIDOR)
n_clear_set_in = np.sum(clear_set_gc <= CORRIDOR)
ambig_in_corridor = ambig_gc <= CORRIDOR

D_trials = np.zeros(N_TRIALS)
t0 = time.time()

for trial in range(N_TRIALS):
    # Randomly assign each ambiguous site to monument or settlement (50/50)
    assignment = rng.randint(0, 2, n_ambig)  # 0 = monument, 1 = settlement

    n_mon = n_clear_mon + np.sum(assignment == 0)
    n_set = n_clear_set + np.sum(assignment == 1)

    n_mon_in = n_clear_mon_in + np.sum((assignment == 0) & ambig_in_corridor)
    n_set_in = n_clear_set_in + np.sum((assignment == 1) & ambig_in_corridor)

    mon_frac = n_mon_in / n_mon if n_mon > 0 else 0
    set_frac = n_set_in / n_set if n_set > 0 else 0

    D_trials[trial] = mon_frac / set_frac if set_frac > 0 else float('inf')

print(f"Completed {N_TRIALS} trials in {time.time()-t0:.1f}s")
print(f"\nD distribution across random reassignments:")
print(f"  Mean:   {np.mean(D_trials):.4f}")
print(f"  Median: {np.median(D_trials):.4f}")
print(f"  Std:    {np.std(D_trials):.4f}")
print(f"  Min:    {np.min(D_trials):.4f}")
print(f"  Max:    {np.max(D_trials):.4f}")
print(f"  95% CI: [{np.percentile(D_trials, 2.5):.4f}, {np.percentile(D_trials, 97.5):.4f}]")

# Key question: does D ever drop below 1.0 (no divergence)?
n_below_1 = np.sum(D_trials < 1.0)
print(f"\n  Trials with D < 1.0: {n_below_1}/{N_TRIALS} ({100*n_below_1/N_TRIALS:.1f}%)")
print(f"  Trials with D > baseline ({D_baseline:.3f}): "
      f"{np.sum(D_trials > D_baseline)}/{N_TRIALS}")

# Coefficient of variation
cv = np.std(D_trials) / np.mean(D_trials)
print(f"  Coefficient of variation: {cv:.4f} ({100*cv:.1f}%)")

if n_below_1 == 0:
    print("\n→ ROBUST: Divergence D > 1.0 in ALL random reassignments")
    print("  Classification ambiguity does not eliminate the signal")
else:
    print(f"\n→ CAUTION: D < 1.0 in {100*n_below_1/N_TRIALS:.1f}% of reassignments")

# ============================================================
# SAVE
# ============================================================
results = {
    'analysis': 'Classification Sensitivity',
    'n_trials': N_TRIALS,
    'site_counts': {
        'clear_monuments': int(n_clear_mon),
        'clear_settlements': int(n_clear_set),
        'ambiguous': int(n_ambig),
        'ambiguous_fraction': float(n_ambig / n_total),
        'ambiguous_breakdown': dict(orig_counts)
    },
    'baseline': {
        'D_clear_only': float(D_baseline),
        'D_with_ambiguous': float(D_with_ambig)
    },
    'reassignment_distribution': {
        'mean': float(np.mean(D_trials)),
        'median': float(np.median(D_trials)),
        'std': float(np.std(D_trials)),
        'min': float(np.min(D_trials)),
        'max': float(np.max(D_trials)),
        'ci_95': [float(np.percentile(D_trials, 2.5)), float(np.percentile(D_trials, 97.5))],
        'cv': float(cv),
        'n_below_1': int(n_below_1),
        'fraction_below_1': float(n_below_1 / N_TRIALS)
    },
    'conclusion': 'Robust' if n_below_1 == 0 else 'Sensitive to classification'
}

with open(os.path.join(OUT, 'classification_sensitivity.json'), 'w') as f:
    json.dump(results, f, indent=2)

# ============================================================
# PLOT
# ============================================================
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.hist(D_trials, bins=50, color='steelblue', alpha=0.7, edgecolor='none')
    ax.axvline(D_baseline, color='red', linewidth=2, linestyle='--',
               label=f'Baseline D={D_baseline:.3f}')
    ax.axvline(1.0, color='black', linewidth=1, linestyle=':',
               label='D=1.0 (no divergence)')
    ax.set_xlabel('Divergence D')
    ax.set_ylabel('Count')
    ax.set_title(f'Classification Sensitivity ({N_TRIALS} reassignments)\n'
                 f'{n_ambig} ambiguous sites ({100*n_ambig/n_total:.1f}%)')
    ax.legend()

    ax2 = axes[1]
    ax2.plot(np.sort(D_trials), np.arange(1, N_TRIALS+1)/N_TRIALS, 'b-')
    ax2.axvline(D_baseline, color='red', linewidth=2, linestyle='--')
    ax2.axvline(1.0, color='black', linewidth=1, linestyle=':')
    ax2.set_xlabel('Divergence D')
    ax2.set_ylabel('Cumulative Fraction')
    ax2.set_title('CDF of D under random reassignment')
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG, '13_classification_sensitivity.png'), dpi=150)
    print(f"\nPlot saved: {os.path.join(FIG, '13_classification_sensitivity.png')}")
except ImportError:
    print("matplotlib not available")

print(f"\nResults saved: {os.path.join(OUT, 'classification_sensitivity.json')}")
print("Done.")
