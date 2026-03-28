#!/usr/bin/env python3
"""
Exhaustive Pole Scan — FILTERED version
========================================
Same as 07_exhaustive_pole_scan.py but only includes poles where
the corridor contains ≥30 monuments AND ≥30 settlements.
This eliminates small-number artifacts.
"""

import csv, json, os, sys, time
import numpy as np

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE = os.path.expanduser("~/megalith_site_research")
OUT = os.path.join(BASE, "paper3/results")
FIG = os.path.join(BASE, "paper3/figures")

POLE_LAT = 59.682122
POLE_LON = -138.646087
R = 6371.0
QC = R * np.pi / 2
CORRIDOR = 50.0
MIN_MON = 30
MIN_SET = 30

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

# Load Pleiades
print("Loading Pleiades...")
mon_lats, mon_lons, set_lats, set_lons = [], [], [], []
with open(os.path.join(BASE, "data/pleiades/pleiades-places-latest.csv")) as f:
    for row in csv.DictReader(f):
        try:
            lat, lon = float(row['reprLat']), float(row['reprLong'])
        except (ValueError, TypeError):
            continue
        if abs(lat) < 0.001 and abs(lon) < 0.001: continue
        fts = set(ft.strip() for ft in row['featureTypes'].split(','))
        if (fts & MONUMENTAL) and not (fts & SETTLEMENT):
            mon_lats.append(lat); mon_lons.append(lon)
        elif (fts & SETTLEMENT) and not (fts & MONUMENTAL):
            set_lats.append(lat); set_lons.append(lon)

mon_lats = np.array(mon_lats); mon_lons = np.array(mon_lons)
set_lats = np.array(set_lats); set_lons = np.array(set_lons)
mon_lats_r = np.radians(mon_lats); mon_lons_r = np.radians(mon_lons)
set_lats_r = np.radians(set_lats); set_lons_r = np.radians(set_lons)
print(f"  Monuments: {len(mon_lats)}, Settlements: {len(set_lats)}")

def compute_divergence(pole_lat, pole_lon):
    la1 = np.radians(pole_lat)
    lo1 = np.radians(pole_lon)
    dl = mon_lats_r - la1; dn = mon_lons_r - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(mon_lats_r)*np.sin(dn/2)**2
    d = R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    mon_in = int(np.sum(np.abs(d - QC) <= CORRIDOR))

    dl = set_lats_r - la1; dn = set_lons_r - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(set_lats_r)*np.sin(dn/2)**2
    d = R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    set_in = int(np.sum(np.abs(d - QC) <= CORRIDOR))

    mon_frac = mon_in / len(mon_lats)
    set_frac = set_in / len(set_lats)
    D = mon_frac / set_frac if set_frac > 0 else float('inf')
    return D, mon_in, set_in

# Alison's circle
D_alison, mon_a, set_a = compute_divergence(POLE_LAT, POLE_LON)
print(f"\nAlison: D={D_alison:.4f}, mon={mon_a}, set={set_a}")

# Scan
print("\nScanning 64,800 poles (1° grid)...")
pole_lats = np.arange(-89.5, 90.5, 1.0)
pole_lons = np.arange(-179.5, 180.5, 1.0)
t0 = time.time()

filtered_D = []
filtered_poles = []
all_D_filtered = []

for i, plat in enumerate(pole_lats):
    for j, plon in enumerate(pole_lons):
        D, mi, si = compute_divergence(plat, plon)
        if mi >= MIN_MON and si >= MIN_SET:
            filtered_D.append(D)
            filtered_poles.append((plat, plon, D, mi, si))

    if (i+1) % 30 == 0:
        print(f"  Row {i+1}/{len(pole_lats)}, filtered so far: {len(filtered_D)}")

total_time = time.time() - t0
print(f"\nDone in {total_time:.0f}s")

filtered_D = np.array(filtered_D)
n_filtered = len(filtered_D)

print(f"\n{'='*60}")
print(f"FILTERED RESULTS (≥{MIN_MON} monuments AND ≥{MIN_SET} settlements)")
print(f"{'='*60}")
print(f"Poles passing filter: {n_filtered} / 64,800 ({100*n_filtered/64800:.1f}%)")

n_ge = np.sum(filtered_D >= D_alison)
percentile = 100 * (1 - n_ge / n_filtered)
rank = int(np.sum(filtered_D > D_alison)) + 1

print(f"Alison D = {D_alison:.4f}")
print(f"Rank: {rank} / {n_filtered}")
print(f"Percentile: {percentile:.2f}th")
print(f"Poles with D ≥ Alison: {n_ge} ({100*n_ge/n_filtered:.2f}%)")

print(f"\nD distribution (filtered):")
print(f"  Mean: {np.mean(filtered_D):.4f}")
print(f"  Median: {np.median(filtered_D):.4f}")
print(f"  p95: {np.percentile(filtered_D, 95):.4f}")
print(f"  p99: {np.percentile(filtered_D, 99):.4f}")
print(f"  Max: {np.max(filtered_D):.4f}")

# Top 10
sorted_poles = sorted(filtered_poles, key=lambda x: x[2], reverse=True)
top10 = []
print(f"\nTop 10 filtered poles:")
for k, (plat, plon, D, mi, si) in enumerate(sorted_poles[:10]):
    print(f"  #{k+1}: ({plat:.1f}, {plon:.1f}) D={D:.4f} mon={mi} set={si}")
    top10.append({'rank': k+1, 'pole_lat': float(plat), 'pole_lon': float(plon),
                  'D': float(D), 'mon': mi, 'set': si})

# Is Alison in the filtered set?
alison_passes = mon_a >= MIN_MON and set_a >= MIN_SET
print(f"\nAlison passes filter: {alison_passes} (mon={mon_a}, set={set_a})")

results = {
    'analysis': 'Exhaustive Pole Scan (Filtered)',
    'filter': f'≥{MIN_MON} monuments AND ≥{MIN_SET} settlements in 50km corridor',
    'total_poles': 64800,
    'filtered_poles': n_filtered,
    'filter_rate': float(n_filtered / 64800),
    'alison': {
        'D': float(D_alison), 'mon': mon_a, 'set': set_a,
        'passes_filter': alison_passes,
        'rank': rank, 'percentile': float(percentile),
        'n_ge_D': int(n_ge)
    },
    'D_distribution': {
        'mean': float(np.mean(filtered_D)),
        'median': float(np.median(filtered_D)),
        'std': float(np.std(filtered_D)),
        'p95': float(np.percentile(filtered_D, 95)),
        'p99': float(np.percentile(filtered_D, 99)),
        'max': float(np.max(filtered_D))
    },
    'top_10': top10
}

with open(os.path.join(OUT, 'exhaustive_pole_scan_filtered.json'), 'w') as f:
    json.dump(results, f, indent=2)

# Plot
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    ax.hist(filtered_D, bins=50, color='steelblue', alpha=0.7, edgecolor='none')
    ax.axvline(D_alison, color='red', linewidth=2, linestyle='--',
               label=f'Alison D={D_alison:.3f} ({percentile:.1f}th %ile)')
    ax.set_xlabel('Divergence D')
    ax.set_ylabel('Count')
    ax.set_title(f'Filtered Pole Scan (≥{MIN_MON} mon, ≥{MIN_SET} set)\n'
                 f'{n_filtered:,} poles pass filter; Alison rank {rank}/{n_filtered}')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, '12b_exhaustive_pole_scan_filtered.png'), dpi=150)
    print(f"\nPlot saved")
except ImportError:
    pass

print(f"\nResults saved: {os.path.join(OUT, 'exhaustive_pole_scan_filtered.json')}")
print("Done.")
