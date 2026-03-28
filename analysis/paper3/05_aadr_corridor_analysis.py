#!/usr/bin/env python3
"""
Paper 3 — Task 4: AADR Corridor Analysis
==========================================
Test whether ancient DNA samples show enrichment near the Great Circle.
Standard MC enrichment test + temporal breakdown.
"""

import csv, json, os, sys, time
import numpy as np

np.random.seed(42)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE = os.path.expanduser("~/megalith_site_research")
OUT = os.path.join(BASE, "paper3/results")
FIG = os.path.join(BASE, "paper3/figures")
POLE_LAT, POLE_LON = 59.682122, -138.646087
R, QC = 6371.0, 6371.0 * np.pi / 2
THRESHOLD = 50
N_MC = 10000

def gc_dist_vec(lats, lons):
    la1, lo1 = np.radians(POLE_LAT), np.radians(POLE_LON)
    la2, lo2 = np.radians(lats), np.radians(lons)
    dl, dn = la2 - la1, lo2 - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(la2)*np.sin(dn/2)**2
    d = R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QC)

# ============================================================
# LOAD AADR
# ============================================================
print("Loading AADR v62.0...")
anno_file = os.path.join(BASE, "paper3/data/aadr/v62.0_1240k_public.anno")

samples = []
with open(anno_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    for row in reader:
        if len(row) < 18:
            continue
        try:
            lat = float(row[16].strip())
            lon = float(row[17].strip())
        except (ValueError, TypeError):
            continue
        if abs(lat) < 0.001 and abs(lon) < 0.001:
            continue

        date_str = row[9].strip()
        try:
            date_bp = float(date_str)
        except (ValueError, TypeError):
            date_bp = None

        gen_id = row[0].strip()
        is_modern = '.DG' in gen_id or date_str == '0' or (date_bp is not None and date_bp <= 0)

        if is_modern:
            continue
        if date_bp is None or date_bp <= 0:
            continue

        samples.append({
            'lat': lat, 'lon': lon,
            'date_bp': date_bp,
            'country': row[15].strip() if len(row) > 15 else '',
            'locality': row[14].strip() if len(row) > 14 else ''
        })

lats = np.array([s['lat'] for s in samples])
lons = np.array([s['lon'] for s in samples])
dates = np.array([s['date_bp'] for s in samples])
n = len(samples)
print(f"Ancient samples with coords and dates: {n}")

# ============================================================
# ENRICHMENT TEST
# ============================================================
print(f"\nRunning {THRESHOLD} km corridor enrichment test ({N_MC} MC trials)...")

dists = gc_dist_vec(lats, lons)
observed = int(np.sum(dists <= THRESHOLD))
print(f"Observed within {THRESHOLD} km: {observed}")

# MC baseline: random great circles
t0 = time.time()
mc_counts = np.zeros(N_MC)
rng = np.random.RandomState(42)
for i in range(N_MC):
    plat = rng.uniform(-90, 90)
    plon = rng.uniform(-180, 180)
    la1, lo1 = np.radians(plat), np.radians(plon)
    la2, lo2 = np.radians(lats), np.radians(lons)
    dl, dn = la2 - la1, lo2 - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(la2)*np.sin(dn/2)**2
    d = R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    rd = np.abs(d - QC)
    mc_counts[i] = np.sum(rd <= THRESHOLD)
    if (i+1) % 2000 == 0:
        print(f"  MC {i+1}/{N_MC}")

mc_mean = mc_counts.mean()
mc_std = mc_counts.std()
z = (observed - mc_mean) / (mc_std + 1e-10)
p_value = np.sum(mc_counts >= observed) / N_MC
enrichment = observed / (mc_mean + 1e-10)
percentile = 100 * np.sum(mc_counts < observed) / N_MC

print(f"\nResults:")
print(f"  Observed: {observed}")
print(f"  MC mean: {mc_mean:.1f} ± {mc_std:.1f}")
print(f"  Z-score: {z:.2f}")
print(f"  Enrichment: {enrichment:.2f}×")
print(f"  p-value: {p_value:.4f}")
print(f"  Percentile: {percentile:.1f}")
print(f"  Time: {time.time()-t0:.1f}s")

# ============================================================
# TEMPORAL BREAKDOWN
# ============================================================
print("\n" + "=" * 70)
print("TEMPORAL BREAKDOWN")
print("=" * 70)

bins_bp = [(0, 1000), (1000, 2000), (2000, 3000), (3000, 4000),
           (4000, 5000), (5000, 7500), (7500, 10000), (10000, 15000),
           (15000, 25000), (25000, 50000), (50000, 120000)]

temporal_results = []
print(f"{'Period':<20} {'n':>6} {'On-circ':>8} {'Expected':>9} {'Z':>6} {'Enrich':>8}")
print("-" * 60)

for lo, hi in bins_bp:
    mask = (dates >= lo) & (dates < hi)
    if np.sum(mask) < 10:
        temporal_results.append({
            'period': f'{lo}-{hi} BP', 'n': int(np.sum(mask)),
            'too_few': True
        })
        print(f"{lo}-{hi} BP{'':<8} {int(np.sum(mask)):>6} — too few")
        continue

    t_lats = lats[mask]
    t_lons = lons[mask]
    t_dists = gc_dist_vec(t_lats, t_lons)
    t_obs = int(np.sum(t_dists <= THRESHOLD))

    # Quick MC (200 trials for temporal bins)
    t_mc = []
    for _ in range(200):
        plat = rng.uniform(-90, 90)
        plon = rng.uniform(-180, 180)
        la1, lo1 = np.radians(plat), np.radians(plon)
        la2, lo2 = np.radians(t_lats), np.radians(t_lons)
        dl, dn = la2 - la1, lo2 - lo1
        a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(la2)*np.sin(dn/2)**2
        d = R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
        rd = np.abs(d - QC)
        t_mc.append(int(np.sum(rd <= THRESHOLD)))

    t_mc = np.array(t_mc)
    t_z = (t_obs - t_mc.mean()) / (t_mc.std() + 1e-10)
    t_enrich = t_obs / (t_mc.mean() + 1e-10)

    temporal_results.append({
        'period': f'{lo}-{hi} BP',
        'n': int(np.sum(mask)),
        'observed': t_obs,
        'expected': float(t_mc.mean()),
        'z': float(t_z),
        'enrichment': float(t_enrich),
        'too_few': False
    })

    flag = " ***" if t_z > 2.0 else " **" if t_z > 1.5 else ""
    print(f"{lo}-{hi} BP{'':<8} {int(np.sum(mask)):>6} {t_obs:>8} {t_mc.mean():>9.1f} {t_z:>6.2f} {t_enrich:>7.2f}×{flag}")

# ============================================================
# ON-CIRCLE SAMPLES
# ============================================================
print("\n" + "=" * 70)
print(f"ON-CIRCLE SAMPLES (within {THRESHOLD} km)")
print("=" * 70)

oncircle_idx = np.where(dists <= THRESHOLD)[0]
oncircle_samples = []
for idx in sorted(oncircle_idx, key=lambda i: dists[i]):
    s = samples[idx]
    d = dists[idx]
    oncircle_samples.append({
        'dist_km': round(float(d), 1),
        'date_bp': s['date_bp'],
        'country': s['country'],
        'locality': s['locality'][:50]
    })
    if len(oncircle_samples) <= 30:
        print(f"  {d:5.1f} km | {s['date_bp']:>8.0f} BP | {s['country']:<15} | {s['locality'][:40]}")

if len(oncircle_samples) > 30:
    print(f"  ... ({len(oncircle_samples)} total)")

from collections import Counter
oc_countries = Counter(s['country'] for s in oncircle_samples)
print(f"\nOn-circle by country:")
for c, cnt in oc_countries.most_common(10):
    print(f"  {c}: {cnt}")

# ============================================================
# SAVE
# ============================================================
results = {
    'n_ancient_samples': n,
    'threshold_km': THRESHOLD,
    'n_mc': N_MC,
    'observed': observed,
    'mc_mean': float(mc_mean),
    'mc_std': float(mc_std),
    'z_score': float(z),
    'p_value': float(p_value),
    'enrichment': float(enrichment),
    'percentile': float(percentile),
    'temporal': temporal_results,
    'on_circle_countries': dict(oc_countries.most_common(20)),
    'n_on_circle': len(oncircle_samples)
}

with open(os.path.join(OUT, 'aadr_corridor_analysis.json'), 'w') as f:
    json.dump(results, f, indent=2)

# ============================================================
# PLOT
# ============================================================
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Left: MC histogram
    ax = axes[0]
    ax.hist(mc_counts, bins=50, color='#3498db', alpha=0.7, edgecolor='black')
    ax.axvline(x=observed, color='red', linewidth=2, linestyle='--',
               label=f'Observed = {observed}')
    ax.axvline(x=mc_mean, color='green', linewidth=1, linestyle=':',
               label=f'Expected = {mc_mean:.0f}')
    ax.set_xlabel(f'Samples within {THRESHOLD} km')
    ax.set_ylabel('MC Count')
    ax.set_title(f'AADR Enrichment Test (Z={z:.2f})')
    ax.legend()
    ax.grid(alpha=0.3)

    # Middle: Temporal breakdown
    ax2 = axes[1]
    valid = [t for t in temporal_results if not t.get('too_few', True)]
    if valid:
        periods = [t['period'] for t in valid]
        zs = [t['z'] for t in valid]
        enrichments = [t['enrichment'] for t in valid]
        x = range(len(periods))
        colors = ['#c0392b' if z > 2 else '#e67e22' if z > 1 else '#3498db' for z in zs]
        ax2.bar(x, zs, color=colors, alpha=0.8, edgecolor='black')
        ax2.set_xticks(list(x))
        ax2.set_xticklabels(periods, rotation=45, ha='right', fontsize=7)
        ax2.axhline(y=0, color='black', linewidth=0.5)
        ax2.axhline(y=1.96, color='red', linestyle='--', alpha=0.5, label='Z=1.96')
        ax2.set_ylabel('Z-score')
        ax2.set_title('AADR Enrichment by Time Period')
        ax2.legend()
        ax2.grid(axis='y', alpha=0.3)

    # Right: Distance histogram
    ax3 = axes[2]
    ax3.hist(dists, bins=np.arange(0, 2001, 50), color='#27ae60', alpha=0.7, edgecolor='black')
    ax3.axvline(x=THRESHOLD, color='red', linewidth=2, linestyle='--', label=f'{THRESHOLD} km')
    ax3.set_xlabel('Distance from Great Circle (km)')
    ax3.set_ylabel('Number of Samples')
    ax3.set_title('AADR Sample Distance Distribution')
    ax3.legend()
    ax3.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG, '10_aadr_corridor.png'), dpi=150)
    print(f"\nPlot saved: {FIG}/10_aadr_corridor.png")
except ImportError:
    print("matplotlib not available")

print(f"Results saved: {OUT}/aadr_corridor_analysis.json")
print("Done.")
