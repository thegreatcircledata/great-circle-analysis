#!/usr/bin/env python3
"""
Paper 3 Robustness Battery — Task 1: Exhaustive Pole Search
============================================================
Generate all pole positions at 1° resolution (360 × 180 = 64,800 poles).
For each, compute monument-settlement divergence D using Pleiades at 50 km.
Report: how many poles produce D ≥ Alison's D? What percentile is Alison?
Where are the top 10 poles geographically?
"""

import csv, json, math, os, sys, time
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
CORRIDOR = 50.0  # km

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

# ============================================================
# LOAD PLEIADES
# ============================================================
print("Loading Pleiades data...")
t0 = time.time()
mon_lats, mon_lons = [], []
set_lats, set_lons = [], []

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
        if is_m and not is_s:
            mon_lats.append(lat)
            mon_lons.append(lon)
        elif is_s and not is_m:
            set_lats.append(lat)
            set_lons.append(lon)

mon_lats = np.array(mon_lats)
mon_lons = np.array(mon_lons)
set_lats = np.array(set_lats)
set_lons = np.array(set_lons)
print(f"  Monuments: {len(mon_lats)}, Settlements: {len(set_lats)} ({time.time()-t0:.1f}s)")

# Precompute radians
mon_lats_r = np.radians(mon_lats)
mon_lons_r = np.radians(mon_lons)
set_lats_r = np.radians(set_lats)
set_lons_r = np.radians(set_lons)

# ============================================================
# DIVERGENCE COMPUTATION
# ============================================================
def compute_divergence(pole_lat, pole_lon):
    """
    For a given pole, compute D = (mon_fraction / set_fraction) - 1
    where fraction = sites within CORRIDOR km of the great circle.
    D > 0 means monuments are enriched relative to settlements.
    """
    la1 = np.radians(pole_lat)
    lo1 = np.radians(pole_lon)

    # Monument distances
    dl = mon_lats_r - la1
    dn = mon_lons_r - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(mon_lats_r)*np.sin(dn/2)**2
    d = R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    mon_gc_dist = np.abs(d - QC)
    mon_in = np.sum(mon_gc_dist <= CORRIDOR)

    # Settlement distances
    dl = set_lats_r - la1
    dn = set_lons_r - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(set_lats_r)*np.sin(dn/2)**2
    d = R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    set_gc_dist = np.abs(d - QC)
    set_in = np.sum(set_gc_dist <= CORRIDOR)

    mon_frac = mon_in / len(mon_lats) if len(mon_lats) > 0 else 0
    set_frac = set_in / len(set_lats) if len(set_lats) > 0 else 0

    if set_frac > 0:
        D = mon_frac / set_frac
    else:
        D = float('inf') if mon_frac > 0 else 1.0

    return D, mon_in, set_in, mon_frac, set_frac

# ============================================================
# ALISON'S CIRCLE FIRST
# ============================================================
print("\nComputing Alison's circle divergence...")
D_alison, mon_in_a, set_in_a, mon_f_a, set_f_a = compute_divergence(POLE_LAT, POLE_LON)
print(f"  Alison pole: ({POLE_LAT:.4f}, {POLE_LON:.4f})")
print(f"  Monuments in corridor: {mon_in_a} ({100*mon_f_a:.2f}%)")
print(f"  Settlements in corridor: {set_in_a} ({100*set_f_a:.2f}%)")
print(f"  D (ratio) = {D_alison:.4f}")

# ============================================================
# EXHAUSTIVE POLE SCAN (1° grid)
# ============================================================
print("\n" + "=" * 70)
print("EXHAUSTIVE POLE SCAN — 1° resolution")
print("=" * 70)

# Due to symmetry, pole at (lat, lon) gives same circle as (-lat, lon+180).
# So we only need to scan one hemisphere of poles.
# But for completeness and to produce the full heatmap, scan all.
# Actually, for efficiency, scan lat 0-90 and mirror.
# A pole at (lat, lon) produces the same great circle as (-lat, (lon+180)%360 - 180).

pole_lats = np.arange(-89.5, 90.5, 1.0)  # -89.5 to 89.5
pole_lons = np.arange(-179.5, 180.5, 1.0)  # -179.5 to 179.5

n_poles = len(pole_lats) * len(pole_lons)
print(f"Total poles: {n_poles}")

D_values = np.zeros((len(pole_lats), len(pole_lons)))
mon_in_arr = np.zeros((len(pole_lats), len(pole_lons)), dtype=int)
set_in_arr = np.zeros((len(pole_lats), len(pole_lons)), dtype=int)

t0 = time.time()
n_done = 0
for i, plat in enumerate(pole_lats):
    for j, plon in enumerate(pole_lons):
        D, mi, si, _, _ = compute_divergence(plat, plon)
        D_values[i, j] = D
        mon_in_arr[i, j] = mi
        set_in_arr[i, j] = si
        n_done += 1

    if (i + 1) % 10 == 0:
        elapsed = time.time() - t0
        rate = n_done / elapsed
        remaining = (n_poles - n_done) / rate
        print(f"  Row {i+1}/{len(pole_lats)} ({n_done}/{n_poles} poles, "
              f"{rate:.0f} poles/s, ~{remaining/60:.0f} min remaining)")

total_time = time.time() - t0
print(f"\nScan complete in {total_time/60:.1f} minutes ({total_time:.0f}s)")

# ============================================================
# ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("RESULTS")
print("=" * 70)

# Flatten for analysis
D_flat = D_values.ravel()

# Filter out inf and nan
valid_mask = np.isfinite(D_flat) & (D_flat > 0)
D_valid = D_flat[valid_mask]

n_ge_alison = np.sum(D_valid >= D_alison)
percentile = 100 * (1 - n_ge_alison / len(D_valid))

print(f"\nAlison's D = {D_alison:.4f}")
print(f"Valid poles (finite D > 0): {len(D_valid)} / {n_poles}")
print(f"Poles with D ≥ Alison's: {n_ge_alison} ({100*n_ge_alison/len(D_valid):.2f}%)")
print(f"Alison's percentile: {percentile:.2f}th")

# Top 10 poles
top_idx = np.argsort(D_valid)[::-1][:100]  # get more to find unique circles
D_sorted = D_valid[top_idx]

# Map back to lat/lon
valid_indices = np.where(valid_mask)[0]
top_poles = []
seen_circles = set()
for k in range(len(top_idx)):
    flat_idx = valid_indices[top_idx[k]]
    i = flat_idx // len(pole_lons)
    j = flat_idx % len(pole_lons)
    plat = pole_lats[i]
    plon = pole_lons[j]
    D_val = D_valid[top_idx[k]]

    # Check if this is a duplicate circle (antipodal pole)
    key1 = (round(plat, 0), round(plon, 0))
    key2 = (round(-plat, 0), round((plon + 180) % 360 - 180, 0))
    if key1 in seen_circles or key2 in seen_circles:
        continue
    seen_circles.add(key1)
    seen_circles.add(key2)

    top_poles.append({
        'rank': len(top_poles) + 1,
        'pole_lat': float(plat),
        'pole_lon': float(plon),
        'D': float(D_val),
        'mon_in_corridor': int(mon_in_arr[i, j]),
        'set_in_corridor': int(set_in_arr[i, j])
    })
    if len(top_poles) >= 10:
        break

print(f"\nTop 10 poles (unique circles):")
print(f"{'Rank':>4} {'Lat':>7} {'Lon':>8} {'D':>8} {'Mon':>5} {'Set':>5}")
print("-" * 42)
for p in top_poles:
    print(f"{p['rank']:>4d} {p['pole_lat']:>7.1f} {p['pole_lon']:>8.1f} {p['D']:>8.3f} {p['mon_in_corridor']:>5d} {p['set_in_corridor']:>5d}")

# Find Alison's rank
alison_rank = np.sum(D_valid > D_alison) + 1
print(f"\nAlison's rank: {alison_rank} out of {len(D_valid)} valid poles")
print(f"Alison's rank among unique circles: ~{alison_rank // 2}")

# ============================================================
# SAVE
# ============================================================
results = {
    'analysis': 'Exhaustive Pole Search',
    'grid_resolution_deg': 1.0,
    'total_poles': int(n_poles),
    'valid_poles': int(len(D_valid)),
    'corridor_km': CORRIDOR,
    'alison': {
        'pole_lat': POLE_LAT,
        'pole_lon': POLE_LON,
        'D': float(D_alison),
        'mon_in_corridor': int(mon_in_a),
        'set_in_corridor': int(set_in_a),
        'mon_fraction': float(mon_f_a),
        'set_fraction': float(set_f_a),
        'rank': int(alison_rank),
        'percentile': float(percentile),
        'n_poles_ge_D': int(n_ge_alison)
    },
    'top_10_poles': top_poles,
    'computation_time_sec': float(total_time),
    'D_distribution': {
        'mean': float(np.mean(D_valid)),
        'median': float(np.median(D_valid)),
        'std': float(np.std(D_valid)),
        'min': float(np.min(D_valid)),
        'max': float(np.max(D_valid)),
        'p5': float(np.percentile(D_valid, 5)),
        'p25': float(np.percentile(D_valid, 25)),
        'p75': float(np.percentile(D_valid, 75)),
        'p95': float(np.percentile(D_valid, 95)),
        'p99': float(np.percentile(D_valid, 99))
    }
}

with open(os.path.join(OUT, 'exhaustive_pole_scan.json'), 'w') as f:
    json.dump(results, f, indent=2)

# ============================================================
# HEATMAP PLOT
# ============================================================
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))

    # Top left: D heatmap (capped for visualization)
    ax = axes[0, 0]
    D_plot = np.clip(D_values, 0, np.percentile(D_valid, 99))
    im = ax.imshow(D_plot, origin='lower', aspect='auto',
                   extent=[-180, 180, -90, 90],
                   cmap='RdBu_r', vmin=0.5, vmax=2.0)
    ax.plot(POLE_LON, POLE_LAT, 'k*', markersize=15, label='Alison pole')
    for p in top_poles[:5]:
        ax.plot(p['pole_lon'], p['pole_lat'], 'gv', markersize=8)
    ax.set_xlabel('Pole Longitude')
    ax.set_ylabel('Pole Latitude')
    ax.set_title(f'Monument/Settlement Divergence D\n(50 km corridor, 1° grid)')
    plt.colorbar(im, ax=ax, label='D (mon_frac / set_frac)')
    ax.legend(loc='lower left', fontsize=8)

    # Top right: D histogram
    ax2 = axes[0, 1]
    ax2.hist(D_valid, bins=100, range=(0, 3), color='steelblue', alpha=0.7, edgecolor='none')
    ax2.axvline(D_alison, color='red', linewidth=2, linestyle='--', label=f'Alison D={D_alison:.3f}')
    ax2.set_xlabel('Divergence D')
    ax2.set_ylabel('Count')
    ax2.set_title(f'Distribution of D across {len(D_valid):,} poles\n'
                  f'Alison: {percentile:.1f}th percentile (rank {alison_rank})')
    ax2.legend()

    # Bottom left: Monuments in corridor heatmap
    ax3 = axes[1, 0]
    im3 = ax3.imshow(mon_in_arr, origin='lower', aspect='auto',
                     extent=[-180, 180, -90, 90], cmap='Reds')
    ax3.plot(POLE_LON, POLE_LAT, 'k*', markersize=12)
    ax3.set_xlabel('Pole Longitude')
    ax3.set_ylabel('Pole Latitude')
    ax3.set_title('Monuments in 50 km corridor')
    plt.colorbar(im3, ax=ax3, label='Count')

    # Bottom right: CDF
    ax4 = axes[1, 1]
    D_sort = np.sort(D_valid)
    cdf = np.arange(1, len(D_sort)+1) / len(D_sort)
    ax4.plot(D_sort, cdf, 'b-', linewidth=1)
    ax4.axvline(D_alison, color='red', linewidth=2, linestyle='--', label=f'Alison D={D_alison:.3f}')
    ax4.axhline(percentile/100, color='red', linewidth=0.5, linestyle=':')
    ax4.set_xlabel('Divergence D')
    ax4.set_ylabel('Cumulative Fraction')
    ax4.set_title('CDF of D values')
    ax4.legend()
    ax4.set_xlim(0, 3)
    ax4.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG, '12_exhaustive_pole_scan.png'), dpi=150)
    print(f"\nPlot saved: {os.path.join(FIG, '12_exhaustive_pole_scan.png')}")
except ImportError:
    print("matplotlib not available")

print(f"\nResults saved: {os.path.join(OUT, 'exhaustive_pole_scan.json')}")
print("Done.")
