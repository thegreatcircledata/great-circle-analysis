#!/usr/bin/env python3
"""
Paper 3 Robustness Battery — Task 6: Independent Origins Test
==============================================================
Compile independently originated civilizations.
For all great circles passing through exactly 4+ of 6-7 origins,
compute mean divergence D. Where does Alison's circle rank?
"""

import json, os, sys, time
import numpy as np
from itertools import combinations

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

# ============================================================
# INDEPENDENT ORIGIN CIVILIZATIONS
# ============================================================
# Standard list from anthropology (Trigger 2003, Marcus & Feinman 1998)
# Using representative monumental center for each

civilizations = {
    'Mesopotamia': {'lat': 31.33, 'lon': 45.64, 'center': 'Ur/Uruk', 'date': '3500 BCE'},
    'Egypt': {'lat': 29.87, 'lon': 31.22, 'center': 'Memphis/Saqqara', 'date': '3100 BCE'},
    'Indus': {'lat': 27.32, 'lon': 68.14, 'center': 'Mohenjo-daro', 'date': '2600 BCE'},
    'China': {'lat': 34.27, 'lon': 108.94, 'center': 'Erlitou/Xian', 'date': '1900 BCE'},
    'Mesoamerica': {'lat': 19.69, 'lon': -98.84, 'center': 'Teotihuacan', 'date': '200 BCE'},
    'Andes': {'lat': -13.16, 'lon': -72.55, 'center': 'Cusco/Caral', 'date': '3000 BCE'},
    'West_Africa': {'lat': 11.50, 'lon': -4.00, 'center': 'Jenne-jeno', 'date': '250 BCE'},
}

civ_names = list(civilizations.keys())
civ_lats = np.array([civilizations[c]['lat'] for c in civ_names])
civ_lons = np.array([civilizations[c]['lon'] for c in civ_names])

print("Independent Origin Civilizations:")
for name, info in civilizations.items():
    print(f"  {name}: {info['center']} ({info['lat']:.2f}, {info['lon']:.2f}), {info['date']}")

# ============================================================
# DISTANCE FROM ALISON'S CIRCLE
# ============================================================
print("\n" + "=" * 70)
print("DISTANCES FROM ALISON'S GREAT CIRCLE")
print("=" * 70)

def gc_distance_to_circle(lats, lons, pole_lat, pole_lon):
    la1 = np.radians(pole_lat)
    lo1 = np.radians(pole_lon)
    la2 = np.radians(lats)
    lo2 = np.radians(lons)
    dl = la2 - la1
    dn = lo2 - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(la2)*np.sin(dn/2)**2
    d = R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QC)

alison_dists = gc_distance_to_circle(civ_lats, civ_lons, POLE_LAT, POLE_LON)

for i, name in enumerate(civ_names):
    marker = "✓" if alison_dists[i] <= 200 else ""
    print(f"  {name:<15} {alison_dists[i]:>8.1f} km  {marker}")

n_within_200 = np.sum(alison_dists <= 200)
n_within_500 = np.sum(alison_dists <= 500)
print(f"\nAlison's circle passes within 200 km of {n_within_200}/7 origins")
print(f"Alison's circle passes within 500 km of {n_within_500}/7 origins")

# ============================================================
# FIND GREAT CIRCLES THROUGH K OF N ORIGINS
# ============================================================
print("\n" + "=" * 70)
print("SEARCHING FOR CIRCLES THROUGH 4+ ORIGINS (200 km corridor)")
print("=" * 70)

THRESHOLD = 200  # km
MIN_ORIGINS = 4

# Strategy: for each pair of civilizations, compute the great circle
# pole and check how many other civilizations it passes near.
# A great circle through 2 points has a unique pole (the cross product).

def latlon_to_xyz(lat, lon):
    """Convert lat/lon to unit vector on sphere."""
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    x = np.cos(lat_r) * np.cos(lon_r)
    y = np.cos(lat_r) * np.sin(lon_r)
    z = np.sin(lat_r)
    return np.array([x, y, z])

def xyz_to_latlon(xyz):
    """Convert unit vector to lat/lon."""
    lat = np.degrees(np.arcsin(xyz[2]))
    lon = np.degrees(np.arctan2(xyz[1], xyz[0]))
    return lat, lon

circles_through_origins = []

for i, j in combinations(range(len(civ_names)), 2):
    # Compute pole of great circle through civ i and civ j
    v1 = latlon_to_xyz(civ_lats[i], civ_lons[i])
    v2 = latlon_to_xyz(civ_lats[j], civ_lons[j])

    pole = np.cross(v1, v2)
    norm = np.linalg.norm(pole)
    if norm < 1e-10:
        continue
    pole = pole / norm

    pole_lat, pole_lon = xyz_to_latlon(pole)

    # Distance from this circle to all civilizations
    dists = gc_distance_to_circle(civ_lats, civ_lons, pole_lat, pole_lon)
    within = np.sum(dists <= THRESHOLD)

    if within >= MIN_ORIGINS:
        which = [civ_names[k] for k in range(len(civ_names)) if dists[k] <= THRESHOLD]
        circles_through_origins.append({
            'defined_by': f"{civ_names[i]} — {civ_names[j]}",
            'pole_lat': float(pole_lat),
            'pole_lon': float(pole_lon),
            'n_origins_within': int(within),
            'origins': which,
            'distances': {civ_names[k]: float(dists[k]) for k in range(len(civ_names))},
            'mean_dist': float(np.mean(dists[dists <= THRESHOLD]))
        })

circles_through_origins.sort(key=lambda c: c['n_origins_within'], reverse=True)

print(f"\nCircles passing within {THRESHOLD} km of 4+ origins: {len(circles_through_origins)}")
for c in circles_through_origins[:10]:
    print(f"  {c['n_origins_within']} origins: {c['defined_by']} → {', '.join(c['origins'])}")
    print(f"    Pole: ({c['pole_lat']:.1f}, {c['pole_lon']:.1f})")

# ============================================================
# COMPUTE DIVERGENCE D FOR EACH ORIGIN-PASSING CIRCLE
# ============================================================
print("\n" + "=" * 70)
print("COMPUTING DIVERGENCE D FOR ORIGIN-PASSING CIRCLES")
print("=" * 70)

# Load Pleiades for D computation
import csv

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

mon_lats_p, mon_lons_p, set_lats_p, set_lons_p = [], [], [], []
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
        if (fts & MONUMENTAL) and not (fts & SETTLEMENT):
            mon_lats_p.append(lat)
            mon_lons_p.append(lon)
        elif (fts & SETTLEMENT) and not (fts & MONUMENTAL):
            set_lats_p.append(lat)
            set_lons_p.append(lon)

mon_lats_p = np.array(mon_lats_p)
mon_lons_p = np.array(mon_lons_p)
set_lats_p = np.array(set_lats_p)
set_lons_p = np.array(set_lons_p)
mon_lats_pr = np.radians(mon_lats_p)
mon_lons_pr = np.radians(mon_lons_p)
set_lats_pr = np.radians(set_lats_p)
set_lons_pr = np.radians(set_lons_p)

CORRIDOR_D = 50.0

def compute_D(pole_lat, pole_lon):
    la1 = np.radians(pole_lat)
    lo1 = np.radians(pole_lon)

    dl = mon_lats_pr - la1
    dn = mon_lons_pr - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(mon_lats_pr)*np.sin(dn/2)**2
    d = R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    mon_in = np.sum(np.abs(d - QC) <= CORRIDOR_D)

    dl = set_lats_pr - la1
    dn = set_lons_pr - lo1
    a = np.sin(dl/2)**2 + np.cos(la1)*np.cos(set_lats_pr)*np.sin(dn/2)**2
    d = R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    set_in = np.sum(np.abs(d - QC) <= CORRIDOR_D)

    mon_frac = mon_in / len(mon_lats_p) if len(mon_lats_p) > 0 else 0
    set_frac = set_in / len(set_lats_p) if len(set_lats_p) > 0 else 0
    D = mon_frac / set_frac if set_frac > 0 else float('inf')
    return D, int(mon_in), int(set_in)

# Alison's D
D_alison, mon_in_a, set_in_a = compute_D(POLE_LAT, POLE_LON)
print(f"\nAlison's D = {D_alison:.4f} (mon={mon_in_a}, set={set_in_a})")

# D for each origin-passing circle
for c in circles_through_origins:
    D, mi, si = compute_D(c['pole_lat'], c['pole_lon'])
    c['D'] = float(D)
    c['mon_in_50km'] = mi
    c['set_in_50km'] = si

# Sort by D
circles_through_origins.sort(key=lambda c: c.get('D', 0), reverse=True)

print(f"\n{'Origins':>7} {'D':>8} {'Mon':>5} {'Set':>5}  Defined by → Through")
print("-" * 80)
for c in circles_through_origins[:15]:
    print(f"{c['n_origins_within']:>7d} {c['D']:>8.3f} {c['mon_in_50km']:>5d} "
          f"{c['set_in_50km']:>5d}  {c['defined_by'][:25]}")

# Where does Alison rank?
all_D = [c['D'] for c in circles_through_origins if np.isfinite(c['D'])]
alison_rank = sum(1 for d in all_D if d > D_alison) + 1
print(f"\nAlison's D rank among origin-passing circles: {alison_rank}/{len(all_D)}")
print(f"Alison's circle passes through {n_within_200}/7 origins within 200km")

# ============================================================
# ALSO TEST: RANDOM CIRCLES FOR BASELINE D
# ============================================================
print("\n" + "=" * 70)
print("RANDOM CIRCLE BASELINE")
print("=" * 70)

N_RAND = 5000
rng = np.random.RandomState(42)
rand_D = np.zeros(N_RAND)

for trial in range(N_RAND):
    plat = np.degrees(np.arcsin(rng.uniform(-1, 1)))
    plon = rng.uniform(-180, 180)
    rand_D[trial], _, _ = compute_D(plat, plon)

print(f"Random circle D: mean={np.mean(rand_D):.4f}, median={np.median(rand_D):.4f}")
print(f"Alison D={D_alison:.4f}, random percentile: "
      f"{100*(1-np.mean(rand_D >= D_alison)):.1f}th")

# ============================================================
# SAVE
# ============================================================
results = {
    'analysis': 'Independent Origins Test',
    'civilizations': civilizations,
    'alison_distances_km': {civ_names[i]: float(alison_dists[i]) for i in range(len(civ_names))},
    'alison_within_200km': int(n_within_200),
    'alison_within_500km': int(n_within_500),
    'origin_threshold_km': THRESHOLD,
    'corridor_km': CORRIDOR_D,
    'circles_through_4plus_origins': len(circles_through_origins),
    'alison': {
        'D': float(D_alison),
        'rank_among_origin_circles': alison_rank,
        'total_origin_circles': len(all_D)
    },
    'random_baseline': {
        'n_trials': N_RAND,
        'mean_D': float(np.mean(rand_D)),
        'median_D': float(np.median(rand_D)),
        'alison_percentile': float(100*(1-np.mean(rand_D >= D_alison)))
    },
    'top_circles': circles_through_origins[:20]
}

with open(os.path.join(OUT, 'independent_origins_test.json'), 'w') as f:
    json.dump(results, f, indent=2)

# ============================================================
# PLOT
# ============================================================
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: D for origin-passing circles vs Alison
    ax = axes[0]
    origin_Ds = [c['D'] for c in circles_through_origins if np.isfinite(c['D'])]
    ax.hist(origin_Ds, bins=30, color='steelblue', alpha=0.7, edgecolor='none',
            label=f'Origin-passing circles (n={len(origin_Ds)})')
    ax.axvline(D_alison, color='red', linewidth=2, linestyle='--',
               label=f'Alison D={D_alison:.3f}')
    ax.axvline(1.0, color='black', linewidth=1, linestyle=':')
    ax.set_xlabel('Divergence D')
    ax.set_ylabel('Count')
    ax.set_title(f'D for circles through 4+ independent origins\n'
                 f'Alison rank: {alison_rank}/{len(all_D)}')
    ax.legend(fontsize=8)

    # Right: Civilization distances from Alison's circle
    ax2 = axes[1]
    y_pos = np.arange(len(civ_names))
    colors = ['green' if d <= 200 else 'orange' if d <= 500 else 'red'
              for d in alison_dists]
    ax2.barh(y_pos, alison_dists, color=colors, alpha=0.7)
    ax2.axvline(200, color='green', linewidth=1, linestyle='--', alpha=0.5, label='200 km')
    ax2.axvline(500, color='orange', linewidth=1, linestyle='--', alpha=0.5, label='500 km')
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([f"{n}\n({civilizations[n]['center']})" for n in civ_names], fontsize=8)
    ax2.set_xlabel('Distance from Alison\'s Great Circle (km)')
    ax2.set_title('Distance to Independent Origin Civilizations')
    ax2.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG, '14_independent_origins.png'), dpi=150)
    print(f"\nPlot saved: {os.path.join(FIG, '14_independent_origins.png')}")
except ImportError:
    print("matplotlib not available")

print(f"\nResults saved: {os.path.join(OUT, 'independent_origins_test.json')}")
print("Done.")
