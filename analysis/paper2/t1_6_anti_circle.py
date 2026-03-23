#!/usr/bin/env python3
"""
Directive T1-6: Anti-Circle Negative Control
=============================================
Compute monument-settlement divergence D for circles with poles
perpendicular or opposite to the Alison pole. These should show D ≈ 0.
"""

import json, math, os, sys, csv
import numpy as np

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "anti_circle")
os.makedirs(OUT_DIR, exist_ok=True)

POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
THRESHOLD_KM = 50
N_TRIALS = 200

MONUMENTAL_TYPES = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus", "shrine"
}
SETTLEMENT_TYPES = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement", "townhouse", "production"
}

def haversine_km(lat1, lon1, lat2, lon2):
    R = EARTH_R_KM
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(min(1.0, a)))

def load_pleiades():
    path = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")
    monuments, settlements = [], []
    with open(path, encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row['reprLat'])
                lon = float(row['reprLong'])
            except (ValueError, KeyError, TypeError):
                continue
            if lat == 0 and lon == 0:
                continue
            ftypes = {t.strip() for t in row.get('featureTypes', '').split(',')}
            if ftypes & MONUMENTAL_TYPES:
                monuments.append((lat, lon))
            elif ftypes & SETTLEMENT_TYPES:
                settlements.append((lat, lon))
    return np.array(monuments), np.array(settlements)

def compute_divergence(pole_lat, pole_lon, mon_arr, set_arr, threshold=THRESHOLD_KM, n_trials=N_TRIALS):
    """Compute monument-settlement divergence D for a given pole."""
    # Distance from great circle for each site
    def dist_from_gc_arr(sites):
        dists = np.empty(len(sites))
        for i, (lat, lon) in enumerate(sites):
            d = haversine_km(lat, lon, pole_lat, pole_lon)
            dists[i] = abs(d - QUARTER_CIRC)
        return dists

    mon_dists = dist_from_gc_arr(mon_arr)
    set_dists = dist_from_gc_arr(set_arr)

    # Observed counts within threshold
    obs_mon = np.sum(mon_dists < threshold)
    obs_set = np.sum(set_dists < threshold)

    # Monte Carlo baseline (distribution-matched)
    mc_mon_counts = []
    mc_set_counts = []
    all_sites = np.vstack([mon_arr, set_arr])

    for _ in range(n_trials):
        # Jitter all sites
        jitter_lat = np.random.normal(0, 2, len(all_sites))
        jitter_lon = np.random.normal(0, 2, len(all_sites))
        jittered = all_sites.copy()
        jittered[:, 0] += jitter_lat
        jittered[:, 1] += jitter_lon

        # Resample monuments and settlements
        idx_mon = np.random.choice(len(jittered), len(mon_arr), replace=True)
        idx_set = np.random.choice(len(jittered), len(set_arr), replace=True)

        jit_mon_dists = dist_from_gc_arr(jittered[idx_mon])
        jit_set_dists = dist_from_gc_arr(jittered[idx_set])

        mc_mon_counts.append(np.sum(jit_mon_dists < threshold))
        mc_set_counts.append(np.sum(jit_set_dists < threshold))

    mc_mon = np.array(mc_mon_counts, dtype=float)
    mc_set = np.array(mc_set_counts, dtype=float)

    # Z-scores
    z_mon = (obs_mon - mc_mon.mean()) / max(mc_mon.std(), 1) if mc_mon.std() > 0 else 0
    z_set = (obs_set - mc_set.mean()) / max(mc_set.std(), 1) if mc_set.std() > 0 else 0

    D = z_mon - z_set

    return {
        'pole_lat': pole_lat,
        'pole_lon': pole_lon,
        'obs_monuments': int(obs_mon),
        'obs_settlements': int(obs_set),
        'mc_mon_mean': float(mc_mon.mean()),
        'mc_set_mean': float(mc_set.mean()),
        'z_monument': float(z_mon),
        'z_settlement': float(z_set),
        'D': float(D)
    }

# ============================================================
# DEFINE TEST CIRCLES
# ============================================================
test_circles = {
    'Alison (reference)': (POLE_LAT, POLE_LON),
    'Perpendicular 1 (lon+90)': (POLE_LAT, POLE_LON + 90),
    'Perpendicular 2 (lat rotated)': (90 - POLE_LAT, POLE_LON + 180),
    'Opposite hemisphere': (0.0, POLE_LON - 90),
    'Antipodal pole': (-POLE_LAT, POLE_LON + 180),
}

# Normalize longitudes
for name in test_circles:
    lat, lon = test_circles[name]
    lon = ((lon + 180) % 360) - 180
    lat = max(-90, min(90, lat))
    test_circles[name] = (lat, lon)

# ============================================================
# RUN ANALYSIS
# ============================================================
print("Loading Pleiades data...")
monuments, settlements = load_pleiades()
print(f"Loaded {len(monuments)} monuments, {len(settlements)} settlements")

results = {}
for name, (plat, plon) in test_circles.items():
    print(f"\nComputing D for: {name} (pole: {plat:.2f}°, {plon:.2f}°)...")
    result = compute_divergence(plat, plon, monuments, settlements)
    results[name] = result
    print(f"  obs_mon={result['obs_monuments']}, obs_set={result['obs_settlements']}")
    print(f"  Z_mon={result['z_monument']:.2f}, Z_set={result['z_settlement']:.2f}")
    print(f"  D = {result['D']:.2f}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("ANTI-CIRCLE NEGATIVE CONTROL RESULTS")
print(f"{'='*70}")
print(f"{'Circle':<35} {'Pole':>20} {'D':>8} {'Z_mon':>8} {'Z_set':>8}")
print("-" * 80)
for name, result in results.items():
    print(f"{name:<35} ({result['pole_lat']:.1f}, {result['pole_lon']:.1f})  "
          f"{result['D']:>7.2f} {result['z_monument']:>7.2f} {result['z_settlement']:>7.2f}")

alison_D = results['Alison (reference)']['D']
print(f"\nAlison D = {alison_D:.2f}")
other_Ds = [v['D'] for k, v in results.items() if k != 'Alison (reference)']
print(f"Perpendicular/opposite D values: {[f'{d:.2f}' for d in other_Ds]}")
print(f"Max |D| of controls: {max(abs(d) for d in other_Ds):.2f}")
if abs(alison_D) > 2 * max(abs(d) for d in other_Ds):
    print("RESULT: Alison D is significantly larger than all controls ✓")
else:
    print("RESULT: Controls show non-trivial D — investigate further")

# ============================================================
# SAVE
# ============================================================
with open(os.path.join(OUT_DIR, "perpendicular_results.json"), 'w') as f:
    json.dump(results, f, indent=2)

summary = []
for name, result in results.items():
    summary.append(f"## {name}\n- Pole: ({result['pole_lat']:.2f}°, {result['pole_lon']:.2f}°)\n"
                   f"- D = {result['D']:.2f} (Z_mon={result['z_monument']:.2f}, Z_set={result['z_settlement']:.2f})\n"
                   f"- Observed: {result['obs_monuments']} monuments, {result['obs_settlements']} settlements within {THRESHOLD_KM}km\n")

summary_md = f"""# Anti-Circle Negative Control Results

## Purpose
Test whether the Alison Great Circle's monument-settlement divergence D is specific to its geometry,
or whether any great circle would produce similar results.

## Method
Computed D using Pleiades monument/settlement classification with {N_TRIALS}-trial Monte Carlo baseline
at {THRESHOLD_KM}km threshold. Tested perpendicular circles and an opposite-hemisphere circle.

## Results
{"".join(summary)}

## Conclusion
Alison D = {alison_D:.2f}
Maximum control |D| = {max(abs(d) for d in other_Ds):.2f}
{"The Alison circle shows significantly higher divergence than all control circles, confirming geometric specificity." if abs(alison_D) > 2 * max(abs(d) for d in other_Ds) else "Some control circles show non-trivial D — further investigation needed."}
"""

with open(os.path.join(OUT_DIR, "summary.md"), 'w') as f:
    f.write(summary_md)

print(f"\nDone! Outputs saved to {OUT_DIR}")
