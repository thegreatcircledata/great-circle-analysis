#!/usr/bin/env python3
"""
Fertile Crescent Control for Deep-Time Enrichment
===================================================
For each of the 100 random circles (same seed=42 as younger_dryas_circle_test.py):
  1. Compute early Holocene enrichment (10,500–8,500 BP in 14C years)
  2. Determine if the circle passes through the Fertile Crescent (28–38°N, 35–50°E)
  3. Compare Alison's 4–6× enrichment to Fertile Crescent–passing circles

Output: results/fertile_crescent_control.json
"""

import csv, math, json, os, time
import numpy as np

EARTH_R_KM = 6371
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
THRESHOLD_KM = 200
N_RANDOM_CIRCLES = 100

# Fertile Crescent bounding box
FC_LAT_MIN, FC_LAT_MAX = 28.0, 38.0
FC_LON_MIN, FC_LON_MAX = 35.0, 50.0

# Early Holocene window: 10,500–8,500 BP (14C years)
EH_START = 8500   # low end
EH_END = 10500    # high end
# Bins that fall within this: [8500,9000), [9000,9500), [9500,10000), [10000,10500)

BASE_DIR = os.path.expanduser("~/megalith_site_research")
BIN_EDGES_BP = list(range(5000, 20500, 500))
BINS_BP = [(BIN_EDGES_BP[i], BIN_EDGES_BP[i+1]) for i in range(len(BIN_EDGES_BP)-1)]

# Early Holocene bin indices
EH_BINS = [(lo, hi) for lo, hi in BINS_BP if lo >= EH_START and hi <= EH_END]
print(f"Early Holocene bins: {EH_BINS}")


def haversine_vec(lat1, lon1, lats, lons):
    lat1_r, lon1_r = np.radians(lat1), np.radians(lon1)
    lat2_r = np.radians(lats)
    lon2_r = np.radians(lons)
    dlat = lat2_r - lat1_r
    dlon = lon2_r - lon1_r
    a = np.sin(dlat/2)**2 + np.cos(lat1_r)*np.cos(lat2_r)*np.sin(dlon/2)**2
    return EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))


def gc_dist_from_pole(pole_lat, pole_lon, lats, lons):
    d = haversine_vec(pole_lat, pole_lon, lats, lons)
    return np.abs(d - QUARTER_CIRC)


def random_pole():
    z = np.random.uniform(-1, 1)
    theta = np.random.uniform(0, 2 * np.pi)
    lat = np.degrees(np.arcsin(z))
    lon = np.degrees(theta) - 180
    return lat, lon


def circle_passes_through_box(pole_lat, pole_lon, lat_min, lat_max, lon_min, lon_max, n_points=3600):
    """Check if a great circle (defined by pole) passes through a lat/lon box.
    Sample n_points around the circle and check if any fall in the box."""
    # The great circle is the set of points exactly QUARTER_CIRC from the pole.
    # Generate points on the circle by rotating around the pole.
    pole_lat_r = np.radians(pole_lat)
    pole_lon_r = np.radians(pole_lon)

    # Create an orthonormal basis. The pole defines the normal to the great circle plane.
    # Normal vector (pole direction):
    nx = np.cos(pole_lat_r) * np.cos(pole_lon_r)
    ny = np.cos(pole_lat_r) * np.sin(pole_lon_r)
    nz = np.sin(pole_lat_r)

    # Find two vectors perpendicular to the normal
    if abs(nz) < 0.9:
        u = np.array([-ny, nx, 0.0])
    else:
        u = np.array([1.0, 0.0, -nx/nz]) if nz != 0 else np.array([1.0, 0.0, 0.0])
    u = u / np.linalg.norm(u)
    v = np.cross(np.array([nx, ny, nz]), u)
    v = v / np.linalg.norm(v)

    # Points on the great circle: cos(t)*u + sin(t)*v
    t = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    xs = np.cos(t) * u[0] + np.sin(t) * v[0]
    ys = np.cos(t) * u[1] + np.sin(t) * v[1]
    zs = np.cos(t) * u[2] + np.sin(t) * v[2]

    lats = np.degrees(np.arcsin(np.clip(zs, -1, 1)))
    lons = np.degrees(np.arctan2(ys, xs))

    in_box = (lats >= lat_min) & (lats <= lat_max) & (lons >= lon_min) & (lons <= lon_max)
    return bool(np.any(in_box))


# Load p3k14c data
print("Loading p3k14c data...")
p3k_file = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")

all_lats = []
all_lons = []
all_ages_bp = []
skipped = 0

with open(p3k_file, encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            lat = float(row['Lat'])
            lon = float(row['Long'])
            age_bp = float(row['Age']) if row['Age'] else None
        except (ValueError, TypeError):
            skipped += 1
            continue
        if age_bp is None or age_bp <= 0:
            skipped += 1
            continue
        if lat == 0 and lon == 0:
            skipped += 1
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            skipped += 1
            continue
        all_lats.append(lat)
        all_lons.append(lon)
        all_ages_bp.append(age_bp)

all_lats = np.array(all_lats)
all_lons = np.array(all_lons)
all_ages_bp = np.array(all_ages_bp)
print(f"  Loaded: {len(all_lats)} dates")

# Precompute global counts per bin
global_counts = {}
for lo, hi in BINS_BP:
    global_counts[(lo, hi)] = int(np.sum((all_ages_bp >= lo) & (all_ages_bp < hi)))


def compute_enrichment(pole_lat, pole_lon):
    """Compute per-bin and early Holocene enrichment for a given great circle."""
    dists = gc_dist_from_pole(pole_lat, pole_lon, all_lats, all_lons)
    near = dists <= THRESHOLD_KM
    n_near_total = int(np.sum(near))
    near_frac = n_near_total / len(all_lats)

    if near_frac == 0:
        return None

    bin_enrichments = {}
    for lo, hi in BINS_BP:
        age_mask = (all_ages_bp >= lo) & (all_ages_bp < hi)
        g = global_counts[(lo, hi)]
        n = int(np.sum(age_mask & near))
        ratio = n / g if g > 0 else 0
        e = ratio / near_frac if near_frac > 0 and g > 0 else 0
        bin_enrichments[(lo, hi)] = {
            'enrichment': round(e, 4),
            'global': g,
            'near': n,
        }

    # Early Holocene enrichment (mean of EH bins with >=50 global dates)
    eh_enrichments = []
    eh_details = []
    for lo, hi in EH_BINS:
        b = bin_enrichments[(lo, hi)]
        if b['global'] >= 50:
            eh_enrichments.append(b['enrichment'])
            eh_details.append({
                'bin': f"{hi}-{lo} BP",
                'enrichment': b['enrichment'],
                'global': b['global'],
                'near': b['near'],
            })

    eh_mean = float(np.mean(eh_enrichments)) if eh_enrichments else 0.0

    return {
        'n_near': n_near_total,
        'near_fraction': round(near_frac, 6),
        'eh_mean_enrichment': round(eh_mean, 4),
        'eh_bin_details': eh_details,
        'n_eh_bins_used': len(eh_enrichments),
    }


# Compute Alison circle
print("\nComputing Alison circle enrichment...")
alison = compute_enrichment(POLE_LAT, POLE_LON)
alison_passes_fc = circle_passes_through_box(POLE_LAT, POLE_LON, FC_LAT_MIN, FC_LAT_MAX, FC_LON_MIN, FC_LON_MAX)
print(f"  Alison EH enrichment: {alison['eh_mean_enrichment']:.3f}x")
print(f"  Alison passes through Fertile Crescent: {alison_passes_fc}")
print(f"  Per-bin: {alison['eh_bin_details']}")

# Compute 100 random circles (same seed as original)
print("\nComputing 100 random circles...")
np.random.seed(42)

random_results = []
t0 = time.time()

for i in range(N_RANDOM_CIRCLES):
    plat, plon = random_pole()
    result = compute_enrichment(plat, plon)
    passes_fc = circle_passes_through_box(plat, plon, FC_LAT_MIN, FC_LAT_MAX, FC_LON_MIN, FC_LON_MAX)

    entry = {
        'index': i,
        'pole_lat': round(plat, 4),
        'pole_lon': round(plon, 4),
        'passes_fertile_crescent': passes_fc,
    }

    if result:
        entry.update(result)
    else:
        entry['eh_mean_enrichment'] = 0.0
        entry['n_near'] = 0
        entry['eh_bin_details'] = []

    random_results.append(entry)

    if (i + 1) % 25 == 0:
        print(f"  {i+1}/100 ({time.time()-t0:.1f}s)")

# Analysis
print("\n" + "=" * 70)
print("RESULTS")
print("=" * 70)

all_eh = [r['eh_mean_enrichment'] for r in random_results if r.get('n_eh_bins_used', 0) > 0]
fc_circles = [r for r in random_results if r['passes_fertile_crescent'] and r.get('n_eh_bins_used', 0) > 0]
non_fc_circles = [r for r in random_results if not r['passes_fertile_crescent'] and r.get('n_eh_bins_used', 0) > 0]

fc_eh = [r['eh_mean_enrichment'] for r in fc_circles]
non_fc_eh = [r['eh_mean_enrichment'] for r in non_fc_circles]

# How many show 4-6× enrichment?
n_4_6x = sum(1 for e in all_eh if 4.0 <= e <= 6.0)
n_gt_4x = sum(1 for e in all_eh if e >= 4.0)

print(f"\nAll 100 random circles:")
print(f"  Valid circles (with EH data): {len(all_eh)}")
print(f"  EH enrichment: mean={np.mean(all_eh):.3f}, median={np.median(all_eh):.3f}, std={np.std(all_eh):.3f}")
print(f"  Range: [{min(all_eh):.3f}, {max(all_eh):.3f}]")
print(f"  Circles with 4-6× EH enrichment: {n_4_6x}")
print(f"  Circles with ≥4× EH enrichment: {n_gt_4x}")

print(f"\nFertile Crescent-passing circles: {len(fc_circles)}")
if fc_eh:
    print(f"  EH enrichment: mean={np.mean(fc_eh):.3f}, median={np.median(fc_eh):.3f}, std={np.std(fc_eh):.3f}")
    print(f"  Range: [{min(fc_eh):.3f}, {max(fc_eh):.3f}]")
    fc_4_6x = sum(1 for e in fc_eh if 4.0 <= e <= 6.0)
    fc_gt_4x = sum(1 for e in fc_eh if e >= 4.0)
    print(f"  FC circles with 4-6× EH enrichment: {fc_4_6x}")
    print(f"  FC circles with ≥4× EH enrichment: {fc_gt_4x}")

    # Where does Alison rank among FC circles?
    alison_eh = alison['eh_mean_enrichment']
    fc_rank = sum(1 for e in fc_eh if e >= alison_eh) + 1
    print(f"\n  Alison EH enrichment: {alison_eh:.3f}x")
    print(f"  Alison rank among FC circles: {fc_rank}/{len(fc_eh)} (percentile: {100*(1-fc_rank/len(fc_eh)):.1f}%)")

    # List all FC circles sorted by EH enrichment
    print(f"\n  All FC-passing circles (sorted by EH enrichment):")
    fc_sorted = sorted(fc_circles, key=lambda x: x['eh_mean_enrichment'], reverse=True)
    for r in fc_sorted:
        print(f"    pole=({r['pole_lat']:.1f}, {r['pole_lon']:.1f})  EH={r['eh_mean_enrichment']:.3f}x  n_near={r['n_near']}")

print(f"\nNon-FC circles: {len(non_fc_circles)}")
if non_fc_eh:
    print(f"  EH enrichment: mean={np.mean(non_fc_eh):.3f}, median={np.median(non_fc_eh):.3f}")

# Is Alison's enrichment unusual among FC circles?
alison_eh = alison['eh_mean_enrichment']
if fc_eh:
    fc_mean = np.mean(fc_eh)
    fc_std = np.std(fc_eh)
    fc_z = (alison_eh - fc_mean) / fc_std if fc_std > 0 else 0
    print(f"\n  Alison vs FC circles:")
    print(f"    Alison EH: {alison_eh:.3f}x")
    print(f"    FC mean: {fc_mean:.3f}x, FC std: {fc_std:.3f}")
    print(f"    Z-score: {fc_z:.3f}")
    print(f"    Interpretation: {'UNUSUAL (Z>2)' if abs(fc_z) > 2 else 'TYPICAL' if abs(fc_z) < 1 else 'SOMEWHAT UNUSUAL (1<Z<2)'}")

# Save results
output = {
    "meta": {
        "analysis": "Fertile Crescent Control for Deep-Time Enrichment",
        "date": "2026-03-20",
        "description": (
            "Tests whether the Alison circle's 4-6x early Holocene enrichment "
            "(10,500-8,500 BP) is unusual among random circles, particularly those "
            "passing through the Fertile Crescent (28-38°N, 35-50°E)."
        ),
        "methodology": {
            "database": "p3k14c",
            "n_dates": int(len(all_lats)),
            "threshold_km": THRESHOLD_KM,
            "early_holocene_window": "8,500-10,500 14C BP",
            "early_holocene_bins": [f"{hi}-{lo} BP" for lo, hi in EH_BINS],
            "fertile_crescent_box": {
                "lat_range": [FC_LAT_MIN, FC_LAT_MAX],
                "lon_range": [FC_LON_MIN, FC_LON_MAX],
            },
            "n_random_circles": N_RANDOM_CIRCLES,
            "random_seed": 42,
            "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
        },
    },
    "alison_circle": {
        "eh_mean_enrichment": alison_eh,
        "passes_fertile_crescent": alison_passes_fc,
        "eh_bin_details": alison['eh_bin_details'],
        "n_near_circle": alison['n_near'],
    },
    "all_100_random": {
        "n_valid": len(all_eh),
        "eh_mean": round(float(np.mean(all_eh)), 4),
        "eh_median": round(float(np.median(all_eh)), 4),
        "eh_std": round(float(np.std(all_eh)), 4),
        "eh_min": round(float(min(all_eh)), 4),
        "eh_max": round(float(max(all_eh)), 4),
        "n_with_4_6x": n_4_6x,
        "n_with_gte_4x": n_gt_4x,
        "percentile_of_alison": round(100 * sum(1 for e in all_eh if e < alison_eh) / len(all_eh), 2),
    },
    "fertile_crescent_circles": {
        "n_passing": len(fc_circles),
        "n_total_random": N_RANDOM_CIRCLES,
        "pct_passing": round(100 * len(fc_circles) / N_RANDOM_CIRCLES, 1),
        "eh_enrichments": sorted([round(e, 4) for e in fc_eh], reverse=True) if fc_eh else [],
        "eh_mean": round(float(np.mean(fc_eh)), 4) if fc_eh else None,
        "eh_median": round(float(np.median(fc_eh)), 4) if fc_eh else None,
        "eh_std": round(float(np.std(fc_eh)), 4) if fc_eh else None,
        "n_with_4_6x": sum(1 for e in fc_eh if 4.0 <= e <= 6.0) if fc_eh else 0,
        "n_with_gte_4x": sum(1 for e in fc_eh if e >= 4.0) if fc_eh else 0,
        "alison_rank_among_fc": fc_rank if fc_eh else None,
        "alison_z_vs_fc": round(float(fc_z), 4) if fc_eh and fc_std > 0 else None,
        "circles": [
            {
                "pole_lat": r['pole_lat'],
                "pole_lon": r['pole_lon'],
                "eh_mean_enrichment": r['eh_mean_enrichment'],
                "n_near": r['n_near'],
                "eh_bin_details": r.get('eh_bin_details', []),
            }
            for r in sorted(fc_circles, key=lambda x: x['eh_mean_enrichment'], reverse=True)
        ],
    },
    "non_fc_circles": {
        "n_circles": len(non_fc_circles),
        "eh_mean": round(float(np.mean(non_fc_eh)), 4) if non_fc_eh else None,
        "eh_median": round(float(np.median(non_fc_eh)), 4) if non_fc_eh else None,
    },
    "interpretation": {
        "alison_eh_enrichment": alison_eh,
        "is_alison_unusual_globally": alison_eh > np.mean(all_eh) + 2*np.std(all_eh) if all_eh else None,
        "is_alison_unusual_among_fc": (abs(fc_z) > 2) if fc_eh and fc_std > 0 else None,
        "verdict": (
            f"Alison's {alison_eh:.2f}x early Holocene enrichment is "
            + ("TYPICAL" if fc_eh and abs(fc_z) < 1 else
               "SOMEWHAT UNUSUAL" if fc_eh and abs(fc_z) < 2 else
               "UNUSUAL")
            + f" among the {len(fc_circles)} Fertile Crescent-passing random circles "
            + f"(FC mean: {np.mean(fc_eh):.2f}x, Z={fc_z:.2f}). "
            + f"{n_4_6x}/100 random circles show 4-6x EH enrichment overall."
        ) if fc_eh else "Insufficient FC circles for comparison.",
    },
}

out_path = os.path.join(BASE_DIR, "results", "fertile_crescent_control.json")
with open(out_path, "w") as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nSaved: {out_path}")
print("\nDONE.")
