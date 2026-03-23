#!/usr/bin/env python3
"""
Alternative Site Triplets: Monument-Settlement Divergence Comparison
=====================================================================
For each triplet of famous archaeological sites, find the best-fit great
circle, compute monument-settlement divergence D using Pleiades, rank
against 100 random circles, and compare to the Alison family (D > 7).

Methodology:
- For each triplet: find the pole that minimizes the max distance from
  the three sites to the great circle (best-fit circle).
- Compute D = Z_monument - Z_settlement (50 km threshold, 200 MC trials)
- Generate 100 distribution-matched random circles and compute D for each
  to establish a percentile rank for the triplet's circle.
"""

import csv, math, os, sys, time, json
import numpy as np
from scipy.optimize import minimize

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
np.random.seed(42)

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * np.pi / 2
THRESHOLD_KM = 50
N_MC = 200
N_RANDOM = 100  # random circles for ranking

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_FILE = os.path.join(BASE_DIR, "results", "triplet_divergence_comparison.json")

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

# ============================================================
# SITE COORDINATES
# ============================================================
SITES = {
    "Giza":           (29.9792,   31.1342),
    "Nazca":          (-14.735,  -75.130),
    "Easter_Island":  (-27.116, -109.350),
    "Angkor_Wat":     (13.4125,  103.867),
    "Machu_Picchu":   (-13.163,  -72.546),
    "Gobekli_Tepe":   (37.2233,   38.9224),
    "Mohenjo_daro":   (27.3242,   68.1375),
    "Stonehenge":     (51.1789,   -1.8262),
    "Persepolis":     (29.9353,   52.8917),
    "Teotihuacan":    (19.6925,  -98.8438),
    "Athens":         (37.9715,   23.7267),  # Parthenon
    "Varanasi":       (25.3109,   83.0107),
    "Cahokia":        (38.6605,  -90.0622),
    "Carnac":         (47.5847,   -3.0780),
    "Luxor":          (25.7188,   32.6573),  # Karnak Temple
    "Tiwanaku":       (-16.5546, -68.6731),
    "Petra":          (30.3285,   35.4444),
    "Borobudur":      (-7.6079,  110.2038),
}

TRIPLETS = [
    ("Giza + Angkor Wat + Machu Picchu",       ["Giza", "Angkor_Wat", "Machu_Picchu"]),
    ("Giza + Gobekli Tepe + Mohenjo-daro",     ["Giza", "Gobekli_Tepe", "Mohenjo_daro"]),
    ("Stonehenge + Giza + Persepolis",          ["Stonehenge", "Giza", "Persepolis"]),
    ("Teotihuacan + Giza + Angkor Wat",        ["Teotihuacan", "Giza", "Angkor_Wat"]),
    ("Machu Picchu + Stonehenge + Angkor Wat",  ["Machu_Picchu", "Stonehenge", "Angkor_Wat"]),
    ("Athens + Persepolis + Varanasi",          ["Athens", "Persepolis", "Varanasi"]),
    ("Cahokia + Carnac + Luxor",               ["Cahokia", "Carnac", "Luxor"]),
    ("Tiwanaku + Petra + Borobudur",           ["Tiwanaku", "Petra", "Borobudur"]),
]

# Alison circle for reference
ALISON_POLE = (59.682122, -138.646087)

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def haversine_km(lat1, lon1, lat2, lon2):
    lat1r, lon1r = np.radians(lat1), np.radians(lon1)
    lat2r, lon2r = np.radians(lat2), np.radians(lon2)
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat/2)**2 + np.cos(lat1r)*np.cos(lat2r)*np.sin(dlon/2)**2
    return EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))

def dist_from_gc_scalar(site_lat, site_lon, pole_lat, pole_lon):
    d = haversine_km(site_lat, site_lon, pole_lat, pole_lon)
    return abs(d - QUARTER_CIRC)

def dist_from_gc_vec(site_lats, site_lons, pole_lat, pole_lon):
    d = haversine_km(pole_lat, pole_lon, site_lats, site_lons)
    return np.abs(d - QUARTER_CIRC)

def count_within_vec(site_lats, site_lons, pole_lat, pole_lon):
    dists = dist_from_gc_vec(site_lats, site_lons, pole_lat, pole_lon)
    return int(np.sum(dists <= THRESHOLD_KM))

def mc_zscore_vec(site_lats, site_lons, pole_lat, pole_lon, n_trials):
    n = len(site_lats)
    observed = count_within_vec(site_lats, site_lons, pole_lat, pole_lon)
    rand_counts = np.empty(n_trials)
    for t in range(n_trials):
        idx = np.random.randint(0, n, size=n)
        rlats = np.clip(site_lats[idx] + np.random.normal(0, 2, n), -90, 90)
        rlons = np.clip(site_lons[idx] + np.random.normal(0, 2, n), -180, 180)
        rand_counts[t] = count_within_vec(rlats, rlons, pole_lat, pole_lon)
    mu = np.mean(rand_counts)
    sigma = np.std(rand_counts)
    return float((observed - mu) / sigma) if sigma > 0 else 0.0

def latlon_to_xyz(lat, lon):
    """Convert lat/lon (degrees) to unit vector on sphere."""
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    return np.array([
        np.cos(lat_r) * np.cos(lon_r),
        np.cos(lat_r) * np.sin(lon_r),
        np.sin(lat_r)
    ])

def xyz_to_latlon(x, y, z):
    """Convert unit vector to lat/lon (degrees)."""
    lat = np.degrees(np.arcsin(np.clip(z, -1, 1)))
    lon = np.degrees(np.arctan2(y, x))
    return lat, lon

def best_fit_pole(site_keys):
    """Find the great circle pole that minimizes max distance to all three sites.

    Strategy: Start from the cross-product of the first two sites' position vectors
    (the exact pole for a 2-site great circle), then optimize to minimize the max
    miss across all three sites.
    """
    coords = [SITES[k] for k in site_keys]

    # Initial guess: pole of great circle through first two sites
    v0 = latlon_to_xyz(coords[0][0], coords[0][1])
    v1 = latlon_to_xyz(coords[1][0], coords[1][1])
    pole_init = np.cross(v0, v1)
    norm = np.linalg.norm(pole_init)
    if norm < 1e-10:
        # Sites are antipodal or identical, use arbitrary perpendicular
        pole_init = np.array([1, 0, 0])
    else:
        pole_init = pole_init / norm

    init_lat, init_lon = xyz_to_latlon(*pole_init)

    def objective(params):
        plat, plon = params
        dists = [dist_from_gc_scalar(c[0], c[1], plat, plon) for c in coords]
        return max(dists)  # minimax

    # Try both the pole and its antipode as starting points
    antipode_lat, antipode_lon = -init_lat, (init_lon + 180) % 360 - 180

    results = []
    for start_lat, start_lon in [(init_lat, init_lon), (antipode_lat, antipode_lon)]:
        res = minimize(objective, [start_lat, start_lon],
                      method='Nelder-Mead',
                      options={'xatol': 0.001, 'fatol': 0.01, 'maxiter': 5000})
        results.append(res)

    best = min(results, key=lambda r: r.fun)
    pole_lat, pole_lon = best.x
    # Normalize longitude
    pole_lon = ((pole_lon + 180) % 360) - 180

    # Compute per-site distances
    site_dists = {k: round(dist_from_gc_scalar(SITES[k][0], SITES[k][1], pole_lat, pole_lon), 1)
                  for k in site_keys}
    max_miss = round(best.fun, 1)

    return pole_lat, pole_lon, site_dists, max_miss

def random_pole_matched(all_lats, all_lons):
    """Generate a distribution-matched random pole."""
    lat = np.random.choice(all_lats) + np.random.normal(0, 10)
    lon = np.random.choice(all_lons) + np.random.normal(0, 10)
    return max(-90.0, min(90.0, lat)), ((lon + 180) % 360) - 180

# ============================================================
# LOAD PLEIADES DATA
# ============================================================
print("=" * 70)
print("LOADING PLEIADES DATA")
print("=" * 70)

pleiades_path = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")
if not os.path.exists(pleiades_path):
    pleiades_path = os.path.join(BASE_DIR, "pleiades-places-latest.csv")

mon_list = []
set_list = []

with open(pleiades_path, encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            lat = float(row["reprLat"])
            lon = float(row["reprLong"])
            if lat == 0 and lon == 0:
                continue
            feature_types = {t.strip() for t in row.get("featureTypes", "").split(",")}
            if feature_types & MONUMENTAL_TYPES:
                mon_list.append((lat, lon))
            elif feature_types & SETTLEMENT_TYPES:
                set_list.append((lat, lon))
        except (ValueError, KeyError, TypeError):
            pass

mon_arr = np.array(mon_list)
set_arr = np.array(set_list)
mon_lats, mon_lons = mon_arr[:, 0], mon_arr[:, 1]
set_lats, set_lons = set_arr[:, 0], set_arr[:, 1]
all_lats = np.concatenate([mon_lats, set_lats])
all_lons = np.concatenate([mon_lons, set_lons])

print(f"Pleiades monuments: {len(mon_arr)}")
print(f"Pleiades settlements: {len(set_arr)}")

# ============================================================
# COMPUTE ALISON CIRCLE REFERENCE D
# ============================================================
print("\n" + "=" * 70)
print("ALISON CIRCLE REFERENCE")
print("=" * 70)

alison_mon_z = mc_zscore_vec(mon_lats, mon_lons, ALISON_POLE[0], ALISON_POLE[1], N_MC)
alison_set_z = mc_zscore_vec(set_lats, set_lons, ALISON_POLE[0], ALISON_POLE[1], N_MC)
alison_D = alison_mon_z - alison_set_z
print(f"Alison pole: ({ALISON_POLE[0]:.4f}, {ALISON_POLE[1]:.4f})")
print(f"  Z_monument = {alison_mon_z:.2f}, Z_settlement = {alison_set_z:.2f}")
print(f"  D = {alison_D:.2f}")

# ============================================================
# PROCESS EACH TRIPLET
# ============================================================
print("\n" + "=" * 70)
print(f"PROCESSING {len(TRIPLETS)} TRIPLETS")
print("=" * 70)

triplet_results = []
t_total_start = time.time()

for ti, (triplet_name, site_keys) in enumerate(TRIPLETS):
    print(f"\n{'─' * 70}")
    print(f"[{ti+1}/{len(TRIPLETS)}] {triplet_name}")
    print(f"{'─' * 70}")

    t_start = time.time()

    # Step 1: Find best-fit great circle
    pole_lat, pole_lon, site_dists, max_miss = best_fit_pole(site_keys)
    print(f"  Best-fit pole: ({pole_lat:.4f}, {pole_lon:.4f})")
    for sk, sd in site_dists.items():
        print(f"    {sk}: {sd:.1f} km from circle")
    print(f"  Max miss: {max_miss:.1f} km")

    # Step 2: Compute D for the triplet's circle
    mon_z = mc_zscore_vec(mon_lats, mon_lons, pole_lat, pole_lon, N_MC)
    set_z = mc_zscore_vec(set_lats, set_lons, pole_lat, pole_lon, N_MC)
    D = mon_z - set_z
    mon_count = count_within_vec(mon_lats, mon_lons, pole_lat, pole_lon)
    set_count = count_within_vec(set_lats, set_lons, pole_lat, pole_lon)
    print(f"  Z_monument = {mon_z:.2f} ({mon_count} within 50km)")
    print(f"  Z_settlement = {set_z:.2f} ({set_count} within 50km)")
    print(f"  D = {D:.2f}")

    # Step 3: Rank against 100 random circles
    print(f"  Computing {N_RANDOM} random circle baselines...")
    random_Ds = []
    for ri in range(N_RANDOM):
        rp_lat, rp_lon = random_pole_matched(all_lats, all_lons)
        r_mon_z = mc_zscore_vec(mon_lats, mon_lons, rp_lat, rp_lon, 50)  # 50 MC for speed
        r_set_z = mc_zscore_vec(set_lats, set_lons, rp_lat, rp_lon, 50)
        random_Ds.append(r_mon_z - r_set_z)
        if (ri + 1) % 25 == 0:
            elapsed = time.time() - t_start
            print(f"    {ri+1}/{N_RANDOM} random circles done ({elapsed:.0f}s)")

    random_Ds = np.array(random_Ds)
    percentile = float(100 * np.sum(random_Ds < D) / len(random_Ds))
    random_mean = float(np.mean(random_Ds))
    random_std = float(np.std(random_Ds))
    sigma_above = (D - random_mean) / random_std if random_std > 0 else 0

    exceeds_alison_family = D > 7.0

    elapsed = time.time() - t_start
    print(f"  Percentile rank: {percentile:.0f}th among random circles")
    print(f"  Random baseline: mean={random_mean:.2f}, std={random_std:.2f}")
    print(f"  Sigma above random mean: {sigma_above:.2f}")
    print(f"  Comparable to Alison family (D > 7)? {'YES' if exceeds_alison_family else 'NO'}")
    print(f"  Time: {elapsed:.0f}s")

    triplet_results.append({
        "name": triplet_name,
        "sites": {k: {"lat": SITES[k][0], "lon": SITES[k][1]} for k in site_keys},
        "best_fit_pole": {"lat": round(pole_lat, 4), "lon": round(pole_lon, 4)},
        "site_distances_km": site_dists,
        "max_miss_km": max_miss,
        "monument_z": round(mon_z, 3),
        "settlement_z": round(set_z, 3),
        "monuments_within_50km": mon_count,
        "settlements_within_50km": set_count,
        "divergence_D": round(D, 3),
        "random_baseline": {
            "n_random": N_RANDOM,
            "mc_trials_per_random": 50,
            "mean_D": round(random_mean, 3),
            "std_D": round(random_std, 3),
            "percentile": round(percentile, 1),
            "sigma_above_mean": round(sigma_above, 2),
        },
        "exceeds_alison_family_threshold": exceeds_alison_family,
        "elapsed_seconds": round(elapsed, 1),
    })

total_elapsed = time.time() - t_total_start

# ============================================================
# SUMMARY TABLE
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY TABLE")
print("=" * 70)

print(f"\n{'Triplet':<45} {'MaxMiss':>7} {'D':>7} {'%ile':>5} {'σ':>6} {'D>7?':>5}")
print("─" * 80)
for tr in triplet_results:
    print(f"{tr['name']:<45} {tr['max_miss_km']:>6.0f}km {tr['divergence_D']:>7.2f} "
          f"{tr['random_baseline']['percentile']:>4.0f}% {tr['random_baseline']['sigma_above_mean']:>5.1f}σ "
          f"{'YES' if tr['exceeds_alison_family_threshold'] else ' no':>5}")

print(f"\n{'Alison circle (reference)':<45} {'6.5km':>7} {alison_D:>7.2f} {'—':>5} {'—':>6} {'YES':>5}")

# Count
n_exceeds = sum(1 for tr in triplet_results if tr["exceeds_alison_family_threshold"])
n_positive = sum(1 for tr in triplet_results if tr["divergence_D"] > 0)
print(f"\n{n_exceeds}/{len(triplet_results)} triplets exceed D > 7 (Alison family threshold)")
print(f"{n_positive}/{len(triplet_results)} triplets have positive D (monument-favoring)")
print(f"Total time: {total_elapsed:.0f}s")

# ============================================================
# SAVE RESULTS
# ============================================================
results = {
    "meta": {
        "description": "Monument-settlement divergence D for alternative archaeological site triplets, "
                       "compared to 100 random circles and the Alison family threshold (D > 7)",
        "date": "2026-03-20",
        "divergence_threshold_km": THRESHOLD_KM,
        "mc_trials_triplet": N_MC,
        "mc_trials_random": 50,
        "n_random_circles": N_RANDOM,
        "pleiades_monuments": len(mon_arr),
        "pleiades_settlements": len(set_arr),
        "alison_family_threshold": 7.0,
    },
    "alison_reference": {
        "pole": {"lat": ALISON_POLE[0], "lon": ALISON_POLE[1]},
        "monument_z": round(alison_mon_z, 3),
        "settlement_z": round(alison_set_z, 3),
        "divergence_D": round(alison_D, 3),
    },
    "triplets": triplet_results,
    "summary": {
        "n_triplets": len(triplet_results),
        "n_exceeding_alison_threshold": n_exceeds,
        "n_positive_D": n_positive,
        "D_values": [tr["divergence_D"] for tr in triplet_results],
        "max_D": max(tr["divergence_D"] for tr in triplet_results),
        "min_D": min(tr["divergence_D"] for tr in triplet_results),
        "mean_D": round(np.mean([tr["divergence_D"] for tr in triplet_results]), 3),
    },
    "elapsed_seconds": round(total_elapsed, 1),
}

os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
with open(OUT_FILE, "w") as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to {OUT_FILE}")
