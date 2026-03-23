#!/usr/bin/env python3
"""
Directive G: Systematic Great Circle Search
=============================================
Scans all possible great circles (2° grid → 0.5° refinement) to find
which pole positions produce the highest monument-settlement divergence D.
Answers: is the Alison Great Circle unique, or are there others?

Phase 1: Coarse grid (2° resolution, Poisson Z-scores) — ~8,236 poles
Phase 2: Fine refinement (0.5°, rigorous MC Z-scores) — around hot spots
Phase 3: Analysis, heatmaps, rankings
"""

import csv, json, math, os, sys, time
import numpy as np
from multiprocessing import Pool, cpu_count

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
np.random.seed(42)

# ============================================================
# CONFIGURATION
# ============================================================
ALISON_POLE_LAT = 59.682122
ALISON_POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * np.pi / 2
THRESHOLD_KM = 50
N_MC_FINE = 100        # MC trials for Phase 2 refinement
COARSE_STEP = 2        # degrees
FINE_STEP = 0.5        # degrees
FINE_BOX = 5           # degrees around each hot spot
HOT_SPOT_SEP = 10      # degrees minimum separation between peaks

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "systematic_gc_search")
os.makedirs(OUT_DIR, exist_ok=True)

# Pleiades type classifications (matching existing codebase)
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
# VECTORIZED HELPER FUNCTIONS
# ============================================================
def dist_from_gc_vec(site_lats, site_lons, pole_lat, pole_lon):
    """Vectorized distance from each site to the great circle defined by pole."""
    lat1r = np.radians(pole_lat)
    lon1r = np.radians(pole_lon)
    lat2r = np.radians(site_lats)
    lon2r = np.radians(site_lons)
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat/2)**2 + np.cos(lat1r)*np.cos(lat2r)*np.sin(dlon/2)**2
    d = EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QUARTER_CIRC)

def count_within_vec(site_lats, site_lons, pole_lat, pole_lon):
    dists = dist_from_gc_vec(site_lats, site_lons, pole_lat, pole_lon)
    return int(np.sum(dists <= THRESHOLD_KM))

def mc_zscore_vec(site_lats, site_lons, pole_lat, pole_lon, n_trials):
    """Per-circle MC Z-score with distribution-matched null."""
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
    return float((observed - mu) / sigma) if sigma > 0 else 0.0, observed, float(mu), float(sigma)

def haversine_deg(lat1, lon1, lat2, lon2):
    """Angular distance in degrees between two points."""
    lat1r, lon1r = np.radians(lat1), np.radians(lon1)
    lat2r, lon2r = np.radians(lat2), np.radians(lon2)
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat/2)**2 + np.cos(lat1r)*np.cos(lat2r)*np.sin(dlon/2)**2
    return np.degrees(2 * np.arcsin(np.sqrt(min(1.0, a))))

# ============================================================
# POISSON Z-SCORE (Option A — fast, approximate)
# ============================================================
# The great circle strip at ±threshold_km covers a fraction of Earth's surface.
# For a 50km threshold, the strip width is 100km total.
# Strip area = 2π·R · 2·threshold = circumference × strip_width
# Earth area = 4π·R²
# Fraction = 2·threshold / (2·R) = threshold / R ... wait, let me be more precise.
#
# The "strip" around a great circle at distance ≤ d_km consists of all points
# within d_km of the circle. On a sphere of radius R, the fraction of the
# sphere's surface within angular distance α of a great circle is:
#   f = sin(α) where α = d_km / R (in radians)
# For d=50km, α = 50/6371 = 0.00785 rad → sin(α) ≈ 0.00785
# So expected_count = N_sites × sin(d/R)

STRIP_FRACTION = np.sin(THRESHOLD_KM / EARTH_R_KM)

def poisson_z(observed, expected):
    """Poisson-approximation Z-score."""
    if expected <= 0:
        return 0.0
    return (observed - expected) / np.sqrt(expected)

# ============================================================
# LOAD PLEIADES DATA
# ============================================================
print("=" * 70)
print("LOADING PLEIADES DATA")
print("=" * 70)

pleiades_path = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")
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
            is_mon = bool(feature_types & MONUMENTAL_TYPES)
            is_set = bool(feature_types & SETTLEMENT_TYPES)
            if is_mon:
                mon_list.append((lat, lon))
            elif is_set:
                set_list.append((lat, lon))
        except (ValueError, KeyError, TypeError):
            pass

mon_arr = np.array(mon_list)
set_arr = np.array(set_list)
mon_lats, mon_lons = mon_arr[:, 0], mon_arr[:, 1]
set_lats, set_lons = set_arr[:, 0], set_arr[:, 1]

N_MON = len(mon_arr)
N_SET = len(set_arr)
EXPECTED_MON = N_MON * STRIP_FRACTION
EXPECTED_SET = N_SET * STRIP_FRACTION

print(f"Monumental sites: {N_MON}")
print(f"Settlement sites: {N_SET}")
print(f"Strip fraction (50km): {STRIP_FRACTION:.6f}")
print(f"Expected monuments per circle: {EXPECTED_MON:.1f}")
print(f"Expected settlements per circle: {EXPECTED_SET:.1f}")

# ============================================================
# PHASE 1: COARSE GRID SCAN
# ============================================================
print(f"\n{'='*70}")
print(f"PHASE 1: COARSE GRID SCAN ({COARSE_STEP}° resolution)")
print(f"{'='*70}")

# Generate pole grid — one hemisphere only (antipodal pole = same circle)
# For lat > 0, scan all longitudes. For lat = 0, scan lon >= 0 only.
# For lat < 0, skip (handled by antipodal symmetry).
poles = []
for lat in np.arange(-90, 90 + COARSE_STEP, COARSE_STEP):
    for lon in np.arange(-180, 180, COARSE_STEP):
        # Antipodal symmetry: (lat, lon) and (-lat, lon+180) define the same GC
        # Keep only one: require lat > 0, or lat == 0 and lon >= 0
        if lat > 0:
            poles.append((lat, lon))
        elif lat == 0 and lon >= 0:
            poles.append((lat, lon))
        # lat < 0: skip (antipodal partner at lat > 0 covers this)

# Special case: lat = 90 (north pole) — only one point needed
# lat = -90 is antipodal to lat=90, so skip it
# Actually the loop already handles this: lat=90 scans all lons but they're
# all the same circle (equator). Let's deduplicate: lat=±90 → one entry.
poles_dedup = []
seen_polar = False
for lat, lon in poles:
    if abs(lat) == 90:
        if not seen_polar:
            poles_dedup.append((lat, lon))
            seen_polar = True
    else:
        poles_dedup.append((lat, lon))
poles = poles_dedup

print(f"Total pole positions: {len(poles)}")

def scan_one_pole(args):
    """Worker function for parallel coarse scan."""
    pole_lat, pole_lon, mon_lats, mon_lons, set_lats, set_lons = args
    mon_count = count_within_vec(mon_lats, mon_lons, pole_lat, pole_lon)
    set_count = count_within_vec(set_lats, set_lons, pole_lat, pole_lon)
    z_mon = poisson_z(mon_count, EXPECTED_MON)
    z_set = poisson_z(set_count, EXPECTED_SET)
    D = z_mon - z_set
    return {
        'pole_lat': float(pole_lat),
        'pole_lon': float(pole_lon),
        'mon_count': int(mon_count),
        'set_count': int(set_count),
        'z_mon': float(z_mon),
        'z_set': float(z_set),
        'D': float(D)
    }

# Run parallel scan
t0 = time.time()
n_workers = min(cpu_count(), 10)
print(f"Using {n_workers} workers")

# Prepare args — pass arrays as shared data via global (pickle-friendly)
# Actually multiprocessing needs picklable args. Let's use a simpler approach:
# process in chunks sequentially but vectorize within each pole.
# Given ~8K poles and vectorized ops, sequential should be fast enough.

results_phase1 = []
chunk_size = 100
for i in range(0, len(poles), chunk_size):
    chunk = poles[i:i+chunk_size]
    for pole_lat, pole_lon in chunk:
        mon_count = count_within_vec(mon_lats, mon_lons, pole_lat, pole_lon)
        set_count = count_within_vec(set_lats, set_lons, pole_lat, pole_lon)
        z_mon = poisson_z(mon_count, EXPECTED_MON)
        z_set = poisson_z(set_count, EXPECTED_SET)
        D = z_mon - z_set
        results_phase1.append({
            'pole_lat': float(pole_lat),
            'pole_lon': float(pole_lon),
            'mon_count': int(mon_count),
            'set_count': int(set_count),
            'z_mon': float(z_mon),
            'z_set': float(z_set),
            'D': float(D)
        })
    elapsed = time.time() - t0
    done = min(i + chunk_size, len(poles))
    rate = done / elapsed if elapsed > 0 else 0
    eta = (len(poles) - done) / rate if rate > 0 else 0
    print(f"  [{done:5d}/{len(poles)}] elapsed={elapsed:.0f}s  ETA={eta:.0f}s", flush=True)

elapsed = time.time() - t0
print(f"\nPhase 1 complete: {len(results_phase1)} poles scanned in {elapsed:.0f}s")

# Save coarse scan
with open(os.path.join(OUT_DIR, "coarse_scan.json"), 'w') as f:
    json.dump({
        'config': {
            'step_deg': COARSE_STEP,
            'threshold_km': THRESHOLD_KM,
            'n_monuments': N_MON,
            'n_settlements': N_SET,
            'strip_fraction': float(STRIP_FRACTION),
            'expected_mon': float(EXPECTED_MON),
            'expected_set': float(EXPECTED_SET),
            'method': 'poisson_approximation'
        },
        'results': results_phase1
    }, f, indent=1)

# ============================================================
# FIND ALISON POLE IN COARSE GRID
# ============================================================
# Find nearest grid point to Alison pole
alison_dists = [haversine_deg(r['pole_lat'], r['pole_lon'], ALISON_POLE_LAT, ALISON_POLE_LON)
                for r in results_phase1]
alison_idx = int(np.argmin(alison_dists))
alison_coarse = results_phase1[alison_idx]
print(f"\nAlison pole nearest grid point: ({alison_coarse['pole_lat']}, {alison_coarse['pole_lon']})")
print(f"  D_approx = {alison_coarse['D']:.2f}  (Z_mon={alison_coarse['z_mon']:.2f}, Z_set={alison_coarse['z_set']:.2f})")
print(f"  Mon count = {alison_coarse['mon_count']}, Set count = {alison_coarse['set_count']}")

# ============================================================
# PHASE 2: FIND HOT SPOTS AND REFINE
# ============================================================
print(f"\n{'='*70}")
print(f"PHASE 2: FINE GRID REFINEMENT ({FINE_STEP}° resolution)")
print(f"{'='*70}")

# Sort by D descending
sorted_results = sorted(results_phase1, key=lambda r: r['D'], reverse=True)

# Find threshold: top 1% or D > 3, whichever captures more
d_values = np.array([r['D'] for r in results_phase1])
d_p99 = np.percentile(d_values, 99)
d_threshold = min(d_p99, 3.0)
print(f"D percentile 99: {d_p99:.2f}")
print(f"Using D threshold: {d_threshold:.2f}")

# Filter: D above threshold AND at least some monuments (exclude empty ocean circles)
# Empty circles have mon_count=0, set_count=0 → artificial D from Poisson formula
hot_candidates = [r for r in results_phase1
                  if r['D'] >= d_threshold and (r['mon_count'] > 0 or r['set_count'] > 0)]
print(f"Candidates above threshold (with sites): {len(hot_candidates)}")
print(f"  (Excluded {sum(1 for r in results_phase1 if r['D'] >= d_threshold and r['mon_count'] == 0 and r['set_count'] == 0)} empty-ocean circles)")

# Cluster into distinct hot spots (>HOT_SPOT_SEP degrees apart)
hot_spots = []
for cand in sorted(hot_candidates, key=lambda r: r['D'], reverse=True):
    is_new = True
    for hs in hot_spots:
        if haversine_deg(cand['pole_lat'], cand['pole_lon'],
                        hs['pole_lat'], hs['pole_lon']) < HOT_SPOT_SEP:
            is_new = False
            break
    if is_new:
        hot_spots.append(cand)

print(f"Distinct hot spots (>{HOT_SPOT_SEP}° apart): {len(hot_spots)}")
for i, hs in enumerate(hot_spots):
    ang = haversine_deg(hs['pole_lat'], hs['pole_lon'], ALISON_POLE_LAT, ALISON_POLE_LON)
    print(f"  #{i+1}: pole=({hs['pole_lat']:.1f}, {hs['pole_lon']:.1f})  D={hs['D']:.2f}  "
          f"ang_from_Alison={ang:.1f}°")

# Refine each hot spot with MC Z-scores
print(f"\nRefining {len(hot_spots)} hot spots with MC Z-scores ({N_MC_FINE} trials)...")
fine_results = []
t0_fine = time.time()

for hi, hs in enumerate(hot_spots):
    center_lat, center_lon = hs['pole_lat'], hs['pole_lon']
    best_D = -999
    best_result = None
    sub_results = []

    lat_range = np.arange(center_lat - FINE_BOX, center_lat + FINE_BOX + FINE_STEP, FINE_STEP)
    lon_range = np.arange(center_lon - FINE_BOX, center_lon + FINE_BOX + FINE_STEP, FINE_STEP)
    # Clip latitude
    lat_range = lat_range[(lat_range >= -90) & (lat_range <= 90)]

    n_sub = len(lat_range) * len(lon_range)
    done_sub = 0

    for plat in lat_range:
        for plon in lon_range:
            # Normalize longitude
            plon_n = ((plon + 180) % 360) - 180
            z_mon, obs_mon, mu_mon, sig_mon = mc_zscore_vec(mon_lats, mon_lons, plat, plon_n, N_MC_FINE)
            z_set, obs_set, mu_set, sig_set = mc_zscore_vec(set_lats, set_lons, plat, plon_n, N_MC_FINE)
            D = z_mon - z_set
            result = {
                'pole_lat': float(plat),
                'pole_lon': float(plon_n),
                'z_mon': float(z_mon), 'z_set': float(z_set), 'D': float(D),
                'mon_count': int(obs_mon), 'set_count': int(obs_set),
                'mon_baseline_mean': float(mu_mon), 'set_baseline_mean': float(mu_set),
                'hot_spot_id': hi
            }
            sub_results.append(result)
            if D > best_D:
                best_D = D
                best_result = result
            done_sub += 1

        # Progress per lat row
        elapsed_fine = time.time() - t0_fine
        print(f"  Hot spot #{hi+1}: [{done_sub}/{n_sub}] sub-poles  "
              f"best D so far={best_D:.2f}  elapsed={elapsed_fine:.0f}s", flush=True)

    fine_results.extend(sub_results)
    ang = haversine_deg(best_result['pole_lat'], best_result['pole_lon'],
                        ALISON_POLE_LAT, ALISON_POLE_LON)
    print(f"  Hot spot #{hi+1} peak: pole=({best_result['pole_lat']:.2f}, {best_result['pole_lon']:.2f})  "
          f"D={best_D:.2f}  ang_from_Alison={ang:.1f}°")

elapsed_fine = time.time() - t0_fine
print(f"\nPhase 2 complete in {elapsed_fine:.0f}s  ({len(fine_results)} sub-poles evaluated)")

# Save fine scan
with open(os.path.join(OUT_DIR, "fine_scan.json"), 'w') as f:
    json.dump({
        'config': {
            'fine_step_deg': FINE_STEP,
            'fine_box_deg': FINE_BOX,
            'n_mc_trials': N_MC_FINE,
            'threshold_km': THRESHOLD_KM,
            'n_hot_spots': len(hot_spots),
            'd_threshold': float(d_threshold),
            'method': 'distribution_matched_mc'
        },
        'hot_spots': hot_spots,
        'results': fine_results
    }, f, indent=1)

# ============================================================
# PHASE 3: ANALYSIS
# ============================================================
print(f"\n{'='*70}")
print(f"PHASE 3: ANALYSIS")
print(f"{'='*70}")

# --- Alison ranking ---
# Compute Alison's exact D with MC
print("Computing Alison pole exact MC D...")
z_mon_alison, obs_mon_a, mu_mon_a, sig_mon_a = mc_zscore_vec(mon_lats, mon_lons,
    ALISON_POLE_LAT, ALISON_POLE_LON, N_MC_FINE)
z_set_alison, obs_set_a, mu_set_a, sig_set_a = mc_zscore_vec(set_lats, set_lons,
    ALISON_POLE_LAT, ALISON_POLE_LON, N_MC_FINE)
D_alison = z_mon_alison - z_set_alison
print(f"  Alison D = {D_alison:.2f}  (Z_mon={z_mon_alison:.2f}, Z_set={z_set_alison:.2f})")
print(f"  Mon: obs={obs_mon_a}, baseline={mu_mon_a:.1f}±{sig_mon_a:.1f}")
print(f"  Set: obs={obs_set_a}, baseline={mu_set_a:.1f}±{sig_set_a:.1f}")

# Rank in coarse scan (Poisson D)
d_values_sorted = sorted(d_values, reverse=True)
alison_coarse_D = alison_coarse['D']
alison_rank_coarse = int(np.sum(d_values >= alison_coarse_D))
print(f"\nAlison coarse rank: #{alison_rank_coarse} out of {len(d_values)}")
print(f"  Percentile: {100*(1 - alison_rank_coarse/len(d_values)):.1f}%")

# Rank among fine-scan peaks
fine_peaks = []
for hi in range(len(hot_spots)):
    hs_results = [r for r in fine_results if r['hot_spot_id'] == hi]
    if hs_results:
        peak = max(hs_results, key=lambda r: r['D'])
        fine_peaks.append(peak)

fine_peaks_sorted = sorted(fine_peaks, key=lambda r: r['D'], reverse=True)

# Top 10 circles
top10 = fine_peaks_sorted[:10]
print(f"\nTop 10 circles by refined D:")
for i, p in enumerate(top10):
    ang = haversine_deg(p['pole_lat'], p['pole_lon'], ALISON_POLE_LAT, ALISON_POLE_LON)
    print(f"  #{i+1}: pole=({p['pole_lat']:.2f}, {p['pole_lon']:.2f})  D={p['D']:.2f}  "
          f"(Zm={p['z_mon']:.1f}, Zs={p['z_set']:.1f})  ang_from_Alison={ang:.1f}°")

with open(os.path.join(OUT_DIR, "top10_circles.json"), 'w') as f:
    json.dump(top10, f, indent=2)

# Alison rank info
alison_rank_info = {
    'alison_pole': {'lat': ALISON_POLE_LAT, 'lon': ALISON_POLE_LON},
    'alison_mc_D': D_alison,
    'alison_mc_z_mon': z_mon_alison,
    'alison_mc_z_set': z_set_alison,
    'alison_mon_count': obs_mon_a,
    'alison_set_count': obs_set_a,
    'coarse_rank': alison_rank_coarse,
    'coarse_total': len(d_values),
    'coarse_percentile': float(100*(1 - alison_rank_coarse/len(d_values))),
    'coarse_nearest_D': alison_coarse_D,
    'n_fine_peaks': len(fine_peaks),
    'fine_peaks_summary': [
        {'rank': i+1, 'pole_lat': p['pole_lat'], 'pole_lon': p['pole_lon'],
         'D': p['D'], 'z_mon': p['z_mon'], 'z_set': p['z_set'],
         'ang_from_alison': haversine_deg(p['pole_lat'], p['pole_lon'],
                                          ALISON_POLE_LAT, ALISON_POLE_LON)}
        for i, p in enumerate(fine_peaks_sorted[:20])
    ]
}
with open(os.path.join(OUT_DIR, "alison_rank.json"), 'w') as f:
    json.dump(alison_rank_info, f, indent=2)

# --- European bias test ---
print("\n--- European Bias Test ---")
EUROPE_LAT = (35, 70)
EUROPE_LON = (-10, 40)

def gc_passes_through_region(pole_lat, pole_lon, lat_range, lon_range, n_samples=360):
    """Check if a great circle passes through a lat/lon box."""
    # Parametrize the great circle
    pole = np.array([
        np.cos(np.radians(pole_lat)) * np.cos(np.radians(pole_lon)),
        np.cos(np.radians(pole_lat)) * np.sin(np.radians(pole_lon)),
        np.sin(np.radians(pole_lat))
    ])
    # Find two orthonormal vectors in the GC plane
    if abs(pole[2]) < 0.9:
        u = np.cross(pole, [0, 0, 1])
    else:
        u = np.cross(pole, [1, 0, 0])
    u = u / np.linalg.norm(u)
    v = np.cross(pole, u)
    v = v / np.linalg.norm(v)

    thetas = np.linspace(0, 2*np.pi, n_samples, endpoint=False)
    points = np.outer(np.cos(thetas), u) + np.outer(np.sin(thetas), v)
    lats = np.degrees(np.arcsin(np.clip(points[:, 2], -1, 1)))
    lons = np.degrees(np.arctan2(points[:, 1], points[:, 0]))

    in_region = ((lats >= lat_range[0]) & (lats <= lat_range[1]) &
                 (lons >= lon_range[0]) & (lons <= lon_range[1]))
    return bool(np.any(in_region))

# Check each coarse result
europe_results = []
non_europe_results = []
for r in results_phase1:
    if gc_passes_through_region(r['pole_lat'], r['pole_lon'], EUROPE_LAT, EUROPE_LON):
        europe_results.append(r)
    else:
        non_europe_results.append(r)

eu_D = [r['D'] for r in europe_results]
non_eu_D = [r['D'] for r in non_europe_results]
alison_through_europe = gc_passes_through_region(ALISON_POLE_LAT, ALISON_POLE_LON,
                                                   EUROPE_LAT, EUROPE_LON)

print(f"Circles through Europe: {len(europe_results)}  mean D={np.mean(eu_D):.2f}  max D={np.max(eu_D):.2f}")
print(f"Circles NOT through Europe: {len(non_europe_results)}  mean D={np.mean(non_eu_D):.2f}  max D={np.max(non_eu_D) if non_eu_D else 'N/A':.2f}")
print(f"Alison circle through Europe: {alison_through_europe}")

# ============================================================
# VISUALIZATIONS
# ============================================================
print(f"\n{'='*70}")
print("GENERATING VISUALIZATIONS")
print(f"{'='*70}")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# --- Heatmap helper ---
def plot_heatmap(results, value_key, title, filename, cmap='RdBu_r', mark_alison=True):
    """Plot a Mollweide projection heatmap of pole positions."""
    fig, ax = plt.subplots(figsize=(14, 7), subplot_kw={'projection': 'mollweide'})
    ax.set_title(title, fontsize=14, pad=20)

    lats = np.array([r['pole_lat'] for r in results])
    lons = np.array([r['pole_lon'] for r in results])
    vals = np.array([r[value_key] for r in results])

    # Convert to radians for Mollweide
    lats_r = np.radians(lats)
    lons_r = np.radians(lons)

    # Color normalization centered at 0
    vmax = max(abs(vals.min()), abs(vals.max()))
    if vals.min() < 0 and vals.max() > 0:
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    else:
        norm = None

    sc = ax.scatter(lons_r, lats_r, c=vals, cmap=cmap, s=3, alpha=0.7, norm=norm)
    plt.colorbar(sc, ax=ax, shrink=0.6, label=value_key)

    if mark_alison:
        ax.plot(np.radians(ALISON_POLE_LON), np.radians(ALISON_POLE_LAT),
                'k*', markersize=15, markeredgecolor='white', markeredgewidth=0.5,
                label='Alison pole', zorder=10)
        ax.legend(loc='lower left', fontsize=10)

    # Mark top 5 peaks
    top5_idx = np.argsort(vals)[-5:]
    for idx in top5_idx:
        ax.plot(lons_r[idx], lats_r[idx], 'o', markersize=8,
                markerfacecolor='none', markeredgecolor='black', markeredgewidth=1.5)

    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {filename}")

# Plot all three heatmaps
plot_heatmap(results_phase1, 'D', 'Monument-Settlement Divergence D (Poisson, 50km)',
             'global_d_heatmap.png')
plot_heatmap(results_phase1, 'z_mon', 'Monument Z-score (Poisson, 50km)',
             'monument_z_heatmap.png', cmap='Reds')
plot_heatmap(results_phase1, 'z_set', 'Settlement Z-score (Poisson, 50km)',
             'settlement_z_heatmap.png', cmap='Blues')

# --- Top 10 globe plot ---
def trace_great_circle(pole_lat, pole_lon, n_points=360):
    """Return (lats, lons) arrays tracing the great circle."""
    pole = np.array([
        np.cos(np.radians(pole_lat)) * np.cos(np.radians(pole_lon)),
        np.cos(np.radians(pole_lat)) * np.sin(np.radians(pole_lon)),
        np.sin(np.radians(pole_lat))
    ])
    if abs(pole[2]) < 0.9:
        u = np.cross(pole, [0, 0, 1])
    else:
        u = np.cross(pole, [1, 0, 0])
    u = u / np.linalg.norm(u)
    v = np.cross(pole, u)
    v = v / np.linalg.norm(v)
    thetas = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    points = np.outer(np.cos(thetas), u) + np.outer(np.sin(thetas), v)
    lats = np.degrees(np.arcsin(np.clip(points[:, 2], -1, 1)))
    lons = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    return lats, lons

fig, ax = plt.subplots(figsize=(14, 7), subplot_kw={'projection': 'mollweide'})
ax.set_title('Top 10 Great Circles by Divergence D', fontsize=14, pad=20)
ax.grid(True, alpha=0.3)

colors = plt.cm.tab10(np.linspace(0, 1, 10))
for i, p in enumerate(top10[:10]):
    gc_lats, gc_lons = trace_great_circle(p['pole_lat'], p['pole_lon'])
    # Sort by longitude for cleaner plotting, split at discontinuities
    gc_lons_r = np.radians(gc_lons)
    gc_lats_r = np.radians(gc_lats)
    # Plot as scatter to avoid wrapping artifacts
    ax.scatter(gc_lons_r, gc_lats_r, c=[colors[i]], s=0.5, alpha=0.6,
               label=f"#{i+1} D={p['D']:.1f}")

# Plot Alison GC
gc_lats, gc_lons = trace_great_circle(ALISON_POLE_LAT, ALISON_POLE_LON)
ax.scatter(np.radians(gc_lons), np.radians(gc_lats), c='black', s=1.5, alpha=0.9,
           label=f"Alison D={D_alison:.1f}", zorder=5)

ax.legend(loc='lower left', fontsize=7, ncol=2, markerscale=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "top10_globe.png"), dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"  Saved top10_globe.png")

# ============================================================
# D DISTRIBUTION HISTOGRAM
# ============================================================
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(d_values, bins=100, color='steelblue', alpha=0.7, edgecolor='black', linewidth=0.3)
ax.axvline(alison_coarse_D, color='red', linewidth=2, linestyle='--',
           label=f'Alison D={alison_coarse_D:.2f}')
if len(fine_peaks_sorted) > 0:
    ax.axvline(fine_peaks_sorted[0]['D'], color='green', linewidth=2, linestyle=':',
               label=f'Global max D={fine_peaks_sorted[0]["D"]:.2f}')
ax.set_xlabel('Divergence D', fontsize=12)
ax.set_ylabel('Count', fontsize=12)
ax.set_title('Distribution of D across all pole positions (2° grid)', fontsize=13)
ax.legend(fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "d_distribution.png"), dpi=150)
plt.close(fig)
print(f"  Saved d_distribution.png")

# ============================================================
# WRITE RESULTS.md
# ============================================================
print(f"\n{'='*70}")
print("WRITING RESULTS.md")
print(f"{'='*70}")

# Count how many have D > various thresholds
d5_count = int(np.sum(d_values > 5))
d3_count = int(np.sum(d_values > 3))
d_max = float(np.max(d_values))

# Determine verdict
if alison_rank_coarse <= 3:
    verdict = "EXCEPTIONAL — The Alison pole is in the top 3 globally. The geometry is genuinely unusual."
elif alison_rank_coarse <= 10:
    verdict = "NOTABLE — The Alison pole is in the top 10 but not unique. Other hot spots exist."
elif d5_count > 20:
    verdict = "GEOGRAPHIC — Many circles show high D. The divergence may reflect Earth's geography."
else:
    verdict = "UNIQUE — Few or no other circles match the Alison divergence."

results_md = f"""# Directive G: Systematic Great Circle Search — Results

## Summary
Scanned **{len(poles)}** pole positions at {COARSE_STEP}° resolution, computing monument-settlement
divergence D for each great circle using Pleiades data ({N_MON} monuments, {N_SET} settlements).

## Key Findings

### Alison Pole Ranking
- **Coarse grid (Poisson Z):** Rank #{alison_rank_coarse} of {len(d_values)} ({100*(1-alison_rank_coarse/len(d_values)):.1f}th percentile)
- **Nearest grid point D:** {alison_coarse_D:.2f}
- **Exact MC D:** {D_alison:.2f} (Z_mon={z_mon_alison:.2f}, Z_set={z_set_alison:.2f})
- Monument count: {obs_mon_a} (baseline: {mu_mon_a:.1f}±{sig_mon_a:.1f})
- Settlement count: {obs_set_a} (baseline: {mu_set_a:.1f}±{sig_set_a:.1f})

### Global Distribution of D
- Mean D: {np.mean(d_values):.2f}
- Std D: {np.std(d_values):.2f}
- Max D: {d_max:.2f}
- Poles with D > 5: {d5_count}
- Poles with D > 3: {d3_count}

### Hot Spots Found
{len(hot_spots)} distinct peaks identified (>{HOT_SPOT_SEP}° separation):

| Rank | Pole Lat | Pole Lon | D (refined) | Z_mon | Z_set | Ang from Alison |
|------|----------|----------|-------------|-------|-------|-----------------|
"""

for i, p in enumerate(fine_peaks_sorted[:20]):
    ang = haversine_deg(p['pole_lat'], p['pole_lon'], ALISON_POLE_LAT, ALISON_POLE_LON)
    results_md += f"| {i+1} | {p['pole_lat']:.2f} | {p['pole_lon']:.2f} | {p['D']:.2f} | {p['z_mon']:.1f} | {p['z_set']:.1f} | {ang:.1f}° |\n"

results_md += f"""
### European Bias Test
- Circles through Europe: {len(europe_results)} (mean D={np.mean(eu_D):.2f}, max D={np.max(eu_D):.2f})
- Circles NOT through Europe: {len(non_europe_results)} (mean D={np.mean(non_eu_D):.2f}, max D={np.max(non_eu_D) if non_eu_D else 'N/A':.2f})
- Alison circle passes through Europe: {alison_through_europe}

### Verdict
**{verdict}**

## Files
- `coarse_scan.json` — D values for all {len(poles)} poles (Poisson Z)
- `fine_scan.json` — Refined D values around {len(hot_spots)} hot spots (MC Z, {N_MC_FINE} trials)
- `global_d_heatmap.png` — Global D heatmap (Mollweide)
- `monument_z_heatmap.png` — Monument Z-score heatmap
- `settlement_z_heatmap.png` — Settlement Z-score heatmap
- `top10_circles.json` — Properties of top 10 circles
- `top10_globe.png` — Top 10 circles on globe
- `alison_rank.json` — Alison ranking details
- `d_distribution.png` — Histogram of D values

## Methodology
- **Phase 1**: {COARSE_STEP}° grid, Poisson-approximation Z-scores
  - Strip fraction = sin({THRESHOLD_KM}/{EARTH_R_KM:.0f}) = {STRIP_FRACTION:.6f}
  - Expected monuments per circle: {EXPECTED_MON:.1f}
  - Expected settlements per circle: {EXPECTED_SET:.1f}
- **Phase 2**: {FINE_STEP}° refinement in ±{FINE_BOX}° boxes around hot spots
  - Distribution-matched MC, {N_MC_FINE} trials per pole
- **Data**: Pleiades Gazetteer (all periods), {N_MON} monumental + {N_SET} settlement sites
"""

with open(os.path.join(OUT_DIR, "RESULTS.md"), 'w') as f:
    f.write(results_md)
print("  Saved RESULTS.md")

print(f"\n{'='*70}")
print("DIRECTIVE G COMPLETE")
print(f"{'='*70}")
