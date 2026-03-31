#!/usr/bin/env python3
"""
Directive 16: Coastal Bias Correction for Distribution-Matched Monte Carlo
===========================================================================
Implements land-constrained MC baseline and runs all 5 analyses.

Analyses:
  1. Land-constrained distribution-matched baseline (N=1000 and N=10000)
  2. Random circles cross-validation (Approach B)
  3. Easter Island isolation test
  4. MC sample size convergence
  5. Corrected summary for paper revision
"""

import csv
import json
import math
import os
import random
import time
import sys

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from shapely.prepared import prep

# ============================================================
# CONSTANTS
# ============================================================
POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2  # ~10,007.5 km

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "coastal_bias_correction")
os.makedirs(OUT_DIR, exist_ok=True)

PLEIADES_PATH = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")
LAND_MASK_PATH = os.path.join(BASE_DIR, "data", "land_mask", "ne_10m_land.shp")

# Pleiades type classifications (from Paper 1)
PLEIADES_MONUMENTAL = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus", "shrine",
}

PLEIADES_SETTLEMENT = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement", "townhouse", "production",
}

SEED = 42

# ============================================================
# LAND MASK — Direct prepared geometry (no rasterization needed)
# ============================================================
print("Loading Natural Earth land mask...")
t0 = time.time()
_land_gdf = gpd.read_file(LAND_MASK_PATH)
_land_union = _land_gdf.union_all()
_land_prepared = prep(_land_union)
print(f"  Land mask loaded and prepared in {time.time()-t0:.1f}s")

# Verify with known points
for label, lat, lon, expected in [
    ("Giza", 29.98, 31.13, True),
    ("Pacific Ocean", 0.0, -150.0, False),
    ("Easter Island", -27.12, -109.35, True),
    ("London", 51.51, -0.12, True),
]:
    result = _land_prepared.contains(Point(lon, lat))
    status = "OK" if result == expected else "MISMATCH"
    print(f"  Verify {label}: {result} (expected {expected}) [{status}]")

# Benchmark speed
t0 = time.time()
_test_pts = [(random.uniform(-90, 90), random.uniform(-180, 180)) for _ in range(10000)]
for lat, lon in _test_pts:
    _land_prepared.contains(Point(lon, lat))
rate = 10000 / (time.time() - t0)
print(f"  Land check speed: {rate:.0f} points/sec")
# For 185 sites × 10k trials × ~3 attempts avg = ~5.5M lookups
# At ~50k/sec this is ~2 minutes per MC run. Acceptable.
est_lookups = 185 * 10000 * 3
est_time = est_lookups / rate
print(f"  Estimated time for 10k-trial MC with 185 sites: {est_time/60:.1f} min")


def is_on_land(lat, lon):
    """Check if a point is on land using prepared geometry."""
    return _land_prepared.contains(Point(lon, lat))


# ============================================================
# GEODESIC MATH
# ============================================================
def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))


def gc_distance(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    return abs(haversine_km(lat, lon, pole_lat, pole_lon) - QUARTER_CIRC)


def gc_distances_vectorized(lats, lons, pole_lat, pole_lon):
    """Vectorized great circle distance computation."""
    lats_r = np.radians(lats)
    lons_r = np.radians(lons)
    pole_lat_r = math.radians(pole_lat)
    pole_lon_r = math.radians(pole_lon)
    px = math.cos(pole_lat_r) * math.cos(pole_lon_r)
    py = math.cos(pole_lat_r) * math.sin(pole_lon_r)
    pz = math.sin(pole_lat_r)
    sx = np.cos(lats_r) * np.cos(lons_r)
    sy = np.cos(lats_r) * np.sin(lons_r)
    sz = np.sin(lats_r)
    dots = np.clip(sx*px + sy*py + sz*pz, -1.0, 1.0)
    angular_dist = np.arccos(dots)
    return np.abs(angular_dist - math.pi/2) * EARTH_R_KM


def count_within_corridor(lats, lons, pole_lat, pole_lon, threshold_km=50):
    dists = gc_distances_vectorized(np.array(lats), np.array(lons), pole_lat, pole_lon)
    return int(np.sum(dists <= threshold_km))


def random_pole():
    z = random.uniform(-1, 1)
    lat = math.degrees(math.asin(z))
    lon = random.uniform(-180, 180)
    return lat, lon


# ============================================================
# JITTER FUNCTIONS
# ============================================================
def jitter_unconstrained(lat, lon, sigma_deg=2.0):
    """Original jitter (no land constraint)."""
    new_lat = lat + random.gauss(0, sigma_deg)
    new_lon = lon + random.gauss(0, sigma_deg)
    new_lat = max(-90, min(90, new_lat))
    new_lon = ((new_lon + 180) % 360) - 180
    return new_lat, new_lon


def jitter_to_land(lat, lon, sigma_deg=2.0, max_attempts=100):
    """Jitter with land constraint — reject ocean placements."""
    for attempt in range(max_attempts):
        new_lat = lat + random.gauss(0, sigma_deg)
        new_lon = lon + random.gauss(0, sigma_deg)
        new_lat = max(-90, min(90, new_lat))
        new_lon = ((new_lon + 180) % 360) - 180
        if is_on_land(new_lat, new_lon):
            return new_lat, new_lon, False  # (lat, lon, fallback_used)
    # Fallback: return original (conservative — does NOT inflate z-score)
    return lat, lon, True


# ============================================================
# DATA LOADING
# ============================================================
def load_pleiades():
    print(f"Loading Pleiades from {PLEIADES_PATH}...")
    sites = []
    with open(PLEIADES_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row.get('reprLat', ''))
                lon = float(row.get('reprLong', ''))
                if -90 <= lat <= 90 and -180 <= lon <= 180:
                    min_date = None
                    try:
                        min_date = int(row.get('minDate', ''))
                    except (ValueError, TypeError):
                        pass
                    feature_types = row.get('featureTypes', '')
                    sites.append({
                        'lat': lat, 'lon': lon,
                        'name': row.get('title', ''),
                        'feature_types': feature_types,
                        'min_date': min_date,
                    })
            except (ValueError, KeyError):
                continue
    print(f"  Loaded {len(sites)} Pleiades sites")
    return sites


def classify_pleiades(feature_types_str):
    if not feature_types_str:
        return 'other'
    types = set(t.strip().lower() for t in feature_types_str.split(','))
    if types & PLEIADES_MONUMENTAL:
        return 'monumental'
    if types & PLEIADES_SETTLEMENT:
        return 'settlement'
    return 'other'


def split_by_category(sites):
    monumental = [s for s in sites if classify_pleiades(s['feature_types']) == 'monumental']
    settlement = [s for s in sites if classify_pleiades(s['feature_types']) == 'settlement']
    ancient_mon = [s for s in monumental if s.get('min_date') is not None and s['min_date'] <= -2000]
    ancient_set = [s for s in settlement if s.get('min_date') is not None and s['min_date'] <= -2000]
    return {
        'monumental': monumental,
        'settlement': settlement,
        'monumental_ancient': ancient_mon,
        'settlement_ancient': ancient_set,
    }


# ============================================================
# MC BASELINES
# ============================================================
def run_mc_unconstrained(site_lats, site_lons, n_trials, threshold_km=50, sigma_deg=2.0):
    """Original distribution-matched baseline (no land constraint)."""
    n = len(site_lats)
    counts = []
    for _ in range(n_trials):
        syn_lats = []
        syn_lons = []
        for i in range(n):
            lat = random.choice(site_lats) + random.gauss(0, sigma_deg)
            lon = random.choice(site_lons) + random.gauss(0, sigma_deg)
            lat = max(-90, min(90, lat))
            lon = ((lon + 180) % 360) - 180
            syn_lats.append(lat)
            syn_lons.append(lon)
        c = count_within_corridor(syn_lats, syn_lons, POLE_LAT, POLE_LON, threshold_km)
        counts.append(c)
    return counts


def run_mc_land_constrained(site_lats, site_lons, n_trials, threshold_km=50,
                            sigma_deg=2.0, track_fallbacks=False):
    """Land-constrained distribution-matched baseline."""
    n = len(site_lats)
    counts = []
    total_fallbacks = 0
    fallback_sites = {}  # site_index -> count of fallbacks

    for trial in range(n_trials):
        syn_lats = []
        syn_lons = []
        for i in range(n):
            new_lat, new_lon, fell_back = jitter_to_land(
                site_lats[i], site_lons[i], sigma_deg=sigma_deg
            )
            syn_lats.append(new_lat)
            syn_lons.append(new_lon)
            if fell_back:
                total_fallbacks += 1
                if track_fallbacks:
                    fallback_sites[i] = fallback_sites.get(i, 0) + 1
        c = count_within_corridor(syn_lats, syn_lons, POLE_LAT, POLE_LON, threshold_km)
        counts.append(c)

    fb_info = {
        'total_fallbacks': total_fallbacks,
        'total_jitters': n * n_trials,
        'fallback_rate': total_fallbacks / (n * n_trials) if n_trials > 0 else 0,
        'fallback_sites': fallback_sites,
    }
    return counts, fb_info


def run_mc_random_circles(site_lats, site_lons, n_trials, threshold_km=50):
    """Approach B: Random circles, fixed sites."""
    n_sites = len(site_lats)
    lats = np.array(site_lats)
    lons = np.array(site_lons)

    # Compute latitude profile of real great circle sites for filtering
    # Get observed lat distribution in 10° bands
    real_dists = gc_distances_vectorized(lats, lons, POLE_LAT, POLE_LON)
    on_corridor_mask = real_dists <= threshold_km

    counts = []
    for _ in range(n_trials):
        # Generate random pole (uniform on sphere)
        z = random.uniform(-1, 1)
        pole_lat = math.degrees(math.asin(z))
        pole_lon = random.uniform(-180, 180)
        c = count_within_corridor(site_lats, site_lons, pole_lat, pole_lon, threshold_km)
        counts.append(c)
    return counts


def compute_z(observed, mc_counts):
    """Compute z-score from observed count and MC distribution."""
    arr = np.array(mc_counts, dtype=float)
    mu = arr.mean()
    sigma = arr.std()
    z = (observed - mu) / sigma if sigma > 0 else 0.0
    p = np.mean(arr >= observed)
    enrichment = observed / mu if mu > 0 else float('inf')
    return {
        'observed': int(observed),
        'mc_mean': round(float(mu), 2),
        'mc_std': round(float(sigma), 2),
        'z_score': round(float(z), 2),
        'enrichment': round(float(enrichment), 2),
        'p_value': round(float(p), 4),
    }


# ============================================================
# ANALYSIS 1: Land-Constrained Distribution-Matched Baseline
# ============================================================
def analysis_1(categories, threshold_km=50):
    print("\n" + "="*70)
    print("ANALYSIS 1: Land-Constrained Distribution-Matched Baseline")
    print("="*70)

    results = {}
    fallback_reports = []

    # Process ancient subsets first (critical), then all-period
    cat_order = ['monumental_ancient', 'settlement_ancient', 'monumental', 'settlement']
    for cat_name in cat_order:
        sites = categories.get(cat_name)
        if not sites:
            print(f"  Skipping {cat_name} (no sites)")
            continue

        lats = [s['lat'] for s in sites]
        lons = [s['lon'] for s in sites]
        observed = count_within_corridor(lats, lons, POLE_LAT, POLE_LON, threshold_km)
        is_ancient = 'ancient' in cat_name
        # Use fewer trials for large (all-period) categories to keep runtime manageable
        n_mc_base = 1000
        n_mc_high = 10000 if is_ancient else 1000
        print(f"\n--- {cat_name}: {len(sites)} sites, {observed} observed within {threshold_km} km ---")

        # Original unconstrained
        random.seed(SEED)
        np.random.seed(SEED)
        print(f"  Running unconstrained MC (N={n_mc_base})...")
        t0 = time.time()
        unc_counts = run_mc_unconstrained(lats, lons, n_mc_base, threshold_km)
        unc_stats = compute_z(observed, unc_counts)
        print(f"    Done in {time.time()-t0:.1f}s — Z={unc_stats['z_score']}, "
              f"enrichment={unc_stats['enrichment']}x")

        # Land-constrained (N=1000)
        random.seed(SEED)
        np.random.seed(SEED)
        print(f"  Running land-constrained MC (N={n_mc_base})...")
        t0 = time.time()
        lc_counts, fb_info = run_mc_land_constrained(
            lats, lons, n_mc_base, threshold_km, track_fallbacks=True
        )
        lc_stats = compute_z(observed, lc_counts)
        print(f"    Done in {time.time()-t0:.1f}s — Z={lc_stats['z_score']}, "
              f"enrichment={lc_stats['enrichment']}x")
        print(f"    Fallback rate: {fb_info['fallback_rate']*100:.2f}% "
              f"({fb_info['total_fallbacks']}/{fb_info['total_jitters']})")

        # Land-constrained (N=10000 for ancient, skip for all-period large categories)
        random.seed(SEED)
        np.random.seed(SEED)
        print(f"  Running land-constrained MC (N={n_mc_high})...")
        t0 = time.time()
        lc10k_counts, fb10k_info = run_mc_land_constrained(
            lats, lons, n_mc_high, threshold_km, track_fallbacks=False
        )
        lc10k_stats = compute_z(observed, lc10k_counts)
        print(f"    Done in {time.time()-t0:.1f}s — Z={lc10k_stats['z_score']}, "
              f"enrichment={lc10k_stats['enrichment']}x")

        results[cat_name] = {
            'n_sites': len(sites),
            'observed': observed,
            'unconstrained_1k': unc_stats,
            'land_constrained_1k': lc_stats,
            'land_constrained_10k': lc10k_stats,
            'fallback_rate': round(fb_info['fallback_rate'], 6),
        }

        # Track fallback sites
        for idx, count in fb_info.get('fallback_sites', {}).items():
            s = sites[idx]
            fallback_reports.append({
                'category': cat_name,
                'name': s['name'],
                'lat': s['lat'],
                'lon': s['lon'],
                'fallback_count': count,
                'fallback_rate': round(count / 1000, 4),
            })

    # Save results
    with open(os.path.join(OUT_DIR, "01_land_constrained_enrichment.json"), 'w') as f:
        json.dump(results, f, indent=2)

    # Save fallback report
    fallback_reports.sort(key=lambda x: -x['fallback_rate'])
    with open(os.path.join(OUT_DIR, "01_fallback_report.csv"), 'w', newline='') as f:
        if fallback_reports:
            writer = csv.DictWriter(f, fieldnames=fallback_reports[0].keys())
            writer.writeheader()
            writer.writerows(fallback_reports)

    # CSV comparison table
    with open(os.path.join(OUT_DIR, "01_land_constrained_enrichment.csv"), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['category', 'n_sites', 'observed',
                         'mc_mean_original', 'z_original', 'enrichment_original',
                         'mc_mean_corrected_1k', 'z_corrected_1k', 'enrichment_corrected_1k',
                         'mc_mean_corrected_10k', 'z_corrected_10k', 'enrichment_corrected_10k',
                         'fallback_rate'])
        for cat_name, r in results.items():
            writer.writerow([
                cat_name, r['n_sites'], r['observed'],
                r['unconstrained_1k']['mc_mean'], r['unconstrained_1k']['z_score'],
                r['unconstrained_1k']['enrichment'],
                r['land_constrained_1k']['mc_mean'], r['land_constrained_1k']['z_score'],
                r['land_constrained_1k']['enrichment'],
                r['land_constrained_10k']['mc_mean'], r['land_constrained_10k']['z_score'],
                r['land_constrained_10k']['enrichment'],
                r['fallback_rate'],
            ])

    print(f"\n  Saved to {OUT_DIR}/01_*")
    return results


# ============================================================
# ANALYSIS 2: Random Circles Cross-Validation (Approach B)
# ============================================================
def analysis_2(categories, threshold_km=50):
    print("\n" + "="*70)
    print("ANALYSIS 2: Random Circles Cross-Validation (Approach B)")
    print("="*70)

    results = {}

    for cat_name, sites in categories.items():
        if len(sites) == 0:
            continue

        lats = [s['lat'] for s in sites]
        lons = [s['lon'] for s in sites]
        observed = count_within_corridor(lats, lons, POLE_LAT, POLE_LON, threshold_km)
        print(f"\n--- {cat_name}: {len(sites)} sites, {observed} observed ---")

        random.seed(SEED)
        np.random.seed(SEED)
        print(f"  Running random circles MC (N=10000)...")
        t0 = time.time()
        rc_counts = run_mc_random_circles(lats, lons, 10000, threshold_km)
        rc_stats = compute_z(observed, rc_counts)
        print(f"    Done in {time.time()-t0:.1f}s — Z={rc_stats['z_score']}, "
              f"enrichment={rc_stats['enrichment']}x")

        results[cat_name] = {
            'n_sites': len(sites),
            'observed': observed,
            'random_circles_10k': rc_stats,
        }

    with open(os.path.join(OUT_DIR, "02_random_circles_enrichment.json"), 'w') as f:
        json.dump(results, f, indent=2)

    with open(os.path.join(OUT_DIR, "02_random_circles_enrichment.csv"), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['category', 'n_sites', 'observed',
                         'mc_mean', 'mc_std', 'z_score', 'enrichment', 'p_value'])
        for cat_name, r in results.items():
            rc = r['random_circles_10k']
            writer.writerow([
                cat_name, r['n_sites'], r['observed'],
                rc['mc_mean'], rc['mc_std'], rc['z_score'],
                rc['enrichment'], rc['p_value'],
            ])

    print(f"\n  Saved to {OUT_DIR}/02_*")
    return results


# ============================================================
# ANALYSIS 3: Easter Island Isolation Test
# ============================================================
def analysis_3(categories, threshold_km=50):
    print("\n" + "="*70)
    print("ANALYSIS 3: Easter Island Isolation Test")
    print("="*70)

    def is_easter_island(site):
        return -28.0 < site['lat'] < -26.5 and -110.0 < site['lon'] < -108.5

    def is_small_island(site):
        """Easter Island + other small Pacific/Mediterranean islands."""
        if is_easter_island(site):
            return True
        # Crete, Cyprus, Malta, Sardinia, Sicily handled as larger landmasses
        # Focus on genuinely small islands where jitter pushes into ocean
        lat, lon = site['lat'], site['lon']
        # Check if site is on a small island by testing if jitter mostly fails
        land_hits = 0
        for _ in range(20):
            test_lat = lat + random.gauss(0, 2.0)
            test_lon = lon + random.gauss(0, 2.0)
            if is_on_land(test_lat, test_lon):
                land_hits += 1
        return land_hits < 5  # Less than 25% of jitters hit land = small island

    def is_egyptian_coastal(site):
        """Within ~50 km of Mediterranean or Red Sea coast in Egypt."""
        lat, lon = site['lat'], site['lon']
        if not (22 < lat < 32 and 25 < lon < 37):
            return False
        # Check if many jitters would go into water
        land_hits = 0
        for _ in range(20):
            test_lat = lat + random.gauss(0, 2.0)
            test_lon = lon + random.gauss(0, 2.0)
            if is_on_land(test_lat, test_lon):
                land_hits += 1
        return land_hits < 15  # Less than 75% hit rate = significantly coastal

    results = {}

    for cat_name in ['monumental_ancient', 'settlement_ancient']:
        sites = categories[cat_name]
        if len(sites) == 0:
            continue

        lats = [s['lat'] for s in sites]
        lons = [s['lon'] for s in sites]
        observed_full = count_within_corridor(lats, lons, POLE_LAT, POLE_LON, threshold_km)

        # Identify groups
        random.seed(SEED)
        easter_idx = [i for i, s in enumerate(sites) if is_easter_island(s)]
        random.seed(SEED)
        island_idx = [i for i, s in enumerate(sites) if is_small_island(s)]
        random.seed(SEED)
        coastal_egypt_idx = [i for i, s in enumerate(sites) if is_egyptian_coastal(s)]

        print(f"\n--- {cat_name} ({len(sites)} sites) ---")
        print(f"  Easter Island sites: {len(easter_idx)}")
        print(f"  Small island sites: {len(island_idx)}")
        print(f"  Coastal Egypt sites: {len(coastal_egypt_idx)}")

        cat_results = {'n_sites': len(sites), 'observed': observed_full, 'tests': {}}

        # Test each exclusion group
        exclusion_groups = {
            'full': set(),
            'no_easter_island': set(easter_idx),
            'no_small_islands': set(island_idx),
            'no_coastal_egypt': set(coastal_egypt_idx),
            'no_all_coastal': set(island_idx) | set(coastal_egypt_idx),
        }

        for group_name, exclude_set in exclusion_groups.items():
            subset = [s for i, s in enumerate(sites) if i not in exclude_set]
            if len(subset) == 0:
                continue
            sub_lats = [s['lat'] for s in subset]
            sub_lons = [s['lon'] for s in subset]
            sub_obs = count_within_corridor(sub_lats, sub_lons, POLE_LAT, POLE_LON, threshold_km)

            # Run land-constrained MC
            random.seed(SEED)
            np.random.seed(SEED)
            lc_counts, _ = run_mc_land_constrained(sub_lats, sub_lons, 1000, threshold_km)
            lc_stats = compute_z(sub_obs, lc_counts)

            cat_results['tests'][group_name] = {
                'n_sites': len(subset),
                'excluded': len(exclude_set),
                'observed': sub_obs,
                **lc_stats,
            }
            print(f"  {group_name}: {len(subset)} sites, obs={sub_obs}, "
                  f"Z={lc_stats['z_score']}, enrichment={lc_stats['enrichment']}x")

        results[cat_name] = cat_results

    with open(os.path.join(OUT_DIR, "03_island_isolation.json"), 'w') as f:
        json.dump(results, f, indent=2)

    with open(os.path.join(OUT_DIR, "03_island_isolation.csv"), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['category', 'exclusion_group', 'n_sites', 'excluded',
                         'observed', 'mc_mean', 'z_score', 'enrichment', 'p_value'])
        for cat_name, r in results.items():
            for group_name, t in r['tests'].items():
                writer.writerow([
                    cat_name, group_name, t['n_sites'], t['excluded'],
                    t['observed'], t['mc_mean'], t['z_score'],
                    t['enrichment'], t['p_value'],
                ])

    print(f"\n  Saved to {OUT_DIR}/03_*")
    return results


# ============================================================
# ANALYSIS 4: MC Sample Size Convergence
# ============================================================
def analysis_4(categories, threshold_km=50):
    print("\n" + "="*70)
    print("ANALYSIS 4: MC Sample Size Convergence")
    print("="*70)

    sample_sizes = [200, 500, 1000, 5000, 10000]
    results = {}

    for cat_name in ['monumental_ancient', 'settlement_ancient']:
        sites = categories[cat_name]
        if len(sites) == 0:
            continue

        lats = [s['lat'] for s in sites]
        lons = [s['lon'] for s in sites]
        observed = count_within_corridor(lats, lons, POLE_LAT, POLE_LON, threshold_km)
        print(f"\n--- {cat_name}: {len(sites)} sites, {observed} observed ---")

        cat_results = {'n_sites': len(sites), 'observed': observed, 'by_n': {}}

        for n_mc in sample_sizes:
            random.seed(SEED)
            np.random.seed(SEED)
            print(f"  N={n_mc}...", end=" ", flush=True)
            t0 = time.time()
            lc_counts, _ = run_mc_land_constrained(lats, lons, n_mc, threshold_km)
            stats = compute_z(observed, lc_counts)
            elapsed = time.time() - t0
            print(f"Z={stats['z_score']}, mean={stats['mc_mean']}, "
                  f"std={stats['mc_std']} ({elapsed:.1f}s)")
            cat_results['by_n'][str(n_mc)] = stats

        results[cat_name] = cat_results

    with open(os.path.join(OUT_DIR, "04_mc_convergence.json"), 'w') as f:
        json.dump(results, f, indent=2)

    with open(os.path.join(OUT_DIR, "04_mc_convergence.csv"), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['category', 'mc_trials', 'observed', 'mc_mean', 'mc_std',
                         'z_score', 'enrichment', 'p_value'])
        for cat_name, r in results.items():
            for n_str, stats in r['by_n'].items():
                writer.writerow([
                    cat_name, n_str, r['observed'],
                    stats['mc_mean'], stats['mc_std'], stats['z_score'],
                    stats['enrichment'], stats['p_value'],
                ])

    print(f"\n  Saved to {OUT_DIR}/04_*")
    return results


# ============================================================
# ANALYSIS 5: Corrected Summary
# ============================================================
def analysis_5(a1_results, a2_results, a3_results, a4_results):
    print("\n" + "="*70)
    print("ANALYSIS 5: Corrected Summary for Paper Revision")
    print("="*70)

    # Build comparison table
    summary = {
        'comparison_table': {},
        'key_findings': [],
    }

    for cat_name in ['monumental_ancient', 'settlement_ancient']:
        if cat_name not in a1_results:
            continue

        a1 = a1_results[cat_name]
        a2 = a2_results.get(cat_name, {})
        a4 = a4_results.get(cat_name, {})

        row = {
            'n_sites': a1['n_sites'],
            'observed': a1['observed'],
            'original_unconstrained': a1['unconstrained_1k'],
            'land_constrained_1k': a1['land_constrained_1k'],
            'land_constrained_10k': a1['land_constrained_10k'],
            'random_circles_10k': a2.get('random_circles_10k', {}),
            'fallback_rate': a1['fallback_rate'],
        }
        summary['comparison_table'][cat_name] = row

    # Compute divergences
    if ('monumental_ancient' in summary['comparison_table'] and
            'settlement_ancient' in summary['comparison_table']):
        mon = summary['comparison_table']['monumental_ancient']
        sett = summary['comparison_table']['settlement_ancient']

        divergences = {}
        for method in ['original_unconstrained', 'land_constrained_1k',
                       'land_constrained_10k', 'random_circles_10k']:
            mon_z = mon.get(method, {}).get('z_score', 0)
            set_z = sett.get(method, {}).get('z_score', 0)
            divergences[method] = round(mon_z - set_z, 2)

        summary['divergences'] = divergences

    # Key findings
    if 'monumental_ancient' in a1_results:
        ma = a1_results['monumental_ancient']
        orig_z = ma['unconstrained_1k']['z_score']
        corr_z = ma['land_constrained_10k']['z_score']
        change_pct = round((1 - corr_z/orig_z) * 100, 1) if orig_z != 0 else 0

        rc_z = a2_results.get('monumental_ancient', {}).get('random_circles_10k', {}).get('z_score', 0)

        summary['key_findings'] = [
            f"Monument Z-score drops from {orig_z} (unconstrained) to {corr_z} "
            f"(land-constrained, N=10k): {change_pct}% reduction",
            f"Random circles (Approach B) Z-score: {rc_z}",
            f"Fallback rate: {ma['fallback_rate']*100:.2f}%",
        ]

        if 'settlement_ancient' in a1_results:
            sa = a1_results['settlement_ancient']
            set_corr_z = sa['land_constrained_10k']['z_score']
            summary['key_findings'].append(
                f"Settlement Z-score (land-constrained): {set_corr_z}"
            )
            summary['key_findings'].append(
                f"Monument-settlement divergence persists: "
                f"{round(corr_z - set_corr_z, 2)} (land-constrained 10k)"
            )

    with open(os.path.join(OUT_DIR, "05_all_results.json"), 'w') as f:
        json.dump(summary, f, indent=2)

    # CSV comparison table
    with open(os.path.join(OUT_DIR, "05_corrected_statistics_table.csv"), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Statistic',
            'Original (Paper 1, N=200)',
            'Unconstrained (N=1000)',
            'Land-Constrained (N=1000)',
            'Land-Constrained (N=10000)',
            'Random Circles (N=10000)',
        ])
        for cat_name in ['monumental_ancient', 'settlement_ancient']:
            if cat_name not in summary['comparison_table']:
                continue
            r = summary['comparison_table'][cat_name]
            label = 'Monument (ancient)' if 'monumental' in cat_name else 'Settlement (ancient)'

            paper1_z = 11.83 if 'monumental' in cat_name else -0.95
            paper1_e = 5.05 if 'monumental' in cat_name else 0.78

            writer.writerow([
                f'{label} Z-score',
                paper1_z,
                r['original_unconstrained'].get('z_score', ''),
                r['land_constrained_1k'].get('z_score', ''),
                r['land_constrained_10k'].get('z_score', ''),
                r['random_circles_10k'].get('z_score', ''),
            ])
            writer.writerow([
                f'{label} enrichment',
                paper1_e,
                r['original_unconstrained'].get('enrichment', ''),
                r['land_constrained_1k'].get('enrichment', ''),
                r['land_constrained_10k'].get('enrichment', ''),
                r['random_circles_10k'].get('enrichment', ''),
            ])

        if 'divergences' in summary:
            writer.writerow([])
            writer.writerow(['Mon-Set Z divergence', '12.78',
                             summary['divergences'].get('original_unconstrained', ''),
                             summary['divergences'].get('land_constrained_1k', ''),
                             summary['divergences'].get('land_constrained_10k', ''),
                             summary['divergences'].get('random_circles_10k', '')])

    # Write narrative
    narrative = _generate_narrative(summary, a1_results, a2_results, a3_results, a4_results)
    with open(os.path.join(OUT_DIR, "05_revision_narrative.md"), 'w') as f:
        f.write(narrative)

    print(f"\n  Saved to {OUT_DIR}/05_*")
    return summary


def _generate_narrative(summary, a1, a2, a3, a4):
    lines = ["# Coastal Bias Correction: Corrected Statistics for Paper Revision\n"]
    lines.append(f"**Date:** {time.strftime('%Y-%m-%d')}\n")
    lines.append("## Summary of Corrections\n")

    if 'monumental_ancient' in a1:
        ma = a1['monumental_ancient']
        orig = ma['unconstrained_1k']
        corr = ma['land_constrained_10k']
        lines.append(
            f"The original Paper 1 distribution-matched Monte Carlo baseline (Approach A) "
            f"did not constrain synthetic site coordinates to land. When jittered with "
            f"sigma=2 degrees (~220 km), coastal and island sites frequently landed in the "
            f"ocean, deflating the expected on-corridor count and inflating z-scores.\n"
        )
        lines.append(
            f"**Monument (ancient, pre-2000 BCE) — 50 km corridor:**\n"
            f"- Original unconstrained (N=1000): Z = {orig['z_score']}, "
            f"enrichment = {orig['enrichment']}x\n"
            f"- Land-constrained (N=10,000): Z = {corr['z_score']}, "
            f"enrichment = {corr['enrichment']}x\n"
            f"- Fallback rate: {ma['fallback_rate']*100:.2f}%\n"
        )

    if 'settlement_ancient' in a1:
        sa = a1['settlement_ancient']
        corr_s = sa['land_constrained_10k']
        lines.append(
            f"**Settlement (ancient, pre-2000 BCE) — 50 km corridor:**\n"
            f"- Land-constrained (N=10,000): Z = {corr_s['z_score']}, "
            f"enrichment = {corr_s['enrichment']}x\n"
        )

    if 'monumental_ancient' in a2:
        rc = a2['monumental_ancient']['random_circles_10k']
        lines.append(
            f"\n## Random Circles Cross-Validation (Approach B)\n\n"
            f"Approach B (random great circles, fixed sites) completely sidesteps the "
            f"coastal bias issue because real sites stay at their actual coordinates.\n\n"
            f"- Monument (ancient) random circles Z = {rc['z_score']}, "
            f"enrichment = {rc['enrichment']}x\n"
        )

    if a3 and 'monumental_ancient' in a3:
        lines.append(f"\n## Easter Island and Coastal Isolation\n\n")
        tests = a3['monumental_ancient']['tests']
        for group, stats in tests.items():
            lines.append(
                f"- {group}: {stats['n_sites']} sites, Z = {stats['z_score']}, "
                f"enrichment = {stats['enrichment']}x\n"
            )

    if a4 and 'monumental_ancient' in a4:
        lines.append(f"\n## MC Convergence\n\n")
        for n_str, stats in a4['monumental_ancient']['by_n'].items():
            lines.append(f"- N = {n_str}: Z = {stats['z_score']}\n")

    # Conclusion
    lines.append(f"\n## Recommendations for Paper Revision\n\n")
    if summary.get('divergences'):
        div_lc = summary['divergences'].get('land_constrained_10k', 0)
        div_rc = summary['divergences'].get('random_circles_10k', 0)
        lines.append(
            f"Monument-settlement divergence under land-constrained correction: {div_lc}\n"
            f"Monument-settlement divergence under random circles: {div_rc}\n\n"
        )

    lines.append(
        "These corrected values should replace the original Paper 1 statistics. "
        "Paper 2 should cite the corrected values and note the methodological "
        "improvement. The core finding — that monumental sites cluster near the "
        "Great Circle significantly more than contemporaneous settlements — should "
        "be evaluated based on the corrected z-scores above.\n"
    )

    return "\n".join(lines)


# ============================================================
# PLOTTING
# ============================================================
def generate_plots(a1_results, a2_results, a4_results):
    """Generate comparison charts."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available, skipping plots")
        return

    # Plot 1: Three-way comparison bar chart
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    categories_to_plot = ['monumental_ancient', 'settlement_ancient']
    methods = ['Original\n(unconstrained)', 'Land-\nconstrained', 'Random\ncircles']
    x = np.arange(len(categories_to_plot))
    width = 0.25

    for i, method_key in enumerate(['unconstrained_1k', 'land_constrained_10k', None]):
        vals = []
        for cat in categories_to_plot:
            if method_key is None:
                # Random circles
                val = a2_results.get(cat, {}).get('random_circles_10k', {}).get('z_score', 0)
            elif cat in a1_results:
                val = a1_results[cat].get(method_key, {}).get('z_score', 0)
            else:
                val = 0
            vals.append(val)
        bars = ax.bar(x + i*width, vals, width, label=methods[i])
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
                    f'Z={val}', ha='center', va='bottom', fontsize=9)

    ax.set_ylabel('Z-score')
    ax.set_title('Coastal Bias Correction: Three-Way Comparison\n'
                 '(Ancient sites, 50 km corridor)')
    ax.set_xticks(x + width)
    ax.set_xticklabels(['Monument\n(ancient)', 'Settlement\n(ancient)'])
    ax.legend()
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.axhline(y=3.0, color='red', linewidth=0.5, linestyle='--', alpha=0.5)
    ax.text(ax.get_xlim()[1], 3.0, ' Z=3 threshold', va='center', fontsize=8, color='red')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "02_approach_comparison.png"), dpi=150)
    plt.close()
    print("  Saved 02_approach_comparison.png")

    # Plot 2: MC convergence
    if a4_results and 'monumental_ancient' in a4_results:
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        for cat_name in ['monumental_ancient', 'settlement_ancient']:
            if cat_name not in a4_results:
                continue
            ns = []
            zs = []
            for n_str, stats in sorted(a4_results[cat_name]['by_n'].items(), key=lambda x: int(x[0])):
                ns.append(int(n_str))
                zs.append(stats['z_score'])
            label = 'Monument (ancient)' if 'monumental' in cat_name else 'Settlement (ancient)'
            ax.plot(ns, zs, 'o-', label=label, markersize=6)
            for n, z in zip(ns, zs):
                ax.annotate(f'Z={z}', (n, z), textcoords="offset points",
                            xytext=(0, 10), ha='center', fontsize=8)

        ax.set_xlabel('MC Sample Size (N)')
        ax.set_ylabel('Z-score (land-constrained)')
        ax.set_title('MC Sample Size Convergence')
        ax.set_xscale('log')
        ax.legend()
        ax.axhline(y=0, color='black', linewidth=0.5)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "04_mc_convergence.png"), dpi=150)
        plt.close()
        print("  Saved 04_mc_convergence.png")

    # Plot 3: Land-constrained comparison
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    cats = ['monumental', 'settlement', 'monumental_ancient', 'settlement_ancient']
    cat_labels = ['Monument\n(all)', 'Settlement\n(all)', 'Monument\n(ancient)', 'Settlement\n(ancient)']
    x = np.arange(len(cats))
    width = 0.35

    orig_vals = []
    corr_vals = []
    for cat in cats:
        if cat in a1_results:
            orig_vals.append(a1_results[cat]['unconstrained_1k']['z_score'])
            corr_vals.append(a1_results[cat]['land_constrained_10k']['z_score'])
        else:
            orig_vals.append(0)
            corr_vals.append(0)

    bars1 = ax.bar(x - width/2, orig_vals, width, label='Unconstrained', color='#e74c3c', alpha=0.8)
    bars2 = ax.bar(x + width/2, corr_vals, width, label='Land-constrained', color='#2ecc71', alpha=0.8)

    for bars in [bars1, bars2]:
        for bar in bars:
            h = bar.get_height()
            if h != 0:
                ax.text(bar.get_x() + bar.get_width()/2, h + 0.15 * (1 if h >= 0 else -1),
                        f'{h:.1f}', ha='center', va='bottom' if h >= 0 else 'top', fontsize=8)

    ax.set_ylabel('Z-score')
    ax.set_title('Original vs Land-Constrained Z-scores\n(50 km corridor)')
    ax.set_xticks(x)
    ax.set_xticklabels(cat_labels)
    ax.legend()
    ax.axhline(y=0, color='black', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "01_land_constrained_comparison.png"), dpi=150)
    plt.close()
    print("  Saved 01_land_constrained_comparison.png")


# ============================================================
# MAIN
# ============================================================
def main():
    print("="*70)
    print("DIRECTIVE 16: COASTAL BIAS CORRECTION")
    print("="*70)
    t_start = time.time()

    # Load data
    all_sites = load_pleiades()
    categories = split_by_category(all_sites)
    for name, sites in categories.items():
        print(f"  {name}: {len(sites)} sites")

    # Run analyses
    a1 = analysis_1(categories)
    a2 = analysis_2(categories)
    a3 = analysis_3(categories)
    a4 = analysis_4(categories)
    a5 = analysis_5(a1, a2, a3, a4)

    # Generate plots
    print("\nGenerating plots...")
    generate_plots(a1, a2, a4)

    elapsed = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"ALL ANALYSES COMPLETE in {elapsed/60:.1f} minutes")
    print(f"Output directory: {OUT_DIR}")
    print(f"{'='*70}")

    # Print key results
    if 'monumental_ancient' in a1:
        ma = a1['monumental_ancient']
        print(f"\n*** KEY RESULTS ***")
        print(f"Monument (ancient) Z-score:")
        print(f"  Paper 1 original:        Z = 11.83  (5.05x enrichment, N=200)")
        print(f"  Unconstrained (N=1000):  Z = {ma['unconstrained_1k']['z_score']}  "
              f"({ma['unconstrained_1k']['enrichment']}x)")
        print(f"  Land-constrained (N=1k): Z = {ma['land_constrained_1k']['z_score']}  "
              f"({ma['land_constrained_1k']['enrichment']}x)")
        print(f"  Land-constrained (N=10k):Z = {ma['land_constrained_10k']['z_score']}  "
              f"({ma['land_constrained_10k']['enrichment']}x)")

    if 'monumental_ancient' in a2:
        rc = a2['monumental_ancient']['random_circles_10k']
        print(f"  Random circles (N=10k):  Z = {rc['z_score']}  ({rc['enrichment']}x)")

    if 'settlement_ancient' in a1:
        sa = a1['settlement_ancient']
        print(f"\nSettlement (ancient) Z-score:")
        print(f"  Land-constrained (N=10k):Z = {sa['land_constrained_10k']['z_score']}  "
              f"({sa['land_constrained_10k']['enrichment']}x)")

    print(f"\nFallback rate: {a1.get('monumental_ancient', {}).get('fallback_rate', 0)*100:.2f}%")


if __name__ == "__main__":
    main()
