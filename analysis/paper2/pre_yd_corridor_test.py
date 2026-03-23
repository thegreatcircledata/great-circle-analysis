#!/usr/bin/env python3
"""
Pre-Younger Dryas Corridor Test
=================================
Tests whether the Great Circle corridor shows anomalous human activity
before the Younger Dryas (>12,800 BP) compared to matched corridors.

Seven tests building in complexity:
  1. Temporal density profile (on-corridor vs global)
  2. Pre-YD corridor enrichment vs matched random corridors (THE key test)
  3. YD disruption profile (corridor vs global)
  4. Regional decomposition
  5. Continuity vs discontinuity across YD
  6. Atlantis region test
  7. Summed probability distribution (simplified uncalibrated version)

Uses three databases: p3k14c, XRONOS, ROAD v32
Output: outputs/pre_yd_corridor_test/

Author: Claude (for Ell)
Date: 2026-03-20
"""

import csv, math, json, os, sys, time, warnings
import numpy as np
from collections import defaultdict

warnings.filterwarnings('ignore')
np.random.seed(42)

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087

# Corridor thresholds
PRIMARY_THRESHOLD_KM = 50
ALT_THRESHOLDS_KM = [25, 100]

# Monte Carlo
N_RANDOM_CIRCLES = 1000

# Temporal bins
BIN_WIDTH = 500
BIN_START = 5000
BIN_END = 25000
BIN_EDGES = list(range(BIN_START, BIN_END + BIN_WIDTH, BIN_WIDTH))
BINS = [(BIN_EDGES[i], BIN_EDGES[i+1]) for i in range(len(BIN_EDGES)-1)]

# Key time periods
PRE_YD_CUTOFF = 12800  # BP
YD_START = 12800
YD_END = 11600
POST_YD_START = 11600
EARLY_HOLOCENE_END = 8000

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "pre_yd_corridor_test")
FIG_DIR = os.path.join(OUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ============================================================
# VECTORIZED HELPERS
# ============================================================
def haversine_vec(lat1, lon1, lats, lons):
    """Haversine distance from single point to array of points (km)."""
    lat1_r, lon1_r = np.radians(lat1), np.radians(lon1)
    lat2_r = np.radians(lats)
    lon2_r = np.radians(lons)
    dlat = lat2_r - lat1_r
    dlon = lon2_r - lon1_r
    a = np.sin(dlat/2)**2 + np.cos(lat1_r)*np.cos(lat2_r)*np.sin(dlon/2)**2
    return EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))


def gc_dist_from_pole(pole_lat, pole_lon, lats, lons):
    """Distance from each point to the great circle defined by pole (km)."""
    d = haversine_vec(pole_lat, pole_lon, lats, lons)
    return np.abs(d - QUARTER_CIRC)


def random_pole():
    """Generate a uniformly random point on the sphere."""
    z = np.random.uniform(-1, 1)
    theta = np.random.uniform(0, 2 * np.pi)
    lat = np.degrees(np.arcsin(z))
    lon = np.degrees(theta) - 180
    return lat, lon


def classify_region(lat, lon):
    """Classify a coordinate into a Great Circle regional segment."""
    if 25 <= lat <= 35 and 25 <= lon <= 40:
        return "Egypt/Levant"
    elif 25 <= lat <= 40 and 40 <= lon <= 65:
        return "Iran/Central Asia"
    elif 5 <= lat <= 35 and 65 <= lon <= 95:
        return "South Asia"
    elif -15 <= lat <= 25 and 95 <= lon <= 140:
        return "SE Asia"
    elif -50 <= lat <= 25 and 140 <= lon <= 180:
        return "Pacific/Oceania"
    elif -20 <= lat <= 5 and -85 <= lon <= -65:
        return "Peru/Amazon"
    elif 15 <= lat <= 35 and -20 <= lon <= 25:
        return "North Africa/Mediterranean"
    elif -10 <= lat <= 15 and -20 <= lon <= 45:
        return "Sub-Saharan Africa"
    else:
        return "Other"


# ============================================================
# DATA LOADING
# ============================================================
def load_p3k14c():
    """Load p3k14c database."""
    filepath = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")
    dates = []
    with open(filepath, encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                age = int(float(row['Age']))
                lat = float(row['Lat'])
                lon = float(row['Long'])
                if 0 < age <= 50000 and -90 <= lat <= 90 and -180 <= lon <= 180:
                    if lat == 0 and lon == 0:
                        continue
                    dates.append({
                        'age': age, 'lat': lat, 'lon': lon,
                        'site_id': row.get('SiteID', ''),
                        'source': 'p3k14c'
                    })
            except (ValueError, KeyError, TypeError):
                continue
    return dates


def load_xronos():
    """Load XRONOS raw JSON database."""
    filepath = os.path.join(BASE_DIR, "data", "xronos", "xronos_raw.json")
    import json as _json
    with open(filepath, 'r') as f:
        data = _json.load(f)
    dates = []
    for rec in data:
        m = rec.get('measurement', {})
        try:
            bp = int(m.get('bp', 0))
            lat = float(m.get('lat', 0))
            lng = float(m.get('lng', 0))
            if 0 < bp <= 50000 and -90 <= lat <= 90 and -180 <= lng <= 180:
                if lat == 0 and lng == 0:
                    continue
                dates.append({
                    'age': bp, 'lat': lat, 'lon': lng,
                    'site_id': m.get('site', ''),
                    'source': 'xronos'
                })
        except (ValueError, TypeError):
            continue
    return dates


def load_road():
    """Load ROAD v32 database (All sheet)."""
    try:
        import openpyxl
    except ImportError:
        print("  WARNING: openpyxl not available, skipping ROAD database")
        return []
    filepath = os.path.join(BASE_DIR, "data", "road", "road_v32.xlsx")
    wb = openpyxl.load_workbook(filepath, read_only=True)
    ws = wb['All']
    dates = []
    for row in ws.iter_rows(min_row=2, values_only=True):
        try:
            age = int(row[8]) if row[8] else None
            lat = float(row[5]) if row[5] else None
            lon = float(row[4]) if row[4] else None
            if age and lat and lon and 0 < age <= 60000:
                if lat == 0 and lon == 0:
                    continue
                if -90 <= lat <= 90 and -180 <= lon <= 180:
                    dates.append({
                        'age': age, 'lat': lat, 'lon': lon,
                        'site_id': str(row[1] or ''),
                        'source': 'road'
                    })
        except (ValueError, TypeError):
            continue
    wb.close()
    return dates


def merge_and_deduplicate(db_list):
    """Merge databases, deduplicating by location (~0.01°) and age (~200 yr)."""
    all_dates = []
    for db in db_list:
        all_dates.extend(db)
    print(f"  Total before dedup: {len(all_dates)}")

    # Sort by rounded lat, lon, age
    all_dates.sort(key=lambda d: (round(d['lat'], 2), round(d['lon'], 2), d['age']))

    deduped = [all_dates[0]]
    for d in all_dates[1:]:
        prev = deduped[-1]
        if (abs(d['lat'] - prev['lat']) < 0.01 and
            abs(d['lon'] - prev['lon']) < 0.01 and
            abs(d['age'] - prev['age']) < 200):
            continue
        deduped.append(d)

    print(f"  After dedup: {len(deduped)}")
    return deduped


# ============================================================
# MAIN ANALYSIS
# ============================================================
def main():
    t0 = time.time()
    print("=" * 70)
    print("PRE-YOUNGER DRYAS CORRIDOR TEST")
    print("=" * 70)
    print(f"Date: 2026-03-20")
    print(f"Pole: {POLE_LAT}°N, {POLE_LON}°E")
    print(f"Primary threshold: {PRIMARY_THRESHOLD_KM} km")
    print(f"Random circles: {N_RANDOM_CIRCLES}")
    print(f"Temporal bins: {BIN_WIDTH}-year from {BIN_START} to {BIN_END} BP")

    # ---- Load databases ----
    print("\n" + "=" * 70)
    print("LOADING DATABASES")
    print("=" * 70)

    print("\n1. Loading p3k14c...")
    p3k = load_p3k14c()
    print(f"   Loaded: {len(p3k)} dates")

    print("\n2. Loading XRONOS (this takes a moment)...")
    xronos = load_xronos()
    print(f"   Loaded: {len(xronos)} dates")

    print("\n3. Loading ROAD v32...")
    road = load_road()
    print(f"   Loaded: {len(road)} dates")

    print("\n4. Merging and deduplicating...")
    all_dates = merge_and_deduplicate([p3k, xronos, road])

    # Convert to numpy arrays
    lats = np.array([d['lat'] for d in all_dates])
    lons = np.array([d['lon'] for d in all_dates])
    ages = np.array([d['age'] for d in all_dates])
    site_ids = np.array([d['site_id'] for d in all_dates])
    sources = np.array([d['source'] for d in all_dates])

    # Data summary
    pre_yd_mask = ages > PRE_YD_CUTOFF
    print(f"\n   Total merged dates: {len(all_dates)}")
    print(f"   Pre-YD (>{PRE_YD_CUTOFF} BP): {np.sum(pre_yd_mask)}")
    print(f"   Post-YD: {np.sum(~pre_yd_mask)}")
    print(f"   Source breakdown:")
    for src in ['p3k14c', 'xronos', 'road']:
        n = np.sum(sources == src)
        n_pre = np.sum((sources == src) & pre_yd_mask)
        print(f"     {src}: {n} total, {n_pre} pre-YD")

    # Compute distance to Great Circle for all dates
    print("\nComputing distances to Great Circle...")
    gc_dists = gc_dist_from_pole(POLE_LAT, POLE_LON, lats, lons)
    on_corridor = gc_dists <= PRIMARY_THRESHOLD_KM
    print(f"  On-corridor (≤{PRIMARY_THRESHOLD_KM}km): {np.sum(on_corridor)}")
    print(f"  Off-corridor: {np.sum(~on_corridor)}")
    print(f"  On-corridor pre-YD: {np.sum(on_corridor & pre_yd_mask)}")

    # ============================================================
    # TEST 1: TEMPORAL DENSITY PROFILE
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 1: TEMPORAL DENSITY PROFILE — On-Corridor vs Global")
    print("=" * 70)

    test1_results = {
        'test': 'Temporal Density Profile',
        'description': 'Enrichment of on-corridor dates relative to global, by 500-year bin',
        'threshold_km': PRIMARY_THRESHOLD_KM,
        'bins': []
    }

    total_on = np.sum(on_corridor)
    total_off = np.sum(~on_corridor)

    print(f"\n{'Bin (BP)':>16s}  {'On':>6s}  {'Off':>8s}  {'Enrichment':>11s}  {'Period':>20s}")
    print("-" * 70)

    for low_bp, high_bp in BINS:
        age_mask = (ages >= low_bp) & (ages < high_bp)
        on_count = int(np.sum(age_mask & on_corridor))
        off_count = int(np.sum(age_mask & ~on_corridor))

        # Enrichment ratio (normalized)
        if off_count > 0 and total_on > 0 and total_off > 0:
            enrichment = (on_count / total_on) / (off_count / total_off)
        else:
            enrichment = 0.0

        # Label key periods
        mid = (low_bp + high_bp) / 2
        if 12800 <= mid <= 13300:
            period = "** YD ONSET **"
        elif 11600 <= mid <= 12800:
            period = "YD duration"
        elif 10000 <= mid <= 11600:
            period = "Post-YD recovery"
        elif 8000 <= mid <= 10000:
            period = "Early Holocene peak"
        elif mid > 13300:
            period = "Pre-YD"
        else:
            period = ""

        test1_results['bins'].append({
            'bin_start': low_bp,
            'bin_end': high_bp,
            'on_corridor_count': on_count,
            'off_corridor_count': off_count,
            'enrichment': round(enrichment, 4),
            'period_label': period
        })

        enr_str = f"{enrichment:.3f}" if enrichment > 0 else "N/A"
        print(f"  {low_bp:>5d}-{high_bp:>5d} BP  {on_count:>6d}  {off_count:>8d}  {enr_str:>11s}  {period:>20s}")

    # Key period summaries
    print("\n  KEY PERIOD SUMMARIES:")
    for period_name, age_lo, age_hi in [
        ("Deep pre-YD (25000-15000)", 15000, 25000),
        ("Late Glacial (15000-12800)", 12800, 15000),
        ("YD onset (12800-12300)", 12300, 12800),
        ("YD duration (12300-11600)", 11600, 12300),
        ("Post-YD recovery (11600-10000)", 10000, 11600),
        ("Early Holocene peak (10000-8000)", 8000, 10000),
        ("Mid Holocene (8000-5000)", 5000, 8000),
    ]:
        mask = (ages >= age_lo) & (ages < age_hi)
        on_n = int(np.sum(mask & on_corridor))
        off_n = int(np.sum(mask & ~on_corridor))
        if off_n > 0 and total_on > 0 and total_off > 0:
            enr = (on_n / total_on) / (off_n / total_off)
        else:
            enr = 0
        print(f"    {period_name}: on={on_n}, off={off_n}, enrichment={enr:.3f}")

    # Save Test 1
    with open(os.path.join(OUT_DIR, "test1_temporal_profile.json"), 'w') as f:
        json.dump(test1_results, f, indent=2)
    print(f"\n  Saved: test1_temporal_profile.json")

    # ============================================================
    # TEST 2: PRE-YD CORRIDOR ENRICHMENT VS MATCHED CORRIDORS
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 2: PRE-YD CORRIDOR ENRICHMENT VS MATCHED CORRIDORS")
    print("  *** THIS IS THE KEY TEST ***")
    print("=" * 70)

    # Observed counts
    obs_pre_yd_on = int(np.sum(pre_yd_mask & on_corridor))
    obs_total_on = int(np.sum(on_corridor))

    # Sub-periods
    deep_pre_yd = ages > 15000  # 25000-15000 BP
    late_glacial = (ages > 12800) & (ages <= 15000)

    obs_deep_on = int(np.sum(deep_pre_yd & on_corridor))
    obs_late_glacial_on = int(np.sum(late_glacial & on_corridor))

    print(f"\n  Observed (Alison Great Circle, {PRIMARY_THRESHOLD_KM}km):")
    print(f"    Total on-corridor: {obs_total_on}")
    print(f"    Pre-YD on-corridor (>12800 BP): {obs_pre_yd_on}")
    print(f"    Deep pre-YD on-corridor (>15000 BP): {obs_deep_on}")
    print(f"    Late Glacial on-corridor (12800-15000 BP): {obs_late_glacial_on}")

    # Monte Carlo: random great circles
    print(f"\n  Running {N_RANDOM_CIRCLES} random great circles...")
    random_pre_yd = []
    random_deep = []
    random_late_glacial = []
    random_total = []

    # Pre-compute pre-YD subset for speed
    pre_yd_lats = lats[pre_yd_mask]
    pre_yd_lons = lons[pre_yd_mask]
    pre_yd_ages = ages[pre_yd_mask]

    deep_mask_sub = pre_yd_ages > 15000
    lg_mask_sub = pre_yd_ages <= 15000  # 12800-15000 (already filtered to >12800)

    for i in range(N_RANDOM_CIRCLES):
        if (i+1) % 100 == 0:
            print(f"    Circle {i+1}/{N_RANDOM_CIRCLES}...")

        plat, plon = random_pole()

        # Total dates on this random circle
        rc_dists = gc_dist_from_pole(plat, plon, lats, lons)
        rc_on = rc_dists <= PRIMARY_THRESHOLD_KM
        random_total.append(int(np.sum(rc_on)))

        # Pre-YD dates on this random circle
        rc_dists_pre = gc_dist_from_pole(plat, plon, pre_yd_lats, pre_yd_lons)
        rc_pre_on = rc_dists_pre <= PRIMARY_THRESHOLD_KM
        random_pre_yd.append(int(np.sum(rc_pre_on)))

        # Deep pre-YD
        random_deep.append(int(np.sum(rc_pre_on & deep_mask_sub)))
        random_late_glacial.append(int(np.sum(rc_pre_on & lg_mask_sub)))

    random_pre_yd = np.array(random_pre_yd)
    random_deep = np.array(random_deep)
    random_late_glacial = np.array(random_late_glacial)
    random_total = np.array(random_total)

    # Z-scores
    def zscore(obs, dist):
        mu = np.mean(dist)
        sigma = np.std(dist)
        if sigma == 0:
            return 0.0, mu, sigma
        return (obs - mu) / sigma, mu, sigma

    z_pre_yd, mu_pre, std_pre = zscore(obs_pre_yd_on, random_pre_yd)
    z_deep, mu_deep, std_deep = zscore(obs_deep_on, random_deep)
    z_lg, mu_lg, std_lg = zscore(obs_late_glacial_on, random_late_glacial)
    z_total, mu_total, std_total = zscore(obs_total_on, random_total)

    # Percentile rank
    pctile_pre = np.mean(random_pre_yd < obs_pre_yd_on) * 100
    pctile_deep = np.mean(random_deep < obs_deep_on) * 100
    pctile_lg = np.mean(random_late_glacial < obs_late_glacial_on) * 100

    # Habitability-matched comparison: circles with ±20% of Alison's total count
    hab_lo = obs_total_on * 0.8
    hab_hi = obs_total_on * 1.2
    hab_matched = (random_total >= hab_lo) & (random_total <= hab_hi)
    n_hab_matched = int(np.sum(hab_matched))
    if n_hab_matched > 10:
        hab_pre_yd = random_pre_yd[hab_matched]
        z_hab, mu_hab, std_hab = zscore(obs_pre_yd_on, hab_pre_yd)
        pctile_hab = np.mean(hab_pre_yd < obs_pre_yd_on) * 100
    else:
        z_hab, mu_hab, std_hab = 0.0, 0.0, 0.0
        pctile_hab = 0.0

    print(f"\n  RESULTS:")
    print(f"    {'Metric':>35s}  {'Observed':>8s}  {'Mean(rand)':>10s}  {'Std':>8s}  {'Z-score':>8s}  {'%ile':>6s}")
    print(f"    {'-'*80}")
    print(f"    {'Total on-corridor':>35s}  {obs_total_on:>8d}  {mu_total:>10.1f}  {std_total:>8.1f}  {z_total:>8.2f}  {'':>6s}")
    print(f"    {'Pre-YD (>12800 BP)':>35s}  {obs_pre_yd_on:>8d}  {mu_pre:>10.1f}  {std_pre:>8.1f}  {z_pre_yd:>8.2f}  {pctile_pre:>5.1f}%")
    print(f"    {'Deep pre-YD (>15000 BP)':>35s}  {obs_deep_on:>8d}  {mu_deep:>10.1f}  {std_deep:>8.1f}  {z_deep:>8.2f}  {pctile_deep:>5.1f}%")
    print(f"    {'Late Glacial (12800-15000 BP)':>35s}  {obs_late_glacial_on:>8d}  {mu_lg:>10.1f}  {std_lg:>8.1f}  {z_lg:>8.2f}  {pctile_lg:>5.1f}%")
    if n_hab_matched > 10:
        print(f"    {'Hab-matched pre-YD':>35s}  {obs_pre_yd_on:>8d}  {mu_hab:>10.1f}  {std_hab:>8.1f}  {z_hab:>8.2f}  {pctile_hab:>5.1f}%")
        print(f"    (Habitability-matched: {n_hab_matched} circles with total count ±20% of Alison)")

    print(f"\n  INTERPRETATION:")
    if z_pre_yd < 2:
        print(f"    Z = {z_pre_yd:.2f} → NO pre-YD anomaly. Corridor is unremarkable before YD.")
        print(f"    Supports conventional narrative.")
    elif z_pre_yd < 3:
        print(f"    Z = {z_pre_yd:.2f} → MILD enrichment. May be habitability-driven.")
    elif z_pre_yd < 5:
        print(f"    Z = {z_pre_yd:.2f} → SIGNIFICANT pre-YD anomaly. Needs scrutiny.")
    else:
        print(f"    Z = {z_pre_yd:.2f} → EXTRAORDINARY. Triple-check everything.")

    test2_results = {
        'test': 'Pre-YD Corridor Enrichment vs Matched Corridors',
        'threshold_km': PRIMARY_THRESHOLD_KM,
        'n_random_circles': N_RANDOM_CIRCLES,
        'observed': {
            'total_on_corridor': obs_total_on,
            'pre_yd_on_corridor': obs_pre_yd_on,
            'deep_pre_yd_on_corridor': obs_deep_on,
            'late_glacial_on_corridor': obs_late_glacial_on,
        },
        'random_distribution': {
            'pre_yd': {'mean': round(mu_pre, 2), 'std': round(std_pre, 2),
                       'z_score': round(z_pre_yd, 4), 'percentile': round(pctile_pre, 2),
                       'min': int(np.min(random_pre_yd)), 'max': int(np.max(random_pre_yd)),
                       'median': int(np.median(random_pre_yd))},
            'deep_pre_yd': {'mean': round(mu_deep, 2), 'std': round(std_deep, 2),
                            'z_score': round(z_deep, 4), 'percentile': round(pctile_deep, 2)},
            'late_glacial': {'mean': round(mu_lg, 2), 'std': round(std_lg, 2),
                             'z_score': round(z_lg, 4), 'percentile': round(pctile_lg, 2)},
            'total': {'mean': round(mu_total, 2), 'std': round(std_total, 2),
                      'z_score': round(z_total, 4)},
        },
        'habitability_matched': {
            'n_matched_circles': n_hab_matched,
            'z_score': round(z_hab, 4) if n_hab_matched > 10 else None,
            'percentile': round(pctile_hab, 2) if n_hab_matched > 10 else None,
        },
        'random_pre_yd_counts': random_pre_yd.tolist(),
    }

    with open(os.path.join(OUT_DIR, "test2_pre_yd_enrichment.json"), 'w') as f:
        json.dump(test2_results, f, indent=2)
    print(f"\n  Saved: test2_pre_yd_enrichment.json")

    # ============================================================
    # TEST 2b: Z-SCORE OVER TIME (extended to all bins)
    # ============================================================
    print("\n  Computing Z-scores by temporal bin (extended Test 2)...")

    # We already have random circle poles, reuse them
    # For each bin, compute on-corridor count for Alison vs distribution of randoms
    bin_zscores = []

    # Store random circle distances for reuse
    print("    Pre-computing random circle distances...")
    random_poles = []
    for i in range(N_RANDOM_CIRCLES):
        plat, plon = random_pole()
        random_poles.append((plat, plon))

    for low_bp, high_bp in BINS:
        age_mask = (ages >= low_bp) & (ages < high_bp)
        bin_lats = lats[age_mask]
        bin_lons = lons[age_mask]

        if len(bin_lats) < 10:
            bin_zscores.append({
                'bin_start': low_bp, 'bin_end': high_bp,
                'z_score': None, 'observed': 0, 'n_dates_bin': int(np.sum(age_mask))
            })
            continue

        obs_bin = int(np.sum(age_mask & on_corridor))

        rand_counts = []
        for plat, plon in random_poles:
            rd = gc_dist_from_pole(plat, plon, bin_lats, bin_lons)
            rand_counts.append(int(np.sum(rd <= PRIMARY_THRESHOLD_KM)))

        rand_counts = np.array(rand_counts)
        z, mu, sigma = zscore(obs_bin, rand_counts)

        bin_zscores.append({
            'bin_start': low_bp, 'bin_end': high_bp,
            'z_score': round(z, 4), 'observed': obs_bin,
            'mean_random': round(mu, 2), 'std_random': round(sigma, 2),
            'n_dates_bin': int(np.sum(age_mask))
        })

    test2_results['bin_zscores'] = bin_zscores

    # Re-save with bin z-scores
    with open(os.path.join(OUT_DIR, "test2_pre_yd_enrichment.json"), 'w') as f:
        json.dump(test2_results, f, indent=2)

    print(f"\n    {'Bin (BP)':>16s}  {'Obs':>5s}  {'Mean':>7s}  {'Z':>7s}")
    print(f"    {'-'*40}")
    for bz in bin_zscores:
        z_str = f"{bz['z_score']:.2f}" if bz['z_score'] is not None else "N/A"
        m_str = f"{bz.get('mean_random', 0):.1f}" if bz.get('mean_random') else "N/A"
        print(f"    {bz['bin_start']:>5d}-{bz['bin_end']:>5d} BP  {bz['observed']:>5d}  {m_str:>7s}  {z_str:>7s}")

    # ============================================================
    # TEST 3: YD DISRUPTION PROFILE
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 3: YD DISRUPTION PROFILE — Corridor vs Global")
    print("=" * 70)

    # Dates in 500 years before YD onset (13300-12800 BP) vs 500 years after (12800-12300 BP)
    pre_yd_window = (ages >= 12800) & (ages < 13300)
    post_yd_onset = (ages >= 12300) & (ages < 12800)

    pre_on = int(np.sum(pre_yd_window & on_corridor))
    post_on = int(np.sum(post_yd_onset & on_corridor))
    pre_off = int(np.sum(pre_yd_window & ~on_corridor))
    post_off = int(np.sum(post_yd_onset & ~on_corridor))

    crash_corridor = post_on / pre_on if pre_on > 0 else 0
    crash_global = post_off / pre_off if pre_off > 0 else 0

    print(f"\n  On-corridor: pre-YD window = {pre_on}, post-YD onset = {post_on}")
    print(f"  Off-corridor: pre-YD window = {pre_off}, post-YD onset = {post_off}")
    print(f"  Crash ratio (corridor): {crash_corridor:.3f}")
    print(f"  Crash ratio (global): {crash_global:.3f}")

    if crash_corridor < crash_global:
        print(f"  → Corridor crashed HARDER than global ({crash_corridor:.3f} vs {crash_global:.3f})")
    elif crash_corridor > crash_global:
        print(f"  → Corridor was MORE RESILIENT than global ({crash_corridor:.3f} vs {crash_global:.3f})")
    else:
        print(f"  → Corridor matches global pattern")

    # Compare to random corridors
    random_crash = []
    for plat, plon in random_poles:
        rd = gc_dist_from_pole(plat, plon, lats, lons)
        rc_on = rd <= PRIMARY_THRESHOLD_KM

        rpre = int(np.sum(pre_yd_window & rc_on))
        rpost = int(np.sum(post_yd_onset & rc_on))

        if rpre > 0:
            random_crash.append(rpost / rpre)

    random_crash = np.array(random_crash)
    if len(random_crash) > 10:
        z_crash, mu_crash, std_crash = zscore(crash_corridor, random_crash)
        print(f"\n  Random corridor crash ratios: mean={mu_crash:.3f}, std={std_crash:.3f}")
        print(f"  Alison crash Z-score: {z_crash:.2f}")
    else:
        z_crash, mu_crash, std_crash = 0, 0, 0
        print(f"\n  Insufficient random corridors with pre-YD dates for comparison")

    # Also check recovery: post-YD (11600-10000) vs pre-YD (13300-12800)
    recovery_window = (ages >= 10000) & (ages < 11600)
    rec_on = int(np.sum(recovery_window & on_corridor))
    rec_off = int(np.sum(recovery_window & ~on_corridor))
    recovery_corridor = rec_on / pre_on if pre_on > 0 else 0
    recovery_global = rec_off / pre_off if pre_off > 0 else 0

    print(f"\n  Recovery ratio (post-YD 11600-10000 / pre-YD 13300-12800):")
    print(f"    Corridor: {recovery_corridor:.3f}")
    print(f"    Global: {recovery_global:.3f}")

    test3_results = {
        'test': 'YD Disruption Profile',
        'corridor': {
            'pre_yd_window_count': pre_on,
            'post_yd_onset_count': post_on,
            'crash_ratio': round(crash_corridor, 4),
            'recovery_count': rec_on,
            'recovery_ratio': round(recovery_corridor, 4),
        },
        'global': {
            'pre_yd_window_count': pre_off,
            'post_yd_onset_count': post_off,
            'crash_ratio': round(crash_global, 4),
            'recovery_count': rec_off,
            'recovery_ratio': round(recovery_global, 4),
        },
        'z_score_crash': round(z_crash, 4) if len(random_crash) > 10 else None,
        'random_crash_mean': round(mu_crash, 4) if len(random_crash) > 10 else None,
        'random_crash_std': round(std_crash, 4) if len(random_crash) > 10 else None,
        'n_random_with_data': len(random_crash),
    }

    with open(os.path.join(OUT_DIR, "test3_yd_disruption.json"), 'w') as f:
        json.dump(test3_results, f, indent=2)
    print(f"\n  Saved: test3_yd_disruption.json")

    # ============================================================
    # TEST 4: REGIONAL DECOMPOSITION
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 4: REGIONAL DECOMPOSITION")
    print("=" * 70)

    # Classify all on-corridor dates into regions
    on_indices = np.where(on_corridor)[0]
    region_map = defaultdict(lambda: {'pre_yd': 0, 'post_yd': 0, 'total': 0})

    for idx in on_indices:
        region = classify_region(lats[idx], lons[idx])
        age = ages[idx]
        region_map[region]['total'] += 1
        if age > PRE_YD_CUTOFF:
            region_map[region]['pre_yd'] += 1
        else:
            region_map[region]['post_yd'] += 1

    print(f"\n  {'Region':>30s}  {'Pre-YD':>7s}  {'Post-YD':>8s}  {'Total':>7s}  {'Pre-YD%':>8s}")
    print(f"  {'-'*65}")
    for region in sorted(region_map.keys(), key=lambda r: region_map[r]['pre_yd'], reverse=True):
        d = region_map[region]
        pct = d['pre_yd'] / d['total'] * 100 if d['total'] > 0 else 0
        print(f"  {region:>30s}  {d['pre_yd']:>7d}  {d['post_yd']:>8d}  {d['total']:>7d}  {pct:>7.1f}%")

    # Leave-one-out: remove each region and recompute pre-YD Z-score
    print(f"\n  Leave-one-out robustness (removing each region):")
    loo_results = {}

    for exclude_region in sorted(region_map.keys()):
        if region_map[exclude_region]['pre_yd'] == 0:
            continue

        # Identify which on-corridor dates to exclude
        exclude_mask = np.zeros(len(lats), dtype=bool)
        for idx in on_indices:
            if classify_region(lats[idx], lons[idx]) == exclude_region:
                exclude_mask[idx] = True

        # Recount pre-YD on-corridor excluding this region
        obs_excl = int(np.sum(pre_yd_mask & on_corridor & ~exclude_mask))

        # For random circles, we also need to exclude dates in that geographic region
        # Simpler: just report the observed count change
        removed = region_map[exclude_region]['pre_yd']
        pct_of_total = removed / obs_pre_yd_on * 100 if obs_pre_yd_on > 0 else 0

        loo_results[exclude_region] = {
            'removed_pre_yd': removed,
            'remaining_pre_yd': obs_excl,
            'pct_of_total_pre_yd': round(pct_of_total, 1)
        }
        print(f"    Exclude {exclude_region}: removed {removed} pre-YD dates ({pct_of_total:.1f}% of total), {obs_excl} remaining")

    test4_results = {
        'test': 'Regional Decomposition',
        'regions': {r: dict(v) for r, v in region_map.items()},
        'leave_one_out': loo_results,
    }

    with open(os.path.join(OUT_DIR, "test4_regional_decomposition.json"), 'w') as f:
        json.dump(test4_results, f, indent=2)
    print(f"\n  Saved: test4_regional_decomposition.json")

    # ============================================================
    # TEST 5: CONTINUITY VS DISCONTINUITY ACROSS YD
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 5: CONTINUITY VS DISCONTINUITY ACROSS YD")
    print("=" * 70)

    # Cluster sites by coordinates (within 5km)
    # Use a simple approach: round to ~0.05° (roughly 5km)
    CLUSTER_RES = 0.05  # degrees (~5km)

    site_clusters = defaultdict(lambda: {'ages': [], 'on_corridor': False})
    for i in range(len(lats)):
        key = (round(lats[i] / CLUSTER_RES) * CLUSTER_RES,
               round(lons[i] / CLUSTER_RES) * CLUSTER_RES)
        site_clusters[key]['ages'].append(ages[i])
        if on_corridor[i]:
            site_clusters[key]['on_corridor'] = True

    # Find continuity sites: have dates both >12800 and <11600
    continuity_on = 0
    continuity_off = 0
    total_sites_on = 0
    total_sites_off = 0

    for key, info in site_clusters.items():
        age_arr = np.array(info['ages'])
        has_pre = np.any(age_arr > YD_START)
        has_post = np.any(age_arr < YD_END)

        if info['on_corridor']:
            total_sites_on += 1
            if has_pre and has_post:
                continuity_on += 1
        else:
            total_sites_off += 1
            if has_pre and has_post:
                continuity_off += 1

    cont_rate_on = continuity_on / total_sites_on if total_sites_on > 0 else 0
    cont_rate_off = continuity_off / total_sites_off if total_sites_off > 0 else 0

    print(f"\n  Site clusters (5km resolution):")
    print(f"    On-corridor: {total_sites_on} clusters, {continuity_on} span the YD ({cont_rate_on*100:.2f}%)")
    print(f"    Off-corridor: {total_sites_off} clusters, {continuity_off} span the YD ({cont_rate_off*100:.2f}%)")

    if cont_rate_on > 0 and cont_rate_off > 0:
        ratio = cont_rate_on / cont_rate_off
        print(f"    Continuity ratio (on/off): {ratio:.2f}×")

    # Compare to random corridors
    random_continuity_rates = []
    for plat, plon in random_poles[:200]:  # Use 200 for speed
        rd = gc_dist_from_pole(plat, plon, lats, lons)
        rc_on = rd <= PRIMARY_THRESHOLD_KM

        # Quick cluster-based continuity
        rc_clusters = defaultdict(lambda: {'ages': [], 'on': False})
        for i in np.where(rc_on)[0]:
            key = (round(lats[i] / CLUSTER_RES) * CLUSTER_RES,
                   round(lons[i] / CLUSTER_RES) * CLUSTER_RES)
            rc_clusters[key]['ages'].append(ages[i])
            rc_clusters[key]['on'] = True

        n_total = len(rc_clusters)
        n_cont = 0
        for k, info in rc_clusters.items():
            aa = np.array(info['ages'])
            if np.any(aa > YD_START) and np.any(aa < YD_END):
                n_cont += 1

        if n_total > 0:
            random_continuity_rates.append(n_cont / n_total)

    random_continuity_rates = np.array(random_continuity_rates)
    if len(random_continuity_rates) > 10:
        z_cont, mu_cont, std_cont = zscore(cont_rate_on, random_continuity_rates)
        print(f"\n  Random corridor continuity: mean={mu_cont*100:.2f}%, std={std_cont*100:.2f}%")
        print(f"  Alison continuity Z-score: {z_cont:.2f}")
    else:
        z_cont = 0

    test5_results = {
        'test': 'Continuity vs Discontinuity Across YD',
        'on_corridor': {
            'total_site_clusters': total_sites_on,
            'continuity_sites': continuity_on,
            'continuity_rate': round(cont_rate_on, 6),
        },
        'off_corridor': {
            'total_site_clusters': total_sites_off,
            'continuity_sites': continuity_off,
            'continuity_rate': round(cont_rate_off, 6),
        },
        'continuity_ratio_on_off': round(cont_rate_on / cont_rate_off, 4) if cont_rate_off > 0 else None,
        'z_score_continuity': round(z_cont, 4) if len(random_continuity_rates) > 10 else None,
    }

    with open(os.path.join(OUT_DIR, "test5_continuity.json"), 'w') as f:
        json.dump(test5_results, f, indent=2)
    print(f"\n  Saved: test5_continuity.json")

    # ============================================================
    # TEST 6: ATLANTIS REGION TEST
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 6: ATLANTIS REGION TEST")
    print("=" * 70)

    atlantis_regions = [
        {"name": "Atlantic west of Gibraltar", "lat": 35.0, "lon": -10.0, "radius_km": 500,
         "rationale": "Plato's description"},
        {"name": "Richat Structure, Mauritania", "lat": 21.1, "lon": -11.4, "radius_km": 200,
         "rationale": "Alt history community claim"},
        {"name": "Persian Gulf basin", "lat": 27.0, "lon": 51.0, "radius_km": 300,
         "rationale": "On Great Circle, dry land during LGM"},
        {"name": "Sunda Shelf, SE Asia", "lat": -5.0, "lon": 110.0, "radius_km": 500,
         "rationale": "On Great Circle, dry land during LGM"},
        {"name": "Doggerland, North Sea", "lat": 55.0, "lon": 3.0, "radius_km": 300,
         "rationale": "Known submerged Mesolithic landscape"},
    ]

    # Global pre-YD/post-YD ratio as baseline
    global_pre = int(np.sum(pre_yd_mask))
    global_post = int(np.sum(~pre_yd_mask))
    global_ratio = global_pre / global_post if global_post > 0 else 0

    print(f"\n  Global pre-YD/post-YD ratio: {global_ratio:.4f} ({global_pre}/{global_post})")

    region_results = []
    for reg in atlantis_regions:
        # Distance from region center to all dates
        dists = haversine_vec(reg['lat'], reg['lon'], lats, lons)
        in_region = dists <= reg['radius_km']

        n_total = int(np.sum(in_region))
        n_pre = int(np.sum(in_region & pre_yd_mask))
        n_post = int(np.sum(in_region & ~pre_yd_mask))

        if n_post > 0:
            region_ratio = n_pre / n_post
            ratio_vs_global = region_ratio / global_ratio if global_ratio > 0 else 0
        else:
            region_ratio = float('inf') if n_pre > 0 else 0
            ratio_vs_global = 0

        # Temporal breakdown
        temporal = {}
        for label, lo, hi in [("25000-15000", 15000, 25000), ("15000-12800", 12800, 15000),
                               ("12800-11600", 11600, 12800), ("11600-8000", 8000, 11600),
                               ("8000-5000", 5000, 8000)]:
            mask = in_region & (ages >= lo) & (ages < hi)
            temporal[label] = int(np.sum(mask))

        result = {
            'name': reg['name'],
            'center_lat': reg['lat'],
            'center_lon': reg['lon'],
            'radius_km': reg['radius_km'],
            'rationale': reg['rationale'],
            'total_dates': n_total,
            'pre_yd_dates': n_pre,
            'post_yd_dates': n_post,
            'pre_post_ratio': round(region_ratio, 4) if region_ratio != float('inf') else 'inf',
            'ratio_vs_global': round(ratio_vs_global, 4),
            'temporal_breakdown': temporal,
        }
        region_results.append(result)

        anomaly = "ANOMALOUS" if ratio_vs_global > 2 else "normal" if n_total > 0 else "NO DATA"
        print(f"\n  {reg['name']}:")
        print(f"    Total: {n_total}, Pre-YD: {n_pre}, Post-YD: {n_post}")
        print(f"    Pre/Post ratio: {region_ratio:.4f} (vs global {global_ratio:.4f}) → {ratio_vs_global:.2f}× global")
        print(f"    Assessment: {anomaly}")
        if n_total < 20:
            print(f"    ⚠ CAVEAT: Very low data coverage ({n_total} dates). Most of this region may be underwater.")

    test6_results = {
        'test': 'Atlantis Region Test',
        'global_pre_post_ratio': round(global_ratio, 6),
        'regions': region_results,
        'caveat': 'Most candidate regions are now partially or fully submerged. Low data counts reflect absence of excavation, NOT necessarily absence of activity.',
    }

    with open(os.path.join(OUT_DIR, "test6_atlantis_regions.json"), 'w') as f:
        json.dump(test6_results, f, indent=2)
    print(f"\n  Saved: test6_atlantis_regions.json")

    # ============================================================
    # TEST 7: SIMPLIFIED SPD (uncalibrated bins)
    # ============================================================
    print("\n" + "=" * 70)
    print("TEST 7: SUMMED PROBABILITY DISTRIBUTION (simplified, uncalibrated)")
    print("=" * 70)
    print("  NOTE: Using uncalibrated ages binned at 200-year intervals.")
    print("  This is a limitation — proper SPD requires IntCal20 calibration.")

    SPD_BIN_WIDTH = 200
    SPD_EDGES = list(range(5000, 25200, SPD_BIN_WIDTH))
    SPD_BINS = [(SPD_EDGES[i], SPD_EDGES[i+1]) for i in range(len(SPD_EDGES)-1)]

    # For each date, create a Gaussian kernel (age ± error)
    # Since we don't have error for all merged dates, use bin counting
    # with a smoothing window instead

    spd_on = []
    spd_off = []

    for low_bp, high_bp in SPD_BINS:
        age_mask = (ages >= low_bp) & (ages < high_bp)
        on_n = int(np.sum(age_mask & on_corridor))
        off_n = int(np.sum(age_mask & ~on_corridor))
        spd_on.append(on_n)
        spd_off.append(off_n)

    spd_on = np.array(spd_on, dtype=float)
    spd_off = np.array(spd_off, dtype=float)

    # Normalize to density
    if np.sum(spd_on) > 0:
        spd_on_norm = spd_on / np.sum(spd_on)
    else:
        spd_on_norm = spd_on
    if np.sum(spd_off) > 0:
        spd_off_norm = spd_off / np.sum(spd_off)
    else:
        spd_off_norm = spd_off

    # Compute divergence in pre-YD bins
    pre_yd_bins = [(lo, hi) for lo, hi in SPD_BINS if lo >= 12800]
    pre_yd_on_total = sum(spd_on[i] for i, (lo, hi) in enumerate(SPD_BINS) if lo >= 12800)
    pre_yd_off_total = sum(spd_off[i] for i, (lo, hi) in enumerate(SPD_BINS) if lo >= 12800)
    post_yd_on_total = sum(spd_on[i] for i, (lo, hi) in enumerate(SPD_BINS) if lo < 11600)
    post_yd_off_total = sum(spd_off[i] for i, (lo, hi) in enumerate(SPD_BINS) if lo < 11600)

    print(f"\n  SPD Summary (200-year bins):")
    print(f"    On-corridor: {int(pre_yd_on_total)} pre-YD, {int(post_yd_on_total)} post-YD bins")
    print(f"    Off-corridor: {int(pre_yd_off_total)} pre-YD, {int(post_yd_off_total)} post-YD bins")

    # Divergence metric: KL-like divergence in pre-YD period
    # Compare on-corridor SPD shape to off-corridor SPD shape
    divergence_pre_yd = []
    divergence_post_yd = []

    for i, (lo, hi) in enumerate(SPD_BINS):
        if spd_off_norm[i] > 0 and spd_on_norm[i] > 0:
            ratio = spd_on_norm[i] / spd_off_norm[i]
            if lo >= 12800:
                divergence_pre_yd.append(ratio)
            elif lo < 11600:
                divergence_post_yd.append(ratio)

    if divergence_pre_yd:
        mean_div_pre = np.mean(divergence_pre_yd)
        print(f"    Mean on/off density ratio (pre-YD bins): {mean_div_pre:.3f}")
    if divergence_post_yd:
        mean_div_post = np.mean(divergence_post_yd)
        print(f"    Mean on/off density ratio (post-YD bins): {mean_div_post:.3f}")

    test7_results = {
        'test': 'Summed Probability Distribution (simplified uncalibrated)',
        'limitation': 'Uses uncalibrated ages in 200-year bins. Proper SPD requires IntCal20 calibration.',
        'spd_bins': [{'bin_start': lo, 'bin_end': hi,
                      'on_corridor': int(spd_on[i]),
                      'off_corridor': int(spd_off[i]),
                      'on_normalized': round(float(spd_on_norm[i]), 8),
                      'off_normalized': round(float(spd_off_norm[i]), 8)}
                     for i, (lo, hi) in enumerate(SPD_BINS)],
        'pre_yd_mean_density_ratio': round(float(np.mean(divergence_pre_yd)), 4) if divergence_pre_yd else None,
        'post_yd_mean_density_ratio': round(float(np.mean(divergence_post_yd)), 4) if divergence_post_yd else None,
    }

    with open(os.path.join(OUT_DIR, "test7_spd.json"), 'w') as f:
        json.dump(test7_results, f, indent=2)
    print(f"\n  Saved: test7_spd.json")

    # ============================================================
    # SENSITIVITY: ALTERNATE THRESHOLDS
    # ============================================================
    print("\n" + "=" * 70)
    print("SENSITIVITY: ALTERNATE CORRIDOR WIDTHS")
    print("=" * 70)

    sensitivity_results = {}
    for thresh in ALT_THRESHOLDS_KM:
        alt_on = gc_dists <= thresh
        alt_pre_on = int(np.sum(pre_yd_mask & alt_on))
        alt_total_on = int(np.sum(alt_on))

        # Quick MC comparison (100 circles for speed)
        rand_pre = []
        for plat, plon in random_poles[:200]:
            rd = gc_dist_from_pole(plat, plon, pre_yd_lats, pre_yd_lons)
            rand_pre.append(int(np.sum(rd <= thresh)))
        rand_pre = np.array(rand_pre)
        z_alt, mu_alt, std_alt = zscore(alt_pre_on, rand_pre)

        print(f"\n  {thresh}km threshold:")
        print(f"    On-corridor total: {alt_total_on}, pre-YD: {alt_pre_on}")
        print(f"    Random mean pre-YD: {mu_alt:.1f} ± {std_alt:.1f}")
        print(f"    Z-score: {z_alt:.2f}")

        sensitivity_results[f"{thresh}km"] = {
            'total_on': alt_total_on,
            'pre_yd_on': alt_pre_on,
            'z_score': round(z_alt, 4),
            'random_mean': round(mu_alt, 2),
            'random_std': round(std_alt, 2),
        }

    # SENSITIVITY: Exclude North America (coordinate obfuscation in p3k14c)
    print("\n  Excluding North America (p3k14c coordinate obfuscation control):")
    na_mask = (lats >= 15) & (lats <= 75) & (lons >= -170) & (lons <= -50)
    non_na_mask = ~na_mask

    pre_yd_non_na = pre_yd_mask & non_na_mask
    on_non_na = on_corridor & non_na_mask
    obs_pre_non_na = int(np.sum(pre_yd_non_na & on_non_na))

    pre_yd_lats_nna = lats[pre_yd_non_na]
    pre_yd_lons_nna = lons[pre_yd_non_na]

    rand_pre_nna = []
    for plat, plon in random_poles[:200]:
        rd = gc_dist_from_pole(plat, plon, pre_yd_lats_nna, pre_yd_lons_nna)
        rand_pre_nna.append(int(np.sum(rd <= PRIMARY_THRESHOLD_KM)))
    rand_pre_nna = np.array(rand_pre_nna)
    z_nna, mu_nna, std_nna = zscore(obs_pre_non_na, rand_pre_nna)

    print(f"    Pre-YD on-corridor (excl. NA): {obs_pre_non_na}")
    print(f"    Z-score (excl. NA): {z_nna:.2f}")

    sensitivity_results['exclude_north_america'] = {
        'pre_yd_on': obs_pre_non_na,
        'z_score': round(z_nna, 4),
    }

    # ============================================================
    # GENERATE FIGURES
    # ============================================================
    print("\n" + "=" * 70)
    print("GENERATING FIGURES")
    print("=" * 70)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        # Figure 1: Temporal Enrichment Profile (THE MONEY PLOT)
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

        bin_mids = [(b['bin_start'] + b['bin_end']) / 2 for b in test1_results['bins']]
        enrichments = [b['enrichment'] for b in test1_results['bins']]
        on_counts = [b['on_corridor_count'] for b in test1_results['bins']]
        off_counts = [b['off_corridor_count'] for b in test1_results['bins']]

        # Top: enrichment ratio
        ax1.bar(bin_mids, enrichments, width=450, color='steelblue', alpha=0.7, edgecolor='navy')
        ax1.axhline(y=1.0, color='black', linestyle='--', linewidth=1, label='Expected (chance)')
        ax1.axvspan(YD_END, YD_START, alpha=0.2, color='cyan', label='Younger Dryas')
        ax1.set_ylabel('Enrichment Ratio\n(on-corridor / global, normalized)', fontsize=11)
        ax1.set_title('Pre-Younger Dryas Corridor Test — Temporal Enrichment Profile\n'
                      f'Great Circle corridor (±{PRIMARY_THRESHOLD_KM}km) | '
                      f'p3k14c + XRONOS + ROAD merged ({len(all_dates):,} dates)',
                      fontsize=13, fontweight='bold')
        ax1.legend(loc='upper left')
        ax1.set_ylim(bottom=0)

        # Bottom: raw counts
        ax2.bar(bin_mids, on_counts, width=450, color='crimson', alpha=0.7, label='On-corridor')
        ax2.axvspan(YD_END, YD_START, alpha=0.2, color='cyan')
        ax2.set_ylabel('Radiocarbon Date Count\n(on-corridor)', fontsize=11)
        ax2.set_xlabel('Age (14C years BP)', fontsize=11)
        ax2.legend(loc='upper left')
        ax2.invert_xaxis()

        plt.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "test1_temporal_enrichment.png"), dpi=150, bbox_inches='tight')
        plt.close(fig)
        print("  Saved: test1_temporal_enrichment.png")

        # Figure 2: Z-score over time
        fig, ax = plt.subplots(figsize=(14, 6))
        bz_mids = [(b['bin_start'] + b['bin_end']) / 2 for b in bin_zscores if b['z_score'] is not None]
        bz_vals = [b['z_score'] for b in bin_zscores if b['z_score'] is not None]

        colors = ['crimson' if z > 2 else 'orange' if z > 1 else 'steelblue' for z in bz_vals]
        ax.bar(bz_mids, bz_vals, width=450, color=colors, alpha=0.7, edgecolor='navy')
        ax.axhline(y=2.0, color='red', linestyle='--', linewidth=1, alpha=0.7, label='Z=2 threshold')
        ax.axhline(y=0, color='black', linewidth=0.5)
        ax.axvspan(YD_END, YD_START, alpha=0.2, color='cyan', label='Younger Dryas')
        ax.set_xlabel('Age (14C years BP)', fontsize=11)
        ax.set_ylabel('Z-score\n(vs 1000 random great circles)', fontsize=11)
        ax.set_title('Corridor Enrichment Z-Score by Temporal Bin\n'
                      f'Great Circle (±{PRIMARY_THRESHOLD_KM}km) vs {N_RANDOM_CIRCLES} random circles',
                      fontsize=13, fontweight='bold')
        ax.legend(loc='upper left')
        ax.invert_xaxis()

        plt.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "test2_zscore_over_time.png"), dpi=150, bbox_inches='tight')
        plt.close(fig)
        print("  Saved: test2_zscore_over_time.png")

        # Figure 3: Pre-YD Z-score histogram
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(random_pre_yd, bins=40, color='steelblue', alpha=0.7, edgecolor='navy', label='Random circles')
        ax.axvline(x=obs_pre_yd_on, color='crimson', linewidth=2, linestyle='--',
                   label=f'Alison GC: {obs_pre_yd_on} (Z={z_pre_yd:.2f})')
        ax.set_xlabel(f'Pre-YD dates within {PRIMARY_THRESHOLD_KM}km of circle', fontsize=11)
        ax.set_ylabel('Count (of 1000 random circles)', fontsize=11)
        ax.set_title(f'Test 2: Pre-YD Corridor Enrichment\nAlison Great Circle vs {N_RANDOM_CIRCLES} Random Great Circles',
                     fontsize=13, fontweight='bold')
        ax.legend(fontsize=11)

        plt.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "test2_pre_yd_histogram.png"), dpi=150, bbox_inches='tight')
        plt.close(fig)
        print("  Saved: test2_pre_yd_histogram.png")

        # Figure 4: YD Crash comparison
        fig, ax = plt.subplots(figsize=(8, 6))
        bars = ax.bar(['Alison\nCorridor', 'Global\n(off-corridor)', f'Random\n(mean±1σ)'],
                      [crash_corridor, crash_global, mu_crash],
                      yerr=[0, 0, std_crash],
                      color=['crimson', 'steelblue', 'gray'], alpha=0.7, edgecolor='navy',
                      capsize=5)
        ax.set_ylabel('Crash Ratio\n(post-YD onset / pre-YD window)', fontsize=11)
        ax.set_title('Test 3: YD Disruption — Crash Ratio Comparison\n'
                     '(13300-12800 BP → 12800-12300 BP)', fontsize=13, fontweight='bold')
        ax.axhline(y=1.0, color='black', linestyle=':', linewidth=1, alpha=0.5)

        plt.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "test3_yd_crash.png"), dpi=150, bbox_inches='tight')
        plt.close(fig)
        print("  Saved: test3_yd_crash.png")

        # Figure 5: Regional decomposition
        regions_sorted = sorted(region_map.keys(), key=lambda r: region_map[r]['pre_yd'], reverse=True)
        regions_sorted = [r for r in regions_sorted if region_map[r]['total'] > 0]

        fig, ax = plt.subplots(figsize=(12, 6))
        x = range(len(regions_sorted))
        pre_vals = [region_map[r]['pre_yd'] for r in regions_sorted]
        post_vals = [region_map[r]['post_yd'] for r in regions_sorted]

        ax.bar(x, pre_vals, color='crimson', alpha=0.7, label='Pre-YD (>12800 BP)')
        ax.bar(x, post_vals, bottom=pre_vals, color='steelblue', alpha=0.7, label='Post-YD (≤12800 BP)')
        ax.set_xticks(x)
        ax.set_xticklabels(regions_sorted, rotation=45, ha='right')
        ax.set_ylabel('Radiocarbon Date Count', fontsize=11)
        ax.set_title('Test 4: Regional Decomposition of On-Corridor Dates', fontsize=13, fontweight='bold')
        ax.legend()

        plt.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "test4_regional_decomposition.png"), dpi=150, bbox_inches='tight')
        plt.close(fig)
        print("  Saved: test4_regional_decomposition.png")

        # Figure 6: SPD comparison
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

        spd_mids = [(lo + hi) / 2 for lo, hi in SPD_BINS]

        ax1.fill_between(spd_mids, spd_on_norm, alpha=0.5, color='crimson', label='On-corridor')
        ax1.fill_between(spd_mids, spd_off_norm, alpha=0.3, color='steelblue', label='Off-corridor')
        ax1.axvspan(YD_END, YD_START, alpha=0.15, color='cyan')
        ax1.set_ylabel('Normalized Density', fontsize=11)
        ax1.set_title('Test 7: Simplified SPD — On-Corridor vs Off-Corridor\n'
                      '(uncalibrated 200-year bins)', fontsize=13, fontweight='bold')
        ax1.legend()

        # Ratio plot
        ratio_spd = np.where(spd_off_norm > 0, spd_on_norm / spd_off_norm, 0)
        ax2.bar(spd_mids, ratio_spd, width=180, color='purple', alpha=0.5)
        ax2.axhline(y=1.0, color='black', linestyle='--')
        ax2.axvspan(YD_END, YD_START, alpha=0.15, color='cyan')
        ax2.set_xlabel('Age (14C years BP)', fontsize=11)
        ax2.set_ylabel('On/Off Density Ratio', fontsize=11)
        ax2.set_ylim(0, max(3, np.max(ratio_spd[ratio_spd < 10]) * 1.1) if np.any(ratio_spd < 10) else 3)
        ax2.invert_xaxis()

        plt.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "test7_spd_comparison.png"), dpi=150, bbox_inches='tight')
        plt.close(fig)
        print("  Saved: test7_spd_comparison.png")

    except ImportError as e:
        print(f"  WARNING: matplotlib not available ({e}). Skipping figures.")

    # ============================================================
    # FINAL SUMMARY
    # ============================================================
    elapsed = time.time() - t0

    print("\n" + "=" * 70)
    print("FINAL SUMMARY — PRE-YOUNGER DRYAS CORRIDOR TEST")
    print("=" * 70)

    summary = {
        'meta': {
            'date': '2026-03-20',
            'databases': ['p3k14c', 'XRONOS', 'ROAD v32'],
            'total_merged_dates': len(all_dates),
            'pre_yd_dates': int(np.sum(pre_yd_mask)),
            'on_corridor_dates': int(np.sum(on_corridor)),
            'threshold_km': PRIMARY_THRESHOLD_KM,
            'n_random_circles': N_RANDOM_CIRCLES,
            'runtime_seconds': round(elapsed, 1),
        },
        'test1_temporal_profile': {
            'description': 'Enrichment of on-corridor dates by 500-year bin',
            'file': 'test1_temporal_profile.json',
        },
        'test2_pre_yd_enrichment': {
            'description': 'THE KEY TEST: Pre-YD corridor enrichment vs random circles',
            'pre_yd_z_score': round(z_pre_yd, 4),
            'deep_pre_yd_z_score': round(z_deep, 4),
            'late_glacial_z_score': round(z_lg, 4),
            'habitability_matched_z': round(z_hab, 4) if n_hab_matched > 10 else None,
            'observed_pre_yd_on_corridor': obs_pre_yd_on,
            'random_mean': round(mu_pre, 2),
            'random_std': round(std_pre, 2),
            'percentile': round(pctile_pre, 2),
            'file': 'test2_pre_yd_enrichment.json',
        },
        'test3_yd_disruption': {
            'description': 'YD crash ratio: corridor vs global',
            'crash_ratio_corridor': round(crash_corridor, 4),
            'crash_ratio_global': round(crash_global, 4),
            'z_score': round(z_crash, 4) if len(random_crash) > 10 else None,
            'file': 'test3_yd_disruption.json',
        },
        'test4_regional_decomposition': {
            'description': 'Regional breakdown of on-corridor pre-YD dates',
            'top_region': max(region_map.keys(), key=lambda r: region_map[r]['pre_yd']),
            'top_region_pre_yd': max(region_map.values(), key=lambda v: v['pre_yd'])['pre_yd'],
            'file': 'test4_regional_decomposition.json',
        },
        'test5_continuity': {
            'description': 'Site continuity across YD boundary',
            'continuity_rate_on': round(cont_rate_on, 6),
            'continuity_rate_off': round(cont_rate_off, 6),
            'file': 'test5_continuity.json',
        },
        'test6_atlantis_regions': {
            'description': 'Pre-YD activity in candidate "Atlantis" regions',
            'file': 'test6_atlantis_regions.json',
        },
        'test7_spd': {
            'description': 'Simplified summed probability distribution',
            'file': 'test7_spd.json',
        },
        'sensitivity': sensitivity_results,
    }

    # Verdict
    print(f"\n  Test 2 (THE KEY TEST): Pre-YD Z = {z_pre_yd:.2f}")
    if z_pre_yd < 2:
        verdict = "CONVENTIONAL — No pre-YD anomaly detected"
        print(f"  VERDICT: {verdict}")
        print(f"  The corridor shows no statistically significant enrichment of human")
        print(f"  activity before the Younger Dryas compared to random great circles.")
    elif z_pre_yd < 3:
        verdict = "MILD ENRICHMENT — Needs careful habitability control"
        print(f"  VERDICT: {verdict}")
    elif z_pre_yd < 5:
        verdict = "SIGNIFICANT — Requires scrutiny and replication"
        print(f"  VERDICT: {verdict}")
    else:
        verdict = "EXTRAORDINARY — Triple-check everything"
        print(f"  VERDICT: {verdict}")

    summary['verdict'] = verdict

    print(f"\n  Test 3 (YD disruption): crash Z = {z_crash:.2f}" if len(random_crash) > 10 else "")
    print(f"  Test 4 (regions): Top contributor = {summary['test4_regional_decomposition']['top_region']}")
    print(f"  Test 5 (continuity): on={cont_rate_on*100:.2f}%, off={cont_rate_off*100:.2f}%")

    print(f"\n  Runtime: {elapsed:.1f} seconds")

    # Save master summary
    with open(os.path.join(OUT_DIR, "RESULTS_SUMMARY.json"), 'w') as f:
        json.dump(summary, f, indent=2)

    # Write RESULTS.md
    results_md = f"""# Pre-Younger Dryas Corridor Test — Results

**Date:** 2026-03-20
**Databases:** p3k14c ({len(p3k):,}) + XRONOS ({len(xronos):,}) + ROAD v32 ({len(road):,})
**After deduplication:** {len(all_dates):,} dates ({int(np.sum(pre_yd_mask)):,} pre-YD)
**Corridor threshold:** {PRIMARY_THRESHOLD_KM}km | Random circles: {N_RANDOM_CIRCLES}
**Runtime:** {elapsed:.1f}s

---

## Verdict: {verdict}

---

## Test 2 — THE KEY RESULT: Pre-YD Corridor Enrichment

| Metric | Observed | Random Mean ± SD | Z-score | Percentile |
|--------|----------|-----------------|---------|------------|
| Pre-YD (>12800 BP) | {obs_pre_yd_on} | {mu_pre:.1f} ± {std_pre:.1f} | **{z_pre_yd:.2f}** | {pctile_pre:.1f}% |
| Deep pre-YD (>15000 BP) | {obs_deep_on} | {mu_deep:.1f} ± {std_deep:.1f} | {z_deep:.2f} | {pctile_deep:.1f}% |
| Late Glacial (12800-15000) | {obs_late_glacial_on} | {mu_lg:.1f} ± {std_lg:.1f} | {z_lg:.2f} | {pctile_lg:.1f}% |
| Total (all periods) | {obs_total_on} | {mu_total:.1f} ± {std_total:.1f} | {z_total:.2f} | — |
"""
    if n_hab_matched > 10:
        results_md += f"| Habitability-matched pre-YD | {obs_pre_yd_on} | {mu_hab:.1f} ± {std_hab:.1f} | {z_hab:.2f} | {pctile_hab:.1f}% |\n"

    results_md += f"""
## Test 3 — YD Disruption

| Metric | Corridor | Global | Z-score |
|--------|----------|--------|---------|
| Crash ratio (post-YD onset / pre-YD) | {crash_corridor:.3f} | {crash_global:.3f} | {z_crash:.2f} |
| Recovery ratio | {recovery_corridor:.3f} | {recovery_global:.3f} | — |

## Test 4 — Regional Decomposition

"""
    for region in sorted(region_map.keys(), key=lambda r: region_map[r]['pre_yd'], reverse=True):
        d = region_map[region]
        if d['total'] > 0:
            pct = d['pre_yd'] / d['total'] * 100
            results_md += f"| {region} | {d['pre_yd']} pre-YD | {d['post_yd']} post-YD | {pct:.1f}% pre-YD |\n"

    results_md += f"""
## Test 5 — Continuity

- On-corridor continuity rate: {cont_rate_on*100:.2f}%
- Off-corridor continuity rate: {cont_rate_off*100:.2f}%
- Continuity sites spanning the YD: {continuity_on} on-corridor, {continuity_off} off-corridor

## Test 6 — Atlantis Regions

"""
    for r in region_results:
        results_md += f"- **{r['name']}**: {r['total_dates']} dates, {r['pre_yd_dates']} pre-YD, ratio vs global = {r['ratio_vs_global']:.2f}×\n"

    results_md += f"""
## Sensitivity Checks

"""
    for k, v in sensitivity_results.items():
        if isinstance(v, dict):
            results_md += f"- **{k}**: Z = {v.get('z_score', 'N/A')}\n"

    results_md += f"""
## Caveats

1. **Calibration:** Tests 1-6 use uncalibrated 14C ages. The YD-era calibration plateau may compress some dates.
2. **p3k14c North America coordinates:** Fuzzed to county centroids. Sensitivity test excluding NA: Z = {z_nna:.2f}
3. **Submerged regions:** Persian Gulf, Sunda Shelf, Doggerland have minimal data due to submersion, not necessarily absence of activity.
4. **Sampling bias:** The Natufian, PPNA, and Clovis are heavily sampled. Matched-corridor comparison controls for this globally but not regionally.
5. **Multiple comparisons:** {len(BINS)} temporal bins tested; apply Bonferroni or FDR correction for individual bin significance.

## Files

- `test1_temporal_profile.json` — Enrichment time series
- `test2_pre_yd_enrichment.json` — Pre-YD Z-scores and random distributions
- `test3_yd_disruption.json` — YD crash ratio analysis
- `test4_regional_decomposition.json` — Regional breakdown
- `test5_continuity.json` — Site continuity across YD
- `test6_atlantis_regions.json` — Candidate region analysis
- `test7_spd.json` — Simplified SPD
- `figures/` — Visualizations
"""

    with open(os.path.join(OUT_DIR, "RESULTS.md"), 'w') as f:
        f.write(results_md)

    print(f"\n  Saved: RESULTS.md")
    print(f"\n  ALL OUTPUTS WRITTEN TO: {OUT_DIR}")
    print(f"\n{'='*70}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
