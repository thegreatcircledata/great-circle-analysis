#!/usr/bin/env python3
"""
Cross-Continental Temporal Divergence Test
============================================
Tests whether the monument-settlement divergence appears independently on
different continents along the Great Circle, and whether temporal peaks
are synchronized across continents with zero known contact.

Steps:
  1. Regional temporal divergence (5 regions, 500-year bins)
  2. Cross-continental synchrony test (correlation of D time series)
  3. South American deep dive (Caral/Norte Chico contemporaneity)
  4. Caral-specific test (distance from GC, temporal spike)
  5. Null expectation (100 random great circles)
"""

import csv, math, random, json, os, sys, time
import numpy as np
from collections import defaultdict, Counter

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
THRESHOLD_KM = 50
REGIONAL_THRESHOLD_KM = 200  # for filtering sites near the circle
N_TRIALS = 200
MIN_SITES = 10  # minimum sites per bin to compute Z

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(OUT_DIR, exist_ok=True)

# 500-year bins from 5000 BCE to 500 CE
BIN_EDGES = list(range(-5000, 1000, 500))
BINS = [(BIN_EDGES[i], BIN_EDGES[i+1]) for i in range(len(BIN_EDGES)-1)]

# Region definitions: name -> (lat_min, lat_max, lon_min, lon_max, use_gc_filter)
REGIONS = {
    "egypt_levant":    {"lat": (25, 35), "lon": (28, 38), "gc_filter": False,
                        "label": "Egypt/Levant"},
    "iran_mesopotamia":{"lat": (28, 38), "lon": (44, 60), "gc_filter": False,
                        "label": "Iran/Mesopotamia"},
    "south_america":   {"lat": (-56, 5), "lon": (-82, -34), "gc_filter": True,
                        "label": "South America"},
    "western_europe":  {"lat": (36, 60), "lon": (-10, 5), "gc_filter": True,
                        "label": "Western Europe"},
    "south_se_asia":   {"lat": (0, 30),  "lon": (70, 110), "gc_filter": True,
                        "label": "South/SE Asia"},
}

# Classification keywords
MONUMENTAL_KW = ['temple', 'pyramid', 'tomb', 'mound', 'monument', 'cairn', 'megalith',
                  'henge', 'huaca', 'ceremonial', 'ritual', 'geoglyph', 'ahu', 'moai',
                  'stela', 'ziggurat', 'palace', 'citadel', 'fortress', 'dolmen',
                  'stone circle', 'standing stone', 'barrow', 'tumulus', 'sanctuary',
                  'necropolis', 'acropolis', 'amphitheatre', 'colosseum']

DOMESTIC_KW = ['village', 'settlement', 'farm', 'camp', 'midden', 'shell mound',
               'quarry', 'mine', 'shelter', 'rockshelter', 'workshop', 'kiln', 'pit',
               'habitation', 'dwelling', 'house', 'homestead', 'terrace', 'field',
               'granary', 'storehouse', 'well', 'cistern', 'fishweir']

PLEIADES_MONUMENTAL = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe",
    "tumulus", "shrine",
}

PLEIADES_SETTLEMENT = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement",
    "townhouse", "production",
}

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def gc_dist_vec(pole_lat, pole_lon, site_lats, site_lons):
    lat1_r, lon1_r = np.radians(pole_lat), np.radians(pole_lon)
    lat2_r, lon2_r = np.radians(site_lats), np.radians(site_lons)
    dlat = lat2_r - lat1_r; dlon = lon2_r - lon1_r
    a = np.sin(dlat/2)**2 + np.cos(lat1_r)*np.cos(lat2_r)*np.sin(dlon/2)**2
    d = EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QUARTER_CIRC)

def gc_dist_single(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    return float(gc_dist_vec(pole_lat, pole_lon, np.array([lat]), np.array([lon]))[0])

def rand_matched_batch(site_lats, site_lons, n):
    idx = np.random.randint(0, len(site_lats), size=n)
    lats = site_lats[idx] + np.random.normal(0, 2, n)
    lons = site_lons[idx] + np.random.normal(0, 2, n)
    return np.clip(lats, -90, 90), np.clip(lons, -180, 180)

def compute_z(pole_lat, pole_lon, lats, lons, threshold_km=THRESHOLD_KM, n_trials=N_TRIALS):
    n = len(lats)
    if n < MIN_SITES:
        return {"n": int(n), "z_score": None, "too_few": True}
    dists = gc_dist_vec(pole_lat, pole_lon, lats, lons)
    observed = int(np.sum(dists <= threshold_km))
    rand_counts = []
    for _ in range(n_trials):
        rl, rn = rand_matched_batch(lats, lons, n)
        rd = gc_dist_vec(pole_lat, pole_lon, rl, rn)
        rand_counts.append(int(np.sum(rd <= threshold_km)))
    mu = np.mean(rand_counts); sigma = np.std(rand_counts)
    z = float((observed - mu) / sigma) if sigma > 0 else 0.0
    enrich = float(observed / mu) if mu > 0 else 0.0
    p_val = float(np.sum(np.array(rand_counts) >= observed) / n_trials)
    return {"n": int(n), "observed": int(observed), "expected": round(float(mu), 2),
            "std": round(float(sigma), 2), "z_score": round(z, 2),
            "enrichment": round(enrich, 3), "p_value": round(p_val, 4), "too_few": False}

def c14_to_cal(age_bp):
    if age_bp <= 0:
        return 1950
    if age_bp < 2500:
        cal_bp = age_bp * 1.0
    elif age_bp < 5000:
        cal_bp = age_bp * 1.05 + 50
    elif age_bp < 8000:
        cal_bp = age_bp * 1.08 + 100
    elif age_bp < 12000:
        cal_bp = age_bp * 1.12 + 200
    else:
        cal_bp = age_bp * 1.15 + 300
    return 1950 - cal_bp

def classify_p3k(site_name):
    if not site_name:
        return "unclassified"
    name_lower = site_name.lower()
    for kw in MONUMENTAL_KW:
        if kw in name_lower:
            return "monumental"
    for kw in DOMESTIC_KW:
        if kw in name_lower:
            return "domestic"
    return "unclassified"

def bin_label(start, end):
    def yr(y):
        return f"{abs(y)} BCE" if y < 0 else f"{y} CE"
    return f"{yr(start)} – {yr(end)}"

def in_region(lat, lon, region):
    lat_min, lat_max = region["lat"]
    lon_min, lon_max = region["lon"]
    return lat_min <= lat <= lat_max and lon_min <= lon <= lon_max

def random_pole():
    """Generate a random great circle pole on the sphere."""
    lat = np.degrees(np.arcsin(2 * random.random() - 1))
    lon = random.uniform(-180, 180)
    return lat, lon

def pearson_corr(x, y):
    """Pearson correlation for two arrays, ignoring NaN pairs."""
    mask = ~(np.isnan(x) | np.isnan(y))
    x_clean = x[mask]; y_clean = y[mask]
    if len(x_clean) < 3:
        return None, 0
    mx = np.mean(x_clean); my = np.mean(y_clean)
    sx = np.std(x_clean); sy = np.std(y_clean)
    if sx == 0 or sy == 0:
        return 0.0, len(x_clean)
    r = np.mean((x_clean - mx) * (y_clean - my)) / (sx * sy)
    return round(float(r), 3), int(len(x_clean))

# ============================================================
# MAIN
# ============================================================
print("=" * 70)
print("CROSS-CONTINENTAL TEMPORAL DIVERGENCE TEST")
print("=" * 70)

t0 = time.time()
np.random.seed(42)
random.seed(42)

# ============================================================
# LOAD DATA
# ============================================================
print("\n--- Loading data ---")

# Load Pleiades
pleiades_file = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")
pleiades_sites = []
with open(pleiades_file, encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            lat = float(row['reprLat']); lon = float(row['reprLong'])
        except (ValueError, KeyError, TypeError):
            continue
        if lat == 0 and lon == 0: continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180): continue
        try:
            min_date = float(row.get('minDate', ''))
        except (ValueError, TypeError):
            continue
        feature_types = {t.strip() for t in row.get('featureTypes', '').split(',') if t.strip()}
        is_mono = bool(feature_types & PLEIADES_MONUMENTAL)
        is_settle = bool(feature_types & PLEIADES_SETTLEMENT)
        if not is_mono and not is_settle: continue
        pleiades_sites.append({
            'lat': lat, 'lon': lon, 'min_date': min_date,
            'classification': 'monumental' if is_mono else 'settlement',
            'source': 'pleiades'
        })
print(f"  Pleiades: {len(pleiades_sites)} classified sites with dates")

# Load p3k14c
p3k_file = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")
site_groups = defaultdict(list)
with open(p3k_file, encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            lat = float(row['Lat']); lon = float(row['Long'])
            age_bp = float(row['Age'])
        except (ValueError, TypeError, KeyError):
            continue
        if lat == 0 and lon == 0: continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180): continue
        site_id = row.get('SiteID', '') or f"anon_{lat:.4f}_{lon:.4f}"
        site_groups[site_id].append({
            'lat': lat, 'lon': lon, 'age_bp': age_bp,
            'site_name': row.get('SiteName', ''),
        })

p3k_sites = []
for sid, rows in site_groups.items():
    ages = [r['age_bp'] for r in rows]
    cls = classify_p3k(rows[0]['site_name'])
    if cls == "unclassified": continue
    cal_year = c14_to_cal(max(ages))
    p3k_sites.append({
        'lat': np.mean([r['lat'] for r in rows]),
        'lon': np.mean([r['lon'] for r in rows]),
        'cal_year': cal_year,
        'site_name': rows[0]['site_name'],
        'classification': 'monumental' if cls == 'monumental' else 'settlement',
        'source': 'p3k14c'
    })
print(f"  p3k14c: {len(p3k_sites)} classified sites")

# Combine all sites
all_sites = []
for s in pleiades_sites:
    all_sites.append({**s, 'cal_year': s['min_date']})
for s in p3k_sites:
    all_sites.append(s)
print(f"  Combined: {len(all_sites)} sites")

# ============================================================
# STEP 1: REGIONAL TEMPORAL DIVERGENCE
# ============================================================
print("\n" + "=" * 70)
print("STEP 1: REGIONAL TEMPORAL DIVERGENCE")
print("=" * 70)

def compute_regional_bins(sites, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    """Compute D time series for a set of sites."""
    results = []
    for start, end in BINS:
        mono_in_bin = [(s['lat'], s['lon']) for s in sites
                       if s['classification'] == 'monumental' and start <= s['cal_year'] < end]
        settle_in_bin = [(s['lat'], s['lon']) for s in sites
                         if s['classification'] == 'settlement' and start <= s['cal_year'] < end]

        m_z = None; s_z = None; D = None
        m_res = None; s_res = None

        if len(mono_in_bin) >= MIN_SITES:
            m_lats = np.array([p[0] for p in mono_in_bin])
            m_lons = np.array([p[1] for p in mono_in_bin])
            m_res = compute_z(pole_lat, pole_lon, m_lats, m_lons)
            if not m_res.get('too_few'):
                m_z = m_res['z_score']

        if len(settle_in_bin) >= MIN_SITES:
            s_lats = np.array([p[0] for p in settle_in_bin])
            s_lons = np.array([p[1] for p in settle_in_bin])
            s_res = compute_z(pole_lat, pole_lon, s_lats, s_lons)
            if not s_res.get('too_few'):
                s_z = s_res['z_score']

        if m_z is not None and s_z is not None:
            D = round(m_z - s_z, 2)

        results.append({
            "bin_start": start, "bin_end": end,
            "label": bin_label(start, end),
            "monument_n": len(mono_in_bin),
            "settlement_n": len(settle_in_bin),
            "monument_z": m_z, "settlement_z": s_z,
            "divergence_D": D,
            "monument_detail": m_res,
            "settlement_detail": s_res,
        })
    return results

regional_results = {}

for region_key, region_def in REGIONS.items():
    print(f"\n  --- {region_def['label']} ---")

    # Filter sites to region
    region_sites = []
    for s in all_sites:
        if not in_region(s['lat'], s['lon'], region_def):
            continue
        # If gc_filter, only include sites within 200km of the circle
        if region_def['gc_filter']:
            d = gc_dist_single(s['lat'], s['lon'])
            if d > REGIONAL_THRESHOLD_KM:
                continue
        region_sites.append(s)

    n_mono = sum(1 for s in region_sites if s['classification'] == 'monumental')
    n_settle = sum(1 for s in region_sites if s['classification'] == 'settlement')
    print(f"    Sites: {len(region_sites)} (mono={n_mono}, settle={n_settle})")

    bins = compute_regional_bins(region_sites)
    regional_results[region_key] = {
        "label": region_def["label"],
        "n_sites": len(region_sites),
        "n_monumental": n_mono,
        "n_settlement": n_settle,
        "bins": bins,
    }

    # Print summary of D values
    for b in bins:
        if b['divergence_D'] is not None:
            print(f"    {b['label']}: D={b['divergence_D']:+.1f}  "
                  f"(mono={b['monument_n']}, settle={b['settlement_n']})")

# ============================================================
# STEP 2: CROSS-CONTINENTAL SYNCHRONY TEST
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: CROSS-CONTINENTAL SYNCHRONY TEST")
print("=" * 70)

# Extract D time series for each region
def get_d_series(region_key):
    bins = regional_results[region_key]["bins"]
    return np.array([b['divergence_D'] if b['divergence_D'] is not None else np.nan
                     for b in bins])

region_d_series = {k: get_d_series(k) for k in REGIONS}
bin_centers = np.array([(s + e) / 2 for s, e in BINS])

# Find peak D bin for each region
print("\n  Peak divergence bins:")
peak_bins = {}
for rk, label in [(k, REGIONS[k]['label']) for k in REGIONS]:
    d_arr = region_d_series[rk]
    valid = ~np.isnan(d_arr)
    if np.any(valid):
        peak_idx = np.nanargmax(d_arr)
        peak_d = d_arr[peak_idx]
        peak_bin = BINS[peak_idx]
        peak_bins[rk] = {"bin": peak_bin, "D": round(float(peak_d), 2),
                         "label": bin_label(*peak_bin)}
        print(f"    {label}: peak D = {peak_d:.2f} at {bin_label(*peak_bin)}")
    else:
        peak_bins[rk] = {"bin": None, "D": None, "label": "N/A"}
        print(f"    {label}: no valid D values")

# Pairwise correlations
print("\n  Pairwise Pearson correlations of D(t):")
region_keys = list(REGIONS.keys())
correlation_matrix = {}
for i in range(len(region_keys)):
    for j in range(i+1, len(region_keys)):
        rk1, rk2 = region_keys[i], region_keys[j]
        r, n = pearson_corr(region_d_series[rk1], region_d_series[rk2])
        pair_key = f"{rk1}_vs_{rk2}"
        correlation_matrix[pair_key] = {"r": r, "n_bins": n,
                                         "labels": [REGIONS[rk1]['label'], REGIONS[rk2]['label']]}
        if r is not None:
            print(f"    {REGIONS[rk1]['label']} vs {REGIONS[rk2]['label']}: "
                  f"r = {r:.3f} (n={n} bins)")
        else:
            print(f"    {REGIONS[rk1]['label']} vs {REGIONS[rk2]['label']}: "
                  f"insufficient data")

# ============================================================
# STEP 3: SOUTH AMERICAN DEEP DIVE
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: SOUTH AMERICAN DEEP DIVE")
print("=" * 70)

sa_sites = [s for s in all_sites
            if in_region(s['lat'], s['lon'], REGIONS['south_america'])
            and gc_dist_single(s['lat'], s['lon']) <= REGIONAL_THRESHOLD_KM]

sa_mono = [s for s in sa_sites if s['classification'] == 'monumental']
sa_settle = [s for s in sa_sites if s['classification'] == 'settlement']

print(f"  SA sites within 200km of GC: {len(sa_sites)} (mono={len(sa_mono)}, settle={len(sa_settle)})")

# List monumental sites with dates
print("\n  SA monumental sites (near GC):")
sa_mono_sorted = sorted(sa_mono, key=lambda s: s['cal_year'])
for s in sa_mono_sorted[:30]:
    yr = s['cal_year']
    yr_str = f"{abs(yr):.0f} BCE" if yr < 0 else f"{yr:.0f} CE"
    d = gc_dist_single(s['lat'], s['lon'])
    print(f"    {yr_str:>10s}  {d:6.1f}km  {s.get('site_name', s.get('name', 'unknown'))}")

# Check for Caral-era monument construction (3000-2000 BCE)
caral_era_mono = [s for s in sa_mono if -3000 <= s['cal_year'] < -2000]
caral_era_settle = [s for s in sa_settle if -3000 <= s['cal_year'] < -2000]
print(f"\n  Caral era (3000-2000 BCE):")
print(f"    Monumental: {len(caral_era_mono)}")
print(f"    Settlement: {len(caral_era_settle)}")

# Temporal distribution of SA monuments near GC
sa_mono_bins = defaultdict(int)
sa_settle_bins = defaultdict(int)
for s in sa_mono:
    for start, end in BINS:
        if start <= s['cal_year'] < end:
            sa_mono_bins[start] += 1
            break
for s in sa_settle:
    for start, end in BINS:
        if start <= s['cal_year'] < end:
            sa_settle_bins[start] += 1
            break

print("\n  SA temporal distribution (near GC):")
for start, end in BINS:
    m = sa_mono_bins.get(start, 0)
    se = sa_settle_bins.get(start, 0)
    if m > 0 or se > 0:
        print(f"    {bin_label(start, end)}: mono={m}, settle={se}")

# ============================================================
# STEP 4: CARAL-SPECIFIC TEST
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: CARAL-SPECIFIC TEST")
print("=" * 70)

CARAL_LAT = -10.89
CARAL_LON = -77.52
caral_gc_dist = gc_dist_single(CARAL_LAT, CARAL_LON)
print(f"  Caral coordinates: {CARAL_LAT}°S, {CARAL_LON}°W")
print(f"  Distance from Great Circle: {caral_gc_dist:.1f} km")
print(f"  Within 50km: {'YES' if caral_gc_dist <= 50 else 'NO'}")

# Coastal Peru sites near GC (within 100km of coast, within 200km of GC)
coastal_peru = [s for s in all_sites
                if -18 <= s['lat'] <= -3 and -82 <= s['lon'] <= -75
                and gc_dist_single(s['lat'], s['lon']) <= 200]

print(f"\n  Coastal Peru sites within 200km of GC: {len(coastal_peru)}")

# Temporal spike analysis
peru_3k_2k = [s for s in coastal_peru if -3000 <= s['cal_year'] < -2000]
peru_other = [s for s in coastal_peru if not (-3000 <= s['cal_year'] < -2000)]
print(f"  Sites dated 3000-2000 BCE: {len(peru_3k_2k)}")
print(f"  Sites dated other periods: {len(peru_other)}")

# Bin distribution
peru_bins = defaultdict(int)
for s in coastal_peru:
    for start, end in BINS:
        if start <= s['cal_year'] < end:
            peru_bins[start] += 1
            break

print("\n  Coastal Peru temporal distribution (near GC):")
for start, end in BINS:
    n = peru_bins.get(start, 0)
    if n > 0:
        print(f"    {bin_label(start, end)}: {n} sites")

# Monument/settlement split for 3000-2000 BCE coastal Peru
peru_3k_mono = [s for s in peru_3k_2k if s['classification'] == 'monumental']
peru_3k_settle = [s for s in peru_3k_2k if s['classification'] == 'settlement']
print(f"\n  3000-2000 BCE coastal Peru near GC:")
print(f"    Monumental: {len(peru_3k_mono)}")
print(f"    Settlement: {len(peru_3k_settle)}")

caral_test = {
    "caral_coords": {"lat": CARAL_LAT, "lon": CARAL_LON},
    "gc_distance_km": round(caral_gc_dist, 1),
    "within_50km": caral_gc_dist <= 50,
    "coastal_peru_near_gc": len(coastal_peru),
    "sites_3000_2000_bce": len(peru_3k_2k),
    "sites_other_periods": len(peru_other),
    "mono_3000_2000_bce": len(peru_3k_mono),
    "settle_3000_2000_bce": len(peru_3k_settle),
    "temporal_distribution": {bin_label(s, e): peru_bins.get(s, 0) for s, e in BINS if peru_bins.get(s, 0) > 0},
}

# ============================================================
# STEP 5: NULL EXPECTATION (100 random great circles)
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: NULL EXPECTATION — 100 random great circles")
print("=" * 70)

N_RANDOM = 100
random_sync_count = 0
random_corrs = []

# For speed, pre-filter sites to the regions once
region_site_pools = {}
for rk, rdef in REGIONS.items():
    pool = [s for s in all_sites if in_region(s['lat'], s['lon'], rdef)]
    region_site_pools[rk] = pool

for trial in range(N_RANDOM):
    if (trial + 1) % 10 == 0:
        print(f"  Random circle {trial+1}/{N_RANDOM}...")

    rpole_lat, rpole_lon = random_pole()

    # For each region, filter by 200km from this random circle and compute D series
    trial_d_series = {}
    for rk, rdef in REGIONS.items():
        pool = region_site_pools[rk]
        if rdef['gc_filter']:
            # Filter by distance from random circle
            filtered = []
            for s in pool:
                d = gc_dist_single(s['lat'], s['lon'], rpole_lat, rpole_lon)
                if d <= REGIONAL_THRESHOLD_KM:
                    filtered.append(s)
        else:
            filtered = pool

        # Compute D time series for this region with the random circle
        d_series = []
        for start, end in BINS:
            mono_bin = [(s['lat'], s['lon']) for s in filtered
                        if s['classification'] == 'monumental' and start <= s['cal_year'] < end]
            settle_bin = [(s['lat'], s['lon']) for s in filtered
                          if s['classification'] == 'settlement' and start <= s['cal_year'] < end]

            D = np.nan
            if len(mono_bin) >= MIN_SITES and len(settle_bin) >= MIN_SITES:
                m_lats = np.array([p[0] for p in mono_bin])
                m_lons = np.array([p[1] for p in mono_bin])
                s_lats = np.array([p[0] for p in settle_bin])
                s_lons = np.array([p[1] for p in settle_bin])
                m_res = compute_z(rpole_lat, rpole_lon, m_lats, m_lons, n_trials=50)
                s_res = compute_z(rpole_lat, rpole_lon, s_lats, s_lons, n_trials=50)
                if not m_res.get('too_few') and not s_res.get('too_few'):
                    if m_res['z_score'] is not None and s_res['z_score'] is not None:
                        D = m_res['z_score'] - s_res['z_score']
            d_series.append(D)

        trial_d_series[rk] = np.array(d_series)

    # Check for synchronized peaks across any 2+ regions
    # Find peak bin for each region (if valid)
    trial_peaks = {}
    for rk in REGIONS:
        d_arr = trial_d_series[rk]
        valid = ~np.isnan(d_arr)
        if np.any(valid):
            peak_idx = int(np.nanargmax(d_arr))
            trial_peaks[rk] = peak_idx

    # Count how many regions share the same peak bin
    if len(trial_peaks) >= 2:
        peak_bins_count = Counter(trial_peaks.values())
        max_shared = max(peak_bins_count.values())
        if max_shared >= 2:
            random_sync_count += 1

    # Compute pairwise correlations for this random circle
    trial_corrs = []
    for i in range(len(region_keys)):
        for j in range(i+1, len(region_keys)):
            rk1, rk2 = region_keys[i], region_keys[j]
            r, n = pearson_corr(trial_d_series[rk1], trial_d_series[rk2])
            if r is not None:
                trial_corrs.append(r)
    if trial_corrs:
        random_corrs.append(np.mean(trial_corrs))

random_sync_pct = random_sync_count / N_RANDOM * 100
random_corr_mean = float(np.mean(random_corrs)) if random_corrs else 0
random_corr_std = float(np.std(random_corrs)) if random_corrs else 0

print(f"\n  Random circles with synchronized D peaks (2+ regions): "
      f"{random_sync_count}/{N_RANDOM} ({random_sync_pct:.1f}%)")
print(f"  Mean pairwise correlation across random circles: "
      f"{random_corr_mean:.3f} ± {random_corr_std:.3f}")

# Compute the actual circle's mean pairwise correlation for comparison
actual_corrs = []
for pk, pv in correlation_matrix.items():
    if pv['r'] is not None:
        actual_corrs.append(pv['r'])
actual_mean_corr = float(np.mean(actual_corrs)) if actual_corrs else 0

# Check if actual peak bins are synchronized
actual_peak_indices = {}
for rk in REGIONS:
    d_arr = region_d_series[rk]
    valid = ~np.isnan(d_arr)
    if np.any(valid):
        actual_peak_indices[rk] = int(np.nanargmax(d_arr))

actual_peak_counter = Counter(actual_peak_indices.values())
actual_max_shared = max(actual_peak_counter.values()) if actual_peak_counter else 0

print(f"\n  Actual circle mean pairwise correlation: {actual_mean_corr:.3f}")
print(f"  Actual circle max shared peak bin: {actual_max_shared} regions")
print(f"  Random circles that match or exceed this: "
      f"{sum(1 for c in random_corrs if c >= actual_mean_corr)}/{N_RANDOM}")

# ============================================================
# PLOTS
# ============================================================
print("\n--- Generating plots ---")

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # Plot 1: Regional D time series overlaid
    fig, ax = plt.subplots(1, 1, figsize=(14, 7))
    colors = {'egypt_levant': '#D4380D', 'iran_mesopotamia': '#FA8C16',
              'south_america': '#1890FF', 'western_europe': '#52C41A',
              'south_se_asia': '#722ED1'}
    markers = {'egypt_levant': 'o', 'iran_mesopotamia': 's',
               'south_america': '^', 'western_europe': 'D',
               'south_se_asia': 'v'}

    for rk in REGIONS:
        d_arr = region_d_series[rk]
        valid = ~np.isnan(d_arr)
        label = REGIONS[rk]['label']
        x = bin_centers[valid]
        y = d_arr[valid]
        if len(x) > 0:
            ax.plot(x, y, f'{markers[rk]}-', color=colors[rk], label=label,
                    linewidth=2, markersize=6, alpha=0.8)

    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.axhline(y=2, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
    ax.axhline(y=-2, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
    ax.axvspan(-3000, -2500, alpha=0.08, color='red', label='3000-2500 BCE')
    ax.set_xlabel('Calendar Year', fontsize=12)
    ax.set_ylabel('Divergence D (monument Z − settlement Z)', fontsize=12)
    ax.set_title('Cross-Continental Monument-Settlement Divergence (500-year bins)', fontsize=14)
    ax.legend(loc='upper left', fontsize=10)
    ax.set_xlim(-5200, 700)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    plot1_path = os.path.join(OUT_DIR, "cross_continental_divergence_timeseries.png")
    plt.savefig(plot1_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {plot1_path}")
    plt.close()

    # Plot 2: Correlation matrix heatmap
    n_regions = len(region_keys)
    corr_mat = np.full((n_regions, n_regions), np.nan)
    for i in range(n_regions):
        corr_mat[i, i] = 1.0
        for j in range(i+1, n_regions):
            pair_key = f"{region_keys[i]}_vs_{region_keys[j]}"
            if pair_key in correlation_matrix and correlation_matrix[pair_key]['r'] is not None:
                corr_mat[i, j] = correlation_matrix[pair_key]['r']
                corr_mat[j, i] = correlation_matrix[pair_key]['r']

    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    labels = [REGIONS[k]['label'] for k in region_keys]
    im = ax.imshow(corr_mat, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    ax.set_xticks(range(n_regions))
    ax.set_yticks(range(n_regions))
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=10)
    ax.set_yticklabels(labels, fontsize=10)
    for i in range(n_regions):
        for j in range(n_regions):
            if not np.isnan(corr_mat[i, j]):
                ax.text(j, i, f'{corr_mat[i, j]:.2f}', ha='center', va='center', fontsize=11,
                        color='white' if abs(corr_mat[i, j]) > 0.5 else 'black')
    plt.colorbar(im, ax=ax, label='Pearson r')
    ax.set_title('Cross-Continental D(t) Correlation Matrix', fontsize=13)
    plt.tight_layout()

    plot2_path = os.path.join(OUT_DIR, "cross_continental_correlation_matrix.png")
    plt.savefig(plot2_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {plot2_path}")
    plt.close()

except ImportError:
    print("  matplotlib not available — skipping plots")
    plot1_path = None
    plot2_path = None

# ============================================================
# SAVE RESULTS
# ============================================================
elapsed = time.time() - t0
print(f"\n{'=' * 70}")
print(f"COMPLETE — {elapsed:.0f}s elapsed")
print(f"{'=' * 70}")

# Convert numpy types for JSON
def clean_for_json(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj) if not np.isnan(obj) else None
    if isinstance(obj, np.ndarray):
        return [clean_for_json(x) for x in obj]
    if isinstance(obj, dict):
        return {k: clean_for_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [clean_for_json(x) for x in obj]
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj

output = clean_for_json({
    "meta": {
        "date": "2026-03-19",
        "description": "Cross-continental temporal divergence test",
        "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
        "threshold_km": THRESHOLD_KM,
        "regional_threshold_km": REGIONAL_THRESHOLD_KM,
        "n_mc_trials": N_TRIALS,
        "n_random_circles": N_RANDOM,
        "min_sites_per_bin": MIN_SITES,
        "bins": [bin_label(s, e) for s, e in BINS],
    },
    "step1_regional_temporal": regional_results,
    "step2_synchrony": {
        "peak_bins": peak_bins,
        "correlation_matrix": correlation_matrix,
    },
    "step3_south_america_deep_dive": {
        "n_sites_near_gc": len(sa_sites),
        "n_monumental": len(sa_mono),
        "n_settlement": len(sa_settle),
        "caral_era_3000_2000_bce": {
            "monumental": len(caral_era_mono),
            "settlement": len(caral_era_settle),
        },
        "sa_temporal_bins": {
            "monumental": {bin_label(s, e): sa_mono_bins.get(s, 0) for s, e in BINS},
            "settlement": {bin_label(s, e): sa_settle_bins.get(s, 0) for s, e in BINS},
        },
    },
    "step4_caral_test": caral_test,
    "step5_null_expectation": {
        "n_random_circles": N_RANDOM,
        "random_sync_count": random_sync_count,
        "random_sync_pct": round(random_sync_pct, 1),
        "random_corr_mean": round(random_corr_mean, 3),
        "random_corr_std": round(random_corr_std, 3),
        "actual_mean_corr": round(actual_mean_corr, 3),
        "actual_max_shared_peak": actual_max_shared,
        "n_random_exceeding_actual_corr": sum(1 for c in random_corrs if c >= actual_mean_corr),
    },
    "plots": {
        "timeseries": plot1_path,
        "correlation_matrix": plot2_path,
    }
})

out_path = os.path.join(OUT_DIR, "cross_continental_temporal_test.json")
with open(out_path, 'w') as f:
    json.dump(output, f, indent=2)
print(f"Saved: {out_path}")
