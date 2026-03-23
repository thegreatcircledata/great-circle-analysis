#!/usr/bin/env python3
"""
Directive C: 250-Year Temporal Resolution Analysis
====================================================
Refines the 500-year bin analysis to 250-year windows from 4000 BCE to 500 CE.

Key questions:
  - Does the onset happen at 3000 BCE, 2750 BCE, or 3250 BCE?
  - Is the collapse gradual or sharp?
  - Is the secondary peak at 500 BCE-0 CE real at finer resolution?

Pleiades: bin by minDate field into 250-year windows.
p3k14c:   bin by Age field (14C years BP -> calendar years) into same windows.

For each bin: compute monument Z, settlement Z, divergence D = monument_Z - settlement_Z.
Reliability flags: Pleiades bins with <20 sites, p3k14c bins with <30 dates.
Monte Carlo: 1000 iterations (upgraded from 200 in the 500yr version).

Output: temporal_divergence_250yr.json + plot.
"""

import csv, math, random, json, os, sys
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
N_TRIALS = 1000

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "temporal_divergence")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(os.path.join(BASE_DIR, "results"), exist_ok=True)

# 250-year bins from 4000 BCE to 500 CE
BIN_EDGES = list(range(-4000, 750, 250))  # -4000, -3750, ..., 250, 500
BINS = [(BIN_EDGES[i], BIN_EDGES[i+1]) for i in range(len(BIN_EDGES)-1)]

# Reliability thresholds
PLEIADES_MIN_SITES = 20
P3K_MIN_DATES = 30

# Pleiades type classifications (from settlement_baseline_test.py)
MONUMENTAL_TYPES = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe",
    "tumulus", "shrine",
}

SETTLEMENT_TYPES = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement",
    "townhouse", "production",
}

# p3k14c keyword classifications (from nazca_followup.py)
MONUMENTAL_KW = ['temple', 'pyramid', 'tomb', 'mound', 'monument', 'cairn', 'megalith',
                  'henge', 'huaca', 'ceremonial', 'ritual', 'geoglyph', 'ahu', 'moai',
                  'stela', 'ziggurat', 'palace', 'citadel', 'fortress', 'dolmen',
                  'stone circle', 'standing stone', 'barrow', 'tumulus', 'sanctuary',
                  'necropolis', 'acropolis', 'amphitheatre', 'colosseum']

DOMESTIC_KW = ['village', 'settlement', 'farm', 'camp', 'midden', 'shell mound',
               'quarry', 'mine', 'shelter', 'rockshelter', 'workshop', 'kiln', 'pit',
               'habitation', 'dwelling', 'house', 'homestead', 'terrace', 'field',
               'granary', 'storehouse', 'well', 'cistern', 'fishweir']


# ============================================================
# VECTORIZED HELPERS
# ============================================================
def gc_dist_vec(pole_lat, pole_lon, site_lats, site_lons):
    lat1_r, lon1_r = np.radians(pole_lat), np.radians(pole_lon)
    lat2_r, lon2_r = np.radians(site_lats), np.radians(site_lons)
    dlat = lat2_r - lat1_r
    dlon = lon2_r - lon1_r
    a = np.sin(dlat/2)**2 + np.cos(lat1_r)*np.cos(lat2_r)*np.sin(dlon/2)**2
    d = EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QUARTER_CIRC)


def rand_matched_batch(site_lats, site_lons, n):
    idx = np.random.randint(0, len(site_lats), size=n)
    lats = site_lats[idx] + np.random.normal(0, 2, n)
    lons = site_lons[idx] + np.random.normal(0, 2, n)
    return np.clip(lats, -90, 90), np.clip(lons, -180, 180)


def compute_z(pole_lat, pole_lon, lats, lons, threshold_km, n_trials):
    """Compute Z-score for a set of sites against the great circle."""
    n = len(lats)
    if n < 10:
        return {"n": n, "observed": 0, "expected": 0, "z_score": float('nan'),
                "enrichment": float('nan'), "too_few": True}

    dists = gc_dist_vec(pole_lat, pole_lon, lats, lons)
    observed = int(np.sum(dists <= threshold_km))

    rand_counts = []
    for _ in range(n_trials):
        rl, rn = rand_matched_batch(lats, lons, n)
        rd = gc_dist_vec(pole_lat, pole_lon, rl, rn)
        rand_counts.append(int(np.sum(rd <= threshold_km)))

    mu = np.mean(rand_counts)
    sigma = np.std(rand_counts)
    z = float((observed - mu) / sigma) if sigma > 0 else 0.0
    enrich = float(observed / mu) if mu > 0 else 0.0

    return {"n": n, "observed": observed, "expected": round(float(mu), 2),
            "std": round(float(sigma), 2), "z_score": round(z, 2),
            "enrichment": round(enrich, 3), "too_few": False}


def bin_label(start, end):
    """Human-readable bin label."""
    def yr(y):
        return f"{abs(y)} BCE" if y < 0 else f"{y} CE"
    return f"{yr(start)} – {yr(end)}"


def c14_to_cal(age_bp):
    """
    Rough conversion from radiocarbon years BP to calendar years BCE/CE.
    Uses a simplified IntCal20-based polynomial approximation.
    For ages < ~10000 BP this is reasonably accurate (+/-100-200 years).
    """
    if age_bp <= 0:
        return 1950  # modern
    # Piecewise linear approximation of IntCal20 cal offset
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

    cal_year = 1950 - cal_bp
    return cal_year


# ============================================================
# LOAD PLEIADES DATA
# ============================================================
print("=" * 70)
print("DIRECTIVE C: 250-YEAR TEMPORAL RESOLUTION ANALYSIS")
print("=" * 70)

print("\nLoading Pleiades data...")
pleiades_file = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")

pleiades_monument_bins = defaultdict(list)
pleiades_settlement_bins = defaultdict(list)
pleiades_skipped = 0
pleiades_total = 0
pleiades_no_date = 0

with open(pleiades_file, encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            lat = float(row['reprLat'])
            lon = float(row['reprLong'])
        except (ValueError, KeyError, TypeError):
            pleiades_skipped += 1
            continue

        if lat == 0 and lon == 0:
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue

        # Get minDate (calendar year, negative = BCE)
        try:
            min_date = float(row.get('minDate', ''))
        except (ValueError, TypeError):
            pleiades_no_date += 1
            continue

        # Classify
        feature_types = {t.strip() for t in row.get('featureTypes', '').split(',') if t.strip()}
        is_monumental = bool(feature_types & MONUMENTAL_TYPES)
        is_settlement = bool(feature_types & SETTLEMENT_TYPES)

        if not is_monumental and not is_settlement:
            continue

        pleiades_total += 1

        # Find the bin
        for start, end in BINS:
            if start <= min_date < end:
                if is_monumental:
                    pleiades_monument_bins[start].append((lat, lon))
                elif is_settlement:
                    pleiades_settlement_bins[start].append((lat, lon))
                break

print(f"  Total classified: {pleiades_total}, No date: {pleiades_no_date}, Skipped: {pleiades_skipped}")

# ============================================================
# LOAD p3k14c DATA
# ============================================================
print("\nLoading p3k14c data...")
p3k_file = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")

# Deduplicate by SiteID, keeping oldest age
site_groups = defaultdict(list)
p3k_raw = 0
with open(p3k_file, encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            lat = float(row['Lat'])
            lon = float(row['Long'])
            age_bp = float(row['Age']) if row['Age'] else None
        except (ValueError, TypeError):
            continue
        if lat == 0 and lon == 0:
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue
        if age_bp is None:
            continue
        p3k_raw += 1
        site_id = row.get('SiteID', '') or f"anon_{lat:.4f}_{lon:.4f}"
        site_groups[site_id].append({
            'lat': lat, 'lon': lon, 'age_bp': age_bp,
            'site_name': row.get('SiteName', ''),
        })

# Aggregate sites: use mean location, oldest age
p3k_sites = []
for sid, rows in site_groups.items():
    ages = [r['age_bp'] for r in rows]
    p3k_sites.append({
        'lat': np.mean([r['lat'] for r in rows]),
        'lon': np.mean([r['lon'] for r in rows]),
        'age_bp': max(ages),
        'site_name': rows[0]['site_name'],
    })

print(f"  Raw rows: {p3k_raw}, Unique sites: {len(p3k_sites)}")

# Classify p3k14c sites
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

p3k_monument_bins = defaultdict(list)
p3k_settlement_bins = defaultdict(list)
p3k_classified = 0
p3k_unclassified = 0

for site in p3k_sites:
    cls = classify_p3k(site['site_name'])
    if cls == "unclassified":
        p3k_unclassified += 1
        continue

    cal_year = c14_to_cal(site['age_bp'])
    p3k_classified += 1

    for start, end in BINS:
        if start <= cal_year < end:
            if cls == "monumental":
                p3k_monument_bins[start].append((site['lat'], site['lon']))
            else:
                p3k_settlement_bins[start].append((site['lat'], site['lon']))
            break

print(f"  Classified: {p3k_classified}, Unclassified: {p3k_unclassified}")

# ============================================================
# COMPUTE Z-SCORES PER BIN
# ============================================================
np.random.seed(42)
random.seed(42)

print("\n" + "=" * 70)
print("PLEIADES — MONUMENT vs SETTLEMENT BY 250-YEAR BIN")
print("=" * 70)

pleiades_results = []
for start, end in BINS:
    label = bin_label(start, end)
    m_sites = pleiades_monument_bins.get(start, [])
    s_sites = pleiades_settlement_bins.get(start, [])

    total_in_bin = len(m_sites) + len(s_sites)
    unreliable = total_in_bin < PLEIADES_MIN_SITES

    print(f"\n  {label}: monuments={len(m_sites)}, settlements={len(s_sites)}"
          f"{' [UNRELIABLE <20 total]' if unreliable else ''}")

    if len(m_sites) >= 10:
        m_lats = np.array([s[0] for s in m_sites])
        m_lons = np.array([s[1] for s in m_sites])
        m_res = compute_z(POLE_LAT, POLE_LON, m_lats, m_lons, THRESHOLD_KM, N_TRIALS)
        print(f"    Monument Z = {m_res['z_score']}")
    else:
        m_res = {"n": len(m_sites), "z_score": float('nan'), "too_few": True}
        print(f"    Monument: too few ({len(m_sites)})")

    if len(s_sites) >= 10:
        s_lats = np.array([s[0] for s in s_sites])
        s_lons = np.array([s[1] for s in s_sites])
        s_res = compute_z(POLE_LAT, POLE_LON, s_lats, s_lons, THRESHOLD_KM, N_TRIALS)
        print(f"    Settlement Z = {s_res['z_score']}")
    else:
        s_res = {"n": len(s_sites), "z_score": float('nan'), "too_few": True}
        print(f"    Settlement: too few ({len(s_sites)})")

    m_z = m_res['z_score'] if not m_res.get('too_few') else None
    s_z = s_res['z_score'] if not s_res.get('too_few') else None

    if m_z is not None and s_z is not None and not (np.isnan(m_z) or np.isnan(s_z)):
        D = round(m_z - s_z, 2)
    else:
        D = None

    pleiades_results.append({
        "bin_start": start,
        "bin_end": end,
        "label": label,
        "monument_n": len(m_sites),
        "settlement_n": len(s_sites),
        "monument_z": m_z,
        "settlement_z": s_z,
        "divergence_D": D,
        "unreliable": unreliable,
        "monument_detail": m_res if not m_res.get('too_few') else None,
        "settlement_detail": s_res if not s_res.get('too_few') else None,
    })

print("\n" + "=" * 70)
print("p3k14c — MONUMENT vs SETTLEMENT BY 250-YEAR BIN")
print("=" * 70)

p3k_results = []
for start, end in BINS:
    label = bin_label(start, end)
    m_sites = p3k_monument_bins.get(start, [])
    s_sites = p3k_settlement_bins.get(start, [])

    total_in_bin = len(m_sites) + len(s_sites)
    unreliable = total_in_bin < P3K_MIN_DATES

    print(f"\n  {label}: monuments={len(m_sites)}, settlements={len(s_sites)}"
          f"{' [UNRELIABLE <30 total]' if unreliable else ''}")

    if len(m_sites) >= 10:
        m_lats = np.array([s[0] for s in m_sites])
        m_lons = np.array([s[1] for s in m_sites])
        m_res = compute_z(POLE_LAT, POLE_LON, m_lats, m_lons, THRESHOLD_KM, N_TRIALS)
        print(f"    Monument Z = {m_res['z_score']}")
    else:
        m_res = {"n": len(m_sites), "z_score": float('nan'), "too_few": True}
        print(f"    Monument: too few ({len(m_sites)})")

    if len(s_sites) >= 10:
        s_lats = np.array([s[0] for s in s_sites])
        s_lons = np.array([s[1] for s in s_sites])
        s_res = compute_z(POLE_LAT, POLE_LON, s_lats, s_lons, THRESHOLD_KM, N_TRIALS)
        print(f"    Settlement Z = {s_res['z_score']}")
    else:
        s_res = {"n": len(s_sites), "z_score": float('nan'), "too_few": True}
        print(f"    Settlement: too few ({len(s_sites)})")

    m_z = m_res['z_score'] if not m_res.get('too_few') else None
    s_z = s_res['z_score'] if not s_res.get('too_few') else None

    if m_z is not None and s_z is not None and not (np.isnan(m_z) or np.isnan(s_z)):
        D = round(m_z - s_z, 2)
    else:
        D = None

    p3k_results.append({
        "bin_start": start,
        "bin_end": end,
        "label": label,
        "monument_n": len(m_sites),
        "settlement_n": len(s_sites),
        "monument_z": m_z,
        "settlement_z": s_z,
        "divergence_D": D,
        "unreliable": unreliable,
        "monument_detail": m_res if not m_res.get('too_few') else None,
        "settlement_detail": s_res if not s_res.get('too_few') else None,
    })


# ============================================================
# SUMMARY TABLE
# ============================================================
print("\n" + "=" * 70)
print("PLEIADES DIVERGENCE TABLE (250-YEAR BINS)")
print("=" * 70)
print(f"{'Bin':<28} {'Mon_N':>6} {'Set_N':>6} {'Mon_Z':>7} {'Set_Z':>7} {'D':>7} {'Flag':>10}")
print("-" * 75)
for r in pleiades_results:
    mn = r['monument_n']
    sn = r['settlement_n']
    mz = f"{r['monument_z']:.2f}" if r['monument_z'] is not None else "---"
    sz = f"{r['settlement_z']:.2f}" if r['settlement_z'] is not None else "---"
    d = f"{r['divergence_D']:.2f}" if r['divergence_D'] is not None else "---"
    flag = "UNRELIABLE" if r['unreliable'] else ""
    print(f"  {r['label']:<26} {mn:>6} {sn:>6} {mz:>7} {sz:>7} {d:>7} {flag:>10}")

print("\n" + "=" * 70)
print("p3k14c DIVERGENCE TABLE (250-YEAR BINS)")
print("=" * 70)
print(f"{'Bin':<28} {'Mon_N':>6} {'Set_N':>6} {'Mon_Z':>7} {'Set_Z':>7} {'D':>7} {'Flag':>10}")
print("-" * 75)
for r in p3k_results:
    mn = r['monument_n']
    sn = r['settlement_n']
    mz = f"{r['monument_z']:.2f}" if r['monument_z'] is not None else "---"
    sz = f"{r['settlement_z']:.2f}" if r['settlement_z'] is not None else "---"
    d = f"{r['divergence_D']:.2f}" if r['divergence_D'] is not None else "---"
    flag = "UNRELIABLE" if r['unreliable'] else ""
    print(f"  {r['label']:<26} {mn:>6} {sn:>6} {mz:>7} {sz:>7} {d:>7} {flag:>10}")

# ============================================================
# INTERPRETATION
# ============================================================
print("\n" + "=" * 70)
print("INTERPRETATION")
print("=" * 70)

pleiades_valid = [r for r in pleiades_results if r['divergence_D'] is not None]
p3k_valid = [r for r in p3k_results if r['divergence_D'] is not None]
pleiades_reliable = [r for r in pleiades_valid if not r['unreliable']]
p3k_reliable = [r for r in p3k_valid if not r['unreliable']]

interpretation = {
    "bin_width": "250 years",
    "bins_range": "4000 BCE to 500 CE",
    "n_bins": len(BINS),
    "monte_carlo_trials": N_TRIALS,
}

# --- Pleiades interpretation ---
if pleiades_valid:
    peak = max(pleiades_valid, key=lambda r: r['divergence_D'])
    interpretation['pleiades_peak_divergence'] = {
        "bin": peak['label'],
        "D": peak['divergence_D'],
        "monument_z": peak['monument_z'],
        "settlement_z": peak['settlement_z'],
        "unreliable": peak['unreliable'],
    }

    d_values = [(r['bin_start'], r['divergence_D'], r['unreliable']) for r in pleiades_valid]
    d_values.sort()

    # Q1: Onset timing — find first reliable bin where D > 2
    onset_bins_reliable = [(b, d) for b, d, u in d_values if d > 2.0 and not u]
    onset_bins_all = [(b, d) for b, d, u in d_values if d > 2.0]
    if onset_bins_reliable:
        b = onset_bins_reliable[0][0]
        interpretation['pleiades_onset'] = {
            "first_bin_D_gt_2": bin_label(b, b + 250),
            "bin_start": b,
            "D_value": onset_bins_reliable[0][1],
            "reliable": True,
        }
        print(f"\n  Pleiades onset (D > 2.0): {bin_label(b, b+250)} (D = {onset_bins_reliable[0][1]:.2f})")
    elif onset_bins_all:
        b = onset_bins_all[0][0]
        interpretation['pleiades_onset'] = {
            "first_bin_D_gt_2": bin_label(b, b + 250),
            "bin_start": b,
            "D_value": onset_bins_all[0][1],
            "reliable": False,
            "note": "bin flagged as unreliable (<20 sites)",
        }
        print(f"\n  Pleiades onset (D > 2.0): {bin_label(b, b+250)} (D = {onset_bins_all[0][1]:.2f}) [UNRELIABLE]")
    else:
        interpretation['pleiades_onset'] = "D never exceeds 2.0"
        print("\n  Pleiades: D never exceeds 2.0")

    # Q2: Collapse — is it gradual or sharp?
    # Look at D values after peak
    peak_start = peak['bin_start']
    post_peak = [(b, d) for b, d, u in d_values if b > peak_start]
    if len(post_peak) >= 2:
        # Check how many bins it takes for D to drop below 1.0
        collapse_bins = 0
        for b, d in post_peak:
            collapse_bins += 1
            if d < 1.0:
                break
        interpretation['pleiades_collapse'] = {
            "bins_to_D_lt_1": collapse_bins,
            "years_to_D_lt_1": collapse_bins * 250,
            "pattern": "sharp" if collapse_bins <= 2 else "gradual",
            "post_peak_D_values": [(bin_label(b, b+250), round(d, 2)) for b, d in post_peak],
        }
        collapse_type = "SHARP" if collapse_bins <= 2 else "GRADUAL"
        print(f"  Pleiades collapse: {collapse_type} ({collapse_bins * 250} years to D < 1.0)")

    # Q3: Secondary peak at 500 BCE-0 CE
    secondary_bins = [(b, d, u) for b, d, u in d_values if -500 <= b <= 0]
    if secondary_bins:
        max_secondary = max(secondary_bins, key=lambda x: x[1])
        interpretation['pleiades_secondary_peak'] = {
            "bins_in_range": [(bin_label(b, b+250), round(d, 2), u) for b, d, u in secondary_bins],
            "peak_bin": bin_label(max_secondary[0], max_secondary[0]+250),
            "peak_D": max_secondary[1],
            "unreliable": max_secondary[2],
            "appears_real": max_secondary[1] > 1.5 and not max_secondary[2],
        }
        status = "APPEARS REAL" if max_secondary[1] > 1.5 and not max_secondary[2] else "WEAK/UNRELIABLE"
        print(f"  Pleiades secondary peak (500 BCE-0 CE): D = {max_secondary[1]:.2f} [{status}]")

    # Chronology match at finer resolution
    pre_3000 = [d for b, d, u in d_values if b < -3000 and not u]
    peak_period = [d for b, d, u in d_values if -3250 <= b <= -2000 and not u]
    post_2000 = [d for b, d, u in d_values if b > -2000 and not u]

    if pre_3000 and peak_period:
        avg_pre = sum(pre_3000) / len(pre_3000)
        avg_peak = sum(peak_period) / len(peak_period)
        avg_post = sum(post_2000) / len(post_2000) if post_2000 else 0
        interpretation['pleiades_chronology_match'] = {
            "avg_D_pre_3000BCE": round(avg_pre, 2),
            "avg_D_3250_2000BCE": round(avg_peak, 2),
            "avg_D_post_2000BCE": round(avg_post, 2),
            "pattern": "rise_and_fall" if avg_peak > avg_pre and avg_peak > avg_post else
                       "monotonic_rise" if avg_peak > avg_pre else
                       "no_clear_pattern",
            "matches_known_chronology": avg_peak > avg_pre,
        }

# --- p3k14c interpretation ---
if p3k_valid:
    peak = max(p3k_valid, key=lambda r: r['divergence_D'])
    interpretation['p3k14c_peak_divergence'] = {
        "bin": peak['label'],
        "D": peak['divergence_D'],
        "monument_z": peak['monument_z'],
        "settlement_z": peak['settlement_z'],
        "unreliable": peak['unreliable'],
    }

    d_values_p3k = [(r['bin_start'], r['divergence_D'], r['unreliable']) for r in p3k_valid]
    d_values_p3k.sort()

    onset_bins_p3k = [(b, d) for b, d, u in d_values_p3k if d > 2.0 and not u]
    if onset_bins_p3k:
        b = onset_bins_p3k[0][0]
        interpretation['p3k14c_onset'] = {
            "first_bin_D_gt_2": bin_label(b, b + 250),
            "bin_start": b,
            "D_value": onset_bins_p3k[0][1],
        }
        print(f"\n  p3k14c onset (D > 2.0): {bin_label(b, b+250)} (D = {onset_bins_p3k[0][1]:.2f})")

    # Secondary peak
    secondary_bins_p3k = [(b, d, u) for b, d, u in d_values_p3k if -500 <= b <= 0]
    if secondary_bins_p3k:
        max_sec = max(secondary_bins_p3k, key=lambda x: x[1])
        interpretation['p3k14c_secondary_peak'] = {
            "bins_in_range": [(bin_label(b, b+250), round(d, 2), u) for b, d, u in secondary_bins_p3k],
            "peak_bin": bin_label(max_sec[0], max_sec[0]+250),
            "peak_D": max_sec[1],
            "unreliable": max_sec[2],
        }

# Summary of key questions
key_questions = {}
if 'pleiades_onset' in interpretation and isinstance(interpretation['pleiades_onset'], dict):
    key_questions['onset_timing'] = interpretation['pleiades_onset'].get('first_bin_D_gt_2', 'unknown')
if 'pleiades_collapse' in interpretation:
    key_questions['collapse_pattern'] = interpretation['pleiades_collapse'].get('pattern', 'unknown')
if 'pleiades_secondary_peak' in interpretation:
    sp = interpretation['pleiades_secondary_peak']
    key_questions['secondary_peak_real'] = sp.get('appears_real', False)
    key_questions['secondary_peak_D'] = sp.get('peak_D')

interpretation['key_questions'] = key_questions

print("\n  Key questions summary:")
for k, v in key_questions.items():
    print(f"    {k}: {v}")

# ============================================================
# RELIABILITY SUMMARY
# ============================================================
pleiades_unreliable_bins = [r['label'] for r in pleiades_results if r['unreliable']]
p3k_unreliable_bins = [r['label'] for r in p3k_results if r['unreliable']]

interpretation['reliability'] = {
    "pleiades_unreliable_bins": pleiades_unreliable_bins,
    "pleiades_unreliable_count": len(pleiades_unreliable_bins),
    "p3k14c_unreliable_bins": p3k_unreliable_bins,
    "p3k14c_unreliable_count": len(p3k_unreliable_bins),
    "pleiades_threshold": f"<{PLEIADES_MIN_SITES} total sites in bin",
    "p3k14c_threshold": f"<{P3K_MIN_DATES} total dates in bin",
}

print(f"\n  Pleiades unreliable bins: {len(pleiades_unreliable_bins)} of {len(BINS)}")
print(f"  p3k14c unreliable bins:  {len(p3k_unreliable_bins)} of {len(BINS)}")

# ============================================================
# SAVE RESULTS
# ============================================================
output = {
    "meta": {
        "analysis": "Directive C: 250-Year Temporal Resolution Analysis",
        "date": "2026-03-19",
        "methodology": {
            "bins": "250-year windows from 4000 BCE to 500 CE",
            "n_bins": len(BINS),
            "pleiades_date_field": "minDate (calendar year)",
            "p3k14c_date_field": "Age (14C years BP, converted to cal years via simplified IntCal20 approximation)",
            "p3k14c_classification": "keyword-based (monumental vs domestic from site names)",
            "pleiades_classification": "featureTypes field (same types as settlement_baseline_test)",
            "threshold_km": THRESHOLD_KM,
            "monte_carlo_trials": N_TRIALS,
            "pleiades_reliability_threshold": PLEIADES_MIN_SITES,
            "p3k14c_reliability_threshold": P3K_MIN_DATES,
            "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
        },
        "key_questions": [
            "Does the onset happen at 3000 BCE or 2750 BCE or 3250 BCE?",
            "Is the collapse gradual or sharp?",
            "Is the secondary peak at 500 BCE-0 CE real at finer resolution?",
        ],
    },
    "pleiades_bins": pleiades_results,
    "p3k14c_bins": p3k_results,
    "interpretation": interpretation,
}

out_path = os.path.join(OUT_DIR, "temporal_divergence_250yr.json")
with open(out_path, "w") as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nSaved: {out_path}")

results_path = os.path.join(BASE_DIR, "results", "temporal_divergence_250yr.json")
with open(results_path, "w") as f:
    json.dump(output, f, indent=2, default=str)
print(f"Saved: {results_path}")

# ============================================================
# GENERATE PLOT
# ============================================================
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(16, 11), sharex=True)

    # --- Pleiades subplot ---
    ax1 = axes[0]
    p_bins_valid = [r for r in pleiades_results if r['divergence_D'] is not None]
    if p_bins_valid:
        x = [(r['bin_start'] + r['bin_end']) / 2 for r in p_bins_valid]
        d_vals = [r['divergence_D'] for r in p_bins_valid]
        m_vals = [r['monument_z'] for r in p_bins_valid]
        s_vals = [r['settlement_z'] for r in p_bins_valid]
        unreliable_mask = [r['unreliable'] for r in p_bins_valid]

        # Bar chart with unreliable bins in a different color
        bar_colors = ['#d4a0a0' if u else '#9467bd' for u in unreliable_mask]
        ax1.bar(x, d_vals, width=200, alpha=0.3, color=bar_colors)
        ax1.plot(x, d_vals, 'o-', color='purple', linewidth=2, markersize=7, label='D (Mon Z - Set Z)')
        ax1.plot(x, m_vals, 's--', color='red', linewidth=1.5, markersize=5, alpha=0.7, label='Monument Z')
        ax1.plot(x, s_vals, '^--', color='blue', linewidth=1.5, markersize=5, alpha=0.7, label='Settlement Z')
        ax1.axhline(y=0, color='black', linewidth=0.5, linestyle='-')
        ax1.axhline(y=2, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
        ax1.axhline(y=-2, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)

        # Shade known monument construction era
        ax1.axvspan(-3000, -2000, alpha=0.08, color='red', label='Peak monument era')

        # Mark unreliable bins
        unreliable_x = [xi for xi, u in zip(x, unreliable_mask) if u]
        unreliable_d = [di for di, u in zip(d_vals, unreliable_mask) if u]
        if unreliable_x:
            ax1.scatter(unreliable_x, unreliable_d, marker='x', color='gray', s=100,
                       zorder=5, label='Unreliable (<20 sites)')

    ax1.set_ylabel('Z-score / Divergence D', fontsize=12)
    ax1.set_title('PLEIADES — Monument-Settlement Divergence by 250-Year Bin (Directive C)',
                   fontsize=13, fontweight='bold')
    ax1.legend(loc='upper left', fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Site count annotations
    for r in p_bins_valid:
        mid = (r['bin_start'] + r['bin_end']) / 2
        ax1.annotate(f"m={r['monument_n']}\ns={r['settlement_n']}",
                     xy=(mid, min(r['divergence_D'], 0) - 0.5),
                     fontsize=5, ha='center', alpha=0.5)

    # --- p3k14c subplot ---
    ax2 = axes[1]
    p3k_bins_valid = [r for r in p3k_results if r['divergence_D'] is not None]
    if p3k_bins_valid:
        x = [(r['bin_start'] + r['bin_end']) / 2 for r in p3k_bins_valid]
        d_vals = [r['divergence_D'] for r in p3k_bins_valid]
        m_vals = [r['monument_z'] for r in p3k_bins_valid]
        s_vals = [r['settlement_z'] for r in p3k_bins_valid]
        unreliable_mask = [r['unreliable'] for r in p3k_bins_valid]

        bar_colors = ['#d4a0a0' if u else '#9467bd' for u in unreliable_mask]
        ax2.bar(x, d_vals, width=200, alpha=0.3, color=bar_colors)
        ax2.plot(x, d_vals, 'o-', color='purple', linewidth=2, markersize=7, label='D (Mon Z - Set Z)')
        ax2.plot(x, m_vals, 's--', color='red', linewidth=1.5, markersize=5, alpha=0.7, label='Monument Z')
        ax2.plot(x, s_vals, '^--', color='blue', linewidth=1.5, markersize=5, alpha=0.7, label='Settlement Z')
        ax2.axhline(y=0, color='black', linewidth=0.5, linestyle='-')
        ax2.axhline(y=2, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
        ax2.axhline(y=-2, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)

        ax2.axvspan(-3000, -2000, alpha=0.08, color='red', label='Peak monument era')

        unreliable_x = [xi for xi, u in zip(x, unreliable_mask) if u]
        unreliable_d = [di for di, u in zip(d_vals, unreliable_mask) if u]
        if unreliable_x:
            ax2.scatter(unreliable_x, unreliable_d, marker='x', color='gray', s=100,
                       zorder=5, label='Unreliable (<30 dates)')

    ax2.set_xlabel('Calendar Year (negative = BCE)', fontsize=12)
    ax2.set_ylabel('Z-score / Divergence D', fontsize=12)
    ax2.set_title('p3k14c — Monument-Settlement Divergence by 250-Year Bin (Directive C)',
                   fontsize=13, fontweight='bold')
    ax2.legend(loc='upper left', fontsize=8)
    ax2.grid(True, alpha=0.3)

    for r in p3k_bins_valid:
        mid = (r['bin_start'] + r['bin_end']) / 2
        ax2.annotate(f"m={r['monument_n']}\ns={r['settlement_n']}",
                     xy=(mid, min(r['divergence_D'], 0) - 0.5),
                     fontsize=5, ha='center', alpha=0.5)

    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, "temporal_divergence_250yr_plot.png")
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {plot_path}")

    plot_path2 = os.path.join(BASE_DIR, "results", "temporal_divergence_250yr_plot.png")
    plt.savefig(plot_path2, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {plot_path2}")

except ImportError as e:
    print(f"\nCould not generate plot: {e}")
    print("Install matplotlib: pip install matplotlib")

print("\nDONE.")
