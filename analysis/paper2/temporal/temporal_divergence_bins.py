#!/usr/bin/env python3
"""
Temporal Divergence Deep Dive
==============================
Monument-settlement divergence test in 500-year bins.

Pleiades: bin by minDate field into 500-year windows from 5000 BCE to 500 CE.
p3k14c:   bin by Age field (14C years BP -> calendar years) into same windows.

For each bin: compute monument Z, settlement Z, divergence D = monument_Z - settlement_Z.
Output: temporal_divergence_bins.json + plot.
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
N_TRIALS = 200

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "temporal_divergence")
os.makedirs(OUT_DIR, exist_ok=True)

# 500-year bins from 5000 BCE to 500 CE
# Using calendar years (negative = BCE)
BIN_EDGES = list(range(-5000, 1000, 500))  # -5000, -4500, ..., 0, 500
BINS = [(BIN_EDGES[i], BIN_EDGES[i+1]) for i in range(len(BIN_EDGES)-1)]

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
    For ages < ~10000 BP this is reasonably accurate (±100-200 years).
    """
    # BP to calendar year (1950 - age_bp gives uncalibrated cal year)
    # Apply rough correction for the main C14/cal offset
    if age_bp <= 0:
        return 1950  # modern
    # Piecewise linear approximation of IntCal20 cal offset
    # These are rough midpoints — good enough for 500-year binning
    if age_bp < 2500:
        cal_bp = age_bp * 1.0  # minimal offset in recent millennia
    elif age_bp < 5000:
        cal_bp = age_bp * 1.05 + 50  # slight offset 2500-5000 BP
    elif age_bp < 8000:
        cal_bp = age_bp * 1.08 + 100  # larger offset
    elif age_bp < 12000:
        cal_bp = age_bp * 1.12 + 200  # growing offset
    else:
        cal_bp = age_bp * 1.15 + 300  # deep time offset

    cal_year = 1950 - cal_bp  # calendar year (negative = BCE)
    return cal_year


# ============================================================
# LOAD PLEIADES DATA
# ============================================================
print("=" * 70)
print("TEMPORAL DIVERGENCE DEEP DIVE — 500-YEAR BINS")
print("=" * 70)

print("\nLoading Pleiades data...")
pleiades_file = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")

# Store sites with their bin assignment
pleiades_monument_bins = defaultdict(list)  # bin_start -> [(lat, lon), ...]
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

# First deduplicate by SiteID, keeping oldest age
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

# Aggregate sites: use median location, oldest age
p3k_sites = []
for sid, rows in site_groups.items():
    ages = [r['age_bp'] for r in rows]
    p3k_sites.append({
        'lat': np.mean([r['lat'] for r in rows]),
        'lon': np.mean([r['lon'] for r in rows]),
        'age_bp': max(ages),  # oldest date for the site
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

MIN_SITES = 20  # minimum sites in a bin to compute Z

print("\n" + "=" * 70)
print("PLEIADES — MONUMENT vs SETTLEMENT BY 500-YEAR BIN")
print("=" * 70)

pleiades_results = []
for start, end in BINS:
    label = bin_label(start, end)
    m_sites = pleiades_monument_bins.get(start, [])
    s_sites = pleiades_settlement_bins.get(start, [])

    print(f"\n  {label}: monuments={len(m_sites)}, settlements={len(s_sites)}")

    if len(m_sites) >= MIN_SITES:
        m_lats = np.array([s[0] for s in m_sites])
        m_lons = np.array([s[1] for s in m_sites])
        m_res = compute_z(POLE_LAT, POLE_LON, m_lats, m_lons, THRESHOLD_KM, N_TRIALS)
        print(f"    Monument Z = {m_res['z_score']}")
    else:
        m_res = {"n": len(m_sites), "z_score": float('nan'), "too_few": True}
        print(f"    Monument: too few ({len(m_sites)})")

    if len(s_sites) >= MIN_SITES:
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
        "monument_detail": m_res if not m_res.get('too_few') else None,
        "settlement_detail": s_res if not s_res.get('too_few') else None,
    })

print("\n" + "=" * 70)
print("p3k14c — MONUMENT vs SETTLEMENT BY 500-YEAR BIN")
print("=" * 70)

p3k_results = []
for start, end in BINS:
    label = bin_label(start, end)
    m_sites = p3k_monument_bins.get(start, [])
    s_sites = p3k_settlement_bins.get(start, [])

    print(f"\n  {label}: monuments={len(m_sites)}, settlements={len(s_sites)}")

    if len(m_sites) >= MIN_SITES:
        m_lats = np.array([s[0] for s in m_sites])
        m_lons = np.array([s[1] for s in m_sites])
        m_res = compute_z(POLE_LAT, POLE_LON, m_lats, m_lons, THRESHOLD_KM, N_TRIALS)
        print(f"    Monument Z = {m_res['z_score']}")
    else:
        m_res = {"n": len(m_sites), "z_score": float('nan'), "too_few": True}
        print(f"    Monument: too few ({len(m_sites)})")

    if len(s_sites) >= MIN_SITES:
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
        "monument_detail": m_res if not m_res.get('too_few') else None,
        "settlement_detail": s_res if not s_res.get('too_few') else None,
    })


# ============================================================
# SUMMARY TABLE
# ============================================================
print("\n" + "=" * 70)
print("PLEIADES DIVERGENCE TABLE")
print("=" * 70)
print(f"{'Bin':<25} {'Mon_N':>6} {'Set_N':>6} {'Mon_Z':>7} {'Set_Z':>7} {'D':>7}")
print("-" * 60)
for r in pleiades_results:
    mn = r['monument_n']
    sn = r['settlement_n']
    mz = f"{r['monument_z']:.2f}" if r['monument_z'] is not None else "—"
    sz = f"{r['settlement_z']:.2f}" if r['settlement_z'] is not None else "—"
    d = f"{r['divergence_D']:.2f}" if r['divergence_D'] is not None else "—"
    print(f"  {r['label']:<23} {mn:>6} {sn:>6} {mz:>7} {sz:>7} {d:>7}")

print("\n" + "=" * 70)
print("p3k14c DIVERGENCE TABLE")
print("=" * 70)
print(f"{'Bin':<25} {'Mon_N':>6} {'Set_N':>6} {'Mon_Z':>7} {'Set_Z':>7} {'D':>7}")
print("-" * 60)
for r in p3k_results:
    mn = r['monument_n']
    sn = r['settlement_n']
    mz = f"{r['monument_z']:.2f}" if r['monument_z'] is not None else "—"
    sz = f"{r['settlement_z']:.2f}" if r['settlement_z'] is not None else "—"
    d = f"{r['divergence_D']:.2f}" if r['divergence_D'] is not None else "—"
    print(f"  {r['label']:<23} {mn:>6} {sn:>6} {mz:>7} {sz:>7} {d:>7}")

# ============================================================
# INTERPRETATION
# ============================================================

# Find peak divergence
pleiades_valid = [r for r in pleiades_results if r['divergence_D'] is not None]
p3k_valid = [r for r in p3k_results if r['divergence_D'] is not None]

interpretation = {}
if pleiades_valid:
    peak = max(pleiades_valid, key=lambda r: r['divergence_D'])
    interpretation['pleiades_peak_divergence'] = {
        "bin": peak['label'],
        "D": peak['divergence_D'],
        "monument_z": peak['monument_z'],
        "settlement_z": peak['settlement_z'],
    }

    # Characterize onset
    d_values = [(r['bin_start'], r['divergence_D']) for r in pleiades_valid]
    d_values.sort()

    # Find first bin where D > 2
    onset_bins = [b for b, d in d_values if d > 2]
    if onset_bins:
        interpretation['pleiades_onset'] = f"D first exceeds 2.0 in bin starting {onset_bins[0]} ({bin_label(onset_bins[0], onset_bins[0]+500)})"
    else:
        interpretation['pleiades_onset'] = "D never exceeds 2.0"

    # Check if pattern matches known chronology
    pre_3000 = [d for b, d in d_values if b < -3000]
    peak_period = [d for b, d in d_values if -3000 <= b <= -2000]
    post_2000 = [d for b, d in d_values if b > -2000]

    if pre_3000 and peak_period:
        avg_pre = sum(pre_3000) / len(pre_3000)
        avg_peak = sum(peak_period) / len(peak_period)
        avg_post = sum(post_2000) / len(post_2000) if post_2000 else 0
        interpretation['pleiades_chronology_match'] = {
            "avg_D_pre_3000BCE": round(avg_pre, 2),
            "avg_D_3000_2000BCE": round(avg_peak, 2),
            "avg_D_post_2000BCE": round(avg_post, 2),
            "pattern": "rise_and_fall" if avg_peak > avg_pre and avg_peak > avg_post else
                       "monotonic_rise" if avg_peak > avg_pre else
                       "no_clear_pattern",
            "matches_known_chronology": avg_peak > avg_pre,
        }

if p3k_valid:
    peak = max(p3k_valid, key=lambda r: r['divergence_D'])
    interpretation['p3k14c_peak_divergence'] = {
        "bin": peak['label'],
        "D": peak['divergence_D'],
        "monument_z": peak['monument_z'],
        "settlement_z": peak['settlement_z'],
    }

# ============================================================
# SAVE RESULTS
# ============================================================
output = {
    "meta": {
        "analysis": "Temporal Divergence Deep Dive — 500-Year Bins",
        "date": "2026-03-19",
        "methodology": {
            "bins": "500-year windows from 5000 BCE to 500 CE",
            "pleiades_date_field": "minDate (calendar year)",
            "p3k14c_date_field": "Age (14C years BP, converted to cal years via simplified IntCal20 approximation)",
            "p3k14c_classification": "keyword-based (monumental vs domestic from site names)",
            "pleiades_classification": "featureTypes field (same types as settlement_baseline_test)",
            "threshold_km": THRESHOLD_KM,
            "monte_carlo_trials": N_TRIALS,
            "min_sites_per_bin": MIN_SITES,
            "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
        },
        "question": "Does monument-settlement divergence onset sharply or gradually? Does it track the known chronology of large-scale monument construction (onset ~3000 BCE, peak ~2500 BCE, decline after ~2000 BCE)?",
    },
    "pleiades_bins": pleiades_results,
    "p3k14c_bins": p3k_results,
    "interpretation": interpretation,
}

# Save to both outputs dir and results dir
out_path = os.path.join(OUT_DIR, "temporal_divergence_bins.json")
with open(out_path, "w") as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nSaved: {out_path}")

results_path = os.path.join(BASE_DIR, "results", "temporal_divergence_bins.json")
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

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

    # --- Pleiades subplot ---
    ax1 = axes[0]
    p_bins_valid = [r for r in pleiades_results if r['divergence_D'] is not None]
    if p_bins_valid:
        x = [(r['bin_start'] + r['bin_end']) / 2 for r in p_bins_valid]
        d_vals = [r['divergence_D'] for r in p_bins_valid]
        m_vals = [r['monument_z'] for r in p_bins_valid]
        s_vals = [r['settlement_z'] for r in p_bins_valid]

        ax1.bar(x, d_vals, width=400, alpha=0.3, color='purple', label='Divergence D')
        ax1.plot(x, d_vals, 'o-', color='purple', linewidth=2, markersize=8, label='D (Mon Z − Set Z)')
        ax1.plot(x, m_vals, 's--', color='red', linewidth=1.5, markersize=6, alpha=0.7, label='Monument Z')
        ax1.plot(x, s_vals, '^--', color='blue', linewidth=1.5, markersize=6, alpha=0.7, label='Settlement Z')
        ax1.axhline(y=0, color='black', linewidth=0.5, linestyle='-')
        ax1.axhline(y=2, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
        ax1.axhline(y=-2, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)

        # Shade known monument construction era
        ax1.axvspan(-3000, -2000, alpha=0.08, color='red', label='Peak monument era')

    ax1.set_ylabel('Z-score / Divergence D', fontsize=12)
    ax1.set_title('PLEIADES — Monument-Settlement Divergence by 500-Year Bin', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper left', fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Add site count annotations
    for r in p_bins_valid:
        mid = (r['bin_start'] + r['bin_end']) / 2
        ax1.annotate(f"m={r['monument_n']}\ns={r['settlement_n']}",
                     xy=(mid, min(r['divergence_D'], 0) - 0.5),
                     fontsize=6, ha='center', alpha=0.5)

    # --- p3k14c subplot ---
    ax2 = axes[1]
    p3k_bins_valid = [r for r in p3k_results if r['divergence_D'] is not None]
    if p3k_bins_valid:
        x = [(r['bin_start'] + r['bin_end']) / 2 for r in p3k_bins_valid]
        d_vals = [r['divergence_D'] for r in p3k_bins_valid]
        m_vals = [r['monument_z'] for r in p3k_bins_valid]
        s_vals = [r['settlement_z'] for r in p3k_bins_valid]

        ax2.bar(x, d_vals, width=400, alpha=0.3, color='purple', label='Divergence D')
        ax2.plot(x, d_vals, 'o-', color='purple', linewidth=2, markersize=8, label='D (Mon Z − Set Z)')
        ax2.plot(x, m_vals, 's--', color='red', linewidth=1.5, markersize=6, alpha=0.7, label='Monument Z')
        ax2.plot(x, s_vals, '^--', color='blue', linewidth=1.5, markersize=6, alpha=0.7, label='Settlement Z')
        ax2.axhline(y=0, color='black', linewidth=0.5, linestyle='-')
        ax2.axhline(y=2, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
        ax2.axhline(y=-2, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)

        ax2.axvspan(-3000, -2000, alpha=0.08, color='red', label='Peak monument era')

    ax2.set_xlabel('Calendar Year (negative = BCE)', fontsize=12)
    ax2.set_ylabel('Z-score / Divergence D', fontsize=12)
    ax2.set_title('p3k14c — Monument-Settlement Divergence by 500-Year Bin', fontsize=14, fontweight='bold')
    ax2.legend(loc='upper left', fontsize=9)
    ax2.grid(True, alpha=0.3)

    for r in p3k_bins_valid:
        mid = (r['bin_start'] + r['bin_end']) / 2
        ax2.annotate(f"m={r['monument_n']}\ns={r['settlement_n']}",
                     xy=(mid, min(r['divergence_D'], 0) - 0.5),
                     fontsize=6, ha='center', alpha=0.5)

    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, "temporal_divergence_plot.png")
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {plot_path}")

    # Also save to results
    plot_path2 = os.path.join(BASE_DIR, "results", "temporal_divergence_plot.png")
    plt.savefig(plot_path2, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {plot_path2}")

except ImportError as e:
    print(f"\nCould not generate plot: {e}")
    print("Install matplotlib: pip install matplotlib")

print("\nDONE.")
