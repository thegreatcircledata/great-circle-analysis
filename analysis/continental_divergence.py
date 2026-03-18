#!/usr/bin/env python3
"""
Continental Divergence Directive
=================================
Determines whether the p3k14c monument-settlement divergence survives
without African/Egyptian sites. Five studies addressing the criticism
that the divergence is an artifact of Egyptian recording practices.

Studies:
  1. p3k14c Divergence by Continent
  2. p3k14c Divergence Without Africa
  3. p3k14c Divergence Within South America Only
  4. p3k14c Egypt-Only Divergence
  5. Keyword Classification Audit
"""

import csv, math, random, json, os, time
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
OUT_DIR = os.path.join(BASE_DIR, "outputs", "continental_divergence")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# KEYWORDS (same as nazca_followup.py Study 2)
# ============================================================
MONUMENTAL_KW = ['temple', 'pyramid', 'tomb', 'mound', 'monument', 'cairn', 'megalith',
                  'henge', 'huaca', 'ceremonial', 'ritual', 'geoglyph', 'ahu', 'moai',
                  'stela', 'ziggurat', 'palace', 'citadel', 'fortress', 'dolmen',
                  'stone circle', 'standing stone', 'barrow', 'tumulus', 'sanctuary',
                  'necropolis', 'acropolis', 'amphitheatre', 'colosseum']

DOMESTIC_KW = ['village', 'settlement', 'farm', 'camp', 'midden', 'shell mound',
               'quarry', 'mine', 'shelter', 'rockshelter', 'workshop', 'kiln', 'pit',
               'habitation', 'dwelling', 'house', 'homestead', 'terrace', 'field',
               'granary', 'storehouse', 'well', 'cistern', 'fishweir']

def classify_p3k_site(site_name):
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

# ============================================================
# VECTORIZED HELPERS
# ============================================================
def haversine_vec(lat1, lon1, lat2_arr, lon2_arr):
    lat1_r, lon1_r = np.radians(lat1), np.radians(lon1)
    lat2_r, lon2_r = np.radians(lat2_arr), np.radians(lon2_arr)
    dlat = lat2_r - lat1_r; dlon = lon2_r - lon1_r
    a = np.sin(dlat/2)**2 + np.cos(lat1_r)*np.cos(lat2_r)*np.sin(dlon/2)**2
    return EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))

def gc_dist_vec(pole_lat, pole_lon, site_lats, site_lons):
    d = haversine_vec(pole_lat, pole_lon, site_lats, site_lons)
    return np.abs(d - QUARTER_CIRC)

def rand_matched_batch(site_lats, site_lons, n):
    idx = np.random.randint(0, len(site_lats), size=n)
    lats = site_lats[idx] + np.random.normal(0, 2, n)
    lons = site_lons[idx] + np.random.normal(0, 2, n)
    return np.clip(lats, -90, 90), np.clip(lons, -180, 180)

def run_mc(site_lats, site_lons, threshold=THRESHOLD_KM, n_trials=N_TRIALS):
    """Run distribution-matched MC and return results dict."""
    n = len(site_lats)
    if n < 5:
        return {"n_sites": int(n), "observed": 0, "expected": 0, "std": 0,
                "z_score": 0, "enrichment": 0, "p_value": 1.0, "too_few": True}
    dists = gc_dist_vec(POLE_LAT, POLE_LON, site_lats, site_lons)
    observed = int(np.sum(dists <= threshold))
    rand_counts = []
    for _ in range(n_trials):
        rl, rn = rand_matched_batch(site_lats, site_lons, n)
        rd = gc_dist_vec(POLE_LAT, POLE_LON, rl, rn)
        rand_counts.append(int(np.sum(rd <= threshold)))
    mu = np.mean(rand_counts); sigma = np.std(rand_counts)
    z = float((observed - mu) / sigma) if sigma > 0 else 0.0
    enrich = float(observed / mu) if mu > 0 else 0.0
    p_val = float(np.sum(np.array(rand_counts) >= observed) / n_trials)
    return {
        "n_sites": int(n), "observed": int(observed),
        "expected": round(float(mu), 2), "std": round(float(sigma), 2),
        "z_score": round(z, 2), "enrichment": round(enrich, 3),
        "p_value": round(p_val, 4)
    }

def save_json(data, filename):
    path = os.path.join(OUT_DIR, filename)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  -> Saved: {path}")

# ============================================================
# LOAD DATA
# ============================================================
print("=" * 70)
print("CONTINENTAL DIVERGENCE DIRECTIVE")
print("=" * 70)

t0 = time.time()
print("\nLoading p3k14c data...")

p3k_all_rows = []
with open(os.path.join(BASE_DIR, "p3k14c_data.csv"), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            lat, lon = float(row['Lat']), float(row['Long'])
            if lat == 0 and lon == 0: continue
            if not (-90 <= lat <= 90 and -180 <= lon <= 180): continue
            p3k_all_rows.append({
                'lat': lat, 'lon': lon,
                'site_id': row.get('SiteID', ''),
                'site_name': row.get('SiteName', ''),
                'country': row.get('Country', ''),
                'continent': row.get('Continent', ''),
            })
        except: pass

# Deduplicate by SiteID (average coords for multi-date sites)
site_groups = defaultdict(list)
for r in p3k_all_rows:
    key = r['site_id'] if r['site_id'] else f"anon_{r['lat']:.4f}_{r['lon']:.4f}"
    site_groups[key].append(r)

p3k_sites = []
for sid, rows in site_groups.items():
    p3k_sites.append({
        'lat': np.mean([r['lat'] for r in rows]),
        'lon': np.mean([r['lon'] for r in rows]),
        'site_name': rows[0]['site_name'],
        'country': rows[0]['country'],
        'continent': rows[0]['continent'],
        'site_id': sid,
        'n_dates': len(rows),
        'classification': classify_p3k_site(rows[0]['site_name']),
    })

p3k_lats = np.array([s['lat'] for s in p3k_sites], dtype=np.float64)
p3k_lons = np.array([s['lon'] for s in p3k_sites], dtype=np.float64)

class_counts = Counter(s['classification'] for s in p3k_sites)
print(f"  Total unique sites: {len(p3k_sites)}")
print(f"  Classification: {dict(class_counts)}")

# ============================================================
# STUDY 1: p3k14c DIVERGENCE BY CONTINENT
# ============================================================
print("\n" + "=" * 70)
print("STUDY 1: p3k14c DIVERGENCE BY CONTINENT")
print("=" * 70)

np.random.seed(42)

continents_of_interest = ["Africa", "Asia", "South America"]
continent_results = {}

for cont in continents_of_interest:
    mask = np.array([s['continent'] == cont for s in p3k_sites])
    cont_sites = [s for s in p3k_sites if s['continent'] == cont]

    mono_mask = np.array([s['continent'] == cont and s['classification'] == 'monumental' for s in p3k_sites])
    dome_mask = np.array([s['continent'] == cont and s['classification'] == 'domestic' for s in p3k_sites])

    n_mono = int(np.sum(mono_mask))
    n_dome = int(np.sum(dome_mask))
    print(f"\n  --- {cont} (total={int(np.sum(mask))}, monumental={n_mono}, domestic={n_dome}) ---")

    if n_mono >= 5:
        mono_res = run_mc(p3k_lats[mono_mask], p3k_lons[mono_mask])
        print(f"  Monumental: Z={mono_res['z_score']}, obs={mono_res['observed']}")
    else:
        mono_res = {"n_sites": n_mono, "z_score": 0, "too_few": True}
        print(f"  Monumental: too few ({n_mono})")

    if n_dome >= 5:
        dome_res = run_mc(p3k_lats[dome_mask], p3k_lons[dome_mask])
        print(f"  Domestic: Z={dome_res['z_score']}, obs={dome_res['observed']}")
    else:
        dome_res = {"n_sites": n_dome, "z_score": 0, "too_few": True}
        print(f"  Domestic: too few ({n_dome})")

    if not mono_res.get("too_few") and not dome_res.get("too_few"):
        divergence = round(mono_res['z_score'] - dome_res['z_score'], 2)
    else:
        divergence = "N/A"

    continent_results[cont] = {
        "monumental": mono_res,
        "domestic": dome_res,
        "divergence": divergence,
    }
    print(f"  Divergence: {divergence}")

# "All without Africa"
print(f"\n  --- All without Africa ---")
no_africa_mono = np.array([s['continent'] != 'Africa' and s['classification'] == 'monumental' for s in p3k_sites])
no_africa_dome = np.array([s['continent'] != 'Africa' and s['classification'] == 'domestic' for s in p3k_sites])
n_mono_noaf = int(np.sum(no_africa_mono)); n_dome_noaf = int(np.sum(no_africa_dome))
print(f"  Monumental (no Africa): {n_mono_noaf}, Domestic (no Africa): {n_dome_noaf}")

mono_noaf_res = run_mc(p3k_lats[no_africa_mono], p3k_lons[no_africa_mono])
dome_noaf_res = run_mc(p3k_lats[no_africa_dome], p3k_lons[no_africa_dome])
div_noaf = round(mono_noaf_res['z_score'] - dome_noaf_res['z_score'], 2)
print(f"  Monumental Z={mono_noaf_res['z_score']}, Domestic Z={dome_noaf_res['z_score']}, Divergence={div_noaf}")

continent_results["all_without_africa"] = {
    "monumental": mono_noaf_res, "domestic": dome_noaf_res, "divergence": div_noaf,
}

# "All without South America"
print(f"\n  --- All without South America ---")
no_sa_mono = np.array([s['continent'] != 'South America' and s['classification'] == 'monumental' for s in p3k_sites])
no_sa_dome = np.array([s['continent'] != 'South America' and s['classification'] == 'domestic' for s in p3k_sites])
n_mono_nosa = int(np.sum(no_sa_mono)); n_dome_nosa = int(np.sum(no_sa_dome))
print(f"  Monumental (no SA): {n_mono_nosa}, Domestic (no SA): {n_dome_nosa}")

mono_nosa_res = run_mc(p3k_lats[no_sa_mono], p3k_lons[no_sa_mono])
dome_nosa_res = run_mc(p3k_lats[no_sa_dome], p3k_lons[no_sa_dome])
div_nosa = round(mono_nosa_res['z_score'] - dome_nosa_res['z_score'], 2)
print(f"  Monumental Z={mono_nosa_res['z_score']}, Domestic Z={dome_nosa_res['z_score']}, Divergence={div_nosa}")

continent_results["all_without_south_america"] = {
    "monumental": mono_nosa_res, "domestic": dome_nosa_res, "divergence": div_nosa,
}

save_json({
    "meta": {"date": "2026-03-18", "dataset": "p3k14c", "threshold_km": THRESHOLD_KM,
             "n_trials": N_TRIALS, "classification_counts": dict(class_counts)},
    "continents": continent_results,
}, "continental_divergence.json")


# ============================================================
# STUDY 2: p3k14c DIVERGENCE WITHOUT AFRICA (detailed)
# ============================================================
print("\n" + "=" * 70)
print("STUDY 2: p3k14c DIVERGENCE WITHOUT AFRICA (CRITICAL TEST)")
print("=" * 70)

# Already computed in Study 1, but provide more detail
# Also compute "without Egypt specifically" (Africa minus Egypt)
egypt_mask = np.array([
    s['continent'] == 'Africa' and 22 <= s['lat'] <= 32 and 25 <= s['lon'] <= 36
    for s in p3k_sites
])
no_egypt_mono = np.array([
    s['classification'] == 'monumental' and not (22 <= s['lat'] <= 32 and 25 <= s['lon'] <= 36)
    for s in p3k_sites
])
no_egypt_dome = np.array([
    s['classification'] == 'domestic' and not (22 <= s['lat'] <= 32 and 25 <= s['lon'] <= 36)
    for s in p3k_sites
])

n_egypt = int(np.sum(egypt_mask))
print(f"  Egyptian sites removed: {n_egypt}")
print(f"  Remaining monumental: {int(np.sum(no_egypt_mono))}, domestic: {int(np.sum(no_egypt_dome))}")

mono_noe_res = run_mc(p3k_lats[no_egypt_mono], p3k_lons[no_egypt_mono])
dome_noe_res = run_mc(p3k_lats[no_egypt_dome], p3k_lons[no_egypt_dome])
div_noe = round(mono_noe_res['z_score'] - dome_noe_res['z_score'], 2)

print(f"  Without Egypt: Mon Z={mono_noe_res['z_score']}, Dom Z={dome_noe_res['z_score']}, Div={div_noe}")
print(f"  Without Africa: Mon Z={mono_noaf_res['z_score']}, Dom Z={dome_noaf_res['z_score']}, Div={div_noaf}")

# Compute the full-dataset divergence for comparison
all_mono = np.array([s['classification'] == 'monumental' for s in p3k_sites])
all_dome = np.array([s['classification'] == 'domestic' for s in p3k_sites])
mono_all_res = run_mc(p3k_lats[all_mono], p3k_lons[all_mono])
dome_all_res = run_mc(p3k_lats[all_dome], p3k_lons[all_dome])
div_all = round(mono_all_res['z_score'] - dome_all_res['z_score'], 2)
print(f"  Full dataset: Mon Z={mono_all_res['z_score']}, Dom Z={dome_all_res['z_score']}, Div={div_all}")

save_json({
    "meta": {"date": "2026-03-18", "dataset": "p3k14c", "threshold_km": THRESHOLD_KM,
             "n_trials": N_TRIALS, "description": "Critical test: does divergence survive without Africa/Egypt?"},
    "full_dataset": {
        "monumental": mono_all_res, "domestic": dome_all_res, "divergence": div_all
    },
    "without_egypt_only": {
        "n_egypt_removed": n_egypt,
        "monumental": mono_noe_res, "domestic": dome_noe_res, "divergence": div_noe
    },
    "without_all_africa": {
        "monumental": mono_noaf_res, "domestic": dome_noaf_res, "divergence": div_noaf
    },
    "interpretation": {
        "egypt_bias_dead": bool(div_noaf > 3),
        "egypt_bias_possible": bool(div_noaf < 2),
        "egypt_bias_marginal": bool(2 <= div_noaf <= 3),
    }
}, "without_africa.json")


# ============================================================
# STUDY 3: p3k14c DIVERGENCE WITHIN SOUTH AMERICA ONLY
# ============================================================
print("\n" + "=" * 70)
print("STUDY 3: p3k14c DIVERGENCE WITHIN SOUTH AMERICA ONLY")
print("=" * 70)

sa_mono_mask = np.array([s['continent'] == 'South America' and s['classification'] == 'monumental' for s in p3k_sites])
sa_dome_mask = np.array([s['continent'] == 'South America' and s['classification'] == 'domestic' for s in p3k_sites])
sa_all_mask = np.array([s['continent'] == 'South America' for s in p3k_sites])

n_sa_mono = int(np.sum(sa_mono_mask))
n_sa_dome = int(np.sum(sa_dome_mask))
n_sa_all = int(np.sum(sa_all_mask))

print(f"  SA total: {n_sa_all}, monumental: {n_sa_mono}, domestic: {n_sa_dome}")

# List SA monumental site names
sa_mono_names = Counter(s['site_name'] for s in p3k_sites if s['continent'] == 'South America' and s['classification'] == 'monumental')
print(f"  SA monumental site examples: {[n for n, _ in sa_mono_names.most_common(10)]}")

sa_dome_names = Counter(s['site_name'] for s in p3k_sites if s['continent'] == 'South America' and s['classification'] == 'domestic')
print(f"  SA domestic site examples: {[n for n, _ in sa_dome_names.most_common(10)]}")

# Run SA overall
sa_all_res = run_mc(p3k_lats[sa_all_mask], p3k_lons[sa_all_mask])
print(f"  SA overall: Z={sa_all_res['z_score']}")

# Run SA monumental
if n_sa_mono >= 5:
    sa_mono_res = run_mc(p3k_lats[sa_mono_mask], p3k_lons[sa_mono_mask])
    print(f"  SA monumental: Z={sa_mono_res['z_score']}, obs={sa_mono_res['observed']}")
else:
    sa_mono_res = {"n_sites": n_sa_mono, "z_score": 0, "too_few": True}
    print(f"  SA monumental: too few ({n_sa_mono})")

# Run SA domestic
if n_sa_dome >= 5:
    sa_dome_res = run_mc(p3k_lats[sa_dome_mask], p3k_lons[sa_dome_mask])
    print(f"  SA domestic: Z={sa_dome_res['z_score']}, obs={sa_dome_res['observed']}")
else:
    sa_dome_res = {"n_sites": n_sa_dome, "z_score": 0, "too_few": True}
    print(f"  SA domestic: too few ({n_sa_dome})")

if not sa_mono_res.get("too_few") and not sa_dome_res.get("too_few"):
    sa_div = round(sa_mono_res['z_score'] - sa_dome_res['z_score'], 2)
else:
    sa_div = "N/A"
print(f"  SA divergence: {sa_div}")

save_json({
    "meta": {"date": "2026-03-18", "dataset": "p3k14c (South America only)",
             "threshold_km": THRESHOLD_KM, "n_trials": N_TRIALS},
    "overall": sa_all_res,
    "monumental": sa_mono_res,
    "domestic": sa_dome_res,
    "divergence": sa_div,
    "monumental_site_examples": [n for n, _ in sa_mono_names.most_common(20)],
    "domestic_site_examples": [n for n, _ in sa_dome_names.most_common(20)],
    "interpretation": {
        "overall_enrichment": sa_all_res.get('z_score', 0),
        "type_divergence_despite_no_overall": bool(sa_div != "N/A" and isinstance(sa_div, (int, float)) and abs(sa_div) > 2 and sa_all_res.get('z_score', 0) < 2),
    }
}, "south_america_divergence.json")


# ============================================================
# STUDY 4: p3k14c EGYPT-ONLY DIVERGENCE
# ============================================================
print("\n" + "=" * 70)
print("STUDY 4: p3k14c EGYPT-ONLY DIVERGENCE")
print("=" * 70)

egypt_sites = [s for s in p3k_sites if 22 <= s['lat'] <= 32 and 25 <= s['lon'] <= 36]
egypt_mono = [s for s in egypt_sites if s['classification'] == 'monumental']
egypt_dome = [s for s in egypt_sites if s['classification'] == 'domestic']
egypt_uncl = [s for s in egypt_sites if s['classification'] == 'unclassified']

print(f"  Egypt total: {len(egypt_sites)}")
print(f"  Egypt monumental: {len(egypt_mono)}, domestic: {len(egypt_dome)}, unclassified: {len(egypt_uncl)}")

# Dates-per-site ratio analysis (recording bias indicator)
egypt_mono_dates = [s['n_dates'] for s in egypt_mono]
egypt_dome_dates = [s['n_dates'] for s in egypt_dome]

mono_mean_dates = np.mean(egypt_mono_dates) if egypt_mono_dates else 0
dome_mean_dates = np.mean(egypt_dome_dates) if egypt_dome_dates else 0
mono_median_dates = np.median(egypt_mono_dates) if egypt_mono_dates else 0
dome_median_dates = np.median(egypt_dome_dates) if egypt_dome_dates else 0

print(f"  Egyptian monumental dates/site: mean={mono_mean_dates:.2f}, median={mono_median_dates:.1f}")
print(f"  Egyptian domestic dates/site: mean={dome_mean_dates:.2f}, median={dome_median_dates:.1f}")

dates_ratio = mono_mean_dates / dome_mean_dates if dome_mean_dates > 0 else float('inf')
print(f"  Dates ratio (mono/dome): {dates_ratio:.2f}")

# Run MC for Egypt
egypt_lats = np.array([s['lat'] for s in egypt_sites])
egypt_lons = np.array([s['lon'] for s in egypt_sites])

egypt_mono_lats = np.array([s['lat'] for s in egypt_mono])
egypt_mono_lons = np.array([s['lon'] for s in egypt_mono])
egypt_dome_lats = np.array([s['lat'] for s in egypt_dome])
egypt_dome_lons = np.array([s['lon'] for s in egypt_dome])

egypt_all_res = run_mc(egypt_lats, egypt_lons) if len(egypt_sites) >= 5 else {"n_sites": len(egypt_sites), "too_few": True}
egypt_mono_res = run_mc(egypt_mono_lats, egypt_mono_lons) if len(egypt_mono) >= 5 else {"n_sites": len(egypt_mono), "too_few": True}
egypt_dome_res = run_mc(egypt_dome_lats, egypt_dome_lons) if len(egypt_dome) >= 5 else {"n_sites": len(egypt_dome), "too_few": True}

if not egypt_mono_res.get("too_few") and not egypt_dome_res.get("too_few"):
    egypt_div = round(egypt_mono_res['z_score'] - egypt_dome_res['z_score'], 2)
else:
    egypt_div = "N/A"

print(f"  Egypt overall Z: {egypt_all_res.get('z_score', 'N/A')}")
print(f"  Egypt monumental Z: {egypt_mono_res.get('z_score', 'N/A')}")
print(f"  Egypt domestic Z: {egypt_dome_res.get('z_score', 'N/A')}")
print(f"  Egypt divergence: {egypt_div}")

# Show monumental site names from Egypt
egypt_mono_names = Counter(s['site_name'] for s in egypt_mono)
egypt_dome_names = Counter(s['site_name'] for s in egypt_dome)
print(f"  Egypt monumental sites: {[n for n, _ in egypt_mono_names.most_common(15)]}")
print(f"  Egypt domestic sites: {[n for n, _ in egypt_dome_names.most_common(15)]}")

save_json({
    "meta": {"date": "2026-03-18", "dataset": "p3k14c (Egypt: lat 22-32, lon 25-36)",
             "threshold_km": THRESHOLD_KM, "n_trials": N_TRIALS},
    "n_egypt_sites": len(egypt_sites),
    "n_monumental": len(egypt_mono),
    "n_domestic": len(egypt_dome),
    "n_unclassified": len(egypt_uncl),
    "overall": egypt_all_res,
    "monumental": egypt_mono_res,
    "domestic": egypt_dome_res,
    "divergence": egypt_div,
    "recording_bias_indicator": {
        "monumental_mean_dates_per_site": round(float(mono_mean_dates), 2),
        "domestic_mean_dates_per_site": round(float(dome_mean_dates), 2),
        "monumental_median_dates_per_site": round(float(mono_median_dates), 1),
        "domestic_median_dates_per_site": round(float(dome_median_dates), 1),
        "ratio_mono_to_dome": round(float(dates_ratio), 2),
        "bias_likely": bool(dates_ratio > 3),
        "bias_marginal": bool(1.5 < dates_ratio <= 3),
        "bias_unlikely": bool(dates_ratio <= 1.5),
    },
    "monumental_site_names": [n for n, _ in egypt_mono_names.most_common(20)],
    "domestic_site_names": [n for n, _ in egypt_dome_names.most_common(20)],
}, "egypt_only.json")


# ============================================================
# STUDY 5: KEYWORD CLASSIFICATION AUDIT
# ============================================================
print("\n" + "=" * 70)
print("STUDY 5: KEYWORD CLASSIFICATION AUDIT")
print("=" * 70)

# Sites within 50 km of the great circle
dists = gc_dist_vec(POLE_LAT, POLE_LON, p3k_lats, p3k_lons)
on_line_mask = dists <= THRESHOLD_KM
on_line_sites = [s for s, m in zip(p3k_sites, on_line_mask) if m]

print(f"  Sites within {THRESHOLD_KM}km of great circle: {len(on_line_sites)}")

# Top 20 most common SiteNames on the line
on_line_names = Counter(s['site_name'] for s in on_line_sites)
print(f"\n  Top 20 most common SiteNames on the line:")
top_names = on_line_names.most_common(20)
audit_entries = []
for name, count in top_names:
    classification = classify_p3k_site(name)
    # Find which keyword matched
    matched_kw = None
    if classification == "monumental":
        for kw in MONUMENTAL_KW:
            if kw in name.lower():
                matched_kw = kw; break
    elif classification == "domestic":
        for kw in DOMESTIC_KW:
            if kw in name.lower():
                matched_kw = kw; break
    entry = {
        "site_name": name,
        "count_on_line": count,
        "classification": classification,
        "matched_keyword": matched_kw,
    }
    audit_entries.append(entry)
    flag = " ***AMBIGUOUS***" if classification == "unclassified" else ""
    print(f"    {count:3d}x  [{classification:12s}]  {name}{flag}")

# Classification breakdown for on-line sites
on_line_class = Counter(s['classification'] for s in on_line_sites)
total_on = len(on_line_sites)
pct_uncl = on_line_class.get('unclassified', 0) / total_on * 100 if total_on > 0 else 0
pct_mono = on_line_class.get('monumental', 0) / total_on * 100 if total_on > 0 else 0
pct_dome = on_line_class.get('domestic', 0) / total_on * 100 if total_on > 0 else 0

print(f"\n  On-line classification breakdown:")
print(f"    Monumental: {on_line_class.get('monumental', 0)} ({pct_mono:.1f}%)")
print(f"    Domestic: {on_line_class.get('domestic', 0)} ({pct_dome:.1f}%)")
print(f"    Unclassified: {on_line_class.get('unclassified', 0)} ({pct_uncl:.1f}%)")

# Compare to off-line breakdown
off_line_sites = [s for s, m in zip(p3k_sites, on_line_mask) if not m]
off_line_class = Counter(s['classification'] for s in off_line_sites)
total_off = len(off_line_sites)
pct_mono_off = off_line_class.get('monumental', 0) / total_off * 100 if total_off > 0 else 0
pct_dome_off = off_line_class.get('domestic', 0) / total_off * 100 if total_off > 0 else 0

print(f"\n  Off-line classification breakdown for comparison:")
print(f"    Monumental: {off_line_class.get('monumental', 0)} ({pct_mono_off:.1f}%)")
print(f"    Domestic: {off_line_class.get('domestic', 0)} ({pct_dome_off:.1f}%)")

# On-line by continent
on_line_by_cont = defaultdict(lambda: {"monumental": 0, "domestic": 0, "unclassified": 0, "total": 0})
for s in on_line_sites:
    on_line_by_cont[s['continent']][s['classification']] += 1
    on_line_by_cont[s['continent']]["total"] += 1

print(f"\n  On-line sites by continent and classification:")
for cont in sorted(on_line_by_cont.keys()):
    c = on_line_by_cont[cont]
    print(f"    {cont}: mono={c['monumental']}, dome={c['domestic']}, uncl={c['unclassified']}, total={c['total']}")

save_json({
    "meta": {"date": "2026-03-18", "dataset": "p3k14c", "threshold_km": THRESHOLD_KM},
    "n_on_line": total_on,
    "on_line_classification": {
        "monumental": on_line_class.get('monumental', 0),
        "domestic": on_line_class.get('domestic', 0),
        "unclassified": on_line_class.get('unclassified', 0),
        "pct_monumental": round(pct_mono, 1),
        "pct_domestic": round(pct_dome, 1),
        "pct_unclassified": round(pct_uncl, 1),
    },
    "off_line_classification": {
        "monumental": off_line_class.get('monumental', 0),
        "domestic": off_line_class.get('domestic', 0),
        "pct_monumental": round(pct_mono_off, 1),
        "pct_domestic": round(pct_dome_off, 1),
    },
    "top_20_site_names": audit_entries,
    "on_line_by_continent": dict(on_line_by_cont),
    "monumental_keywords": MONUMENTAL_KW,
    "domestic_keywords": DOMESTIC_KW,
}, "classification_audit.json")


# ============================================================
# FINAL SUMMARY
# ============================================================
elapsed = time.time() - t0
print(f"\n{'=' * 70}")
print(f"ALL STUDIES COMPLETE — {elapsed:.0f}s elapsed")
print(f"{'=' * 70}")
print(f"Output directory: {OUT_DIR}")

print(f"\n  CRITICAL RESULTS:")
print(f"  Full dataset divergence:      {div_all}")
print(f"  Without Egypt divergence:     {div_noe}")
print(f"  Without Africa divergence:    {div_noaf}")
print(f"  South America divergence:     {sa_div}")
print(f"  Egypt-only divergence:        {egypt_div}")
print(f"  Egypt dates ratio (mono/dom): {dates_ratio:.2f}")
