#!/usr/bin/env python3
"""
New Database Analysis
=====================
Great Circle tests on three newly downloaded databases:
1. DARE (Digital Atlas of the Roman Empire) — 29,760 ancient places
2. Historic England Scheduled Monuments — 20,026 monuments
3. EAMENA Sistan Heritage Places — 1,668 sites (Iran/Afghanistan)

For each: overall enrichment, monument-settlement divergence where applicable.
"""

import csv, math, json, os, re, time
import numpy as np
from collections import Counter, defaultdict

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
THRESHOLD_KM = 50
N_TRIALS = 1000  # directive says 1000-trial MC

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "new_databases")
os.makedirs(OUT_DIR, exist_ok=True)

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
# DATABASE 1: DARE (Digital Atlas of the Roman Empire)
# ============================================================
print("=" * 70)
print("DATABASE 1: DARE (Digital Atlas of the Roman Empire)")
print("=" * 70)

t0 = time.time()
np.random.seed(42)

# Load type mapping
with open(os.path.join(BASE_DIR, "dare_data", "type_mapping.json")) as f:
    DARE_TYPE_MAP = json.load(f)

# Load GeoJSON
with open(os.path.join(BASE_DIR, "dare_data", "places2.geojson")) as f:
    dare_data = json.load(f)

dare_sites = []
for feat in dare_data["features"]:
    props = feat["properties"]
    geom = feat["geometry"]
    if geom["type"] != "Point":
        continue
    lon, lat = geom["coordinates"]
    if lat == 0 and lon == 0:
        continue
    if not (-90 <= lat <= 90 and -180 <= lon <= 180):
        continue
    type_code = str(props.get("type", ""))
    type_name = DARE_TYPE_MAP.get(type_code, "unknown")
    dare_sites.append({
        "lat": lat, "lon": lon,
        "name": props.get("latin", "") or props.get("modern", ""),
        "type_code": type_code,
        "type_name": type_name,
        "major": props.get("major", 0),
    })

print(f"  Loaded {len(dare_sites)} DARE sites")

# Classify DARE types as monumental vs settlement
DARE_MONUMENTAL = {"temple", "amphitheater", "theatre", "odeon", "stadium", "circus",
                   "monument", "necropolis", "tumulus", "aqueduct/dam/cistern/spring",
                   "basilica", "forum", "bath", "lighthouse", "bouleuterion"}
DARE_SETTLEMENT = {"city", "town", "civitas", "villa", "settlement",
                   "Late-Roman oppidum", "Iron-Age oppidum", "oasis"}
DARE_MILITARY = {"fort", "fortress", "camp", "fortlet/tower", "wall"}
DARE_INFRASTRUCTURE = {"station", "bridge", "port", "road/milestone", "mine/quarry",
                       "production", "well"}

def classify_dare(type_name):
    if type_name in DARE_MONUMENTAL:
        return "monumental"
    elif type_name in DARE_SETTLEMENT:
        return "settlement"
    elif type_name in DARE_MILITARY:
        return "military"
    elif type_name in DARE_INFRASTRUCTURE:
        return "infrastructure"
    else:
        return "other"

for s in dare_sites:
    s["category"] = classify_dare(s["type_name"])

cat_counts = Counter(s["category"] for s in dare_sites)
type_counts = Counter(s["type_name"] for s in dare_sites)
print(f"  Category breakdown: {dict(cat_counts)}")
print(f"  Top types: {type_counts.most_common(10)}")

dare_lats = np.array([s["lat"] for s in dare_sites])
dare_lons = np.array([s["lon"] for s in dare_sites])

# Overall GC test
print(f"\n  Running overall GC test ({len(dare_sites)} sites, {N_TRIALS} trials)...")
dare_overall = run_mc(dare_lats, dare_lons)
print(f"  Overall: Z={dare_overall['z_score']}, obs={dare_overall['observed']}, "
      f"exp={dare_overall['expected']}, enrich={dare_overall['enrichment']}x")

# Monument-settlement divergence
dare_mono_mask = np.array([s["category"] == "monumental" for s in dare_sites])
dare_sett_mask = np.array([s["category"] == "settlement" for s in dare_sites])
dare_mil_mask = np.array([s["category"] == "military" for s in dare_sites])

print(f"\n  Running monumental MC ({int(np.sum(dare_mono_mask))} sites)...")
dare_mono_res = run_mc(dare_lats[dare_mono_mask], dare_lons[dare_mono_mask])
print(f"  Monumental: Z={dare_mono_res['z_score']}, obs={dare_mono_res['observed']}")

print(f"  Running settlement MC ({int(np.sum(dare_sett_mask))} sites)...")
dare_sett_res = run_mc(dare_lats[dare_sett_mask], dare_lons[dare_sett_mask])
print(f"  Settlement: Z={dare_sett_res['z_score']}, obs={dare_sett_res['observed']}")

print(f"  Running military MC ({int(np.sum(dare_mil_mask))} sites)...")
dare_mil_res = run_mc(dare_lats[dare_mil_mask], dare_lons[dare_mil_mask])
print(f"  Military: Z={dare_mil_res['z_score']}, obs={dare_mil_res['observed']}")

dare_div = round(dare_mono_res["z_score"] - dare_sett_res["z_score"], 2)
print(f"\n  DARE Divergence (monumental - settlement): {dare_div}")

# On-line site breakdown
dists = gc_dist_vec(POLE_LAT, POLE_LON, dare_lats, dare_lons)
on_line_mask = dists <= THRESHOLD_KM
on_line_cats = Counter(s["category"] for s, m in zip(dare_sites, on_line_mask) if m)
on_line_types = Counter(s["type_name"] for s, m in zip(dare_sites, on_line_mask) if m)
print(f"\n  On-line sites ({int(np.sum(on_line_mask))}): {dict(on_line_cats)}")
print(f"  On-line types: {on_line_types.most_common(10)}")

# On-line site names (for audit)
on_line_names = [(s["name"], s["type_name"]) for s, m in zip(dare_sites, on_line_mask) if m]

save_json({
    "meta": {"date": "2026-03-18", "dataset": "DARE", "n_sites": len(dare_sites),
             "threshold_km": THRESHOLD_KM, "n_trials": N_TRIALS,
             "source": "github.com/johaahlf/dare/places2.geojson"},
    "category_breakdown": dict(cat_counts),
    "type_breakdown": dict(type_counts.most_common(20)),
    "overall": dare_overall,
    "monumental": dare_mono_res,
    "settlement": dare_sett_res,
    "military": dare_mil_res,
    "divergence": dare_div,
    "on_line_category_breakdown": dict(on_line_cats),
    "on_line_type_breakdown": dict(on_line_types.most_common(15)),
    "on_line_sites": [{"name": n, "type": t} for n, t in on_line_names[:50]],
    "type_mapping": DARE_TYPE_MAP,
}, "dare_analysis.json")

print(f"  DARE analysis done: {time.time()-t0:.1f}s")


# ============================================================
# DATABASE 2: HISTORIC ENGLAND SCHEDULED MONUMENTS
# ============================================================
print(f"\n{'=' * 70}")
print("DATABASE 2: HISTORIC ENGLAND SCHEDULED MONUMENTS")
print("=" * 70)

t1 = time.time()

# Parse the CSV — coordinates are in "POINT (lon lat)" format
he_sites = []
with open(os.path.join(BASE_DIR, "historic_england_data", "scheduled-monument.csv"), encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        point = row.get("point", "")
        name = row.get("name", "")
        if not point:
            continue
        # Parse "POINT (lon lat)"
        m = re.match(r'POINT\s*\(\s*([-\d.]+)\s+([-\d.]+)\s*\)', point)
        if not m:
            continue
        lon, lat = float(m.group(1)), float(m.group(2))
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue
        he_sites.append({"lat": lat, "lon": lon, "name": name})

print(f"  Loaded {len(he_sites)} Historic England sites")

# Classify HE monuments using name keywords
HE_MONUMENTAL_KW = ["stone circle", "henge", "standing stone", "menhir",
                     "stone row", "dolmen", "burial chamber", "barrow",
                     "cairn", "tumulus", "passage grave", "tomb", "pyramid",
                     "abbey", "priory", "church", "cathedral", "chapel",
                     "monastery", "cross", "temple", "amphitheatre",
                     "stone setting", "cursus"]
HE_SETTLEMENT_KW = ["settlement", "village", "farmstead", "homestead",
                     "enclosure", "hillfort", "camp", "castle", "fort",
                     "moated site", "manor", "house"]

def classify_he(name):
    nl = name.lower()
    for kw in HE_MONUMENTAL_KW:
        if kw in nl:
            return "monumental"
    for kw in HE_SETTLEMENT_KW:
        if kw in nl:
            return "settlement"
    return "other"

for s in he_sites:
    s["category"] = classify_he(s["name"])

he_cats = Counter(s["category"] for s in he_sites)
print(f"  Category breakdown: {dict(he_cats)}")

he_lats = np.array([s["lat"] for s in he_sites])
he_lons = np.array([s["lon"] for s in he_sites])

# The Great Circle barely touches England — this should be a negative control
# Let's first check how many sites are within various thresholds
for t in [50, 100, 200, 500]:
    dists = gc_dist_vec(POLE_LAT, POLE_LON, he_lats, he_lons)
    n_within = int(np.sum(dists <= t))
    print(f"  Sites within {t}km of GC: {n_within}")

# Overall GC test
print(f"\n  Running overall GC test ({len(he_sites)} sites, {N_TRIALS} trials)...")
he_overall = run_mc(he_lats, he_lons)
print(f"  Overall: Z={he_overall['z_score']}, obs={he_overall['observed']}, "
      f"exp={he_overall['expected']}, enrich={he_overall['enrichment']}x")

# Monument-settlement divergence
he_mono_mask = np.array([s["category"] == "monumental" for s in he_sites])
he_sett_mask = np.array([s["category"] == "settlement" for s in he_sites])

print(f"\n  Running monumental MC ({int(np.sum(he_mono_mask))} sites)...")
he_mono_res = run_mc(he_lats[he_mono_mask], he_lons[he_mono_mask])
print(f"  Monumental: Z={he_mono_res['z_score']}, obs={he_mono_res['observed']}")

print(f"  Running settlement MC ({int(np.sum(he_sett_mask))} sites)...")
he_sett_res = run_mc(he_lats[he_sett_mask], he_lons[he_sett_mask])
print(f"  Settlement: Z={he_sett_res['z_score']}, obs={he_sett_res['observed']}")

he_div = round(he_mono_res["z_score"] - he_sett_res["z_score"], 2)
print(f"\n  HE Divergence (monumental - settlement): {he_div}")

# Check what's nearest the line
dists_he = gc_dist_vec(POLE_LAT, POLE_LON, he_lats, he_lons)
nearest_idx = np.argsort(dists_he)[:10]
print(f"\n  10 nearest HE sites to the Great Circle:")
for i in nearest_idx:
    print(f"    {dists_he[i]:.1f}km: {he_sites[i]['name']} ({he_sites[i]['category']})")

save_json({
    "meta": {"date": "2026-03-18", "dataset": "Historic England Scheduled Monuments",
             "n_sites": len(he_sites), "threshold_km": THRESHOLD_KM, "n_trials": N_TRIALS,
             "source": "files.planning.data.gov.uk/dataset/scheduled-monument.csv",
             "note": "Expected negative control - Great Circle avoids England"},
    "category_breakdown": dict(he_cats),
    "overall": he_overall,
    "monumental": he_mono_res,
    "settlement": he_sett_res,
    "divergence": he_div,
    "nearest_to_gc": [{"distance_km": round(float(dists_he[i]), 1),
                        "name": he_sites[i]["name"],
                        "category": he_sites[i]["category"]}
                       for i in nearest_idx],
    "geographic_note": "All sites in England (lat ~50-56N, lon ~-6 to 2E). "
                       "Great Circle passes through southern England peripherally.",
}, "historic_england_analysis.json")

print(f"  HE analysis done: {time.time()-t1:.1f}s")


# ============================================================
# DATABASE 3: EAMENA SISTAN HERITAGE PLACES
# ============================================================
print(f"\n{'=' * 70}")
print("DATABASE 3: EAMENA SISTAN HERITAGE PLACES")
print("=" * 70)

t2 = time.time()

with open(os.path.join(BASE_DIR, "eamena_data", "sistan_part1_hps.geojson")) as f:
    eamena_data = json.load(f)

eamena_sites = []
for feat in eamena_data["features"]:
    geom = feat["geometry"]
    if geom is None or geom["type"] != "Point":
        continue
    coords = geom["coordinates"]
    lon, lat = coords[0], coords[1]
    if lat == 0 and lon == 0:
        continue
    if not (-90 <= lat <= 90 and -180 <= lon <= 180):
        continue
    props = feat["properties"]
    site_type = props.get("Site Feature Form Type", "")
    country = props.get("Country Type", "")
    eamena_id = props.get("EAMENA ID", "")
    eamena_sites.append({
        "lat": lat, "lon": lon,
        "eamena_id": eamena_id,
        "site_type": site_type,
        "country": country,
    })

print(f"  Loaded {len(eamena_sites)} EAMENA Sistan sites")

# Check geographic extent
eamena_lats = np.array([s["lat"] for s in eamena_sites])
eamena_lons = np.array([s["lon"] for s in eamena_sites])
print(f"  Lat range: {eamena_lats.min():.2f} to {eamena_lats.max():.2f}")
print(f"  Lon range: {eamena_lons.min():.2f} to {eamena_lons.max():.2f}")

# Country breakdown
countries = Counter(s["country"] for s in eamena_sites)
print(f"  Countries: {dict(countries)}")

# Type breakdown
types = Counter(s["site_type"] for s in eamena_sites)
print(f"  Types: {dict(types)}")

# How far from the Great Circle?
dists_eamena = gc_dist_vec(POLE_LAT, POLE_LON, eamena_lats, eamena_lons)
print(f"\n  Distance to GC: min={dists_eamena.min():.1f}km, "
      f"median={np.median(dists_eamena):.1f}km, max={dists_eamena.max():.1f}km")
for t in [50, 100, 200, 500]:
    n_within = int(np.sum(dists_eamena <= t))
    print(f"  Sites within {t}km: {n_within}")

# Run GC test
print(f"\n  Running overall GC test ({len(eamena_sites)} sites, {N_TRIALS} trials)...")
eamena_overall = run_mc(eamena_lats, eamena_lons)
print(f"  Overall: Z={eamena_overall['z_score']}, obs={eamena_overall['observed']}, "
      f"exp={eamena_overall['expected']}, enrich={eamena_overall['enrichment']}x")

save_json({
    "meta": {"date": "2026-03-18", "dataset": "EAMENA Sistan Heritage Places",
             "n_sites": len(eamena_sites), "threshold_km": THRESHOLD_KM, "n_trials": N_TRIALS,
             "source": "zenodo.org/records/10375902",
             "coverage": "Eastern Iran / SW Afghanistan (Sistan region)"},
    "geographic_extent": {
        "lat_min": round(float(eamena_lats.min()), 2),
        "lat_max": round(float(eamena_lats.max()), 2),
        "lon_min": round(float(eamena_lons.min()), 2),
        "lon_max": round(float(eamena_lons.max()), 2),
    },
    "country_breakdown": dict(countries),
    "type_breakdown": dict(types),
    "overall": eamena_overall,
    "note": "Small regional dataset (Sistan only). Full EAMENA database requires contributor access.",
}, "eamena_analysis.json")

print(f"  EAMENA analysis done: {time.time()-t2:.1f}s")


# ============================================================
# SUMMARY
# ============================================================
elapsed = time.time() - t0
print(f"\n{'=' * 70}")
print(f"ALL ANALYSES COMPLETE — {elapsed:.0f}s elapsed")
print(f"{'=' * 70}")
print(f"Output directory: {OUT_DIR}")

print(f"\n  SUMMARY:")
print(f"  {'Database':<40} {'N':>8} {'Z':>8} {'Obs':>6} {'Exp':>8} {'Enrich':>8} {'Div':>8}")
print(f"  {'-'*40} {'-'*8} {'-'*8} {'-'*6} {'-'*8} {'-'*8} {'-'*8}")
print(f"  {'DARE (all)':<40} {dare_overall['n_sites']:>8} {dare_overall['z_score']:>8} "
      f"{dare_overall['observed']:>6} {dare_overall['expected']:>8} {dare_overall['enrichment']:>8} {'':>8}")
print(f"  {'DARE monumental':<40} {dare_mono_res['n_sites']:>8} {dare_mono_res['z_score']:>8} "
      f"{dare_mono_res['observed']:>6} {dare_mono_res['expected']:>8} {dare_mono_res['enrichment']:>8} {'':>8}")
print(f"  {'DARE settlement':<40} {dare_sett_res['n_sites']:>8} {dare_sett_res['z_score']:>8} "
      f"{dare_sett_res['observed']:>6} {dare_sett_res['expected']:>8} {dare_sett_res['enrichment']:>8} {'':>8}")
print(f"  {'DARE divergence':<40} {'':>8} {'':>8} {'':>6} {'':>8} {'':>8} {dare_div:>8}")
print(f"  {'Historic England (all)':<40} {he_overall['n_sites']:>8} {he_overall['z_score']:>8} "
      f"{he_overall['observed']:>6} {he_overall['expected']:>8} {he_overall['enrichment']:>8} {'':>8}")
print(f"  {'HE monumental':<40} {he_mono_res['n_sites']:>8} {he_mono_res['z_score']:>8} "
      f"{he_mono_res['observed']:>6} {he_mono_res['expected']:>8} {he_mono_res['enrichment']:>8} {'':>8}")
print(f"  {'HE settlement':<40} {he_sett_res['n_sites']:>8} {he_sett_res['z_score']:>8} "
      f"{he_sett_res['observed']:>6} {he_sett_res['expected']:>8} {he_sett_res['enrichment']:>8} {'':>8}")
print(f"  {'HE divergence':<40} {'':>8} {'':>8} {'':>6} {'':>8} {'':>8} {he_div:>8}")
print(f"  {'EAMENA Sistan (all)':<40} {eamena_overall['n_sites']:>8} {eamena_overall['z_score']:>8} "
      f"{eamena_overall['observed']:>6} {eamena_overall['expected']:>8} {eamena_overall['enrichment']:>8} {'':>8}")
