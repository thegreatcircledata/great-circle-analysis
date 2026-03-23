#!/usr/bin/env python3
"""
Egypt Preservation Test — 1000 MC Trials
==========================================
Re-run of egypt_preservation_test.py with N_MC_TRIALS=1000 for baseline
Monte Carlo Z-scores, matching the main paper's methodology.

Per Perplexity methodological note: the original test used 200 MC trials
for each Z-score computation. This re-run uses 1000 for consistency
with the main paper's methodology.

The +350 estimated test still uses 1000 outer iterations.
"""

import csv, math, random, json, os, time
import numpy as np
from collections import Counter

# ============================================================
# CONFIGURATION — ONLY CHANGE: N_MC_TRIALS = 1000
# ============================================================
EARTH_R_KM = 6371
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
THRESHOLD_KM = 50
N_MC_TRIALS = 1000        # <<< UPGRADED from 200 to 1000
N_MONTE_CARLO_ITER = 1000  # outer MC for +350 test (unchanged)

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# 58 BURIED SETTLEMENT SITES (from Perplexity reports)
# ============================================================

VALLEY_SITES = [
    {"name": "Memphis / Mit Rahina",              "lat": 29.851, "lon": 31.253},
    {"name": "Thebes East Bank (Luxor)",          "lat": 25.687, "lon": 32.639},
    {"name": "Thebes West Bank (Malkata/Aten)",   "lat": 25.720, "lon": 32.601},
    {"name": "Amarna / Akhetaten periphery",      "lat": 27.645, "lon": 30.899},
    {"name": "Hermopolis / el-Ashmunein",         "lat": 27.781, "lon": 30.805},
    {"name": "Asyut / Lycopolis",                 "lat": 27.167, "lon": 31.133},
    {"name": "Abydos / Wah-Sut",                  "lat": 26.184, "lon": 31.920},
    {"name": "Dendara / Tentyra",                 "lat": 26.142, "lon": 32.670},
    {"name": "Tell Edfu / Djeba",                 "lat": 24.978, "lon": 32.874},
    {"name": "Esna / Latopolis",                  "lat": 25.293, "lon": 32.553},
    {"name": "Elephantine / Abu",                 "lat": 24.088, "lon": 32.888},
    {"name": "el-Kab / Nekheb",                   "lat": 25.115, "lon": 32.797},
    {"name": "Hierakonpolis / Nekhen",             "lat": 25.099, "lon": 32.779},
    {"name": "Kom Ombo / Nubt",                   "lat": 24.451, "lon": 32.928},
    {"name": "Shashotep / Shutb",                 "lat": 27.117, "lon": 31.150},
    {"name": "Qus / Apollinopolis Parva",         "lat": 26.078, "lon": 32.754},
    {"name": "Naqada / Nubt",                     "lat": 25.900, "lon": 32.696},
    {"name": "Armant / Hermonthis",               "lat": 25.617, "lon": 32.531},
    {"name": "Tod / Djerty",                      "lat": 25.595, "lon": 32.534},
    {"name": "Gebelein / Inerty",                 "lat": 25.481, "lon": 32.471},
    {"name": "el-Mo'alla / Hefat",                "lat": 25.346, "lon": 32.535},
    {"name": "Medamud / Madu",                    "lat": 25.734, "lon": 32.705},
]

DELTA_SITES = [
    {"name": "Tell el-Dab'a / Avaris",            "lat": 30.792, "lon": 31.833},
    {"name": "Qantir / Pi-Ramesse",               "lat": 30.800, "lon": 31.833},
    {"name": "Tell el-Fara'in / Buto",            "lat": 31.200, "lon": 30.733},
    {"name": "Sa el-Hagar / Sais",                "lat": 30.967, "lon": 30.767},
    {"name": "Tell Basta / Bubastis",             "lat": 30.569, "lon": 31.517},
    {"name": "Tell el-Rub'a / Mendes",            "lat": 30.958, "lon": 31.516},
    {"name": "San el-Hagar / Tanis",              "lat": 30.975, "lon": 31.883},
    {"name": "Tell Atrib / Athribis",             "lat": 30.471, "lon": 31.188},
    {"name": "Abusir Bana / Busiris",             "lat": 30.912, "lon": 31.243},
    {"name": "Samannud / Sebennytos",             "lat": 30.950, "lon": 31.250},
    {"name": "Tell el-Moqdam / Leontopolis",      "lat": 30.683, "lon": 31.358},
    {"name": "Tell el-Balamun",                   "lat": 31.259, "lon": 31.517},
    {"name": "Tell Nabasha / Imet",               "lat": 30.783, "lon": 31.900},
    {"name": "Naukratis / Kom Ge'if",             "lat": 30.900, "lon": 30.583},
    {"name": "Tell el-Yahudiya",                  "lat": 30.293, "lon": 31.333},
    {"name": "Heliopolis / Iunu",                 "lat": 30.131, "lon": 31.313},
    {"name": "Tell Belim / Sethroe",              "lat": 30.978, "lon": 32.099},
    {"name": "Tell Tebilla / Onuphis",            "lat": 31.057, "lon": 31.774},
    {"name": "Tell Baqliya / Hermopolis Parva",   "lat": 30.950, "lon": 31.417},
    {"name": "Tell Ibrahim Awad",                 "lat": 30.717, "lon": 31.800},
    {"name": "Tell el-Farkha",                    "lat": 30.871, "lon": 31.603},
    {"name": "Tell el-Iswid",                     "lat": 30.833, "lon": 31.767},
    {"name": "Tell el-Samara",                    "lat": 30.850, "lon": 31.620},
    {"name": "Merimda Beni Salama",               "lat": 30.317, "lon": 30.850},
    {"name": "Tell el-Maskhuta / Pithom",         "lat": 30.550, "lon": 32.083},
    {"name": "Tell el-Retaba",                    "lat": 30.550, "lon": 31.850},
    {"name": "Behbeit el-Hagar / Per-Hebit",      "lat": 31.017, "lon": 31.517},
    {"name": "Tell Defenna / Daphnae",            "lat": 30.883, "lon": 32.167},
    {"name": "Itjtawy / near Lisht",              "lat": 29.570, "lon": 31.230},
    {"name": "Thinis / near Girga",               "lat": 26.330, "lon": 31.880},
    {"name": "Herakleopolis / Ihnasya el-Medina", "lat": 29.080, "lon": 30.930},
    {"name": "Coptos / Qift",                     "lat": 25.980, "lon": 32.820},
    {"name": "Hu / Diospolis Parva",              "lat": 26.020, "lon": 32.350},
    {"name": "Pelusium / Tell el-Farama",         "lat": 31.050, "lon": 32.550},
    {"name": "Kafr Hassan Dawood",                "lat": 30.500, "lon": 32.000},
    {"name": "Kom el-Khawaled",                   "lat": 31.100, "lon": 30.500},
]

ALL_BURIED = VALLEY_SITES + DELTA_SITES

print(f"Buried settlement sites: {len(VALLEY_SITES)} Valley + {len(DELTA_SITES)} Delta = {len(ALL_BURIED)} total")

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

def gc_dist_single(lat, lon):
    d_pole = haversine_vec(POLE_LAT, POLE_LON, np.array([lat]), np.array([lon]))[0]
    return abs(d_pole - QUARTER_CIRC)

# ============================================================
# CLASSIFICATION KEYWORDS
# ============================================================
EGYPT_MONUMENTAL_KW = ['pyramid', 'mastaba', 'temple', 'tomb', 'burial', 'cemetery',
                        'mortuary', 'funerary', 'necropolis', 'shrine', 'chapel',
                        'obelisk', 'fortress', 'fort', 'rock art', 'rock cut',
                        'megalith', 'standing stone', 'monument', 'sanctuary',
                        'enclosure', 'sacred', 'ritual']

EGYPT_SETTLEMENT_KW = ['settlement', 'village', 'town', 'city', 'habitation',
                        'domestic', 'residential', 'camp', 'farm', 'urban',
                        'tell', 'kom', 'mound']

def classify_pleiades_site(feat_types, periods=""):
    s = (str(feat_types) + " " + str(periods)).lower()
    for kw in EGYPT_MONUMENTAL_KW:
        if kw in s:
            return "monumental"
    for kw in EGYPT_SETTLEMENT_KW:
        if kw in s:
            return "settlement"
    return "unknown"

# ============================================================
# MONTE CARLO ENGINE — NOW 1000 TRIALS
# ============================================================
def run_mc(site_lats, site_lons, threshold=THRESHOLD_KM, n_trials=N_MC_TRIALS):
    n = len(site_lats)
    if n < 5:
        return {"n_sites": int(n), "too_few": True}
    dists = gc_dist_vec(POLE_LAT, POLE_LON, site_lats, site_lons)
    observed = int(np.sum(dists <= threshold))
    rand_counts = []
    for _ in range(n_trials):
        idx = np.random.randint(0, n, size=n)
        rl = site_lats[idx] + np.random.normal(0, 2, n)
        rn = site_lons[idx] + np.random.normal(0, 2, n)
        rl = np.clip(rl, -90, 90)
        rn = np.clip(rn, -180, 180)
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

def compute_divergence(mono_res, settle_res):
    if mono_res.get("too_few") or settle_res.get("too_few"):
        return "N/A"
    return round(mono_res["z_score"] - settle_res["z_score"], 2)

# ============================================================
# FLOODPLAIN RANDOM POINT GENERATOR
# ============================================================
def rand_nile_floodplain():
    if random.random() < 0.70:
        lat = random.uniform(30.0, 31.5)
        lon = random.uniform(29.5, 32.5)
        center_lon = 31.0
        half_width = 0.5 + (lat - 30.0) * 1.3
        if abs(lon - center_lon) > half_width:
            lon = center_lon + random.uniform(-half_width, half_width)
    else:
        lat = random.uniform(24.0, 30.0)
        if lat < 25.5:
            nile_lon = 32.8 + (lat - 24.0) * 0.05
        elif lat < 26.5:
            nile_lon = 32.6 - (lat - 25.5) * 0.8
        elif lat < 27.5:
            nile_lon = 31.8 - (lat - 26.5) * 0.5
        elif lat < 28.5:
            nile_lon = 31.3 - (lat - 27.5) * 0.2
        else:
            nile_lon = 31.1 - (lat - 28.5) * 0.1
        lon = nile_lon + random.uniform(-0.15, 0.15)
    return lat, lon

# ============================================================
# MAIN ANALYSIS
# ============================================================
print("=" * 70)
print("EGYPT PRESERVATION TEST — 1000 MC TRIALS")
print("=" * 70)

t0 = time.time()
np.random.seed(42)
random.seed(42)

# Load Pleiades Egypt data
print("\n--- Loading Pleiades Egypt data ---")
pleiades_csv = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")

egypt_sites = []
with open(pleiades_csv, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            lat = float(row.get('reprLat', '') or '0')
            lon = float(row.get('reprLong', '') or '0')
        except ValueError:
            continue
        if not (22.0 <= lat <= 32.0 and 25.0 <= lon <= 36.0):
            continue
        if lat == 0 and lon == 0:
            continue
        name = row.get('title', '')
        feat_types = row.get('featureTypes', '')
        periods = row.get('timePeriods', '')
        classification = classify_pleiades_site(feat_types, periods)
        egypt_sites.append({
            'name': name, 'lat': lat, 'lon': lon,
            'type': feat_types, 'classification': classification, 'source': 'pleiades'
        })

class_counts = Counter(s['classification'] for s in egypt_sites)
print(f"  Pleiades Egypt sites: {len(egypt_sites)}")
print(f"  Classification: {dict(class_counts)}")

# Compute GC distance for buried settlements
buried_with_dist = []
for s in ALL_BURIED:
    d = gc_dist_single(s['lat'], s['lon'])
    buried_with_dist.append({
        **s, 'gc_distance_km': round(d, 2),
        'within_50km': bool(d <= 50), 'classification': 'settlement',
        'source': 'perplexity_buried'
    })

n_within_50 = sum(1 for s in buried_with_dist if s['within_50km'])
mean_dist = np.mean([s['gc_distance_km'] for s in buried_with_dist])
print(f"\n  Buried sites within 50km of GC: {n_within_50}/{len(buried_with_dist)}")
print(f"  Mean GC distance: {mean_dist:.1f} km")

# ============================================================
# TEST 1: BASELINE
# ============================================================
print("\n" + "=" * 70)
print("TEST 1: BASELINE — Current Pleiades Egypt data (1000 MC trials)")
print("=" * 70)

pleiades_mono = [s for s in egypt_sites if s['classification'] == 'monumental']
pleiades_settle = [s for s in egypt_sites if s['classification'] == 'settlement']
pleiades_unknown = [s for s in egypt_sites if s['classification'] == 'unknown']

print(f"  Monumental: {len(pleiades_mono)}, Settlement: {len(pleiades_settle)}, Unknown: {len(pleiades_unknown)}")

mono_lats = np.array([s['lat'] for s in pleiades_mono])
mono_lons = np.array([s['lon'] for s in pleiades_mono])
settle_lats = np.array([s['lat'] for s in pleiades_settle])
settle_lons = np.array([s['lon'] for s in pleiades_settle])

baseline_mono = run_mc(mono_lats, mono_lons)
baseline_settle = run_mc(settle_lats, settle_lons)
baseline_D = compute_divergence(baseline_mono, baseline_settle)

print(f"  Monumental: Z={baseline_mono.get('z_score', 'N/A')}, obs={baseline_mono.get('observed', 'N/A')}/{baseline_mono.get('n_sites', 0)}")
print(f"  Settlement: Z={baseline_settle.get('z_score', 'N/A')}, obs={baseline_settle.get('observed', 'N/A')}/{baseline_settle.get('n_sites', 0)}")
print(f"  Divergence D = {baseline_D}")

# ============================================================
# TEST 2: +58 BURIED SETTLEMENTS
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: +58 BURIED SETTLEMENTS (1000 MC trials)")
print("=" * 70)

augmented_settle = pleiades_settle + buried_with_dist
aug_settle_lats = np.array([s['lat'] for s in augmented_settle])
aug_settle_lons = np.array([s['lon'] for s in augmented_settle])

print(f"  Original: {len(pleiades_settle)}, Added: {len(buried_with_dist)}, Total: {len(augmented_settle)}")

aug_mono = run_mc(mono_lats, mono_lons)
aug_settle = run_mc(aug_settle_lats, aug_settle_lons)
aug_D = compute_divergence(aug_mono, aug_settle)

print(f"  Monumental: Z={aug_mono.get('z_score', 'N/A')}")
print(f"  Settlement: Z={aug_settle.get('z_score', 'N/A')}")
print(f"  Divergence D = {aug_D}")

# ============================================================
# TEST 3: +350 ESTIMATED (Monte Carlo floodplain)
# ============================================================
print("\n" + "=" * 70)
print(f"TEST 3: +350 ESTIMATED — Monte Carlo ({N_MONTE_CARLO_ITER} iterations, {N_MC_TRIALS} MC trials each)")
print("=" * 70)

mc350_divergences = []
mc350_settle_z = []
mc350_mono_z = []

for i in range(N_MONTE_CARLO_ITER):
    if (i + 1) % 100 == 0:
        elapsed = time.time() - t0
        print(f"  Iteration {i+1}/{N_MONTE_CARLO_ITER}... ({elapsed:.0f}s elapsed)")

    rand_settlements = [rand_nile_floodplain() for _ in range(350)]
    rand_lats = np.array([p[0] for p in rand_settlements])
    rand_lons = np.array([p[1] for p in rand_settlements])

    combined_lats = np.concatenate([settle_lats, rand_lats])
    combined_lons = np.concatenate([settle_lons, rand_lons])

    settle_res = run_mc(combined_lats, combined_lons)
    mono_res = run_mc(mono_lats, mono_lons)

    if not settle_res.get("too_few") and not mono_res.get("too_few"):
        D = mono_res['z_score'] - settle_res['z_score']
        mc350_divergences.append(D)
        mc350_settle_z.append(settle_res['z_score'])
        mc350_mono_z.append(mono_res['z_score'])

mc350_divergences = np.array(mc350_divergences)
mc350_settle_z = np.array(mc350_settle_z)
mc350_mono_z = np.array(mc350_mono_z)

mc350_mean_D = float(np.mean(mc350_divergences))
mc350_std_D = float(np.std(mc350_divergences))
mc350_median_D = float(np.median(mc350_divergences))
mc350_pct_positive = float(np.sum(mc350_divergences > 0) / len(mc350_divergences) * 100)
mc350_pct_gt2 = float(np.sum(mc350_divergences > 2) / len(mc350_divergences) * 100)

print(f"\n  Results across {len(mc350_divergences)} valid iterations:")
print(f"  Mean D: {mc350_mean_D:.2f} ± {mc350_std_D:.2f}")
print(f"  Median D: {mc350_median_D:.2f}")
print(f"  D > 0 in {mc350_pct_positive:.1f}% of iterations")
print(f"  D > 2 in {mc350_pct_gt2:.1f}% of iterations")

# ============================================================
# SUMMARY
# ============================================================
elapsed = time.time() - t0
print(f"\n{'=' * 70}")
print(f"EGYPT PRESERVATION TEST — 1000 MC TRIALS SUMMARY ({elapsed:.0f}s)")
print(f"{'=' * 70}")
print(f"  Test 1 (Baseline):       D = {baseline_D}")
print(f"  Test 2 (+58 buried):     D = {aug_D}")
print(f"  Test 3 (+350 estimated): D = {mc350_mean_D:.2f} ± {mc350_std_D:.2f}")
print(f"                           D > 2 in {mc350_pct_gt2:.1f}% of iterations")

change_56 = aug_D - baseline_D if isinstance(baseline_D, (int, float)) and isinstance(aug_D, (int, float)) else "N/A"
survives_56 = isinstance(aug_D, (int, float)) and aug_D > 2
survives_350 = mc350_pct_gt2 > 50

# ============================================================
# SAVE RESULTS
# ============================================================
output = {
    "meta": {
        "date": "2026-03-19",
        "description": "Egypt preservation test — RE-RUN with 1000 MC trials (was 200). Matches main paper methodology.",
        "upgrade_note": "N_MC_TRIALS increased from 200 to 1000 per Perplexity methodological recommendation",
        "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
        "threshold_km": THRESHOLD_KM,
        "n_mc_trials": N_MC_TRIALS,
        "n_monte_carlo_iterations": N_MONTE_CARLO_ITER,
    },
    "buried_sites": {
        "n_valley": len(VALLEY_SITES),
        "n_delta": len(DELTA_SITES),
        "n_total": len(ALL_BURIED),
        "n_within_50km_gc": n_within_50,
        "mean_gc_distance_km": round(float(mean_dist), 1),
    },
    "test1_baseline": {
        "label": "Pleiades Egypt as-is (1000 MC)",
        "n_monumental": len(pleiades_mono),
        "n_settlement": len(pleiades_settle),
        "monumental": baseline_mono,
        "settlement": baseline_settle,
        "divergence_D": baseline_D
    },
    "test2_plus_58": {
        "label": "+58 named buried settlement sites (1000 MC)",
        "n_monumental": len(pleiades_mono),
        "n_settlement": len(augmented_settle),
        "n_added": len(buried_with_dist),
        "monumental": aug_mono,
        "settlement": aug_settle,
        "divergence_D": aug_D,
        "change_from_baseline": change_56
    },
    "test3_plus_350_mc": {
        "label": "+350 random floodplain settlements (1000 MC trials per Z-score)",
        "n_iterations": N_MONTE_CARLO_ITER,
        "n_valid_iterations": len(mc350_divergences),
        "n_random_settlements_per_iter": 350,
        "mean_D": round(mc350_mean_D, 2),
        "std_D": round(mc350_std_D, 2),
        "median_D": round(mc350_median_D, 2),
        "pct_D_positive": round(mc350_pct_positive, 1),
        "pct_D_gt_2": round(mc350_pct_gt2, 1),
        "mean_settlement_z": round(float(np.mean(mc350_settle_z)), 2),
        "mean_monument_z": round(float(np.mean(mc350_mono_z)), 2),
        "D_percentiles": {
            "5th": round(float(np.percentile(mc350_divergences, 5)), 2),
            "25th": round(float(np.percentile(mc350_divergences, 25)), 2),
            "50th": round(float(np.percentile(mc350_divergences, 50)), 2),
            "75th": round(float(np.percentile(mc350_divergences, 75)), 2),
            "95th": round(float(np.percentile(mc350_divergences, 95)), 2),
        }
    },
    "comparison_with_200mc": {
        "note": "Compare with egypt_preservation_test.json (200 MC trials) — results should be consistent but with tighter confidence intervals",
        "original_baseline_D": "see egypt_preservation_test.json",
        "original_mc_trials": 200,
        "upgraded_mc_trials": 1000,
    },
    "verdict": {
        "survives_plus_58": survives_56,
        "survives_plus_350": survives_350,
        "interpretation": (
            "Divergence SURVIVES preservation correction (1000 MC)" if (survives_56 and survives_350)
            else "Divergence WEAKENED but present" if (survives_56 or survives_350)
            else "Divergence DOES NOT SURVIVE preservation correction"
        )
    }
}

out_path = os.path.join(OUT_DIR, "egypt_preservation_test_1000mc.json")
with open(out_path, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved: {out_path}")

print("\nDONE.")
