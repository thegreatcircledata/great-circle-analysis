#!/usr/bin/env python3
"""
Statistical Polish Directive
==============================
Two standard spatial statistics analyses for PLOS ONE submission:
1. Ripley's K-function at multiple scales (Megalithic Portal, Pleiades, p3k14c)
2. Benjamini-Hochberg FDR correction across all tests
3. Summary statistics for paper text
"""

import csv, math, json, os, time, re
import numpy as np
from collections import Counter, defaultdict
from scipy import stats as sp_stats
import xml.etree.ElementTree as ET
import glob as glob_mod

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
N_TRIALS = 1000

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "statistical_polish")
os.makedirs(OUT_DIR, exist_ok=True)

R_VALUES = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250, 300, 400, 500]

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

def save_json(data, filename):
    path = os.path.join(OUT_DIR, filename)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  -> Saved: {path}")


# ============================================================
# K-FUNCTION ENGINE
# ============================================================
def k_function(site_lats, site_lons, r_values=R_VALUES, n_trials=N_TRIALS,
               label="dataset"):
    """Compute K-function analog for sites relative to the Great Circle.

    For each distance threshold r, counts sites within r km of the circle
    and compares to distribution-matched MC baseline.
    """
    n = len(site_lats)
    print(f"  Computing K-function for {label} ({n} sites, {len(r_values)} thresholds, {n_trials} trials)...")

    # Compute observed distances to GC
    dists = gc_dist_vec(POLE_LAT, POLE_LON, site_lats, site_lons)

    # Observed counts at each r
    observed = {r: int(np.sum(dists <= r)) for r in r_values}

    # MC baseline: for each trial, generate n random points and count at all r
    mc_counts = {r: [] for r in r_values}
    t0 = time.time()
    for trial in range(n_trials):
        rl, rn = rand_matched_batch(site_lats, site_lons, n)
        rd = gc_dist_vec(POLE_LAT, POLE_LON, rl, rn)
        for r in r_values:
            mc_counts[r].append(int(np.sum(rd <= r)))
        if (trial + 1) % 200 == 0:
            elapsed = time.time() - t0
            print(f"    {trial+1}/{n_trials} trials ({elapsed:.1f}s)")

    results = {}
    for r in r_values:
        obs = observed[r]
        mc = np.array(mc_counts[r])
        mu = float(np.mean(mc))
        sigma = float(np.std(mc))
        L_r = (obs / mu) - 1 if mu > 0 else 0
        Z_r = (obs - mu) / sigma if sigma > 0 else 0
        ci_lower = float(np.percentile(mc, 2.5))
        ci_upper = float(np.percentile(mc, 97.5))
        p_val = float(np.sum(mc >= obs) / n_trials)

        results[str(r)] = {
            "observed": int(obs),
            "expected": round(mu, 2),
            "std": round(sigma, 2),
            "L_r": round(L_r, 4),
            "Z_r": round(Z_r, 2),
            "ci_lower": round(ci_lower, 1),
            "ci_upper": round(ci_upper, 1),
            "p_value": round(p_val, 4),
        }

    # Find peak significance and significant range
    z_vals = {r: results[str(r)]["Z_r"] for r in r_values}
    peak_r = max(z_vals, key=z_vals.get)
    peak_z = z_vals[peak_r]
    sig_range = [r for r in r_values if z_vals[r] > 3]
    sig_min = min(sig_range) if sig_range else None
    sig_max = max(sig_range) if sig_range else None

    print(f"  Peak Z = {peak_z} at r = {peak_r} km")
    if sig_range:
        print(f"  Significant (Z>3) range: {sig_min}-{sig_max} km")
    else:
        print(f"  No distance shows Z > 3")

    return {
        "n_sites": int(n),
        "r_values": r_values,
        "results": results,
        "peak_significance_km": int(peak_r),
        "peak_z": round(peak_z, 2),
        "significant_range_km": [sig_min, sig_max] if sig_range else None,
    }


# ============================================================
# DATA LOADING
# ============================================================
print("=" * 70)
print("STATISTICAL POLISH DIRECTIVE")
print("=" * 70)

t_start = time.time()
np.random.seed(42)

# --- Megalithic Portal ---
print("\nLoading Megalithic Portal data...")
kml_pattern = os.path.join(BASE_DIR, "MegP_*.kml")
kml_files = sorted(glob_mod.glob(kml_pattern))
mp_sites = []
for path in kml_files:
    fname = os.path.basename(path)
    site_type = fname.replace("MegP_", "").replace(".kml", "").replace("_", " ")
    try:
        tree = ET.parse(path)
    except ET.ParseError:
        continue
    root = tree.getroot()
    # Try both KML namespaces
    ns_options = [
        {'kml': 'http://www.opengis.net/kml/2.2'},
        {'kml': 'http://earth.google.com/kml/2.0'},
        {'kml': 'http://earth.google.com/kml/2.1'},
    ]
    # Detect namespace from root tag
    root_tag = root.tag
    if 'earth.google.com/kml/2.0' in root_tag:
        ns = ns_options[1]
        ns_uri = 'http://earth.google.com/kml/2.0'
    elif 'earth.google.com/kml/2.1' in root_tag:
        ns = ns_options[2]
        ns_uri = 'http://earth.google.com/kml/2.1'
    else:
        ns = ns_options[0]
        ns_uri = 'http://www.opengis.net/kml/2.2'
    for pm in root.iter(f'{{{ns_uri}}}Placemark'):
        coords_el = pm.find(f'.//{{{ns_uri}}}coordinates')
        if coords_el is None or not coords_el.text:
            continue
        try:
            parts = coords_el.text.strip().split(",")
            lon, lat = float(parts[0]), float(parts[1])
            if -90 <= lat <= 90 and -180 <= lon <= 180:
                mp_sites.append({"lat": lat, "lon": lon, "type": site_type})
        except (ValueError, IndexError):
            continue

# Dedup
seen = set()
mp_dedup = []
for s in mp_sites:
    key = (round(s["lat"], 3), round(s["lon"], 3))
    if key not in seen:
        seen.add(key)
        mp_dedup.append(s)
mp_sites = mp_dedup
print(f"  Megalithic Portal: {len(mp_sites)} unique sites")

# --- Pleiades ---
print("Loading Pleiades data...")
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

pleiades_mon = []
pleiades_set = []
pleiades_path = os.path.join(BASE_DIR, "pleiades-places-latest.csv")
with open(pleiades_path, encoding="utf-8") as f:
    for row in csv.DictReader(f):
        try:
            lat = float(row["reprLat"]); lon = float(row["reprLong"])
        except (ValueError, KeyError):
            continue
        if lat == 0 and lon == 0: continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180): continue
        feature_types = {t.strip() for t in row.get("featureTypes", "").split(",") if t.strip()}
        # Filter to ancient (pre-2000 BCE where available)
        min_date = None
        try:
            min_date = int(row.get("minDate", ""))
        except (ValueError, TypeError):
            pass
        if min_date is not None and min_date > -2000:
            continue
        if feature_types & PLEIADES_MONUMENTAL:
            pleiades_mon.append((lat, lon))
        elif feature_types & PLEIADES_SETTLEMENT:
            pleiades_set.append((lat, lon))
print(f"  Pleiades ancient: {len(pleiades_mon)} monumental, {len(pleiades_set)} settlement")

# --- p3k14c ---
print("Loading p3k14c data...")
p3k_all = []
with open(os.path.join(BASE_DIR, "p3k14c_data.csv"), encoding="utf-8") as f:
    for row in csv.DictReader(f):
        try:
            lat, lon = float(row["Lat"]), float(row["Long"])
            if lat == 0 and lon == 0: continue
            if not (-90 <= lat <= 90 and -180 <= lon <= 180): continue
            p3k_all.append({"lat": lat, "lon": lon,
                            "site_id": row.get("SiteID", ""),
                            "site_name": row.get("SiteName", "")})
        except: pass
# Dedup by SiteID
p3k_groups = defaultdict(list)
for r in p3k_all:
    key = r["site_id"] if r["site_id"] else f"anon_{r['lat']:.4f}_{r['lon']:.4f}"
    p3k_groups[key].append(r)
p3k_sites = []
for sid, rows in p3k_groups.items():
    p3k_sites.append({
        "lat": np.mean([r["lat"] for r in rows]),
        "lon": np.mean([r["lon"] for r in rows]),
    })
print(f"  p3k14c: {len(p3k_sites)} unique sites")


# ============================================================
# TEST 2: BENJAMINI-HOCHBERG FDR CORRECTION (Priority 1)
# ============================================================
print(f"\n{'=' * 70}")
print("TEST 2: BENJAMINI-HOCHBERG FDR CORRECTION")
print("=" * 70)

def z_to_p(z, two_tailed=False):
    """Convert Z-score to p-value (one-tailed by default)."""
    p = 1 - sp_stats.norm.cdf(abs(z))
    if two_tailed:
        p *= 2
    return min(p, 1.0)

# Collect all tests with p-values
# Format: (name, category, p_value, source)
tests = [
    # PRIMARY FINDINGS
    ("MP overall (dist-matched, 50km)", "Primary", 0.000999, "peer_review_response"),
    ("Pleiades ancient monuments (50km)", "Primary", 0.000999, "peer_review_response"),
    ("Pleiades ancient settlements (50km)", "Primary", 0.9371, "peer_review_response"),
    ("p3k14c overall (50km)", "Primary", z_to_p(7.21), "nazca_solution_validation"),
    ("Habitability percentile", "Primary", 0.43, "peer_review_response"),
    ("10,000-circle divergence (0/10,000)", "Primary", 0.0001, "divergence_10k"),

    # ROBUSTNESS / VALIDATION
    ("p3k14c LocAccuracy>=2 (50km)", "Robustness", z_to_p(7.86), "peer_review_response"),
    ("p3k14c without SA (50km)", "Robustness", z_to_p(5.09), "p3k14c_sampling_control"),
    ("Split-sample validation (100% > Z=5)", "Validation", z_to_p(9.45), "split_sample_validation"),
    ("Leave-one-out: Egypt removed", "Validation", z_to_p(4.51), "spatial_block_validation"),
    ("Leave-one-out: Peru removed", "Validation", z_to_p(10.41), "spatial_block_validation"),
    ("Leave-one-out: Easter Island removed", "Validation", z_to_p(12.55), "spatial_block_validation"),

    # CONTINENTAL DECOMPOSITION
    ("p3k14c Africa (50km)", "Continental", z_to_p(3.53), "hemisphere_decomposition"),
    ("p3k14c Asia (50km)", "Continental", 0.05310, "hemisphere_decomposition"),
    ("p3k14c South America (50km)", "Continental", z_to_p(9.25), "hemisphere_decomposition"),
    ("p3k14c Australia (50km)", "Continental", 1.0, "hemisphere_decomposition"),
    ("p3k14c Europe (50km)", "Continental", 1.0, "hemisphere_decomposition"),
    ("p3k14c North America (50km)", "Continental", 1.0, "hemisphere_decomposition"),

    # NEW DATABASE VALIDATION
    ("DARE overall (50km)", "Validation", z_to_p(4.77), "new_databases"),
    ("DARE monuments (50km)", "Validation", 0.329, "new_databases"),
    ("DARE settlements (50km)", "Validation", 1.0, "new_databases"),
    ("Historic England (50km)", "Control", 1.0, "new_databases"),

    # MONUMENT-SETTLEMENT DIVERGENCE
    ("p3k14c monument Z (50km)", "Divergence", z_to_p(8.09), "continental_divergence"),
    ("p3k14c domestic Z (50km)", "Divergence", z_to_p(1.91), "continental_divergence"),
    ("p3k14c divergence without Africa", "Divergence", z_to_p(4.59), "continental_divergence"),
    ("p3k14c divergence without SA", "Divergence", z_to_p(3.23), "continental_divergence"),
    ("DARE divergence (monumental - settlement)", "Divergence", z_to_p(6.15), "new_databases"),
    ("SA-only divergence (50km)", "Divergence", z_to_p(5.68), "continental_divergence"),

    # TEMPORAL BINS (5 bins)
    ("Temporal: >10000 BP", "Exploratory", z_to_p(0.91), "nazca_solution_validation"),
    ("Temporal: 5000-10000 BP", "Exploratory", 0.01265, "nazca_solution_validation"),
    ("Temporal: 3000-5000 BP", "Exploratory", z_to_p(4.05), "nazca_solution_validation"),
    ("Temporal: 1000-3000 BP", "Exploratory", 0.001795, "nazca_solution_validation"),
    ("Temporal: <1000 BP", "Exploratory", z_to_p(9.03), "nazca_solution_validation"),

    # TYPE ENRICHMENT (selected types with enough data)
    ("Type: Geoglyph enrichment", "Exploratory", z_to_p(5.84), "hemisphere_decomposition"),
    ("Type: Sculptured Stone enrichment", "Exploratory", z_to_p(6.39), "hemisphere_decomposition"),
    ("Type: Holy Well enrichment", "Exploratory", 0.00216, "hemisphere_decomposition"),
    ("Type: Temple enrichment", "Exploratory", 0.00917, "hemisphere_decomposition"),
    ("Type: Standing Stones enrichment", "Exploratory", 0.04876, "hemisphere_decomposition"),

    # STRATIFIED MC
    ("Stratified MC (alison)", "Robustness", z_to_p(4.42), "data_quality_checks"),

    # OLD WORLD / NEW WORLD
    ("Old World non-Euro monumental Z", "Exploratory", z_to_p(11.03), "hemisphere_decomposition"),
    ("New World monumental Z", "Exploratory", z_to_p(8.1), "hemisphere_decomposition"),
    ("New World settlement Z", "Exploratory", z_to_p(7.97), "hemisphere_decomposition"),
]

# Extract p-values
test_names = [t[0] for t in tests]
categories = [t[1] for t in tests]
p_values = np.array([t[2] for t in tests])

# Replace exact 0 with a small value (can't take log of 0)
p_values = np.maximum(p_values, 1e-16)

# Apply Benjamini-Hochberg correction manually (avoid statsmodels dependency)
n_tests = len(p_values)
sorted_indices = np.argsort(p_values)
sorted_p = p_values[sorted_indices]

# BH adjusted p-values
bh_adjusted = np.zeros(n_tests)
for i in range(n_tests):
    rank = i + 1
    bh_adjusted[i] = sorted_p[i] * n_tests / rank

# Enforce monotonicity (step-up)
for i in range(n_tests - 2, -1, -1):
    bh_adjusted[i] = min(bh_adjusted[i], bh_adjusted[i + 1])
bh_adjusted = np.minimum(bh_adjusted, 1.0)

# Map back to original order
p_adjusted = np.zeros(n_tests)
for i, idx in enumerate(sorted_indices):
    p_adjusted[idx] = bh_adjusted[i]

# Determine significance
alpha = 0.05
significant = p_adjusted <= alpha

# Report
print(f"\n  Total tests: {n_tests}")
print(f"  Significant before correction: {np.sum(p_values <= alpha)}")
print(f"  Significant after BH correction (q=0.05): {int(np.sum(significant))}")

print(f"\n  {'Test':<50} {'Cat':<12} {'p_orig':>10} {'p_adj':>10} {'Sig':>5}")
print(f"  {'-'*50} {'-'*12} {'-'*10} {'-'*10} {'-'*5}")
for i in range(n_tests):
    sig_marker = "YES" if significant[i] else "no"
    p_o = f"{p_values[i]:.6f}" if p_values[i] > 1e-10 else f"{p_values[i]:.2e}"
    p_a = f"{p_adjusted[i]:.6f}" if p_adjusted[i] > 1e-10 else f"{p_adjusted[i]:.2e}"
    print(f"  {test_names[i]:<50} {categories[i]:<12} {p_o:>10} {p_a:>10} {sig_marker:>5}")

# Find what was lost
lost = []
for i in range(n_tests):
    if p_values[i] <= alpha and not significant[i]:
        lost.append(test_names[i])

# Primary findings survival check
primary_tests = [i for i in range(n_tests) if categories[i] == "Primary"]
primary_survive = all(significant[i] or p_values[i] > alpha for i in primary_tests
                      if p_values[i] <= alpha)

# Divergence tests survival
div_tests = [i for i in range(n_tests) if categories[i] == "Divergence"]
div_survive = all(significant[i] or p_values[i] > alpha for i in div_tests
                  if p_values[i] <= alpha)

print(f"\n  Lost after correction: {lost if lost else 'None'}")
print(f"  Primary findings survive: {primary_survive}")
print(f"  Divergence findings survive: {div_survive}")

bh_results = {
    "meta": {"date": "2026-03-18", "method": "Benjamini-Hochberg FDR",
             "alpha": 0.05, "total_tests": int(n_tests)},
    "tests": [
        {"name": test_names[i], "category": categories[i],
         "p_original": round(float(p_values[i]), 10),
         "p_adjusted": round(float(p_adjusted[i]), 10),
         "significant_at_005": bool(significant[i])}
        for i in range(n_tests)
    ],
    "n_significant_original": int(np.sum(p_values <= alpha)),
    "n_significant_corrected": int(np.sum(significant)),
    "primary_findings_survive": bool(primary_survive),
    "divergence_findings_survive": bool(div_survive),
    "lost_after_correction": lost,
    "category_summary": {
        cat: {
            "n_tests": int(np.sum(np.array(categories) == cat)),
            "n_significant": int(np.sum(significant[np.array(categories) == cat])),
        }
        for cat in sorted(set(categories))
    },
}
save_json(bh_results, "benjamini_hochberg.json")


# ============================================================
# TEST 1A: K-FUNCTION FOR MEGALITHIC PORTAL
# ============================================================
print(f"\n{'=' * 70}")
print("TEST 1A: K-FUNCTION — MEGALITHIC PORTAL")
print("=" * 70)

mp_lats = np.array([s["lat"] for s in mp_sites])
mp_lons = np.array([s["lon"] for s in mp_sites])
mp_kfunc = k_function(mp_lats, mp_lons, R_VALUES, N_TRIALS, "Megalithic Portal")

save_json({
    "meta": {"date": "2026-03-18", "dataset": "Megalithic Portal",
             "n_trials": N_TRIALS},
    **mp_kfunc,
}, "k_function_mp.json")


# ============================================================
# TEST 1B: K-FUNCTION FOR PLEIADES (MONUMENTS vs SETTLEMENTS)
# ============================================================
print(f"\n{'=' * 70}")
print("TEST 1B: K-FUNCTION — PLEIADES MONUMENTS vs SETTLEMENTS")
print("=" * 70)

mon_lats = np.array([s[0] for s in pleiades_mon])
mon_lons = np.array([s[1] for s in pleiades_mon])
set_lats = np.array([s[0] for s in pleiades_set])
set_lons = np.array([s[1] for s in pleiades_set])

print("\n  --- Monuments ---")
pleiades_mon_kfunc = k_function(mon_lats, mon_lons, R_VALUES, N_TRIALS, "Pleiades Monuments")

print("\n  --- Settlements ---")
pleiades_set_kfunc = k_function(set_lats, set_lons, R_VALUES, N_TRIALS, "Pleiades Settlements")

# Compute divergence at each scale
divergence_by_r = {}
for r in R_VALUES:
    mon_z = pleiades_mon_kfunc["results"][str(r)]["Z_r"]
    set_z = pleiades_set_kfunc["results"][str(r)]["Z_r"]
    divergence_by_r[str(r)] = round(mon_z - set_z, 2)

save_json({
    "meta": {"date": "2026-03-18", "dataset": "Pleiades (ancient, pre-2000 BCE)",
             "n_trials": N_TRIALS},
    "monuments": pleiades_mon_kfunc,
    "settlements": pleiades_set_kfunc,
    "divergence_by_r": divergence_by_r,
    "peak_divergence_km": int(max(divergence_by_r, key=lambda r: divergence_by_r[r])),
    "peak_divergence_z": max(divergence_by_r.values()),
}, "k_function_pleiades.json")


# ============================================================
# TEST 1C: K-FUNCTION FOR P3K14C
# ============================================================
print(f"\n{'=' * 70}")
print("TEST 1C: K-FUNCTION — P3K14C")
print("=" * 70)

p3k_lats = np.array([s["lat"] for s in p3k_sites])
p3k_lons = np.array([s["lon"] for s in p3k_sites])
p3k_kfunc = k_function(p3k_lats, p3k_lons, R_VALUES, N_TRIALS, "p3k14c")

save_json({
    "meta": {"date": "2026-03-18", "dataset": "p3k14c",
             "n_trials": N_TRIALS},
    **p3k_kfunc,
}, "k_function_p3k14c.json")


# ============================================================
# FINAL SUMMARY
# ============================================================
elapsed = time.time() - t_start
print(f"\n{'=' * 70}")
print(f"ALL ANALYSES COMPLETE — {elapsed:.0f}s elapsed")
print(f"{'=' * 70}")
print(f"Output directory: {OUT_DIR}")

print(f"\n  K-FUNCTION PEAKS:")
print(f"    Megalithic Portal: peak Z={mp_kfunc['peak_z']} at {mp_kfunc['peak_significance_km']}km, "
      f"significant range: {mp_kfunc['significant_range_km']}")
print(f"    Pleiades Monuments: peak Z={pleiades_mon_kfunc['peak_z']} at {pleiades_mon_kfunc['peak_significance_km']}km, "
      f"significant range: {pleiades_mon_kfunc['significant_range_km']}")
print(f"    Pleiades Settlements: peak Z={pleiades_set_kfunc['peak_z']} at {pleiades_set_kfunc['peak_significance_km']}km, "
      f"significant range: {pleiades_set_kfunc['significant_range_km']}")
print(f"    p3k14c: peak Z={p3k_kfunc['peak_z']} at {p3k_kfunc['peak_significance_km']}km, "
      f"significant range: {p3k_kfunc['significant_range_km']}")

peak_div_r = max(divergence_by_r, key=lambda r: divergence_by_r[r])
print(f"\n  PLEIADES DIVERGENCE PEAK: {divergence_by_r[peak_div_r]} Z-units at {peak_div_r}km")

print(f"\n  BH CORRECTION SUMMARY:")
print(f"    {bh_results['n_significant_original']}/{n_tests} significant before correction")
print(f"    {bh_results['n_significant_corrected']}/{n_tests} significant after BH (q=0.05)")
print(f"    Primary findings survive: {bh_results['primary_findings_survive']}")
print(f"    Lost after correction: {bh_results['lost_after_correction']}")
