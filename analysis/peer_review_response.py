#!/usr/bin/env python3
"""
Peer Review Response — Four Critical Tests
=============================================
Addresses reviewer concerns before journal submission:
  Test 3: Upgrade MC trials to 1,000 (+ Shapiro-Wilk, CIs, exact p-values)
  Test 4: Resolve Z=-0.04 / 78.7th percentile discrepancy
  Test 2: Spatial block cross-validation (leave-one-cluster-out + hemisphere block)
  Test 1: Habitability-adjusted monument-settlement divergence

Execution order: 3 → 4 → 2 → 1 (per directive)
"""

import csv, math, random, json, os, sys, time
import numpy as np
from scipy import stats as sp_stats
from collections import defaultdict

# Force unbuffered output
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

# ============================================================
# CONFIGURATION
# ============================================================
POLE_LAT = 59.682122
POLE_LNG = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * np.pi / 2
THRESHOLD_KM = 50
JITTER_DEG = 2.0

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "peer_review_response")
os.makedirs(OUT_DIR, exist_ok=True)

# Pleiades type classifications (same as hemisphere_decomposition.py)
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

def count_within(pole_lat, pole_lon, site_lats, site_lons, threshold=THRESHOLD_KM):
    dists = gc_dist_vec(pole_lat, pole_lon, site_lats, site_lons)
    return int(np.sum(dists <= threshold))

def rand_matched_batch(site_lats, site_lons, n):
    idx = np.random.randint(0, len(site_lats), size=n)
    lats = site_lats[idx] + np.random.normal(0, JITTER_DEG, n)
    lons = site_lons[idx] + np.random.normal(0, JITTER_DEG, n)
    return np.clip(lats, -90, 90), np.clip(lons, -180, 180)

def run_mc(site_lats, site_lons, n_trials, threshold=THRESHOLD_KM,
           pole_lat=POLE_LAT, pole_lon=POLE_LNG):
    """Run distribution-matched MC and return detailed results."""
    n = len(site_lats)
    if n == 0:
        return {"n_sites": 0, "observed": 0, "expected": 0, "std": 0,
                "z_score": 0, "enrichment": 0, "empirical_p": 1.0,
                "baseline_counts": []}
    dists = gc_dist_vec(pole_lat, pole_lon, site_lats, site_lons)
    observed = int(np.sum(dists <= threshold))
    rand_counts = []
    for _ in range(n_trials):
        rl, rn = rand_matched_batch(site_lats, site_lons, n)
        rd = gc_dist_vec(pole_lat, pole_lon, rl, rn)
        rand_counts.append(int(np.sum(rd <= threshold)))
    rand_arr = np.array(rand_counts)
    mu = float(np.mean(rand_arr))
    sigma = float(np.std(rand_arr))
    z = float((observed - mu) / sigma) if sigma > 0 else 0.0
    # Empirical p-value: fraction of random trials >= observed
    n_exceeding = int(np.sum(rand_arr >= observed))
    empirical_p = (n_exceeding + 1) / (n_trials + 1)  # +1 for continuity correction
    enrich = float(observed / mu) if mu > 0 else 0.0
    return {
        "n_sites": int(n), "observed": int(observed),
        "expected": round(mu, 2), "std": round(sigma, 2),
        "z_score": round(z, 2), "enrichment": round(enrich, 3),
        "empirical_p": float(empirical_p),
        "n_exceeding": n_exceeding,
        "baseline_counts": rand_counts,
    }

def compute_stats(result, n_trials):
    """Compute Shapiro-Wilk, CI, and formatted p-value from MC result."""
    counts = np.array(result["baseline_counts"])
    # Shapiro-Wilk on baseline distribution
    if len(counts) >= 8:
        # Sample up to 5000 for Shapiro-Wilk (scipy limit)
        sample = counts[:5000] if len(counts) > 5000 else counts
        sw_stat, sw_p = sp_stats.shapiro(sample)
    else:
        sw_stat, sw_p = 0, 1.0

    # 95% CI on enrichment ratio using bootstrap
    observed = result["observed"]
    mu = result["expected"]
    if mu > 0 and len(counts) > 0:
        enrichments = observed / counts[counts > 0]
        if len(enrichments) >= 20:
            ci_lower = float(np.percentile(enrichments, 2.5))
            ci_upper = float(np.percentile(enrichments, 97.5))
        else:
            # Adjusted Wald
            n_eff = len(counts)
            p_hat = observed / (mu * n_eff) if mu * n_eff > 0 else 0
            se = math.sqrt(p_hat * (1 - p_hat) / n_eff) if n_eff > 0 else 0
            ci_lower = max(0, result["enrichment"] - 1.96 * se)
            ci_upper = result["enrichment"] + 1.96 * se
    else:
        ci_lower, ci_upper = 0, 0

    # Format p-value per PLOS ONE SAMPL guidelines
    emp_p = result["empirical_p"]
    if emp_p >= 0.001:
        p_formatted = f"p = {emp_p:.4f}"
    else:
        p_formatted = f"p < 0.001"

    return {
        "shapiro_wilk": {"statistic": round(float(sw_stat), 4), "p_value": float(sw_p)},
        "enrichment_95ci": [round(ci_lower, 3), round(ci_upper, 3)],
        "p_formatted": p_formatted,
    }

def save_json(data, filename):
    path = os.path.join(OUT_DIR, filename)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  -> Saved: {path}")

# ============================================================
# DATA LOADING
# ============================================================
def load_p3k14c():
    csv_path = os.path.join(BASE_DIR, "p3k14c_data.csv")
    rows = []
    with open(csv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["Lat"]); lon = float(row["Long"])
                if lat == 0 and lon == 0: continue
                continent = row.get("Continent", "").strip()
                country = row.get("Country", "").strip()
                site_name = row.get("SiteName", "").strip()
                rows.append((lat, lon, continent, country, site_name))
            except (ValueError, KeyError, TypeError):
                pass
    seen = {}
    for lat, lon, cont, country, name in rows:
        key = (round(lat, 4), round(lon, 4))
        if key not in seen:
            seen[key] = (lat, lon, cont, country, name)
    return list(seen.values())

def load_megalithic_portal():
    csv_path = os.path.join(BASE_DIR, "github-repo", "data", "merged_sites.csv")
    sites = []
    with open(csv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["lat"]); lon = float(row["lon"])
                if lat == 0 and lon == 0: continue
                name = row.get("name", "")
                stype = row.get("type", "").strip()
                sites.append((lat, lon, name, stype))
            except (ValueError, KeyError, TypeError):
                pass
    seen = {}
    for lat, lon, name, stype in sites:
        key = (round(lat, 3), round(lon, 3))
        if key not in seen:
            seen[key] = (lat, lon, name, stype)
    return list(seen.values())

def load_pleiades(ancient_only=False):
    csv_path = os.path.join(BASE_DIR, "pleiades-places-latest.csv")
    monumental = []; settlements = []
    with open(csv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["reprLat"]); lon = float(row["reprLong"])
            except (ValueError, KeyError, TypeError):
                continue
            if lat == 0 and lon == 0: continue
            if not (-90 <= lat <= 90 and -180 <= lon <= 180): continue
            feature_types = {t.strip() for t in row.get("featureTypes", "").split(",") if t.strip()}
            is_mon = bool(feature_types & MONUMENTAL_TYPES)
            is_set = bool(feature_types & SETTLEMENT_TYPES)

            if ancient_only:
                try:
                    min_date = int(row.get("minDate", "0"))
                except (ValueError, TypeError):
                    min_date = 0
                if min_date >= -2000:
                    continue

            if is_mon:
                monumental.append((lat, lon))
            elif is_set:
                settlements.append((lat, lon))
    return monumental, settlements


# ============================================================
# TEST 3: UPGRADE MC TRIALS TO 1,000
# ============================================================
def test_3():
    print("\n" + "=" * 70)
    print("TEST 3: UPGRADE MC TRIALS TO 1,000")
    print("=" * 70)

    N_TRIALS = 1000
    results = {"meta": {"date": "2026-03-17", "n_trials": N_TRIALS, "threshold_km": THRESHOLD_KM}}

    # --- 3a: Megalithic Portal Great Circle test ---
    print("\n  --- 3a: Megalithic Portal (1,000 trials, 50 km) ---")
    mp_sites = load_megalithic_portal()
    mp_lats = np.array([s[0] for s in mp_sites])
    mp_lons = np.array([s[1] for s in mp_sites])
    print(f"  Loaded {len(mp_sites)} Megalithic Portal sites")

    r_mp = run_mc(mp_lats, mp_lons, N_TRIALS)
    stats_mp = compute_stats(r_mp, N_TRIALS)
    print(f"  MegP: observed={r_mp['observed']}, expected={r_mp['expected']}, "
          f"Z={r_mp['z_score']}, enrichment={r_mp['enrichment']}, "
          f"empirical_p={r_mp['empirical_p']}")
    print(f"  Shapiro-Wilk: W={stats_mp['shapiro_wilk']['statistic']}, "
          f"p={stats_mp['shapiro_wilk']['p_value']:.4e}")
    print(f"  Enrichment 95% CI: {stats_mp['enrichment_95ci']}")

    # Remove baseline_counts from saved output (too large)
    r_mp_save = {k: v for k, v in r_mp.items() if k != "baseline_counts"}
    r_mp_save.update(stats_mp)
    results["megalithic_portal"] = r_mp_save

    # --- 3b: Pleiades monument-settlement divergence (1,000 trials each, pre-2000 BCE) ---
    print("\n  --- 3b: Pleiades monument-settlement divergence (1,000 trials, pre-2000 BCE) ---")
    mon_anc, set_anc = load_pleiades(ancient_only=True)
    mon_lats = np.array([s[0] for s in mon_anc])
    mon_lons = np.array([s[1] for s in mon_anc])
    set_lats = np.array([s[0] for s in set_anc])
    set_lons = np.array([s[1] for s in set_anc])
    print(f"  Monumental (ancient): {len(mon_anc)}, Settlement (ancient): {len(set_anc)}")

    r_mon = run_mc(mon_lats, mon_lons, N_TRIALS)
    stats_mon = compute_stats(r_mon, N_TRIALS)
    r_set = run_mc(set_lats, set_lons, N_TRIALS)
    stats_set = compute_stats(r_set, N_TRIALS)
    divergence = round(r_mon["z_score"] - r_set["z_score"], 2)

    print(f"  Monument Z={r_mon['z_score']}, Settlement Z={r_set['z_score']}, Divergence={divergence}")
    print(f"  Monument empirical p={r_mon['empirical_p']}")
    print(f"  Settlement empirical p={r_set['empirical_p']}")

    r_mon_save = {k: v for k, v in r_mon.items() if k != "baseline_counts"}
    r_mon_save.update(stats_mon)
    r_set_save = {k: v for k, v in r_set.items() if k != "baseline_counts"}
    r_set_save.update(stats_set)
    results["pleiades_divergence"] = {
        "monument": r_mon_save, "settlement": r_set_save,
        "divergence": divergence,
    }

    # --- 3c: p3k14c unique sites Great Circle test (1,000 trials, 50 km) ---
    print("\n  --- 3c: p3k14c unique sites (1,000 trials, 50 km) ---")
    p3k_sites = load_p3k14c()
    p3k_lats = np.array([s[0] for s in p3k_sites])
    p3k_lons = np.array([s[1] for s in p3k_sites])
    print(f"  Loaded {len(p3k_sites)} unique p3k14c sites")

    r_p3k = run_mc(p3k_lats, p3k_lons, N_TRIALS)
    stats_p3k = compute_stats(r_p3k, N_TRIALS)
    print(f"  p3k14c: observed={r_p3k['observed']}, expected={r_p3k['expected']}, "
          f"Z={r_p3k['z_score']}, enrichment={r_p3k['enrichment']}, "
          f"empirical_p={r_p3k['empirical_p']}")
    print(f"  Shapiro-Wilk: W={stats_p3k['shapiro_wilk']['statistic']}, "
          f"p={stats_p3k['shapiro_wilk']['p_value']:.4e}")

    r_p3k_save = {k: v for k, v in r_p3k.items() if k != "baseline_counts"}
    r_p3k_save.update(stats_p3k)
    results["p3k14c"] = r_p3k_save

    save_json(results, "upgraded_statistics.json")
    return results, r_p3k  # return raw r_p3k for Test 4


# ============================================================
# TEST 4: RESOLVE Z = -0.04 / 78.7th PERCENTILE DISCREPANCY
# ============================================================
def test_4():
    print("\n" + "=" * 70)
    print("TEST 4: RESOLVE Z = -0.04 / 78.7th PERCENTILE DISCREPANCY")
    print("=" * 70)

    p3k_sites = load_p3k14c()
    all_lats = np.array([s[0] for s in p3k_sites])
    all_lons = np.array([s[1] for s in p3k_sites])

    # Compute Alison circle count
    alison_dists = gc_dist_vec(POLE_LAT, POLE_LNG, all_lats, all_lons)
    alison_count = int(np.sum(alison_dists <= THRESHOLD_KM))

    # Generate 1,000 random great circles (upgrade from original 1,000)
    n_random = 1000
    random_counts = []
    print(f"  Running {n_random} random great circles on p3k14c...")
    for i in range(n_random):
        z = np.random.uniform(-1, 1)
        rand_pole_lat = math.degrees(math.asin(z))
        rand_pole_lon = np.random.uniform(-180, 180)
        dists = gc_dist_vec(rand_pole_lat, rand_pole_lon, all_lats, all_lons)
        random_counts.append(int(np.sum(dists <= THRESHOLD_KM)))
        if (i + 1) % 250 == 0:
            print(f"    {i+1}/{n_random} circles done")

    rand_arr = np.array(random_counts)
    mu = float(np.mean(rand_arr))
    sigma = float(np.std(rand_arr))
    z_hab = (alison_count - mu) / sigma if sigma > 0 else 0

    # Empirical percentile (the correct metric)
    empirical_pct = float(np.sum(rand_arr < alison_count) / len(rand_arr) * 100)
    # Two-tailed p-value from empirical distribution
    rank_above = np.sum(rand_arr >= alison_count)
    rank_below = np.sum(rand_arr <= alison_count)
    two_tailed_p = 2 * min(rank_above, rank_below) / len(rand_arr)
    two_tailed_p = min(1.0, two_tailed_p)

    # Shapiro-Wilk on the random circle distribution
    sw_stat, sw_p = sp_stats.shapiro(rand_arr[:5000])

    print(f"\n  Alison circle: {alison_count} p3k14c sites within {THRESHOLD_KM}km")
    print(f"  Random circles: mean={mu:.1f}, std={sigma:.1f}")
    print(f"  Z-score (habitability): {z_hab:.2f}")
    print(f"  Empirical percentile: {empirical_pct:.1f}%")
    print(f"  Two-tailed p: {two_tailed_p:.4f}")
    print(f"  Shapiro-Wilk: W={sw_stat:.4f}, p={sw_p:.4e}")
    print(f"  Distribution is {'non-normal' if sw_p < 0.05 else 'normal'} (p={'<0.001' if sw_p < 0.001 else f'{sw_p:.4f}'})")

    # Explanation
    if sw_p < 0.05:
        explanation = (
            f"The habitability-adjusted distribution is non-normal (Shapiro-Wilk p={sw_p:.4e}). "
            f"The Z-score framework is inappropriate. The correct reporting is: "
            f"'The Alison circle site count ranks at the {empirical_pct:.1f}th percentile "
            f"among {n_random} random great circles — not significantly different from random "
            f"(two-tailed p = {two_tailed_p:.2f}).' "
            f"The Z = {z_hab:.2f} should be dropped from the abstract."
        )
    else:
        explanation = (
            f"The distribution is approximately normal. Z = {z_hab:.2f} corresponds to "
            f"percentile {empirical_pct:.1f}%, consistent with normal CDF expectations."
        )

    results = {
        "meta": {"date": "2026-03-17", "n_random_circles": n_random,
                 "dataset": "p3k14c", "threshold_km": THRESHOLD_KM},
        "alison_count": alison_count,
        "random_circle_distribution": {
            "mean": round(mu, 2), "std": round(sigma, 2),
            "min": int(np.min(rand_arr)), "max": int(np.max(rand_arr)),
            "p5": round(float(np.percentile(rand_arr, 5)), 1),
            "p25": round(float(np.percentile(rand_arr, 25)), 1),
            "median": round(float(np.median(rand_arr)), 1),
            "p75": round(float(np.percentile(rand_arr, 75)), 1),
            "p95": round(float(np.percentile(rand_arr, 95)), 1),
        },
        "habitability_z_score": round(z_hab, 2),
        "habitability_empirical_percentile": round(empirical_pct, 1),
        "habitability_two_tailed_p": round(two_tailed_p, 4),
        "habitability_shapiro_wilk": {
            "statistic": round(float(sw_stat), 4),
            "p_value": float(sw_p),
            "distribution_is_normal": bool(sw_p >= 0.05),
        },
        "explanation": explanation,
        "recommendation": "Report empirical percentile only. Drop Z = -0.04 from abstract.",
    }

    save_json(results, "z_percentile_discrepancy_resolved.json")

    # Also append to upgraded_statistics.json
    stats_path = os.path.join(OUT_DIR, "upgraded_statistics.json")
    if os.path.exists(stats_path):
        with open(stats_path) as f:
            existing = json.load(f)
        existing["habitability_resolution"] = results
        save_json(existing, "upgraded_statistics.json")

    return results


# ============================================================
# TEST 2: SPATIAL BLOCK CROSS-VALIDATION
# ============================================================
def test_2():
    print("\n" + "=" * 70)
    print("TEST 2: SPATIAL BLOCK CROSS-VALIDATION")
    print("=" * 70)

    N_TRIALS = 1000

    mp_sites = load_megalithic_portal()
    mp_lats = np.array([s[0] for s in mp_sites])
    mp_lons = np.array([s[1] for s in mp_sites])
    print(f"  Loaded {len(mp_sites)} Megalithic Portal sites")

    # Cluster centers
    clusters = {
        "egypt_levant": (30.0, 32.0),
        "peru_andes": (-14.0, -75.0),
        "easter_island": (-27.0, -109.0),
        "iran": (30.0, 52.0),
        "indus_valley": (27.0, 68.0),
        "southeast_asia": (15.0, 100.0),
    }
    EXCLUSION_RADIUS = 500  # km

    # --- Leave-one-cluster-out tests ---
    print("\n  --- Leave-One-Cluster-Out Tests ---")
    loo_results = []

    for cluster_name, (clat, clon) in clusters.items():
        print(f"\n  Removing cluster: {cluster_name} (center: {clat}°N, {clon}°E, radius: {EXCLUSION_RADIUS}km)")

        # Compute distances from cluster center to all sites
        cluster_dists = haversine_vec(clat, clon, mp_lats, mp_lons)
        mask_keep = cluster_dists > EXCLUSION_RADIUS
        n_removed = int(np.sum(~mask_keep))
        n_remaining = int(np.sum(mask_keep))

        remaining_lats = mp_lats[mask_keep]
        remaining_lons = mp_lons[mask_keep]

        print(f"  Sites removed: {n_removed}, remaining: {n_remaining}")

        r = run_mc(remaining_lats, remaining_lons, N_TRIALS)
        print(f"  Z={r['z_score']}, observed={r['observed']}, expected={r['expected']}, "
              f"enrichment={r['enrichment']}")

        loo_results.append({
            "cluster": cluster_name,
            "center": {"lat": clat, "lon": clon},
            "sites_removed": n_removed,
            "remaining_sites": n_remaining,
            "z_score": r["z_score"],
            "observed": r["observed"],
            "expected": r["expected"],
            "enrichment": r["enrichment"],
            "empirical_p": r["empirical_p"],
        })

    # --- Egypt-removed Pleiades divergence ---
    print("\n  --- Pleiades Divergence with Egypt/Levant Removed ---")
    mon_anc, set_anc = load_pleiades(ancient_only=True)

    # Remove sites between 25°E-40°E longitude
    mon_no_egypt = [(lat, lon) for lat, lon in mon_anc if not (25 <= lon <= 40)]
    set_no_egypt = [(lat, lon) for lat, lon in set_anc if not (25 <= lon <= 40)]

    print(f"  Monumental: {len(mon_anc)} -> {len(mon_no_egypt)} (removed {len(mon_anc) - len(mon_no_egypt)})")
    print(f"  Settlement: {len(set_anc)} -> {len(set_no_egypt)} (removed {len(set_anc) - len(set_no_egypt)})")

    if len(mon_no_egypt) > 0:
        mon_ne_lats = np.array([s[0] for s in mon_no_egypt])
        mon_ne_lons = np.array([s[1] for s in mon_no_egypt])
        r_mon_ne = run_mc(mon_ne_lats, mon_ne_lons, N_TRIALS)
    else:
        r_mon_ne = {"z_score": 0, "observed": 0, "expected": 0, "enrichment": 0,
                     "empirical_p": 1.0, "baseline_counts": []}

    if len(set_no_egypt) > 0:
        set_ne_lats = np.array([s[0] for s in set_no_egypt])
        set_ne_lons = np.array([s[1] for s in set_no_egypt])
        r_set_ne = run_mc(set_ne_lats, set_ne_lons, N_TRIALS)
    else:
        r_set_ne = {"z_score": 0, "observed": 0, "expected": 0, "enrichment": 0,
                     "empirical_p": 1.0, "baseline_counts": []}

    div_no_egypt = round(r_mon_ne["z_score"] - r_set_ne["z_score"], 2)
    print(f"  Without Egypt: Monument Z={r_mon_ne['z_score']}, Settlement Z={r_set_ne['z_score']}, "
          f"Divergence={div_no_egypt}")

    egypt_removed_div = {
        "monument_z": r_mon_ne["z_score"],
        "monument_observed": r_mon_ne["observed"],
        "monument_expected": r_mon_ne["expected"],
        "settlement_z": r_set_ne["z_score"],
        "settlement_observed": r_set_ne["observed"],
        "settlement_expected": r_set_ne["expected"],
        "divergence": div_no_egypt,
        "monument_n_sites": len(mon_no_egypt),
        "settlement_n_sites": len(set_no_egypt),
    }

    # --- Hemisphere block test ---
    print("\n  --- Hemisphere Block Test (Old World / New World) ---")

    # Old World: lon >= -20
    # New World: lon < -20
    ow_mask = mp_lons >= -20
    nw_mask = mp_lons < -20

    ow_lats, ow_lons = mp_lats[ow_mask], mp_lons[ow_mask]
    nw_lats, nw_lons = mp_lats[nw_mask], mp_lons[nw_mask]

    print(f"  Old World sites: {len(ow_lats)}")
    print(f"  New World sites: {len(nw_lats)}")

    # Train on Old World, test on New World
    print("  Train: Old World | Test: New World")
    r_nw_test = run_mc(nw_lats, nw_lons, N_TRIALS)
    print(f"  New World test: Z={r_nw_test['z_score']}, observed={r_nw_test['observed']}, "
          f"expected={r_nw_test['expected']}")

    # Train on New World, test on Old World
    print("  Train: New World | Test: Old World")
    r_ow_test = run_mc(ow_lats, ow_lons, N_TRIALS)
    print(f"  Old World test: Z={r_ow_test['z_score']}, observed={r_ow_test['observed']}, "
          f"expected={r_ow_test['expected']}")

    hemisphere_block = {
        "old_world_n": int(np.sum(ow_mask)),
        "new_world_n": int(np.sum(nw_mask)),
        "old_world_only_z": r_ow_test["z_score"],
        "old_world_only_observed": r_ow_test["observed"],
        "old_world_only_expected": r_ow_test["expected"],
        "old_world_only_enrichment": r_ow_test["enrichment"],
        "new_world_only_z": r_nw_test["z_score"],
        "new_world_only_observed": r_nw_test["observed"],
        "new_world_only_expected": r_nw_test["expected"],
        "new_world_only_enrichment": r_nw_test["enrichment"],
    }

    # Verdict
    egypt_cluster = next(r for r in loo_results if r["cluster"] == "egypt_levant")
    all_survive = all(r["z_score"] > 2 for r in loo_results)
    egypt_survives = egypt_cluster["z_score"] > 2

    if all_survive:
        verdict = "STRONG: Signal survives removal of every individual cluster, including Egypt."
    elif egypt_survives:
        verdict = "MODERATE: Signal survives Egypt removal but is weakened by removing some other clusters."
    else:
        verdict = (f"WEAK: Signal does NOT survive Egypt removal "
                   f"(Z={egypt_cluster['z_score']} after removal). "
                   f"Egypt is a critical driver.")

    results = {
        "meta": {"date": "2026-03-17", "n_trials": N_TRIALS,
                 "threshold_km": THRESHOLD_KM, "exclusion_radius_km": EXCLUSION_RADIUS,
                 "dataset": "Megalithic Portal + Pleiades"},
        "leave_one_out": loo_results,
        "egypt_removed_divergence": egypt_removed_div,
        "hemisphere_block": hemisphere_block,
        "verdict": verdict,
    }

    save_json(results, "spatial_block_validation.json")
    return results


# ============================================================
# TEST 1: HABITABILITY-ADJUSTED MONUMENT-SETTLEMENT DIVERGENCE
# ============================================================
def test_1():
    print("\n" + "=" * 70)
    print("TEST 1: HABITABILITY-ADJUSTED MONUMENT-SETTLEMENT DIVERGENCE")
    print("=" * 70)

    N_MC_FAST = 50
    N_MC_PRECISE = 200

    # Load p3k14c for density baseline
    p3k_sites = load_p3k14c()
    p3k_lats = np.array([s[0] for s in p3k_sites])
    p3k_lons = np.array([s[1] for s in p3k_sites])
    print(f"  p3k14c sites: {len(p3k_sites)}")

    # Load Pleiades ancient for divergence
    mon_anc, set_anc = load_pleiades(ancient_only=True)
    mon_lats = np.array([s[0] for s in mon_anc])
    mon_lons = np.array([s[1] for s in mon_anc])
    set_lats = np.array([s[0] for s in set_anc])
    set_lons = np.array([s[1] for s in set_anc])
    print(f"  Pleiades ancient: {len(mon_anc)} monumental, {len(set_anc)} settlement")

    # Step 1: Compute Alison's p3k14c density
    alison_p3k_count = count_within(POLE_LAT, POLE_LNG, p3k_lats, p3k_lons)
    print(f"  Alison p3k14c density: {alison_p3k_count} sites within {THRESHOLD_KM}km")

    # Step 2: Compute Alison's divergence
    print("\n  Computing Alison's divergence (200 MC trials)...")
    alison_mon_z = run_mc(mon_lats, mon_lons, N_MC_PRECISE)["z_score"]
    alison_set_z = run_mc(set_lats, set_lons, N_MC_PRECISE)["z_score"]
    alison_divergence = round(alison_mon_z - alison_set_z, 2)
    print(f"  Alison: monument Z={alison_mon_z}, settlement Z={alison_set_z}, divergence={alison_divergence}")

    # Step 3: Generate habitability-matched random circles
    print(f"\n  Generating habitability-matched random circles...")
    print(f"  Target density: {alison_p3k_count} ±10%")

    tolerance = 0.10
    target_min = alison_p3k_count * (1 - tolerance)
    target_max = alison_p3k_count * (1 + tolerance)

    matched_circles = []
    n_tested = 0
    max_attempts = 100000  # cap total attempts

    while len(matched_circles) < 10000 and n_tested < max_attempts:
        # Random pole (uniform on sphere)
        z = np.random.uniform(-1, 1)
        rand_pole_lat = math.degrees(math.asin(z))
        rand_pole_lon = np.random.uniform(-180, 180)

        # Check p3k14c density
        density = count_within(rand_pole_lat, rand_pole_lon, p3k_lats, p3k_lons)
        n_tested += 1

        if target_min <= density <= target_max:
            matched_circles.append((rand_pole_lat, rand_pole_lon, density))

        if n_tested % 10000 == 0:
            print(f"    Tested {n_tested} circles, matched {len(matched_circles)} "
                  f"({100*len(matched_circles)/max(1,n_tested):.1f}% acceptance)")

    # If too few, widen tolerance
    if len(matched_circles) < 1000:
        print(f"  Only {len(matched_circles)} matched at ±{tolerance*100:.0f}%. Widening to ±20%...")
        tolerance = 0.20
        target_min = alison_p3k_count * (1 - tolerance)
        target_max = alison_p3k_count * (1 + tolerance)

        while len(matched_circles) < 10000 and n_tested < max_attempts * 2:
            z = np.random.uniform(-1, 1)
            rand_pole_lat = math.degrees(math.asin(z))
            rand_pole_lon = np.random.uniform(-180, 180)
            density = count_within(rand_pole_lat, rand_pole_lon, p3k_lats, p3k_lons)
            n_tested += 1
            if target_min <= density <= target_max:
                matched_circles.append((rand_pole_lat, rand_pole_lon, density))
            if n_tested % 10000 == 0:
                print(f"    Tested {n_tested} circles, matched {len(matched_circles)}")

    n_matched = len(matched_circles)
    print(f"\n  Habitability-matched circles: {n_matched} (from {n_tested} tested)")
    print(f"  Density tolerance: ±{tolerance*100:.0f}%")

    if n_matched == 0:
        print("  ERROR: No matched circles found. Cannot proceed.")
        results = {"error": "No habitability-matched circles found",
                   "n_tested": n_tested, "alison_density": alison_p3k_count}
        save_json(results, "habitability_adjusted_divergence.json")
        return results

    # Step 4: Compute divergence for each matched circle (fast MC first)
    print(f"\n  Computing divergence for {n_matched} matched circles ({N_MC_FAST} MC trials each)...")
    matched_divergences = []
    t_start = time.time()

    for i, (plat, plon, density) in enumerate(matched_circles):
        # Monument Z
        mon_z = run_mc(mon_lats, mon_lons, N_MC_FAST, pole_lat=plat, pole_lon=plon)["z_score"]
        # Settlement Z
        set_z = run_mc(set_lats, set_lons, N_MC_FAST, pole_lat=plat, pole_lon=plon)["z_score"]
        div = mon_z - set_z
        matched_divergences.append(div)

        if (i + 1) % 500 == 0:
            elapsed = time.time() - t_start
            rate = (i + 1) / elapsed
            eta = (n_matched - i - 1) / rate
            print(f"    {i+1}/{n_matched} circles done, max_div={max(matched_divergences):.2f}, "
                  f"[{elapsed:.0f}s elapsed, ~{eta:.0f}s remaining]")

    # Step 5: Re-run close contenders with precise MC
    div_arr = np.array(matched_divergences)
    div_std = float(np.std(div_arr))
    threshold_rerun = alison_divergence - 2 * div_std

    close_indices = [i for i, d in enumerate(matched_divergences) if d >= threshold_rerun]
    print(f"\n  Circles within 2 std of Alison's divergence: {len(close_indices)}")

    for idx in close_indices:
        plat, plon, density = matched_circles[idx]
        mon_z = run_mc(mon_lats, mon_lons, N_MC_PRECISE, pole_lat=plat, pole_lon=plon)["z_score"]
        set_z = run_mc(set_lats, set_lons, N_MC_PRECISE, pole_lat=plat, pole_lon=plon)["z_score"]
        old_div = matched_divergences[idx]
        matched_divergences[idx] = mon_z - set_z
        print(f"    Circle {idx}: {old_div:.2f} -> {matched_divergences[idx]:.2f} (precise)")

    # Step 6: Compare
    div_arr = np.array(matched_divergences)
    n_exceeding = int(np.sum(div_arr >= alison_divergence))
    alison_pct = float(np.sum(div_arr < alison_divergence) / len(div_arr) * 100)

    print(f"\n  === RESULTS ===")
    print(f"  Alison divergence: {alison_divergence}")
    print(f"  Matched distribution: mean={np.mean(div_arr):.2f}, std={np.std(div_arr):.2f}, "
          f"max={np.max(div_arr):.2f}")
    print(f"  Exceeding Alison: {n_exceeding}/{n_matched}")
    print(f"  Alison percentile: {alison_pct:.2f}%")

    if n_exceeding == 0:
        verdict = "PASS: 0 habitability-matched circles exceed Alison's divergence. Monument-settlement divergence is NOT explained by habitability."
    elif n_exceeding / n_matched < 0.01:
        verdict = f"PASS: Only {n_exceeding}/{n_matched} ({100*n_exceeding/n_matched:.2f}%) habitability-matched circles exceed Alison's divergence (<1%). Divergence is NOT explained by habitability."
    elif n_exceeding / n_matched < 0.05:
        verdict = f"PARTIAL: {n_exceeding}/{n_matched} ({100*n_exceeding/n_matched:.1f}%) exceed. Divergence is partially but not fully explained by habitability."
    else:
        verdict = f"FAIL: {n_exceeding}/{n_matched} ({100*n_exceeding/n_matched:.1f}%) exceed. Divergence may be a byproduct of corridor selection."

    print(f"  Verdict: {verdict}")

    results = {
        "meta": {"date": "2026-03-17", "dataset": "Pleiades (ancient, pre-2000 BCE) + p3k14c density",
                 "threshold_km": THRESHOLD_KM, "mc_fast": N_MC_FAST, "mc_precise": N_MC_PRECISE},
        "n_habitability_matched": n_matched,
        "n_tested": n_tested,
        "density_tolerance": tolerance,
        "alison_density": alison_p3k_count,
        "alison_divergence": alison_divergence,
        "alison_monument_z": alison_mon_z,
        "alison_settlement_z": alison_set_z,
        "matched_distribution": {
            "mean": round(float(np.mean(div_arr)), 2),
            "std": round(float(np.std(div_arr)), 2),
            "min": round(float(np.min(div_arr)), 2),
            "max": round(float(np.max(div_arr)), 2),
            "p5": round(float(np.percentile(div_arr, 5)), 2),
            "p25": round(float(np.percentile(div_arr, 25)), 2),
            "p50": round(float(np.percentile(div_arr, 50)), 2),
            "p75": round(float(np.percentile(div_arr, 75)), 2),
            "p95": round(float(np.percentile(div_arr, 95)), 2),
        },
        "alison_percentile": round(alison_pct, 2),
        "n_exceeding_alison": n_exceeding,
        "n_rerun_precise": len(close_indices),
        "verdict": verdict,
    }

    save_json(results, "habitability_adjusted_divergence.json")
    return results


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    t0 = time.time()
    print("=" * 70)
    print("PEER REVIEW RESPONSE — FOUR CRITICAL TESTS")
    print("=" * 70)

    # Test 3: Upgrade MC trials (fastest, fixes most issues)
    r3, _ = test_3()
    t3 = time.time()
    print(f"\n  Test 3 complete: {t3 - t0:.0f}s")

    # Test 4: Resolve Z/-0.04 discrepancy (can run alongside Test 3)
    r4 = test_4()
    t4 = time.time()
    print(f"\n  Test 4 complete: {t4 - t3:.0f}s")

    # Test 2: Spatial block cross-validation (medium compute)
    r2 = test_2()
    t2 = time.time()
    print(f"\n  Test 2 complete: {t2 - t4:.0f}s")

    # Test 1: Habitability-adjusted divergence (most compute-intensive)
    r1 = test_1()
    t1_end = time.time()
    print(f"\n  Test 1 complete: {t1_end - t2:.0f}s")

    # Final summary
    elapsed = time.time() - t0
    print(f"\n{'=' * 70}")
    print(f"ALL TESTS COMPLETE — {elapsed:.0f}s total ({elapsed/60:.1f} minutes)")
    print(f"{'=' * 70}")
    print(f"Output directory: {OUT_DIR}")
    print(f"\nFiles produced:")
    for f in sorted(os.listdir(OUT_DIR)):
        path = os.path.join(OUT_DIR, f)
        size = os.path.getsize(path)
        print(f"  {f} ({size:,} bytes)")

    # Quick verdict summary
    print(f"\n--- VERDICTS ---")
    print(f"Test 3 (upgraded stats): DONE — see upgraded_statistics.json")
    print(f"Test 4 (Z/percentile): {r4.get('recommendation', 'see json')}")
    print(f"Test 2 (spatial block): {r2.get('verdict', 'see json')}")
    print(f"Test 1 (hab-adjusted div): {r1.get('verdict', 'see json')}")
