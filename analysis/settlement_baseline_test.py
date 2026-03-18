#!/usr/bin/env python3
"""
Settlement Baseline Test
========================
Determines whether MONUMENTAL sites cluster on the Great Circle more than
ordinary settlements in the same regions. Uses the Pleiades Gazetteer.

Distribution-matched Monte Carlo: 200 trials, ±2° jitter, 50 km threshold.
"""

import csv, math, random, json, os, sys
from collections import Counter

# ============================================================
# CONFIGURATION
# ============================================================
POLE_LAT = 59.682122
POLE_LNG = -138.646087
EARTH_R_KM = 6371
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
THRESHOLD_KM = 50
N_TRIALS = 200
JITTER_DEG = 2.0

DATA_FILE = os.path.join(os.path.dirname(__file__), "pleiades-places-latest.csv")
OUT_FILE = os.path.join(os.path.dirname(__file__), "results", "settlement_baseline_test.json")

# Feature type classifications
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


# ============================================================
# HELPERS
# ============================================================
def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))


def dist_from_great_circle(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LNG)
    return abs(d - QUARTER_CIRC)


def load_pleiades():
    """Load and classify Pleiades sites into monumental vs settlement."""
    monumental = []
    settlements = []
    monumental_ancient = []
    settlements_ancient = []

    with open(DATA_FILE, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["reprLat"])
                lon = float(row["reprLong"])
            except (ValueError, KeyError, TypeError):
                continue

            if lat == 0 and lon == 0:
                continue
            if not (-90 <= lat <= 90 and -180 <= lon <= 180):
                continue

            feature_types = {t.strip() for t in row.get("featureTypes", "").split(",") if t.strip()}

            is_monumental = bool(feature_types & MONUMENTAL_TYPES)
            is_settlement = bool(feature_types & SETTLEMENT_TYPES)

            # If a site is both, classify as monumental (temples in settlements)
            if is_monumental:
                monumental.append((lat, lon))
                try:
                    min_date = float(row.get("minDate", ""))
                    if min_date < -2000:
                        monumental_ancient.append((lat, lon))
                except (ValueError, TypeError):
                    pass
            elif is_settlement:
                settlements.append((lat, lon))
                try:
                    min_date = float(row.get("minDate", ""))
                    if min_date < -2000:
                        settlements_ancient.append((lat, lon))
                except (ValueError, TypeError):
                    pass

    return monumental, settlements, monumental_ancient, settlements_ancient


def count_within(sites, threshold_km):
    return sum(1 for lat, lon in sites if dist_from_great_circle(lat, lon) <= threshold_km)


def monte_carlo(sites, n_trials, threshold_km, jitter_deg):
    """Distribution-matched Monte Carlo: pick random real site coords + jitter."""
    lats = [s[0] for s in sites]
    lons = [s[1] for s in sites]
    n = len(sites)
    counts = []

    for trial in range(n_trials):
        c = 0
        for _ in range(n):
            lat = random.choice(lats) + random.gauss(0, jitter_deg)
            lon = random.choice(lons) + random.gauss(0, jitter_deg)
            lat = max(-90, min(90, lat))
            lon = max(-180, min(180, lon))
            if dist_from_great_circle(lat, lon) <= threshold_km:
                c += 1
        counts.append(c)
        if (trial + 1) % 50 == 0:
            print(f"    {trial+1}/{n_trials} done")

    return counts


def analyze(label, sites, threshold_km, n_trials, jitter_deg):
    """Run full analysis on a group of sites."""
    n = len(sites)
    if n == 0:
        return {"n_sites": 0, "within_50km": 0, "expected": 0, "z_score": 0, "enrichment": 0}

    observed = count_within(sites, threshold_km)
    print(f"\n  {label}: {n} sites, {observed} within {threshold_km} km")

    print(f"  Running {n_trials}-trial Monte Carlo...")
    mc_counts = monte_carlo(sites, n_trials, threshold_km, jitter_deg)

    mu = sum(mc_counts) / len(mc_counts)
    sigma = (sum((x - mu)**2 for x in mc_counts) / len(mc_counts)) ** 0.5
    enrichment = observed / mu if mu > 0 else float('inf')
    z = (observed - mu) / sigma if sigma > 0 else 0

    print(f"  Expected: {mu:.1f} ± {sigma:.1f}")
    print(f"  Enrichment: {enrichment:.2f}x | Z-score: {z:.2f}")

    return {
        "n_sites": n,
        "within_50km": observed,
        "expected": round(mu, 1),
        "z_score": round(z, 2),
        "enrichment": round(enrichment, 2),
    }


# ============================================================
# MAIN
# ============================================================
def main():
    random.seed(42)

    print("=" * 60)
    print("SETTLEMENT BASELINE TEST")
    print("Pleiades Gazetteer — Monumental vs Ordinary Settlements")
    print("=" * 60)

    print("\nLoading Pleiades data...")
    monumental, settlements, monumental_ancient, settlements_ancient = load_pleiades()
    print(f"  Monumental sites: {len(monumental)}")
    print(f"  Settlement sites: {len(settlements)}")
    print(f"  Monumental (pre-2000 BCE): {len(monumental_ancient)}")
    print(f"  Settlements (pre-2000 BCE): {len(settlements_ancient)}")

    # Run all four analyses
    print("\n" + "=" * 60)
    print("ALL PERIODS")
    print("=" * 60)
    r_monumental = analyze("MONUMENTAL", monumental, THRESHOLD_KM, N_TRIALS, JITTER_DEG)
    r_settlements = analyze("SETTLEMENTS", settlements, THRESHOLD_KM, N_TRIALS, JITTER_DEG)

    print("\n" + "=" * 60)
    print("ANCIENT ONLY (pre-2000 BCE)")
    print("=" * 60)
    r_monumental_ancient = analyze("MONUMENTAL (ancient)", monumental_ancient, THRESHOLD_KM, N_TRIALS, JITTER_DEG)
    r_settlements_ancient = analyze("SETTLEMENTS (ancient)", settlements_ancient, THRESHOLD_KM, N_TRIALS, JITTER_DEG)

    # Determine verdict
    z_m = r_monumental["z_score"]
    z_s = r_settlements["z_score"]
    z_diff = z_m - z_s

    if z_diff > 2:
        verdict = "MONUMENTAL ENRICHED"
    elif z_diff < -2:
        verdict = "SETTLEMENTS ENRICHED"
    else:
        verdict = "EQUAL"

    # Ancient verdict
    z_ma = r_monumental_ancient["z_score"]
    z_sa = r_settlements_ancient["z_score"]
    z_diff_ancient = z_ma - z_sa

    if z_diff > 2 and z_diff_ancient > 2:
        verdict_detail = "STRONGEST: Monumental enrichment holds across all periods"
    elif z_diff > 2:
        verdict_detail = "Monumental enrichment in all periods, weaker in ancient subset"
    elif abs(z_diff) <= 2:
        verdict_detail = "Geographic coincidence — both groups cluster equally"
    else:
        verdict_detail = "Settlements cluster more — unexpected result"

    results = {
        "monumental": r_monumental,
        "settlements": r_settlements,
        "monumental_ancient": r_monumental_ancient,
        "settlements_ancient": r_settlements_ancient,
        "z_difference_all": round(z_diff, 2),
        "z_difference_ancient": round(z_diff_ancient, 2),
        "verdict": verdict,
        "verdict_detail": verdict_detail,
        "methodology": {
            "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LNG},
            "threshold_km": THRESHOLD_KM,
            "monte_carlo_trials": N_TRIALS,
            "jitter_degrees": JITTER_DEG,
            "data_source": "Pleiades Gazetteer (pleiades-places-latest.csv)",
            "monumental_types": sorted(MONUMENTAL_TYPES),
            "settlement_types": sorted(SETTLEMENT_TYPES),
        }
    }

    # Print summary
    print("\n" + "=" * 60)
    print("FINAL RESULTS")
    print("=" * 60)
    print(f"\n{'Group':<25} {'N':>6} {'Within':>7} {'Expected':>9} {'Enrich':>7} {'Z':>7}")
    print("-" * 62)
    for label, r in [("Monumental", r_monumental), ("Settlements", r_settlements),
                      ("Monumental (ancient)", r_monumental_ancient),
                      ("Settlements (ancient)", r_settlements_ancient)]:
        print(f"{label:<25} {r['n_sites']:>6} {r['within_50km']:>7} {r['expected']:>9} "
              f"{r.get('enrichment',''):>6}x {r['z_score']:>7}")

    print(f"\nZ-score difference (all periods): {z_diff:.2f}")
    print(f"Z-score difference (ancient):     {z_diff_ancient:.2f}")
    print(f"\nVERDICT: {verdict}")
    print(f"DETAIL:  {verdict_detail}")

    os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
    with open(OUT_FILE, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {OUT_FILE}")


if __name__ == "__main__":
    main()
