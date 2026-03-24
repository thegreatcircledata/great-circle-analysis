#!/usr/bin/env python3
"""
Oral Tradition & Mythology Spatial Mapping — Great Circle Analysis
==================================================================
Tests whether specific mythological motifs cluster along the Great Circle
more than expected by chance, using Berezkin's Analytical Catalogue.

Data source: macleginn/mythology-queries (GitHub), which provides a parsed
version of Berezkin & Duvakin's Electronic Analytical Catalogue.

Method: For each motif category, compute enrichment ratio (on-corridor
prevalence / off-corridor prevalence) and Monte Carlo significance using
10,000 random great circles.
"""

import csv
import json
import math
import os
import random
import sys
from collections import defaultdict

# ============================================================
# CONSTANTS (from great-circle-analysis/utils.py)
# ============================================================
POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2  # ~10,007.5 km
CORRIDOR_WIDTH_KM = 200  # wider band for ethnic territories
N_MONTE_CARLO = 10000
SEED = 42

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ============================================================
# TARGET MOTIF CATEGORIES
# ============================================================
# Each category maps to a list of Berezkin motif codes present in the dataset.
# Selected based on the directive's criteria: deep antiquity, wide distribution,
# and relevance to corridor/travel hypotheses.

TARGET_MOTIFS = {
    "flood_deluge": {
        "description": "Flood and deluge myths — universal cataclysmic water narratives",
        "codes": [
            "c2_3",    # Deluge and conflagration combined
            "c4_3",    # The flood: fruits fall from a tree
            "c7_3",    # The flood: breaking the dam
            "c10_7",   # The flood: the wet tails
            "c10a_7",  # The flood: birds cling to the sky
            "c9a_7",   # The flood: people turn into toads
            "c5a_3",   # Bird-scouts (post-flood)
            "c5b_3",   # Animal-scouts (post-flood)
            "c8_5",    # The primeval couple of siblings (repopulation)
            "c8a_5",   # Grinding stones match
            "c8b_5",   # Siblings change their looks and marry
            "L28a_10", # Eating of snake meat triggers flood
            "g9a_3",   # Restored forest and flood
        ],
    },
    "earth_diver": {
        "description": "Earth-diver creation myths — being dives into primordial waters to bring up earth",
        "codes": [
            "c6_10",   # The diver
            "c6a_3",   # The diver is turtle or frog
            "c6b_3",   # The diver is muskrat or beaver
            "c6c_3",   # The diver is a bird
            "c6c1_3",  # Birds: successful and unsuccessful divers
            "c6d_3",   # The earth-diver
            "c6e_3",   # The diver is a crustacean
            "c6g_3",   # The diver is a wild boar
            "c6h_3",   # An insect brings the earth
            "c6I_3",   # Dirt stuck to body turns into the earth
            "b4_3",    # The fished out earth
        ],
    },
    "world_axis_sky_support": {
        "description": "World axis, sky support, world tree — axis mundi concepts",
        "codes": [
            "i12_3",    # The world axis
            "b21_3",    # Destruction of the world tree
            "b77_3",    # Primeval sky close to earth
            "b77a_3",   # Giant pushed the sky up
            "b77b_3",   # Sky touched with a long object
            "b77b1_3",  # Sky touched with a pestle
            "b77b2_3",  # Sky touched with a broom
            "b77c_3",   # Serpent pushes sky up
            "c24_3",    # The world in danger of falling down
            "b1a_3",    # Three worlds: the universe divided between gods
            "b2d_3",    # Marriage of the sky and the earth
        ],
    },
    "travelling_transformer": {
        "description": "Travelling transformer / culture hero who journeys and reshapes the landscape",
        "codes": [
            "b28_3",    # Travelling transformer
            "b28a_8",   # (sub-variant)
            "b28b_7",   # (sub-variant)
            "b28c_7",   # (sub-variant)
            "b28d_7",   # (sub-variant)
            "b80_3",    # Measuring of the world
        ],
    },
    "cosmic_hunt": {
        "description": "Cosmic hunt — d'Huy's best-dated myth, possibly Paleolithic",
        "codes": [
            "b42_2",    # Cosmic hunt
            "b42a_2",   # Blood or fat drops to earth
            "b42b_2",   # Sky hunters pursue an ungulate
            "b42c_2",   # Sky hunters pursue a bear
            "b42d_2",   # Sky hunters pursue a tapir
            "b42e_2",   # Sky hunters pursue a rhea
            "b42f_2",   # (variant)
            "b42h_2",   # Object of cosmic hunt is Belt of Orion
            "b42h1_2",  # (variant)
            "b42i_2",   # Object is Lady in the Chair
            "b42k_2",   # Cosmic hunt and the Pleiades
            "b42l_2",   # (variant)
            "b42m_2",   # (variant)
            "B42m1_2",  # (variant)
            "b42mm_2",  # (variant)
            "b42mn_2",  # (variant)
            "b42n_2",   # Orion is one person
            "b42o_2",   # (variant)
            "B42o1_2",  # Ursa major is a fisher
            "b42p_2",   # Ursa major is a bear
            "b42q_2",   # Ursa major is a carriage
            "b42r_2",   # Belt of Orion: one runs after another
            "b42s_2",   # (variant)
            "b42t_2",   # Ursa major is a big mammal
        ],
    },
    "star_celestial": {
        "description": "Star myths and celestial orientation — Pleiades, Milky Way, Polaris, constellations",
        "codes": [
            "i55_2",    # Stars are openings
            "i72_2",    # Stars are people
            "i73_2",    # Stars are sparks
            "i85_2",    # Polaris is a pole, a nail
            "i85a_2",   # Horses around Polaris
            "i85b_2",   # Polaris is a man
            "i100_2",   # The Pleiades are girls
            "i100a_2",  # The Pleiades are mother with children
            "i100b_2",  # The Pleiades are a group of people
            "i102_2",   # Milky Way is a tree
            "i109_2",   # Milky Way is the path of the Sun
            "i58_2",    # Milky Way is the way of birds
            "i62_2",    # Milky Way is a river
            "i65_2",    # Milky Way of the dead
            "b46_2",    # Ursa major as seven men
            "b46a_2",   # Stolen star of the Pleiades
            "b47_2",    # The Pleiades bring cold
            "i99_2",    # The Pleiades are boys or men
            "i101_2",   # Big Dipper is poles, a nailed skin
            "i128_2",   # Ursa major is a dipper
        ],
    },
    "giants_earlier_race": {
        "description": "Giants, earlier races, and primordial beings",
        "codes": [
            "i20a_8",   # The sky giants
            "k54_10",   # Two giants
            "b68_7",    # The giant grouse
            "e1a_5",    # First people of unstable materials
            "i16_5",    # Body anomalies of the first people
            "c13_3",    # The objects' revolt (animated objects)
            "c14_3",    # Monsters destroy people
        ],
    },
}


# ============================================================
# GEODESIC MATH
# ============================================================

def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))


def gc_distance(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    return abs(haversine_km(lat, lon, pole_lat, pole_lon) - QUARTER_CIRC)


def random_pole():
    z = random.uniform(-1, 1)
    lat = math.degrees(math.asin(z))
    lon = random.uniform(-180, 180)
    return lat, lon


def latitude_matched_random_pole(target_lats, n_bins=18):
    """Generate a random pole whose great circle has a similar latitude profile
    to the real Great Circle. This controls for the fact that motif distributions
    are latitude-dependent."""
    # Accept-reject: generate random pole, check that the resulting great circle
    # passes through a similar range of latitudes
    target_min = min(target_lats)
    target_max = max(target_lats)
    target_range = target_max - target_min
    for _ in range(1000):
        plat, plon = random_pole()
        # Quick check: the great circle defined by this pole spans latitudes
        # from (90 - |pole_lat|) in both hemispheres
        gc_max_lat = 90.0 - abs(plat)
        gc_min_lat = -gc_max_lat
        # Check that the latitude range of this circle overlaps substantially
        # with the latitude range of ethnic groups on the real corridor
        overlap_min = max(gc_min_lat, target_min)
        overlap_max = min(gc_max_lat, target_max)
        if overlap_max > overlap_min:
            overlap = overlap_max - overlap_min
            if overlap >= 0.5 * target_range:
                return plat, plon
    # Fallback: just return random
    return random_pole()


# ============================================================
# DATA LOADING
# ============================================================

def load_data():
    """Load Berezkin CSV and coords, return list of ethnic group dicts."""
    # Load the main CSV
    csv_path = os.path.join(DATA_DIR, "berezkin_new.csv")
    with open(csv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames
        rows = list(reader)

    # Identify motif columns (everything after the metadata columns)
    meta_cols = {"idn", "groups", "latit", "longit", "sum", "cosm_mot",
                 "areas1", "AREA2", "areas3", "lang1", "lang2", "lang3"}
    motif_cols = [h for h in headers if h not in meta_cols]

    groups = []
    for row in rows:
        try:
            lat = float(row["latit"])
            lon = float(row["longit"])
        except (ValueError, KeyError):
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue

        motifs = {}
        for mc in motif_cols:
            try:
                motifs[mc] = int(row[mc])
            except (ValueError, KeyError):
                motifs[mc] = 0

        groups.append({
            "name": row.get("groups", ""),
            "lat": lat,
            "lon": lon,
            "area": row.get("areas1", ""),
            "lang1": row.get("lang1", ""),
            "motifs": motifs,
            "gc_dist": gc_distance(lat, lon),
        })

    return groups, motif_cols


def validate_motif_codes(groups, motif_cols):
    """Check which target motif codes actually exist in the dataset."""
    available = set(motif_cols)
    report = {}
    for cat_name, cat_info in TARGET_MOTIFS.items():
        found = [c for c in cat_info["codes"] if c in available]
        missing = [c for c in cat_info["codes"] if c not in available]
        report[cat_name] = {"found": found, "missing": missing}
        if missing:
            print(f"  WARNING: {cat_name} — {len(missing)} codes not in dataset: {missing}")
    return report


# ============================================================
# ANALYSIS: PER-MOTIF ENRICHMENT TEST
# ============================================================

def compute_enrichment(groups, motif_codes, corridor_km=CORRIDOR_WIDTH_KM,
                       pole_lat=POLE_LAT, pole_lon=POLE_LON):
    """Compute enrichment ratio for a set of motif codes.

    A group "has" the category if ANY of the motif codes = 1.
    """
    on_corridor_has = 0
    on_corridor_total = 0
    off_corridor_has = 0
    off_corridor_total = 0

    for g in groups:
        dist = abs(haversine_km(g["lat"], g["lon"], pole_lat, pole_lon) - QUARTER_CIRC)
        has_motif = any(g["motifs"].get(code, 0) == 1 for code in motif_codes)

        if dist <= corridor_km:
            on_corridor_total += 1
            if has_motif:
                on_corridor_has += 1
        else:
            off_corridor_total += 1
            if has_motif:
                off_corridor_has += 1

    on_prev = on_corridor_has / on_corridor_total if on_corridor_total > 0 else 0
    off_prev = off_corridor_has / off_corridor_total if off_corridor_total > 0 else 0
    enrichment = on_prev / off_prev if off_prev > 0 else float("inf")

    return {
        "on_corridor_has": on_corridor_has,
        "on_corridor_total": on_corridor_total,
        "on_corridor_prevalence": round(on_prev, 4),
        "off_corridor_has": off_corridor_has,
        "off_corridor_total": off_corridor_total,
        "off_corridor_prevalence": round(off_prev, 4),
        "enrichment_ratio": round(enrichment, 4),
    }


def monte_carlo_enrichment(groups, motif_codes, n_trials=N_MONTE_CARLO,
                           corridor_km=CORRIDOR_WIDTH_KM, seed=SEED,
                           latitude_matched=True):
    """Monte Carlo test: generate random great circles and compute enrichment
    for each, then rank the real Great Circle."""
    random.seed(seed)

    # Real enrichment
    real = compute_enrichment(groups, motif_codes, corridor_km)

    # Get latitude profile of on-corridor groups for matching
    on_lats = [g["lat"] for g in groups if g["gc_dist"] <= corridor_km]

    # Monte Carlo
    mc_enrichments = []
    for i in range(n_trials):
        if latitude_matched and on_lats:
            plat, plon = latitude_matched_random_pole(on_lats)
        else:
            plat, plon = random_pole()
        mc_result = compute_enrichment(groups, motif_codes, corridor_km, plat, plon)
        mc_enrichments.append(mc_result["enrichment_ratio"])

    # Statistics
    mc_enrichments.sort()
    percentile = sum(1 for e in mc_enrichments if e <= real["enrichment_ratio"]) / n_trials
    p_value = 1.0 - percentile
    mean_mc = sum(mc_enrichments) / len(mc_enrichments)
    std_mc = (sum((e - mean_mc) ** 2 for e in mc_enrichments) / len(mc_enrichments)) ** 0.5
    z_score = (real["enrichment_ratio"] - mean_mc) / std_mc if std_mc > 0 else 0

    real["mc_mean_enrichment"] = round(mean_mc, 4)
    real["mc_std_enrichment"] = round(std_mc, 4)
    real["mc_percentile"] = round(percentile * 100, 1)
    real["mc_p_value"] = round(p_value, 4)
    real["mc_z_score"] = round(z_score, 2)
    real["n_trials"] = n_trials
    real["latitude_matched"] = latitude_matched

    return real


# ============================================================
# ANALYSIS: MULTI-MOTIF COMPOSITE TEST
# ============================================================

def composite_motif_score(groups, corridor_km=CORRIDOR_WIDTH_KM):
    """For each group, count how many target categories are present.
    Compare mean score on-corridor vs off-corridor."""
    scores = []
    for g in groups:
        score = 0
        for cat_name, cat_info in TARGET_MOTIFS.items():
            codes = cat_info["codes"]
            if any(g["motifs"].get(c, 0) == 1 for c in codes):
                score += 1
        scores.append(score)

    on_scores = []
    off_scores = []
    for g, s in zip(groups, scores):
        if g["gc_dist"] <= corridor_km:
            on_scores.append(s)
        else:
            off_scores.append(s)

    on_mean = sum(on_scores) / len(on_scores) if on_scores else 0
    off_mean = sum(off_scores) / len(off_scores) if off_scores else 0

    return {
        "n_categories": len(TARGET_MOTIFS),
        "on_corridor_n": len(on_scores),
        "on_corridor_mean_score": round(on_mean, 3),
        "off_corridor_n": len(off_scores),
        "off_corridor_mean_score": round(off_mean, 3),
        "difference": round(on_mean - off_mean, 3),
        "scores_on": on_scores,
        "scores_off": off_scores,
    }


def monte_carlo_composite(groups, n_trials=N_MONTE_CARLO, corridor_km=CORRIDOR_WIDTH_KM,
                          seed=SEED):
    """Monte Carlo for the composite motif score test."""
    random.seed(seed)

    # Real composite
    real = composite_motif_score(groups, corridor_km)
    real_diff = real["difference"]

    # Precompute per-group scores
    group_scores = []
    for g in groups:
        score = 0
        for cat_info in TARGET_MOTIFS.values():
            if any(g["motifs"].get(c, 0) == 1 for c in cat_info["codes"]):
                score += 1
        group_scores.append(score)

    on_lats = [g["lat"] for g in groups if g["gc_dist"] <= corridor_km]

    mc_diffs = []
    for _ in range(n_trials):
        plat, plon = latitude_matched_random_pole(on_lats) if on_lats else random_pole()
        on_s = []
        off_s = []
        for g, s in zip(groups, group_scores):
            dist = abs(haversine_km(g["lat"], g["lon"], plat, plon) - QUARTER_CIRC)
            if dist <= corridor_km:
                on_s.append(s)
            else:
                off_s.append(s)
        on_m = sum(on_s) / len(on_s) if on_s else 0
        off_m = sum(off_s) / len(off_s) if off_s else 0
        mc_diffs.append(on_m - off_m)

    mc_diffs.sort()
    percentile = sum(1 for d in mc_diffs if d <= real_diff) / n_trials
    p_value = 1.0 - percentile
    mean_mc = sum(mc_diffs) / len(mc_diffs)
    std_mc = (sum((d - mean_mc) ** 2 for d in mc_diffs) / len(mc_diffs)) ** 0.5
    z_score = (real_diff - mean_mc) / std_mc if std_mc > 0 else 0

    # Remove raw score lists for JSON output
    del real["scores_on"]
    del real["scores_off"]

    real["mc_mean_diff"] = round(mean_mc, 4)
    real["mc_std_diff"] = round(std_mc, 4)
    real["mc_percentile"] = round(percentile * 100, 1)
    real["mc_p_value"] = round(p_value, 4)
    real["mc_z_score"] = round(z_score, 2)
    real["n_trials"] = n_trials

    return real


# ============================================================
# OUTPUT: MASTER TABLE
# ============================================================

def write_motif_distributions(groups, motif_cols):
    """Write the master CSV: motif_code | ethnic_group | lat | lon | on_corridor | present"""
    path = os.path.join(OUTPUT_DIR, "motif_distributions.csv")
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ethnic_group", "lat", "lon", "gc_distance_km", "on_corridor",
                          "area", "lang_family", "category", "motif_code", "motif_present"])
        for g in groups:
            on_corr = 1 if g["gc_dist"] <= CORRIDOR_WIDTH_KM else 0
            for cat_name, cat_info in TARGET_MOTIFS.items():
                for code in cat_info["codes"]:
                    if code in g["motifs"]:
                        writer.writerow([
                            g["name"], g["lat"], g["lon"],
                            round(g["gc_dist"], 1), on_corr,
                            g["area"], g["lang1"],
                            cat_name, code, g["motifs"].get(code, 0),
                        ])
    print(f"  Written: {path}")


def save_json(data, filename):
    path = os.path.join(OUTPUT_DIR, filename)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    print(f"  Written: {path}")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("ORAL TRADITION & MYTHOLOGY SPATIAL MAPPING")
    print("Great Circle Corridor Analysis")
    print("=" * 70)
    print()

    # Load data
    print("[1/6] Loading Berezkin data...")
    groups, motif_cols = load_data()
    print(f"  Loaded {len(groups)} ethnic groups with {len(motif_cols)} motif columns")

    # Validate motif codes
    print("\n[2/6] Validating target motif codes...")
    code_report = validate_motif_codes(groups, motif_cols)
    for cat_name, info in code_report.items():
        print(f"  {cat_name}: {len(info['found'])}/{len(info['found']) + len(info['missing'])} codes found")

    # Corridor stats
    on_corridor = [g for g in groups if g["gc_dist"] <= CORRIDOR_WIDTH_KM]
    print(f"\n  Groups within {CORRIDOR_WIDTH_KM}km of Great Circle: {len(on_corridor)}/{len(groups)}")
    print(f"  Corridor fraction: {len(on_corridor)/len(groups)*100:.1f}%")

    # Write master table
    print("\n[3/6] Writing motif distributions CSV...")
    write_motif_distributions(groups, motif_cols)

    # Per-motif enrichment with Monte Carlo
    print(f"\n[4/6] Running per-motif enrichment tests ({N_MONTE_CARLO} MC trials each)...")
    enrichment_results = {}
    for cat_name, cat_info in TARGET_MOTIFS.items():
        valid_codes = [c for c in cat_info["codes"] if c in set(motif_cols)]
        if not valid_codes:
            print(f"  SKIP {cat_name}: no valid codes")
            continue
        print(f"  Testing {cat_name} ({len(valid_codes)} motifs)...", end=" ", flush=True)
        result = monte_carlo_enrichment(groups, valid_codes)
        result["description"] = cat_info["description"]
        result["n_motifs_tested"] = len(valid_codes)
        result["motif_codes"] = valid_codes
        enrichment_results[cat_name] = result
        sig = "***" if result["mc_p_value"] < 0.001 else "**" if result["mc_p_value"] < 0.01 else "*" if result["mc_p_value"] < 0.05 else "ns"
        print(f"enrichment={result['enrichment_ratio']:.2f}, p={result['mc_p_value']:.4f} {sig}")

    save_json(enrichment_results, "enrichment_by_motif.json")

    # Composite test
    print(f"\n[5/6] Running composite motif score test ({N_MONTE_CARLO} MC trials)...")
    composite_result = monte_carlo_composite(groups)
    save_json(composite_result, "composite_motif_score.json")
    print(f"  On-corridor mean: {composite_result['on_corridor_mean_score']:.3f}")
    print(f"  Off-corridor mean: {composite_result['off_corridor_mean_score']:.3f}")
    print(f"  Difference: {composite_result['difference']:.3f}, p={composite_result['mc_p_value']:.4f}")

    # Summary
    print("\n[6/6] Generating summary...")
    summary = {
        "analysis": "Oral Tradition & Mythology Spatial Mapping",
        "data_source": "Berezkin Analytical Catalogue (via macleginn/mythology-queries)",
        "n_ethnic_groups": len(groups),
        "n_on_corridor": len(on_corridor),
        "corridor_width_km": CORRIDOR_WIDTH_KM,
        "n_monte_carlo_trials": N_MONTE_CARLO,
        "pole": {"lat": POLE_LAT, "lon": POLE_LON},
        "categories_tested": len(enrichment_results),
        "significant_at_005": sum(1 for r in enrichment_results.values() if r["mc_p_value"] < 0.05),
        "significant_at_001": sum(1 for r in enrichment_results.values() if r["mc_p_value"] < 0.01),
        "per_motif_results": {
            k: {
                "enrichment": v["enrichment_ratio"],
                "p_value": v["mc_p_value"],
                "z_score": v["mc_z_score"],
                "on_corridor_prevalence": v["on_corridor_prevalence"],
                "off_corridor_prevalence": v["off_corridor_prevalence"],
            }
            for k, v in enrichment_results.items()
        },
        "composite_test": {
            "on_mean": composite_result["on_corridor_mean_score"],
            "off_mean": composite_result["off_corridor_mean_score"],
            "p_value": composite_result["mc_p_value"],
            "z_score": composite_result["mc_z_score"],
        },
    }
    save_json(summary, "analysis_summary.json")

    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"\nEthnic groups: {len(groups)} total, {len(on_corridor)} on-corridor ({len(on_corridor)/len(groups)*100:.1f}%)")
    print(f"\nPer-category enrichment (latitude-matched MC, {N_MONTE_CARLO} trials):")
    print(f"{'Category':<30} {'Enrich':>8} {'On-prev':>8} {'Off-prev':>9} {'p-value':>8} {'Z':>6}")
    print("-" * 75)
    for cat, r in sorted(enrichment_results.items(), key=lambda x: x[1]["mc_p_value"]):
        sig = "***" if r["mc_p_value"] < 0.001 else "**" if r["mc_p_value"] < 0.01 else "*" if r["mc_p_value"] < 0.05 else ""
        print(f"{cat:<30} {r['enrichment_ratio']:>8.3f} {r['on_corridor_prevalence']:>8.3f} "
              f"{r['off_corridor_prevalence']:>9.3f} {r['mc_p_value']:>8.4f} {r['mc_z_score']:>5.1f} {sig}")

    print(f"\nComposite test: on={composite_result['on_corridor_mean_score']:.3f}, "
          f"off={composite_result['off_corridor_mean_score']:.3f}, "
          f"p={composite_result['mc_p_value']:.4f}, z={composite_result['mc_z_score']:.1f}")

    print("\nDone. Outputs in:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
