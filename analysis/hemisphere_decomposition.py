#!/usr/bin/env python3
"""
Hemisphere Decomposition Directive
====================================
Determines whether the Great Circle signal is a global phenomenon or primarily
an Old World corridor effect. Seven studies decomposing the signal by hemisphere,
continent, site type, and habitability.

Studies:
  1. Pleiades Settlement Baseline by Hemisphere (geographic segments)
  2. Megalithic Portal Old World vs New World Split
  3. p3k14c Continental Signal Decomposition
  4. Old World Corridor Enrichment Profile
  5. New World Type Enrichment
  6. Ice Age (LGM) Connectivity Test
  7. "Boring Explanation" Quantification (habitability-adjusted signal)
"""

import csv, math, random, json, os, time, glob, re
import xml.etree.ElementTree as ET
import numpy as np
from collections import defaultdict, Counter
from scipy import stats as sp_stats

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

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "hemisphere_decomposition")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# VECTORIZED HELPERS (from nazca_followup.py pattern)
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

def scalar_gc_dist(lat, lon):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat, lon, POLE_LAT, POLE_LNG])
    dlat = lat2 - lat1; dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    d = EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))
    return abs(d - QUARTER_CIRC)

def rand_matched_batch(site_lats, site_lons, n):
    idx = np.random.randint(0, len(site_lats), size=n)
    lats = site_lats[idx] + np.random.normal(0, JITTER_DEG, n)
    lons = site_lons[idx] + np.random.normal(0, JITTER_DEG, n)
    return np.clip(lats, -90, 90), np.clip(lons, -180, 180)

def run_mc(site_lats, site_lons, threshold=THRESHOLD_KM, n_trials=N_TRIALS,
           pole_lat=POLE_LAT, pole_lon=POLE_LNG):
    """Run distribution-matched MC and return results dict."""
    n = len(site_lats)
    if n == 0:
        return {"n_sites": 0, "observed": 0, "expected": 0, "std": 0, "z_score": 0,
                "enrichment": 0, "p_value": 1.0}
    dists = gc_dist_vec(pole_lat, pole_lon, site_lats, site_lons)
    observed = int(np.sum(dists <= threshold))
    rand_counts = []
    for _ in range(n_trials):
        rl, rn = rand_matched_batch(site_lats, site_lons, n)
        rd = gc_dist_vec(pole_lat, pole_lon, rl, rn)
        rand_counts.append(int(np.sum(rd <= threshold)))
    mu = np.mean(rand_counts); sigma = np.std(rand_counts)
    z = float((observed - mu) / sigma) if sigma > 0 else 0.0
    p_val = float(1 - sp_stats.norm.cdf(z)) if z > 0 else 1.0
    enrich = float(observed / mu) if mu > 0 else 0.0
    return {
        "n_sites": int(n), "observed": int(observed),
        "expected": round(float(mu), 2), "std": round(float(sigma), 2),
        "z_score": round(z, 2), "enrichment": round(enrich, 3),
        "p_value": float(p_val)
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
    """Load p3k14c data, deduplicate by coordinate."""
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
                rows.append((lat, lon, continent, country))
            except (ValueError, KeyError, TypeError):
                pass
    # Deduplicate by coordinate (4 decimal places)
    seen = {}
    for lat, lon, cont, country in rows:
        key = (round(lat, 4), round(lon, 4))
        if key not in seen:
            seen[key] = (lat, lon, cont, country)
    return list(seen.values())

def load_megalithic_portal():
    """Load merged_sites.csv from the Megalithic Portal."""
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
                source = row.get("source", "").strip()
                sites.append((lat, lon, name, stype, source))
            except (ValueError, KeyError, TypeError):
                pass
    # Deduplicate by coordinate (3 decimal places)
    seen = {}
    for lat, lon, name, stype, source in sites:
        key = (round(lat, 3), round(lon, 3))
        if key not in seen:
            seen[key] = (lat, lon, name, stype, source)
    return list(seen.values())

def load_pleiades():
    """Load Pleiades data with feature type classification."""
    csv_path = os.path.join(BASE_DIR, "pleiades-places-latest.csv")
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
            if is_mon:
                monumental.append((lat, lon))
            elif is_set:
                settlements.append((lat, lon))
    return monumental, settlements

# ============================================================
# STUDY 1: PLEIADES SETTLEMENT BASELINE BY HEMISPHERE
# ============================================================
def study_1():
    print("\n" + "=" * 70)
    print("STUDY 1: PLEIADES SETTLEMENT BASELINE BY GEOGRAPHIC SEGMENT")
    print("=" * 70)

    monumental, settlements = load_pleiades()
    print(f"  Loaded {len(monumental)} monumental, {len(settlements)} settlement sites")

    segments = {
        "egypt_levant": {"lon_min": 25, "lon_max": 40, "label": "Egypt/Levant (25-40°E)"},
        "iran_mesopotamia": {"lon_min": 40, "lon_max": 60, "label": "Iran/Mesopotamia (40-60°E)"},
        "south_asia": {"lon_min": 60, "lon_max": 80, "label": "South Asia (60-80°E)"},
    }

    results = {"meta": {"date": "2026-03-17", "dataset": "Pleiades Gazetteer",
                         "threshold_km": THRESHOLD_KM, "n_trials": N_TRIALS},
               "segments": {}}

    for seg_key, seg in segments.items():
        print(f"\n  --- {seg['label']} ---")
        mon_seg = [(lat, lon) for lat, lon in monumental
                   if seg["lon_min"] <= lon < seg["lon_max"]]
        set_seg = [(lat, lon) for lat, lon in settlements
                   if seg["lon_min"] <= lon < seg["lon_max"]]

        print(f"  Monumental: {len(mon_seg)} | Settlement: {len(set_seg)}")

        if len(mon_seg) > 0:
            mon_lats = np.array([s[0] for s in mon_seg])
            mon_lons = np.array([s[1] for s in mon_seg])
            r_mon = run_mc(mon_lats, mon_lons)
            print(f"  Monumental Z = {r_mon['z_score']}")
        else:
            r_mon = {"n_sites": 0, "z_score": 0}

        if len(set_seg) > 0:
            set_lats = np.array([s[0] for s in set_seg])
            set_lons = np.array([s[1] for s in set_seg])
            r_set = run_mc(set_lats, set_lons)
            print(f"  Settlement Z = {r_set['z_score']}")
        else:
            r_set = {"n_sites": 0, "z_score": 0}

        divergence = round(r_mon["z_score"] - r_set["z_score"], 2)
        results["segments"][seg_key] = {
            "label": seg["label"],
            "monumental": r_mon, "settlement": r_set,
            "divergence": divergence
        }

    # Summary
    divs = {k: v["divergence"] for k, v in results["segments"].items()}
    peak_seg = max(divs, key=divs.get)
    results["summary"] = {
        "segment_divergences": divs,
        "peak_segment": peak_seg,
        "gradient_exists": divs.get("egypt_levant", 0) > divs.get("iran_mesopotamia", 0) > divs.get("south_asia", 0)
    }

    print(f"\n  Divergences: {divs}")
    print(f"  Peak segment: {peak_seg}")
    save_json(results, "pleiades_segment_settlement.json")
    return results


# ============================================================
# STUDY 2: MEGALITHIC PORTAL OLD WORLD vs NEW WORLD SPLIT
# ============================================================
def study_2():
    print("\n" + "=" * 70)
    print("STUDY 2: MEGALITHIC PORTAL OLD WORLD vs NEW WORLD SPLIT")
    print("=" * 70)

    sites = load_megalithic_portal()
    print(f"  Loaded {len(sites)} unique sites from Megalithic Portal")

    # Classify site types as monumental vs settlement
    MONUMENTAL_KEYWORDS = {
        "stone circle", "henge", "standing stone", "menhir", "dolmen",
        "passage grave", "burial chamber", "cairn", "chambered cairn",
        "chambered tomb", "nuraghe", "temple", "pyramid", "tumulus",
        "tomb", "sanctuary", "cursus", "stone row", "ancient cross",
        "marker stone", "holy well", "geoglyph", "carving", "petroglyph",
        "pictish symbol", "ogham stone", "round tower", "broch",
    }
    SETTLEMENT_KEYWORDS = {
        "settlement", "village", "hillfort", "fortified settlement",
        "crannog", "souterrain", "fogou", "castro",
    }

    def classify(stype):
        st = stype.lower()
        for kw in MONUMENTAL_KEYWORDS:
            if kw in st:
                return "monumental"
        for kw in SETTLEMENT_KEYWORDS:
            if kw in st:
                return "settlement"
        return "other"

    # Split by hemisphere
    old_world_noneuro = []  # lon 20E-180E OR (Africa: lon -20 to 55, lat < 35N)
    new_world = []  # Americas: lon -180 to -20

    for lat, lon, name, stype, source in sites:
        cat = classify(stype)
        if -180 <= lon < -20:
            new_world.append((lat, lon, stype, cat))
        elif lon >= 20 or (-20 <= lon < 20 and lat < 35):
            # Old World non-European: lon >= 20E, or Africa (lon -20 to 20, lat < 35N)
            if lon >= 20:
                old_world_noneuro.append((lat, lon, stype, cat))
            elif lat < 35:  # Africa between -20 and 20 longitude
                old_world_noneuro.append((lat, lon, stype, cat))

    print(f"  Old World non-European: {len(old_world_noneuro)} sites")
    print(f"  New World (Americas): {len(new_world)} sites")

    results = {"meta": {"date": "2026-03-17", "dataset": "Megalithic Portal",
                         "threshold_km": THRESHOLD_KM, "n_trials": N_TRIALS}}

    # Old World non-European
    ow_mon = [(lat, lon) for lat, lon, _, cat in old_world_noneuro if cat == "monumental"]
    ow_set = [(lat, lon) for lat, lon, _, cat in old_world_noneuro if cat == "settlement"]
    ow_all = [(lat, lon) for lat, lon, _, _ in old_world_noneuro]

    print(f"\n  Old World NE — Monumental: {len(ow_mon)}, Settlement: {len(ow_set)}, Total: {len(ow_all)}")

    if len(ow_mon) > 0:
        r_ow_mon = run_mc(np.array([s[0] for s in ow_mon]), np.array([s[1] for s in ow_mon]))
    else:
        r_ow_mon = {"n_sites": 0, "z_score": 0}
    if len(ow_set) > 0:
        r_ow_set = run_mc(np.array([s[0] for s in ow_set]), np.array([s[1] for s in ow_set]))
    else:
        r_ow_set = {"n_sites": 0, "z_score": 0}
    r_ow_all = run_mc(np.array([s[0] for s in ow_all]), np.array([s[1] for s in ow_all]))

    print(f"  Old World NE — Mon Z={r_ow_mon['z_score']}, Set Z={r_ow_set['z_score']}, All Z={r_ow_all['z_score']}")

    # New World
    nw_mon = [(lat, lon) for lat, lon, _, cat in new_world if cat == "monumental"]
    nw_set = [(lat, lon) for lat, lon, _, cat in new_world if cat == "settlement"]
    nw_all = [(lat, lon) for lat, lon, _, _ in new_world]

    print(f"\n  New World — Monumental: {len(nw_mon)}, Settlement: {len(nw_set)}, Total: {len(nw_all)}")

    if len(nw_mon) > 0:
        r_nw_mon = run_mc(np.array([s[0] for s in nw_mon]), np.array([s[1] for s in nw_mon]))
    else:
        r_nw_mon = {"n_sites": 0, "z_score": 0}
    if len(nw_set) > 0:
        r_nw_set = run_mc(np.array([s[0] for s in nw_set]), np.array([s[1] for s in nw_set]))
    else:
        r_nw_set = {"n_sites": 0, "z_score": 0}
    if len(nw_all) > 0:
        r_nw_all = run_mc(np.array([s[0] for s in nw_all]), np.array([s[1] for s in nw_all]))
    else:
        r_nw_all = {"n_sites": 0, "z_score": 0}

    print(f"  New World — Mon Z={r_nw_mon['z_score']}, Set Z={r_nw_set['z_score']}, All Z={r_nw_all['z_score']}")

    ow_div = round(r_ow_mon["z_score"] - r_ow_set["z_score"], 2)
    nw_div = round(r_nw_mon["z_score"] - r_nw_set["z_score"], 2)

    results["old_world_non_european"] = {
        "monumental": r_ow_mon, "settlement": r_ow_set, "all": r_ow_all,
        "divergence": ow_div
    }
    results["new_world"] = {
        "monumental": r_nw_mon, "settlement": r_nw_set, "all": r_nw_all,
        "divergence": nw_div
    }
    results["comparison"] = {
        "old_world_divergence": ow_div, "new_world_divergence": nw_div,
        "divergence_ratio": round(ow_div / nw_div, 2) if nw_div != 0 else "inf"
    }

    print(f"\n  OW divergence: {ow_div} | NW divergence: {nw_div}")
    save_json(results, "portal_hemisphere_settlement.json")
    return results


# ============================================================
# STUDY 3: P3K14C CONTINENTAL SIGNAL DECOMPOSITION
# ============================================================
def study_3():
    print("\n" + "=" * 70)
    print("STUDY 3: P3K14C CONTINENTAL SIGNAL DECOMPOSITION")
    print("=" * 70)

    p3k_sites = load_p3k14c()
    print(f"  Loaded {len(p3k_sites)} unique p3k14c sites")

    # Group by continent
    by_continent = defaultdict(list)
    for lat, lon, cont, country in p3k_sites:
        by_continent[cont].append((lat, lon, country))

    target_continents = ["Africa", "Asia", "South America"]
    results = {"meta": {"date": "2026-03-17", "dataset": "p3k14c",
                         "threshold_km": THRESHOLD_KM, "n_trials": N_TRIALS},
               "continents": {}}

    for cont in target_continents:
        cdata = by_continent.get(cont, [])
        print(f"\n  --- {cont} ({len(cdata)} sites) ---")
        if len(cdata) == 0:
            results["continents"][cont] = {"n_sites": 0, "z_score": 0}
            continue
        lats = np.array([s[0] for s in cdata])
        lons = np.array([s[1] for s in cdata])
        r = run_mc(lats, lons)
        print(f"  {cont}: Z = {r['z_score']}, obs = {r['observed']}, exp = {r['expected']}")
        results["continents"][cont] = r

    # Africa sub-regions
    print(f"\n  --- Africa Sub-Regions ---")
    africa = by_continent.get("Africa", [])
    sub_regions = {
        "north_africa_egypt": [(lat, lon) for lat, lon, c in africa
                                if lat > 20 and 25 <= lon <= 35],
        "negev_levant": [(lat, lon) for lat, lon, c in africa
                          if lat > 29 and 34 <= lon <= 36],
        "sub_saharan": [(lat, lon) for lat, lon, c in africa
                         if lat <= 20],
        "rest_of_north_africa": [(lat, lon) for lat, lon, c in africa
                                  if lat > 20 and not (25 <= lon <= 35)],
    }

    results["africa_subregions"] = {}
    for region_key, region_sites in sub_regions.items():
        print(f"  {region_key}: {len(region_sites)} sites")
        if len(region_sites) > 0:
            lats = np.array([s[0] for s in region_sites])
            lons = np.array([s[1] for s in region_sites])
            r = run_mc(lats, lons)
            print(f"    Z = {r['z_score']}")
        else:
            r = {"n_sites": 0, "z_score": 0}
        results["africa_subregions"][region_key] = r

    save_json(results, "p3k14c_continental_decomposition.json")
    return results


# ============================================================
# STUDY 4: OLD WORLD CORRIDOR ENRICHMENT PROFILE
# ============================================================
def study_4():
    print("\n" + "=" * 70)
    print("STUDY 4: OLD WORLD CORRIDOR ENRICHMENT PROFILE")
    print("=" * 70)

    sites = load_megalithic_portal()

    # Sample points along the great circle in the Old World corridor (25°E to 110°E)
    # First, compute great circle points
    pole_lat_r = math.radians(POLE_LAT)
    pole_lon_r = math.radians(POLE_LNG)

    # Generate points along the great circle for every degree of arc
    gc_points = []
    for deg in range(360):
        angle = math.radians(deg)
        # Point on great circle = rotate pole by 90° around the equatorial plane
        # Using the parametric form:
        # The great circle with pole P is the set of points at distance QUARTER_CIRC from P
        lat = math.asin(
            math.sin(pole_lat_r) * math.cos(math.pi/2) +
            math.cos(pole_lat_r) * math.sin(math.pi/2) * math.cos(angle)
        )
        lon = pole_lon_r + math.atan2(
            math.sin(angle) * math.sin(math.pi/2) * math.cos(pole_lat_r),
            math.cos(math.pi/2) - math.sin(pole_lat_r) * math.sin(lat)
        )
        lat_d = math.degrees(lat)
        lon_d = math.degrees(lon)
        # Normalize longitude
        if lon_d > 180: lon_d -= 360
        if lon_d < -180: lon_d += 360
        gc_points.append((lat_d, lon_d, deg))

    # Filter to Old World corridor: longitude 25°E to 110°E
    ow_gc_points = [(lat, lon, deg) for lat, lon, deg in gc_points if 25 <= lon <= 110]
    ow_gc_points.sort(key=lambda x: x[1])  # sort by longitude
    print(f"  Great circle points in Old World corridor (25-110°E): {len(ow_gc_points)}")

    # For each point, count Megalithic Portal sites within 50 km
    site_lats = np.array([s[0] for s in sites])
    site_lons = np.array([s[1] for s in sites])

    profile = []
    for gc_lat, gc_lon, deg in ow_gc_points:
        # Count sites within 50 km of this GC point
        dists = haversine_vec(gc_lat, gc_lon, site_lats, site_lons)
        nearby = int(np.sum(dists <= THRESHOLD_KM))

        # Also compute local background: sites within 500 km for density estimation
        bg_count = int(np.sum(dists <= 500))

        profile.append({
            "gc_lat": round(gc_lat, 3), "gc_lon": round(gc_lon, 3),
            "arc_degree": deg,
            "sites_within_50km": nearby,
            "sites_within_500km": bg_count,
        })

    # Compute residuals: expected = (sites_500km / area_500km) * area_50km
    area_50 = math.pi * 50**2
    area_500 = math.pi * 500**2
    for p in profile:
        bg_density = p["sites_within_500km"] / area_500
        p["expected_50km"] = round(bg_density * area_50, 3)
        p["residual"] = round(p["sites_within_50km"] - p["expected_50km"], 3)

    # Identify peaks
    peak_threshold = 2  # residual > 2 sites
    peaks = [p for p in profile if p["residual"] > peak_threshold]

    # Cluster peaks by longitude proximity
    clusters = []
    if peaks:
        current_cluster = [peaks[0]]
        for i in range(1, len(peaks)):
            if peaks[i]["gc_lon"] - peaks[i-1]["gc_lon"] < 5:  # within 5° longitude
                current_cluster.append(peaks[i])
            else:
                clusters.append(current_cluster)
                current_cluster = [peaks[i]]
        clusters.append(current_cluster)

    cluster_summary = []
    for cl in clusters:
        lon_range = (cl[0]["gc_lon"], cl[-1]["gc_lon"])
        max_residual = max(p["residual"] for p in cl)
        total_sites = sum(p["sites_within_50km"] for p in cl)
        cluster_summary.append({
            "lon_range": [round(lon_range[0], 1), round(lon_range[1], 1)],
            "n_points": len(cl),
            "max_residual": round(max_residual, 2),
            "total_sites_50km": total_sites,
        })

    results = {
        "meta": {"date": "2026-03-17", "dataset": "Megalithic Portal",
                 "corridor": "25E to 110E longitude", "threshold_km": THRESHOLD_KM},
        "profile": profile,
        "peaks_above_threshold": len(peaks),
        "peak_clusters": cluster_summary,
        "n_gc_points_in_corridor": len(ow_gc_points),
    }

    print(f"  Profile computed: {len(profile)} points")
    print(f"  Peaks (residual > {peak_threshold}): {len(peaks)}")
    print(f"  Peak clusters: {len(cluster_summary)}")
    for i, cl in enumerate(cluster_summary):
        print(f"    Cluster {i+1}: lon {cl['lon_range']}, max_residual={cl['max_residual']}, sites={cl['total_sites_50km']}")

    save_json(results, "old_world_enrichment_profile.json")
    return results


# ============================================================
# STUDY 5: NEW WORLD TYPE ENRICHMENT
# ============================================================
def study_5():
    print("\n" + "=" * 70)
    print("STUDY 5: NEW WORLD TYPE ENRICHMENT")
    print("=" * 70)

    sites = load_megalithic_portal()

    # New World sites: Americas + Easter Island (lon -180 to -20, OR Easter Island ~-109)
    nw_sites = [(lat, lon, name, stype) for lat, lon, name, stype, source in sites
                if lon < -20 or (abs(lat - (-27.1)) < 1 and abs(lon - (-109.3)) < 1)]

    print(f"  New World sites: {len(nw_sites)}")

    # Compute which are on the line
    on_line = []
    off_line = []
    for lat, lon, name, stype in nw_sites:
        d = scalar_gc_dist(lat, lon)
        if d <= THRESHOLD_KM:
            on_line.append((lat, lon, name, stype))
        else:
            off_line.append((lat, lon, name, stype))

    print(f"  On-line (within {THRESHOLD_KM}km): {len(on_line)}")
    print(f"  Off-line: {len(off_line)}")

    # Count by type
    all_types_on = Counter(stype for _, _, _, stype in on_line)
    all_types_off = Counter(stype for _, _, _, stype in off_line)
    all_types_total = Counter(stype for _, _, _, stype in nw_sites)

    # Overall fraction on the line
    overall_frac = len(on_line) / len(nw_sites) if len(nw_sites) > 0 else 0

    # Per-type enrichment
    type_results = {}
    for stype in sorted(all_types_total.keys()):
        total = all_types_total[stype]
        on = all_types_on.get(stype, 0)
        off = all_types_off.get(stype, 0)
        frac_on = on / total if total > 0 else 0
        expected_on = overall_frac * total
        enrichment = frac_on / overall_frac if overall_frac > 0 else 0

        type_results[stype] = {
            "total_new_world": total,
            "on_line": on,
            "off_line": off,
            "fraction_on_line": round(frac_on, 4),
            "expected_on_line": round(expected_on, 2),
            "enrichment_vs_background": round(enrichment, 3),
        }

    # Sort by enrichment
    enriched_types = sorted(type_results.items(),
                            key=lambda x: x[1]["enrichment_vs_background"], reverse=True)

    # Run MC for specific types with enough data
    mc_by_type = {}
    for stype, info in enriched_types[:10]:  # top 10 enriched types
        type_sites = [(lat, lon) for lat, lon, _, st in nw_sites if st == stype]
        if len(type_sites) >= 5:  # need at least 5 sites for meaningful MC
            lats = np.array([s[0] for s in type_sites])
            lons = np.array([s[1] for s in type_sites])
            r = run_mc(lats, lons)
            mc_by_type[stype] = r
            print(f"  {stype} ({info['total_new_world']} sites): Z={r['z_score']}, on_line={info['on_line']}")

    results = {
        "meta": {"date": "2026-03-17", "dataset": "Megalithic Portal (New World)",
                 "threshold_km": THRESHOLD_KM, "n_trials": N_TRIALS},
        "overall": {
            "total_new_world": len(nw_sites),
            "on_line": len(on_line), "off_line": len(off_line),
            "fraction_on_line": round(overall_frac, 4),
        },
        "type_enrichment": type_results,
        "mc_by_type": mc_by_type,
        "top_enriched": [{"type": stype, **info} for stype, info in enriched_types[:15]],
    }

    save_json(results, "new_world_type_enrichment.json")
    return results


# ============================================================
# STUDY 6: ICE AGE (LGM) CONNECTIVITY TEST
# ============================================================
def study_6():
    print("\n" + "=" * 70)
    print("STUDY 6: LGM CONNECTIVITY TEST")
    print("=" * 70)

    # --- Megalithic Portal ---
    mp_sites = load_megalithic_portal()
    p3k_sites = load_p3k14c()

    # Define LGM-connected vs ocean-separated regions
    lgm_connected_bounds = {
        "egypt_levant": {"lat": (20, 35), "lon": (25, 40)},
        "iran": {"lat": (25, 40), "lon": (40, 65)},
        "indus_valley": {"lat": (20, 35), "lon": (65, 80)},
        "se_asia": {"lat": (-10, 25), "lon": (80, 120)},
    }
    ocean_separated_bounds = {
        "peru_andes": {"lat": (-25, 5), "lon": (-85, -65)},
        "easter_island": {"lat": (-30, -25), "lon": (-112, -107)},
    }

    def in_region(lat, lon, bounds):
        return (bounds["lat"][0] <= lat <= bounds["lat"][1] and
                bounds["lon"][0] <= lon <= bounds["lon"][1])

    def classify_connectivity(lat, lon):
        for name, bounds in lgm_connected_bounds.items():
            if in_region(lat, lon, bounds):
                return "lgm_connected", name
        for name, bounds in ocean_separated_bounds.items():
            if in_region(lat, lon, bounds):
                return "ocean_separated", name
        return None, None

    # Megalithic Portal split
    mp_lgm = []; mp_ocean = []
    mp_lgm_by_region = defaultdict(list)
    mp_ocean_by_region = defaultdict(list)

    for lat, lon, name, stype, source in mp_sites:
        conn, region = classify_connectivity(lat, lon)
        if conn == "lgm_connected":
            mp_lgm.append((lat, lon))
            mp_lgm_by_region[region].append((lat, lon))
        elif conn == "ocean_separated":
            mp_ocean.append((lat, lon))
            mp_ocean_by_region[region].append((lat, lon))

    print(f"\n  Megalithic Portal:")
    print(f"  LGM-connected: {len(mp_lgm)} sites")
    for r, s in mp_lgm_by_region.items():
        print(f"    {r}: {len(s)}")
    print(f"  Ocean-separated: {len(mp_ocean)} sites")
    for r, s in mp_ocean_by_region.items():
        print(f"    {r}: {len(s)}")

    results = {"meta": {"date": "2026-03-17", "threshold_km": THRESHOLD_KM, "n_trials": N_TRIALS},
               "megalithic_portal": {}, "p3k14c": {}}

    # MC for each group (Megalithic Portal)
    for label, group_sites in [("lgm_connected", mp_lgm), ("ocean_separated", mp_ocean)]:
        if len(group_sites) > 0:
            lats = np.array([s[0] for s in group_sites])
            lons = np.array([s[1] for s in group_sites])
            r = run_mc(lats, lons)
            print(f"  MP {label}: Z={r['z_score']}")
        else:
            r = {"n_sites": 0, "z_score": 0}
        results["megalithic_portal"][label] = r

    # Per-region MC (Megalithic Portal)
    results["megalithic_portal"]["by_region"] = {}
    for region in list(lgm_connected_bounds.keys()) + list(ocean_separated_bounds.keys()):
        rsites = mp_lgm_by_region.get(region, []) + mp_ocean_by_region.get(region, [])
        if len(rsites) > 0:
            lats = np.array([s[0] for s in rsites])
            lons = np.array([s[1] for s in rsites])
            r = run_mc(lats, lons)
            print(f"    MP {region}: Z={r['z_score']} (n={len(rsites)})")
        else:
            r = {"n_sites": 0, "z_score": 0}
        results["megalithic_portal"]["by_region"][region] = r

    # p3k14c split
    p3k_lgm = defaultdict(list); p3k_ocean = defaultdict(list)
    for lat, lon, cont, country in p3k_sites:
        conn, region = classify_connectivity(lat, lon)
        if conn == "lgm_connected":
            p3k_lgm[region].append((lat, lon))
        elif conn == "ocean_separated":
            p3k_ocean[region].append((lat, lon))

    all_p3k_lgm = [s for sites in p3k_lgm.values() for s in sites]
    all_p3k_ocean = [s for sites in p3k_ocean.values() for s in sites]

    print(f"\n  p3k14c:")
    print(f"  LGM-connected: {len(all_p3k_lgm)} sites")
    print(f"  Ocean-separated: {len(all_p3k_ocean)} sites")

    for label, group_sites in [("lgm_connected", all_p3k_lgm), ("ocean_separated", all_p3k_ocean)]:
        if len(group_sites) > 0:
            lats = np.array([s[0] for s in group_sites])
            lons = np.array([s[1] for s in group_sites])
            r = run_mc(lats, lons)
            print(f"  p3k {label}: Z={r['z_score']}")
        else:
            r = {"n_sites": 0, "z_score": 0}
        results["p3k14c"][label] = r

    # Distance gradient test (LGM-connected regions ordered by distance from Egypt)
    region_order = ["egypt_levant", "iran", "indus_valley", "se_asia"]
    gradient = []
    for region in region_order:
        mp_r = results["megalithic_portal"]["by_region"].get(region, {})
        gradient.append({"region": region, "z_score": mp_r.get("z_score", 0),
                          "n_sites": mp_r.get("n_sites", 0)})

    z_values = [g["z_score"] for g in gradient if g["n_sites"] > 0]
    gradient_exists = len(z_values) >= 2 and z_values == sorted(z_values, reverse=True)

    results["distance_gradient"] = {
        "region_order": gradient,
        "gradient_exists": gradient_exists,
        "description": "Regions ordered by distance from Egypt (west to east)"
    }

    # Comparison
    lgm_z = results["megalithic_portal"]["lgm_connected"].get("z_score", 0)
    ocean_z = results["megalithic_portal"]["ocean_separated"].get("z_score", 0)
    results["comparison"] = {
        "lgm_z": lgm_z, "ocean_z": ocean_z,
        "lgm_stronger": lgm_z > ocean_z,
        "z_difference": round(lgm_z - ocean_z, 2)
    }

    print(f"\n  LGM Z={lgm_z} vs Ocean Z={ocean_z}")
    print(f"  Gradient test: {gradient_exists}")
    save_json(results, "lgm_connectivity_test.json")
    return results


# ============================================================
# STUDY 7: "BORING EXPLANATION" QUANTIFICATION
# ============================================================
def study_7():
    print("\n" + "=" * 70)
    print("STUDY 7: HABITABILITY-ADJUSTED SIGNAL")
    print("=" * 70)

    p3k_sites = load_p3k14c()
    print(f"  Loaded {len(p3k_sites)} p3k14c sites")

    all_lats = np.array([s[0] for s in p3k_sites])
    all_lons = np.array([s[1] for s in p3k_sites])

    # Step 1: Does the Alison circle pass through higher-density human occupation?
    print("\n  Step 1: Human occupation density along the Alison circle vs random circles")

    # Compute density along the Alison circle
    alison_dists = gc_dist_vec(POLE_LAT, POLE_LNG, all_lats, all_lons)
    alison_count = int(np.sum(alison_dists <= THRESHOLD_KM))

    # Generate 1000 random great circles and count sites within 50 km
    n_random_circles = 1000
    random_counts = []
    print(f"  Running {n_random_circles} random circles...")

    for i in range(n_random_circles):
        # Generate random pole (uniform on sphere)
        z = np.random.uniform(-1, 1)
        rand_pole_lat = math.degrees(math.asin(z))
        rand_pole_lon = np.random.uniform(-180, 180)

        dists = gc_dist_vec(rand_pole_lat, rand_pole_lon, all_lats, all_lons)
        random_counts.append(int(np.sum(dists <= THRESHOLD_KM)))

        if (i + 1) % 200 == 0:
            print(f"    {i+1}/{n_random_circles} circles done")

    mu_rand = np.mean(random_counts)
    sigma_rand = np.std(random_counts)
    z_habitability = (alison_count - mu_rand) / sigma_rand if sigma_rand > 0 else 0
    percentile = np.sum(np.array(random_counts) < alison_count) / len(random_counts) * 100

    print(f"  Alison circle: {alison_count} sites within {THRESHOLD_KM}km")
    print(f"  Random circles: {mu_rand:.1f} ± {sigma_rand:.1f}")
    print(f"  Habitability Z: {z_habitability:.2f}")
    print(f"  Percentile: {percentile:.1f}%")

    # Step 2: Habitability-adjusted Z-score
    # Use density-weighted random circles as baseline
    # Each random circle is itself a sample of "where the line passes through populated areas"
    # So the distribution of random circle counts IS the habitability baseline
    # The Z_adjusted = (alison_count - mu_random_circles) / sigma_random_circles
    # This is exactly z_habitability computed above!

    # But we should also compare to distribution-matched MC
    # Distribution-matched MC accounts for the geographic DISTRIBUTION of sites
    # Random circles account for the DENSITY along any circle
    dm_result = run_mc(all_lats, all_lons)
    print(f"\n  Distribution-matched MC: Z={dm_result['z_score']}")
    print(f"  Random-circle baseline:  Z={z_habitability:.2f}")
    print(f"  Difference: {dm_result['z_score'] - z_habitability:.2f}")

    # Step 3: Residual signal
    z_adjusted = round(z_habitability, 2)
    z_standard = dm_result["z_score"]

    if z_adjusted > 3:
        diagnosis = "SIGNIFICANT: Signal survives habitability correction. Something beyond population density drives the pattern."
    elif z_adjusted > 2:
        diagnosis = "MARGINAL: Signal weakened but not eliminated by habitability correction."
    elif z_adjusted > 0:
        diagnosis = "WEAK: Most signal explained by habitability corridor, small residual remains."
    else:
        diagnosis = "NULL: Pattern fully explained by habitable corridors. Circle passes through where people lived."

    results = {
        "meta": {"date": "2026-03-17", "dataset": "p3k14c",
                 "threshold_km": THRESHOLD_KM, "n_random_circles": n_random_circles,
                 "n_mc_trials": N_TRIALS},
        "alison_circle": {
            "sites_within_threshold": alison_count,
        },
        "random_circle_baseline": {
            "mean": round(float(mu_rand), 2),
            "std": round(float(sigma_rand), 2),
            "z_habitability": round(float(z_habitability), 2),
            "percentile": round(float(percentile), 1),
            "n_circles": n_random_circles,
        },
        "distribution_matched_baseline": dm_result,
        "comparison": {
            "z_standard": z_standard,
            "z_habitability_adjusted": z_adjusted,
            "z_reduction": round(z_standard - z_adjusted, 2),
            "fraction_explained_by_habitability": round(
                (z_standard - z_adjusted) / z_standard, 3) if z_standard > 0 else 0,
        },
        "diagnosis": diagnosis,
    }

    print(f"\n  Z standard: {z_standard}")
    print(f"  Z adjusted: {z_adjusted}")
    print(f"  Diagnosis: {diagnosis}")

    save_json(results, "habitability_adjusted_signal.json")
    return results


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    t0 = time.time()
    print("=" * 70)
    print("HEMISPHERE DECOMPOSITION DIRECTIVE")
    print("=" * 70)

    r1 = study_1()
    r2 = study_2()
    r3 = study_3()
    r4 = study_4()
    r5 = study_5()
    r6 = study_6()
    r7 = study_7()

    elapsed = time.time() - t0
    print(f"\n{'=' * 70}")
    print(f"ALL STUDIES COMPLETE — {elapsed:.0f}s elapsed")
    print(f"{'=' * 70}")
    print(f"Output directory: {OUT_DIR}")
