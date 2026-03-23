#!/usr/bin/env python3
"""
Collinearity Resolution Tests
==============================
Seven tests probing the MECHANISM behind Giza-Nazca-Easter Island collinearity.

1. Constraint Propagation — given two anchors, where must the third be?
2. Continental Position — do Africa-Peru-Pacific zones force collinearity?
3. Random Civilization Simulation — monument-weighted triplet sampling
4. Habitable Island Scan — is Easter Island the only option on the circle?
5. Geological Source Alignment — are tectonic features collinear too?
6. Atmospheric Circulation — does the circle track the Hadley cell boundary?
7. Equivalent Site Substitution — do culturally equivalent sites also align?
"""

import json, math, os, sys, time
import numpy as np
import netCDF4 as nc
from pathlib import Path
from itertools import product as iterproduct
from scipy.stats import percentileofscore

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

# ============================================================
# CONFIGURATION
# ============================================================
R = 6371.0
POLE_LAT = 59.682122
POLE_LON = -138.646087

GIZA = (29.979, 31.134)
NAZCA = (-14.700, -75.100)
EASTER = (-27.100, -109.333)

BASE = Path(os.path.expanduser("~/megalith_site_research"))
OUT = BASE / "outputs" / "collinearity_resolution"
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# LOAD ETOPO1
# ============================================================
print("=" * 70)
print("LOADING ETOPO1")
print("=" * 70)

etopo_path = BASE / "data" / "geophysical" / "etopo" / "ETOPO1_Ice_g_gmt4.grd"
ds = nc.Dataset(str(etopo_path))
etopo_lons = ds.variables['x'][:]
etopo_lats = ds.variables['y'][:]
etopo_z = ds.variables['z']

lon_min_e, lon_step = float(etopo_lons[0]), float(etopo_lons[1] - etopo_lons[0])
lat_min_e, lat_step = float(etopo_lats[0]), float(etopo_lats[1] - etopo_lats[0])
n_lats = len(etopo_lats)
n_lons = len(etopo_lons)
print(f"ETOPO1: {n_lats} x {n_lons}, step={lat_step:.4f}°")


# ============================================================
# UTILITY FUNCTIONS
# ============================================================
def get_elev(lat, lon):
    li = int(round((lat - lat_min_e) / lat_step))
    lo = int(round((lon - lon_min_e) / lon_step))
    li = max(0, min(li, n_lats - 1))
    lo = max(0, min(lo, n_lons - 1))
    return float(etopo_z[li, lo])


def get_elev_batch(lats, lons):
    lat_idx = np.round((np.asarray(lats) - lat_min_e) / lat_step).astype(int)
    lon_idx = np.round((np.asarray(lons) - lon_min_e) / lon_step).astype(int)
    lat_idx = np.clip(lat_idx, 0, n_lats - 1)
    lon_idx = np.clip(lon_idx, 0, n_lons - 1)
    elevs = np.empty(len(lat_idx))
    for i in range(len(lat_idx)):
        elevs[i] = float(etopo_z[lat_idx[i], lon_idx[i]])
    return elevs


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    return R * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def latlon_to_cart(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.array([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])


def cart_to_latlon(v):
    v = v / np.linalg.norm(v)
    lat = np.degrees(np.arcsin(np.clip(v[2], -1, 1)))
    lon = np.degrees(np.arctan2(v[1], v[0]))
    return lat, lon


def pole_from_two_sites(lat1, lon1, lat2, lon2):
    """Compute the great circle pole from two sites."""
    v1 = latlon_to_cart(lat1, lon1)
    v2 = latlon_to_cart(lat2, lon2)
    pole = np.cross(v1, v2)
    n = np.linalg.norm(pole)
    if n < 1e-12:
        return None
    return pole / n


def distance_to_gc(lat, lon, pole):
    """Distance in km from a point to the great circle defined by pole (unit vector)."""
    v = latlon_to_cart(lat, lon)
    dot = np.clip(np.dot(v, pole), -1, 1)
    ang = np.arccos(dot)
    return abs(ang - np.pi/2) * R


def trace_great_circle(pole, step_deg=0.1):
    """Trace GC at given step. Returns list of (lat, lon) tuples."""
    if abs(pole[2]) < 0.9:
        u = np.cross(pole, [0, 0, 1])
    else:
        u = np.cross(pole, [1, 0, 0])
    u = u / np.linalg.norm(u)
    v = np.cross(pole, u)
    v = v / np.linalg.norm(v)

    n_pts = int(360.0 / step_deg)
    thetas = np.linspace(0, 2*np.pi, n_pts, endpoint=False)
    pts = np.outer(np.cos(thetas), u) + np.outer(np.sin(thetas), v)
    lats = np.degrees(np.arcsin(np.clip(pts[:, 2], -1, 1)))
    lons = np.degrees(np.arctan2(pts[:, 1], pts[:, 0]))
    return list(zip(lats.tolist(), lons.tolist()))


def random_pole():
    z = np.random.uniform(-1, 1)
    theta = np.random.uniform(0, 2*np.pi)
    lat = np.degrees(np.arcsin(z))
    lon = np.degrees(theta) - 180
    return latlon_to_cart(lat, lon)


def compute_triplet_fit(lat1, lon1, lat2, lon2, lat3, lon3):
    """Compute mean distance (km) of three sites from their best-fit great circle."""
    v1 = latlon_to_cart(lat1, lon1)
    v2 = latlon_to_cart(lat2, lon2)
    v3 = latlon_to_cart(lat3, lon3)
    # SVD to find best-fit plane
    M = np.vstack([v1, v2, v3])
    M_centered = M - M.mean(axis=0)
    _, s, Vt = np.linalg.svd(M_centered)
    pole = Vt[-1]  # normal to best-fit plane
    pole = pole / np.linalg.norm(pole)
    # Mean distance of three points from this great circle
    dists = []
    for v in [v1, v2, v3]:
        dot = np.clip(np.dot(v, pole), -1, 1)
        ang = np.arccos(dot)
        dists.append(abs(ang - np.pi/2) * R)
    return np.mean(dists)


# ============================================================
# PACIFIC ISLANDS DATABASE
# ============================================================
PACIFIC_ISLANDS = [
    {"name": "Easter Island (Rapa Nui)", "lat": -27.100, "lon": -109.333, "population": 7750},
    {"name": "Pitcairn Island", "lat": -25.07, "lon": -130.10, "population": 47},
    {"name": "Mangareva (Gambier)", "lat": -23.12, "lon": -134.97, "population": 1239},
    {"name": "Raivavae", "lat": -23.87, "lon": -147.67, "population": 940},
    {"name": "Rapa Iti", "lat": -27.62, "lon": -144.34, "population": 507},
    {"name": "Tahiti", "lat": -17.68, "lon": -149.45, "population": 189517},
    {"name": "Moorea", "lat": -17.53, "lon": -149.83, "population": 17816},
    {"name": "Bora Bora", "lat": -16.50, "lon": -151.74, "population": 10605},
    {"name": "Huahine", "lat": -16.75, "lon": -150.99, "population": 6303},
    {"name": "Raiatea", "lat": -16.83, "lon": -151.43, "population": 12237},
    {"name": "Rangiroa", "lat": -15.13, "lon": -147.65, "population": 2567},
    {"name": "Hao", "lat": -18.07, "lon": -140.95, "population": 1344},
    {"name": "Nuku Hiva (Marquesas)", "lat": -8.86, "lon": -140.10, "population": 3120},
    {"name": "Hiva Oa (Marquesas)", "lat": -9.77, "lon": -139.00, "population": 2190},
    {"name": "Fatu Hiva (Marquesas)", "lat": -10.47, "lon": -138.67, "population": 611},
    {"name": "Rarotonga (Cook Is.)", "lat": -21.23, "lon": -159.78, "population": 13007},
    {"name": "Aitutaki (Cook Is.)", "lat": -18.86, "lon": -159.79, "population": 2000},
    {"name": "Niue", "lat": -19.05, "lon": -169.87, "population": 1620},
    {"name": "Upolu (Samoa)", "lat": -13.83, "lon": -171.76, "population": 143418},
    {"name": "Savaii (Samoa)", "lat": -13.63, "lon": -172.45, "population": 43142},
    {"name": "Tongatapu (Tonga)", "lat": -21.21, "lon": -175.15, "population": 75416},
    {"name": "Vavau (Tonga)", "lat": -18.62, "lon": -173.98, "population": 14922},
    {"name": "Viti Levu (Fiji)", "lat": -17.77, "lon": 178.07, "population": 600000},
    {"name": "Vanua Levu (Fiji)", "lat": -16.58, "lon": 179.22, "population": 135961},
    {"name": "Hawaii (Big Island)", "lat": 19.50, "lon": -155.50, "population": 200629},
    {"name": "Maui", "lat": 20.80, "lon": -156.32, "population": 164221},
    {"name": "Oahu", "lat": 21.47, "lon": -157.98, "population": 1016508},
    {"name": "Kauai", "lat": 22.08, "lon": -159.52, "population": 72293},
    {"name": "Galapagos (Santa Cruz)", "lat": -0.74, "lon": -90.30, "population": 15000},
    {"name": "Juan Fernandez (Robinson Crusoe)", "lat": -33.64, "lon": -78.83, "population": 926},
    {"name": "Sala y Gomez", "lat": -26.47, "lon": -105.37, "population": 0},
    {"name": "Norfolk Island", "lat": -29.04, "lon": 167.95, "population": 1748},
    {"name": "Lord Howe Island", "lat": -31.55, "lon": 159.08, "population": 382},
    {"name": "New Caledonia (Noumea)", "lat": -22.28, "lon": 166.46, "population": 271407},
    {"name": "Wallis Island", "lat": -13.30, "lon": -176.20, "population": 8333},
    {"name": "Futuna", "lat": -14.27, "lon": -178.15, "population": 3613},
    {"name": "Tokelau (Atafu)", "lat": -8.58, "lon": -172.50, "population": 524},
    {"name": "Tuvalu (Funafuti)", "lat": -8.52, "lon": 179.20, "population": 6320},
    {"name": "Kiribati (Tarawa)", "lat": 1.45, "lon": 173.00, "population": 63439},
    {"name": "Marshall Islands (Majuro)", "lat": 7.09, "lon": 171.38, "population": 27797},
    {"name": "Palau (Koror)", "lat": 7.34, "lon": 134.47, "population": 11754},
    {"name": "Guam", "lat": 13.44, "lon": 144.79, "population": 168783},
    {"name": "Saipan (CNMI)", "lat": 15.19, "lon": 145.75, "population": 43385},
    {"name": "Pohnpei (Micronesia)", "lat": 6.88, "lon": 158.23, "population": 36196},
    {"name": "Kosrae (Micronesia)", "lat": 5.32, "lon": 162.98, "population": 6616},
    {"name": "Yap (Micronesia)", "lat": 9.54, "lon": 138.13, "population": 11377},
    {"name": "Chuuk (Micronesia)", "lat": 7.45, "lon": 151.85, "population": 48654},
    {"name": "Nauru", "lat": -0.52, "lon": 166.93, "population": 10876},
    {"name": "Henderson Island", "lat": -24.37, "lon": -128.33, "population": 0},
    {"name": "Ducie Island", "lat": -24.67, "lon": -124.78, "population": 0},
    {"name": "Clipperton Island", "lat": 10.30, "lon": -109.22, "population": 0},
    {"name": "Howland Island", "lat": 0.81, "lon": -176.62, "population": 0},
    {"name": "Baker Island", "lat": 0.19, "lon": -176.48, "population": 0},
    {"name": "Jarvis Island", "lat": -0.37, "lon": -160.01, "population": 0},
    {"name": "Palmyra Atoll", "lat": 5.88, "lon": -162.08, "population": 0},
    {"name": "Johnston Atoll", "lat": 16.73, "lon": -169.53, "population": 0},
    {"name": "Wake Island", "lat": 19.28, "lon": 166.63, "population": 0},
    {"name": "Midway Atoll", "lat": 28.21, "lon": -177.38, "population": 0},
    {"name": "Tubuai", "lat": -23.35, "lon": -149.48, "population": 2170},
    {"name": "Rurutu", "lat": -22.47, "lon": -151.35, "population": 2466},
    {"name": "Rimatara", "lat": -22.64, "lon": -152.81, "population": 873},
    {"name": "Maupiti", "lat": -16.45, "lon": -152.27, "population": 1194},
    {"name": "Tikehau", "lat": -15.00, "lon": -148.17, "population": 500},
    {"name": "Makatea", "lat": -15.83, "lon": -148.25, "population": 84},
    {"name": "Fakarava", "lat": -16.05, "lon": -145.65, "population": 855},
    {"name": "Anaa", "lat": -17.35, "lon": -145.50, "population": 411},
    {"name": "Tureia", "lat": -20.78, "lon": -138.57, "population": 313},
    {"name": "Reao", "lat": -18.47, "lon": -136.42, "population": 367},
    {"name": "Pukarua", "lat": -18.30, "lon": -137.02, "population": 168},
    {"name": "Tatakoto", "lat": -17.35, "lon": -138.45, "population": 264},
    {"name": "San Felix Island", "lat": -26.28, "lon": -80.13, "population": 0},
    {"name": "San Ambrosio Island", "lat": -26.35, "lon": -79.88, "population": 0},
    {"name": "Desventuradas Islands", "lat": -26.35, "lon": -80.05, "population": 0},
]


# ============================================================
# TEST 1: CONSTRAINT PROPAGATION
# ============================================================
def test_1_constraint_propagation():
    print("\n" + "=" * 70)
    print("TEST 1: CONSTRAINT PROPAGATION")
    print("Given two anchors, how much habitable land does their GC cross")
    print("in the region where the third must be?")
    print("=" * 70)

    pairs = [
        ("Giza+Nazca→Pacific", GIZA, NAZCA, "Pacific",
         {"lon_min": -180, "lon_max": -70, "lat_min": -50, "lat_max": 30}),
        ("Giza+Easter→South_America", GIZA, EASTER, "South America",
         {"lon_min": -85, "lon_max": -65, "lat_min": -30, "lat_max": 5}),
        ("Nazca+Easter→NorthAfrica_ME", NAZCA, EASTER, "North Africa/ME",
         {"lon_min": -10, "lon_max": 70, "lat_min": 15, "lat_max": 40}),
    ]

    results = {}
    for name, site1, site2, target, bounds in pairs:
        print(f"\n--- {name} ---")
        pole = pole_from_two_sites(site1[0], site1[1], site2[0], site2[1])
        if pole is None:
            print("  ERROR: degenerate pole")
            continue

        # Trace circle at 0.1° resolution
        circle = trace_great_circle(pole, step_deg=0.1)

        # Filter to region
        region_pts = [(lat, lon) for lat, lon in circle
                      if bounds["lon_min"] <= lon <= bounds["lon_max"]
                      and bounds["lat_min"] <= lat <= bounds["lat_max"]]

        # Check land
        if region_pts:
            rp_lats = [p[0] for p in region_pts]
            rp_lons = [p[1] for p in region_pts]
            elevs = get_elev_batch(rp_lats, rp_lons)
            land_pts = [(lat, lon) for (lat, lon), e in zip(region_pts, elevs) if e > 0]
        else:
            land_pts = []

        land_frac = len(land_pts) / max(len(region_pts), 1)
        print(f"  Circle points in region: {len(region_pts)}")
        print(f"  Land points: {len(land_pts)}")
        print(f"  Land fraction: {land_frac:.1%}")

        pair_result = {
            "pair": name,
            "target_region": target,
            "circle_points_in_region": len(region_pts),
            "land_points": len(land_pts),
            "land_fraction": round(land_frac, 4),
        }

        # For Pacific: find inhabited islands near the circle
        if target == "Pacific":
            print("\n  Scanning Pacific islands near circle...")
            island_results = []
            for isl in PACIFIC_ISLANDS:
                # Distance from island to great circle
                d = distance_to_gc(isl["lat"], isl["lon"], pole)
                island_results.append({
                    "name": isl["name"],
                    "lat": isl["lat"],
                    "lon": isl["lon"],
                    "population": isl["population"],
                    "distance_km": round(d, 1),
                })
            island_results.sort(key=lambda x: x["distance_km"])

            within_50 = [r for r in island_results if r["distance_km"] <= 50]
            within_100 = [r for r in island_results if r["distance_km"] <= 100]
            within_200 = [r for r in island_results if r["distance_km"] <= 200]
            inhabited_50 = [r for r in within_50 if r["population"] > 0]
            inhabited_100 = [r for r in within_100 if r["population"] > 0]
            inhabited_200 = [r for r in within_200 if r["population"] > 0]

            print(f"  Islands within 50km: {len(within_50)} ({len(inhabited_50)} inhabited)")
            print(f"  Islands within 100km: {len(within_100)} ({len(inhabited_100)} inhabited)")
            print(f"  Islands within 200km: {len(within_200)} ({len(inhabited_200)} inhabited)")
            print(f"\n  Top 20 closest islands:")
            for r in island_results[:20]:
                pop_str = f"pop={r['population']}" if r['population'] > 0 else "uninhabited"
                print(f"    {r['name']}: {r['distance_km']:.0f} km ({pop_str})")

            pair_result["islands_within_50km"] = within_50
            pair_result["islands_within_100km"] = within_100
            pair_result["islands_within_200km"] = within_200
            pair_result["inhabited_within_50km"] = len(inhabited_50)
            pair_result["inhabited_within_100km"] = len(inhabited_100)
            pair_result["inhabited_within_200km"] = len(inhabited_200)
            pair_result["all_islands_ranked"] = island_results

        results[name] = pair_result

    # Verdict
    pacific = results.get("Giza+Nazca→Pacific", {})
    n_inh_100 = pacific.get("inhabited_within_100km", -1)
    if n_inh_100 == 1:
        verdict = "FORCED — Easter Island is the ONLY inhabited island within 100km of the Giza-Nazca great circle in the Pacific"
    elif n_inh_100 == 0:
        verdict = "NOT EVEN EASTER ISLAND is within 100km — check distance calculation"
    else:
        verdict = f"NOT FORCED — {n_inh_100} inhabited islands within 100km; Easter Island is one of several options"

    results["verdict"] = verdict
    print(f"\n  VERDICT: {verdict}")

    out_path = OUT / "constraint_propagation.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Saved: {out_path}")
    return results


# ============================================================
# TEST 2: CONTINENTAL POSITION
# ============================================================
def test_2_continental_position():
    print("\n" + "=" * 70)
    print("TEST 2: CONTINENTAL POSITION")
    print("Among random Africa-Peru-Pacific triplets, how many are as collinear?")
    print("=" * 70)

    # Actual fit
    actual_fit = compute_triplet_fit(GIZA[0], GIZA[1], NAZCA[0], NAZCA[1],
                                     EASTER[0], EASTER[1])
    print(f"  Actual Giza-Nazca-Easter fit: {actual_fit:.2f} km")

    np.random.seed(42)
    n_trials = 100000

    # Pacific islands (inhabited only)
    pacific_inhabited = [i for i in PACIFIC_ISLANDS if i["population"] > 0]

    configs = {
        "narrow_zones": {
            "africa": {"lat_min": 25, "lat_max": 32, "lon_min": 28, "lon_max": 35},
            "peru": {"lat_min": -20, "lat_max": -5, "lon_min": -80, "lon_max": -70},
            "description": "Narrow: Nile Valley + Peru coast + Pacific islands"
        },
        "expanded_zones": {
            "africa": {"lat_min": 25, "lat_max": 35, "lon_min": 25, "lon_max": 50},
            "peru": {"lat_min": -30, "lat_max": -5, "lon_min": -80, "lon_max": -65},
            "description": "Expanded: Egypt/Levant + Peru/Chile coast + Pacific islands"
        },
    }

    results = {"actual_fit_km": round(actual_fit, 2)}

    for config_name, cfg in configs.items():
        print(f"\n  --- {cfg['description']} ---")
        af = cfg["africa"]
        pe = cfg["peru"]

        fits = np.empty(n_trials)
        for i in range(n_trials):
            lat1 = np.random.uniform(af["lat_min"], af["lat_max"])
            lon1 = np.random.uniform(af["lon_min"], af["lon_max"])
            lat2 = np.random.uniform(pe["lat_min"], pe["lat_max"])
            lon2 = np.random.uniform(pe["lon_min"], pe["lon_max"])
            isl = pacific_inhabited[np.random.randint(len(pacific_inhabited))]
            fits[i] = compute_triplet_fit(lat1, lon1, lat2, lon2, isl["lat"], isl["lon"])

        count_better = int(np.sum(fits <= actual_fit))
        p_val = count_better / n_trials
        pct = percentileofscore(fits, actual_fit)

        print(f"  Trials: {n_trials}")
        print(f"  Fits <= {actual_fit:.2f} km: {count_better} ({p_val:.6f})")
        print(f"  Percentile of actual: {pct:.1f}th")
        print(f"  Median fit: {np.median(fits):.1f} km")
        print(f"  Mean fit: {np.mean(fits):.1f} km")
        print(f"  5th/25th/50th/75th/95th: {np.percentile(fits, [5,25,50,75,95])}")

        results[config_name] = {
            "description": cfg["description"],
            "n_trials": n_trials,
            "count_better": count_better,
            "p_value": round(p_val, 8),
            "percentile": round(pct, 2),
            "median_fit_km": round(float(np.median(fits)), 2),
            "mean_fit_km": round(float(np.mean(fits)), 2),
            "percentiles": {str(p): round(float(np.percentile(fits, p)), 2)
                           for p in [1, 5, 10, 25, 50, 75, 90, 95, 99]},
        }

    # Verdict
    p_narrow = results["narrow_zones"]["p_value"]
    p_expanded = results["expanded_zones"]["p_value"]
    if p_narrow > 0.05:
        verdict = f"Continental geometry FORCES collinearity (p={p_narrow:.4f} for narrow zones)"
    elif p_expanded > 0.05:
        verdict = f"Broad continental zones explain it (p_expanded={p_expanded:.4f}) but specific positions add some (p_narrow={p_narrow:.4f})"
    else:
        verdict = f"Specific positions matter beyond continental geometry (p_narrow={p_narrow:.6f}, p_expanded={p_expanded:.6f})"

    results["verdict"] = verdict
    print(f"\n  VERDICT: {verdict}")

    out_path = OUT / "continental_position.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Saved: {out_path}")
    return results


# ============================================================
# TEST 3: RANDOM CIVILIZATION SIMULATION
# ============================================================
def test_3_random_civilization():
    print("\n" + "=" * 70)
    print("TEST 3: RANDOM CIVILIZATION SIMULATION")
    print("Monument-weighted sampling from each region")
    print("=" * 70)

    # Load Pleiades monuments as proxy for "where civilizations build"
    import csv
    pleiades_path = BASE / "data" / "pleiades" / "pleiades-places-latest.csv"

    monument_types = {"temple", "temple-2", "pyramid", "amphitheatre",
                      "sanctuary", "sanctuary-2", "fort", "fortress",
                      "church", "mosque", "tumulus", "tomb", "mausoleum",
                      "monument", "megalithic-monument"}

    africa_monuments = []
    sa_monuments = []

    with open(pleiades_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            lat_str = row.get("reprLat", "")
            lon_str = row.get("reprLong", "")
            if not lat_str or not lon_str:
                continue
            try:
                lat, lon = float(lat_str), float(lon_str)
            except ValueError:
                continue
            ft = set(row.get("featureTypes", "").lower().replace(",", " ").split())
            if not (ft & monument_types):
                continue
            # Africa/ME zone
            if 20 <= lat <= 35 and 25 <= lon <= 50:
                africa_monuments.append((lat, lon))
            # South America zone (Pleiades is weak here, supplement)
            if -25 <= lat <= 0 and -85 <= lon <= -65:
                sa_monuments.append((lat, lon))

    print(f"  Pleiades monuments in Africa/ME zone: {len(africa_monuments)}")
    print(f"  Pleiades monuments in South America zone: {len(sa_monuments)}")

    # Supplement South America with known sites (Pleiades is Old World-focused)
    sa_supplement = [
        (-14.700, -75.100),   # Nazca
        (-10.893, -77.520),   # Caral
        (-13.163, -72.545),   # Machu Picchu
        (-9.593, -77.177),    # Chavin de Huantar
        (-16.554, -68.673),   # Tiwanaku
        (-13.516, -71.978),   # Cusco/Sacsayhuaman
        (-6.780, -79.940),    # Sipan/Huaca Rajada
        (-8.109, -79.074),    # Chan Chan
        (-7.619, -79.424),    # Huaca del Sol
        (-12.015, -76.831),   # Pachacamac
        (-15.500, -75.125),   # Palpa
        (-17.040, -71.604),   # Inca site near Arequipa
        (-15.640, -69.890),   # Sillustani
        (-11.150, -77.350),   # Maranga (Lima)
        (-8.110, -79.050),    # Huanchaco
    ]
    sa_all = sa_monuments + sa_supplement
    print(f"  South America total (with supplement): {len(sa_all)}")

    if len(africa_monuments) < 5:
        print("  WARNING: too few Africa monuments in Pleiades — expanding zone")
        # Fallback: use broader zone
        africa_monuments = [
            (29.979, 31.134),  # Giza
            (29.871, 31.216),  # Saqqara
            (29.790, 31.209),  # Dahshur
            (25.719, 32.657),  # Luxor/Karnak
            (26.185, 31.919),  # Abydos
            (22.337, 31.626),  # Abu Simbel
            (31.244, 29.962),  # Alexandria
            (30.045, 31.233),  # Memphis
            (27.182, 31.168),  # Akhetaten (Amarna)
            (25.740, 32.602),  # Thebes/Valley of Kings
            (24.468, 32.900),  # Kom Ombo
            (24.088, 32.884),  # Aswan
            (33.515, 36.310),  # Baalbek
            (32.220, 35.173),  # Megiddo
            (34.836, 36.282),  # Homs/Palmyra approach
        ]

    actual_fit = compute_triplet_fit(GIZA[0], GIZA[1], NAZCA[0], NAZCA[1],
                                     EASTER[0], EASTER[1])
    pacific_inhabited = [i for i in PACIFIC_ISLANDS if i["population"] > 0]

    np.random.seed(123)
    n_trials = 100000
    fits = np.empty(n_trials)

    for i in range(n_trials):
        a = africa_monuments[np.random.randint(len(africa_monuments))]
        s = sa_all[np.random.randint(len(sa_all))]
        p = pacific_inhabited[np.random.randint(len(pacific_inhabited))]
        fits[i] = compute_triplet_fit(a[0], a[1], s[0], s[1], p["lat"], p["lon"])

    count_better = int(np.sum(fits <= actual_fit))
    p_val = count_better / n_trials
    pct = percentileofscore(fits, actual_fit)

    print(f"\n  Actual fit: {actual_fit:.2f} km")
    print(f"  Trials: {n_trials}")
    print(f"  Better or equal: {count_better} (p={p_val:.6f})")
    print(f"  Percentile: {pct:.1f}th")
    print(f"  Median: {np.median(fits):.1f} km, Mean: {np.mean(fits):.1f} km")

    results = {
        "actual_fit_km": round(actual_fit, 2),
        "n_africa_monuments": len(africa_monuments),
        "n_sa_monuments": len(sa_all),
        "n_pacific_islands": len(pacific_inhabited),
        "n_trials": n_trials,
        "count_better": count_better,
        "p_value": round(p_val, 8),
        "percentile": round(pct, 2),
        "median_fit_km": round(float(np.median(fits)), 2),
        "mean_fit_km": round(float(np.mean(fits)), 2),
        "percentiles": {str(p): round(float(np.percentile(fits, p)), 2)
                       for p in [1, 5, 10, 25, 50, 75, 90, 95, 99]},
    }

    if p_val > 0.05:
        verdict = f"Civilizational geography FORCES collinearity (p={p_val:.4f})"
    elif p_val > 0.01:
        verdict = f"Marginally forced by civilizational geography (p={p_val:.4f})"
    else:
        verdict = f"Specific monument positions matter (p={p_val:.6f})"

    results["verdict"] = verdict
    print(f"  VERDICT: {verdict}")

    out_path = OUT / "random_civilization.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Saved: {out_path}")
    return results


# ============================================================
# TEST 4: HABITABLE ISLAND SCAN
# ============================================================
def test_4_island_scan():
    print("\n" + "=" * 70)
    print("TEST 4: HABITABLE ISLAND SCAN")
    print("Among ALL Pacific islands, how many fall near the Giza-Nazca circle?")
    print("=" * 70)

    pole = pole_from_two_sites(GIZA[0], GIZA[1], NAZCA[0], NAZCA[1])

    results_list = []
    for isl in PACIFIC_ISLANDS:
        d = distance_to_gc(isl["lat"], isl["lon"], pole)
        results_list.append({
            "name": isl["name"],
            "lat": isl["lat"],
            "lon": isl["lon"],
            "population": isl["population"],
            "distance_km": round(d, 1),
        })

    results_list.sort(key=lambda x: x["distance_km"])

    print("\n  Top 25 closest Pacific islands to Giza-Nazca great circle:")
    for r in results_list[:25]:
        pop_str = f"pop={r['population']:,}" if r['population'] > 0 else "uninhabited"
        print(f"    {r['name']}: {r['distance_km']:.0f} km ({pop_str})")

    thresholds = [10, 25, 50, 100, 200, 500]
    summary = {}
    for t in thresholds:
        within = [r for r in results_list if r["distance_km"] <= t]
        inhabited = [r for r in within if r["population"] > 0]
        print(f"\n  Within {t}km: {len(within)} total, {len(inhabited)} inhabited")
        for isl in inhabited:
            print(f"    - {isl['name']} ({isl['distance_km']:.0f} km, pop={isl['population']:,})")
        summary[f"within_{t}km"] = {
            "total": len(within),
            "inhabited": len(inhabited),
            "islands": [r["name"] for r in inhabited]
        }

    # Easter Island distance
    easter_dist = next(r["distance_km"] for r in results_list
                       if "Easter" in r["name"])
    print(f"\n  Easter Island distance: {easter_dist:.1f} km")

    # Verdict
    n_inh_100 = summary["within_100km"]["inhabited"]
    if n_inh_100 <= 1:
        verdict = f"Easter Island is UNIQUE — only {n_inh_100} inhabited island(s) within 100km of the Giza-Nazca circle in the entire Pacific"
    else:
        names = summary["within_100km"]["islands"]
        verdict = f"Easter Island is one of {n_inh_100} inhabited islands within 100km: {', '.join(names)}"

    output = {
        "easter_island_distance_km": easter_dist,
        "all_islands_ranked": results_list,
        "threshold_summary": summary,
        "verdict": verdict,
    }

    print(f"\n  VERDICT: {verdict}")

    out_path = OUT / "pacific_island_scan.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"  Saved: {out_path}")
    return output


# ============================================================
# TEST 5: GEOLOGICAL SOURCE ALIGNMENT
# ============================================================
def test_5_geological():
    print("\n" + "=" * 70)
    print("TEST 5: GEOLOGICAL SOURCE ALIGNMENT")
    print("Are the tectonic features that created these sites also collinear?")
    print("=" * 70)

    geological_sources = [
        {"name": "Gulf of Suez Rift (East African Rift → Nile Valley)", "lat": 29.5, "lon": 32.5},
        {"name": "Peru-Chile Trench (at Nazca latitude)", "lat": -15.0, "lon": -76.0},
        {"name": "Easter Hotspot", "lat": -26.5, "lon": -106.0},
    ]

    # Human sites
    human_fit = compute_triplet_fit(GIZA[0], GIZA[1], NAZCA[0], NAZCA[1],
                                    EASTER[0], EASTER[1])

    # Geological sources
    geo = geological_sources
    geo_fit = compute_triplet_fit(geo[0]["lat"], geo[0]["lon"],
                                  geo[1]["lat"], geo[1]["lon"],
                                  geo[2]["lat"], geo[2]["lon"])

    print(f"  Human site collinearity: {human_fit:.2f} km")
    print(f"  Geological source collinearity: {geo_fit:.2f} km")

    # Compute geological pole
    v1 = latlon_to_cart(geo[0]["lat"], geo[0]["lon"])
    v2 = latlon_to_cart(geo[1]["lat"], geo[1]["lon"])
    v3 = latlon_to_cart(geo[2]["lat"], geo[2]["lon"])
    M = np.vstack([v1, v2, v3])
    M_c = M - M.mean(axis=0)
    _, s, Vt = np.linalg.svd(M_c)
    geo_pole = Vt[-1]
    geo_pole = geo_pole / np.linalg.norm(geo_pole)
    geo_pole_latlon = cart_to_latlon(geo_pole)

    # Alison pole
    alison_pole = latlon_to_cart(POLE_LAT, POLE_LON)
    alison_pole_latlon = (POLE_LAT, POLE_LON)

    # Angular distance between poles
    dot = np.clip(np.dot(geo_pole, alison_pole), -1, 1)
    ang_dist_deg = np.degrees(np.arccos(abs(dot)))  # abs because poles can flip
    print(f"\n  Geological pole: ({geo_pole_latlon[0]:.2f}, {geo_pole_latlon[1]:.2f})")
    print(f"  Alison pole: ({alison_pole_latlon[0]:.2f}, {alison_pole_latlon[1]:.2f})")
    print(f"  Angular distance between poles: {ang_dist_deg:.1f}°")

    # Compare each human site's distance from geological circle vs human circle
    print("\n  Distance from each site to geological vs human great circle:")
    for name, (lat, lon) in [("Giza", GIZA), ("Nazca", NAZCA), ("Easter", EASTER)]:
        d_geo = distance_to_gc(lat, lon, geo_pole)
        d_human_pole = pole_from_two_sites(GIZA[0], GIZA[1], NAZCA[0], NAZCA[1])
        # Use the SVD-based pole for consistency
        v = latlon_to_cart(GIZA[0], GIZA[1])
        v2_ = latlon_to_cart(NAZCA[0], NAZCA[1])
        v3_ = latlon_to_cart(EASTER[0], EASTER[1])
        Mh = np.vstack([v, v2_, v3_])
        Mh_c = Mh - Mh.mean(axis=0)
        _, _, Vth = np.linalg.svd(Mh_c)
        human_pole = Vth[-1]
        human_pole = human_pole / np.linalg.norm(human_pole)
        d_human = distance_to_gc(lat, lon, human_pole)
        print(f"    {name}: geological={d_geo:.1f} km, human={d_human:.1f} km")

    # Random baseline: how often are three random tectonic features this collinear?
    # Generate random "rift/trench/hotspot" triplets
    # Rifts: roughly 15-45°N/S, various longitudes
    # Trenches: along subduction zones
    # Hotspots: scattered across ocean floors
    # Simplified: just random triplets spanning similar geographic extent
    np.random.seed(77)
    n_trials = 100000
    random_fits = np.empty(n_trials)
    for i in range(n_trials):
        # Random Africa/ME point
        lat1 = np.random.uniform(20, 35)
        lon1 = np.random.uniform(25, 50)
        # Random Pacific trench point (ring of fire west coast SA)
        lat2 = np.random.uniform(-25, -5)
        lon2 = np.random.uniform(-82, -70)
        # Random Pacific hotspot
        lat3 = np.random.uniform(-35, -15)
        lon3 = np.random.uniform(-130, -90)
        random_fits[i] = compute_triplet_fit(lat1, lon1, lat2, lon2, lat3, lon3)

    geo_pct = percentileofscore(random_fits, geo_fit)
    print(f"\n  Geological fit percentile vs random: {geo_pct:.1f}th")
    print(f"  Median random fit: {np.median(random_fits):.1f} km")

    if ang_dist_deg < 5:
        verdict = "Tectonic sources are ALIGNED with human sites — plate tectonics is the ultimate cause"
    elif ang_dist_deg < 15:
        verdict = f"Tectonic sources are MODERATELY aligned ({ang_dist_deg:.1f}° apart) — geology partially explains it"
    else:
        verdict = f"Tectonic sources are NOT aligned ({ang_dist_deg:.1f}° apart) — human positions aren't inherited from geology"

    results = {
        "human_fit_km": round(human_fit, 2),
        "geological_fit_km": round(geo_fit, 2),
        "geological_sources": geological_sources,
        "geological_pole": {"lat": round(geo_pole_latlon[0], 2), "lon": round(geo_pole_latlon[1], 2)},
        "alison_pole": {"lat": POLE_LAT, "lon": POLE_LON},
        "pole_angular_distance_deg": round(ang_dist_deg, 2),
        "geological_fit_percentile_vs_random": round(geo_pct, 2),
        "verdict": verdict,
    }

    print(f"\n  VERDICT: {verdict}")

    out_path = OUT / "geological_alignment.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Saved: {out_path}")
    return results


# ============================================================
# TEST 6: ATMOSPHERIC CIRCULATION
# ============================================================
def test_6_atmospheric():
    print("\n" + "=" * 70)
    print("TEST 6: ATMOSPHERIC CIRCULATION / HADLEY CELL")
    print("Does the Great Circle track the 30th parallel (desert belt)?")
    print("=" * 70)

    # Trace the Alison circle
    alison_pole = latlon_to_cart(POLE_LAT, POLE_LON)
    circle = trace_great_circle(alison_pole, step_deg=0.5)

    # Get elevation for all points
    c_lats = [p[0] for p in circle]
    c_lons = [p[1] for p in circle]
    elevs = get_elev_batch(c_lats, c_lons)

    land_mask = elevs > 0
    land_lats = np.array(c_lats)[land_mask]
    land_lons = np.array(c_lons)[land_mask]

    # Fraction of land points near 30th parallel (within 5°)
    near_30 = np.sum(np.abs(np.abs(land_lats) - 30) <= 5)
    total_land = len(land_lats)
    frac_near_30 = near_30 / max(total_land, 1)

    print(f"  Total circle points: {len(circle)}")
    print(f"  Land points: {total_land}")
    print(f"  Land points within 5° of 30th parallel (N or S): {near_30}")
    print(f"  Fraction near 30°: {frac_near_30:.1%}")

    # Compare to random circles
    print("\n  Computing random circle baseline (1000 random poles)...")
    np.random.seed(99)
    random_fractions = []
    for _ in range(1000):
        rp = random_pole()
        rc = trace_great_circle(rp, step_deg=0.5)
        rc_lats = [p[0] for p in rc]
        rc_lons = [p[1] for p in rc]
        rc_elevs = get_elev_batch(rc_lats, rc_lons)
        rc_land = np.array(rc_lats)[rc_elevs > 0]
        if len(rc_land) > 0:
            rc_near = np.sum(np.abs(np.abs(rc_land) - 30) <= 5)
            random_fractions.append(rc_near / len(rc_land))

    pct = percentileofscore(random_fractions, frac_near_30)
    print(f"  Percentile vs random: {pct:.0f}th")
    print(f"  Mean random fraction: {np.mean(random_fractions):.3f}")
    print(f"  Alison fraction: {frac_near_30:.3f}")

    # Also check: latitude distribution of land points on the circle
    lat_bins = np.arange(-90, 91, 10)
    hist, _ = np.histogram(land_lats, bins=lat_bins)
    print(f"\n  Latitude distribution of land points on Alison circle:")
    for i in range(len(hist)):
        if hist[i] > 0:
            print(f"    {lat_bins[i]:+4.0f}° to {lat_bins[i+1]:+4.0f}°: {hist[i]} points")

    # Fraction in arid latitudes (15-35° N or S)
    arid_mask = ((np.abs(land_lats) >= 15) & (np.abs(land_lats) <= 35))
    frac_arid = np.sum(arid_mask) / max(total_land, 1)
    print(f"\n  Fraction in arid belt (15-35° N/S): {frac_arid:.1%}")

    if pct > 75:
        verdict = f"Circle TRACKS the desert belt ({pct:.0f}th percentile for 30° proximity)"
    elif pct > 50:
        verdict = f"Circle is somewhat aligned with desert belt ({pct:.0f}th percentile)"
    else:
        verdict = f"Circle does NOT especially track the desert belt ({pct:.0f}th percentile)"

    results = {
        "total_circle_points": len(circle),
        "land_points": total_land,
        "near_30th_parallel_count": int(near_30),
        "fraction_near_30": round(frac_near_30, 4),
        "percentile_vs_random": round(pct, 1),
        "mean_random_fraction": round(float(np.mean(random_fractions)), 4),
        "fraction_arid_belt": round(frac_arid, 4),
        "latitude_distribution": {f"{lat_bins[i]:+.0f}_to_{lat_bins[i+1]:+.0f}": int(hist[i])
                                  for i in range(len(hist)) if hist[i] > 0},
        "verdict": verdict,
    }

    print(f"\n  VERDICT: {verdict}")

    out_path = OUT / "atmospheric_circulation.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Saved: {out_path}")
    return results


# ============================================================
# TEST 7: EQUIVALENT SITE SUBSTITUTION
# ============================================================
def test_7_equivalent_substitution():
    print("\n" + "=" * 70)
    print("TEST 7: EQUIVALENT SITE SUBSTITUTION")
    print("Replace each anchor with a culturally equivalent site")
    print("=" * 70)

    equivalents = {
        "Giza": [
            {"name": "Giza", "lat": 29.979, "lon": 31.134},
            {"name": "Saqqara", "lat": 29.871, "lon": 31.216},
            {"name": "Dahshur", "lat": 29.790, "lon": 31.209},
            {"name": "Luxor/Karnak", "lat": 25.719, "lon": 32.657},
            {"name": "Abydos", "lat": 26.185, "lon": 31.919},
            {"name": "Abu Simbel", "lat": 22.337, "lon": 31.626},
            {"name": "Baalbek", "lat": 34.007, "lon": 36.204},
            {"name": "Petra", "lat": 30.329, "lon": 35.444},
        ],
        "Nazca": [
            {"name": "Nazca", "lat": -14.700, "lon": -75.100},
            {"name": "Palpa geoglyphs", "lat": -14.517, "lon": -75.183},
            {"name": "Caral", "lat": -10.893, "lon": -77.520},
            {"name": "Chavin de Huantar", "lat": -9.593, "lon": -77.177},
            {"name": "Tiwanaku", "lat": -16.554, "lon": -68.673},
            {"name": "Machu Picchu", "lat": -13.163, "lon": -72.545},
            {"name": "Cusco/Sacsayhuaman", "lat": -13.516, "lon": -71.978},
            {"name": "Chan Chan", "lat": -8.109, "lon": -79.074},
        ],
        "Easter Island": [
            {"name": "Easter Island", "lat": -27.100, "lon": -109.333},
            {"name": "Pitcairn", "lat": -25.07, "lon": -130.10},
            {"name": "Mangareva", "lat": -23.12, "lon": -134.97},
            {"name": "Marquesas (Nuku Hiva)", "lat": -8.86, "lon": -140.10},
            {"name": "Rarotonga", "lat": -21.23, "lon": -159.78},
            {"name": "Galapagos", "lat": -0.95, "lon": -90.97},
            {"name": "Tahiti", "lat": -17.68, "lon": -149.45},
            {"name": "Samoa (Upolu)", "lat": -13.83, "lon": -171.76},
            {"name": "Tongatapu", "lat": -21.21, "lon": -175.15},
            {"name": "Fiji (Viti Levu)", "lat": -17.77, "lon": 178.07},
        ],
    }

    n_giza = len(equivalents["Giza"])
    n_nazca = len(equivalents["Nazca"])
    n_easter = len(equivalents["Easter Island"])
    total = n_giza * n_nazca * n_easter
    print(f"  Giza equivalents: {n_giza}")
    print(f"  Nazca equivalents: {n_nazca}")
    print(f"  Easter equivalents: {n_easter}")
    print(f"  Total combinations: {total}")

    all_results = []
    for eg in equivalents["Giza"]:
        for na in equivalents["Nazca"]:
            for ei in equivalents["Easter Island"]:
                fit = compute_triplet_fit(eg["lat"], eg["lon"],
                                         na["lat"], na["lon"],
                                         ei["lat"], ei["lon"])
                all_results.append({
                    "giza_equiv": eg["name"],
                    "nazca_equiv": na["name"],
                    "easter_equiv": ei["name"],
                    "fit_km": round(fit, 2),
                })

    all_results.sort(key=lambda x: x["fit_km"])

    # Find the actual triplet
    actual = next(r for r in all_results
                  if r["giza_equiv"] == "Giza"
                  and r["nazca_equiv"] == "Nazca"
                  and r["easter_equiv"] == "Easter Island")
    actual_rank = next(i for i, r in enumerate(all_results)
                       if r["giza_equiv"] == "Giza"
                       and r["nazca_equiv"] == "Nazca"
                       and r["easter_equiv"] == "Easter Island") + 1

    print(f"\n  Actual Giza-Nazca-Easter: {actual['fit_km']:.2f} km (rank {actual_rank}/{total})")
    print(f"\n  Top 20 combinations:")
    for i, r in enumerate(all_results[:20]):
        marker = " <<<" if (r["giza_equiv"] == "Giza" and r["nazca_equiv"] == "Nazca"
                           and r["easter_equiv"] == "Easter Island") else ""
        print(f"    {i+1}. {r['giza_equiv']} + {r['nazca_equiv']} + {r['easter_equiv']}: "
              f"{r['fit_km']:.2f} km{marker}")

    # Threshold analysis
    thresholds = [1, 5, 10, 50, 100, 500]
    threshold_counts = {}
    for t in thresholds:
        count = sum(1 for r in all_results if r["fit_km"] <= t)
        threshold_counts[f"within_{t}km"] = count
        print(f"  Combinations with fit < {t} km: {count}/{total}")

    # How many combinations beat the actual?
    better = sum(1 for r in all_results if r["fit_km"] < actual["fit_km"])
    print(f"\n  Combinations with BETTER fit than actual: {better}/{total}")

    if actual_rank == 1:
        verdict = "The actual Giza-Nazca-Easter triplet is the BEST combination — specific sites matter"
    elif actual_rank <= 5:
        verdict = f"The actual triplet ranks {actual_rank}/{total} — near the top, specific positions are important"
    elif better / total > 0.5:
        verdict = f"The actual triplet ranks {actual_rank}/{total} — many equivalent combos are better, regions matter more than specific sites"
    else:
        verdict = f"The actual triplet ranks {actual_rank}/{total} — moderately special among equivalents"

    results = {
        "n_combinations": total,
        "actual_triplet": actual,
        "actual_rank": actual_rank,
        "better_combinations": better,
        "top_20": all_results[:20],
        "threshold_counts": threshold_counts,
        "all_results": all_results,
        "verdict": verdict,
    }

    print(f"\n  VERDICT: {verdict}")

    out_path = OUT / "equivalent_substitution.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Saved: {out_path}")
    return results


# ============================================================
# SYNTHESIS
# ============================================================
def synthesize(all_results):
    print("\n" + "=" * 70)
    print("SYNTHESIS: COLLINEARITY RESOLUTION")
    print("=" * 70)

    verdicts = {}
    for test_name, result in all_results.items():
        if isinstance(result, dict) and "verdict" in result:
            verdicts[test_name] = result["verdict"]
            print(f"\n  {test_name}: {result['verdict']}")

    # Overall assessment
    forced_count = sum(1 for v in verdicts.values()
                       if "FORCE" in v.upper() or "ONLY" in v.upper())
    special_count = sum(1 for v in verdicts.values()
                        if "SPECIFIC" in v.upper() or "MATTER" in v.upper()
                        or "SPECIAL" in v.upper() or "BEST" in v.upper())

    if forced_count > special_count:
        overall = ("RESOLVED: The Giza-Nazca-Easter Island collinearity is largely explained "
                   "by continental geometry and where habitable land exists. The circle is "
                   "famous because of which sites it connects, not because the connection "
                   "is geometrically remarkable.")
    elif special_count > forced_count:
        overall = ("UNRESOLVED: The specific positions of Giza, Nazca, and Easter Island "
                   "produce unusually good collinearity even among geographic/cultural "
                   "equivalents. The collinearity appears to be a genuine geometric coincidence.")
    else:
        overall = ("MIXED: Some aspects of the collinearity are forced by geography, "
                   "but the specific positions contribute. The truth lies between "
                   "'pure coincidence' and 'forced by continental geometry.'")

    print(f"\n  OVERALL: {overall}")

    synthesis = {
        "verdicts": verdicts,
        "forced_indicators": forced_count,
        "special_indicators": special_count,
        "overall_assessment": overall,
    }

    out_path = OUT / "synthesis.json"
    with open(out_path, "w") as f:
        json.dump(synthesis, f, indent=2)
    print(f"  Saved: {out_path}")
    return synthesis


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    t0 = time.time()
    all_results = {}

    print(f"\nStarting collinearity resolution tests at {time.strftime('%H:%M:%S')}")

    all_results["test_1_constraint_propagation"] = test_1_constraint_propagation()
    all_results["test_2_continental_position"] = test_2_continental_position()
    all_results["test_3_random_civilization"] = test_3_random_civilization()
    all_results["test_4_island_scan"] = test_4_island_scan()
    all_results["test_5_geological"] = test_5_geological()
    all_results["test_6_atmospheric"] = test_6_atmospheric()
    all_results["test_7_equivalent_substitution"] = test_7_equivalent_substitution()

    synthesize(all_results)

    elapsed = time.time() - t0
    print(f"\n{'=' * 70}")
    print(f"ALL TESTS COMPLETE in {elapsed/60:.1f} minutes")
    print(f"{'=' * 70}")
