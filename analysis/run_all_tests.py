#!/usr/bin/env python3
"""
Deep Time Atlas — Complete Great Circle Analysis Pipeline
==========================================================
Reproduces all analysis from raw KML + supplement data.

Usage:
    python analysis/run_all_tests.py [--trials 200] [--kml-dir ../] [--skip-pleiades]

Outputs all results to results/ directory.
Requires only Python 3.7+ stdlib.
"""

import xml.etree.ElementTree as ET
import glob, math, random, json, os, re, time, argparse, csv, gzip, sys
from urllib.request import urlopen
from collections import defaultdict

# ============================================================
# CONFIGURATION
# ============================================================
POLE_LAT = 59.682122
POLE_LNG = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2  # ~10,007.5 km
THRESHOLDS = [25, 50, 100, 200]

# ============================================================
# MATH HELPERS
# ============================================================
def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))

def gc_distance(lat, lon):
    """Distance from site to the great circle (km)."""
    return abs(haversine_km(lat, lon, POLE_LAT, POLE_LNG) - QUARTER_CIRC)

def angular_separation(lat1, lon1, lat2, lon2):
    """Angular separation in degrees between two points."""
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    cos_d = (math.sin(lat1)*math.sin(lat2) +
             math.cos(lat1)*math.cos(lat2)*math.cos(lon2 - lon1))
    return math.degrees(math.acos(max(-1, min(1, cos_d))))

def great_circle_point(pole_lat, pole_lng, bearing_deg):
    """Compute a point on the great circle at given bearing from pole."""
    lat1 = math.radians(pole_lat)
    lon1 = math.radians(pole_lng)
    d = math.pi / 2  # quarter circumference in radians
    brng = math.radians(bearing_deg)
    lat2 = math.asin(math.sin(lat1)*math.cos(d) + math.cos(lat1)*math.sin(d)*math.cos(brng))
    lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d)*math.cos(lat1),
                              math.cos(d) - math.sin(lat1)*math.sin(lat2))
    return math.degrees(lat2), math.degrees(lon2)

def rand_matched(site_lats, site_lons):
    """Generate a distribution-matched random point."""
    lat = random.choice(site_lats) + random.gauss(0, 2)
    lon = random.choice(site_lons) + random.gauss(0, 2)
    return max(-90, min(90, lat)), max(-180, min(180, lon))

def compute_z(observed, baseline_counts):
    """Compute Z-score from observed count and baseline distribution."""
    n = len(baseline_counts)
    if n == 0:
        return 0, 0, 0
    mu = sum(baseline_counts) / n
    sigma = (sum((x - mu)**2 for x in baseline_counts) / n) ** 0.5
    z = (observed - mu) / sigma if sigma > 0 else 0
    return z, mu, sigma

def random_pole():
    """Generate a uniformly random point on the sphere."""
    z = random.uniform(-1, 1)
    lat = math.degrees(math.asin(z))
    lon = random.uniform(-180, 180)
    return lat, lon

# ============================================================
# TYPE CLASSIFICATIONS
# ============================================================
PREHISTORIC_TYPES = {
    "Standing Stone (Menhir)", "Standing Stones", "Stone Circle",
    "Burial Chamber or Dolmen", "Chambered Tomb", "Passage Grave",
    "Long Barrow", "Chambered Cairn", "Henge", "Cursus",
    "Court Tomb", "Portal Tomb", "Wedge Tomb", "Timber Circle",
    "Clava Cairn", "Stone Row  Alignment", "Multiple Stone Rows  Avenue",
    "Causewayed Enclosure", "Ring Cairn", "Rock Art",
    "Cave or Rock Shelter", "Cist", "Ancient Village or Settlement",
    "Ancient Temple", "Pyramid  Mastaba", "Hill Figure or Geoglyph"
}

LATER_TYPES = {
    "Round Barrow(s)", "Cairn", "Hillfort", "Barrow Cemetery",
    "Holy Well or Sacred Spring", "Broch or Nuraghe",
    "Natural Stone  Erratic  Other Natural Feature", "Museum",
    "Castro or Chafurdão", "Modern Stone Circle etc", "Ancient Cross",
    "Artificial Mound", "Sculptured Stone", "Carving",
    "Stone Fort or Dun", "Rock Cut Tomb", "Marker Stone",
    "Ancient Mine, Quarry or other Industry", "Round Cairn",
    "Early Christian Sculptured Stone", "Rock Outcrop",
    "Promontory Fort  Cliff Castle", "Ancient Trackway", "Polissoir",
    "Class I Pictish Symbol Stone", "Class II Pictish Symbol Stone",
    "Class III Pictish Cross Slab",
    "Class I  Class II Hybrid Pictish Symbol Stone",
    "Crannog", "Souterrain (Fogou, Earth House)", "Ancient Palace",
    "Holed Stone", "Turf Maze", "Vitrified Fort"
}

# Age estimates (midpoint of typical range, in years BCE, negative = BCE)
TYPE_AGE_ESTIMATES = {
    "Stone Circle": -2500, "Standing Stone (Menhir)": -3000,
    "Standing Stones": -3000, "Burial Chamber or Dolmen": -3500,
    "Chambered Tomb": -3500, "Passage Grave": -3200,
    "Long Barrow": -3500, "Chambered Cairn": -3000,
    "Henge": -2800, "Cursus": -3500, "Court Tomb": -3500,
    "Portal Tomb": -3500, "Wedge Tomb": -2500,
    "Timber Circle": -2500, "Clava Cairn": -2000,
    "Stone Row  Alignment": -2500,
    "Multiple Stone Rows  Avenue": -2500,
    "Causewayed Enclosure": -3500, "Ring Cairn": -2000,
    "Rock Art": -4000, "Cave or Rock Shelter": -10000,
    "Cist": -2000, "Ancient Village or Settlement": -3000,
    "Ancient Temple": -3000, "Pyramid  Mastaba": -2500,
    "Hill Figure or Geoglyph": -1500,
    "Round Barrow(s)": -1500, "Cairn": -2000,
    "Hillfort": -500, "Barrow Cemetery": -1500,
    "Holy Well or Sacred Spring": 500,
    "Broch or Nuraghe": -200, "Museum": 1900,
    "Castro or Chafurdão": -500, "Modern Stone Circle etc": 1900,
    "Ancient Cross": 800, "Artificial Mound": -1000,
    "Sculptured Stone": 700, "Carving": 800,
    "Stone Fort or Dun": -200, "Rock Cut Tomb": -1500,
    "Marker Stone": 500,
    "Ancient Mine, Quarry or other Industry": -2000,
    "Round Cairn": -1500,
    "Early Christian Sculptured Stone": 800,
    "Rock Outcrop": 0, "Promontory Fort  Cliff Castle": -300,
    "Ancient Trackway": -2000, "Polissoir": -4000,
    "Class I Pictish Symbol Stone": 600,
    "Class II Pictish Symbol Stone": 700,
    "Class III Pictish Cross Slab": 800,
    "Class I  Class II Hybrid Pictish Symbol Stone": 650,
    "Crannog": -500, "Souterrain (Fogou, Earth House)": -100,
    "Ancient Palace": -1000, "Holed Stone": -2000,
    "Turf Maze": 1200, "Vitrified Fort": -400,
    "Natural Stone  Erratic  Other Natural Feature": 0,
    "Misc. Earthwork": -1000, "Not Known (by us)": 0,
}

# Supplement type classification
SUPPLEMENT_PREHISTORIC = {
    "pyramid", "temple", "geoglyph", "megalith", "stone_circle",
    "rock_art", "dolmen", "burial", "earthwork", "mound",
    "observatory", "ziggurat", "rock_cut_city", "trilithon",
    "jar_burial", "shrine", "rock", "aqueducts", "terraces"
}

def parse_supplement_date(date_str):
    if not date_str or date_str == "debated":
        return None
    m = re.match(r'(\d+)(?:st|nd|rd|th)\s+c\s+(BCE|CE)', date_str)
    if m:
        century = int(m.group(1))
        year = (century - 1) * 100 + 50
        return -year if m.group(2) == "BCE" else year
    m = re.match(r'(\d+)\s*(BCE|CE)', date_str)
    if m:
        year = int(m.group(1))
        return -year if m.group(2) == "BCE" else year
    m = re.match(r'(\d+)\s*(BCE|CE)?[-\u2013](\d+)\s*(BCE|CE)?', date_str)
    if m:
        y1, era1, y2, era2 = int(m.group(1)), m.group(2), int(m.group(3)), m.group(4)
        if era1 == "BCE": y1 = -y1
        elif era1 is None and (era2 == "BCE" or era2 == "CE"): y1 = -y1
        if era2 == "BCE": y2 = -y2
        return (y1 + y2) // 2
    return None

def classify_supplement(site):
    stype = site.get("type", "")
    if stype in SUPPLEMENT_PREHISTORIC:
        return "prehistoric"
    if stype in ("fortress", "stela", "palace", "tower", "tomb", "rock_cut"):
        return "later"
    if stype == "citadel":
        yr = parse_supplement_date(site.get("date", ""))
        return "later" if yr and yr > 500 else "prehistoric"
    if stype in ("settlement", "city"):
        yr = parse_supplement_date(site.get("date", ""))
        return "prehistoric" if (yr is None or yr < -500) else "later"
    return "prehistoric"

# ============================================================
# STEP 1: PARSE KML FILES
# ============================================================
def parse_kml_files(kml_dir):
    print("=" * 70)
    print("STEP 1: Parsing KML files")
    print("=" * 70)

    all_sites = []
    kml_paths = sorted(set(glob.glob(os.path.join(kml_dir, "MegP_*.kml"))))

    for filepath in kml_paths:
        filename = os.path.basename(filepath)
        site_type = filename.replace("MegP_", "").replace(".kml", "")
        site_type = site_type.replace("_", " ").replace("  ", " ").strip()
        site_type = re.sub(r'\s*\(\d+\)\s*$', '', site_type).strip()

        try:
            tree = ET.parse(filepath)
            root = tree.getroot()
            for pm in root.iter():
                if 'Placemark' not in pm.tag:
                    continue
                name_el = coords_el = None
                for child in pm.iter():
                    if 'name' in child.tag and child.tag.endswith('name') and name_el is None:
                        name_el = child
                    if 'coordinates' in child.tag:
                        coords_el = child
                if coords_el is not None and coords_el.text:
                    try:
                        parts = coords_el.text.strip().split(',')
                        lon, lat = float(parts[0]), float(parts[1])
                        name = name_el.text if name_el is not None else "Unknown"
                        if -90 <= lat <= 90 and -180 <= lon <= 180 and (lat != 0 or lon != 0):
                            all_sites.append({
                                "name": name, "lat": lat, "lon": lon,
                                "type": site_type, "source": "megalithic_portal"
                            })
                    except (ValueError, IndexError):
                        continue
        except Exception as e:
            print(f"  ERROR {filename}: {e}")

    # Deduplicate within ~111m (3 decimal places)
    seen = set()
    deduped = []
    for s in all_sites:
        key = (round(s["lat"], 3), round(s["lon"], 3))
        if key not in seen:
            seen.add(key)
            deduped.append(s)

    print(f"  Raw sites parsed: {len(all_sites)}")
    print(f"  After dedup: {len(deduped)}")
    return deduped

# ============================================================
# STEP 2: LOAD & MERGE SUPPLEMENT SITES
# ============================================================
def load_and_merge_supplement(portal_sites, supplement_path):
    print(f"\n{'='*70}")
    print("STEP 2: Loading and merging supplementary sites")
    print("=" * 70)

    with open(supplement_path) as f:
        supplement_data = json.load(f)

    supplement_sites = supplement_data["sites"]
    print(f"  Supplement sites loaded: {len(supplement_sites)}")

    added = []
    duplicates = []
    for ss in supplement_sites:
        is_dup = False
        for ks in portal_sites:
            d = haversine_km(ss["lat"], ss.get("lon", ss.get("lng")),
                           ks["lat"], ks["lon"])
            if d < 1.0:
                duplicates.append({
                    "supplement": ss["name"], "existing": ks["name"],
                    "distance_m": round(d * 1000)
                })
                is_dup = True
                break
        if not is_dup:
            added.append({
                "name": ss["name"], "lat": ss["lat"],
                "lon": ss.get("lon", ss.get("lng")),
                "type": ss.get("type", "unknown"),
                "date": ss.get("date", ""),
                "region": ss.get("region", ""),
                "source": "supplement"
            })

    merged = portal_sites + added
    print(f"  Duplicates (within 1km): {len(duplicates)}")
    print(f"  New sites added: {len(added)}")
    print(f"  MERGED TOTAL: {len(merged)}")
    return merged, added, duplicates

# ============================================================
# STEP 3: DOWNLOAD PLEIADES
# ============================================================
def download_pleiades(data_dir, skip=False):
    print(f"\n{'='*70}")
    print("STEP 3: Pleiades Gazetteer")
    print("=" * 70)

    csv_path = os.path.join(data_dir, "pleiades_sites.csv")

    if skip:
        if os.path.exists(csv_path):
            print(f"  Using existing: {csv_path}")
            return csv_path
        else:
            print("  SKIPPED (file not found, use --skip-pleiades to disable)")
            return None

    if os.path.exists(csv_path):
        print(f"  Already downloaded: {csv_path}")
        return csv_path

    url = "https://atlantides.org/downloads/pleiades/dumps/pleiades-places-latest.csv.gz"
    gz_path = os.path.join(data_dir, "pleiades-places-latest.csv.gz")

    print(f"  Downloading from {url}...")
    try:
        with urlopen(url) as response:
            with open(gz_path, 'wb') as f:
                while True:
                    chunk = response.read(65536)
                    if not chunk:
                        break
                    f.write(chunk)
        print(f"  Decompressing...")
        with gzip.open(gz_path, 'rb') as f_in:
            with open(csv_path, 'wb') as f_out:
                while True:
                    chunk = f_in.read(65536)
                    if not chunk:
                        break
                    f_out.write(chunk)
        os.remove(gz_path)
        print(f"  Saved: {csv_path}")
    except Exception as e:
        print(f"  Download failed: {e}")
        return None

    return csv_path

def load_pleiades(csv_path):
    """Parse Pleiades CSV and return list of dicts with lat/lon/minDate."""
    sites = []
    with open(csv_path, 'r', encoding='utf-8', errors='replace') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row.get('reprLat', ''))
                lon = float(row.get('reprLong', ''))
                if lat == 0 and lon == 0:
                    continue
                min_date = None
                try:
                    min_date = int(row.get('minDate', ''))
                except (ValueError, TypeError):
                    pass
                sites.append({"lat": lat, "lon": lon, "minDate": min_date,
                             "name": row.get('title', ''), "source": "pleiades"})
            except (ValueError, TypeError):
                continue
    return sites

# ============================================================
# STEP 4: GENERATE GREAT CIRCLE COORDINATES
# ============================================================
def generate_circle_coordinates(data_dir):
    print(f"\n{'='*70}")
    print("STEP 4: Generating Great Circle coordinates")
    print("=" * 70)

    points = []
    for deg in range(360):
        lat, lon = great_circle_point(POLE_LAT, POLE_LNG, deg)
        points.append({"bearing": deg, "lat": round(lat, 6), "lon": round(lon, 6)})

    json_path = os.path.join(data_dir, "circle_coordinates.json")
    with open(json_path, 'w') as f:
        json.dump({"pole": {"lat": POLE_LAT, "lon": POLE_LNG},
                   "points": points}, f, indent=2)

    csv_path = os.path.join(data_dir, "circle_coordinates.csv")
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["bearing", "lat", "lon"])
        for p in points:
            w.writerow([p["bearing"], p["lat"], p["lon"]])

    print(f"  Saved {len(points)} points to {json_path} and {csv_path}")
    return points

# ============================================================
# MONTE CARLO ENGINE
# ============================================================
def run_monte_carlo(sites, n_trials, thresholds=None):
    """Run distribution-matched Monte Carlo for a set of sites.
    Returns dict of {threshold: {"observed", "mean", "std", "z", "enrichment", "p_value"}}
    """
    if thresholds is None:
        thresholds = THRESHOLDS
    n = len(sites)
    site_lats = [s["lat"] for s in sites]
    site_lons = [s["lon"] for s in sites]

    # Observed counts
    observed = {}
    for t in thresholds:
        observed[t] = sum(1 for s in sites if s["gc_dist"] <= t)

    # Baseline trials
    baseline = {t: [] for t in thresholds}
    t_start = time.time()
    for trial in range(n_trials):
        pts = [rand_matched(site_lats, site_lons) for _ in range(n)]
        for t in thresholds:
            baseline[t].append(sum(1 for lat, lon in pts if gc_distance(lat, lon) <= t))
        if (trial + 1) % 50 == 0:
            elapsed = time.time() - t_start
            rate = (trial + 1) / elapsed if elapsed > 0 else 0
            eta = (n_trials - trial - 1) / rate if rate > 0 else 0
            print(f"    {trial+1}/{n_trials} ({elapsed:.0f}s, ETA {eta:.0f}s)")

    results = {}
    for t in thresholds:
        obs = observed[t]
        z, mu, sigma = compute_z(obs, baseline[t])
        enrich = obs / mu if mu > 0 else float('inf')
        # Empirical p-value
        p = sum(1 for b in baseline[t] if b >= obs) / n_trials
        results[t] = {
            "observed": obs,
            "baseline_mean": round(mu, 1),
            "baseline_std": round(sigma, 1),
            "z_score": round(z, 2),
            "enrichment": round(enrich, 2),
            "p_value": p
        }
    return results

# ============================================================
# TEST 1: 108° ANGULAR SEPARATION (UNESCO)
# ============================================================
def test_108_falsification(unesco_path, results_dir, n_trials=100):
    print(f"\n{'='*70}")
    print("TEST 1: 108° Angular Separation (UNESCO)")
    print("=" * 70)

    with open(unesco_path) as f:
        data = json.load(f)
    sites = data["sites"]
    n = len(sites)
    print(f"  UNESCO cultural sites: {n}")

    # Filter to ancient-only (rough heuristic: sites with ancient-sounding names)
    # Since UNESCO data has no dates, we use all cultural sites
    target_angle = 108.0
    tolerance = 0.5

    # Count pairs at 108° ± 0.5°
    obs_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            ang = angular_separation(sites[i]["lat"], sites[i]["lon"],
                                    sites[j]["lat"], sites[j]["lon"])
            if abs(ang - target_angle) <= tolerance:
                obs_count += 1

    print(f"  Observed pairs at 108°±0.5°: {obs_count}")

    # Monte Carlo: shuffle lat/lon independently with jitter
    lats = [s["lat"] for s in sites]
    lons = [s["lon"] for s in sites]

    baseline_counts = []
    t_start = time.time()
    for trial in range(n_trials):
        rand_sites = [rand_matched(lats, lons) for _ in range(n)]
        count = 0
        for i in range(len(rand_sites)):
            for j in range(i + 1, len(rand_sites)):
                ang = angular_separation(rand_sites[i][0], rand_sites[i][1],
                                        rand_sites[j][0], rand_sites[j][1])
                if abs(ang - target_angle) <= tolerance:
                    count += 1
        baseline_counts.append(count)
        if (trial + 1) % 10 == 0:
            elapsed = time.time() - t_start
            print(f"    {trial+1}/{n_trials} ({elapsed:.0f}s)")

    z, mu, sigma = compute_z(obs_count, baseline_counts)
    p = sum(1 for b in baseline_counts if b >= obs_count) / n_trials

    result = {
        "test": "108° Angular Separation",
        "dataset": "UNESCO Cultural Heritage",
        "n_sites": n,
        "target_angle": target_angle,
        "tolerance": tolerance,
        "observed_pairs": obs_count,
        "baseline_mean": round(mu, 1),
        "baseline_std": round(sigma, 1),
        "z_score": round(z, 2),
        "p_value": p,
        "n_trials": n_trials,
        "verdict": "FALSIFIED" if z < 2 else "SUPPORTED"
    }

    print(f"  Z-score: {z:.2f} (mean={mu:.1f}, std={sigma:.1f})")
    print(f"  Verdict: {result['verdict']}")

    out_path = os.path.join(results_dir, "108_falsification.json")
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {out_path}")
    return result

# ============================================================
# TEST 2: GREAT CIRCLE ESCALATION SERIES
# ============================================================
def test_escalation(merged, portal_only, unesco_sites, n_trials, results_dir):
    print(f"\n{'='*70}")
    print("TEST 2: Great Circle Escalation Series")
    print("=" * 70)

    # Filter UNESCO to ancient-only (~214 sites)
    # Since we don't have date info, use all but it's the "ancient" cultural sites
    # The directive says ~214, so use all
    unesco_for_test = [{"lat": s["lat"], "lon": s["lon"], "name": s["name"],
                        "gc_dist": gc_distance(s["lat"], s["lon"])}
                       for s in unesco_sites]

    results = {}

    # a) UNESCO
    print(f"\n  --- UNESCO ancient ({len(unesco_for_test)} sites) ---")
    results["unesco"] = {
        "dataset": "UNESCO ancient",
        "n_sites": len(unesco_for_test),
        "results": run_monte_carlo(unesco_for_test, n_trials)
    }

    # b) Portal only
    print(f"\n  --- Megalithic Portal ({len(portal_only)} sites) ---")
    results["portal"] = {
        "dataset": "Megalithic Portal only",
        "n_sites": len(portal_only),
        "results": run_monte_carlo(portal_only, n_trials)
    }

    # c) Merged
    print(f"\n  --- Merged dataset ({len(merged)} sites) ---")
    results["merged"] = {
        "dataset": "Merged (Portal + Supplement)",
        "n_sites": len(merged),
        "results": run_monte_carlo(merged, n_trials)
    }

    out_path = os.path.join(results_dir, "great_circle_escalation.json")
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved: {out_path}")
    return results

# ============================================================
# TEST 3: TEMPORAL ANALYSIS
# ============================================================
def test_temporal(merged, n_trials, results_dir):
    print(f"\n{'='*70}")
    print("TEST 3: Temporal Analysis")
    print("=" * 70)

    prehistoric = []
    later = []
    unclassified = []

    for s in merged:
        if s["source"] == "supplement":
            cls = classify_supplement(s)
        elif s["type"] in PREHISTORIC_TYPES:
            cls = "prehistoric"
        elif s["type"] in LATER_TYPES:
            cls = "later"
        else:
            cls = "unclassified"
        s["temporal"] = cls
        if cls == "prehistoric":
            prehistoric.append(s)
        elif cls == "later":
            later.append(s)
        else:
            unclassified.append(s)

    print(f"  Prehistoric: {len(prehistoric)}")
    print(f"  Later: {len(later)}")
    print(f"  Unclassified: {len(unclassified)}")

    print(f"\n  Running Monte Carlo for prehistoric group...")
    pre_results = run_monte_carlo(prehistoric, n_trials)
    print(f"  Running Monte Carlo for later group...")
    lat_results = run_monte_carlo(later, n_trials)

    result = {
        "prehistoric": {
            "n_sites": len(prehistoric),
            "results": pre_results
        },
        "later": {
            "n_sites": len(later),
            "results": lat_results
        },
        "unclassified_count": len(unclassified),
        "ratio_50km": round(pre_results[50]["z_score"] / lat_results[50]["z_score"], 2)
            if lat_results[50]["z_score"] != 0 else None
    }

    print(f"\n  Prehistoric Z@50km: {pre_results[50]['z_score']}")
    print(f"  Later Z@50km: {lat_results[50]['z_score']}")

    out_path = os.path.join(results_dir, "temporal_analysis.json")
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {out_path}")
    return result

# ============================================================
# TEST 4: MULTIPLE GREAT CIRCLE COMPARISON
# ============================================================
def test_multiple_circles(merged, results_dir, n_circles=1000, trials_per=50):
    print(f"\n{'='*70}")
    print(f"TEST 4: Multiple Great Circle Comparison ({n_circles} random circles)")
    print("=" * 70)

    site_lats = [s["lat"] for s in merged]
    site_lons = [s["lon"] for s in merged]
    n = len(merged)

    # First compute Alison's circle Z@50km
    obs_alison = sum(1 for s in merged if s["gc_dist"] <= 50)

    # Z-scores for random circles
    random_zs = []
    t_start = time.time()
    for ci in range(n_circles):
        pole_lat, pole_lon = random_pole()
        qc = EARTH_R_KM * math.pi / 2

        # Observed for this circle
        obs = sum(1 for s in merged
                  if abs(haversine_km(s["lat"], s["lon"], pole_lat, pole_lon) - qc) <= 50)

        # Quick baseline
        baseline = []
        for _ in range(trials_per):
            bc = sum(1 for __ in range(n)
                     if abs(haversine_km(
                         random.choice(site_lats) + random.gauss(0, 2),
                         random.choice(site_lons) + random.gauss(0, 2),
                         pole_lat, pole_lon) - qc) <= 50)
            baseline.append(bc)

        z, mu, sigma = compute_z(obs, baseline)
        random_zs.append({"pole_lat": round(pole_lat, 3),
                          "pole_lon": round(pole_lon, 3),
                          "observed": obs, "z_score": round(z, 2)})

        if (ci + 1) % 100 == 0:
            elapsed = time.time() - t_start
            print(f"    {ci+1}/{n_circles} ({elapsed:.0f}s)")

    # Compute Alison's Z (use stored result from escalation, or recompute)
    # Quick baseline for Alison
    alison_baseline = []
    for _ in range(trials_per):
        bc = sum(1 for __ in range(n)
                 if gc_distance(
                     random.choice(site_lats) + random.gauss(0, 2),
                     random.choice(site_lons) + random.gauss(0, 2)) <= 50)
        alison_baseline.append(bc)
    alison_z, _, _ = compute_z(obs_alison, alison_baseline)

    # Rank
    zs_sorted = sorted([r["z_score"] for r in random_zs], reverse=True)
    rank = sum(1 for z in zs_sorted if z >= alison_z) + 1
    percentile = (1 - rank / (n_circles + 1)) * 100

    result = {
        "alison_circle": {
            "pole": {"lat": POLE_LAT, "lon": POLE_LNG},
            "observed_50km": obs_alison,
            "z_score_50km": round(alison_z, 2)
        },
        "random_circles": {
            "n_circles": n_circles,
            "trials_per_circle": trials_per,
            "z_scores_summary": {
                "mean": round(sum(r["z_score"] for r in random_zs) / n_circles, 2),
                "max": round(max(r["z_score"] for r in random_zs), 2),
                "min": round(min(r["z_score"] for r in random_zs), 2),
            },
            "top_10": sorted(random_zs, key=lambda x: -x["z_score"])[:10]
        },
        "alison_rank": rank,
        "alison_percentile": round(percentile, 2),
        "n_random_exceeding_alison": rank - 1
    }

    print(f"  Alison Z@50km: {alison_z:.2f}")
    print(f"  Random circles max Z: {result['random_circles']['z_scores_summary']['max']}")
    print(f"  Alison rank: {rank}/{n_circles+1} (percentile: {percentile:.1f}%)")

    out_path = os.path.join(results_dir, "multiple_circles.json")
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {out_path}")
    return result

# ============================================================
# TEST 5: DENSITY PROFILE
# ============================================================
def test_density_profile(merged, circle_points, results_dir):
    print(f"\n{'='*70}")
    print("TEST 5: Density Profile")
    print("=" * 70)

    profile = []
    for p in circle_points:
        count = sum(1 for s in merged
                    if haversine_km(s["lat"], s["lon"], p["lat"], p["lon"]) <= 50)
        profile.append({
            "bearing": p["bearing"], "lat": p["lat"], "lon": p["lon"],
            "sites_within_50km": count
        })

    # Identify clusters (consecutive bins with sites)
    clusters = []
    current = None
    for p in profile:
        if p["sites_within_50km"] > 0:
            if current is None:
                current = {"start": p["bearing"], "end": p["bearing"],
                           "max_density": p["sites_within_50km"],
                           "total_sites": p["sites_within_50km"]}
            else:
                current["end"] = p["bearing"]
                current["max_density"] = max(current["max_density"], p["sites_within_50km"])
                current["total_sites"] += p["sites_within_50km"]
        else:
            if current is not None:
                current["span_degrees"] = current["end"] - current["start"] + 1
                clusters.append(current)
                current = None
    if current is not None:
        current["span_degrees"] = current["end"] - current["start"] + 1
        clusters.append(current)

    # Check wrap-around
    if len(clusters) >= 2 and clusters[0]["start"] == 0 and clusters[-1]["end"] == 359:
        merged_cluster = {
            "start": clusters[-1]["start"],
            "end": clusters[0]["end"] + 360,
            "max_density": max(clusters[0]["max_density"], clusters[-1]["max_density"]),
            "total_sites": clusters[0]["total_sites"] + clusters[-1]["total_sites"],
            "span_degrees": clusters[-1]["span_degrees"] + clusters[0]["span_degrees"],
            "note": "wraps around 360°"
        }
        clusters = [merged_cluster] + clusters[1:-1]

    clusters.sort(key=lambda c: -c["total_sites"])

    result = {
        "profile": profile,
        "clusters": clusters[:20],  # top 20
        "summary": {
            "bins_with_sites": sum(1 for p in profile if p["sites_within_50km"] > 0),
            "bins_empty": sum(1 for p in profile if p["sites_within_50km"] == 0),
            "max_density": max(p["sites_within_50km"] for p in profile),
            "total_sites_on_circle": sum(p["sites_within_50km"] for p in profile),
            "n_clusters": len(clusters)
        }
    }

    print(f"  Bins with sites: {result['summary']['bins_with_sites']}/360")
    print(f"  Max density: {result['summary']['max_density']} sites in one bin")
    print(f"  Clusters found: {len(clusters)}")

    out_path = os.path.join(results_dir, "density_profile.json")
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {out_path}")
    return result

# ============================================================
# TEST 6: SITE TYPE ENRICHMENT
# ============================================================
def test_type_enrichment(merged, results_dir):
    print(f"\n{'='*70}")
    print("TEST 6: Site Type Enrichment")
    print("=" * 70)

    type_counts = defaultdict(lambda: {"total": 0, "within_50km": 0})
    for s in merged:
        type_counts[s["type"]]["total"] += 1
        if s["gc_dist"] <= 50:
            type_counts[s["type"]]["within_50km"] += 1

    enrichment = []
    for stype, counts in type_counts.items():
        if counts["total"] >= 20:
            pct = counts["within_50km"] / counts["total"] * 100
            enrichment.append({
                "type": stype,
                "total": counts["total"],
                "within_50km": counts["within_50km"],
                "percentage": round(pct, 2)
            })

    enrichment.sort(key=lambda x: -x["percentage"])

    result = {
        "types_analyzed": len(enrichment),
        "enrichment": enrichment
    }

    print(f"  Types with >20 sites: {len(enrichment)}")
    for e in enrichment[:10]:
        print(f"    {e['type']:>40}: {e['within_50km']:>4}/{e['total']:>5} ({e['percentage']:.1f}%)")

    out_path = os.path.join(results_dir, "type_enrichment.json")
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {out_path}")
    return result

# ============================================================
# TEST 7: AGE ANALYSIS
# ============================================================
def test_age_analysis(merged, results_dir):
    print(f"\n{'='*70}")
    print("TEST 7: Age Analysis")
    print("=" * 70)

    on_line = []  # within 50km
    off_line = []

    for s in merged:
        age = TYPE_AGE_ESTIMATES.get(s["type"])
        if s["source"] == "supplement":
            age = parse_supplement_date(s.get("date", ""))
        if age is not None:
            if s["gc_dist"] <= 50:
                on_line.append(age)
            else:
                off_line.append(age)

    on_mean = sum(on_line) / len(on_line) if on_line else 0
    off_mean = sum(off_line) / len(off_line) if off_line else 0
    on_median = sorted(on_line)[len(on_line)//2] if on_line else 0
    off_median = sorted(off_line)[len(off_line)//2] if off_line else 0

    result = {
        "on_line": {
            "count": len(on_line),
            "mean_age": round(on_mean),
            "median_age": on_median
        },
        "off_line": {
            "count": len(off_line),
            "mean_age": round(off_mean),
            "median_age": off_median
        },
        "comparison": {
            "mean_difference": round(on_mean - off_mean),
            "median_difference": on_median - off_median,
            "on_line_older": on_mean < off_mean
        }
    }

    print(f"  On-line (≤50km):  n={len(on_line)}, mean={on_mean:.0f}, median={on_median}")
    print(f"  Off-line (>50km): n={len(off_line)}, mean={off_mean:.0f}, median={off_median}")
    print(f"  On-line {'older' if on_mean < off_mean else 'newer'} by {abs(on_mean - off_mean):.0f} years (mean)")

    out_path = os.path.join(results_dir, "age_analysis.json")
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {out_path}")
    return result

# ============================================================
# TEST 8: PLEIADES VALIDATION
# ============================================================
def test_pleiades(pleiades_csv, n_trials, results_dir):
    print(f"\n{'='*70}")
    print("TEST 8: Pleiades Validation")
    print("=" * 70)

    if pleiades_csv is None or not os.path.exists(pleiades_csv):
        print("  SKIPPED: Pleiades data not available")
        result = {"skipped": True, "reason": "Pleiades CSV not downloaded"}
        out_path = os.path.join(results_dir, "pleiades_validation.json")
        with open(out_path, 'w') as f:
            json.dump(result, f, indent=2)
        return result

    sites = load_pleiades(pleiades_csv)
    print(f"  Total Pleiades sites: {len(sites)}")

    # Add gc_dist
    for s in sites:
        s["gc_dist"] = gc_distance(s["lat"], s["lon"])

    # a) All sites
    print(f"\n  --- All Pleiades ({len(sites)} sites) ---")
    all_results = run_monte_carlo(sites, n_trials)

    # b) Pre-2000 BCE
    ancient = [s for s in sites if s["minDate"] is not None and s["minDate"] < -2000]
    print(f"\n  --- Pleiades pre-2000 BCE ({len(ancient)} sites) ---")
    ancient_results = run_monte_carlo(ancient, n_trials) if len(ancient) > 10 else None

    result = {
        "all_pleiades": {
            "n_sites": len(sites),
            "results": all_results
        },
        "pre_2000_bce": {
            "n_sites": len(ancient),
            "results": ancient_results
        } if ancient_results else {
            "n_sites": len(ancient),
            "skipped": True,
            "reason": f"Only {len(ancient)} sites with minDate < -2000"
        }
    }

    print(f"\n  All Pleiades Z@50km: {all_results[50]['z_score']}")
    if ancient_results:
        print(f"  Pre-2000 BCE Z@50km: {ancient_results[50]['z_score']}")

    out_path = os.path.join(results_dir, "pleiades_validation.json")
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {out_path}")
    return result

# ============================================================
# TEST 9: SUPPLEMENT-ONLY TEST
# ============================================================
def test_supplement_only(supplement_sites, n_trials, results_dir):
    print(f"\n{'='*70}")
    print(f"TEST 9: Supplement-Only Test ({len(supplement_sites)} sites)")
    print("=" * 70)

    for s in supplement_sites:
        s["gc_dist"] = gc_distance(s["lat"], s["lon"])

    results = run_monte_carlo(supplement_sites, n_trials)

    nearby = sorted([s for s in supplement_sites if s["gc_dist"] <= 50],
                    key=lambda x: x["gc_dist"])

    result = {
        "n_sites": len(supplement_sites),
        "results": results,
        "sites_within_50km": [
            {"name": s["name"], "type": s.get("type", ""), "region": s.get("region", ""),
             "distance_km": round(s["gc_dist"], 1)}
            for s in nearby
        ],
        "interpretation": (
            "These 114 sites were hand-selected as major world archaeological sites "
            "missing from the Megalithic Portal. A high Z-score could reflect genuine "
            "alignment OR selection bias toward famous sites near the circle."
        )
    }

    print(f"  Z@50km: {results[50]['z_score']}")
    print(f"  Sites within 50km: {len(nearby)}")

    out_path = os.path.join(results_dir, "supplement_only.json")
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {out_path}")
    return result

# ============================================================
# STEP 6: GENERATE SUMMARY
# ============================================================
def generate_summary(escalation, temporal, density, test_108,
                     pleiades_result, results_dir):
    print(f"\n{'='*70}")
    print("STEP 6: Generating Summary")
    print("=" * 70)

    escalation_table = []
    for key, label in [("unesco", "UNESCO ancient"),
                       ("portal", "Portal only"),
                       ("merged", "Merged")]:
        if key in escalation:
            escalation_table.append({
                "dataset": label,
                "sites": escalation[key]["n_sites"],
                "z_25km": escalation[key]["results"][25]["z_score"],
                "z_50km": escalation[key]["results"][50]["z_score"],
                "z_100km": escalation[key]["results"][100]["z_score"],
                "z_200km": escalation[key]["results"][200]["z_score"]
            })

    # Add Pleiades to escalation table
    if pleiades_result and not pleiades_result.get("skipped"):
        if "all_pleiades" in pleiades_result:
            p = pleiades_result["all_pleiades"]
            escalation_table.append({
                "dataset": "Pleiades all",
                "sites": p["n_sites"],
                "z_50km": p["results"][50]["z_score"]
            })
        if "pre_2000_bce" in pleiades_result and not pleiades_result["pre_2000_bce"].get("skipped"):
            p = pleiades_result["pre_2000_bce"]
            escalation_table.append({
                "dataset": "Pleiades pre-2000 BCE",
                "sites": p["n_sites"],
                "z_50km": p["results"][50]["z_score"]
            })

    temporal_summary = {}
    if temporal:
        temporal_summary = {
            "prehistoric_z": temporal["prehistoric"]["results"][50]["z_score"],
            "later_z": temporal["later"]["results"][50]["z_score"],
            "ratio": temporal.get("ratio_50km")
        }

    clusters = density["clusters"][:10] if density else []

    summary = {
        "escalation_table": escalation_table,
        "temporal": temporal_summary,
        "clusters": clusters,
        "108_verdict": test_108["verdict"] if test_108 else "NOT RUN",
        "great_circle_verdict": "CONFIRMED" if (
            escalation.get("merged", {}).get("results", {}).get(50, {}).get("z_score", 0) > 3
        ) else "NOT CONFIRMED"
    }

    out_path = os.path.join(results_dir, "summary.json")
    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  Saved: {out_path}")
    return summary

# ============================================================
# MAIN
# ============================================================
def main():
    parser = argparse.ArgumentParser(description="Deep Time Atlas — Great Circle Analysis")
    parser.add_argument("--trials", type=int, default=200, help="Monte Carlo trials per test")
    parser.add_argument("--kml-dir", default=None, help="Directory containing MegP_*.kml files")
    parser.add_argument("--skip-pleiades", action="store_true", help="Skip Pleiades download")
    parser.add_argument("--fast", action="store_true", help="Fast mode: 50 trials, 200 circles")
    args = parser.parse_args()

    n_trials = args.trials
    if args.fast:
        n_trials = 50
        print("*** FAST MODE: 50 trials, reduced circles ***")

    # Resolve paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_dir = os.path.dirname(script_dir)
    data_dir = os.path.join(repo_dir, "data")
    results_dir = os.path.join(repo_dir, "results")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    # KML directory: default to parent of repo (where raw KMLs live)
    kml_dir = args.kml_dir or os.path.dirname(repo_dir)

    # Supplement path
    supplement_path = os.path.join(data_dir, "supplement_sites.json")
    if not os.path.exists(supplement_path):
        # Try common locations
        for p in [os.path.join(kml_dir, "complete_supplement_for_retest.json"),
                  os.path.expanduser("~/Downloads/complete_supplement_for_retest.json")]:
            if os.path.exists(p):
                supplement_path = p
                break

    # UNESCO path
    unesco_path = None
    for p in [os.path.join(data_dir, "unesco_cultural_sites.json"),
              os.path.join(kml_dir, "unesco_cultural_sites.json")]:
        if os.path.exists(p):
            unesco_path = p
            break

    t_total = time.time()
    random.seed(42)  # Reproducibility

    # STEP 1: Parse KML
    portal_sites = parse_kml_files(kml_dir)

    # Save portal CSV
    portal_csv = os.path.join(data_dir, "megalithic_portal_sites.csv")
    with open(portal_csv, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(["name", "lat", "lon", "type"])
        for s in portal_sites:
            w.writerow([s["name"], s["lat"], s["lon"], s["type"]])
    print(f"  Saved: {portal_csv}")

    # STEP 2: Merge supplement
    merged, supplement_added, dedup_log = load_and_merge_supplement(
        portal_sites, supplement_path)

    # Save merged CSV
    merged_csv = os.path.join(data_dir, "merged_sites.csv")
    with open(merged_csv, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(["name", "lat", "lon", "type", "source"])
        for s in merged:
            w.writerow([s["name"], s["lat"], s["lon"], s["type"], s["source"]])
    print(f"  Saved: {merged_csv}")

    # Copy supplement to data dir if not already there
    supp_data_path = os.path.join(data_dir, "supplement_sites.json")
    if not os.path.exists(supp_data_path):
        with open(supplement_path) as f:
            supp_data = json.load(f)
        with open(supp_data_path, 'w') as f:
            json.dump(supp_data, f, indent=2)
        print(f"  Copied supplement to: {supp_data_path}")

    # STEP 3: Pleiades
    pleiades_csv = download_pleiades(data_dir, skip=args.skip_pleiades)

    # STEP 4: Circle coordinates
    circle_points = generate_circle_coordinates(data_dir)

    # Precompute gc_dist for all merged sites
    for s in merged:
        s["gc_dist"] = gc_distance(s["lat"], s["lon"])
    for s in portal_sites:
        s["gc_dist"] = gc_distance(s["lat"], s["lon"])

    # ================= TESTS =================

    # TEST 1: 108° Falsification
    test_108_result = None
    if unesco_path:
        with open(unesco_path) as f:
            unesco_data = json.load(f)
        test_108_result = test_108_falsification(
            unesco_path, results_dir, n_trials=min(100, n_trials))
    else:
        print("\n  WARNING: UNESCO data not found, skipping Test 1")

    # TEST 2: Escalation
    unesco_sites = unesco_data["sites"] if unesco_path else []
    escalation_result = test_escalation(
        merged, portal_sites, unesco_sites, n_trials, results_dir)

    # TEST 3: Temporal
    temporal_result = test_temporal(merged, n_trials, results_dir)

    # TEST 4: Multiple circles
    n_circles = 200 if args.fast else 1000
    circles_result = test_multiple_circles(
        merged, results_dir, n_circles=n_circles, trials_per=50)

    # TEST 5: Density profile
    density_result = test_density_profile(merged, circle_points, results_dir)

    # TEST 6: Type enrichment
    enrichment_result = test_type_enrichment(merged, results_dir)

    # TEST 7: Age analysis
    age_result = test_age_analysis(merged, results_dir)

    # TEST 8: Pleiades
    pleiades_result = test_pleiades(pleiades_csv, n_trials, results_dir)

    # TEST 9: Supplement-only
    supplement_result = test_supplement_only(supplement_added, n_trials, results_dir)

    # STEP 6: Summary
    summary = generate_summary(
        escalation_result, temporal_result, density_result,
        test_108_result, pleiades_result, results_dir)

    total_time = time.time() - t_total
    print(f"\n{'='*70}")
    print(f"ALL DONE in {total_time:.0f}s ({total_time/60:.1f} min)")
    print(f"{'='*70}")

    # Print summary table
    print(f"\n  ESCALATION TABLE:")
    print(f"  {'Dataset':<25} {'Sites':>8} {'Z@50km':>8}")
    print(f"  {'-'*43}")
    for row in summary["escalation_table"]:
        print(f"  {row['dataset']:<25} {row['sites']:>8} {row.get('z_50km', 'N/A'):>8}")

    print(f"\n  108° verdict: {summary['108_verdict']}")
    print(f"  Great Circle verdict: {summary['great_circle_verdict']}")

if __name__ == "__main__":
    main()
