#!/usr/bin/env python3
"""
Directive 12, Script 00: Subterranean Site Database
====================================================
Compiles a global database of subterranean archaeological sites from:
  1. Wikidata (filtered by subterranean keywords in name/type)
  2. Pleiades (filtered by featureTypes containing cave/underground terms)
  3. Megalithic Portal KMLs (Cave_or_Rock_Shelter, Rock_Cut_Tomb, Souterrain)
  4. UNESCO World Heritage subterranean sites (hand-curated)
  5. Regional datasets (Cappadocia, India, Peru, Levant)

Each site is classified into one of 5 subterranean types:
  - natural_cave: Natural cave with human use (paintings, habitation, burial)
  - modified_cave: Cave with carved extensions or architectural additions
  - rock_cut_monumental: Rock-cut temples, tombs, planned chambers
  - artificial_complex: Underground cities, tunnel networks, catacombs
  - utilitarian: Mines, quarries, storage, aqueducts

Output: subterranean_sites_master.csv
"""

import csv
import json
import math
import os
import re
import xml.etree.ElementTree as ET

# ── Constants ──────────────────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
R_EARTH = 6371.0
QUARTER_CIRC = R_EARTH * math.pi / 2  # ~10,007.5 km

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(BASE, "..", "..", "..")  # megalith_site_research/

# Minimum thresholds for "subterranean" (from directive bias section)
# >3m below surface OR >10m tunnel OR >2 connected chambers
# Applied loosely since most databases don't have depth data


# ── Geo utilities ──────────────────────────────────────────────────────
def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlat, dlon = la2 - la1, lo2 - lo1
    a = math.sin(dlat / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin(dlon / 2) ** 2
    return R_EARTH * 2 * math.asin(math.sqrt(min(1.0, a)))


def gc_distance(lat, lon):
    d = haversine_km(POLE_LAT, POLE_LON, lat, lon)
    return abs(d - QUARTER_CIRC)


def assign_cluster(lat, lon):
    if haversine_km(-27.1, -109.3, lat, lon) < 200:
        return "Easter Island"
    if 20 <= lat <= 36 and 25 <= lon <= 40:
        return "Egypt / Levant"
    if -25 <= lat <= 0 and -80 <= lon <= -60:
        return "Peru / Andes"
    if 25 <= lat <= 38 and 44 <= lon <= 60:
        return "Iran / Persia"
    if 20 <= lat <= 35 and 60 <= lon <= 80:
        return "Indus Valley"
    if 0 <= lat <= 25 and 95 <= lon <= 115:
        return "Southeast Asia"
    return "Other"


# ── Type classification ────────────────────────────────────────────────
SUBTERRANEAN_KEYWORDS = [
    "cave", "grotto", "grotta", "cueva", "höhle", "cavern",
    "hypogeum", "hypogea", "hypogée",
    "underground", "subterranean", "souterrain", "fogou",
    "rock-cut", "rock cut", "rupestre",
    "catacomb", "ossuary",
    "tunnel", "gallery", "shaft",
    "lava tube",
    "mine", "quarry",
    "puquio", "qanat", "foggara",
    "cenote",
]

# Keywords that indicate MONUMENTAL rock-cut architecture
MONUMENTAL_KEYWORDS = [
    "temple", "church", "chapel", "mosque", "monastery",
    "tomb", "necropolis", "mausoleum", "sarcophag",
    "palace", "city", "dwelling", "habitat",
    "serapeum", "hypogeum", "catacomb",
    "petra", "ajanta", "ellora", "barabar",
    "derinkuyu", "kaymakli",
]

UTILITARIAN_KEYWORDS = [
    "mine", "quarry", "aqueduct", "puquio", "qanat",
    "cistern", "well", "drain", "sewer", "storage",
    "foggara", "irrigation",
]


def is_subterranean(name, type_str, description=""):
    """Check if a site is subterranean based on textual fields."""
    text = f"{name} {type_str} {description}".lower()
    return any(kw in text for kw in SUBTERRANEAN_KEYWORDS)


def classify_subterranean(name, type_str, description=""):
    """Classify a subterranean site into one of 5 categories."""
    text = f"{name} {type_str} {description}".lower()

    # Check utilitarian first (mines, aqueducts)
    if any(kw in text for kw in UTILITARIAN_KEYWORDS):
        return "utilitarian"

    # Artificial underground complex (cities, multi-level systems)
    if any(kw in text for kw in ["underground city", "derinkuyu", "kaymakli",
                                   "özkonak", "underground dwelling",
                                   "subterranean city", "catacomb"]):
        return "artificial_complex"

    # Rock-cut monumental (temples, tombs carved from rock)
    if "rock-cut" in text or "rock cut" in text or "rupestre" in text:
        if any(kw in text for kw in ["temple", "church", "tomb", "necropolis"]):
            return "rock_cut_monumental"
        return "rock_cut_monumental"  # rock-cut anything is monumental

    # Check if it has monumental indicators
    if any(kw in text for kw in ["temple", "tomb", "necropolis", "serapeum",
                                   "hypogeum", "mausoleum", "palace"]):
        return "rock_cut_monumental"

    # Modified cave (paintings, carvings, architectural additions)
    if any(kw in text for kw in ["painting", "art", "carving", "decorated",
                                   "burial", "ritual", "modified"]):
        return "modified_cave"

    # Default: natural cave with human use
    return "natural_cave"


def estimate_date_bp(name, type_str, description="", date_field=None):
    """Rough date estimate in years BP (before 1950). Returns None if unknown."""
    if date_field is not None:
        try:
            yr = float(date_field)
            if yr < 0:  # BCE
                return int(1950 - yr)
            elif yr < 1950:
                return int(1950 - yr)
            return None
        except (ValueError, TypeError):
            pass

    text = f"{name} {type_str} {description}".lower()

    # Very rough heuristics
    if any(kw in text for kw in ["paleolithic", "palaeolithic"]):
        return 40000
    if any(kw in text for kw in ["mesolithic", "epipaleolithic"]):
        return 12000
    if any(kw in text for kw in ["neolithic"]):
        return 7000
    if any(kw in text for kw in ["bronze age"]):
        return 4000
    if any(kw in text for kw in ["iron age"]):
        return 2800
    if any(kw in text for kw in ["roman", "hellenistic"]):
        return 2200
    if any(kw in text for kw in ["medieval", "mediaeval", "byzantine"]):
        return 1000

    return None


# ── Data Loaders ───────────────────────────────────────────────────────

def load_wikidata_subterranean():
    """Filter Wikidata archaeological sites for subterranean entries."""
    path = os.path.join(ROOT, "data", "wikidata", "wikidata_archaeological_sites.csv")
    sites = []
    total = 0
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            total += 1
            try:
                lat = float(row["lat"])
                lon = float(row["lon"])
            except (ValueError, KeyError):
                continue

            name = row.get("name", "")
            type_str = row.get("type", "")
            inception = row.get("inception", "")

            if not is_subterranean(name, type_str):
                continue

            sites.append({
                "name": name.strip() or "Unnamed",
                "lat": lat,
                "lon": lon,
                "source": "wikidata",
                "type_raw": type_str,
                "subtype": classify_subterranean(name, type_str),
                "date_bp": estimate_date_bp(name, type_str, date_field=inception),
                "gc_distance_km": round(gc_distance(lat, lon), 1),
                "cluster": assign_cluster(lat, lon),
            })

    print(f"Wikidata: {len(sites)} subterranean sites (from {total} total)")
    return sites


def load_pleiades_subterranean():
    """Filter Pleiades for cave/underground feature types."""
    path = os.path.join(ROOT, "data", "pleiades", "pleiades-places-latest.csv")
    sites = []
    total = 0
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            total += 1
            try:
                lat = float(row["reprLat"])
                lon = float(row["reprLong"])
            except (ValueError, KeyError):
                continue

            name = row.get("title", "")
            feature = row.get("featureTypes", "")
            desc = row.get("description", "")
            min_date = row.get("minDate", "")

            if not is_subterranean(name, feature, desc):
                continue

            sites.append({
                "name": name.strip() or "Unnamed",
                "lat": lat,
                "lon": lon,
                "source": "pleiades",
                "type_raw": feature,
                "subtype": classify_subterranean(name, feature, desc),
                "date_bp": estimate_date_bp(name, feature, desc, date_field=min_date),
                "gc_distance_km": round(gc_distance(lat, lon), 1),
                "cluster": assign_cluster(lat, lon),
            })

    print(f"Pleiades: {len(sites)} subterranean sites (from {total} total)")
    return sites


def parse_kml_coordinates(kml_path):
    """Extract placemarks with coordinates from a KML file."""
    sites = []
    try:
        tree = ET.parse(kml_path)
    except Exception:
        return sites

    root = tree.getroot()
    # Handle KML namespace
    ns = ""
    if root.tag.startswith("{"):
        ns = root.tag.split("}")[0] + "}"

    for pm in root.iter(f"{ns}Placemark"):
        name_el = pm.find(f"{ns}name")
        name = name_el.text.strip() if name_el is not None and name_el.text else "Unnamed"

        # Get coordinates from Point
        point = pm.find(f".//{ns}Point/{ns}coordinates")
        if point is None or point.text is None:
            continue

        coords = point.text.strip().split(",")
        if len(coords) < 2:
            continue

        try:
            lon = float(coords[0])
            lat = float(coords[1])
        except ValueError:
            continue

        sites.append({"name": name, "lat": lat, "lon": lon})

    return sites


def load_megalithic_portal_subterranean():
    """Load Cave, Rock Cut Tomb, and Souterrain KMLs from Megalithic Portal."""
    portal_dir = os.path.join(ROOT, "data", "megalithic_portal")
    kml_files = {
        "Cave_or_Rock_Shelter": "natural_cave",
        "Rock_Cut_Tomb": "rock_cut_monumental",
        "Souterrain_(Fogou,_Earth_House)": "artificial_complex",
    }

    all_sites = []
    for kml_name, default_subtype in kml_files.items():
        path = os.path.join(portal_dir, f"MegP_{kml_name}.kml")
        if not os.path.exists(path):
            print(f"  Warning: {path} not found")
            continue

        raw = parse_kml_coordinates(path)
        for s in raw:
            subtype = classify_subterranean(s["name"], kml_name)
            # Override if the classifier returns natural_cave but KML is Rock_Cut_Tomb
            if default_subtype == "rock_cut_monumental" and subtype == "natural_cave":
                subtype = "rock_cut_monumental"
            if default_subtype == "artificial_complex" and subtype == "natural_cave":
                subtype = "artificial_complex"

            all_sites.append({
                "name": s["name"],
                "lat": s["lat"],
                "lon": s["lon"],
                "source": "megalithic_portal",
                "type_raw": kml_name.replace("_", " "),
                "subtype": subtype,
                "date_bp": estimate_date_bp(s["name"], kml_name),
                "gc_distance_km": round(gc_distance(s["lat"], s["lon"]), 1),
                "cluster": assign_cluster(s["lat"], s["lon"]),
            })

    print(f"Megalithic Portal: {len(all_sites)} subterranean sites")
    return all_sites


def load_curated_sites():
    """Hand-curated major subterranean sites from the directive."""
    sites = [
        # ── Egypt ──
        {"name": "Osiris Shaft (Giza)", "lat": 29.976, "lon": 31.131,
         "subtype": "rock_cut_monumental", "date_bp": 4500},
        {"name": "Subterranean Chamber, Great Pyramid", "lat": 29.9792, "lon": 31.1342,
         "subtype": "rock_cut_monumental", "date_bp": 4550},
        {"name": "Saqqara Serapeum", "lat": 29.871, "lon": 31.215,
         "subtype": "artificial_complex", "date_bp": 3350},
        {"name": "Valley of the Kings (KV complex)", "lat": 25.740, "lon": 32.601,
         "subtype": "rock_cut_monumental", "date_bp": 3500},
        {"name": "Osireion at Abydos", "lat": 26.185, "lon": 31.919,
         "subtype": "rock_cut_monumental", "date_bp": 3300},
        {"name": "Ibis Catacombs, Saqqara", "lat": 29.868, "lon": 31.217,
         "subtype": "artificial_complex", "date_bp": 2500},

        # ── Jordan / Levant ──
        {"name": "Petra (rock-cut city)", "lat": 30.329, "lon": 35.444,
         "subtype": "rock_cut_monumental", "date_bp": 2350},
        {"name": "Beit Guvrin-Maresha caves", "lat": 31.614, "lon": 34.896,
         "subtype": "artificial_complex", "date_bp": 2800},
        {"name": "Qumran Caves (Dead Sea Scrolls)", "lat": 31.741, "lon": 35.459,
         "subtype": "modified_cave", "date_bp": 2100},

        # ── Iran ──
        {"name": "Naqsh-e Rostam (rock-cut tombs)", "lat": 29.988, "lon": 52.874,
         "subtype": "rock_cut_monumental", "date_bp": 2500},
        {"name": "Persepolis drainage tunnels", "lat": 29.935, "lon": 52.891,
         "subtype": "utilitarian", "date_bp": 2500},

        # ── India ──
        {"name": "Ajanta Caves", "lat": 20.552, "lon": 75.700,
         "subtype": "rock_cut_monumental", "date_bp": 2150},
        {"name": "Ellora Caves", "lat": 20.027, "lon": 75.179,
         "subtype": "rock_cut_monumental", "date_bp": 1550},
        {"name": "Barabar Caves", "lat": 25.005, "lon": 85.063,
         "subtype": "rock_cut_monumental", "date_bp": 2200},
        {"name": "Elephanta Caves", "lat": 18.963, "lon": 72.932,
         "subtype": "rock_cut_monumental", "date_bp": 1450},
        {"name": "Kanheri Caves", "lat": 19.209, "lon": 72.907,
         "subtype": "rock_cut_monumental", "date_bp": 1950},
        {"name": "Udayagiri Caves", "lat": 23.480, "lon": 75.772,
         "subtype": "rock_cut_monumental", "date_bp": 1600},

        # ── Cappadocia ──
        {"name": "Derinkuyu underground city", "lat": 38.374, "lon": 34.735,
         "subtype": "artificial_complex", "date_bp": 2750},
        {"name": "Kaymakli underground city", "lat": 38.463, "lon": 34.752,
         "subtype": "artificial_complex", "date_bp": 2750},
        {"name": "Özkonak underground city", "lat": 38.680, "lon": 34.730,
         "subtype": "artificial_complex", "date_bp": 2500},
        {"name": "Mazi underground city", "lat": 38.555, "lon": 34.665,
         "subtype": "artificial_complex", "date_bp": 2500},
        {"name": "Tatlarin underground city", "lat": 38.510, "lon": 34.558,
         "subtype": "artificial_complex", "date_bp": 2400},
        {"name": "Gaziemir underground city", "lat": 38.300, "lon": 34.700,
         "subtype": "artificial_complex", "date_bp": 2400},
        {"name": "Saratli underground city", "lat": 38.520, "lon": 34.810,
         "subtype": "artificial_complex", "date_bp": 2400},
        {"name": "Acigöl underground city", "lat": 38.558, "lon": 34.517,
         "subtype": "artificial_complex", "date_bp": 2300},

        # ── Peru ──
        {"name": "Chavín de Huántar galleries", "lat": -9.593, "lon": -77.177,
         "subtype": "artificial_complex", "date_bp": 3000},
        {"name": "Nazca puquio (Cantalloc)", "lat": -14.828, "lon": -75.119,
         "subtype": "utilitarian", "date_bp": 1500},
        {"name": "Nazca puquio (Ocongalla)", "lat": -14.810, "lon": -75.145,
         "subtype": "utilitarian", "date_bp": 1500},

        # ── Easter Island ──
        {"name": "Ana Kai Tangata (painted cave)", "lat": -27.175, "lon": -109.435,
         "subtype": "modified_cave", "date_bp": 500},
        {"name": "Ana Te Pahu (lava tube complex)", "lat": -27.100, "lon": -109.400,
         "subtype": "natural_cave", "date_bp": 700},

        # ── Malta (near but not on circle) ──
        {"name": "Ħal Saflieni Hypogeum", "lat": 35.869, "lon": 14.507,
         "subtype": "artificial_complex", "date_bp": 5000},

        # ── Turkey (non-Cappadocia) ──
        {"name": "Basilica Cistern (Istanbul)", "lat": 41.008, "lon": 28.978,
         "subtype": "artificial_complex", "date_bp": 1420},

        # ── China ──
        {"name": "Mogao Caves (Dunhuang)", "lat": 40.042, "lon": 94.809,
         "subtype": "rock_cut_monumental", "date_bp": 1590},
        {"name": "Longmen Grottoes", "lat": 34.566, "lon": 112.471,
         "subtype": "rock_cut_monumental", "date_bp": 1550},
        {"name": "Yungang Grottoes", "lat": 40.112, "lon": 113.133,
         "subtype": "rock_cut_monumental", "date_bp": 1490},

        # ── Ethiopia ──
        {"name": "Lalibela rock-cut churches", "lat": 12.032, "lon": 39.043,
         "subtype": "rock_cut_monumental", "date_bp": 750},

        # ── Italy ──
        {"name": "Catacombs of Rome (San Callisto)", "lat": 41.861, "lon": 12.514,
         "subtype": "artificial_complex", "date_bp": 1750},
        {"name": "Catacombs of Paris", "lat": 48.834, "lon": 2.332,
         "subtype": "artificial_complex", "date_bp": 230},
        {"name": "Sassi di Matera", "lat": 40.666, "lon": 16.611,
         "subtype": "rock_cut_monumental", "date_bp": 9000},

        # ── France ──
        {"name": "Lascaux Cave", "lat": 45.054, "lon": 1.168,
         "subtype": "modified_cave", "date_bp": 17000},
        {"name": "Chauvet Cave", "lat": 44.387, "lon": 4.418,
         "subtype": "modified_cave", "date_bp": 36000},

        # ── Spain ──
        {"name": "Altamira Cave", "lat": 43.378, "lon": -4.120,
         "subtype": "modified_cave", "date_bp": 15000},

        # ── Indonesia ──
        {"name": "Maros-Pangkep caves (Sulawesi)", "lat": -4.985, "lon": 119.643,
         "subtype": "modified_cave", "date_bp": 45000},

        # ── South Africa ──
        {"name": "Blombos Cave", "lat": -34.412, "lon": 21.223,
         "subtype": "modified_cave", "date_bp": 100000},

        # ── Australia ──
        {"name": "Nawarla Gabarnmang", "lat": -13.029, "lon": 132.909,
         "subtype": "modified_cave", "date_bp": 28000},

        # ── PNG ──
        {"name": "Niah Cave (Borneo)", "lat": 3.812, "lon": 113.768,
         "subtype": "natural_cave", "date_bp": 40000},

        # ── Greece ──
        {"name": "Diros Caves (Mani)", "lat": 36.639, "lon": 22.385,
         "subtype": "natural_cave", "date_bp": 6000},

        # ── Mexico ──
        {"name": "Cenotes of Yucatan (Ik Kil)", "lat": 20.666, "lon": -88.551,
         "subtype": "natural_cave", "date_bp": 1200},

        # ── Afghanistan ──
        {"name": "Bamiyan cave complexes", "lat": 34.833, "lon": 67.827,
         "subtype": "rock_cut_monumental", "date_bp": 1550},

        # ── Control: major sites clearly OFF the circle ──
        {"name": "Cu Chi Tunnels (Vietnam)", "lat": 11.143, "lon": 106.463,
         "subtype": "artificial_complex", "date_bp": 50},
        {"name": "Edinburgh Vaults", "lat": 55.949, "lon": -3.187,
         "subtype": "artificial_complex", "date_bp": 230},
    ]

    output = []
    for s in sites:
        output.append({
            "name": s["name"],
            "lat": s["lat"],
            "lon": s["lon"],
            "source": "curated",
            "type_raw": s["subtype"],
            "subtype": s["subtype"],
            "date_bp": s.get("date_bp"),
            "gc_distance_km": round(gc_distance(s["lat"], s["lon"]), 1),
            "cluster": assign_cluster(s["lat"], s["lon"]),
        })

    print(f"Curated: {len(output)} major subterranean sites")
    return output


# ── Deduplication ──────────────────────────────────────────────────────
def deduplicate(all_sites, radius_km=2.0):
    """Deduplicate within radius_km, preferring curated > portal > pleiades > wikidata."""
    priority = {"curated": 0, "megalithic_portal": 1, "pleiades": 2, "wikidata": 3}

    def score(s):
        has_name = 1 if s["name"] != "Unnamed" else 0
        return (-priority.get(s["source"], 99), has_name)

    all_sites.sort(key=score, reverse=True)

    kept = []
    for site in all_sites:
        is_dup = False
        for existing in kept:
            if existing["source"] == site["source"]:
                continue
            if haversine_km(site["lat"], site["lon"], existing["lat"], existing["lon"]) < radius_km:
                is_dup = True
                break
        if not is_dup:
            kept.append(site)

    return kept


# ── Main ───────────────────────────────────────────────────────────────
def build_database():
    print("=" * 70)
    print("SUBTERRANEAN ARCHAEOLOGY — Building Master Database")
    print("=" * 70)

    wikidata = load_wikidata_subterranean()
    pleiades = load_pleiades_subterranean()
    portal = load_megalithic_portal_subterranean()
    curated = load_curated_sites()

    all_sites = curated + portal + pleiades + wikidata
    print(f"\nTotal before dedup: {len(all_sites)}")

    deduped = deduplicate(all_sites)
    print(f"Total after dedup:  {len(deduped)}")

    # Summary stats
    by_source = {}
    by_subtype = {}
    for s in deduped:
        by_source[s["source"]] = by_source.get(s["source"], 0) + 1
        by_subtype[s["subtype"]] = by_subtype.get(s["subtype"], 0) + 1

    print(f"\n--- By source ---")
    for src, count in sorted(by_source.items()):
        print(f"  {src}: {count}")

    print(f"\n--- By subterranean type ---")
    for st, count in sorted(by_subtype.items(), key=lambda x: -x[1]):
        print(f"  {st}: {count}")

    # Distance distribution
    within_50 = sum(1 for s in deduped if s["gc_distance_km"] <= 50)
    within_100 = sum(1 for s in deduped if s["gc_distance_km"] <= 100)
    within_200 = sum(1 for s in deduped if s["gc_distance_km"] <= 200)
    within_500 = sum(1 for s in deduped if s["gc_distance_km"] <= 500)

    print(f"\n--- Distance to Great Circle ---")
    print(f"  Within  50 km: {within_50}")
    print(f"  Within 100 km: {within_100}")
    print(f"  Within 200 km: {within_200}")
    print(f"  Within 500 km: {within_500}")
    print(f"  Total:         {len(deduped)}")

    # Save CSV
    out_path = os.path.join(BASE, "subterranean_sites_master.csv")
    fieldnames = ["name", "lat", "lon", "source", "type_raw", "subtype",
                  "date_bp", "gc_distance_km", "cluster"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for s in sorted(deduped, key=lambda x: x["gc_distance_km"]):
            writer.writerow(s)

    print(f"\nSaved {out_path} ({len(deduped)} sites)")

    return deduped


# Export for use by other scripts
SUBTERRANEAN_SITES = None

def get_sites():
    global SUBTERRANEAN_SITES
    if SUBTERRANEAN_SITES is None:
        SUBTERRANEAN_SITES = build_database()
    return SUBTERRANEAN_SITES


if __name__ == "__main__":
    get_sites()
