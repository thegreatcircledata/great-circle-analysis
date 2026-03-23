#!/usr/bin/env python3
"""
Regional Onset Analysis: Earliest monumental construction within 200km of the Great Circle
Tests diffusion vs independent invention by comparing onset dates across regions.
"""

import json
import csv
import math
import re
import os
from collections import defaultdict
from pathlib import Path

BASE = Path("/Users/elliotallan/megalith_site_research")

# ── Great Circle definition ──
with open(BASE / "github-repo/data/circle_coordinates.json") as f:
    circle_data = json.load(f)
POLE = (circle_data["pole"]["lat"], circle_data["pole"]["lon"])
CIRCLE_POINTS = [(p["lat"], p["lon"]) for p in circle_data["points"]]

# ── Region definitions ──
REGIONS = {
    "Egypt": {"lat_range": (25, 31), "lon_range": (29, 34)},
    "Levant/Negev": {"lat_range": (29, 33), "lon_range": (34, 36)},
    "Mesopotamia/Iran": {"lat_range": (28, 38), "lon_range": (44, 60)},
    "South Asia/Indus": {"lat_range": (24, 30), "lon_range": (66, 73)},
    "Peru/Andes": {"lat_range": (-25, 5), "lon_range": (-82, -65)},  # broad box for SA
    "Easter Island": {"lat_range": (-28, -26), "lon_range": (-110, -108)},
    "SE Asia": {"lat_range": (0, 20), "lon_range": (95, 110)},
}

# ── Haversine distance ──
def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(a))

def min_distance_to_circle(lat, lon):
    """Minimum distance from point to any point on the great circle."""
    # Normalize longitude
    if lon > 180:
        lon -= 360
    elif lon < -180:
        lon += 360

    min_d = float('inf')
    for clat, clon in CIRCLE_POINTS:
        # Normalize circle longitude too
        norm_clon = clon
        if norm_clon > 180:
            norm_clon -= 360
        elif norm_clon < -180:
            norm_clon += 360
        d = haversine(lat, lon, clat, norm_clon)
        if d < min_d:
            min_d = d
    return min_d

def in_region(lat, lon, region_def):
    lat_min, lat_max = region_def["lat_range"]
    lon_min, lon_max = region_def["lon_range"]
    return lat_min <= lat <= lat_max and lon_min <= lon <= lon_max

# ── Date parsing ──
def parse_date_bce(date_str):
    """Parse various date formats. Returns years BCE (positive = BCE, negative = CE).
    Returns None if unparseable."""
    if date_str is None or str(date_str).strip() == '':
        return None
    s = str(date_str).strip()

    # Pleiades-style: negative = BCE, positive = CE
    try:
        val = float(s)
        if val <= 0:
            return abs(val)  # BCE
        else:
            return -val  # CE (negative in our system)
    except ValueError:
        pass

    # "3000 BCE" or "3000 BC"
    m = re.match(r'(\d+)\s*(BCE|BC)', s, re.IGNORECASE)
    if m:
        return int(m.group(1))

    # "500 CE" or "500 AD"
    m = re.match(r'(\d+)\s*(CE|AD)', s, re.IGNORECASE)
    if m:
        return -int(m.group(1))

    # "3rd millennium BCE" -> 3000 BCE
    m = re.match(r'(\d+)\w*\s*millennium\s*(BCE|BC)', s, re.IGNORECASE)
    if m:
        return int(m.group(1)) * 1000

    # "13th c CE" -> 1300 CE
    m = re.match(r'(\d+)\w*\s*c(?:entury)?\s*(CE|AD)', s, re.IGNORECASE)
    if m:
        return -(int(m.group(1)) * 100)

    # "13th c BCE"
    m = re.match(r'(\d+)\w*\s*c(?:entury)?\s*(BCE|BC)', s, re.IGNORECASE)
    if m:
        return int(m.group(1)) * 100

    # Range: "500 BCE-500 CE" -> take earliest
    m = re.match(r'(\d+)\s*(BCE|BC)', s, re.IGNORECASE)
    if m:
        return int(m.group(1))

    return None

def date_to_display(bce_val):
    """Convert internal BCE value to display string."""
    if bce_val is None:
        return "undated"
    if bce_val > 0:
        return f"{int(bce_val)} BCE"
    elif bce_val == 0:
        return "0 CE/BCE"
    else:
        return f"{int(abs(bce_val))} CE"

# ── Monumental classification ──
MONUMENTAL_KEYWORDS = {
    'pyramid', 'mastaba', 'temple', 'ziggurat', 'megalith', 'dolmen', 'menhir',
    'stone circle', 'henge', 'passage grave', 'chambered', 'cairn', 'barrow',
    'monumental', 'monument', 'palace', 'fortress', 'fort', 'citadel', 'wall',
    'ceremonial', 'ritual', 'tomb', 'mausoleum', 'stupa', 'geoglyph', 'moai',
    'ahu', 'platform', 'huaca', 'ashlar', 'cyclopean', 'sanctuary', 'shrine',
    'enclosure', 'tumulus', 'burial', 'necropolis', 'acropolis', 'standing stone',
    'cursus', 'earthwork', 'mound', 'artificial mound', 'long barrow', 'causeway',
    'broch', 'nuraghe', 'taula', 'naveta', 'hillfort', 'promontory fort',
    'stone row', 'alignment', 'avenue', 'rock cut tomb', 'rock art', 'carving',
}

SETTLEMENT_KEYWORDS = {
    'settlement', 'village', 'town', 'city', 'camp', 'domestic', 'habitation',
    'dwelling', 'house', 'farm', 'field', 'workshop', 'kiln', 'mine', 'quarry',
    'cave', 'rock shelter', 'shell midden', 'midden', 'pit', 'well', 'spring',
}

def classify_site(type_str, name_str=""):
    """Returns 'monumental', 'settlement', or 'other'."""
    combined = (str(type_str) + " " + str(name_str)).lower()
    for kw in MONUMENTAL_KEYWORDS:
        if kw in combined:
            return "monumental"
    for kw in SETTLEMENT_KEYWORDS:
        if kw in combined:
            return "settlement"
    return "other"

# ── Load databases ──
def load_p3k14c():
    """p3k14c radiocarbon database - 36K+ sites with C14 dates."""
    sites = []
    with open(BASE / "data/p3k14c/p3k14c_data.csv", encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row['Lat'])
                lon = float(row['Long'])
                age = float(row['Age'])  # C14 years BP
                # Convert C14 BP to calendar years BCE (rough: BP - 1950 = BCE)
                bce = age - 1950
                if bce < 0:
                    bce_val = -abs(bce)  # CE date
                else:
                    bce_val = bce  # BCE date
                sites.append({
                    'name': row.get('SiteName', ''),
                    'lat': lat, 'lon': lon,
                    'date_bce': bce_val,
                    'type': row.get('Period', ''),
                    'source': 'p3k14c',
                    'classification': classify_site(row.get('Period', ''), row.get('SiteName', ''))
                })
            except (ValueError, KeyError):
                continue
    return sites

def load_pleiades():
    """Pleiades gazetteer - classical world places with dates."""
    sites = []
    with open(BASE / "data/pleiades/pleiades-places-latest.csv", encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row['reprLat'])
                lon = float(row['reprLong'])
                min_date = parse_date_bce(row.get('minDate', ''))
                # Filter out geological/paleontological dates (>50000 BCE)
                if min_date is not None and min_date > 50000:
                    min_date = None
                feature_types = row.get('featureTypes', '')
                name = row.get('title', '')
                classification = classify_site(feature_types, name)
                sites.append({
                    'name': name,
                    'lat': lat, 'lon': lon,
                    'date_bce': min_date,
                    'type': feature_types,
                    'source': 'pleiades',
                    'classification': classification
                })
            except (ValueError, KeyError):
                continue
    return sites

def load_supplement():
    """Hand-curated supplement sites with dates."""
    sites = []
    with open(BASE / "github-repo/data/supplement_sites.json") as f:
        data = json.load(f)
    for s in data['sites']:
        date_bce = parse_date_bce(s.get('date', ''))
        classification = classify_site(s.get('type', ''), s.get('name', ''))
        sites.append({
            'name': s['name'],
            'lat': s['lat'], 'lon': s['lon'],
            'date_bce': date_bce,
            'type': s.get('type', ''),
            'source': 'supplement',
            'classification': classification
        })
    return sites

def load_xronos():
    """Xronos radiocarbon database."""
    sites = []
    with open(BASE / "data/xronos/xronos_sites.csv", encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row['lat'])
                lon = float(row['lon'])
                n_dates = int(row.get('n_dates', 0))
                site_type = row.get('site_type', '') or row.get('feature_type', '')
                name = row.get('site_name', '')
                classification = classify_site(site_type, name)
                # Xronos doesn't have aggregate dates in this export, skip undated
                sites.append({
                    'name': name,
                    'lat': lat, 'lon': lon,
                    'date_bce': None,  # No dates in this export
                    'type': site_type,
                    'source': 'xronos',
                    'classification': classification
                })
            except (ValueError, KeyError):
                continue
    return sites

def load_wikidata():
    """Wikidata archaeological sites."""
    sites = []
    with open(BASE / "data/wikidata/wikidata_archaeological_sites.csv", encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row['lat'])
                lon = float(row['lon'])
                inception = row.get('inception', '')
                date_bce = parse_date_bce(inception) if inception else None
                site_type = row.get('type', '')
                name = row.get('name', '')
                classification = classify_site(site_type, name)
                sites.append({
                    'name': name,
                    'lat': lat, 'lon': lon,
                    'date_bce': date_bce,
                    'type': site_type,
                    'source': 'wikidata',
                    'classification': classification
                })
            except (ValueError, KeyError):
                continue
    return sites

def load_peru():
    """Peru Ministry archaeological sites."""
    sites = []
    with open(BASE / "data/peru/peru_ministry_sites.csv", encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row['lat'])
                lon = float(row['lon'])
                name = row.get('name', '')
                site_type = row.get('type', '')
                classification_field = row.get('classification', '')
                classification = classify_site(site_type + ' ' + classification_field, name)
                if classification_field == 'monumental':
                    classification = 'monumental'
                sites.append({
                    'name': name,
                    'lat': lat, 'lon': lon,
                    'date_bce': None,  # No dates in this dataset
                    'type': site_type,
                    'source': 'peru_ministry',
                    'classification': classification
                })
            except (ValueError, KeyError):
                continue
    return sites

def load_south_asia():
    """South Asia sites."""
    sites = []
    with open(BASE / "data/south_asia/south_asia_sites.csv", encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row['lat'])
                lon = float(row['lon'])
                name = row.get('name', '')
                site_type = row.get('type', '')
                classification_field = row.get('classification', '')
                classification = classify_site(site_type + ' ' + classification_field, name)
                if classification_field == 'monumental':
                    classification = 'monumental'
                sites.append({
                    'name': name,
                    'lat': lat, 'lon': lon,
                    'date_bce': None,
                    'type': site_type,
                    'source': 'south_asia',
                    'classification': classification
                })
            except (ValueError, KeyError):
                continue
    return sites

def load_known_sites():
    """Well-established archaeological sites with literature-sourced dates.
    Fills gaps where databases lack dates for major sites."""
    known = [
        # Egypt - near circle
        {"name": "Nabta Playa stone circle", "lat": 22.51, "lon": 30.72, "date_bce": 5000, "type": "stone circle", "classification": "monumental"},
        {"name": "Step Pyramid of Djoser", "lat": 29.8713, "lon": 31.2164, "date_bce": 2700, "type": "pyramid", "classification": "monumental"},
        {"name": "Great Pyramid of Giza", "lat": 29.9792, "lon": 31.1342, "date_bce": 2560, "type": "pyramid", "classification": "monumental"},
        {"name": "Abydos (Umm el-Qa'ab)", "lat": 26.185, "lon": 31.92, "date_bce": 3100, "type": "royal tomb complex", "classification": "monumental"},

        # Levant
        {"name": "Tower of Jericho", "lat": 31.871, "lon": 35.444, "date_bce": 8000, "type": "stone tower", "classification": "monumental"},
        {"name": "Jericho PPNA settlement", "lat": 31.871, "lon": 35.444, "date_bce": 9000, "type": "settlement", "classification": "settlement"},
        {"name": "Beidha PPNB", "lat": 30.373, "lon": 35.439, "date_bce": 7200, "type": "settlement", "classification": "settlement"},

        # Mesopotamia/Iran
        {"name": "Eridu Temple I", "lat": 30.817, "lon": 45.996, "date_bce": 5400, "type": "temple", "classification": "monumental"},
        {"name": "Chogha Mish", "lat": 32.22, "lon": 48.55, "date_bce": 6800, "type": "settlement", "classification": "settlement"},
        {"name": "Susa Acropolis", "lat": 32.189, "lon": 48.246, "date_bce": 4200, "type": "temple platform", "classification": "monumental"},
        {"name": "Uruk Eanna precinct", "lat": 31.323, "lon": 45.637, "date_bce": 3500, "type": "temple complex", "classification": "monumental"},
        {"name": "Ur ziggurat", "lat": 30.963, "lon": 46.103, "date_bce": 2100, "type": "ziggurat", "classification": "monumental"},
        {"name": "Chogha Zanbil ziggurat", "lat": 32.009, "lon": 48.521, "date_bce": 1250, "type": "ziggurat", "classification": "monumental"},
        {"name": "Tepe Sialk", "lat": 33.968, "lon": 51.407, "date_bce": 5500, "type": "ziggurat platform", "classification": "monumental"},
        {"name": "Godin Tepe", "lat": 34.503, "lon": 48.068, "date_bce": 5000, "type": "fortified settlement", "classification": "monumental"},

        # South Asia/Indus
        {"name": "Mehrgarh (earliest)", "lat": 29.402, "lon": 67.599, "date_bce": 7000, "type": "settlement", "classification": "settlement"},
        {"name": "Mohenjo-daro Great Bath", "lat": 27.329, "lon": 68.139, "date_bce": 2500, "type": "monumental platform/bath", "classification": "monumental"},
        {"name": "Mohenjo-daro citadel", "lat": 27.329, "lon": 68.139, "date_bce": 2600, "type": "citadel", "classification": "monumental"},
        {"name": "Kot Diji fort", "lat": 27.345, "lon": 68.711, "date_bce": 3300, "type": "fortified settlement", "classification": "monumental"},
        {"name": "Amri settlement", "lat": 26.217, "lon": 68.168, "date_bce": 3600, "type": "settlement", "classification": "settlement"},
        {"name": "Dholavira citadel", "lat": 23.887, "lon": 70.211, "date_bce": 2650, "type": "citadel", "classification": "monumental"},
        {"name": "Kalibangan citadel", "lat": 29.473, "lon": 74.130, "date_bce": 2900, "type": "citadel", "classification": "monumental"},
        {"name": "Rehman Dheri fortification", "lat": 31.829, "lon": 70.397, "date_bce": 3200, "type": "fortified settlement", "classification": "monumental"},

        # Peru/Andes
        {"name": "Caral (Norte Chico)", "lat": -10.893, "lon": -77.520, "date_bce": 3000, "type": "pyramid complex", "classification": "monumental"},
        {"name": "Sechin Bajo", "lat": -9.467, "lon": -78.250, "date_bce": 3500, "type": "plaza/platform", "classification": "monumental"},
        {"name": "Huaca Prieta", "lat": -7.926, "lon": -79.315, "date_bce": 3100, "type": "mound", "classification": "monumental"},
        {"name": "Kotosh (Templo de las Manos Cruzadas)", "lat": -9.933, "lon": -76.267, "date_bce": 2300, "type": "temple", "classification": "monumental"},
        {"name": "Aspero", "lat": -10.847, "lon": -77.599, "date_bce": 3000, "type": "pyramid complex", "classification": "monumental"},
        {"name": "Chavín de Huántar", "lat": -9.593, "lon": -77.177, "date_bce": 1200, "type": "temple complex", "classification": "monumental"},
        {"name": "Nazca Lines", "lat": -14.735, "lon": -75.130, "date_bce": 500, "type": "geoglyph", "classification": "monumental"},

        # Easter Island
        {"name": "Ahu Tongariki", "lat": -27.126, "lon": -109.277, "date_bce": -1200, "type": "ahu/platform", "classification": "monumental"},
        {"name": "Ahu Akivi", "lat": -27.115, "lon": -109.395, "date_bce": -1460, "type": "ahu/platform", "classification": "monumental"},
        {"name": "Rano Raraku quarry", "lat": -27.122, "lon": -109.289, "date_bce": -1100, "type": "moai quarry", "classification": "monumental"},
        {"name": "Anakena earliest occupation", "lat": -27.074, "lon": -109.322, "date_bce": -800, "type": "settlement", "classification": "settlement"},

        # SE Asia
        {"name": "Ban Chiang (earliest)", "lat": 17.413, "lon": 102.832, "date_bce": 3600, "type": "settlement", "classification": "settlement"},
        {"name": "Non Nok Tha", "lat": 16.770, "lon": 102.360, "date_bce": 3000, "type": "settlement", "classification": "settlement"},
        {"name": "Khao Sam Kaeo", "lat": 9.250, "lon": 99.217, "date_bce": 400, "type": "fortified settlement", "classification": "monumental"},
        {"name": "Oc Eo", "lat": 10.360, "lon": 105.147, "date_bce": -100, "type": "port/temple", "classification": "monumental"},
        {"name": "Phnom Da", "lat": 10.823, "lon": 104.934, "date_bce": -500, "type": "temple", "classification": "monumental"},
        {"name": "Sambor Prei Kuk", "lat": 12.873, "lon": 105.028, "date_bce": -600, "type": "temple complex", "classification": "monumental"},
        {"name": "Plain of Jars", "lat": 19.428, "lon": 103.167, "date_bce": 500, "type": "megalithic jars", "classification": "monumental"},
    ]
    sites = []
    for s in known:
        sites.append({
            'name': s['name'],
            'lat': s['lat'],
            'lon': s['lon'],
            'date_bce': s['date_bce'],
            'type': s['type'],
            'source': 'literature',
            'classification': s['classification'],
        })
    return sites

def load_deep_time():
    """Deep-time expansion sites with dates."""
    sites = []
    fpath = BASE / "data/supplementary/deep-time-expansion-sites.json"
    if fpath.exists():
        with open(fpath) as f:
            data = json.load(f)
        for s in data if isinstance(data, list) else data.get('sites', []):
            try:
                date_bce = parse_date_bce(s.get('date', s.get('date_bce', '')))
                classification = classify_site(s.get('type', ''), s.get('name', ''))
                sites.append({
                    'name': s.get('name', ''),
                    'lat': float(s['lat']),
                    'lon': float(s['lon']),
                    'date_bce': date_bce,
                    'type': s.get('type', ''),
                    'source': 'deep_time',
                    'classification': classification
                })
            except (ValueError, KeyError):
                continue
    return sites

def load_unesco():
    """UNESCO cultural heritage sites."""
    sites = []
    with open(BASE / "data/unesco/unesco_cultural_sites.json") as f:
        data = json.load(f)
    site_list = data if isinstance(data, list) else data.get('sites', data.get('features', []))
    for s in site_list:
        try:
            if isinstance(s, dict):
                # Handle different JSON structures
                props = s.get('properties', s)
                lat = float(props.get('lat', props.get('latitude', 0)))
                lon = float(props.get('lon', props.get('longitude', 0)))
                if lat == 0 and lon == 0:
                    continue
                name = props.get('name', props.get('title', ''))
                date_bce = parse_date_bce(props.get('date', props.get('inscription_date', '')))
                site_type = props.get('type', props.get('category', ''))
                classification = classify_site(site_type, name)
                sites.append({
                    'name': name,
                    'lat': lat, 'lon': lon,
                    'date_bce': date_bce,
                    'type': site_type,
                    'source': 'unesco',
                    'classification': classification
                })
        except (ValueError, KeyError, TypeError):
            continue
    return sites

# ── Known site dates from archaeological knowledge ──
# For sites we identify but that lack dates in the databases,
# we can fill in well-established dates from the literature
KNOWN_DATES = {
    # Egypt
    'giza': 4550,       # Great Pyramid ~2600 BCE -> 4550 BP... no, let's use BCE
    'saqqara': 4700,    # Step Pyramid ~2700 BCE
    'abydos': 5200,     # Early Dynastic ~3200 BCE
    'hierakonpolis': 5500, # Pre-dynastic ~3500 BCE
    'nabta playa': 7000,  # stone circle ~5000 BCE

    # Levant
    'jericho': 10000,   # Tower of Jericho ~8000 BCE
    'göbekli': 11600,   # ~9600 BCE  (outside region but reference)

    # Mesopotamia
    'eridu': 6500,      # ~4500 BCE
    'uruk': 5500,       # ~3500 BCE
    'ur': 5000,         # ~3000 BCE
    'susa': 6000,       # ~4000 BCE
    'chogha zanbil': 3250, # ~1250 BCE
    'sialk': 7500,      # ~5500 BCE
    'godin tepe': 5000, # ~3000 BCE

    # South Asia
    'mohenjo-daro': 4600, # ~2600 BCE
    'harappa': 5300,    # ~3300 BCE
    'dholavira': 4700,  # ~2700 BCE
    'kalibangan': 4900, # ~2900 BCE
    'mehrgarh': 9000,   # ~7000 BCE
    'rakhigarhi': 4600, # ~2600 BCE

    # Peru
    'caral': 5000,      # ~3000 BCE
    'sechin bajo': 5500, # ~3500 BCE
    'huaca prieta': 7500, # ~5500 BCE
    'kotosh': 4000,     # ~2000 BCE
    'chavín': 3200,     # ~1200 BCE
    'nazca': 2500,      # ~500 BCE

    # Easter Island
    'ahu tongariki': 700, # ~1250 CE -> -1250
    'rano raraku': 750,   # ~1200 CE

    # SE Asia
    'ban chiang': 5500, # ~3500 BCE
    'khao sam kaeo': 2400, # ~400 BCE
    'non nok tha': 5000, # ~3000 BCE
    'ban non wat': 4000, # ~2000 BCE
    'phu lon': 4000,    # ~2000 BCE
}

# ── Main analysis ──
def main():
    print("=" * 90)
    print("REGIONAL ONSET ANALYSIS: Earliest Construction Within 200km of Great Circle")
    print("=" * 90)
    print(f"\nGreat Circle: pole at ({POLE[0]:.2f}°N, {POLE[1]:.2f}°E)")
    print(f"Analysis radius: 200 km from circle\n")

    # Load all databases
    print("Loading databases...")
    all_sites = []

    loaders = [
        ("p3k14c", load_p3k14c),
        ("Pleiades", load_pleiades),
        ("Supplement", load_supplement),
        ("Wikidata", load_wikidata),
        ("Literature (known dates)", load_known_sites),
        ("Deep-time", load_deep_time),
        ("UNESCO", load_unesco),
        ("Peru Ministry", load_peru),
        ("South Asia", load_south_asia),
    ]

    for name, loader in loaders:
        try:
            sites = loader()
            all_sites.extend(sites)
            print(f"  {name}: {len(sites):,} sites")
        except Exception as e:
            print(f"  {name}: FAILED - {e}")

    print(f"\nTotal sites loaded: {len(all_sites):,}")

    # For each region, find sites within 200km of the circle
    print("\n" + "=" * 90)

    results = {}

    for region_name, region_def in REGIONS.items():
        print(f"\n{'─' * 90}")
        print(f"  REGION: {region_name}")
        print(f"  Bounding box: {region_def['lat_range'][0]}-{region_def['lat_range'][1]}°N, "
              f"{region_def['lon_range'][0]}-{region_def['lon_range'][1]}°E")
        print(f"{'─' * 90}")

        # Filter sites in region
        regional_sites = [s for s in all_sites if in_region(s['lat'], s['lon'], region_def)]
        print(f"  Sites in region: {len(regional_sites):,}")

        # Compute distance to circle for each site
        within_200km = []
        for s in regional_sites:
            d = min_distance_to_circle(s['lat'], s['lon'])
            s_copy = dict(s)
            s_copy['gc_distance'] = d
            if d <= 200:
                within_200km.append(s_copy)

        print(f"  Sites within 200km of circle: {len(within_200km):,}")

        # Separate by classification
        monumental = [s for s in within_200km if s['classification'] == 'monumental']
        settlements = [s for s in within_200km if s['classification'] == 'settlement']
        dated_monumental = [s for s in monumental if s['date_bce'] is not None]
        dated_settlements = [s for s in settlements if s['date_bce'] is not None]
        dated_any = [s for s in within_200km if s['date_bce'] is not None]

        print(f"  Monumental: {len(monumental)} ({len(dated_monumental)} dated)")
        print(f"  Settlement: {len(settlements)} ({len(dated_settlements)} dated)")
        print(f"  Other/unclassified: {len(within_200km) - len(monumental) - len(settlements)}")

        # Find earliest monumental
        if dated_monumental:
            earliest_mon = max(dated_monumental, key=lambda s: s['date_bce'])
            print(f"\n  ★ EARLIEST MONUMENTAL (within 200km):")
            print(f"    {earliest_mon['name']}")
            print(f"    Type: {earliest_mon['type']} | Source: {earliest_mon['source']}")
            print(f"    Date: {date_to_display(earliest_mon['date_bce'])}")
            print(f"    Distance to circle: {earliest_mon['gc_distance']:.1f} km")
            print(f"    Location: ({earliest_mon['lat']:.4f}, {earliest_mon['lon']:.4f})")
        else:
            print(f"\n  ★ EARLIEST MONUMENTAL: No dated monumental sites found")
            earliest_mon = None
            # Show undated monumental sites
            if monumental:
                print(f"    (Undated monumental sites exist: {', '.join(s['name'] for s in monumental[:5])})")

        # Find earliest of any type
        if dated_any:
            earliest_any = max(dated_any, key=lambda s: s['date_bce'])
            print(f"\n  ★ EARLIEST ANY TYPE (within 200km):")
            print(f"    {earliest_any['name']}")
            print(f"    Type: {earliest_any['type']} ({earliest_any['classification']}) | Source: {earliest_any['source']}")
            print(f"    Date: {date_to_display(earliest_any['date_bce'])}")
            print(f"    Distance to circle: {earliest_any['gc_distance']:.1f} km")
            print(f"    Location: ({earliest_any['lat']:.4f}, {earliest_any['lon']:.4f})")
        else:
            earliest_any = None
            print(f"\n  ★ EARLIEST ANY TYPE: No dated sites found")

        # Top 5 earliest monumental
        if len(dated_monumental) > 1:
            print(f"\n  Top 5 earliest monumental sites:")
            for i, s in enumerate(sorted(dated_monumental, key=lambda s: -s['date_bce'])[:5]):
                print(f"    {i+1}. {s['name']} — {date_to_display(s['date_bce'])} — "
                      f"{s['gc_distance']:.0f}km — {s['source']}")

        # Top 5 earliest any type
        if len(dated_any) > 1:
            print(f"\n  Top 5 earliest sites (any type):")
            for i, s in enumerate(sorted(dated_any, key=lambda s: -s['date_bce'])[:5]):
                print(f"    {i+1}. {s['name']} — {date_to_display(s['date_bce'])} — "
                      f"{s['type']} ({s['classification']}) — {s['gc_distance']:.0f}km — {s['source']}")

        # Monument density vs settlement density over time
        # Bin into 1000-year intervals
        print(f"\n  Temporal density (dated sites, 1000-year bins):")
        bins = defaultdict(lambda: {'monumental': 0, 'settlement': 0, 'other': 0})
        for s in dated_any:
            bin_val = int(s['date_bce'] // 1000) * 1000
            bins[bin_val][s['classification']] += 1

        onset_date = None
        for bin_val in sorted(bins.keys(), reverse=True):
            counts = bins[bin_val]
            m, s_count = counts['monumental'], counts['settlement']
            marker = ""
            if m > s_count and m > 0 and onset_date is None:
                onset_date = bin_val
                marker = " ← ONSET (monuments > settlements)"
            print(f"    {date_to_display(bin_val)}-{date_to_display(bin_val-999)}: "
                  f"M={m} S={s_count} O={counts['other']}{marker}")

        if onset_date is not None:
            print(f"\n  Monument > Settlement onset: ~{date_to_display(onset_date)}")
        else:
            print(f"\n  Monument > Settlement onset: Never observed in dated sample")

        results[region_name] = {
            'total_in_region': len(regional_sites),
            'within_200km': len(within_200km),
            'monumental': len(monumental),
            'settlements': len(settlements),
            'earliest_monumental': {
                'name': earliest_mon['name'] if earliest_mon else None,
                'date': date_to_display(earliest_mon['date_bce']) if earliest_mon else None,
                'distance_km': earliest_mon['gc_distance'] if earliest_mon else None,
            } if earliest_mon else None,
            'earliest_any': {
                'name': earliest_any['name'] if earliest_any else None,
                'date': date_to_display(earliest_any['date_bce']) if earliest_any else None,
                'distance_km': earliest_any['gc_distance'] if earliest_any else None,
            } if earliest_any else None,
            'monument_onset': date_to_display(onset_date) if onset_date else None,
        }

    # ── Summary comparison ──
    print(f"\n\n{'=' * 90}")
    print("CROSS-REGIONAL COMPARISON: Diffusion vs Independent Invention")
    print("=" * 90)
    print(f"\n{'Region':<20} {'Earliest Monument':<25} {'Earliest Any':<25} {'Mon>Set Onset':<20} {'Dist(km)':<10}")
    print("─" * 100)

    for region_name in REGIONS:
        r = results[region_name]
        mon = r['earliest_monumental']
        any_s = r['earliest_any']
        mon_date = mon['date'] if mon else 'N/A'
        any_date = any_s['date'] if any_s else 'N/A'
        mon_onset = r['monument_onset'] or 'Never'
        mon_dist = f"{mon['distance_km']:.0f}" if mon else 'N/A'
        print(f"{region_name:<20} {mon_date:<25} {any_date:<25} {mon_onset:<20} {mon_dist:<10}")

    print(f"\n{'─' * 100}")
    print("\nINTERPRETATION:")

    # Collect onset dates for comparison
    onset_dates = {}
    for region_name in REGIONS:
        r = results[region_name]
        if r['earliest_monumental']:
            date_str = r['earliest_monumental']['date']
            # Parse back to numeric for comparison
            m = re.match(r'(\d+)\s*BCE', date_str)
            if m:
                onset_dates[region_name] = int(m.group(1))
            m2 = re.match(r'(\d+)\s*CE', date_str)
            if m2:
                onset_dates[region_name] = -int(m2.group(1))

    if len(onset_dates) >= 2:
        values = list(onset_dates.values())
        spread = max(values) - min(values)
        avg = sum(values) / len(values)

        print(f"  Earliest onset dates across regions: {onset_dates}")
        print(f"  Spread: {spread} years")
        print(f"  Mean onset: ~{int(avg)} {'BCE' if avg > 0 else 'CE'}")

        if spread < 1000:
            print(f"  → SUGGESTIVE OF DIFFUSION: All regions onset within ~{spread} years")
        elif spread < 3000:
            print(f"  → AMBIGUOUS: {spread}-year spread could be either diffusion or independent")
        else:
            print(f"  → SUGGESTIVE OF INDEPENDENT DEVELOPMENT: {spread}-year spread across regions")

    # Save results
    output_path = BASE / "results/regional_onset_analysis.json"
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {output_path}")

if __name__ == '__main__':
    main()
