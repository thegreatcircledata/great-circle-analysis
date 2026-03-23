#!/usr/bin/env python3
"""
Directive T1-7: Linguistic/Cultural Connectivity
=================================================
Using Glottolog language family data, test whether languages along the
Great Circle corridor show more family-sharing than off-corridor languages.
"""

import json, math, os, sys, csv
import numpy as np
from collections import Counter, defaultdict

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "linguistic_connectivity")
os.makedirs(OUT_DIR, exist_ok=True)

POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
CORRIDOR_KM = 500

def haversine_km(lat1, lon1, lat2, lon2):
    R = EARTH_R_KM
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(min(1.0, a)))

def dist_from_gc(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QUARTER_CIRC)

# ============================================================
# DOWNLOAD GLOTTOLOG IF NEEDED
# ============================================================
glottolog_dir = os.path.join(BASE_DIR, "data", "glottolog")
languages_csv = os.path.join(glottolog_dir, "languages.csv")

if not os.path.exists(languages_csv):
    print("Glottolog data not found. Attempting download...")
    os.makedirs(glottolog_dir, exist_ok=True)
    import urllib.request, zipfile
    url = "https://glottolog.org/glottolog/glottologinformation?type=languages&sEcho=1"
    # Try the direct CSV from Glottolog's cldf release
    cldf_url = "https://raw.githubusercontent.com/glottolog/glottolog-cldf/master/cldf/languages.csv"
    try:
        urllib.request.urlretrieve(cldf_url, languages_csv)
        print("Downloaded Glottolog languages.csv")
    except Exception as e:
        print(f"Could not download from GitHub: {e}")
        # Fallback: use languoids.csv from Glottolog
        languoid_url = "https://cdstar.eva.mpg.de/bitstreams/EAEA0-E7DE-FA06-8817-0/glottolog_languoid.csv.zip"
        try:
            zip_path = os.path.join(glottolog_dir, "glottolog_languoid.csv.zip")
            urllib.request.urlretrieve(languoid_url, zip_path)
            with zipfile.ZipFile(zip_path, 'r') as zf:
                zf.extractall(glottolog_dir)
            print("Downloaded and extracted Glottolog languoid data")
        except Exception as e2:
            print(f"Download failed: {e2}")
            print("Will use built-in language family data instead.")

# ============================================================
# LOAD LANGUAGE DATA
# ============================================================
print("\nLoading language data...")
languages = []

# Try Glottolog CLDF format first
if os.path.exists(languages_csv):
    # First pass: build ID->Name map for family resolution
    # In Glottolog CLDF, Family_ID references the ID of the top-level family row
    id_to_name = {}
    family_ids_with_coords = {}
    with open(languages_csv, encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rid = row.get('ID', '')
            rname = row.get('Name', '')
            rlevel = row.get('Level', '')
            if rid:
                id_to_name[rid] = rname
            if rlevel == 'family':
                family_ids_with_coords[rid] = rname

    print(f"Built ID->Name map: {len(id_to_name)} entries, {len(family_ids_with_coords)} families")

    # Second pass: load languages with coordinates
    with open(languages_csv, encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            level = row.get('Level', '')
            if level != 'language':
                continue  # Only individual languages, not families/subfamilies

            try:
                lat = float(row.get('Latitude', ''))
                lon = float(row.get('Longitude', ''))
            except (ValueError, TypeError):
                continue
            if lat == 0 and lon == 0:
                continue

            name = row.get('Name', '')
            family_id = row.get('Family_ID', '')
            is_isolate = row.get('Is_Isolate', '')
            macroarea = row.get('Macroarea', '')
            lang_id = row.get('ID', '')

            # Resolve family name from Family_ID
            if family_id and family_id in id_to_name:
                family = id_to_name[family_id]
            elif is_isolate and is_isolate.lower() == 'true':
                family = f'Isolate:{name}'
            else:
                family = 'Unclassified'

            dist = dist_from_gc(lat, lon)
            languages.append({
                'name': name,
                'lat': lat,
                'lon': lon,
                'family': family,
                'subfamily': '',
                'macroarea': macroarea,
                'id': lang_id,
                'dist_from_gc': dist,
                'on_corridor': dist < CORRIDOR_KM
            })
else:
    # Use a comprehensive built-in dataset of major language families with representative coordinates
    print("Using built-in language family reference data...")
    # This is a curated subset — actual analysis would use full Glottolog download
    BUILT_IN_FAMILIES = [
        # Format: (name, lat, lon, family, subfamily)
        # Along the corridor (Nile-Levant-Mesopotamia-Iran-India-SE Asia-Pacific)
        ("Egyptian Arabic", 30.0, 31.2, "Afro-Asiatic", "Semitic"),
        ("Coptic", 28.0, 31.0, "Afro-Asiatic", "Egyptian"),
        ("Hebrew", 31.5, 35.0, "Afro-Asiatic", "Semitic"),
        ("Aramaic", 37.0, 43.0, "Afro-Asiatic", "Semitic"),
        ("Akkadian", 33.0, 44.0, "Afro-Asiatic", "Semitic"),
        ("Sumerian", 31.3, 46.0, "Sumerian", "Isolate"),
        ("Persian", 32.0, 52.0, "Indo-European", "Iranian"),
        ("Kurdish", 36.0, 44.0, "Indo-European", "Iranian"),
        ("Elamite", 32.0, 49.0, "Elamite", "Isolate"),
        ("Hindi", 28.0, 77.0, "Indo-European", "Indo-Aryan"),
        ("Sanskrit", 27.0, 82.0, "Indo-European", "Indo-Aryan"),
        ("Tamil", 11.0, 79.0, "Dravidian", "South Dravidian"),
        ("Brahui", 29.0, 66.0, "Dravidian", "North Dravidian"),
        ("Thai", 14.0, 100.5, "Kra-Dai", "Tai"),
        ("Khmer", 12.0, 105.0, "Austroasiatic", "Khmeric"),
        ("Malay", 3.0, 102.0, "Austronesian", "Malayo-Polynesian"),
        ("Javanese", -7.5, 110.0, "Austronesian", "Malayo-Polynesian"),
        ("Rapa Nui", -27.1, -109.4, "Austronesian", "Polynesian"),
        ("Quechua", -13.5, -72.0, "Quechuan", "Quechua II"),
        ("Aymara", -16.5, -68.0, "Aymaran", "Aymaran"),
        ("Nahuatl", 19.0, -99.0, "Uto-Aztecan", "Aztecan"),
        ("Maya", 17.0, -89.0, "Mayan", "Yucatecan"),
        # Off corridor examples
        ("Finnish", 61.0, 25.0, "Uralic", "Finnic"),
        ("Hungarian", 47.0, 19.0, "Uralic", "Ugric"),
        ("Mandarin", 40.0, 116.0, "Sino-Tibetan", "Sinitic"),
        ("Cantonese", 23.0, 113.0, "Sino-Tibetan", "Sinitic"),
        ("Tibetan", 30.0, 91.0, "Sino-Tibetan", "Tibetic"),
        ("Japanese", 36.0, 140.0, "Japonic", "Japanese"),
        ("Korean", 37.5, 127.0, "Koreanic", "Korean"),
        ("Yoruba", 8.0, 4.0, "Niger-Congo", "Volta-Niger"),
        ("Swahili", -6.0, 37.0, "Niger-Congo", "Bantu"),
        ("Zulu", -29.0, 31.0, "Niger-Congo", "Bantu"),
        ("Berber", 33.0, 2.0, "Afro-Asiatic", "Berber"),
        ("Hausa", 12.0, 8.0, "Afro-Asiatic", "Chadic"),
        ("English", 52.0, -1.0, "Indo-European", "Germanic"),
        ("German", 51.0, 10.0, "Indo-European", "Germanic"),
        ("French", 47.0, 2.0, "Indo-European", "Romance"),
        ("Spanish", 40.0, -4.0, "Indo-European", "Romance"),
        ("Greek", 38.0, 24.0, "Indo-European", "Hellenic"),
        ("Russian", 56.0, 38.0, "Indo-European", "Slavic"),
        ("Turkish", 39.0, 32.0, "Turkic", "Oghuz"),
        ("Mongolian", 48.0, 107.0, "Mongolic", "Mongolic"),
        ("Navajo", 36.0, -109.0, "Na-Dene", "Athabaskan"),
        ("Inuktitut", 64.0, -68.0, "Eskimo-Aleut", "Inuit"),
    ]

    for name, lat, lon, family, subfamily in BUILT_IN_FAMILIES:
        dist = dist_from_gc(lat, lon)
        languages.append({
            'name': name,
            'lat': lat,
            'lon': lon,
            'family': family,
            'subfamily': subfamily,
            'macroarea': '',
            'id': '',
            'dist_from_gc': dist,
            'on_corridor': dist < CORRIDOR_KM
        })

print(f"Total languages loaded: {len(languages)}")
on_corridor = [l for l in languages if l['on_corridor']]
off_corridor = [l for l in languages if not l['on_corridor']]
print(f"On corridor (<{CORRIDOR_KM}km): {len(on_corridor)}")
print(f"Off corridor: {len(off_corridor)}")

# ============================================================
# ANALYSIS
# ============================================================
print(f"\n{'='*70}")
print("LINGUISTIC CONNECTIVITY ANALYSIS")
print(f"{'='*70}")

# Test 1: Family sharing rate
def compute_family_sharing(lang_list, max_pairs=50000):
    """For pairs of languages, what fraction share a family?"""
    n = len(lang_list)
    if n < 2:
        return 0, 0, []

    same_family = 0
    total_pairs = 0
    pair_distances = []

    for i in range(n):
        for j in range(i+1, n):
            d = haversine_km(lang_list[i]['lat'], lang_list[i]['lon'],
                           lang_list[j]['lat'], lang_list[j]['lon'])
            shares = lang_list[i]['family'] == lang_list[j]['family']
            pair_distances.append((d, shares))
            if shares:
                same_family += 1
            total_pairs += 1
            if total_pairs >= max_pairs:
                break
        if total_pairs >= max_pairs:
            break

    rate = same_family / total_pairs if total_pairs > 0 else 0
    return rate, total_pairs, pair_distances

on_rate, on_pairs, on_pair_data = compute_family_sharing(on_corridor)
off_rate, off_pairs, off_pair_data = compute_family_sharing(off_corridor)

print(f"\n1. FAMILY SHARING RATE")
print(f"   On corridor: {on_rate*100:.1f}% ({on_pairs} pairs)")
print(f"   Off corridor: {off_rate*100:.1f}% ({off_pairs} pairs)")

# Distance-matched comparison
def sharing_at_distance(pair_data, d_min, d_max):
    relevant = [(d, s) for d, s in pair_data if d_min <= d <= d_max]
    if not relevant:
        return None, 0
    shares = sum(1 for _, s in relevant if s)
    return shares / len(relevant), len(relevant)

print(f"\n   Distance-matched family sharing:")
for d_min, d_max in [(0, 1000), (1000, 3000), (3000, 5000), (5000, 10000)]:
    on_r, on_n = sharing_at_distance(on_pair_data, d_min, d_max)
    off_r, off_n = sharing_at_distance(off_pair_data, d_min, d_max)
    on_str = f"{on_r*100:.1f}% (n={on_n})" if on_r is not None else "no data"
    off_str = f"{off_r*100:.1f}% (n={off_n})" if off_r is not None else "no data"
    print(f"   {d_min}-{d_max}km: on-corridor={on_str}, off-corridor={off_str}")

# Test 2: Language diversity (families per region)
on_families = set(l['family'] for l in on_corridor)
off_families = set(l['family'] for l in off_corridor)

print(f"\n2. LANGUAGE DIVERSITY")
print(f"   On corridor: {len(on_families)} families for {len(on_corridor)} languages")
print(f"   Off corridor: {len(off_families)} families for {len(off_corridor)} languages")
print(f"   On corridor families: {sorted(on_families)}")

# Diversity per degree of arc (rough)
on_diversity = len(on_families) / max(1, len(on_corridor)) * 100
off_diversity = len(off_families) / max(1, len(off_corridor)) * 100
print(f"   Diversity index (families/100 languages): on={on_diversity:.1f}, off={off_diversity:.1f}")

# Test 3: Cross-family contact (different families on same corridor)
cross_family_on = set()
for i in range(len(on_corridor)):
    for j in range(i+1, len(on_corridor)):
        if on_corridor[i]['family'] != on_corridor[j]['family']:
            d = haversine_km(on_corridor[i]['lat'], on_corridor[i]['lon'],
                           on_corridor[j]['lat'], on_corridor[j]['lon'])
            if d < 1000:  # nearby different families = contact zone
                cross_family_on.add((on_corridor[i]['family'], on_corridor[j]['family']))

print(f"\n3. CROSS-FAMILY CONTACT ZONES (within 1000km)")
print(f"   On corridor: {len(cross_family_on)} unique family contact pairs")
for pair in sorted(cross_family_on):
    print(f"     {pair[0]} ↔ {pair[1]}")

# ============================================================
# RESULTS
# ============================================================
results = {
    'corridor_km': CORRIDOR_KM,
    'n_on_corridor': len(on_corridor),
    'n_off_corridor': len(off_corridor),
    'family_sharing': {
        'on_corridor_rate': on_rate,
        'off_corridor_rate': off_rate,
        'on_pairs': on_pairs,
        'off_pairs': off_pairs
    },
    'diversity': {
        'on_corridor_families': len(on_families),
        'off_corridor_families': len(off_families),
        'on_corridor_family_list': sorted(on_families),
        'off_corridor_family_list': sorted(off_families),
        'on_diversity_index': on_diversity,
        'off_diversity_index': off_diversity
    },
    'cross_family_contacts': {
        'on_corridor_pairs': [list(p) for p in sorted(cross_family_on)]
    },
    'interpretation': (
        'Higher on-corridor diversity suggests the corridor crosses multiple language zones, '
        'consistent with a geographic/cultural boundary rather than a homogenizing route.'
        if on_diversity > off_diversity else
        'Lower on-corridor diversity could suggest a homogenizing cultural route.'
    )
}

with open(os.path.join(OUT_DIR, "family_sharing.json"), 'w') as f:
    json.dump(results, f, indent=2)

with open(os.path.join(OUT_DIR, "diversity_comparison.json"), 'w') as f:
    json.dump({
        'on_corridor': [{'name': l['name'], 'family': l['family'], 'lat': l['lat'], 'lon': l['lon'],
                         'dist_from_gc': l['dist_from_gc']} for l in on_corridor],
        'off_corridor': [{'name': l['name'], 'family': l['family'], 'lat': l['lat'], 'lon': l['lon'],
                         'dist_from_gc': l['dist_from_gc']} for l in off_corridor]
    }, f, indent=2)

print(f"\nDone! Outputs saved to {OUT_DIR}")
