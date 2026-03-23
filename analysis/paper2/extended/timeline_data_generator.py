#!/usr/bin/env python3
"""
Timeline Visualization Data Generator
======================================
Compiles ALL dated sites within 100km of the Great Circle from p3k14c,
XRONOS, and Pleiades. Bins into 9 time periods from 25,000 BP to 300 BP.

Outputs:
  - timeline_data.json        — all sites by period with enrichment
  - timeline_geojson.json     — GeoJSON FeatureCollection for web mapping
  - timeline_summary.json     — counts and enrichment by period
  - timeline_static.png       — multi-panel static image (one panel per period)
"""

import csv, math, json, os, sys
import numpy as np
from collections import defaultdict

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
THRESHOLD_KM = 100  # 100km corridor for timeline
N_TRIALS = 200
SEED = 42

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "extended_analysis")
os.makedirs(OUT_DIR, exist_ok=True)

# Time periods (start_bp, end_bp, label)
# BP = Before Present (1950)
PERIODS = [
    (25000, 12800, "Pre-Younger Dryas (25,000–12,800 BP)"),
    (12800, 11600, "Younger Dryas (12,800–11,600 BP)"),
    (11600, 10000, "Early Holocene (11,600–10,000 BP)"),
    (10000,  8000, "PPNB Era (10,000–8,000 BP)"),
    ( 8000,  5000, "Neolithic/Chalcolithic (8,000–5,000 BP)"),
    ( 5000,  3000, "Bronze Age (5,000–3,000 BP)"),
    ( 3000,  2000, "Iron Age (3,000–2,000 BP)"),
    ( 2000,   800, "Classical/Medieval (2,000–800 BP)"),
    (  800,   300, "Late (800–300 BP)"),
]

# Site type classification keywords
MONUMENTAL_KW = [
    'temple', 'pyramid', 'tomb', 'mound', 'monument', 'cairn', 'megalith',
    'henge', 'huaca', 'ceremonial', 'ritual', 'geoglyph', 'ahu', 'moai',
    'stela', 'ziggurat', 'palace', 'citadel', 'fortress', 'dolmen',
    'stone circle', 'standing stone', 'barrow', 'tumulus', 'sanctuary',
    'necropolis', 'acropolis', 'amphitheatre', 'colosseum', 'cathedral',
    'church', 'mosque', 'shrine', 'pagoda', 'stupa', 'nuraghe',
    'passage grave', 'gallery grave', 'earthwork', 'petroglyph',
    'rock art', 'cave art',
]

DOMESTIC_KW = [
    'village', 'settlement', 'farm', 'camp', 'midden', 'shell mound',
    'quarry', 'mine', 'shelter', 'rockshelter', 'workshop', 'kiln', 'pit',
    'habitation', 'dwelling', 'house', 'homestead', 'terrace', 'field',
    'granary', 'storehouse', 'well', 'cistern', 'fishweir', 'port',
]

PLEIADES_MONUMENTAL = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe",
    "tumulus", "shrine",
}

PLEIADES_SETTLEMENT = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement",
    "townhouse", "production",
}

# ============================================================
# VECTORIZED HELPERS
# ============================================================
def gc_dist_vec(pole_lat, pole_lon, site_lats, site_lons):
    """Distance from each site to the Great Circle (vectorized)."""
    lat1_r, lon1_r = np.radians(pole_lat), np.radians(pole_lon)
    lat2_r, lon2_r = np.radians(site_lats), np.radians(site_lons)
    dlat = lat2_r - lat1_r
    dlon = lon2_r - lon1_r
    a = np.sin(dlat/2)**2 + np.cos(lat1_r)*np.cos(lat2_r)*np.sin(dlon/2)**2
    d = EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QUARTER_CIRC)


def dist_from_circle(lat, lon):
    """Scalar version for single site."""
    lat1 = math.radians(POLE_LAT)
    lon1 = math.radians(POLE_LON)
    lat2 = math.radians(lat)
    lon2 = math.radians(lon)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    d = EARTH_R_KM * 2 * math.asin(min(1.0, math.sqrt(a)))
    return abs(d - QUARTER_CIRC)


def bearing_from_pole(lat, lon):
    """Azimuthal position along the circle (0-360 degrees from pole)."""
    lat1 = math.radians(POLE_LAT)
    lon1 = math.radians(POLE_LON)
    lat2 = math.radians(lat)
    lon2 = math.radians(lon)
    dlon = lon2 - lon1
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(dlon)
    return math.degrees(math.atan2(x, y)) % 360


def c14_to_cal_bp(age_bp):
    """Convert C14 years BP to approximate calendar years BP (IntCal20 approx)."""
    if age_bp <= 0:
        return 0
    if age_bp < 2500:
        return age_bp * 1.0
    elif age_bp < 5000:
        return age_bp * 1.05 + 50
    elif age_bp < 8000:
        return age_bp * 1.08 + 100
    elif age_bp < 12000:
        return age_bp * 1.12 + 200
    else:
        return age_bp * 1.15 + 300


def cal_year_to_bp(cal_year):
    """Convert calendar year (negative=BCE) to BP."""
    return 1950 - cal_year


def bp_to_cal_year(bp):
    """Convert BP to calendar year (negative=BCE)."""
    return 1950 - bp


def classify_name(name):
    """Classify a site by name keywords."""
    if not name:
        return "unknown"
    name_lower = name.lower()
    for kw in MONUMENTAL_KW:
        if kw in name_lower:
            return "monumental"
    for kw in DOMESTIC_KW:
        if kw in name_lower:
            return "settlement"
    return "unknown"


def rand_matched_batch(site_lats, site_lons, n):
    """Generate distribution-matched random sites."""
    idx = np.random.randint(0, len(site_lats), size=n)
    lats = site_lats[idx] + np.random.normal(0, 2, n)
    lons = site_lons[idx] + np.random.normal(0, 2, n)
    return np.clip(lats, -90, 90), np.clip(lons, -180, 180)


def compute_enrichment(near_sites, all_lats, all_lons, threshold_km, n_trials):
    """Compute enrichment ratio: observed / expected near-circle count."""
    n_total = len(all_lats)
    observed = len(near_sites)
    if n_total < 10:
        return {"observed": observed, "expected": 0, "enrichment": float('nan'),
                "z_score": float('nan'), "too_few": True}

    rand_counts = []
    for _ in range(n_trials):
        rl, rn = rand_matched_batch(all_lats, all_lons, n_total)
        rd = gc_dist_vec(POLE_LAT, POLE_LON, rl, rn)
        rand_counts.append(int(np.sum(rd <= threshold_km)))

    mu = np.mean(rand_counts)
    sigma = np.std(rand_counts)
    z = float((observed - mu) / sigma) if sigma > 0 else 0.0
    enrich = float(observed / mu) if mu > 0 else 0.0

    return {
        "observed": observed,
        "expected": round(float(mu), 2),
        "std": round(float(sigma), 2),
        "enrichment": round(enrich, 3),
        "z_score": round(z, 2),
        "too_few": False,
    }


# ============================================================
# LOAD ALL THREE DATABASES
# ============================================================
print("=" * 70)
print("TIMELINE VISUALIZATION DATA GENERATOR")
print("=" * 70)

all_sites = []  # master list: {lat, lon, date_bp, name, type, source, dist_km, arc_deg}

# --- 1. p3k14c ---
print("\n[1/3] Loading p3k14c...")
p3k_file = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")
site_groups = defaultdict(list)
p3k_raw = 0

with open(p3k_file, encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            lat = float(row['Lat'])
            lon = float(row['Long'])
            age_bp = float(row['Age']) if row['Age'] else None
        except (ValueError, TypeError):
            continue
        if lat == 0 and lon == 0:
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue
        if age_bp is None or age_bp <= 0:
            continue
        p3k_raw += 1
        site_id = row.get('SiteID', '') or f"anon_{lat:.4f}_{lon:.4f}"
        site_groups[site_id].append({
            'lat': lat, 'lon': lon, 'age_bp': age_bp,
            'name': row.get('SiteName', ''),
        })

# Deduplicate: one entry per site, use oldest date
p3k_count = 0
for sid, rows in site_groups.items():
    lat = np.mean([r['lat'] for r in rows])
    lon = np.mean([r['lon'] for r in rows])
    age_bp = max(r['age_bp'] for r in rows)
    cal_bp = c14_to_cal_bp(age_bp)
    name = rows[0]['name']
    d = dist_from_circle(lat, lon)
    if d <= THRESHOLD_KM:
        arc = bearing_from_pole(lat, lon)
        site_type = classify_name(name)
        all_sites.append({
            'lat': round(lat, 5),
            'lon': round(lon, 5),
            'date_bp': round(cal_bp),
            'name': name or f"p3k14c_{sid}",
            'type': site_type,
            'source': 'p3k14c',
            'dist_km': round(d, 1),
            'arc_deg': round(arc, 2),
        })
        p3k_count += 1

print(f"  Raw rows: {p3k_raw}, Near circle (≤{THRESHOLD_KM}km): {p3k_count}")

# --- 2. XRONOS ---
print("\n[2/3] Loading XRONOS...")
xronos_file = os.path.join(BASE_DIR, "data", "xronos", "xronos_sites.csv")
xronos_count = 0

if os.path.exists(xronos_file):
    with open(xronos_file, encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                lat = float(row['lat'])
                lon = float(row['lon'])
            except (ValueError, TypeError):
                continue
            if lat == 0 and lon == 0:
                continue
            if not (-90 <= lat <= 90 and -180 <= lon <= 180):
                continue

            d = dist_from_circle(lat, lon)
            if d > THRESHOLD_KM:
                continue

            # XRONOS doesn't have direct dating in the sites CSV;
            # use site_type/feature_type for classification
            name = row.get('site_name', '')
            site_type_raw = row.get('site_type', '') or ''
            feature_type_raw = row.get('feature_type', '') or ''
            combined = f"{name} {site_type_raw} {feature_type_raw}"
            site_type = classify_name(combined)

            # Skip XRONOS sites without dates (n_dates field just counts available dates)
            # We can't assign a time period without a date, so skip undated
            # XRONOS sites will supplement spatial coverage but won't have temporal assignment
            # unless we cross-reference with the raw JSON
            # For now, skip XRONOS from temporal binning (p3k14c covers most of the same C14 dates)
            # but include them as "undated" spatial context
            arc = bearing_from_pole(lat, lon)
            xronos_count += 1
            # Note: not adding to all_sites for temporal analysis since no date field
            # They'll be noted in the summary

    print(f"  Near circle (≤{THRESHOLD_KM}km): {xronos_count} (spatial only, no dates in sites CSV)")
else:
    print(f"  XRONOS file not found: {xronos_file}")

# --- 3. Pleiades ---
print("\n[3/3] Loading Pleiades...")
pleiades_file = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")
pleiades_count = 0
pleiades_no_date = 0

with open(pleiades_file, encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            lat = float(row['reprLat'])
            lon = float(row['reprLong'])
        except (ValueError, KeyError, TypeError):
            continue
        if lat == 0 and lon == 0:
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue

        d = dist_from_circle(lat, lon)
        if d > THRESHOLD_KM:
            continue

        # Get minDate (calendar year, negative = BCE)
        try:
            min_date = float(row.get('minDate', ''))
        except (ValueError, TypeError):
            pleiades_no_date += 1
            continue

        date_bp = cal_year_to_bp(min_date)
        if date_bp < 0:
            continue  # future dates

        # Classify by featureTypes
        feature_types = {t.strip() for t in row.get('featureTypes', '').split(',') if t.strip()}
        if feature_types & PLEIADES_MONUMENTAL:
            site_type = "monumental"
        elif feature_types & PLEIADES_SETTLEMENT:
            site_type = "settlement"
        else:
            site_type = classify_name(row.get('title', ''))

        name = row.get('title', '') or f"pleiades_{row.get('id', '')}"
        arc = bearing_from_pole(lat, lon)

        all_sites.append({
            'lat': round(lat, 5),
            'lon': round(lon, 5),
            'date_bp': round(date_bp),
            'name': name,
            'type': site_type,
            'source': 'pleiades',
            'dist_km': round(d, 1),
            'arc_deg': round(arc, 2),
        })
        pleiades_count += 1

print(f"  Near circle with dates: {pleiades_count}, No date: {pleiades_no_date}")

# --- Deduplicate across sources (within 1km) ---
print(f"\nTotal sites before dedup: {len(all_sites)}")

# Sort by source priority (p3k14c first — has C14 dates, then Pleiades)
all_sites.sort(key=lambda s: (0 if s['source'] == 'p3k14c' else 1, s['lat'], s['lon']))

deduped = []
used_coords = set()
for site in all_sites:
    key = (round(site['lat'], 2), round(site['lon'], 2))  # ~1km grid
    if key not in used_coords:
        used_coords.add(key)
        deduped.append(site)

all_sites = deduped
print(f"After dedup (~1km grid): {len(all_sites)}")

# ============================================================
# BIN INTO TIME PERIODS
# ============================================================
print("\n" + "=" * 70)
print("BINNING SITES INTO TIME PERIODS")
print("=" * 70)

np.random.seed(SEED)

# Collect ALL sites (not just near-circle) for enrichment baseline
# We need the full lat/lon arrays from each source
print("\nBuilding full-dataset arrays for enrichment baseline...")
all_p3k_lats, all_p3k_lons = [], []
with open(p3k_file, encoding='utf-8') as f:
    seen_ids = set()
    for row in csv.DictReader(f):
        try:
            lat = float(row['Lat'])
            lon = float(row['Long'])
            age_bp = float(row['Age']) if row['Age'] else None
        except (ValueError, TypeError):
            continue
        if lat == 0 and lon == 0:
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue
        if age_bp is None or age_bp <= 0:
            continue
        sid = row.get('SiteID', '') or f"anon_{lat:.4f}_{lon:.4f}"
        if sid not in seen_ids:
            seen_ids.add(sid)
            all_p3k_lats.append(lat)
            all_p3k_lons.append(lon)

all_pleiades_lats, all_pleiades_lons = [], []
with open(pleiades_file, encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            lat = float(row['reprLat'])
            lon = float(row['reprLong'])
        except (ValueError, KeyError, TypeError):
            continue
        if lat == 0 and lon == 0:
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue
        all_pleiades_lats.append(lat)
        all_pleiades_lons.append(lon)

# Merge for combined baseline
all_baseline_lats = np.array(all_p3k_lats + all_pleiades_lats)
all_baseline_lons = np.array(all_p3k_lons + all_pleiades_lons)
print(f"  Baseline sites: {len(all_baseline_lats)} (p3k14c unique: {len(all_p3k_lats)}, Pleiades: {len(all_pleiades_lats)})")

timeline_data = []

for start_bp, end_bp, label in PERIODS:
    # Filter sites in this period
    period_sites = [s for s in all_sites if end_bp <= s['date_bp'] < start_bp]

    # Compute enrichment against random circles
    enrichment_result = compute_enrichment(
        period_sites, all_baseline_lats, all_baseline_lons,
        THRESHOLD_KM, N_TRIALS
    )

    period_entry = {
        'period_label': label,
        'start_bp': start_bp,
        'end_bp': end_bp,
        'start_cal': bp_to_cal_year(start_bp),
        'end_cal': bp_to_cal_year(end_bp),
        'sites': [{
            'lat': s['lat'],
            'lon': s['lon'],
            'name': s['name'],
            'type': s['type'],
            'arc_deg': s['arc_deg'],
            'dist_km': s['dist_km'],
            'date_bp': s['date_bp'],
            'source': s['source'],
        } for s in period_sites],
        'count': len(period_sites),
        'by_type': {
            'monumental': sum(1 for s in period_sites if s['type'] == 'monumental'),
            'settlement': sum(1 for s in period_sites if s['type'] == 'settlement'),
            'unknown': sum(1 for s in period_sites if s['type'] == 'unknown'),
        },
        'by_source': {
            'p3k14c': sum(1 for s in period_sites if s['source'] == 'p3k14c'),
            'pleiades': sum(1 for s in period_sites if s['source'] == 'pleiades'),
        },
        'enrichment': enrichment_result,
    }
    timeline_data.append(period_entry)

    print(f"\n  {label}")
    print(f"    Sites: {len(period_sites)} (monumental: {period_entry['by_type']['monumental']}, "
          f"settlement: {period_entry['by_type']['settlement']}, unknown: {period_entry['by_type']['unknown']})")
    print(f"    Enrichment: {enrichment_result['enrichment']:.2f}x, Z={enrichment_result['z_score']:.2f}")

# ============================================================
# OUTPUT 1: timeline_data.json
# ============================================================
print("\n" + "=" * 70)
print("SAVING OUTPUTS")
print("=" * 70)

output_main = {
    "meta": {
        "analysis": "Timeline Visualization Data — Site Appearance Along the Great Circle",
        "date": "2026-03-21",
        "methodology": {
            "threshold_km": THRESHOLD_KM,
            "monte_carlo_trials": N_TRIALS,
            "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
            "sources": ["p3k14c (C14 dates)", "Pleiades (minDate field)"],
            "dedup_method": "~1km grid, p3k14c priority",
            "c14_calibration": "simplified IntCal20 piecewise linear approximation",
            "note": "XRONOS sites excluded from temporal bins (no date in sites CSV)",
        },
        "total_sites": len(all_sites),
        "periods": len(PERIODS),
        "optimized_for": "Mapbox GL / Cesium / Leaflet timeline plugin",
    },
    "periods": timeline_data,
}

path1 = os.path.join(OUT_DIR, "timeline_data.json")
with open(path1, "w") as f:
    json.dump(output_main, f, indent=2, default=str)
print(f"  [1/4] {path1}")

# ============================================================
# OUTPUT 2: timeline_geojson.json — GeoJSON FeatureCollection
# ============================================================
features = []
for period in timeline_data:
    for site in period['sites']:
        features.append({
            "type": "Feature",
            "geometry": {
                "type": "Point",
                "coordinates": [site['lon'], site['lat']],
            },
            "properties": {
                "name": site['name'],
                "type": site['type'],
                "date_bp": site['date_bp'],
                "date_cal": bp_to_cal_year(site['date_bp']),
                "period": period['period_label'],
                "period_start_bp": period['start_bp'],
                "period_end_bp": period['end_bp'],
                "arc_deg": site['arc_deg'],
                "dist_km": site['dist_km'],
                "source": site['source'],
            },
        })

geojson = {
    "type": "FeatureCollection",
    "properties": {
        "title": "Great Circle Timeline — Archaeological Sites by Period",
        "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
        "threshold_km": THRESHOLD_KM,
        "total_features": len(features),
    },
    "features": features,
}

path2 = os.path.join(OUT_DIR, "timeline_geojson.json")
with open(path2, "w") as f:
    json.dump(geojson, f, indent=2, default=str)
print(f"  [2/4] {path2}")

# ============================================================
# OUTPUT 3: timeline_summary.json
# ============================================================
summary = {
    "meta": output_main["meta"],
    "periods": [],
}
for period in timeline_data:
    summary["periods"].append({
        "label": period["period_label"],
        "start_bp": period["start_bp"],
        "end_bp": period["end_bp"],
        "start_cal": period["start_cal"],
        "end_cal": period["end_cal"],
        "count": period["count"],
        "by_type": period["by_type"],
        "by_source": period["by_source"],
        "enrichment": period["enrichment"]["enrichment"],
        "z_score": period["enrichment"]["z_score"],
        "expected": period["enrichment"]["expected"],
    })

# Overall stats
total = sum(p["count"] for p in summary["periods"])
summary["overall"] = {
    "total_dated_sites_near_circle": total,
    "xronos_spatial_only": xronos_count,
    "peak_period": max(summary["periods"], key=lambda p: p["count"])["label"],
    "highest_enrichment_period": max(
        [p for p in summary["periods"] if p["count"] > 0],
        key=lambda p: p["enrichment"] if not math.isnan(p["enrichment"]) else 0
    )["label"] if any(p["count"] > 0 for p in summary["periods"]) else "none",
}

path3 = os.path.join(OUT_DIR, "timeline_summary.json")
with open(path3, "w") as f:
    json.dump(summary, f, indent=2, default=str)
print(f"  [3/4] {path3}")

# ============================================================
# OUTPUT 4: timeline_static.png — multi-panel map
# ============================================================
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    import matplotlib.colors as mcolors

    # Generate great circle path for plotting
    def gc_path(pole_lat, pole_lon, n_points=500):
        """Generate lat/lon points along the great circle."""
        plat = math.radians(pole_lat)
        plon = math.radians(pole_lon)
        lats, lons = [], []
        for i in range(n_points):
            az = 2 * math.pi * i / n_points
            # Point at 90° from pole along azimuth
            lat2 = math.asin(math.sin(plat)*math.cos(math.pi/2) +
                             math.cos(plat)*math.sin(math.pi/2)*math.cos(az))
            lon2 = plon + math.atan2(
                math.sin(az)*math.sin(math.pi/2)*math.cos(plat),
                math.cos(math.pi/2) - math.sin(plat)*math.sin(lat2))
            lats.append(math.degrees(lat2))
            lons.append(math.degrees(lon2))
        return lons, lats

    gc_lons, gc_lats = gc_path(POLE_LAT, POLE_LON)

    # Color map for site types
    type_colors = {
        'monumental': '#e63946',
        'settlement': '#457b9d',
        'unknown': '#999999',
    }

    n_periods = len(PERIODS)
    n_cols = 3
    n_rows = (n_periods + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, n_rows * 5),
                              subplot_kw={'projection': None})
    axes = axes.flatten()

    for i, period in enumerate(timeline_data):
        ax = axes[i]

        # World background
        ax.set_xlim(-180, 180)
        ax.set_ylim(-70, 85)
        ax.set_facecolor('#f0f0f0')
        ax.set_aspect('equal')

        # Great circle
        ax.plot(gc_lons, gc_lats, '.', color='gold', markersize=0.5, alpha=0.5, zorder=1)

        # Sites
        if period['sites']:
            for stype, color in type_colors.items():
                sx = [s['lon'] for s in period['sites'] if s['type'] == stype]
                sy = [s['lat'] for s in period['sites'] if s['type'] == stype]
                if sx:
                    ax.scatter(sx, sy, c=color, s=12, alpha=0.7, zorder=3,
                               edgecolors='black', linewidths=0.3, label=stype)

        # Title with count and enrichment
        enrich_str = f"{period['enrichment']['enrichment']:.1f}x" if not math.isnan(period['enrichment'].get('enrichment', float('nan'))) else "n/a"
        z_str = f"Z={period['enrichment']['z_score']:.1f}" if not math.isnan(period['enrichment'].get('z_score', float('nan'))) else ""
        ax.set_title(f"{period['period_label']}\n{period['count']} sites | {enrich_str} enrichment {z_str}",
                     fontsize=9, fontweight='bold')
        ax.tick_params(labelsize=6)

        if i == 0:
            ax.legend(loc='lower left', fontsize=7, markerscale=1.5)

    # Hide unused panels
    for j in range(len(timeline_data), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f"Archaeological Sites Within {THRESHOLD_KM}km of the Great Circle — By Time Period\n"
                 f"Sources: p3k14c + Pleiades | {len(all_sites)} total dated sites",
                 fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    path4 = os.path.join(OUT_DIR, "timeline_static.png")
    plt.savefig(path4, dpi=150, bbox_inches='tight')
    print(f"  [4/4] {path4}")

except ImportError as e:
    print(f"  [4/4] Could not generate plot: {e}")
    print("  Install matplotlib: pip install matplotlib")

# ============================================================
# PRINT SUMMARY TABLE
# ============================================================
print("\n" + "=" * 70)
print("TIMELINE SUMMARY")
print("=" * 70)
print(f"{'Period':<45} {'Sites':>6} {'Mon':>5} {'Set':>5} {'Enrich':>8} {'Z':>6}")
print("-" * 75)
for p in summary["periods"]:
    e = f"{p['enrichment']:.2f}x" if not math.isnan(p['enrichment']) else "n/a"
    z = f"{p['z_score']:.1f}" if not math.isnan(p['z_score']) else "n/a"
    print(f"  {p['label']:<43} {p['count']:>6} {p['by_type']['monumental']:>5} "
          f"{p['by_type']['settlement']:>5} {e:>8} {z:>6}")
print("-" * 75)
print(f"  {'TOTAL':<43} {total:>6}")

print("\nDONE.")
