#!/usr/bin/env python3
"""
Study 3: The Indian Subcontinent Gap Diagnosis
===============================================
The great circle passes through northern India/Pakistan where the Indus Valley
Civilization should produce strong signal, but results show weak signal.
Diagnose whether this is database coverage, site type distribution, genuine
absence, or the major sites being off-circle.
"""

import csv
import json
import math
import os
import random
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── Shared Constants & Geometry ─────────────────────────────────────────────

POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
THRESHOLD_KM = 50

OUT_DIR = Path("/Users/elliotallan/megalith_site_research/outputs/next_wave/south_asia_gap")
OUT_DIR.mkdir(parents=True, exist_ok=True)

DATA_ROOT = Path("/Users/elliotallan/megalith_site_research/data")

random.seed(42)
np.random.seed(42)


def haversine_km(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return 2 * EARTH_R_KM * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def dist_from_gc(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QUARTER_CIRC)


# ── Type Sets ───────────────────────────────────────────────────────────────

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

# ── Major IVC Sites ─────────────────────────────────────────────────────────

MAJOR_IVC = {
    "Mohenjo-daro":  (27.3242, 68.1362),
    "Harappa":       (30.6310, 72.8651),
    "Lothal":        (22.5218, 72.2498),
    "Dholavira":     (23.8876, 70.2129),
    "Kalibangan":    (29.4736, 74.1300),
    "Rakhigarhi":    (29.2808, 76.1156),
}

REFERENCE_SITES = {
    "Giza":          6.5,
    "Nazca":         1.0,
    "Easter Island": 10.3,
}

# ── Bounding Boxes ──────────────────────────────────────────────────────────

BOXES = {
    "South Asia":  (5, 35, 60, 95),
    "Egypt":       (22, 31, 24, 36),
    "Anatolia":    (36, 42, 26, 45),
    "Italy":       (36, 47, 7, 19),
}

def in_box(lat, lon, box):
    lat_min, lat_max, lon_min, lon_max = box
    return lat_min <= lat <= lat_max and lon_min <= lon <= lon_max


# ── Data Loading ────────────────────────────────────────────────────────────

def load_pleiades():
    """Load Pleiades places CSV."""
    path = DATA_ROOT / "pleiades" / "pleiades-places-latest.csv"
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                lat = float(r["reprLat"])
                lon = float(r["reprLong"])
            except (ValueError, TypeError):
                continue
            ft_raw = r.get("featureTypes", "")
            rows.append({"lat": lat, "lon": lon, "featureType": ft_raw,
                         "title": r.get("title", "")})
    return rows


def classify_pleiades_type(ft_raw):
    """Classify a Pleiades featureTypes string."""
    if not ft_raw:
        return "unknown"
    fts = {t.strip().lower().replace('"', '').replace("'", "")
           for t in ft_raw.replace(",", " ").split()}
    if fts & MONUMENTAL_TYPES:
        return "monumental"
    if fts & SETTLEMENT_TYPES:
        return "settlement"
    return "other"


def load_p3k14c():
    """Load p3k14c radiocarbon data."""
    path = DATA_ROOT / "p3k14c" / "p3k14c_data.csv"
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                lat = float(r["Lat"])
                lon = float(r["Long"])
            except (ValueError, TypeError):
                continue
            rows.append({"lat": lat, "lon": lon,
                         "site_name": r.get("SiteName", ""),
                         "site_id": r.get("SiteID", ""),
                         "age": r.get("Age", "")})
    return rows


def load_harappan_kmz():
    """Parse Harappan Sites KMZ for site coordinates."""
    kmz_path = Path("/Users/elliotallan/megalith_site_research/"
                    "trade-corridor-analysis/data/raw/harappan/Harappan_Sites.kmz")
    if not kmz_path.exists():
        return []

    sites = []
    ns = {"kml": "http://earth.google.com/kml/2.2"}
    with zipfile.ZipFile(kmz_path) as z:
        with z.open("doc.kml") as f:
            tree = ET.parse(f)
    root = tree.getroot()
    for pm in root.findall(".//kml:Placemark", ns):
        name_el = pm.find("kml:name", ns)
        name = name_el.text.strip() if name_el is not None and name_el.text else ""
        coords_el = pm.find(".//kml:coordinates", ns)
        if coords_el is None or not coords_el.text:
            continue
        parts = coords_el.text.strip().split(",")
        if len(parts) < 2:
            continue
        try:
            lon, lat = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        # Determine period from style
        style_el = pm.find("kml:styleUrl", ns)
        period = ""
        if style_el is not None and style_el.text:
            period = style_el.text.replace("#", "").replace("indus_", "")
        sites.append({"name": name, "lat": lat, "lon": lon, "period": period})
    return sites


def load_south_asia_merged():
    """Load the merged South Asia dataset."""
    path = DATA_ROOT / "south_asia" / "south_asia_merged.json"
    if not path.exists():
        return []
    with open(path) as f:
        data = json.load(f)
    return data.get("sites", [])


def load_wikipedia_ivc():
    """Load Wikipedia IVC geocoded sites."""
    path = DATA_ROOT / "wikidata" / "wikipedia_ivc_geocoded.json"
    if not path.exists():
        return []
    with open(path) as f:
        data = json.load(f)
    return data.get("sites", [])


# ── Phase 1: How Close Are the Major Sites? ─────────────────────────────────

def phase1_major_site_distances():
    print("\n" + "=" * 70)
    print("PHASE 1: Major IVC Site Distances from Great Circle")
    print("=" * 70)

    results = {}
    for name, (lat, lon) in MAJOR_IVC.items():
        d = dist_from_gc(lat, lon)
        results[name] = {
            "lat": lat, "lon": lon, "gc_dist_km": round(d, 2),
            "within_50km": d <= 50,
            "within_100km": d <= 100,
            "within_200km": d <= 200,
        }
        status = "NEAR" if d <= 50 else ("MODERATE" if d <= 100 else "FAR")
        print(f"  {name:20s}  {d:8.1f} km  [{status}]")

    print(f"\n  Reference sites:")
    for name, d in REFERENCE_SITES.items():
        print(f"  {name:20s}  {d:8.1f} km")

    within_50 = sum(1 for v in results.values() if v["within_50km"])
    within_100 = sum(1 for v in results.values() if v["within_100km"])
    within_200 = sum(1 for v in results.values() if v["within_200km"])
    print(f"\n  Summary: {within_50}/6 within 50km, {within_100}/6 within 100km, "
          f"{within_200}/6 within 200km")

    return results


# ── Phase 2: Database Coverage Audit ─────────────────────────────────────────

def phase2_coverage_audit(pleiades, p3k14c):
    print("\n" + "=" * 70)
    print("PHASE 2: Database Coverage Audit")
    print("=" * 70)

    results = {}
    for db_name, dataset, lat_key, lon_key in [
        ("Pleiades", pleiades, "lat", "lon"),
        ("p3k14c", p3k14c, "lat", "lon"),
    ]:
        counts = {}
        for region, box in BOXES.items():
            count = sum(1 for s in dataset if in_box(s[lat_key], s[lon_key], box))
            counts[region] = count
        results[db_name] = counts

        print(f"\n  {db_name}:")
        sa = counts["South Asia"]
        for region, count in counts.items():
            ratio = (sa / count * 100) if count > 0 else float("inf")
            print(f"    {region:15s}  {count:7d} sites"
                  + (f"  (SA = {ratio:.0f}% of this)" if region != "South Asia" else ""))

    return results


# ── Phase 3: Type Distribution ───────────────────────────────────────────────

def phase3_type_distribution(pleiades):
    print("\n" + "=" * 70)
    print("PHASE 3: Pleiades Type Distribution by Region")
    print("=" * 70)

    results = {}
    for region, box in BOXES.items():
        region_sites = [s for s in pleiades if in_box(s["lat"], s["lon"], box)]
        total = len(region_sites)
        if total == 0:
            results[region] = {"total": 0, "monumental": 0, "settlement": 0,
                               "other": 0, "mon_frac": 0, "set_frac": 0}
            continue

        type_counts = Counter()
        for s in region_sites:
            cls = classify_pleiades_type(s["featureType"])
            type_counts[cls] += 1

        mon = type_counts.get("monumental", 0)
        sett = type_counts.get("settlement", 0)
        oth = type_counts.get("other", 0) + type_counts.get("unknown", 0)

        results[region] = {
            "total": total,
            "monumental": mon,
            "settlement": sett,
            "other": oth,
            "mon_frac": round(mon / total, 3),
            "set_frac": round(sett / total, 3),
        }
        print(f"  {region:15s}  total={total:5d}  "
              f"mon={mon:4d} ({mon/total*100:5.1f}%)  "
              f"sett={sett:4d} ({sett/total*100:5.1f}%)  "
              f"other={oth:4d} ({oth/total*100:5.1f}%)")

    return results


# ── Phase 4: South Asia Divergence (Monte Carlo) ────────────────────────────

def phase4_divergence(pleiades, p3k14c):
    print("\n" + "=" * 70)
    print("PHASE 4: South Asia Divergence Analysis")
    print("=" * 70)

    sa_box = BOXES["South Asia"]

    # Gather Pleiades SA sites with classification
    sa_pleiades = []
    for s in pleiades:
        if in_box(s["lat"], s["lon"], sa_box):
            cls = classify_pleiades_type(s["featureType"])
            if cls in ("monumental", "settlement"):
                sa_pleiades.append({"lat": s["lat"], "lon": s["lon"],
                                    "type": cls, "source": "pleiades"})

    # Gather p3k14c SA sites (classify as settlement by default)
    sa_p3k = []
    for s in p3k14c:
        if in_box(s["lat"], s["lon"], sa_box):
            sa_p3k.append({"lat": s["lat"], "lon": s["lon"],
                           "type": "settlement", "source": "p3k14c"})

    # Combine
    sa_all = sa_pleiades + sa_p3k

    monuments = [s for s in sa_all if s["type"] == "monumental"]
    settlements = [s for s in sa_all if s["type"] == "settlement"]

    mon_near = sum(1 for s in monuments if dist_from_gc(s["lat"], s["lon"]) <= THRESHOLD_KM)
    sett_near = sum(1 for s in settlements if dist_from_gc(s["lat"], s["lon"]) <= THRESHOLD_KM)

    n_mon = len(monuments)
    n_sett = len(settlements)

    print(f"  Monuments:  {n_mon:5d} total, {mon_near:3d} within {THRESHOLD_KM}km "
          f"({mon_near/n_mon*100:.1f}%)" if n_mon else "  Monuments: 0")
    print(f"  Settlements:{n_sett:5d} total, {sett_near:3d} within {THRESHOLD_KM}km "
          f"({sett_near/n_sett*100:.1f}%)" if n_sett else "  Settlements: 0")

    # MC simulation: shuffle type labels, count near-GC for each type
    N_TRIALS = 500
    mon_null = []
    sett_null = []
    labels = [s["type"] for s in sa_all]
    coords = [(s["lat"], s["lon"]) for s in sa_all]
    near_flags = [dist_from_gc(lat, lon) <= THRESHOLD_KM for lat, lon in coords]

    for _ in range(N_TRIALS):
        shuffled = labels.copy()
        random.shuffle(shuffled)
        mc_mon = sum(1 for i, lbl in enumerate(shuffled)
                     if lbl == "monumental" and near_flags[i])
        mc_sett = sum(1 for i, lbl in enumerate(shuffled)
                      if lbl == "settlement" and near_flags[i])
        mon_null.append(mc_mon)
        sett_null.append(mc_sett)

    mon_mean = np.mean(mon_null)
    mon_std = np.std(mon_null) if np.std(mon_null) > 0 else 1
    sett_mean = np.mean(sett_null)
    sett_std = np.std(sett_null) if np.std(sett_null) > 0 else 1

    mon_z = (mon_near - mon_mean) / mon_std
    sett_z = (sett_near - sett_mean) / sett_std
    divergence = mon_z - sett_z

    print(f"\n  MC ({N_TRIALS} trials):")
    print(f"    Monument Z-score:   {mon_z:+.2f}  (obs={mon_near}, "
          f"null mean={mon_mean:.1f}, std={mon_std:.2f})")
    print(f"    Settlement Z-score: {sett_z:+.2f}  (obs={sett_near}, "
          f"null mean={sett_mean:.1f}, std={sett_std:.2f})")
    print(f"    DIVERGENCE (mon_Z - sett_Z): {divergence:+.2f}")

    return {
        "n_monuments": n_mon,
        "n_settlements": n_sett,
        "monuments_near_gc": mon_near,
        "settlements_near_gc": sett_near,
        "mon_z": round(mon_z, 3),
        "sett_z": round(sett_z, 3),
        "divergence": round(divergence, 3),
        "n_trials": N_TRIALS,
    }


# ── Phase 4b: Divergence using South Asia Merged dataset ────────────────────

def phase4b_divergence_merged():
    print("\n" + "-" * 50)
    print("  Phase 4b: Merged South Asia Dataset Divergence")
    print("-" * 50)

    sa_merged = load_south_asia_merged()
    if not sa_merged:
        print("  [SKIP] south_asia_merged.json not found")
        return None

    monuments = [s for s in sa_merged if s.get("classification") == "monumental"]
    settlements = [s for s in sa_merged if s.get("classification") in
                   ("settlement", "later_modern", "settlement_domestic")]
    # Also include any with gc_distance already computed
    all_classified = monuments + settlements

    mon_near = sum(1 for s in monuments
                   if dist_from_gc(s["lat"], s["lon"]) <= THRESHOLD_KM)
    sett_near = sum(1 for s in settlements
                    if dist_from_gc(s["lat"], s["lon"]) <= THRESHOLD_KM)

    n_mon = len(monuments)
    n_sett = len(settlements)

    print(f"  Monuments:   {n_mon:5d} total, {mon_near:3d} within {THRESHOLD_KM}km")
    print(f"  Settlements: {n_sett:5d} total, {sett_near:3d} within {THRESHOLD_KM}km")

    # MC
    N_TRIALS = 500
    labels = (["monumental"] * n_mon) + (["settlement"] * n_sett)
    coords = [(s["lat"], s["lon"]) for s in monuments + settlements]
    near_flags = [dist_from_gc(lat, lon) <= THRESHOLD_KM for lat, lon in coords]

    mon_null, sett_null = [], []
    for _ in range(N_TRIALS):
        shuffled = labels.copy()
        random.shuffle(shuffled)
        mc_mon = sum(1 for i, lbl in enumerate(shuffled)
                     if lbl == "monumental" and near_flags[i])
        mc_sett = sum(1 for i, lbl in enumerate(shuffled)
                      if lbl == "settlement" and near_flags[i])
        mon_null.append(mc_mon)
        sett_null.append(mc_sett)

    mon_mean = np.mean(mon_null)
    mon_std = max(np.std(mon_null), 1)
    sett_mean = np.mean(sett_null)
    sett_std = max(np.std(sett_null), 1)

    mon_z = (mon_near - mon_mean) / mon_std
    sett_z = (sett_near - sett_mean) / sett_std
    divergence = mon_z - sett_z

    print(f"  MC: mon_Z={mon_z:+.2f}, sett_Z={sett_z:+.2f}, divergence={divergence:+.2f}")

    return {
        "n_monuments": n_mon,
        "n_settlements": n_sett,
        "monuments_near_gc": mon_near,
        "settlements_near_gc": sett_near,
        "mon_z": round(mon_z, 3),
        "sett_z": round(sett_z, 3),
        "divergence": round(divergence, 3),
        "n_trials": N_TRIALS,
        "source": "south_asia_merged.json",
    }


# ── Phase 5: GC Segment Analysis ────────────────────────────────────────────

def gc_point(pole_lat, pole_lon, bearing_deg, dist_km):
    """Given a pole, bearing, and distance, compute the point on the sphere."""
    lat1 = math.radians(pole_lat)
    lon1 = math.radians(pole_lon)
    d = dist_km / EARTH_R_KM
    brng = math.radians(bearing_deg)

    lat2 = math.asin(math.sin(lat1) * math.cos(d)
                     + math.cos(lat1) * math.sin(d) * math.cos(brng))
    lon2 = lon1 + math.atan2(math.sin(brng) * math.sin(d) * math.cos(lat1),
                              math.cos(d) - math.sin(lat1) * math.sin(lat2))
    return math.degrees(lat2), math.degrees(lon2)


def normalize_lon(lon):
    """Normalize longitude to [-180, 180]."""
    return (lon + 180) % 360 - 180


def sample_gc_points(n=360):
    """Sample n points along the great circle."""
    points = []
    for i in range(n):
        bearing = i * (360 / n)
        lat, lon = gc_point(POLE_LAT, POLE_LON, bearing, QUARTER_CIRC)
        lon = normalize_lon(lon)
        points.append((lat, lon))
    return points


def phase5_gc_segment():
    print("\n" + "=" * 70)
    print("PHASE 5: Great Circle Segment Through South Asia")
    print("=" * 70)

    gc_points = sample_gc_points(720)

    # Filter to South Asia longitude range 60-95E
    sa_segment = [(lat, lon) for lat, lon in gc_points
                  if 55 <= lon <= 100]

    if not sa_segment:
        print("  No GC points found in South Asia longitude range!")
        return {"segment_points": []}

    sa_segment.sort(key=lambda p: p[1])  # sort by longitude

    print(f"  GC segment through South Asia ({len(sa_segment)} sample points):")
    # Report key waypoints
    waypoints = []
    for i in range(0, len(sa_segment), max(1, len(sa_segment) // 10)):
        lat, lon = sa_segment[i]
        waypoints.append({"lat": round(lat, 2), "lon": round(lon, 2)})
        terrain = describe_terrain(lat, lon)
        print(f"    {lat:6.2f}N  {lon:6.2f}E  — {terrain}")

    # Lat range tells us what part of India it crosses
    lats = [p[0] for p in sa_segment]
    lons = [p[1] for p in sa_segment]
    print(f"\n  Latitude range in SA: {min(lats):.1f} to {max(lats):.1f}")
    print(f"  Longitude range:      {min(lons):.1f} to {max(lons):.1f}")

    return {
        "segment_points": [{"lat": round(lat, 2), "lon": round(lon, 2)}
                           for lat, lon in sa_segment],
        "lat_range": [round(min(lats), 2), round(max(lats), 2)],
        "lon_range": [round(min(lons), 2), round(max(lons), 2)],
    }


def describe_terrain(lat, lon):
    """Rough geographic description based on lat/lon."""
    if 59 <= lon <= 63:
        return "Eastern Iran / Balochistan border"
    if 63 < lon <= 68:
        if lat > 30:
            return "Northern Balochistan / FATA (Pakistan)"
        elif lat > 26:
            return "Sindh (Pakistan) — Indus Valley heartland"
        else:
            return "Southern Sindh / Thar Desert fringe"
    if 68 < lon <= 72:
        if lat > 30:
            return "Punjab (Pakistan) — Harappa region"
        elif lat > 26:
            return "Eastern Sindh / Thar Desert"
        else:
            return "Kutch / Gujarat coast"
    if 72 < lon <= 77:
        if lat > 30:
            return "Punjab / Haryana (India) — Rakhigarhi region"
        elif lat > 26:
            return "Rajasthan / Gujarat — Dholavira region"
        else:
            return "Gujarat coast / Kathiawar"
    if 77 < lon <= 82:
        if lat > 28:
            return "Gangetic Plain / UP"
        elif lat > 22:
            return "Central India / Malwa Plateau"
        else:
            return "Western Ghats / Deccan"
    if 82 < lon <= 88:
        if lat > 25:
            return "Eastern Gangetic Plain / Bihar"
        else:
            return "Odisha / Eastern Ghats"
    if 88 < lon <= 95:
        if lat > 22:
            return "Bengal / Bangladesh"
        else:
            return "Bay of Bengal coast"
    return "South Asia (unspecified)"


# ── Harappan KMZ Analysis ───────────────────────────────────────────────────

def analyze_harappan_kmz():
    print("\n" + "=" * 70)
    print("BONUS: Harappan KMZ Dataset Analysis")
    print("=" * 70)

    sites = load_harappan_kmz()
    if not sites:
        print("  [SKIP] KMZ not found or empty")
        return None

    print(f"  Total Harappan sites from KMZ: {len(sites)}")

    # Distance distribution
    dists = []
    near_sites = []
    for s in sites:
        d = dist_from_gc(s["lat"], s["lon"])
        dists.append(d)
        s["gc_dist_km"] = round(d, 2)
        if d <= THRESHOLD_KM:
            near_sites.append(s)

    within_50 = sum(1 for d in dists if d <= 50)
    within_100 = sum(1 for d in dists if d <= 100)
    within_200 = sum(1 for d in dists if d <= 200)

    print(f"  Within 50km:  {within_50:4d} ({within_50/len(sites)*100:.1f}%)")
    print(f"  Within 100km: {within_100:4d} ({within_100/len(sites)*100:.1f}%)")
    print(f"  Within 200km: {within_200:4d} ({within_200/len(sites)*100:.1f}%)")
    print(f"  Median distance: {np.median(dists):.0f} km")
    print(f"  Mean distance:   {np.mean(dists):.0f} km")

    # Period breakdown
    periods = Counter(s["period"] for s in sites)
    print(f"\n  Period breakdown:")
    for period, count in periods.most_common():
        near_in_period = sum(1 for s in near_sites if s["period"] == period)
        print(f"    {period:30s}  {count:4d} total, {near_in_period:3d} near GC")

    # What fraction SHOULD be near GC if randomly distributed?
    # SA box area approx: (35-5)*(95-60)*cos(20) ~= 30*35*0.94 = 987 deg^2
    # 50km strip width ~ 0.45 deg lat, spanning ~30 deg lon = 13.5 deg^2
    # Expected fraction ~ 13.5/987 ~ 1.4%
    expected_frac = 0.014
    actual_frac = within_50 / len(sites) if len(sites) > 0 else 0
    enrichment = actual_frac / expected_frac if expected_frac > 0 else 0

    print(f"\n  Random expectation (50km strip): ~{expected_frac*100:.1f}%")
    print(f"  Actual near-GC fraction:         {actual_frac*100:.1f}%")
    print(f"  Enrichment factor:               {enrichment:.2f}x")

    return {
        "total_sites": len(sites),
        "within_50km": within_50,
        "within_100km": within_100,
        "within_200km": within_200,
        "median_dist_km": round(float(np.median(dists)), 1),
        "mean_dist_km": round(float(np.mean(dists)), 1),
        "enrichment_50km": round(enrichment, 2),
        "period_counts": dict(periods),
        "nearest_sites": sorted([{"name": s["name"], "gc_dist_km": s["gc_dist_km"],
                                   "period": s["period"]}
                                  for s in near_sites],
                                 key=lambda x: x["gc_dist_km"])[:20],
    }


# ── Wikipedia IVC Analysis ──────────────────────────────────────────────────

def analyze_wikipedia_ivc():
    print("\n" + "-" * 50)
    print("  Wikipedia IVC Geocoded Sites")
    print("-" * 50)

    sites = load_wikipedia_ivc()
    if not sites:
        print("  [SKIP] Not found")
        return None

    print(f"  Total sites: {len(sites)}")
    dists = []
    for s in sites:
        d = dist_from_gc(s["lat"], s["lon"])
        dists.append(d)
    within_50 = sum(1 for d in dists if d <= 50)
    within_100 = sum(1 for d in dists if d <= 100)
    print(f"  Within 50km:  {within_50}")
    print(f"  Within 100km: {within_100}")
    print(f"  Median distance: {np.median(dists):.0f} km")

    return {
        "total": len(sites),
        "within_50km": within_50,
        "within_100km": within_100,
        "median_dist_km": round(float(np.median(dists)), 1),
    }


# ── Plotting ────────────────────────────────────────────────────────────────

def plot_ivc_distances(phase1_results, harappan_results):
    """Plot major IVC site distances vs reference sites."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: Major IVC sites bar chart
    ax = axes[0]
    names = list(phase1_results.keys())
    dists = [phase1_results[n]["gc_dist_km"] for n in names]
    colors = ["#2ecc71" if d <= 50 else "#f39c12" if d <= 100 else "#e74c3c"
              for d in dists]
    bars = ax.barh(names, dists, color=colors, edgecolor="black", linewidth=0.5)
    ax.axvline(50, color="green", linestyle="--", alpha=0.7, label="50km threshold")
    ax.axvline(100, color="orange", linestyle="--", alpha=0.7, label="100km")
    ax.axvline(200, color="red", linestyle="--", alpha=0.7, label="200km")

    # Add reference lines
    for ref_name, ref_dist in REFERENCE_SITES.items():
        ax.axvline(ref_dist, color="blue", linestyle=":", alpha=0.4)
        ax.text(ref_dist + 2, len(names) - 0.3, ref_name, fontsize=7,
                color="blue", alpha=0.6)

    ax.set_xlabel("Distance from Great Circle (km)")
    ax.set_title("Major IVC Sites: Distance from Great Circle")
    ax.legend(fontsize=8)
    ax.set_xlim(0, max(dists) * 1.2)

    # Right: Harappan KMZ distance histogram
    ax2 = axes[1]
    if harappan_results:
        # Recompute all distances for histogram
        kmz_sites = load_harappan_kmz()
        all_dists = [dist_from_gc(s["lat"], s["lon"]) for s in kmz_sites]
        ax2.hist(all_dists, bins=50, color="#3498db", edgecolor="black",
                 linewidth=0.3, alpha=0.8)
        ax2.axvline(50, color="green", linestyle="--", linewidth=2,
                    label=f"50km (n={harappan_results['within_50km']})")
        ax2.axvline(100, color="orange", linestyle="--", linewidth=2,
                    label=f"100km (n={harappan_results['within_100km']})")
        ax2.set_xlabel("Distance from Great Circle (km)")
        ax2.set_ylabel("Count")
        ax2.set_title(f"All {len(kmz_sites)} Harappan Sites: GC Distance Distribution")
        ax2.legend(fontsize=8)
    else:
        ax2.text(0.5, 0.5, "Harappan KMZ\nnot available", ha="center", va="center",
                 fontsize=14, transform=ax2.transAxes)
        ax2.set_title("Harappan Sites (unavailable)")

    plt.tight_layout()
    path = OUT_DIR / "ivc_distances.png"
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved: {path}")


def plot_coverage_comparison(phase2_results, phase3_results):
    """Plot database coverage comparison across regions."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: Site counts by region/database
    ax = axes[0]
    regions = list(BOXES.keys())
    x = np.arange(len(regions))
    w = 0.35

    pleiades_counts = [phase2_results["Pleiades"][r] for r in regions]
    p3k_counts = [phase2_results["p3k14c"][r] for r in regions]

    ax.bar(x - w / 2, pleiades_counts, w, label="Pleiades", color="#2ecc71",
           edgecolor="black", linewidth=0.5)
    ax.bar(x + w / 2, p3k_counts, w, label="p3k14c", color="#3498db",
           edgecolor="black", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(regions)
    ax.set_ylabel("Number of Sites")
    ax.set_title("Database Coverage by Region")
    ax.legend()
    ax.set_yscale("log")

    # Right: Type distribution
    ax2 = axes[1]
    mon_fracs = [phase3_results[r]["mon_frac"] * 100 for r in regions]
    set_fracs = [phase3_results[r]["set_frac"] * 100 for r in regions]
    oth_fracs = [100 - mon_fracs[i] - set_fracs[i] for i in range(len(regions))]

    ax2.bar(x, mon_fracs, w * 2, label="Monumental", color="#e74c3c",
            edgecolor="black", linewidth=0.5)
    ax2.bar(x, set_fracs, w * 2, bottom=mon_fracs, label="Settlement",
            color="#3498db", edgecolor="black", linewidth=0.5)
    ax2.bar(x, oth_fracs, w * 2,
            bottom=[mon_fracs[i] + set_fracs[i] for i in range(len(regions))],
            label="Other", color="#95a5a6", edgecolor="black", linewidth=0.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(regions)
    ax2.set_ylabel("Percentage")
    ax2.set_title("Pleiades Type Distribution by Region")
    ax2.legend(fontsize=8)

    plt.tight_layout()
    path = OUT_DIR / "coverage_comparison.png"
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved: {path}")


def plot_gc_segment(phase5_results, phase1_results, harappan_results):
    """Plot GC segment through South Asia with IVC sites."""
    fig, ax = plt.subplots(figsize=(12, 8))

    # GC segment
    seg = phase5_results.get("segment_points", [])
    if seg:
        gc_lons = [p["lon"] for p in seg]
        gc_lats = [p["lat"] for p in seg]
        ax.plot(gc_lons, gc_lats, "r-", linewidth=2, label="Great Circle", zorder=3)

    # Major IVC sites
    for name, data in phase1_results.items():
        color = "green" if data["gc_dist_km"] <= 50 else (
            "orange" if data["gc_dist_km"] <= 100 else "red")
        ax.plot(data["lon"], data["lat"], "D", color=color, markersize=10,
                markeredgecolor="black", zorder=5)
        ax.annotate(f"{name}\n({data['gc_dist_km']:.0f}km)",
                    (data["lon"], data["lat"]),
                    textcoords="offset points", xytext=(8, 5),
                    fontsize=7, fontweight="bold")

    # Harappan KMZ sites
    if harappan_results:
        kmz_sites = load_harappan_kmz()
        near = [(s["lon"], s["lat"]) for s in kmz_sites
                if dist_from_gc(s["lat"], s["lon"]) <= 50]
        far = [(s["lon"], s["lat"]) for s in kmz_sites
               if dist_from_gc(s["lat"], s["lon"]) > 50]
        if far:
            ax.scatter([p[0] for p in far], [p[1] for p in far],
                       s=3, c="#aaaaaa", alpha=0.3, label=f"Harappan >50km (n={len(far)})",
                       zorder=1)
        if near:
            ax.scatter([p[0] for p in near], [p[1] for p in near],
                       s=8, c="#2ecc71", alpha=0.6,
                       label=f"Harappan <=50km (n={len(near)})", zorder=2)

    ax.set_xlim(58, 92)
    ax.set_ylim(18, 36)
    ax.set_xlabel("Longitude (E)")
    ax.set_ylabel("Latitude (N)")
    ax.set_title("Great Circle Through South Asia with IVC Sites")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal")

    plt.tight_layout()
    path = OUT_DIR / "gc_south_asia_map.png"
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved: {path}")


# ── Diagnosis ────────────────────────────────────────────────────────────────

def diagnose(phase1, phase2, phase3, phase4, phase4b, harappan, phase5):
    """Synthesize findings into a diagnosis."""
    print("\n" + "=" * 70)
    print("DIAGNOSIS")
    print("=" * 70)

    factors = []

    # Factor 1: Major site distances
    ivc_dists = {n: v["gc_dist_km"] for n, v in phase1.items()}
    near_major = sum(1 for d in ivc_dists.values() if d <= 50)
    moderate_major = sum(1 for d in ivc_dists.values() if d <= 100)
    if near_major >= 2:
        factors.append(("MAJOR_SITES_NEAR", f"{near_major}/6 major IVC sites within 50km"))
    elif moderate_major >= 3:
        factors.append(("MAJOR_SITES_MODERATE",
                        f"{moderate_major}/6 within 100km, only {near_major}/6 within 50km"))
    else:
        factors.append(("MAJOR_SITES_FAR",
                        f"Only {near_major}/6 within 50km, {moderate_major}/6 within 100km"))

    # Factor 2: Database coverage
    sa_pleiades = phase2["Pleiades"]["South Asia"]
    eg_pleiades = phase2["Pleiades"]["Egypt"]
    coverage_ratio = sa_pleiades / eg_pleiades if eg_pleiades > 0 else 0
    if coverage_ratio < 0.1:
        factors.append(("SEVERE_COVERAGE_GAP",
                        f"SA Pleiades = {sa_pleiades}, Egypt = {eg_pleiades} "
                        f"(ratio {coverage_ratio:.2f})"))
    elif coverage_ratio < 0.5:
        factors.append(("MODERATE_COVERAGE_GAP",
                        f"SA Pleiades = {sa_pleiades}, Egypt = {eg_pleiades} "
                        f"(ratio {coverage_ratio:.2f})"))
    else:
        factors.append(("COVERAGE_ADEQUATE",
                        f"SA/Egypt ratio = {coverage_ratio:.2f}"))

    # Factor 3: Type distribution
    sa_mon_frac = phase3["South Asia"]["mon_frac"]
    eg_mon_frac = phase3["Egypt"]["mon_frac"]
    if sa_mon_frac < eg_mon_frac * 0.3:
        factors.append(("LOW_MONUMENT_FRACTION",
                        f"SA monumental = {sa_mon_frac*100:.1f}%, "
                        f"Egypt = {eg_mon_frac*100:.1f}%"))
    else:
        factors.append(("MONUMENT_FRACTION_OK",
                        f"SA monumental = {sa_mon_frac*100:.1f}%, "
                        f"Egypt = {eg_mon_frac*100:.1f}%"))

    # Factor 4: Divergence
    div = phase4["divergence"]
    if div > 1.5:
        factors.append(("STRONG_DIVERGENCE", f"Monument-settlement divergence = {div:+.2f}"))
    elif div > 0:
        factors.append(("WEAK_DIVERGENCE", f"Monument-settlement divergence = {div:+.2f}"))
    elif div > -1.5:
        factors.append(("NO_DIVERGENCE", f"Monument-settlement divergence = {div:+.2f}"))
    else:
        factors.append(("NEGATIVE_DIVERGENCE", f"Monument-settlement divergence = {div:+.2f}"))

    # Factor 5: Harappan enrichment
    if harappan:
        enrich = harappan["enrichment_50km"]
        if enrich > 1.5:
            factors.append(("HARAPPAN_ENRICHED",
                            f"Harappan sites {enrich:.1f}x enriched near GC"))
        elif enrich > 0.8:
            factors.append(("HARAPPAN_BASELINE",
                            f"Harappan sites {enrich:.1f}x — near random expectation"))
        else:
            factors.append(("HARAPPAN_DEPLETED",
                            f"Harappan sites {enrich:.1f}x — depleted near GC"))

    # Overall diagnosis
    factor_names = {f[0] for f in factors}

    if "SEVERE_COVERAGE_GAP" in factor_names:
        verdict = "EXPLAINED"
        explanation = ("The weak signal is primarily explained by severe database coverage gaps. "
                       "Pleiades has minimal South Asian coverage compared to Mediterranean regions.")
    elif "MAJOR_SITES_FAR" in factor_names:
        verdict = "EXPLAINED"
        explanation = ("Major IVC sites are mostly far from the great circle. "
                       "The GC passes through the region but misses the main centers.")
    elif ("MAJOR_SITES_NEAR" in factor_names and
          "NEGATIVE_DIVERGENCE" in factor_names):
        verdict = "PROBLEMATIC"
        explanation = ("Major sites are near the GC and coverage is adequate, but "
                       "divergence is negative — settlements cluster more than monuments.")
    elif ("MAJOR_SITES_NEAR" in factor_names and
          ("NO_DIVERGENCE" in factor_names or "WEAK_DIVERGENCE" in factor_names)):
        verdict = "INTERESTING"
        explanation = ("Major IVC sites are near the GC but divergence is weak/zero. "
                       "The corridor is real but IVC was settlement-dominated, not monument-dominated.")
    elif "STRONG_DIVERGENCE" in factor_names:
        verdict = "UNEXPECTED_POSITIVE"
        explanation = ("Despite coverage gaps, there is a strong monument-settlement divergence "
                       "in South Asia — warranting further investigation with better data.")
    else:
        verdict = "MIXED"
        explanation = "Multiple factors contribute; no single dominant explanation."

    # Refine with harappan data if available
    if harappan and "HARAPPAN_ENRICHED" in factor_names:
        explanation += (" The Harappan KMZ dataset shows enrichment near the GC, "
                        "suggesting the corridor does intersect IVC territory.")
    elif harappan and "HARAPPAN_DEPLETED" in factor_names:
        explanation += (" The Harappan KMZ dataset shows depletion near the GC, "
                        "confirming the main IVC settlement axis lies away from this corridor.")

    print(f"\n  VERDICT: {verdict}")
    print(f"  {explanation}")
    print(f"\n  Contributing factors:")
    for fname, fdesc in factors:
        print(f"    [{fname}] {fdesc}")

    return {
        "verdict": verdict,
        "explanation": explanation,
        "factors": [{"factor": f, "detail": d} for f, d in factors],
    }


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("Study 3: South Asia Gap Diagnosis")
    print("=" * 70)

    # Load data
    print("\nLoading data...")
    pleiades = load_pleiades()
    print(f"  Pleiades: {len(pleiades)} sites")
    p3k14c = load_p3k14c()
    print(f"  p3k14c:   {len(p3k14c)} samples")

    # Phase 1
    phase1 = phase1_major_site_distances()

    # Phase 2
    phase2 = phase2_coverage_audit(pleiades, p3k14c)

    # Phase 3
    phase3 = phase3_type_distribution(pleiades)

    # Phase 4
    phase4 = phase4_divergence(pleiades, p3k14c)
    phase4b = phase4b_divergence_merged()

    # Harappan KMZ
    harappan = analyze_harappan_kmz()

    # Wikipedia IVC
    wiki_ivc = analyze_wikipedia_ivc()

    # Phase 5
    phase5 = phase5_gc_segment()

    # Diagnosis
    diagnosis = diagnose(phase1, phase2, phase3, phase4, phase4b, harappan, phase5)

    # Plots
    print("\n" + "=" * 70)
    print("GENERATING PLOTS")
    print("=" * 70)
    plot_ivc_distances(phase1, harappan)
    plot_coverage_comparison(phase2, phase3)
    plot_gc_segment(phase5, phase1, harappan)

    # Save results
    results = {
        "study": "S3: South Asia Gap Diagnosis",
        "phase1_major_site_distances": phase1,
        "phase2_coverage_audit": phase2,
        "phase3_type_distribution": phase3,
        "phase4_divergence_pleiades_p3k": phase4,
        "phase4b_divergence_merged": phase4b,
        "harappan_kmz_analysis": harappan,
        "wikipedia_ivc": wiki_ivc,
        "phase5_gc_segment": {
            "lat_range": phase5.get("lat_range"),
            "lon_range": phase5.get("lon_range"),
            "n_sample_points": len(phase5.get("segment_points", [])),
        },
        "diagnosis": diagnosis,
    }

    results_path = OUT_DIR / "results.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved: {results_path}")

    # Generate RESULTS.md
    write_results_md(results, phase1, phase2, phase3, phase4, phase4b,
                     harappan, wiki_ivc, phase5, diagnosis)

    print("\nStudy 3 complete.")


def write_results_md(results, phase1, phase2, phase3, phase4, phase4b,
                     harappan, wiki_ivc, phase5, diagnosis):
    """Write summary markdown."""
    md = []
    md.append("# Study 3: South Asia Gap Diagnosis\n")
    md.append("## Question\n")
    md.append("Why does the great circle show weak signal through northern India/Pakistan,")
    md.append("where the Indus Valley Civilization should produce strong alignment?\n")

    md.append("## Phase 1: Major IVC Site Distances\n")
    md.append("| Site | GC Distance (km) | Within 50km | Within 100km |")
    md.append("|------|-------------------|-------------|--------------|")
    for name, data in phase1.items():
        md.append(f"| {name} | {data['gc_dist_km']:.1f} | "
                  f"{'Yes' if data['within_50km'] else 'No'} | "
                  f"{'Yes' if data['within_100km'] else 'No'} |")
    md.append("")
    md.append("Reference: Giza ~6.5km, Nazca ~1.0km, Easter Island ~10.3km\n")

    md.append("## Phase 2: Database Coverage\n")
    md.append("| Region | Pleiades | p3k14c |")
    md.append("|--------|----------|--------|")
    for region in BOXES:
        md.append(f"| {region} | {phase2['Pleiades'][region]:,} | "
                  f"{phase2['p3k14c'][region]:,} |")
    md.append("")

    md.append("## Phase 3: Type Distribution (Pleiades)\n")
    md.append("| Region | Total | Monumental % | Settlement % |")
    md.append("|--------|-------|-------------|-------------|")
    for region in BOXES:
        d = phase3[region]
        md.append(f"| {region} | {d['total']:,} | "
                  f"{d['mon_frac']*100:.1f}% | {d['set_frac']*100:.1f}% |")
    md.append("")

    md.append("## Phase 4: Divergence Analysis\n")
    md.append("### Pleiades + p3k14c\n")
    md.append(f"- Monuments: {phase4['n_monuments']} total, "
              f"{phase4['monuments_near_gc']} within 50km")
    md.append(f"- Settlements: {phase4['n_settlements']} total, "
              f"{phase4['settlements_near_gc']} within 50km")
    md.append(f"- Monument Z-score: {phase4['mon_z']:+.3f}")
    md.append(f"- Settlement Z-score: {phase4['sett_z']:+.3f}")
    md.append(f"- **Divergence: {phase4['divergence']:+.3f}**\n")

    if phase4b:
        md.append("### Merged South Asia Dataset\n")
        md.append(f"- Monuments: {phase4b['n_monuments']} total, "
                  f"{phase4b['monuments_near_gc']} within 50km")
        md.append(f"- Settlements: {phase4b['n_settlements']} total, "
                  f"{phase4b['settlements_near_gc']} within 50km")
        md.append(f"- **Divergence: {phase4b['divergence']:+.3f}**\n")

    if harappan:
        md.append("## Harappan KMZ Analysis\n")
        md.append(f"- Total sites: {harappan['total_sites']:,}")
        md.append(f"- Within 50km of GC: {harappan['within_50km']} "
                  f"({harappan['within_50km']/harappan['total_sites']*100:.1f}%)")
        md.append(f"- Median distance: {harappan['median_dist_km']:.0f} km")
        md.append(f"- Enrichment factor (50km): **{harappan['enrichment_50km']:.2f}x**")
        md.append("")
        if harappan["nearest_sites"]:
            md.append("Nearest Harappan sites to GC:")
            md.append("| Site | Distance (km) | Period |")
            md.append("|------|--------------|--------|")
            for s in harappan["nearest_sites"][:10]:
                md.append(f"| {s['name']} | {s['gc_dist_km']:.1f} | {s['period']} |")
            md.append("")

    if wiki_ivc:
        md.append("## Wikipedia IVC Sites\n")
        md.append(f"- Total: {wiki_ivc['total']}")
        md.append(f"- Within 50km: {wiki_ivc['within_50km']}")
        md.append(f"- Median distance: {wiki_ivc['median_dist_km']:.0f} km\n")

    md.append("## Phase 5: GC Segment Through South Asia\n")
    if phase5.get("lat_range"):
        md.append(f"- Latitude range: {phase5['lat_range'][0]} to {phase5['lat_range'][1]}")
        md.append(f"- Longitude range: {phase5['lon_range'][0]} to {phase5['lon_range'][1]}")
    md.append("")

    md.append("## Diagnosis\n")
    md.append(f"**Verdict: {diagnosis['verdict']}**\n")
    md.append(f"{diagnosis['explanation']}\n")
    md.append("### Contributing Factors\n")
    for f in diagnosis["factors"]:
        md.append(f"- **{f['factor']}**: {f['detail']}")
    md.append("")

    md.append("## Plots\n")
    md.append("- `ivc_distances.png` — Major IVC site distances + Harappan histogram")
    md.append("- `coverage_comparison.png` — Database coverage and type distribution")
    md.append("- `gc_south_asia_map.png` — GC segment through South Asia with sites")

    md_path = OUT_DIR / "RESULTS.md"
    with open(md_path, "w") as f:
        f.write("\n".join(md) + "\n")
    print(f"  Saved: {md_path}")


if __name__ == "__main__":
    main()
