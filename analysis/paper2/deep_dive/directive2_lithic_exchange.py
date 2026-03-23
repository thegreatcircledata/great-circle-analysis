#!/usr/bin/env python3
"""
Deep Dive Directive 2: Stone Tool Sourcing / Lithic Exchange Networks

Question: Do stone tool exchange routes show directionality along the Great Circle?

Uses:
  - Pofatu database (Pacific basalt/obsidian sourcing, 7,759+ samples)
  - Published obsidian trade routes (Near East, Mesoamerica, South America)
  - Egyptian quarry-to-pyramid transport routes
  - Known long-distance exchange networks from literature

Tests:
  - For each trade route: compute bearing and compare to local GC bearing
  - Mean offset (0°=along GC, 90°=perpendicular, 45°=random)
  - Distance-weighted analysis (long routes more meaningful)
  - Regional breakdown
  - Temporal analysis where dates available
"""

import json
import math
import os
import sys
import csv
import numpy as np
from scipy import stats as sp_stats
from pathlib import Path

# ── Constants ───────────────────────────────────────────────────────────────
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
OUT_DIR = PROJECT_ROOT / "outputs" / "deep_dive_tests"
OUT_DIR.mkdir(parents=True, exist_ok=True)
POFATU_DIR = PROJECT_ROOT / "data" / "pofatu" / "pofatu-data-1.3" / "dist"

# ── Geometry helpers ────────────────────────────────────────────────────────
def haversine_km(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi/2)**2 + math.cos(phi1)*math.cos(phi2)*math.sin(dlam/2)**2
    return 2 * EARTH_R_KM * math.atan2(math.sqrt(a), math.sqrt(1-a))

def dist_from_gc(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QUARTER_CIRC)

def forward_azimuth(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dlam = math.radians(lon2 - lon1)
    x = math.sin(dlam) * math.cos(phi2)
    y = math.cos(phi1)*math.sin(phi2) - math.sin(phi1)*math.cos(phi2)*math.cos(dlam)
    return math.degrees(math.atan2(x, y)) % 360

def gc_bearing_at(lat, lon):
    bearing_to_pole = forward_azimuth(lat, lon, POLE_LAT, POLE_LON)
    return (bearing_to_pole + 90) % 360

def midpoint(lat1, lon1, lat2, lon2):
    return ((lat1+lat2)/2, (lon1+lon2)/2)


# ── Published lithic exchange routes ────────────────────────────────────────
# Each route: (name, src_lat, src_lon, dst_lat, dst_lon, material, period, region)

PUBLISHED_ROUTES = [
    # ── Near Eastern Obsidian (Cann & Renfrew framework + updates) ──
    # Anatolian sources to Levantine/Mesopotamian sites
    ("Göllü Dağ → Çatalhöyük", 38.25, 34.07, 37.67, 32.83, "obsidian", "7000-5500 BCE", "Near_East"),
    ("Göllü Dağ → Mersin", 38.25, 34.07, 36.80, 34.63, "obsidian", "6000-4000 BCE", "Near_East"),
    ("Göllü Dağ → Jericho", 38.25, 34.07, 31.87, 35.44, "obsidian", "8000-6000 BCE", "Near_East"),
    ("Göllü Dağ → Byblos", 38.25, 34.07, 34.12, 35.65, "obsidian", "5000-3000 BCE", "Near_East"),
    ("Göllü Dağ → Tell Brak", 38.25, 34.07, 36.67, 41.06, "obsidian", "4000-3000 BCE", "Near_East"),
    ("Nemrut Dağ → Jarmo", 38.64, 42.24, 35.55, 44.95, "obsidian", "7000-5000 BCE", "Near_East"),
    ("Nemrut Dağ → Tell es-Sawwan", 38.64, 42.24, 34.20, 43.80, "obsidian", "6000-5000 BCE", "Near_East"),
    ("Nemrut Dağ → Tepe Gawra", 38.64, 42.24, 36.42, 43.27, "obsidian", "5000-3500 BCE", "Near_East"),
    ("Bingöl → Çayönü", 38.88, 40.50, 38.22, 39.72, "obsidian", "8000-6000 BCE", "Near_East"),
    ("Bingöl → Tell Halaf", 38.88, 40.50, 36.82, 40.03, "obsidian", "6000-5000 BCE", "Near_East"),
    ("Melos → Franchthi Cave", 36.68, 24.43, 37.42, 23.13, "obsidian", "11000-3000 BCE", "Mediterranean"),
    ("Melos → Knossos", 36.68, 24.43, 35.30, 25.16, "obsidian", "6000-2000 BCE", "Mediterranean"),
    ("Lipari → Multiple Neolithic sites", 38.47, 14.95, 40.85, 14.25, "obsidian", "5000-3000 BCE", "Mediterranean"),
    ("Sardinia Monte Arci → Corsica", 39.85, 8.87, 42.15, 9.09, "obsidian", "5000-3000 BCE", "Mediterranean"),
    ("Pantelleria → Sicily", 36.78, 11.99, 37.50, 14.00, "obsidian", "5000-3000 BCE", "Mediterranean"),

    # ── Egyptian quarry-to-monument routes (Klemm & Klemm, Aston et al.) ──
    ("Aswan granite → Giza", 24.09, 32.90, 29.98, 31.13, "granite", "2600-2500 BCE", "Egypt"),
    ("Tura limestone → Giza", 29.94, 31.28, 29.98, 31.13, "limestone", "2600-2500 BCE", "Egypt"),
    ("Hatnub alabaster → Amarna", 27.72, 30.97, 27.65, 30.90, "alabaster", "1350 BCE", "Egypt"),
    ("Wadi Hammamat → Thebes", 25.98, 33.60, 25.72, 32.64, "greywacke", "2000-1500 BCE", "Egypt"),
    ("Gebel el-Silsila → Karnak", 24.63, 32.93, 25.72, 32.66, "sandstone", "1500-1000 BCE", "Egypt"),
    ("Aswan granite → Abu Simbel", 24.09, 32.90, 22.34, 31.63, "granite", "1250 BCE", "Egypt"),
    ("Fayum basalt → Giza", 29.30, 30.60, 29.98, 31.13, "basalt", "2600-2500 BCE", "Egypt"),

    # ── Peruvian obsidian (Burger & Glascock, Tripcevich) ──
    ("Quispisisa → Cahuachi (Nazca)", -14.10, -74.35, -14.82, -75.12, "obsidian", "200 BCE-600 CE", "South_America"),
    ("Quispisisa → Chavín de Huántar", -14.10, -74.35, -9.59, -77.18, "obsidian", "1200-500 BCE", "South_America"),
    ("Quispisisa → Wari (Ayacucho)", -14.10, -74.35, -13.06, -74.22, "obsidian", "600-1000 CE", "South_America"),
    ("Alca → Tiwanaku", -15.71, -72.08, -16.55, -68.67, "obsidian", "400-1000 CE", "South_America"),
    ("Alca → Colca Valley sites", -15.71, -72.08, -15.63, -71.60, "obsidian", "1500-500 BCE", "South_America"),
    ("Alca → Arequipa sites", -15.71, -72.08, -16.40, -71.54, "obsidian", "500 BCE-500 CE", "South_America"),

    # ── Mesoamerican obsidian ──
    ("Pachuca → Teotihuacan", 20.10, -98.73, 19.69, -98.84, "obsidian", "100 BCE-700 CE", "Mesoamerica"),
    ("Pachuca → Monte Albán", 20.10, -98.73, 17.04, -96.77, "obsidian", "500 BCE-500 CE", "Mesoamerica"),
    ("El Chayal → Tikal", 14.54, -90.37, 17.22, -89.62, "obsidian", "600 BCE-900 CE", "Mesoamerica"),
    ("El Chayal → Copán", 14.54, -90.37, 14.84, -89.14, "obsidian", "400-900 CE", "Mesoamerica"),
    ("Ixtepeque → Quiriguá", 14.40, -89.68, 15.27, -89.04, "obsidian", "400-900 CE", "Mesoamerica"),

    # ── Easter Island connections ──
    ("Pitcairn (basalt) → Easter Island", -25.07, -130.10, -27.10, -109.30, "basalt", "1200-1500 CE", "Pacific"),

    # ── British/European flint trade ──
    ("Grimes Graves → Wessex", 52.47, 0.77, 51.18, -1.83, "flint", "3000-2000 BCE", "Europe"),
    ("Grand Pressigny → Britain", 46.92, 0.80, 51.50, -1.00, "flint", "2800-2400 BCE", "Europe"),
    ("Krzemionki → Central Europe", 50.97, 21.50, 50.07, 14.44, "flint", "3500-2000 BCE", "Europe"),

    # ── Jade/greenstone ──
    ("Mont Viso (jadeitite) → British Isles", 44.67, 7.09, 51.18, -1.83, "jadeitite", "4500-3500 BCE", "Europe"),
    ("Motagua Valley → Tikal", 15.00, -89.50, 17.22, -89.62, "jade", "600 BCE-900 CE", "Mesoamerica"),
]


def load_pofatu_routes():
    """Load source-to-artefact-find-site routes from Pofatu database."""
    samples_path = POFATU_DIR / "samples.csv"
    if not samples_path.exists():
        print(f"  Pofatu samples.csv not found at {samples_path}")
        return []

    routes = []
    seen = set()

    with open(samples_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cat = row.get("sample_category", "")
            lat_str = row.get("location_latitude", "")
            lon_str = row.get("location_longitude", "")
            site = row.get("site_name", "")
            region = row.get("location_region", "")
            subregion = row.get("location_subregion", "")
            name = row.get("sample_name", "")

            if not lat_str or not lon_str:
                continue
            try:
                lat = float(lat_str)
                lon = float(lon_str)
            except ValueError:
                continue

            # We need source-artefact pairs. Pofatu categorizes samples as
            # SOURCE or ARTEFACT. We need to match artefacts to their source.
            # For now, extract unique source locations and artefact locations
            # per contribution (study).

            if cat in ("SOURCE", "ARTEFACT"):
                contrib = row.get("contribution_name", "")
                key = (contrib, name, lat, lon)
                if key not in seen:
                    seen.add(key)
                    routes.append({
                        "name": name,
                        "category": cat,
                        "lat": lat,
                        "lon": lon,
                        "region": region,
                        "subregion": subregion,
                        "site": site,
                        "contribution": contrib,
                    })

    # Build source-artefact pairs per contribution
    # Group by contribution
    from collections import defaultdict
    by_contrib = defaultdict(lambda: {"sources": [], "artefacts": []})
    for r in routes:
        if r["category"] == "SOURCE":
            by_contrib[r["contribution"]]["sources"].append(r)
        else:
            by_contrib[r["contribution"]]["artefacts"].append(r)

    # For each contribution with both sources and artefacts,
    # create routes from each source to each artefact location
    pofatu_routes = []
    seen_routes = set()

    for contrib, data in by_contrib.items():
        sources = data["sources"]
        artefacts = data["artefacts"]
        if not sources or not artefacts:
            continue

        # Get unique source locations (within 0.1°)
        unique_sources = {}
        for s in sources:
            key = (round(s["lat"], 1), round(s["lon"], 1))
            if key not in unique_sources:
                unique_sources[key] = s

        # Get unique artefact locations
        unique_artefacts = {}
        for a in artefacts:
            key = (round(a["lat"], 1), round(a["lon"], 1))
            if key not in unique_artefacts:
                unique_artefacts[key] = a

        for sk, src in unique_sources.items():
            for ak, art in unique_artefacts.items():
                if sk == ak:
                    continue  # same location
                dist = haversine_km(src["lat"], src["lon"], art["lat"], art["lon"])
                if dist < 10:
                    continue  # too close, probably same site
                route_key = (sk, ak)
                if route_key in seen_routes:
                    continue
                seen_routes.add(route_key)

                pofatu_routes.append((
                    f"{src['subregion']}→{art['subregion']}",
                    src["lat"], src["lon"],
                    art["lat"], art["lon"],
                    "basalt/obsidian",
                    "prehistoric",
                    "Pacific"
                ))

    print(f"  Pofatu: extracted {len(pofatu_routes)} unique source-artefact routes")
    return pofatu_routes


def analyze_routes(all_routes):
    """Compute GC-alignment metrics for all trade routes."""

    analyzed = []
    for route in all_routes:
        name, src_lat, src_lon, dst_lat, dst_lon, material, period, region = route

        # Route bearing (source to destination)
        route_bearing = forward_azimuth(src_lat, src_lon, dst_lat, dst_lon)

        # Route distance
        route_dist = haversine_km(src_lat, src_lon, dst_lat, dst_lon)

        # Midpoint
        mid_lat, mid_lon = midpoint(src_lat, src_lon, dst_lat, dst_lon)

        # GC bearing at midpoint
        gc_bear = gc_bearing_at(mid_lat, mid_lon)

        # Distance of midpoint from GC
        gc_dist = dist_from_gc(mid_lat, mid_lon)

        # Offset: angular difference between route and GC bearing
        # We use mod 180 because trade is bidirectional
        raw_offset = abs(route_bearing - gc_bear)
        if raw_offset > 180:
            raw_offset = 360 - raw_offset
        if raw_offset > 90:
            offset = 180 - raw_offset  # fold to 0-90° range
        else:
            offset = raw_offset

        analyzed.append({
            "name": name,
            "src": (src_lat, src_lon),
            "dst": (dst_lat, dst_lon),
            "material": material,
            "period": period,
            "region": region,
            "route_bearing": round(route_bearing, 1),
            "route_distance_km": round(route_dist, 1),
            "gc_bearing_at_midpoint": round(gc_bear, 1),
            "gc_dist_midpoint_km": round(gc_dist, 1),
            "offset_from_gc_deg": round(offset, 1),
        })

    return analyzed


def run_analysis(analyzed_routes):
    """Statistical analysis of route-GC alignment."""

    results = {
        "meta": {
            "analysis": "Deep Dive Directive 2: Lithic Exchange Network Directionality",
            "date": "2026-03-21",
            "n_total_routes": len(analyzed_routes),
            "pole": {"lat": POLE_LAT, "lon": POLE_LON},
        },
    }

    offsets = np.array([r["offset_from_gc_deg"] for r in analyzed_routes])
    distances = np.array([r["route_distance_km"] for r in analyzed_routes])

    print(f"\n{'='*70}")
    print(f"DIRECTIVE 2: LITHIC EXCHANGE NETWORK DIRECTIONALITY")
    print(f"{'='*70}")
    print(f"Total routes analyzed: {len(analyzed_routes)}")

    # ── Overall statistics ───────────────────────────────────────────
    print(f"\n{'─'*50}")
    print("Overall Route-GC Alignment")
    print(f"{'─'*50}")

    mean_offset = np.mean(offsets)
    median_offset = np.median(offsets)
    std_offset = np.std(offsets)

    # Under random expectation (uniform on [0, 90°]), mean = 45°
    # Test against null of mean=45°
    t_stat, p_ttest = sp_stats.ttest_1samp(offsets, 45.0)

    # Also use a one-sample KS test against uniform on [0, 90]
    offsets_normalized = offsets / 90.0  # scale to [0, 1]
    ks_stat, p_ks = sp_stats.kstest(offsets_normalized, 'uniform')

    results["overall"] = {
        "mean_offset_deg": round(float(mean_offset), 2),
        "median_offset_deg": round(float(median_offset), 2),
        "std_offset_deg": round(float(std_offset), 2),
        "null_expectation_deg": 45.0,
        "t_test_vs_45": {"t": round(float(t_stat), 4), "p": round(float(p_ttest), 6)},
        "ks_test_vs_uniform": {"D": round(float(ks_stat), 4), "p": round(float(p_ks), 6)},
        "interpretation": (
            f"Mean offset = {mean_offset:.1f}°. "
            f"{'<30° suggests along-GC alignment' if mean_offset < 30 else ''}"
            f"{'~45° consistent with random' if 35 < mean_offset < 55 else ''}"
            f"{'>60° suggests perpendicular to GC' if mean_offset > 60 else ''}"
        ),
    }

    print(f"  Mean offset from GC: {mean_offset:.1f}° (random expectation: 45°)")
    print(f"  Median offset: {median_offset:.1f}°, SD: {std_offset:.1f}°")
    print(f"  t-test vs 45°: t={t_stat:.4f}, p={p_ttest:.6f}")
    print(f"  KS test vs uniform: D={ks_stat:.4f}, p={p_ks:.6f}")

    # ── Distance-weighted analysis ───────────────────────────────────
    print(f"\n{'─'*50}")
    print("Distance-weighted analysis")
    print(f"{'─'*50}")

    weighted_mean = np.average(offsets, weights=distances)
    # Correlation between route length and offset
    r_dist_offset, p_dist_offset = sp_stats.pearsonr(distances, offsets)
    rho_dist_offset, p_rho = sp_stats.spearmanr(distances, offsets)

    results["distance_weighted"] = {
        "weighted_mean_offset": round(float(weighted_mean), 2),
        "correlation_dist_vs_offset": {
            "pearson_r": round(float(r_dist_offset), 4),
            "pearson_p": round(float(p_dist_offset), 6),
            "spearman_rho": round(float(rho_dist_offset), 4),
            "spearman_p": round(float(p_rho), 6),
        },
    }

    print(f"  Distance-weighted mean offset: {weighted_mean:.1f}°")
    print(f"  Correlation (route length vs offset):")
    print(f"    Pearson r={r_dist_offset:.4f} (p={p_dist_offset:.4f})")
    print(f"    Spearman rho={rho_dist_offset:.4f} (p={p_rho:.4f})")

    # Short vs long routes
    dist_median = np.median(distances)
    short_mask = distances < dist_median
    long_mask = distances >= dist_median

    short_offset_mean = np.mean(offsets[short_mask])
    long_offset_mean = np.mean(offsets[long_mask])
    t_short_long, p_short_long = sp_stats.ttest_ind(offsets[short_mask], offsets[long_mask])

    results["distance_weighted"]["short_vs_long"] = {
        "distance_cutoff_km": round(float(dist_median), 1),
        "short_routes_mean_offset": round(float(short_offset_mean), 2),
        "long_routes_mean_offset": round(float(long_offset_mean), 2),
        "t_test": {"t": round(float(t_short_long), 4), "p": round(float(p_short_long), 6)},
    }

    print(f"  Short routes (<{dist_median:.0f}km): mean offset = {short_offset_mean:.1f}°")
    print(f"  Long routes (≥{dist_median:.0f}km): mean offset = {long_offset_mean:.1f}°")
    print(f"  Difference: t={t_short_long:.4f}, p={p_short_long:.4f}")

    # ── Regional breakdown ───────────────────────────────────────────
    print(f"\n{'─'*50}")
    print("Regional breakdown")
    print(f"{'─'*50}")

    regions = set(r["region"] for r in analyzed_routes)
    regional_results = {}

    for region in sorted(regions):
        region_routes = [r for r in analyzed_routes if r["region"] == region]
        if len(region_routes) < 2:
            continue
        region_offsets = np.array([r["offset_from_gc_deg"] for r in region_routes])
        region_dists = np.array([r["route_distance_km"] for r in region_routes])

        mean_off = np.mean(region_offsets)
        mean_dist = np.mean(region_dists)

        t, p = sp_stats.ttest_1samp(region_offsets, 45.0) if len(region_offsets) >= 3 else (np.nan, np.nan)

        regional_results[region] = {
            "n_routes": len(region_routes),
            "mean_offset": round(float(mean_off), 2),
            "mean_distance_km": round(float(mean_dist), 1),
            "t_vs_45": round(float(t), 4) if not np.isnan(t) else None,
            "p_vs_45": round(float(p), 6) if not np.isnan(p) else None,
        }
        sig = f" {'***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'}" if not np.isnan(p) else ""
        print(f"  {region}: n={len(region_routes)}, mean offset={mean_off:.1f}°, mean dist={mean_dist:.0f}km{sig}")

    results["regional"] = regional_results

    # ── GC proximity analysis ────────────────────────────────────────
    print(f"\n{'─'*50}")
    print("Routes near vs far from GC")
    print(f"{'─'*50}")

    gc_dists = np.array([r["gc_dist_midpoint_km"] for r in analyzed_routes])

    for gc_thresh in [100, 500, 1000]:
        near_gc = offsets[gc_dists < gc_thresh]
        far_gc = offsets[gc_dists >= gc_thresh]
        if len(near_gc) >= 3 and len(far_gc) >= 3:
            t, p = sp_stats.ttest_ind(near_gc, far_gc)
            print(f"  GC threshold {gc_thresh}km: near n={len(near_gc)} mean={np.mean(near_gc):.1f}° | far n={len(far_gc)} mean={np.mean(far_gc):.1f}° | t={t:.3f} p={p:.4f}")

    # ── Crossing test ────────────────────────────────────────────────
    print(f"\n{'─'*50}")
    print("Route GC-crossing analysis")
    print(f"{'─'*50}")

    n_cross = 0
    for r in analyzed_routes:
        src_gc_dist = dist_from_gc(r["src"][0], r["src"][1])
        dst_gc_dist = dist_from_gc(r["dst"][0], r["dst"][1])
        # A route crosses the GC if source and destination are on opposite sides
        # Simple heuristic: both within threshold but route is long enough
        src_pole_dist = haversine_km(r["src"][0], r["src"][1], POLE_LAT, POLE_LON)
        dst_pole_dist = haversine_km(r["dst"][0], r["dst"][1], POLE_LAT, POLE_LON)
        if (src_pole_dist - QUARTER_CIRC) * (dst_pole_dist - QUARTER_CIRC) < 0:
            n_cross += 1

    results["gc_crossing"] = {
        "n_cross": n_cross,
        "n_total": len(analyzed_routes),
        "fraction": round(n_cross / len(analyzed_routes), 4) if analyzed_routes else 0,
    }
    print(f"  Routes crossing GC: {n_cross}/{len(analyzed_routes)} ({n_cross/len(analyzed_routes)*100:.1f}%)")

    # ── Key sub-test: Easter Island-Pitcairn ──────────────────────────
    print(f"\n{'─'*50}")
    print("Key sub-test: Easter Island-Pitcairn route")
    print(f"{'─'*50}")

    ei_routes = [r for r in analyzed_routes if "Pitcairn" in r["name"] or "Easter" in r["name"]]
    if ei_routes:
        r = ei_routes[0]
        print(f"  Route: {r['name']}")
        print(f"  Route bearing: {r['route_bearing']:.1f}°")
        print(f"  GC bearing at midpoint: {r['gc_bearing_at_midpoint']:.1f}°")
        print(f"  Offset from GC: {r['offset_from_gc_deg']:.1f}°")
        print(f"  Route distance: {r['route_distance_km']:.0f} km")
        print(f"  {'ALIGNED with GC' if r['offset_from_gc_deg'] < 30 else 'NOT aligned with GC' if r['offset_from_gc_deg'] > 60 else 'Moderate alignment with GC'}")

        results["easter_island_pitcairn"] = r
    else:
        print("  No Easter Island-Pitcairn route found in dataset")

    return results, analyzed_routes


def make_plots(analyzed_routes, results):
    """Generate visualization plots."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    offsets = np.array([r["offset_from_gc_deg"] for r in analyzed_routes])
    distances = np.array([r["route_distance_km"] for r in analyzed_routes])

    # ── Histogram of offsets ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.hist(offsets, bins=18, range=(0, 90), alpha=0.7,
            color='steelblue', edgecolor='navy', linewidth=0.5)
    ax.axvline(x=45, color='red', linestyle='--', linewidth=2, label='Random expectation (45°)')
    ax.axvline(x=np.mean(offsets), color='green', linestyle='-', linewidth=2,
               label=f'Observed mean ({np.mean(offsets):.1f}°)')
    ax.set_xlabel("Offset from GC bearing (degrees, 0-90°)")
    ax.set_ylabel("Count")
    ax.set_title(f"Route-GC Offset Distribution (n={len(offsets)})")
    ax.legend(fontsize=9)

    # ── Scatter: distance vs offset ──
    ax = axes[1]
    regions = set(r["region"] for r in analyzed_routes)
    colors = {"Near_East": "gold", "Mediterranean": "orange", "Egypt": "brown",
              "South_America": "green", "Mesoamerica": "red", "Pacific": "blue",
              "Europe": "purple"}
    for region in sorted(regions):
        mask = np.array([r["region"] == region for r in analyzed_routes])
        ax.scatter(distances[mask], offsets[mask], alpha=0.7, label=region,
                  color=colors.get(region, "gray"), s=40, edgecolors="black", linewidth=0.3)
    ax.axhline(y=45, color='red', linestyle='--', alpha=0.5, label='Random (45°)')
    ax.set_xlabel("Route distance (km)")
    ax.set_ylabel("Offset from GC bearing (degrees)")
    ax.set_title("Route Length vs GC Alignment")
    ax.legend(fontsize=7, loc='upper right')

    plt.tight_layout()
    plt.savefig(OUT_DIR / "route_gc_offsets.png", dpi=150, bbox_inches="tight")
    print(f"  Saved: {OUT_DIR / 'route_gc_offsets.png'}")
    plt.close()

    # ── Map with routes ──
    try:
        fig, ax = plt.subplots(figsize=(16, 8))

        # Plot GC as a band
        gc_lats, gc_lons = [], []
        for lon_deg in np.arange(-180, 181, 1):
            # Points on GC are QUARTER_CIRC from pole
            # Solve for lat given lon such that haversine(lat,lon,pole) = QUARTER_CIRC
            # Use iterative approach
            for lat_test in np.arange(-85, 86, 0.5):
                d = haversine_km(lat_test, lon_deg, POLE_LAT, POLE_LON)
                if abs(d - QUARTER_CIRC) < 50:
                    gc_lats.append(lat_test)
                    gc_lons.append(lon_deg)
                    break

        ax.plot(gc_lons, gc_lats, 'r-', linewidth=1.5, alpha=0.5, label='Great Circle')

        # Plot routes
        for r in analyzed_routes:
            color = colors.get(r["region"], "gray")
            alpha = 0.4 if r["offset_from_gc_deg"] > 45 else 0.8
            lw = 0.5 if r["route_distance_km"] < 200 else 1.5
            ax.plot([r["src"][1], r["dst"][1]], [r["src"][0], r["dst"][0]],
                   color=color, alpha=alpha, linewidth=lw)
            ax.scatter([r["src"][1]], [r["src"][0]], c=color, s=15, marker='o', zorder=5)
            ax.scatter([r["dst"][1]], [r["dst"][0]], c=color, s=15, marker='^', zorder=5)

        ax.set_xlim(-180, 180)
        ax.set_ylim(-60, 70)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title("Lithic Exchange Routes and Great Circle")
        ax.legend(fontsize=8)

        plt.tight_layout()
        plt.savefig(OUT_DIR / "lithic_exchange_map.png", dpi=150, bbox_inches="tight")
        print(f"  Saved: {OUT_DIR / 'lithic_exchange_map.png'}")
        plt.close()
    except Exception as e:
        print(f"  Map plot error: {e}")


def main():
    print("Loading published trade routes...")
    all_routes = list(PUBLISHED_ROUTES)
    print(f"  Published routes: {len(PUBLISHED_ROUTES)}")

    print("Loading Pofatu database routes...")
    pofatu = load_pofatu_routes()
    all_routes.extend(pofatu)
    print(f"  Total routes: {len(all_routes)}")

    print("\nAnalyzing route-GC alignment...")
    analyzed = analyze_routes(all_routes)

    results, analyzed = run_analysis(analyzed)

    # ── Interpretation ──
    print(f"\n{'='*70}")
    print("INTERPRETATION")
    print(f"{'='*70}")

    mean_offset = results["overall"]["mean_offset_deg"]
    interpretation = {
        "key_question": "Do stone tool exchange routes show directionality along the Great Circle?",
        "findings": [],
        "caveats": [
            "Published routes are biased toward well-studied regions (Near East, Mediterranean)",
            "Many 'routes' are inferred from source matching, not documented travel paths",
            "Short-distance routes (<100km) are dominated by local geography",
            "The Pofatu database provides coordinates but source-artefact matching requires interpretation",
            "Route bearing is source-to-destination, actual travel paths may differ significantly",
            "Random expectation assumes uniform bearing distribution, which geography can bias",
        ],
    }

    if mean_offset < 35:
        interpretation["findings"].append(
            f"Mean offset ({mean_offset:.1f}°) is below random expectation (45°), "
            "suggesting some tendency for routes to align with the GC"
        )
    elif mean_offset > 55:
        interpretation["findings"].append(
            f"Mean offset ({mean_offset:.1f}°) is above random expectation (45°), "
            "suggesting routes tend to run perpendicular to the GC"
        )
    else:
        interpretation["findings"].append(
            f"Mean offset ({mean_offset:.1f}°) is consistent with random expectation (45°) — "
            "no evidence for GC-aligned directionality"
        )

    results["interpretation"] = interpretation

    for f in interpretation["findings"]:
        print(f"  • {f}")

    # Generate plots
    print(f"\n{'─'*50}")
    print("Generating plots...")
    try:
        make_plots(analyzed, results)
    except Exception as e:
        print(f"  Plot error: {e}")

    # Save results
    def jsonify(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, tuple):
            return list(obj)
        raise TypeError(f"Not JSON serializable: {type(obj)}")

    with open(OUT_DIR / "lithic_exchange.json", "w") as f:
        json.dump(results, f, indent=2, default=jsonify)
    print(f"  Saved: {OUT_DIR / 'lithic_exchange.json'}")

    with open(OUT_DIR / "lithic_exchange_routes.json", "w") as f:
        json.dump(analyzed, f, indent=2, default=jsonify)
    print(f"  Saved: {OUT_DIR / 'lithic_exchange_routes.json'}")

    return results


if __name__ == "__main__":
    main()
