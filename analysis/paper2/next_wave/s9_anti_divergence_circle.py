#!/usr/bin/env python3
"""
Study 9: Anti-Divergence Circle Search
========================================
Find the great circle that maximizes SETTLEMENT clustering while
minimizing MONUMENT clustering -- the inverse of the Alison pattern.
The geography of this "anti-circle" tells us what makes the Alison
circle different.

Phases:
  1. Pole scan for inverse divergence (settlement / monument ratio)
  2. Characterize the top anti-circle (geography + MC confirmation)
  3. Compare anti-circle to Alison circle
  4. Orthogonality test (angular separation between poles)
"""

import sys
import os
import math
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import OrderedDict

sys.stdout = os.fdopen(sys.stdout.fileno(), "w", buffering=1)

# ============================================================
# CONSTANTS
# ============================================================
POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
THRESHOLD_KM = 50
MC_TRIALS = 200

PLEIADES_CSV = "/Users/elliotallan/megalith_site_research/data/pleiades/pleiades-places-latest.csv"
OUTPUT_DIR = "/Users/elliotallan/megalith_site_research/outputs/next_wave/anti_divergence"

np.random.seed(42)

# ============================================================
# TYPE SETS
# ============================================================
MONUMENTAL_TYPES = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus", "shrine",
}
SETTLEMENT_TYPES = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement", "townhouse", "production",
}

# ============================================================
# GEOMETRY
# ============================================================
def haversine_km(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return 2 * EARTH_R_KM * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def dist_from_gc(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QUARTER_CIRC)


def dist_from_gc_vec(site_lats, site_lons, pole_lat, pole_lon):
    """Vectorized distance from great circle defined by pole."""
    lat1r = np.radians(pole_lat)
    lon1r = np.radians(pole_lon)
    lat2r = np.radians(site_lats)
    lon2r = np.radians(site_lons)
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1r) * np.cos(lat2r) * np.sin(dlon / 2) ** 2
    d = 2 * EARTH_R_KM * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return np.abs(d - QUARTER_CIRC)


def angular_separation(lat1, lon1, lat2, lon2):
    """Angular separation between two points on the sphere (degrees)."""
    return haversine_km(lat1, lon1, lat2, lon2) / EARTH_R_KM * (180 / math.pi)


# ============================================================
# DATA LOADING
# ============================================================
def load_pleiades():
    """Load Pleiades CSV and classify sites."""
    monuments = []
    settlements = []
    all_sites = []
    with open(PLEIADES_CSV, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["reprLat"])
                lon = float(row["reprLong"])
            except (ValueError, TypeError):
                continue
            ft_raw = row.get("featureTypes", "")
            feature_types = {t.strip() for t in ft_raw.split(",") if t.strip()}
            is_monument = bool(feature_types & MONUMENTAL_TYPES)
            is_settlement = bool(feature_types & SETTLEMENT_TYPES)
            site = {"lat": lat, "lon": lon, "title": row.get("title", "")}
            all_sites.append(site)
            if is_monument:
                monuments.append(site)
            if is_settlement:
                settlements.append(site)
    return all_sites, monuments, settlements


# ============================================================
# PHASE 1: Pole Scan
# ============================================================
def pole_scan(mon_lats, mon_lons, set_lats, set_lons):
    """Scan pole grid, compute inverse divergence for each."""
    lat_grid = np.arange(0, 91, 5)
    lon_grid = np.arange(-180, 180, 10)
    total_poles = len(lat_grid) * len(lon_grid)
    print(f"Phase 1: Scanning {total_poles} pole candidates...")

    results = []
    done = 0
    report_interval = max(1, total_poles // 20)

    for plat in lat_grid:
        for plon in lon_grid:
            mon_dists = dist_from_gc_vec(mon_lats, mon_lons, plat, plon)
            set_dists = dist_from_gc_vec(set_lats, set_lons, plat, plon)
            mon_on = int(np.sum(mon_dists < THRESHOLD_KM))
            set_on = int(np.sum(set_dists < THRESHOLD_KM))
            inv_div = set_on / max(mon_on, 1)
            results.append({
                "pole_lat": float(plat),
                "pole_lon": float(plon),
                "monument_count": mon_on,
                "settlement_count": set_on,
                "inverse_divergence": inv_div,
            })
            done += 1
            if done % report_interval == 0:
                pct = done / total_poles * 100
                print(f"  ... {pct:.0f}% complete ({done}/{total_poles})")

    results.sort(key=lambda x: x["inverse_divergence"], reverse=True)
    print(f"Phase 1 complete. Top inverse divergence = {results[0]['inverse_divergence']:.2f}")
    return results


# ============================================================
# PHASE 2: Characterize Top Anti-Circle
# ============================================================
def trace_gc_path(pole_lat, pole_lon, n_points=360):
    """Trace points along the great circle defined by a pole.
    The GC is the set of points exactly QUARTER_CIRC from the pole."""
    points = []
    plat_r = math.radians(pole_lat)
    plon_r = math.radians(pole_lon)
    # Angular distance from pole = 90 degrees
    d_rad = math.pi / 2
    for i in range(n_points):
        bearing = math.radians(i * 360 / n_points)
        lat2 = math.asin(
            math.sin(plat_r) * math.cos(d_rad)
            + math.cos(plat_r) * math.sin(d_rad) * math.cos(bearing)
        )
        lon2 = plon_r + math.atan2(
            math.sin(bearing) * math.sin(d_rad) * math.cos(plat_r),
            math.cos(d_rad) - math.sin(plat_r) * math.sin(lat2),
        )
        points.append((math.degrees(lat2), math.degrees(lon2)))
    return points


def geographic_features_along_gc(gc_points):
    """Describe major geographic regions the GC passes through."""
    regions = []
    for lat, lon in gc_points:
        # Rough region classification
        if 25 <= lat <= 45 and -10 <= lon <= 40:
            regions.append("Mediterranean")
        elif 20 <= lat <= 35 and 40 <= lon <= 80:
            regions.append("Middle East / Central Asia")
        elif 10 <= lat <= 45 and 80 <= lon <= 140:
            regions.append("East Asia")
        elif -10 <= lat <= 15 and 90 <= lon <= 160:
            regions.append("Southeast Asia")
        elif 30 <= lat <= 55 and -130 <= lon <= -60:
            regions.append("North America")
        elif -40 <= lat <= 15 and -80 <= lon <= -30:
            regions.append("South America")
        elif -40 <= lat <= 40 and 10 <= lon <= 55:
            regions.append("Africa")
        elif 45 <= lat <= 72 and -10 <= lon <= 60:
            regions.append("Northern Europe / Russia")
        elif -50 <= lat <= -10 and 110 <= lon <= 180:
            regions.append("Oceania")
    # Deduplicate preserving order
    seen = set()
    unique = []
    for r in regions:
        if r not in seen:
            seen.add(r)
            unique.append(r)
    return unique


def mc_random_pole(site_lats, site_lons, observed_on, n_trials=MC_TRIALS):
    """MC test: compare observed count to counts from random great circles.
    Generates random poles and counts sites within THRESHOLD_KM of each
    random great circle. Returns Z-score and p-value."""
    counts = np.zeros(n_trials)
    for i in range(n_trials):
        # Random pole: uniform on hemisphere (lat 0-90)
        rand_lat = np.degrees(np.arcsin(np.random.uniform(0, 1)))
        rand_lon = np.random.uniform(-180, 180)
        dists = dist_from_gc_vec(site_lats, site_lons, rand_lat, rand_lon)
        counts[i] = np.sum(dists < THRESHOLD_KM)
    mean_mc = float(np.mean(counts))
    std_mc = float(np.std(counts))
    z = (observed_on - mean_mc) / std_mc if std_mc > 0 else 0.0
    p = float(np.mean(counts >= observed_on))
    return z, p, mean_mc, std_mc


# ============================================================
# PHASE 4: Orthogonality
# ============================================================
def compute_orthogonality(alison_lat, alison_lon, anti_lat, anti_lon):
    """Compute angular separation between two poles."""
    sep = angular_separation(alison_lat, alison_lon, anti_lat, anti_lon)
    return sep


# ============================================================
# PLOTS
# ============================================================
def plot_pole_heatmap(results, output_path):
    """Heatmap of inverse divergence across pole grid."""
    lats = sorted(set(r["pole_lat"] for r in results))
    lons = sorted(set(r["pole_lon"] for r in results))
    grid = np.zeros((len(lats), len(lons)))
    lat_idx = {v: i for i, v in enumerate(lats)}
    lon_idx = {v: i for i, v in enumerate(lons)}
    for r in results:
        grid[lat_idx[r["pole_lat"]], lon_idx[r["pole_lon"]]] = r["inverse_divergence"]

    fig, ax = plt.subplots(figsize=(14, 6))
    im = ax.imshow(
        grid, origin="lower", aspect="auto",
        extent=[min(lons), max(lons), min(lats), max(lats)],
        cmap="YlOrRd",
    )
    plt.colorbar(im, ax=ax, label="Inverse Divergence (settlement / monument)")
    ax.set_xlabel("Pole Longitude")
    ax.set_ylabel("Pole Latitude")
    ax.set_title("Study 9: Anti-Divergence Pole Scan\n(Higher = more settlements, fewer monuments on circle)")

    # Mark top pole
    top = results[0]
    ax.plot(top["pole_lon"], top["pole_lat"], "k*", markersize=15, label="Top anti-pole")
    # Mark Alison pole
    ax.plot(POLE_LON, POLE_LAT, "c^", markersize=12, label="Alison pole")
    ax.legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_gc_comparison(anti_pole, output_path):
    """Plot both great circles on a world map projection."""
    alison_pts = trace_gc_path(POLE_LAT, POLE_LON)
    anti_pts = trace_gc_path(anti_pole["pole_lat"], anti_pole["pole_lon"])

    fig, ax = plt.subplots(figsize=(14, 7))
    # Simple coastline approximation -- just plot the circles
    al_lats = [p[0] for p in alison_pts]
    al_lons = [p[1] for p in alison_pts]
    an_lats = [p[0] for p in anti_pts]
    an_lons = [p[1] for p in anti_pts]

    ax.scatter(al_lons, al_lats, s=1, c="blue", alpha=0.5, label="Alison circle")
    ax.scatter(an_lons, an_lats, s=1, c="red", alpha=0.5, label="Anti-divergence circle")

    ax.plot(POLE_LON, POLE_LAT, "b^", markersize=10, label="Alison pole")
    ax.plot(anti_pole["pole_lon"], anti_pole["pole_lat"], "r*", markersize=12, label="Anti pole")

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Study 9: Alison Circle vs Anti-Divergence Circle")
    ax.legend(loc="lower left")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_divergence_comparison(alison_stats, anti_stats, output_path):
    """Bar chart comparing monument/settlement counts for both circles."""
    labels = ["Monuments on circle", "Settlements on circle", "Divergence ratio"]
    alison_vals = [
        alison_stats["monument_count"],
        alison_stats["settlement_count"],
        alison_stats["divergence"],
    ]
    anti_vals = [
        anti_stats["monument_count"],
        anti_stats["settlement_count"],
        anti_stats["inverse_divergence"],
    ]

    x = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))
    bars1 = ax.bar(x - width / 2, alison_vals, width, label="Alison circle", color="steelblue")
    bars2 = ax.bar(x + width / 2, anti_vals, width, label="Anti-divergence circle", color="indianred")

    ax.set_ylabel("Count / Ratio")
    ax.set_title("Study 9: Alison vs Anti-Divergence Circle Comparison")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    # Value labels
    for bar in bars1:
        h = bar.get_height()
        ax.annotate(f"{h:.1f}", xy=(bar.get_x() + bar.get_width() / 2, h),
                    xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=9)
    for bar in bars2:
        h = bar.get_height()
        ax.annotate(f"{h:.1f}", xy=(bar.get_x() + bar.get_width() / 2, h),
                    xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  Saved: {output_path}")


# ============================================================
# MAIN
# ============================================================
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # --- Load data ---
    print("Loading Pleiades data...")
    all_sites, monuments, settlements = load_pleiades()
    print(f"  Total sites: {len(all_sites)}")
    print(f"  Monuments: {len(monuments)}")
    print(f"  Settlements: {len(settlements)}")

    mon_lats = np.array([s["lat"] for s in monuments])
    mon_lons = np.array([s["lon"] for s in monuments])
    set_lats = np.array([s["lat"] for s in settlements])
    set_lons = np.array([s["lon"] for s in settlements])

    # === PHASE 1: Pole Scan ===
    scan_results = pole_scan(mon_lats, mon_lons, set_lats, set_lons)
    top_10 = scan_results[:10]
    print("\nTop 10 Anti-Divergence Poles:")
    print(f"  {'Rank':<5} {'Lat':>6} {'Lon':>7} {'Mon':>5} {'Set':>5} {'InvDiv':>8}")
    for i, r in enumerate(top_10):
        print(f"  {i+1:<5} {r['pole_lat']:>6.1f} {r['pole_lon']:>7.1f} "
              f"{r['monument_count']:>5} {r['settlement_count']:>5} "
              f"{r['inverse_divergence']:>8.2f}")

    # === PHASE 2: Characterize Top Anti-Circle ===
    print("\nPhase 2: Characterizing top anti-circle...")
    anti_pole = top_10[0]
    gc_points = trace_gc_path(anti_pole["pole_lat"], anti_pole["pole_lon"])
    regions = geographic_features_along_gc(gc_points)
    print(f"  Anti-circle pole: ({anti_pole['pole_lat']}, {anti_pole['pole_lon']})")
    print(f"  Regions crossed: {', '.join(regions)}")

    # MC test for settlement enrichment (compare to random great circles)
    anti_set_dists = dist_from_gc_vec(set_lats, set_lons, anti_pole["pole_lat"], anti_pole["pole_lon"])
    observed_set_on = int(np.sum(anti_set_dists < THRESHOLD_KM))
    z_set, p_set, mean_set, std_set = mc_random_pole(
        set_lats, set_lons, observed_set_on, MC_TRIALS
    )
    print(f"  Settlement MC: observed={observed_set_on}, mean={mean_set:.1f}, "
          f"std={std_set:.1f}, Z={z_set:.2f}, p={p_set:.4f}")

    # MC for monument count (compare to random great circles -- should be low)
    anti_mon_dists = dist_from_gc_vec(mon_lats, mon_lons, anti_pole["pole_lat"], anti_pole["pole_lon"])
    observed_mon_on = int(np.sum(anti_mon_dists < THRESHOLD_KM))
    z_mon, p_mon_high, mean_mon, std_mon = mc_random_pole(
        mon_lats, mon_lons, observed_mon_on, MC_TRIALS
    )
    p_mon_low = 1.0 - p_mon_high  # fraction of random circles with FEWER monuments
    print(f"  Monument MC: observed={observed_mon_on}, mean={mean_mon:.1f}, "
          f"std={std_mon:.1f}, Z={z_mon:.2f}, p(low)={p_mon_low:.4f}")

    # === PHASE 3: Compare to Alison Circle ===
    print("\nPhase 3: Comparing to Alison circle...")
    alison_mon_dists = dist_from_gc_vec(mon_lats, mon_lons, POLE_LAT, POLE_LON)
    alison_set_dists = dist_from_gc_vec(set_lats, set_lons, POLE_LAT, POLE_LON)
    alison_mon_on = int(np.sum(alison_mon_dists < THRESHOLD_KM))
    alison_set_on = int(np.sum(alison_set_dists < THRESHOLD_KM))
    alison_divergence = alison_mon_on / max(alison_set_on, 1)

    alison_stats = {
        "monument_count": alison_mon_on,
        "settlement_count": alison_set_on,
        "divergence": round(alison_divergence, 4),
        "pole_lat": POLE_LAT,
        "pole_lon": POLE_LON,
    }

    anti_stats = {
        "monument_count": anti_pole["monument_count"],
        "settlement_count": anti_pole["settlement_count"],
        "inverse_divergence": round(anti_pole["inverse_divergence"], 4),
        "pole_lat": anti_pole["pole_lat"],
        "pole_lon": anti_pole["pole_lon"],
    }

    print(f"  Alison: {alison_mon_on} mon, {alison_set_on} set, div={alison_divergence:.4f}")
    print(f"  Anti:   {anti_pole['monument_count']} mon, {anti_pole['settlement_count']} set, "
          f"inv_div={anti_pole['inverse_divergence']:.4f}")

    # === PHASE 4: Orthogonality ===
    print("\nPhase 4: Orthogonality test...")
    sep = compute_orthogonality(POLE_LAT, POLE_LON, anti_pole["pole_lat"], anti_pole["pole_lon"])
    is_orthogonal = 80 <= sep <= 100
    print(f"  Angular separation: {sep:.2f} degrees")
    print(f"  Orthogonal (80-100 deg): {is_orthogonal}")

    # === PLOTS ===
    print("\nGenerating plots...")
    plot_pole_heatmap(scan_results, os.path.join(OUTPUT_DIR, "pole_heatmap.png"))
    plot_gc_comparison(anti_pole, os.path.join(OUTPUT_DIR, "gc_comparison.png"))
    plot_divergence_comparison(alison_stats, anti_stats, os.path.join(OUTPUT_DIR, "divergence_comparison.png"))

    # === SAVE RESULTS ===
    results = OrderedDict()
    results["study"] = "S9: Anti-Divergence Circle Search"
    results["description"] = (
        "Find the great circle maximizing settlement clustering while "
        "minimizing monument clustering -- the inverse of the Alison pattern."
    )
    results["phase1_pole_scan"] = {
        "total_poles_scanned": len(scan_results),
        "top_10": top_10,
    }
    results["phase2_anti_circle"] = {
        "pole_lat": anti_pole["pole_lat"],
        "pole_lon": anti_pole["pole_lon"],
        "monument_count": anti_pole["monument_count"],
        "settlement_count": anti_pole["settlement_count"],
        "inverse_divergence": round(anti_pole["inverse_divergence"], 4),
        "regions_crossed": regions,
        "settlement_mc": {
            "observed": observed_set_on,
            "mean": round(mean_set, 2),
            "std": round(std_set, 2),
            "z_score": round(z_set, 2),
            "p_value": round(p_set, 4),
            "trials": MC_TRIALS,
        },
        "monument_mc": {
            "observed": observed_mon_on,
            "mean": round(mean_mon, 2),
            "std": round(std_mon, 2),
            "z_score": round(z_mon, 2),
            "p_low": round(p_mon_low, 4),
            "trials": MC_TRIALS,
        },
    }
    results["phase3_comparison"] = {
        "alison_circle": alison_stats,
        "anti_circle": anti_stats,
        "ratio_of_ratios": round(
            anti_pole["inverse_divergence"] / max(alison_divergence, 0.001), 4
        ),
    }
    results["phase4_orthogonality"] = {
        "angular_separation_deg": round(sep, 2),
        "is_orthogonal_80_100": is_orthogonal,
    }

    json_path = os.path.join(OUTPUT_DIR, "results.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Saved: {json_path}")

    # === RESULTS.md ===
    md_path = os.path.join(OUTPUT_DIR, "RESULTS.md")
    with open(md_path, "w") as f:
        f.write("# Study 9: Anti-Divergence Circle Search\n\n")

        f.write("## Overview\n")
        f.write("This study searches for the great circle that maximizes settlement clustering\n")
        f.write("while minimizing monument clustering -- the inverse of the Alison pattern.\n")
        f.write("By identifying the 'anti-circle', we can understand what geographic and\n")
        f.write("archaeological properties distinguish the Alison circle.\n\n")

        f.write("## Phase 1: Pole Scan Results\n\n")
        f.write(f"Scanned {len(scan_results)} pole candidates (lat 0-90, step 5; lon -180 to 180, step 10).\n\n")
        f.write("### Top 10 Anti-Divergence Poles\n\n")
        f.write("| Rank | Lat | Lon | Monuments | Settlements | Inv. Divergence |\n")
        f.write("|------|-----|-----|-----------|-------------|------------------|\n")
        for i, r in enumerate(top_10):
            f.write(f"| {i+1} | {r['pole_lat']:.1f} | {r['pole_lon']:.1f} | "
                    f"{r['monument_count']} | {r['settlement_count']} | "
                    f"{r['inverse_divergence']:.2f} |\n")
        f.write("\n")

        f.write("## Phase 2: Anti-Circle Characterization\n\n")
        f.write(f"**Anti-circle pole:** ({anti_pole['pole_lat']}, {anti_pole['pole_lon']})\n\n")
        f.write(f"**Regions crossed:** {', '.join(regions)}\n\n")
        f.write(f"- Settlements within 50km: {anti_pole['settlement_count']}\n")
        f.write(f"- Monuments within 50km: {anti_pole['monument_count']}\n")
        f.write(f"- Inverse divergence: {anti_pole['inverse_divergence']:.2f}\n\n")

        f.write("### Monte Carlo Validation\n\n")
        f.write(f"Settlement enrichment: Z = {z_set:.2f}, p = {p_set:.4f} ({MC_TRIALS} trials)\n\n")
        f.write(f"Monument depletion: Z = {z_mon:.2f}, p(low) = {p_mon_low:.4f} ({MC_TRIALS} trials)\n\n")

        f.write("## Phase 3: Alison vs Anti-Circle Comparison\n\n")
        f.write("| Metric | Alison Circle | Anti-Circle |\n")
        f.write("|--------|--------------|-------------|\n")
        f.write(f"| Pole (lat, lon) | ({POLE_LAT}, {POLE_LON}) | "
                f"({anti_pole['pole_lat']}, {anti_pole['pole_lon']}) |\n")
        f.write(f"| Monuments within 50km | {alison_mon_on} | {anti_pole['monument_count']} |\n")
        f.write(f"| Settlements within 50km | {alison_set_on} | {anti_pole['settlement_count']} |\n")
        f.write(f"| Mon/Set ratio | {alison_divergence:.4f} | "
                f"{anti_pole['monument_count']/max(anti_pole['settlement_count'],1):.4f} |\n")
        f.write(f"| Set/Mon ratio | {alison_set_on/max(alison_mon_on,1):.4f} | "
                f"{anti_pole['inverse_divergence']:.4f} |\n")
        f.write("\n")

        f.write("## Phase 4: Orthogonality Test\n\n")
        f.write(f"Angular separation between Alison pole and anti-circle pole: **{sep:.2f} degrees**\n\n")
        if is_orthogonal:
            f.write("The poles are approximately orthogonal (80-100 deg). This means the two\n")
            f.write("circles are nearly perpendicular, which would be a geometric constraint\n")
            f.write("difficult to explain by chance alone.\n\n")
        else:
            f.write(f"The poles are NOT orthogonal ({sep:.1f} deg, outside 80-100 range).\n")
            f.write("The two circles are not geometrically constrained to be perpendicular.\n\n")

        f.write("## Plots\n\n")
        f.write("- `pole_heatmap.png` -- Inverse divergence across all scanned poles\n")
        f.write("- `gc_comparison.png` -- Alison circle vs anti-divergence circle\n")
        f.write("- `divergence_comparison.png` -- Bar chart comparison of counts and ratios\n")

    print(f"  Saved: {md_path}")
    print("\nStudy 9 complete.")


if __name__ == "__main__":
    main()
