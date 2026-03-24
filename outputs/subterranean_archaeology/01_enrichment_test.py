#!/usr/bin/env python3
"""
Directive 12, Script 01: Subterranean Site Enrichment Test
==========================================================
Tests whether subterranean archaeological sites cluster along the Great Circle
using distribution-matched Monte Carlo (10,000 random great circles).

Key test: monumental-vs-utilitarian divergence
  - Hypothesis: rock_cut_monumental + artificial_complex should show enrichment
    (like surface monuments), while utilitarian sites should NOT
  - This replicates the surface monument-settlement divergence underground

Output:
  - enrichment_by_type.json
  - underground_divergence.json
  - subterranean_enrichment_map.png
"""

import csv
import json
import math
import os
import sys
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Constants ──────────────────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
R_EARTH = 6371.0
QUARTER_CIRC = R_EARTH * math.pi / 2
N_MC = 10000
THRESHOLDS_KM = [50, 100, 200, 500]

BASE = os.path.dirname(os.path.abspath(__file__))
np.random.seed(42)


# ── Geo utilities ──────────────────────────────────────────────────────
def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlat, dlon = la2 - la1, lo2 - lo1
    a = math.sin(dlat / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin(dlon / 2) ** 2
    return R_EARTH * 2 * math.asin(math.sqrt(min(1.0, a)))


def gc_distance_from_pole(lat, lon, pole_lat, pole_lon):
    d = haversine_km(pole_lat, pole_lon, lat, lon)
    return abs(d - QUARTER_CIRC)


def random_pole():
    """Generate a uniformly random point on the sphere."""
    z = np.random.uniform(-1, 1)
    theta = np.random.uniform(0, 2 * np.pi)
    lat = math.degrees(math.asin(z))
    lon = math.degrees(theta) - 180
    return lat, lon


def get_circle_points(pole_lat, pole_lon, n=360):
    """Get n points on the great circle defined by a pole."""
    points = []
    plat = math.radians(pole_lat)
    plon = math.radians(pole_lon)
    for i in range(n):
        bearing = math.radians(i * 360 / n)
        lat = math.asin(
            math.sin(plat) * math.cos(math.pi / 2)
            + math.cos(plat) * math.sin(math.pi / 2) * math.cos(bearing)
        )
        lon = plon + math.atan2(
            math.sin(bearing) * math.sin(math.pi / 2) * math.cos(plat),
            math.cos(math.pi / 2) - math.sin(plat) * math.sin(lat)
        )
        points.append((math.degrees(lat), math.degrees(lon)))
    return points


# ── Load data ──────────────────────────────────────────────────────────
def load_sites():
    path = os.path.join(BASE, "subterranean_sites_master.csv")
    sites = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sites.append({
                "name": row["name"],
                "lat": float(row["lat"]),
                "lon": float(row["lon"]),
                "subtype": row["subtype"],
                "gc_distance_km": float(row["gc_distance_km"]),
            })
    return sites


# ── Enrichment test ────────────────────────────────────────────────────
def count_within(sites, pole_lat, pole_lon, threshold_km):
    """Count sites within threshold_km of the great circle defined by pole."""
    count = 0
    for s in sites:
        d = gc_distance_from_pole(s["lat"], s["lon"], pole_lat, pole_lon)
        if d <= threshold_km:
            count += 1
    return count


def run_enrichment(sites, label="all"):
    """Run Monte Carlo enrichment test for a set of sites."""
    print(f"\n--- Enrichment: {label} ({len(sites)} sites) ---")

    results = {}
    for thresh in THRESHOLDS_KM:
        actual = count_within(sites, POLE_LAT, POLE_LON, thresh)

        mc_counts = np.zeros(N_MC)
        for i in range(N_MC):
            plat, plon = random_pole()
            mc_counts[i] = count_within(sites, plat, plon, thresh)

        mc_mean = np.mean(mc_counts)
        mc_std = np.std(mc_counts)
        z_score = (actual - mc_mean) / mc_std if mc_std > 0 else 0
        p_value = np.mean(mc_counts >= actual)

        results[f"{thresh}km"] = {
            "actual": int(actual),
            "mc_mean": round(float(mc_mean), 2),
            "mc_std": round(float(mc_std), 2),
            "z_score": round(float(z_score), 2),
            "p_value": round(float(p_value), 4),
            "significant": bool(p_value < 0.05),
        }

        sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
        print(f"  {thresh:>4d} km: actual={actual:>4d}  expected={mc_mean:>7.1f}±{mc_std:.1f}  "
              f"Z={z_score:>6.2f}  p={p_value:.4f} {sig}")

    return results


def main():
    print("=" * 70)
    print("SUBTERRANEAN SITE ENRICHMENT TEST")
    print(f"Monte Carlo: {N_MC} random great circles")
    print("=" * 70)

    sites = load_sites()
    print(f"Loaded {len(sites)} subterranean sites")

    # ── Test 1: All subterranean sites ──
    all_results = run_enrichment(sites, "All subterranean")

    # ── Test 2: By subtype ──
    subtypes = {}
    for s in sites:
        st = s["subtype"]
        if st not in subtypes:
            subtypes[st] = []
        subtypes[st].append(s)

    type_results = {}
    for st in ["natural_cave", "modified_cave", "rock_cut_monumental",
               "artificial_complex", "utilitarian"]:
        if st in subtypes:
            type_results[st] = run_enrichment(subtypes[st], st)
        else:
            print(f"\n  {st}: no sites found")

    # ── Test 3: Monumental vs Utilitarian divergence ──
    print("\n" + "=" * 70)
    print("MONUMENTAL vs UTILITARIAN DIVERGENCE")
    print("=" * 70)

    monumental_sites = [s for s in sites if s["subtype"] in
                        ("rock_cut_monumental", "artificial_complex")]
    utilitarian_sites = [s for s in sites if s["subtype"] == "utilitarian"]

    print(f"\nMonumental (rock-cut + artificial): {len(monumental_sites)} sites")
    print(f"Utilitarian (mines, aqueducts):     {len(utilitarian_sites)} sites")

    mon_results = run_enrichment(monumental_sites, "MONUMENTAL")
    util_results = run_enrichment(utilitarian_sites, "UTILITARIAN")

    divergence = {}
    for thresh in THRESHOLDS_KM:
        key = f"{thresh}km"
        z_mon = mon_results[key]["z_score"]
        z_util = util_results[key]["z_score"]
        d = z_mon - z_util
        divergence[key] = {
            "z_monumental": z_mon,
            "z_utilitarian": z_util,
            "divergence_D": round(d, 2),
            "monumental_enriched": z_mon > 1.96,
            "utilitarian_enriched": z_util > 1.96,
            "pattern_matches_surface": z_mon > z_util,
        }
        print(f"\n  {thresh} km: Z_mon={z_mon:+.2f}  Z_util={z_util:+.2f}  "
              f"D={d:+.2f}  {'✓ matches surface pattern' if z_mon > z_util else '✗ does not match'}")

    # ── Save results ──
    enrichment_output = {
        "analysis": "Subterranean site enrichment by type",
        "n_sites": len(sites),
        "n_monte_carlo": N_MC,
        "all_subterranean": all_results,
        "by_type": type_results,
    }
    with open(os.path.join(BASE, "enrichment_by_type.json"), "w") as f:
        json.dump(enrichment_output, f, indent=2)
    print(f"\nSaved enrichment_by_type.json")

    divergence_output = {
        "analysis": "Monumental vs utilitarian underground divergence",
        "n_monumental": len(monumental_sites),
        "n_utilitarian": len(utilitarian_sites),
        "monumental_enrichment": mon_results,
        "utilitarian_enrichment": util_results,
        "divergence": divergence,
    }
    with open(os.path.join(BASE, "underground_divergence.json"), "w") as f:
        json.dump(divergence_output, f, indent=2)
    print(f"Saved underground_divergence.json")

    # ── Visualization ──
    make_map(sites)


def make_map(sites):
    """World map with subterranean sites colored by type, circle overlaid."""
    fig, ax = plt.subplots(figsize=(18, 9))

    # Draw circle
    circle_pts = get_circle_points(POLE_LAT, POLE_LON, 720)
    # Sort segments to avoid wrapping artifacts
    segments = []
    current = [circle_pts[0]]
    for i in range(1, len(circle_pts)):
        if abs(circle_pts[i][1] - circle_pts[i - 1][1]) > 180:
            segments.append(current)
            current = [circle_pts[i]]
        else:
            current.append(circle_pts[i])
    segments.append(current)
    for seg in segments:
        lats = [p[0] for p in seg]
        lons = [p[1] for p in seg]
        ax.plot(lons, lats, "r-", linewidth=1.5, alpha=0.6, zorder=2)

    # Color map for subtypes
    colors = {
        "natural_cave": "#4CAF50",       # green
        "modified_cave": "#2196F3",      # blue
        "rock_cut_monumental": "#FF9800", # orange
        "artificial_complex": "#F44336",  # red
        "utilitarian": "#9E9E9E",        # gray
    }
    labels = {
        "natural_cave": "Natural cave",
        "modified_cave": "Modified cave",
        "rock_cut_monumental": "Rock-cut monumental",
        "artificial_complex": "Artificial complex",
        "utilitarian": "Utilitarian",
    }

    # Plot sites (utilitarian and natural_cave first so monumental is on top)
    plot_order = ["utilitarian", "natural_cave", "modified_cave",
                  "rock_cut_monumental", "artificial_complex"]
    for st in plot_order:
        subset = [s for s in sites if s["subtype"] == st]
        if not subset:
            continue
        lats = [s["lat"] for s in subset]
        lons = [s["lon"] for s in subset]
        ax.scatter(lons, lats, c=colors.get(st, "#999"),
                   s=8, alpha=0.5, zorder=3, label=f"{labels.get(st, st)} ({len(subset)})")

    # Highlight sites within 200km
    near = [s for s in sites if s["gc_distance_km"] <= 200
            and s["subtype"] in ("rock_cut_monumental", "artificial_complex")]
    if near:
        lats = [s["lat"] for s in near]
        lons = [s["lon"] for s in near]
        ax.scatter(lons, lats, c="none", edgecolors="red", s=60,
                   linewidths=1.5, zorder=4, label=f"Monumental ≤200km ({len(near)})")

    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 75)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Subterranean Archaeological Sites & The Great Circle\n"
                 f"(N={len(sites)} sites, red line = Great Circle)")
    ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.2)

    fig.tight_layout()
    path = os.path.join(BASE, "subterranean_enrichment_map.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
