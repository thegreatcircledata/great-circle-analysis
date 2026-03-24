#!/usr/bin/env python3
"""
Generate visualizations for the oral tradition mapping analysis:
1. mythology_corridor_map.png — world map with motif density and Great Circle
2. dispersal_route_overlay.png — d'Huy's Cosmic Hunt reconstruction + Great Circle
3. enrichment_summary.png — bar chart of enrichment ratios
"""

import csv
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oral_tradition_analysis import (
    load_data, TARGET_MOTIFS, gc_distance, haversine_km,
    POLE_LAT, POLE_LON, EARTH_R_KM, QUARTER_CIRC, CORRIDOR_WIDTH_KM
)

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))


def great_circle_points(pole_lat, pole_lon, n=360):
    """Generate points along the great circle defined by its pole."""
    plat = math.radians(pole_lat)
    plon = math.radians(pole_lon)
    points = []
    for i in range(n):
        bearing = math.radians(i * 360.0 / n)
        d = math.pi / 2
        lat2 = math.asin(math.sin(plat) * math.cos(d) +
                          math.cos(plat) * math.sin(d) * math.cos(bearing))
        lon2 = plon + math.atan2(math.sin(bearing) * math.sin(d) * math.cos(plat),
                                  math.cos(d) - math.sin(plat) * math.sin(lat2))
        points.append((math.degrees(lon2), math.degrees(lat2)))
    return points


def draw_great_circle(ax, pole_lat, pole_lon, color="red", lw=1.5, alpha=0.7, label="Great Circle"):
    """Draw the great circle, handling the antimeridian wrap."""
    pts = great_circle_points(pole_lat, pole_lon, 720)
    # Split at antimeridian
    segments = []
    current = [pts[0]]
    for i in range(1, len(pts)):
        if abs(pts[i][0] - pts[i-1][0]) > 180:
            segments.append(current)
            current = [pts[i]]
        else:
            current.append(pts[i])
    segments.append(current)
    for j, seg in enumerate(segments):
        xs, ys = zip(*seg)
        ax.plot(xs, ys, color=color, linewidth=lw, alpha=alpha,
                label=label if j == 0 else None, zorder=3)


def load_coastlines_simple():
    """Return a very simple coastline approximation using matplotlib's built-in data."""
    # We'll just draw continent outlines from a simple dataset
    # Using a minimal approach since we don't have cartopy
    return None


def fig1_corridor_map(groups):
    """World map with ethnic groups colored by composite motif score, Great Circle overlay."""
    fig, ax = plt.subplots(1, 1, figsize=(16, 8))
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_facecolor("#e8e8e8")
    ax.set_xlabel("Longitude", fontsize=10)
    ax.set_ylabel("Latitude", fontsize=10)
    ax.set_title("Mythology Motif Density Along the Great Circle\n"
                 "(Berezkin Catalogue, 926 ethnic groups)", fontsize=13, fontweight="bold")

    # Grid
    ax.grid(True, alpha=0.3, linestyle="--")

    # Compute composite score for each group
    scores = []
    for g in groups:
        score = 0
        for cat_info in TARGET_MOTIFS.values():
            if any(g["motifs"].get(c, 0) == 1 for c in cat_info["codes"]):
                score += 1
        scores.append(score)

    max_score = max(scores) if scores else 1

    # Separate on/off corridor
    on_lons, on_lats, on_scores = [], [], []
    off_lons, off_lats, off_scores = [], [], []
    for g, s in zip(groups, scores):
        if g["gc_dist"] <= CORRIDOR_WIDTH_KM:
            on_lons.append(g["lon"])
            on_lats.append(g["lat"])
            on_scores.append(s)
        else:
            off_lons.append(g["lon"])
            off_lats.append(g["lat"])
            off_scores.append(s)

    # Plot off-corridor groups
    sc_off = ax.scatter(off_lons, off_lats, c=off_scores, cmap="YlOrRd",
                         s=20, alpha=0.5, vmin=0, vmax=max_score,
                         edgecolors="none", zorder=2, label="Off-corridor")

    # Plot on-corridor groups (larger, with edge)
    sc_on = ax.scatter(on_lons, on_lats, c=on_scores, cmap="YlOrRd",
                        s=60, alpha=0.9, vmin=0, vmax=max_score,
                        edgecolors="black", linewidths=0.8, zorder=4,
                        label="On-corridor (≤200km)")

    # Great Circle
    draw_great_circle(ax, POLE_LAT, POLE_LON)

    # Colorbar
    cbar = plt.colorbar(sc_off, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label(f"Target motif categories present (out of {len(TARGET_MOTIFS)})", fontsize=9)

    ax.legend(loc="lower left", fontsize=9)

    path = os.path.join(OUTPUT_DIR, "mythology_corridor_map.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Written: {path}")


def fig2_enrichment_chart():
    """Bar chart of enrichment ratios with confidence intervals."""
    path_in = os.path.join(OUTPUT_DIR, "enrichment_by_motif.json")
    with open(path_in) as f:
        data = json.load(f)

    cats = sorted(data.keys(), key=lambda k: data[k]["enrichment_ratio"], reverse=True)
    enrichments = [data[c]["enrichment_ratio"] for c in cats]
    mc_means = [data[c]["mc_mean_enrichment"] for c in cats]
    mc_stds = [data[c]["mc_std_enrichment"] for c in cats]
    p_values = [data[c]["mc_p_value"] for c in cats]

    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(cats))
    bars = ax.bar(x, enrichments, color=["#2ca02c" if e > 1 else "#d62728" for e in enrichments],
                  alpha=0.7, zorder=2, edgecolor="black", linewidth=0.5)

    # MC baseline ± 2σ
    ax.errorbar(x, mc_means, yerr=[2*s for s in mc_stds], fmt="none",
                ecolor="gray", capsize=4, zorder=3, label="MC mean ± 2σ")
    ax.scatter(x, mc_means, color="gray", s=30, zorder=4)

    # Reference line at 1.0
    ax.axhline(1.0, color="black", linestyle="--", alpha=0.5, linewidth=1)

    # Labels
    labels = [c.replace("_", "\n") for c in cats]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8, ha="center")
    ax.set_ylabel("Enrichment Ratio (on-corridor / off-corridor)", fontsize=10)
    ax.set_title("Mythological Motif Enrichment Along the Great Circle\n"
                 "(200km corridor, 10K latitude-matched MC trials)", fontsize=12, fontweight="bold")

    # Add p-values
    for i, (e, p) in enumerate(zip(enrichments, p_values)):
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        ax.text(i, e + 0.03, f"p={p:.3f}\n{sig}", ha="center", va="bottom", fontsize=7)

    ax.legend(fontsize=9)
    ax.set_ylim(0, max(enrichments) + 0.4)

    path = os.path.join(OUTPUT_DIR, "enrichment_summary.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Written: {path}")


def fig3_cosmic_hunt_overlay(groups):
    """Overlay d'Huy's reconstructed Cosmic Hunt dispersal with the Great Circle.

    d'Huy (2012, 2016) reconstructed the Cosmic Hunt myth's dispersal from
    Africa through Eurasia and into the Americas. The approximate route
    (from his published maps) follows:
    East Africa → Middle East → Central Asia → Siberia → Beringia → Americas.

    This is a rough digitization from his published figures.
    """
    # Approximate waypoints from d'Huy's Cosmic Hunt dispersal reconstruction
    # (digitized from published maps in d'Huy 2012, 2016)
    cosmic_hunt_route = [
        # Origin in Africa
        (0, 30),      # East Africa
        (35, 30),     # Middle East / Levant
        (50, 35),     # Iran / Central Asia gateway
        (65, 42),     # Central Asia
        (80, 48),     # Kazakhstan / southern Siberia
        (100, 52),    # Central Siberia
        (120, 55),    # Eastern Siberia
        (140, 58),    # Yakutia
        (160, 62),    # Northeast Siberia
        (170, 65),    # Chukotka
        (-170, 65),   # Beringia
        (-155, 62),   # Alaska
        (-140, 58),   # Yukon
        (-120, 50),   # Pacific Northwest
        (-105, 40),   # Great Basin / Plains
        (-85, 30),    # Southeast US / Mexico
        (-70, 10),    # Central America / Caribbean
        (-60, -5),    # Amazon
        (-65, -20),   # Southern Brazil
        (-60, -35),   # Patagonia approach
    ]

    # Also add the European branch
    cosmic_hunt_europe = [
        (35, 30),     # Middle East
        (30, 40),     # Anatolia
        (25, 45),     # Balkans
        (15, 48),     # Central Europe
        (5, 48),      # Western Europe
        (-5, 55),     # Britain / Scandinavia approach
        (20, 60),     # Scandinavia
        (30, 65),     # Northern Europe
    ]

    fig, ax = plt.subplots(figsize=(16, 8))
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_facecolor("#f0f0f0")
    ax.grid(True, alpha=0.3, linestyle="--")

    # Great Circle
    draw_great_circle(ax, POLE_LAT, POLE_LON, color="red", lw=2, alpha=0.6)

    # Cosmic Hunt groups
    ch_codes = TARGET_MOTIFS["cosmic_hunt"]["codes"]
    ch_lons, ch_lats = [], []
    no_ch_lons, no_ch_lats = [], []
    for g in groups:
        has_ch = any(g["motifs"].get(c, 0) == 1 for c in ch_codes)
        if has_ch:
            ch_lons.append(g["lon"])
            ch_lats.append(g["lat"])
        else:
            no_ch_lons.append(g["lon"])
            no_ch_lats.append(g["lat"])

    ax.scatter(no_ch_lons, no_ch_lats, c="lightgray", s=10, alpha=0.4, zorder=1)
    ax.scatter(ch_lons, ch_lats, c="#1f77b4", s=30, alpha=0.7, zorder=2,
               edgecolors="black", linewidths=0.3, label=f"Cosmic Hunt attested (n={len(ch_lons)})")

    # d'Huy route
    route_lons, route_lats = zip(*cosmic_hunt_route)
    ax.plot(route_lons, route_lats, color="#ff7f0e", linewidth=2.5, alpha=0.8,
            zorder=3, label="d'Huy Cosmic Hunt dispersal (approx.)")
    ax.scatter(route_lons, route_lats, c="#ff7f0e", s=15, zorder=4)

    eu_lons, eu_lats = zip(*cosmic_hunt_europe)
    ax.plot(eu_lons, eu_lats, color="#ff7f0e", linewidth=2, alpha=0.6,
            linestyle="--", zorder=3, label="European branch")
    ax.scatter(eu_lons, eu_lats, c="#ff7f0e", s=15, zorder=4)

    ax.set_xlabel("Longitude", fontsize=10)
    ax.set_ylabel("Latitude", fontsize=10)
    ax.set_title("Cosmic Hunt Myth Distribution vs. Great Circle\n"
                 "with d'Huy's reconstructed dispersal route", fontsize=13, fontweight="bold")
    ax.legend(loc="lower left", fontsize=9)

    path = os.path.join(OUTPUT_DIR, "dispersal_route_overlay.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Written: {path}")


def fig4_bandwidth_sensitivity():
    """Plot enrichment vs bandwidth for each motif category."""
    path_in = os.path.join(OUTPUT_DIR, "bandwidth_sensitivity.json")
    if not os.path.exists(path_in):
        print("  Skipping bandwidth plot (no data)")
        return

    with open(path_in) as f:
        data = json.load(f)

    bandwidths = sorted(int(k) for k in data["bandwidths"].keys())
    cats = list(TARGET_MOTIFS.keys())

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = plt.cm.Set2(np.linspace(0, 1, len(cats)))

    for i, cat in enumerate(cats):
        enrichments = []
        for bw in bandwidths:
            r = data["bandwidths"][str(bw)].get(cat, {})
            enrichments.append(r.get("enrichment", None))
        valid = [(bw, e) for bw, e in zip(bandwidths, enrichments) if e is not None]
        if valid:
            bws, es = zip(*valid)
            ax.plot(bws, es, "-o", color=colors[i], label=cat.replace("_", " "),
                    linewidth=1.5, markersize=5)

    ax.axhline(1.0, color="black", linestyle="--", alpha=0.5, linewidth=1)
    ax.set_xlabel("Corridor Width (km)", fontsize=10)
    ax.set_ylabel("Enrichment Ratio", fontsize=10)
    ax.set_title("Motif Enrichment vs. Corridor Width", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_ylim(0, 2.0)

    path = os.path.join(OUTPUT_DIR, "bandwidth_sensitivity.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Written: {path}")


def main():
    print("Generating visualizations...")
    groups, motif_cols = load_data()

    fig1_corridor_map(groups)
    fig2_enrichment_chart()
    fig3_cosmic_hunt_overlay(groups)
    fig4_bandwidth_sensitivity()
    print("Done.")


if __name__ == "__main__":
    main()
