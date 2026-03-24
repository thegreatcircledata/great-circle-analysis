#!/usr/bin/env python3
"""
Directive 12, Script 02: Cave-to-Construction Timeline
======================================================
Tests whether the Great Circle corridor contains sites representing
every stage of the underground architecture transition:
  Stage 1: Natural cave habitation (>40,000 BP)
  Stage 2: Cave art / ritual modification (40,000–10,000 BP)
  Stage 3: Cave burial / mortuary use (10,000–5,000 BP)
  Stage 4: Rock-cut tombs (5,000–3,000 BP)
  Stage 5: Rock-cut temples (3,000–2,000 BP)
  Stage 6: Artificial underground complexes (2,500 BP–present)

Tests:
  A. Stage presence (does corridor have all 6 stages?)
  B. Spearman correlation: complexity vs date for on-corridor sites
  C. Comparison of progression rate on-corridor vs off-corridor

Output:
  - cave_to_construction_timeline.json
  - progression_test.json
  - stage_distribution_map.png
  - underground_evolution_timeline.png
"""

import csv
import json
import math
import os
import numpy as np
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Constants ──────────────────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
R_EARTH = 6371.0
QUARTER_CIRC = R_EARTH * math.pi / 2
CORRIDOR_KM = 200  # on-corridor threshold

BASE = os.path.dirname(os.path.abspath(__file__))

# Stage definitions
STAGES = {
    1: {"name": "Natural cave habitation", "min_bp": 40000, "max_bp": None,
        "subtypes": ["natural_cave"]},
    2: {"name": "Cave art / ritual modification", "min_bp": 10000, "max_bp": 40000,
        "subtypes": ["modified_cave"]},
    3: {"name": "Cave burial / mortuary use", "min_bp": 5000, "max_bp": 10000,
        "subtypes": ["modified_cave", "natural_cave"]},
    4: {"name": "Rock-cut tombs", "min_bp": 3000, "max_bp": 5000,
        "subtypes": ["rock_cut_monumental"]},
    5: {"name": "Rock-cut temples", "min_bp": 2000, "max_bp": 3000,
        "subtypes": ["rock_cut_monumental"]},
    6: {"name": "Artificial underground complexes", "min_bp": None, "max_bp": 2500,
        "subtypes": ["artificial_complex"]},
}


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlat, dlon = la2 - la1, lo2 - lo1
    a = math.sin(dlat / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin(dlon / 2) ** 2
    return R_EARTH * 2 * math.asin(math.sqrt(min(1.0, a)))


def bearing_from_pole(lat, lon):
    """Compute bearing from pole to point (position along circle in degrees)."""
    plat, plon = math.radians(POLE_LAT), math.radians(POLE_LON)
    lat_r, lon_r = math.radians(lat), math.radians(lon)
    dlon = lon_r - plon
    x = math.sin(dlon) * math.cos(lat_r)
    y = math.cos(plat) * math.sin(lat_r) - math.sin(plat) * math.cos(lat_r) * math.cos(dlon)
    bearing = math.degrees(math.atan2(x, y)) % 360
    return bearing


def assign_stage(subtype, date_bp):
    """Assign a complexity stage based on subtype and date."""
    if date_bp is None:
        # Assign based on subtype alone
        if subtype == "natural_cave":
            return 1
        if subtype == "modified_cave":
            return 2
        if subtype == "rock_cut_monumental":
            return 5
        if subtype == "artificial_complex":
            return 6
        if subtype == "utilitarian":
            return None  # not part of the progression
        return None

    if subtype == "natural_cave" and date_bp >= 40000:
        return 1
    if subtype in ("modified_cave",) and 10000 <= date_bp < 40000:
        return 2
    if subtype in ("modified_cave", "natural_cave") and 5000 <= date_bp < 10000:
        return 3
    if subtype == "rock_cut_monumental" and 3000 <= date_bp < 5000:
        return 4
    if subtype == "rock_cut_monumental" and date_bp < 3000:
        return 5
    if subtype == "artificial_complex":
        return 6

    # Fallback to subtype-based
    if subtype == "natural_cave":
        return 1
    if subtype == "modified_cave":
        return 2
    if subtype == "rock_cut_monumental":
        return 5 if (date_bp and date_bp < 3000) else 4
    if subtype == "artificial_complex":
        return 6

    return None


def load_sites():
    path = os.path.join(BASE, "subterranean_sites_master.csv")
    sites = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            date_bp = None
            if row["date_bp"] and row["date_bp"] != "None":
                try:
                    date_bp = float(row["date_bp"])
                except ValueError:
                    pass

            sites.append({
                "name": row["name"],
                "lat": float(row["lat"]),
                "lon": float(row["lon"]),
                "subtype": row["subtype"],
                "date_bp": date_bp,
                "gc_distance_km": float(row["gc_distance_km"]),
                "cluster": row["cluster"],
            })
    return sites


def main():
    print("=" * 70)
    print("CAVE-TO-CONSTRUCTION TIMELINE")
    print("=" * 70)

    sites = load_sites()

    # Exclude utilitarian sites from the progression analysis
    sites = [s for s in sites if s["subtype"] != "utilitarian"]
    print(f"Loaded {len(sites)} non-utilitarian subterranean sites")

    # Assign stages
    for s in sites:
        s["stage"] = assign_stage(s["subtype"], s["date_bp"])
        s["bearing"] = bearing_from_pole(s["lat"], s["lon"])

    staged = [s for s in sites if s["stage"] is not None]
    on_corridor = [s for s in staged if s["gc_distance_km"] <= CORRIDOR_KM]
    off_corridor = [s for s in staged if s["gc_distance_km"] > CORRIDOR_KM]

    print(f"\nSites with assigned stage: {len(staged)}")
    print(f"  On corridor (≤{CORRIDOR_KM}km): {len(on_corridor)}")
    print(f"  Off corridor:                {len(off_corridor)}")

    # ── Test A: Stage presence ──
    print(f"\n--- Stage Presence Test ---")
    on_stages = set(s["stage"] for s in on_corridor)
    off_stages = set(s["stage"] for s in off_corridor)

    print(f"\n{'Stage':>5s}  {'Name':<35s}  {'On-corr':>8s}  {'Off-corr':>8s}")
    print(f"{'-' * 5}  {'-' * 35}  {'-' * 8}  {'-' * 8}")

    timeline_data = []
    for stage_num in range(1, 7):
        stage_info = STAGES[stage_num]
        on_count = sum(1 for s in on_corridor if s["stage"] == stage_num)
        off_count = sum(1 for s in off_corridor if s["stage"] == stage_num)
        present_on = "✓" if stage_num in on_stages else "✗"
        present_off = "✓" if stage_num in off_stages else "✗"
        print(f"  {stage_num:>3d}  {stage_info['name']:<35s}  {present_on} {on_count:>5d}  "
              f"{present_off} {off_count:>5d}")

        # Collect on-corridor sites for this stage
        for s in on_corridor:
            if s["stage"] == stage_num:
                timeline_data.append({
                    "name": s["name"],
                    "stage": stage_num,
                    "stage_name": stage_info["name"],
                    "date_bp": s["date_bp"],
                    "lat": s["lat"],
                    "lon": s["lon"],
                    "bearing": round(s["bearing"], 1),
                    "gc_distance_km": s["gc_distance_km"],
                    "cluster": s["cluster"],
                })

    n_on = len(on_stages)
    n_off = len(off_stages)
    print(f"\n  Stages represented on corridor:  {n_on}/6 → {on_stages}")
    print(f"  Stages represented off corridor: {n_off}/6 → {off_stages}")

    # ── Test B: Spearman correlation (complexity vs date) ──
    print(f"\n--- Temporal Progression Test ---")

    # Use only sites with dates
    on_dated = [s for s in on_corridor if s["date_bp"] is not None]
    off_dated = [s for s in off_corridor if s["date_bp"] is not None]

    def spearman_test(subset, label):
        if len(subset) < 5:
            print(f"  {label}: too few dated sites ({len(subset)})")
            return None
        stages = [s["stage"] for s in subset]
        dates = [s["date_bp"] for s in subset]
        # Higher stage = more complex, higher date_bp = older
        # We expect negative correlation (older → simpler stage)
        rho, p = stats.spearmanr(dates, stages)
        print(f"  {label}: N={len(subset)}, rho={rho:+.3f}, p={p:.4f} "
              f"{'→ Significant' if p < 0.05 else ''}")
        return {"n": len(subset), "rho": round(float(rho), 3),
                "p_value": round(float(p), 4)}

    on_spearman = spearman_test(on_dated, "On-corridor")
    off_spearman = spearman_test(off_dated, "Off-corridor")

    # ── Test C: Stage distribution by region ──
    print(f"\n--- Geographic Distribution of Stages ---")
    regions = {}
    for s in on_corridor:
        r = s["cluster"]
        if r not in regions:
            regions[r] = set()
        regions[r].add(s["stage"])

    for r, stages in sorted(regions.items()):
        print(f"  {r:<20s}: stages {sorted(stages)}")

    # ── Save results ──
    timeline_output = {
        "analysis": "Cave-to-construction timeline along the corridor",
        "corridor_threshold_km": CORRIDOR_KM,
        "stage_definitions": {str(k): v["name"] for k, v in STAGES.items()},
        "stage_presence": {
            "on_corridor": sorted(list(on_stages)),
            "off_corridor": sorted(list(off_stages)),
            "on_corridor_count": n_on,
            "off_corridor_count": n_off,
        },
        "on_corridor_sites": timeline_data,
    }
    with open(os.path.join(BASE, "cave_to_construction_timeline.json"), "w") as f:
        json.dump(timeline_output, f, indent=2)

    progression_output = {
        "analysis": "Temporal progression test (complexity vs age)",
        "on_corridor_spearman": on_spearman,
        "off_corridor_spearman": off_spearman,
        "geographic_distribution": {r: sorted(list(s)) for r, s in regions.items()},
    }
    with open(os.path.join(BASE, "progression_test.json"), "w") as f:
        json.dump(progression_output, f, indent=2)

    print(f"\nSaved cave_to_construction_timeline.json")
    print(f"Saved progression_test.json")

    # ── Visualizations ──
    make_timeline_figure(staged, on_corridor, off_corridor)
    make_stage_map(on_corridor)


def make_timeline_figure(all_staged, on_corr, off_corr):
    """Temporal progression: stage vs date for on/off corridor."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    for ax, subset, title in [(ax1, on_corr, f"On-corridor (≤{CORRIDOR_KM}km)"),
                               (ax2, off_corr, "Off-corridor")]:
        dated = [s for s in subset if s["date_bp"] is not None]
        if not dated:
            ax.set_title(f"{title}\n(no dated sites)")
            continue

        dates = [s["date_bp"] for s in dated]
        stages = [s["stage"] for s in dated]

        colors = {1: "#4CAF50", 2: "#2196F3", 3: "#9C27B0",
                  4: "#FF9800", 5: "#F44336", 6: "#795548"}
        c = [colors.get(s, "#999") for s in stages]

        ax.scatter(dates, stages, c=c, s=30, alpha=0.6, zorder=3)
        ax.set_xlabel("Age (years BP)")
        ax.set_ylabel("Complexity Stage")
        ax.set_title(f"{title}\n(N={len(dated)} dated sites)")
        ax.set_yticks(range(1, 7))
        ax.set_yticklabels([STAGES[i]["name"][:30] for i in range(1, 7)], fontsize=7)
        ax.invert_xaxis()
        ax.grid(True, alpha=0.2)

    fig.suptitle("Underground Architecture: Temporal Progression", fontsize=14)
    fig.tight_layout()
    path = os.path.join(BASE, "underground_evolution_timeline.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved {path}")


def make_stage_map(on_corridor):
    """Map of on-corridor sites colored by stage."""
    fig, ax = plt.subplots(figsize=(18, 9))

    # Draw circle
    plat, plon = math.radians(POLE_LAT), math.radians(POLE_LON)
    pts = []
    for b in range(720):
        bearing = math.radians(b * 0.5)
        lat = math.asin(
            math.sin(plat) * math.cos(math.pi / 2)
            + math.cos(plat) * math.sin(math.pi / 2) * math.cos(bearing))
        lon = plon + math.atan2(
            math.sin(bearing) * math.sin(math.pi / 2) * math.cos(plat),
            math.cos(math.pi / 2) - math.sin(plat) * math.sin(lat))
        pts.append((math.degrees(lat), math.degrees(lon)))

    segs = []
    cur = [pts[0]]
    for i in range(1, len(pts)):
        if abs(pts[i][1] - pts[i - 1][1]) > 180:
            segs.append(cur)
            cur = [pts[i]]
        else:
            cur.append(pts[i])
    segs.append(cur)
    for seg in segs:
        ax.plot([p[1] for p in seg], [p[0] for p in seg], "r-", lw=1.5, alpha=0.5)

    colors = {1: "#4CAF50", 2: "#2196F3", 3: "#9C27B0",
              4: "#FF9800", 5: "#F44336", 6: "#795548"}

    for stage_num in range(1, 7):
        subset = [s for s in on_corridor if s["stage"] == stage_num]
        if not subset:
            continue
        lats = [s["lat"] for s in subset]
        lons = [s["lon"] for s in subset]
        ax.scatter(lons, lats, c=colors[stage_num], s=25, alpha=0.7, zorder=3,
                   label=f"Stage {stage_num}: {STAGES[stage_num]['name']}")

    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 75)
    ax.set_title(f"On-Corridor Subterranean Sites by Stage (≤{CORRIDOR_KM}km)")
    ax.legend(loc="lower left", fontsize=7)
    ax.grid(True, alpha=0.2)

    fig.tight_layout()
    path = os.path.join(BASE, "stage_distribution_map.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
