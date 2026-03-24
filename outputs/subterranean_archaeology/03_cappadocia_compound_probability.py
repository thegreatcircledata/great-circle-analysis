#!/usr/bin/env python3
"""
Directive 12, Script 03: Cappadocia Distance Test & Compound Probability
========================================================================
1. Compute distance from all known Cappadocian underground cities to the Great Circle.
2. Compute distances from other major underground complex clusters (Petra, Beit Guvrin,
   Naqsh-e Rostam, Ajanta, Ellora) to the circle.
3. Monte Carlo compound probability: for 10,000 random great circles, how often does
   one pass within X km of MULTIPLE independent underground complex clusters?

Output:
  - cappadocia_distances.csv
  - underground_complex_compound_probability.json
  - underground_clusters_map.png
"""

import json
import math
import os
import csv
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Constants ──────────────────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
R_EARTH = 6371.0
QUARTER_CIRC = R_EARTH * math.pi / 2
N_MC = 10000

BASE = os.path.dirname(os.path.abspath(__file__))
np.random.seed(42)


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlat, dlon = la2 - la1, lo2 - lo1
    a = math.sin(dlat / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin(dlon / 2) ** 2
    return R_EARTH * 2 * math.asin(math.sqrt(min(1.0, a)))


def gc_distance_from_pole(lat, lon, plat, plon):
    d = haversine_km(plat, plon, lat, lon)
    return abs(d - QUARTER_CIRC)


def random_pole():
    z = np.random.uniform(-1, 1)
    theta = np.random.uniform(0, 2 * np.pi)
    return math.degrees(math.asin(z)), math.degrees(theta) - 180


# ── Cappadocian underground cities ────────────────────────────────────
CAPPADOCIA_CITIES = [
    {"name": "Derinkuyu", "lat": 38.374, "lon": 34.735, "depth_m": 85, "levels": 18},
    {"name": "Kaymakli", "lat": 38.463, "lon": 34.752, "depth_m": 40, "levels": 8},
    {"name": "Özkonak", "lat": 38.680, "lon": 34.730, "depth_m": 40, "levels": 10},
    {"name": "Mazi (Mataza)", "lat": 38.555, "lon": 34.665, "depth_m": 30, "levels": 4},
    {"name": "Tatlarin", "lat": 38.510, "lon": 34.558, "depth_m": 25, "levels": 5},
    {"name": "Gaziemir", "lat": 38.300, "lon": 34.700, "depth_m": 15, "levels": 4},
    {"name": "Saratli", "lat": 38.520, "lon": 34.810, "depth_m": 20, "levels": 5},
    {"name": "Acigöl", "lat": 38.558, "lon": 34.517, "depth_m": 15, "levels": 3},
    {"name": "Güzelyurt (Gelveri)", "lat": 38.305, "lon": 34.323, "depth_m": 20, "levels": 4},
    {"name": "Mucur", "lat": 39.064, "lon": 34.386, "depth_m": 10, "levels": 3},
    {"name": "Avanos", "lat": 38.716, "lon": 34.849, "depth_m": 10, "levels": 3},
    {"name": "Gülşehir", "lat": 38.758, "lon": 34.613, "depth_m": 15, "levels": 3},
    {"name": "Nevşehir (underground city)", "lat": 38.625, "lon": 34.712, "depth_m": 45, "levels": 7},
    {"name": "Çardak", "lat": 38.410, "lon": 34.600, "depth_m": 10, "levels": 3},
    {"name": "Sivasa/Gökçe", "lat": 38.485, "lon": 34.500, "depth_m": 10, "levels": 3},
    {"name": "Tilköy", "lat": 38.430, "lon": 34.780, "depth_m": 10, "levels": 2},
]

# Major underground complex clusters worldwide
UNDERGROUND_CLUSTERS = [
    {"name": "Cappadocia cluster (centroid)", "lat": 38.50, "lon": 34.68,
     "region": "Turkey", "type": "underground cities"},
    {"name": "Beit Guvrin-Maresha", "lat": 31.614, "lon": 34.896,
     "region": "Israel", "type": "cave city (~3500 chambers)"},
    {"name": "Petra", "lat": 30.329, "lon": 35.444,
     "region": "Jordan", "type": "rock-cut city (800+ facades)"},
    {"name": "Naqsh-e Rostam", "lat": 29.988, "lon": 52.874,
     "region": "Iran", "type": "rock-cut royal tombs"},
    {"name": "Ajanta Caves", "lat": 20.552, "lon": 75.700,
     "region": "India", "type": "rock-cut temples (30 caves)"},
    {"name": "Ellora Caves", "lat": 20.027, "lon": 75.179,
     "region": "India", "type": "rock-cut temples (34 caves)"},
    {"name": "Barabar Caves", "lat": 25.005, "lon": 85.063,
     "region": "India", "type": "oldest rock-cut architecture"},
    {"name": "Chavín de Huántar", "lat": -9.593, "lon": -77.177,
     "region": "Peru", "type": "underground galleries"},
    {"name": "Saqqara Serapeum", "lat": 29.871, "lon": 31.215,
     "region": "Egypt", "type": "underground galleries (~400m tunnels)"},
    {"name": "Valley of the Kings", "lat": 25.740, "lon": 32.601,
     "region": "Egypt", "type": "rock-cut tomb complex (63 tombs)"},
    {"name": "Lalibela", "lat": 12.032, "lon": 39.043,
     "region": "Ethiopia", "type": "rock-cut churches"},
    {"name": "Mogao Caves", "lat": 40.042, "lon": 94.809,
     "region": "China", "type": "rock-cut temples (492 caves)"},
]


def main():
    print("=" * 70)
    print("CAPPADOCIA DISTANCE TEST & COMPOUND PROBABILITY")
    print("=" * 70)

    # ── Part 1: Cappadocia distances ──
    print(f"\n--- Cappadocia Underground Cities ---")
    print(f"{'City':<30s}  {'Lat':>8s}  {'Lon':>8s}  {'GC Dist':>8s}  {'Depth':>6s}  {'Levels':>6s}")
    print(f"{'-' * 30}  {'-' * 8}  {'-' * 8}  {'-' * 8}  {'-' * 6}  {'-' * 6}")

    cap_results = []
    for city in CAPPADOCIA_CITIES:
        d = gc_distance_from_pole(city["lat"], city["lon"], POLE_LAT, POLE_LON)
        city["gc_distance_km"] = round(d, 1)
        cap_results.append(city)
        print(f"  {city['name']:<28s}  {city['lat']:>8.3f}  {city['lon']:>8.3f}  "
              f"{d:>7.1f}km  {city['depth_m']:>4d}m  {city['levels']:>4d}")

    dists = [c["gc_distance_km"] for c in cap_results]
    print(f"\n  Closest: {min(dists):.1f} km ({min(cap_results, key=lambda x: x['gc_distance_km'])['name']})")
    print(f"  Farthest: {max(dists):.1f} km")
    print(f"  Mean: {np.mean(dists):.1f} km")
    print(f"  Median: {np.median(dists):.1f} km")

    # Save Cappadocia CSV
    with open(os.path.join(BASE, "cappadocia_distances.csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "lat", "lon", "gc_distance_km",
                                                "depth_m", "levels"])
        writer.writeheader()
        for c in sorted(cap_results, key=lambda x: x["gc_distance_km"]):
            writer.writerow(c)
    print(f"  Saved cappadocia_distances.csv")

    # ── Part 2: All major underground clusters ──
    print(f"\n--- Major Underground Complex Clusters ---")
    print(f"{'Site':<35s}  {'Region':<10s}  {'GC Dist':>8s}")
    print(f"{'-' * 35}  {'-' * 10}  {'-' * 8}")

    for cluster in UNDERGROUND_CLUSTERS:
        d = gc_distance_from_pole(cluster["lat"], cluster["lon"], POLE_LAT, POLE_LON)
        cluster["gc_distance_km"] = round(d, 1)
        print(f"  {cluster['name']:<33s}  {cluster['region']:<10s}  {d:>7.1f}km")

    # ── Part 3: Compound probability ──
    print(f"\n--- Compound Probability Test ---")
    print(f"Testing: for random great circles, how often does one pass near")
    print(f"MULTIPLE underground complex clusters simultaneously?")

    # Use 5 geographically distinct clusters (exclude nearby duplicates)
    test_clusters = [
        {"name": "Cappadocia", "lat": 38.50, "lon": 34.68},
        {"name": "Petra/Levant", "lat": 30.95, "lon": 35.17},  # midpoint Petra/Beit Guvrin
        {"name": "Naqsh-e Rostam/Iran", "lat": 29.988, "lon": 52.874},
        {"name": "Egypt (Saqqara/VoK)", "lat": 27.8, "lon": 31.9},  # midpoint
        {"name": "Chavín/Peru", "lat": -9.593, "lon": -77.177},
    ]

    # Also test a wider set including India
    wide_clusters = test_clusters + [
        {"name": "Ajanta/India", "lat": 20.29, "lon": 75.44},  # midpoint Ajanta/Ellora
    ]

    for threshold in [200, 500, 1000]:
        print(f"\n  Threshold: {threshold} km")

        # Actual counts
        actual_near = []
        actual_near_wide = []
        for cl in test_clusters:
            d = gc_distance_from_pole(cl["lat"], cl["lon"], POLE_LAT, POLE_LON)
            if d <= threshold:
                actual_near.append(cl["name"])
        for cl in wide_clusters:
            d = gc_distance_from_pole(cl["lat"], cl["lon"], POLE_LAT, POLE_LON)
            if d <= threshold:
                actual_near_wide.append(cl["name"])

        print(f"    Actual GC hits (5 clusters): {len(actual_near)}/5 → {actual_near}")
        print(f"    Actual GC hits (6 clusters): {len(actual_near_wide)}/6 → {actual_near_wide}")

        # Monte Carlo
        mc_hits_5 = np.zeros(N_MC)
        mc_hits_6 = np.zeros(N_MC)
        for i in range(N_MC):
            plat, plon = random_pole()
            hits_5 = 0
            hits_6 = 0
            for cl in test_clusters:
                d = gc_distance_from_pole(cl["lat"], cl["lon"], plat, plon)
                if d <= threshold:
                    hits_5 += 1
            for cl in wide_clusters:
                d = gc_distance_from_pole(cl["lat"], cl["lon"], plat, plon)
                if d <= threshold:
                    hits_6 += 1
            mc_hits_5[i] = hits_5
            mc_hits_6[i] = hits_6

        n_actual_5 = len(actual_near)
        n_actual_6 = len(actual_near_wide)
        p5 = np.mean(mc_hits_5 >= n_actual_5) if n_actual_5 > 0 else 1.0
        p6 = np.mean(mc_hits_6 >= n_actual_6) if n_actual_6 > 0 else 1.0

        print(f"    MC: P(≥{n_actual_5} of 5 within {threshold}km) = {p5:.4f}")
        print(f"    MC: P(≥{n_actual_6} of 6 within {threshold}km) = {p6:.4f}")
        print(f"    MC mean hits (5 clusters): {np.mean(mc_hits_5):.2f}±{np.std(mc_hits_5):.2f}")

    # ── Part 4: Specific pairwise distances ──
    # Check if Cappadocia, Petra, and Naqsh-e Rostam are geographically independent
    print(f"\n--- Geographic Independence Check ---")
    pairs = [
        ("Cappadocia", 38.50, 34.68, "Petra", 30.329, 35.444),
        ("Cappadocia", 38.50, 34.68, "Naqsh-e Rostam", 29.988, 52.874),
        ("Petra", 30.329, 35.444, "Naqsh-e Rostam", 29.988, 52.874),
        ("Egypt", 27.8, 31.9, "Petra", 30.329, 35.444),
        ("Cappadocia", 38.50, 34.68, "Egypt", 27.8, 31.9),
    ]
    for n1, la1, lo1, n2, la2, lo2 in pairs:
        d = haversine_km(la1, lo1, la2, lo2)
        print(f"  {n1} ↔ {n2}: {d:.0f} km")

    # ── Save compound probability results ──
    compound_results = {
        "analysis": "Compound probability of proximity to multiple underground clusters",
        "n_monte_carlo": N_MC,
        "cappadocia_summary": {
            "n_cities": len(CAPPADOCIA_CITIES),
            "closest_km": round(min(dists), 1),
            "mean_km": round(float(np.mean(dists)), 1),
            "closest_city": min(cap_results, key=lambda x: x["gc_distance_km"])["name"],
        },
        "all_clusters": [
            {"name": c["name"], "region": c["region"],
             "gc_distance_km": c["gc_distance_km"]}
            for c in sorted(UNDERGROUND_CLUSTERS, key=lambda x: x["gc_distance_km"])
        ],
    }
    with open(os.path.join(BASE, "underground_complex_compound_probability.json"), "w") as f:
        json.dump(compound_results, f, indent=2)
    print(f"\nSaved underground_complex_compound_probability.json")

    # ── Visualization ──
    make_cluster_map()


def make_cluster_map():
    """Map of major underground complex clusters with the Great Circle."""
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
        ax.plot([p[1] for p in seg], [p[0] for p in seg], "r-", lw=2, alpha=0.6)

    # Plot Cappadocia cities
    cap_lats = [c["lat"] for c in CAPPADOCIA_CITIES]
    cap_lons = [c["lon"] for c in CAPPADOCIA_CITIES]
    ax.scatter(cap_lons, cap_lats, c="purple", s=40, zorder=4, marker="^",
               label=f"Cappadocia underground cities ({len(CAPPADOCIA_CITIES)})")

    # Plot other clusters
    for cl in UNDERGROUND_CLUSTERS:
        if "Cappadocia" in cl["name"]:
            continue
        color = "blue" if cl["gc_distance_km"] <= 200 else "gray"
        size = 80 if cl["gc_distance_km"] <= 200 else 50
        ax.scatter(cl["lon"], cl["lat"], c=color, s=size, zorder=5, marker="s",
                   edgecolors="black", linewidths=0.5)
        ax.annotate(f"{cl['name']}\n({cl['gc_distance_km']:.0f}km)",
                    (cl["lon"], cl["lat"]), fontsize=6,
                    xytext=(5, 5), textcoords="offset points")

    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 75)
    ax.set_title("Major Underground Complex Clusters & The Great Circle\n"
                 "(distances in km from circle)")
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(True, alpha=0.2)

    fig.tight_layout()
    path = os.path.join(BASE, "underground_clusters_map.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
