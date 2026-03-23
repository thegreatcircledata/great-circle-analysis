#!/usr/bin/env python3
"""
Directive 6: Off-Circle Cluster Membership Analysis

For each famous off-circle site, compute distance to each of the 12 high-D
great circles from the systematic search, plus the Alison circle.
Report which sites fall within 50km of which circles.
"""

import json
import numpy as np
import os

R_EARTH = 6371.0  # km
THRESHOLD_KM = 50.0

# --- Famous off-circle sites ---
FAMOUS_SITES = [
    {"name": "Göbekli Tepe",    "lat": 37.223,  "lon": 38.922},
    {"name": "Stonehenge",      "lat": 51.179,  "lon": -1.826},
    {"name": "Angkor Wat",      "lat": 13.412,  "lon": 103.867},
    {"name": "Teotihuacan",     "lat": 19.692,  "lon": -98.844},
    {"name": "Great Zimbabwe",  "lat": -20.267, "lon": 30.933},
    {"name": "Borobudur",       "lat": -7.608,  "lon": 110.204},
    {"name": "Machu Picchu",    "lat": -13.163, "lon": -72.546},
    {"name": "Chichen Itza",    "lat": 20.683,  "lon": -88.569},
    {"name": "Cahokia",         "lat": 38.655,  "lon": -90.062},
    {"name": "Karnak/Luxor",    "lat": 25.719,  "lon": 32.657},
]


def to_cart(lat_d, lon_d):
    """Convert lat/lon in degrees to unit 3D vector."""
    lat = np.radians(lat_d)
    lon = np.radians(lon_d)
    return np.array([np.cos(lat) * np.cos(lon),
                     np.cos(lat) * np.sin(lon),
                     np.sin(lat)])


def dist_to_great_circle_km(site_lat, site_lon, pole_lat, pole_lon):
    """
    Distance from a point to a great circle defined by its pole.
    The great circle is the locus of points exactly 90° from the pole.
    Distance = |90° - angular_sep(site, pole)| converted to km.
    """
    v_site = to_cart(site_lat, site_lon)
    v_pole = to_cart(pole_lat, pole_lon)
    dot = np.clip(np.dot(v_site, v_pole), -1, 1)
    ang_sep_deg = np.degrees(np.arccos(dot))
    offset_deg = abs(ang_sep_deg - 90.0)
    return offset_deg * (np.pi / 180.0) * R_EARTH


def main():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Load top 10 circles
    with open(os.path.join(base, "outputs/systematic_gc_search/top10_circles.json")) as f:
        top10 = json.load(f)

    # Load alison_rank for ranks 11-12 and the Alison pole
    with open(os.path.join(base, "outputs/systematic_gc_search/alison_rank.json")) as f:
        alison_data = json.load(f)

    # Build the 12 high-D circles (top 10 + ranks 11, 12)
    fine_peaks = alison_data["fine_peaks_summary"]
    circles = []
    for i, c in enumerate(top10):
        circles.append({
            "id": f"rank_{i+1}",
            "pole_lat": c["pole_lat"],
            "pole_lon": c["pole_lon"],
            "D": c["D"],
        })
    # Add ranks 11 and 12
    for peak in fine_peaks:
        if peak["rank"] in (11, 12):
            circles.append({
                "id": f"rank_{peak['rank']}",
                "pole_lat": peak["pole_lat"],
                "pole_lon": peak["pole_lon"],
                "D": peak["D"],
            })

    # Add Alison circle as reference
    alison_circle = {
        "id": "alison",
        "pole_lat": alison_data["alison_pole"]["lat"],
        "pole_lon": alison_data["alison_pole"]["lon"],
        "D": alison_data["alison_mc_D"],
    }

    all_circles = circles + [alison_circle]  # 13 total

    # --- Compute distances ---
    results = []
    for site in FAMOUS_SITES:
        distances = {}
        for c in all_circles:
            d = dist_to_great_circle_km(site["lat"], site["lon"],
                                        c["pole_lat"], c["pole_lon"])
            distances[c["id"]] = round(d, 1)

        closest_id = min(distances, key=distances.get)
        closest_km = distances[closest_id]
        within_50km = {cid: dkm for cid, dkm in distances.items() if dkm <= THRESHOLD_KM}

        results.append({
            "site": site["name"],
            "lat": site["lat"],
            "lon": site["lon"],
            "distances_km": distances,
            "closest_circle": closest_id,
            "closest_distance_km": closest_km,
            "within_50km": within_50km,
            "n_circles_within_50km": len(within_50km),
        })

    # --- Which circles capture multiple famous sites? ---
    circle_captures = {c["id"]: [] for c in all_circles}
    for r in results:
        for cid in r["within_50km"]:
            circle_captures[cid].append(r["site"])

    multi_capture = {cid: sites for cid, sites in circle_captures.items() if len(sites) >= 2}

    # --- Print table ---
    print("\n" + "=" * 120)
    print("DIRECTIVE 6: Off-Circle Site Membership in 12 High-D Circles + Alison Circle")
    print("=" * 120)

    # Header
    circle_ids = [c["id"] for c in all_circles]
    header = f"{'Site':<20}" + "".join(f"{cid:>10}" for cid in circle_ids) + f"{'Closest':>12}{'Dist(km)':>10}"
    print(header)
    print("-" * len(header))

    for r in results:
        row = f"{r['site']:<20}"
        for cid in circle_ids:
            d = r["distances_km"][cid]
            marker = f"{d:>9.0f}*" if d <= THRESHOLD_KM else f"{d:>10.0f}"
            row += marker
        row += f"{r['closest_circle']:>12}{r['closest_distance_km']:>10.1f}"
        print(row)

    print(f"\n* = within {THRESHOLD_KM:.0f} km threshold")

    # Summary
    print("\n--- Sites within 50 km of at least one circle ---")
    for r in results:
        if r["n_circles_within_50km"] > 0:
            circles_str = ", ".join(f"{cid} ({r['within_50km'][cid]:.1f} km)"
                                   for cid in r["within_50km"])
            print(f"  {r['site']}: {circles_str}")

    if not any(r["n_circles_within_50km"] > 0 for r in results):
        print("  (none)")

    print("\n--- Circles capturing multiple famous sites (within 50 km) ---")
    if multi_capture:
        for cid, sites in multi_capture.items():
            print(f"  {cid}: {', '.join(sites)}")
    else:
        print("  (none)")

    # Also show sites within 100 km for broader view
    print("\n--- Sites within 100 km of at least one circle ---")
    for r in results:
        close = {cid: d for cid, d in r["distances_km"].items() if d <= 100.0}
        if close:
            circles_str = ", ".join(f"{cid} ({close[cid]:.1f} km)" for cid in close)
            print(f"  {r['site']}: {circles_str}")

    # --- Save output ---
    output = {
        "description": "Directive 6: Off-circle site membership in 12 high-D circles + Alison",
        "threshold_km": THRESHOLD_KM,
        "n_circles": len(all_circles),
        "circles": [{
            "id": c["id"],
            "pole_lat": c["pole_lat"],
            "pole_lon": c["pole_lon"],
            "D": c["D"],
        } for c in all_circles],
        "site_results": results,
        "circle_captures_within_50km": circle_captures,
        "multi_capture_circles": multi_capture,
    }

    out_path = os.path.join(base, "outputs/extended_analysis/off_circle_cluster_membership.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
