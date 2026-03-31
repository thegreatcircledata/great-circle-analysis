#!/usr/bin/env python3
"""
Task 3 — Corridor width profile: monument-to-settlement ratio at 5 km
intervals from the Great Circle (pole: 59.682122N, 138.646087W).

Pleiades sites filtered to minDate <= -2000 (pre-2000 BCE).
Dual-classified sites (both monumental + settlement types) counted in both.

Produces:
  A) Global profile (all qualifying Pleiades sites)
  B) Egypt-Iran corridor profile (20-38 N, 25-60 E)

Output: corridor_width_profile.json
"""

import csv
import json
import math

# ---------- constants ----------
POLE_LAT = math.radians(59.682122)
POLE_LON = math.radians(-138.646087)
R_EARTH = 6371.0  # km

MONUMENTAL = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2", "amphitheatre",
    "hippodrome", "aqueduct", "dam", "lighthouse", "arch", "theatre",
    "basilica", "acropolis", "nuraghe", "tumulus", "shrine",
}
SETTLEMENT = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry", "bath",
    "cistern", "bridge", "fortified-settlement", "townhouse", "production",
}

BANDS = [(i, i + 5) for i in range(0, 100, 5)]
CORRIDOR = {"lat": (20, 38), "lon": (25, 60)}

CSV_PATH = "/Users/elliotallan/megalith_site_research/data/pleiades/pleiades-places-latest.csv"
OUT_PATH = "/Users/elliotallan/megalith_site_research/outputs/merged_paper_prep/corridor_width_profile.json"


# ---------- geometry ----------
def distance_to_great_circle_km(lat_rad, lon_rad):
    """Shortest distance (km) from a point to the great circle whose pole
    is (POLE_LAT, POLE_LON).  GC = locus at 90 deg from pole."""
    cos_ang = (math.sin(lat_rad) * math.sin(POLE_LAT)
               + math.cos(lat_rad) * math.cos(POLE_LAT)
               * math.cos(lon_rad - POLE_LON))
    cos_ang = max(-1.0, min(1.0, cos_ang))
    return abs(math.acos(cos_ang) - math.pi / 2) * R_EARTH


# ---------- data loading ----------
def load_sites():
    sites = []
    with open(CSV_PATH, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            # date filter
            md = row.get("minDate", "").strip()
            if not md:
                continue
            try:
                min_date = float(md)
            except ValueError:
                continue
            if min_date > -2000:
                continue

            # coordinates
            rlat, rlon = row.get("reprLat", "").strip(), row.get("reprLong", "").strip()
            if not rlat or not rlon:
                continue
            try:
                lat, lon = float(rlat), float(rlon)
            except ValueError:
                continue

            # feature classification
            types = {t.strip() for t in row.get("featureTypes", "").split(",") if t.strip()}
            is_mon = bool(types & MONUMENTAL)
            is_set = bool(types & SETTLEMENT)
            if not is_mon and not is_set:
                continue

            dist = distance_to_great_circle_km(math.radians(lat), math.radians(lon))
            sites.append(dict(lat=lat, lon=lon, is_mon=is_mon, is_set=is_set, dist_km=dist))
    return sites


# ---------- profile computation ----------
def compute_profile(site_list):
    profile = []
    for lo, hi in BANDS:
        band = [s for s in site_list if lo <= s["dist_km"] < hi]
        mon = sum(1 for s in band if s["is_mon"])
        sett = sum(1 for s in band if s["is_set"])
        ratio = round(mon / sett, 4) if sett > 0 else None
        profile.append(dict(band_km=f"{lo}-{hi}", monuments=mon, settlements=sett, ratio=ratio))
    return profile


def in_corridor(s):
    return (CORRIDOR["lat"][0] <= s["lat"] <= CORRIDOR["lat"][1]
            and CORRIDOR["lon"][0] <= s["lon"] <= CORRIDOR["lon"][1])


# ---------- main ----------
def main():
    sites = load_sites()
    corridor_sites = [s for s in sites if in_corridor(s)]

    global_profile = compute_profile(sites)
    corridor_profile = compute_profile(corridor_sites)

    # Print summary
    print(f"Pre-2000 BCE classified sites: {len(sites)} total, {len(corridor_sites)} in corridor")
    for label, prof in [("Global", global_profile), ("Corridor", corridor_profile)]:
        print(f"\n{label} profile:")
        for b in prof:
            print(f"  {b['band_km']:>7s} km  mon={b['monuments']:4d}  sett={b['settlements']:4d}  ratio={b['ratio']}")

    # Check for 1.43 at ~10 km
    print("\n--- 1.43 ratio check ---")
    finding = "Neither profile produces a 1.43 ratio at ~10 km."
    for label, prof in [("global", global_profile), ("egypt_iran_corridor", corridor_profile)]:
        for b in prof:
            if b["ratio"] is not None and abs(b["ratio"] - 1.43) < 0.05:
                finding = f"{label} produces ratio {b['ratio']} at {b['band_km']} km."
                print(f"  Match: {finding}")

    # The 80-85 km global band at 1.3333 is the closest to 1.43 in any band.
    # With pre-2000 BCE Pleiades data, the near-circle bands are heavily
    # monument-dominated (Giza pyramids cluster at 5-10 km) so the ratio is
    # much higher than 1.43 in those bands.
    print(f"  {finding}")

    # Save
    output = {
        "parameters": {
            "pole_lat_deg": 59.682122,
            "pole_lon_deg": -138.646087,
            "min_date_filter": -2000,
            "band_width_km": 5,
            "max_distance_km": 100,
            "monumental_types": sorted(MONUMENTAL),
            "settlement_types": sorted(SETTLEMENT),
            "classification": "dual-classified sites counted in both categories",
        },
        "global_profile": global_profile,
        "egypt_iran_corridor_profile": corridor_profile,
        "corridor_bounds": {"lat_range": list(CORRIDOR["lat"]),
                            "lon_range": list(CORRIDOR["lon"])},
        "ratio_143_finding": finding,
    }
    with open(OUT_PATH, "w") as fh:
        json.dump(output, fh, indent=2)
    print(f"\nSaved to {OUT_PATH}")


if __name__ == "__main__":
    main()
