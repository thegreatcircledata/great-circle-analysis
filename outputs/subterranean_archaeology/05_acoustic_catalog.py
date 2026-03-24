#!/usr/bin/env python3
"""
Directive 12, Script 05: Acoustic Properties of Subterranean Circle Sites
=========================================================================
EXPLORATORY ANALYSIS — descriptive catalog only (dataset too small for
Monte Carlo testing).

Compiles all published acoustic measurements from underground sites
worldwide and computes their distance to the Great Circle.

Sources:
  - Cook et al. 2008 (Chavín de Huántar acoustics)
  - Abel et al. 2008 (Chavín resonance measurements)
  - Jahn et al. 1996 (Princeton PEAR lab — various sites)
  - Devereux 2001 "Stone Age Soundtracks"
  - Debertolis & Bisconti 2013 (Ħal Saflieni Hypogeum)
  - Reznikoff 2006 (French cave acoustics)

Output:
  - acoustic_sites_catalog.csv
  - acoustic_analysis.json
"""

import csv
import json
import math
import os

BASE = os.path.dirname(os.path.abspath(__file__))

POLE_LAT = 59.682122
POLE_LON = -138.646087
R_EARTH = 6371.0
QUARTER_CIRC = R_EARTH * math.pi / 2


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlat, dlon = la2 - la1, lo2 - lo1
    a = math.sin(dlat / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin(dlon / 2) ** 2
    return R_EARTH * 2 * math.asin(math.sqrt(min(1.0, a)))


def gc_distance(lat, lon):
    d = haversine_km(POLE_LAT, POLE_LON, lat, lon)
    return abs(d - QUARTER_CIRC)


# ── Catalog of sites with published acoustic measurements ─────────────
ACOUSTIC_SITES = [
    # Peer-reviewed measurements
    {"name": "Chavín de Huántar (galleries)",
     "lat": -9.593, "lon": -77.177,
     "resonance_hz": 100, "type": "underground gallery",
     "reference": "Cook et al. 2008; Abel et al. 2008",
     "notes": "Standing waves at ~100 Hz in underground galleries, pututu shell horns"},

    {"name": "Ħal Saflieni Hypogeum (Oracle Chamber)",
     "lat": 35.869, "lon": 14.507,
     "resonance_hz": 110, "type": "underground temple",
     "reference": "Debertolis & Bisconti 2013; Jahn et al. 1996",
     "notes": "110 Hz resonance in the Oracle Chamber, male voice range"},

    {"name": "Great Pyramid King's Chamber",
     "lat": 29.979, "lon": 31.134,
     "resonance_hz": 121, "type": "pyramid chamber",
     "reference": "Reid 2007 (cymatics); Dunn 1998 (contested)",
     "notes": "121 Hz fundamental reported, granite chamber acoustics. Note: some measurements contested"},

    {"name": "Newgrange passage tomb",
     "lat": 53.695, "lon": -6.475,
     "resonance_hz": 110, "type": "passage tomb",
     "reference": "Jahn et al. 1996 (Princeton PEAR lab)",
     "notes": "Sound amplification in the chamber, 110 Hz standing wave"},

    {"name": "Wayland's Smithy chambered tomb",
     "lat": 51.568, "lon": -1.596,
     "resonance_hz": 112, "type": "chambered tomb",
     "reference": "Jahn et al. 1996",
     "notes": "PEAR lab measurement series"},

    {"name": "Cairn L, Loughcrew",
     "lat": 53.741, "lon": -7.117,
     "resonance_hz": 110, "type": "passage tomb",
     "reference": "Jahn et al. 1996; Devereux 2001",
     "notes": "Part of systematic UK/Ireland megalithic acoustic survey"},

    {"name": "Maeshowe chambered cairn",
     "lat": 58.997, "lon": -3.188,
     "resonance_hz": 107, "type": "chambered cairn",
     "reference": "Watson & Keating 1999",
     "notes": "Orkney passage grave, acoustic measurements"},

    {"name": "Le Menec (Carnac) dolmen",
     "lat": 47.591, "lon": -3.081,
     "resonance_hz": 110, "type": "dolmen",
     "reference": "Devereux 2001",
     "notes": "Brittany megalithic site acoustic survey"},

    {"name": "Chichén Itzá (El Castillo chirp)",
     "lat": 20.683, "lon": -88.569,
     "resonance_hz": None, "type": "pyramid (exterior)",
     "reference": "Declercq et al. 2004",
     "notes": "Staircase produces chirped echo mimicking quetzal bird call. Not underground but documented acoustic design"},

    {"name": "Stonehenge (interior)",
     "lat": 51.179, "lon": -1.826,
     "resonance_hz": None, "type": "stone circle",
     "reference": "Till 2010; Fazenda et al. 2013",
     "notes": "Acoustic modeling suggests interior sound amplification. Not subterranean"},

    {"name": "Arcy-sur-Cure caves (France)",
     "lat": 47.592, "lon": 3.768,
     "resonance_hz": None, "type": "painted cave",
     "reference": "Reznikoff & Dauvois 1988; Reznikoff 2006",
     "notes": "Correlation between cave painting locations and acoustic resonance points"},

    {"name": "Le Portel cave (France)",
     "lat": 43.018, "lon": 1.606,
     "resonance_hz": None, "type": "painted cave",
     "reference": "Reznikoff & Dauvois 1988",
     "notes": "Acoustic mapping of painted vs unpainted chambers"},

    {"name": "Niaux cave (France)",
     "lat": 42.830, "lon": 1.583,
     "resonance_hz": None, "type": "painted cave",
     "reference": "Reznikoff 2006",
     "notes": "Paintings concentrate at acoustic hot spots"},

    {"name": "Font-de-Gaume (France)",
     "lat": 44.935, "lon": 1.007,
     "resonance_hz": None, "type": "painted cave",
     "reference": "Reznikoff 2006",
     "notes": "Acoustic-painting correlation documented"},

    {"name": "Diros Caves (Greece, Mani)",
     "lat": 36.639, "lon": 22.385,
     "resonance_hz": None, "type": "natural cave",
     "reference": "Informal measurement only",
     "notes": "Noted acoustic properties but no published measurements"},

    {"name": "Derinkuyu underground city",
     "lat": 38.374, "lon": 34.735,
     "resonance_hz": None, "type": "underground city",
     "reference": "No published acoustic data",
     "notes": "Ventilation shafts may have acoustic properties, no formal study"},

    {"name": "Saqqara Serapeum",
     "lat": 29.871, "lon": 31.215,
     "resonance_hz": None, "type": "underground gallery",
     "reference": "No published acoustic data",
     "notes": "Granite sarcophagi and tunnel geometry suggest acoustic potential, no study"},
]


def main():
    print("=" * 70)
    print("ACOUSTIC PROPERTIES OF SUBTERRANEAN SITES — Exploratory Catalog")
    print("=" * 70)

    # Compute distances
    for site in ACOUSTIC_SITES:
        site["gc_distance_km"] = round(gc_distance(site["lat"], site["lon"]), 1)

    # Print catalog
    print(f"\n{'Site':<40s}  {'Hz':>5s}  {'Type':<20s}  {'GC Dist':>8s}")
    print(f"{'-' * 40}  {'-' * 5}  {'-' * 20}  {'-' * 8}")
    for s in sorted(ACOUSTIC_SITES, key=lambda x: x["gc_distance_km"]):
        hz_str = f"{s['resonance_hz']:>5d}" if s["resonance_hz"] else "  N/A"
        print(f"  {s['name'][:38]:<38s}  {hz_str}  {s['type'][:18]:<20s}  {s['gc_distance_km']:>7.1f}km")

    # Sites with actual resonance measurements
    measured = [s for s in ACOUSTIC_SITES if s["resonance_hz"] is not None]
    print(f"\n--- Sites with published resonance measurements: {len(measured)} ---")

    within_200 = [s for s in measured if s["gc_distance_km"] <= 200]
    within_500 = [s for s in measured if s["gc_distance_km"] <= 500]
    print(f"  Within 200km of circle: {len(within_200)} ({[s['name'][:25] for s in within_200]})")
    print(f"  Within 500km of circle: {len(within_500)}")
    print(f"  Total with measurements: {len(measured)}")

    # Resonance frequency distribution
    freqs = [s["resonance_hz"] for s in measured]
    import numpy as np
    print(f"\n--- Resonance Frequencies ---")
    print(f"  Range: {min(freqs)}–{max(freqs)} Hz")
    print(f"  Mean: {np.mean(freqs):.1f} Hz")
    print(f"  Median: {np.median(freqs):.1f} Hz")
    print(f"  Note: 110 Hz (±10) is the dominant resonance in megalithic chambers")
    print(f"  This corresponds to the adult male chest/skull resonance range")

    # Cave painting–acoustic correlation (Reznikoff)
    reznikoff = [s for s in ACOUSTIC_SITES if "Reznikoff" in s.get("reference", "")]
    print(f"\n--- Cave Painting–Acoustic Correlation (Reznikoff) ---")
    print(f"  Sites tested: {len(reznikoff)}")
    for s in reznikoff:
        print(f"    {s['name']}: {s['gc_distance_km']:.0f}km from circle")
    print(f"  Finding: paintings cluster at locations of maximum acoustic resonance")
    print(f"  Interpretation: cave painters selected locations for acoustic properties")
    print(f"  Caveat: all Reznikoff sites are in France, ~3000+km from circle")

    # ── Save ──
    with open(os.path.join(BASE, "acoustic_sites_catalog.csv"), "w", newline="") as f:
        fields = ["name", "lat", "lon", "resonance_hz", "type",
                   "gc_distance_km", "reference", "notes"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for s in sorted(ACOUSTIC_SITES, key=lambda x: x["gc_distance_km"]):
            writer.writerow({k: s.get(k) for k in fields})

    analysis = {
        "analysis": "Acoustic properties of subterranean sites (exploratory)",
        "status": "DESCRIPTIVE ONLY — insufficient data for statistical testing",
        "total_sites_cataloged": len(ACOUSTIC_SITES),
        "sites_with_measurements": len(measured),
        "measured_within_200km": len(within_200),
        "measured_within_500km": len(within_500),
        "frequency_range_hz": {"min": min(freqs), "max": max(freqs),
                                "mean": round(float(np.mean(freqs)), 1)},
        "key_findings": [
            "Chavín de Huántar (on-circle, 100 Hz) has peer-reviewed acoustic measurements",
            "Great Pyramid (on-circle) has contested but documented measurements",
            "Saqqara Serapeum (on-circle) has no published acoustic study despite obvious potential",
            "Dominant frequency across megalithic chambers is ~110 Hz (±10)",
            "French cave painting-acoustic correlation not on circle (control data)",
        ],
        "caveat": "Dataset too small (<20 measured sites) for Monte Carlo. "
                  "Acoustic measurements are rare and often contested. "
                  "Report as catalog, not statistical test.",
    }
    with open(os.path.join(BASE, "acoustic_analysis.json"), "w") as f:
        json.dump(analysis, f, indent=2)

    print(f"\nSaved acoustic_sites_catalog.csv")
    print(f"Saved acoustic_analysis.json")


if __name__ == "__main__":
    main()
