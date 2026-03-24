#!/usr/bin/env python3
"""
Analysis 1: The Great Circle as a Frozen Equator
=================================================
Tests whether the Great Circle could have been Earth's equator at any past epoch.
Examines axial precession, true polar wander, and obliquity variation.
"""

import json
import math
import os

# Great Circle pole
POLE_LAT = 59.682122  # °N
POLE_LON = -138.646087  # °W

# Angular distance from GC pole to geographic North Pole
ang_dist = 90.0 - POLE_LAT  # = 30.318°

results = {
    "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
    "angular_distance_gc_pole_to_north_pole_deg": round(ang_dist, 3),
}

# ============================================================
# 1. Axial Precession
# ============================================================
# Axial precession changes the direction of Earth's spin axis relative to the
# STARS (the celestial pole traces a circle around the ecliptic pole with period
# ~25,772 years). But it does NOT change the geographic pole's position on Earth's
# surface. The North Pole stays at 90°N geographic regardless of precession.
#
# Therefore, precession alone CANNOT make the Great Circle into the equator.
# The equator is always 90° from the geographic pole, and the GC pole is at 59.68°N,
# which is 30.3° from the geographic pole — precession doesn't change this.

results["axial_precession"] = {
    "conclusion": "CANNOT align Great Circle with equator",
    "reason": "Axial precession changes the celestial pole direction but NOT the geographic pole position on Earth's surface. The geographic North Pole remains at 90°N regardless of precession epoch. The Great Circle pole is 30.3° from the geographic pole at ALL epochs.",
    "period_years": 25772,
}

# ============================================================
# 2. True Polar Wander (TPW)
# ============================================================
# TPW shifts the geographic pole relative to the mantle. Typical rate: ~1°/Myr
# (Steinberger & Torsvik 2008). Some estimates are slightly faster.
#
# To move the pole 30.3° would require ~30 million years at typical rates.
# This is orders of magnitude beyond human timescales.
#
# Rapid TPW events:
# - Kirschvink 1997 proposed "inertial interchange true polar wander" (IITPW)
#   in the Late Precambrian (~800 Ma), potentially ~90° in ~10-15 Myr
# - No convincing evidence for rapid TPW in the Quaternary (<2.6 Ma)
# - The most cited Quaternary TPW rates are 0.5-1.0°/Myr (Besse & Courtillot 2002)
# - Even the fastest proposed Quaternary excursions are <5° over millions of years
# - No mechanism exists for >1° TPW within human history (~300,000 years)

tpw_time_myr = ang_dist / 1.0  # at 1°/Myr
tpw_time_fast_myr = ang_dist / 5.0  # at hypothetical fast rate of 5°/Myr

results["true_polar_wander"] = {
    "conclusion": "Great Circle was NOT Earth's equator at any point in human history",
    "degrees_needed": round(ang_dist, 1),
    "typical_rate_deg_per_myr": 1.0,
    "time_at_typical_rate_myr": round(tpw_time_myr, 1),
    "fastest_proposed_quaternary_rate_deg_per_myr": 5.0,
    "time_at_fastest_rate_myr": round(tpw_time_fast_myr, 1),
    "human_timescale_note": "Even at the fastest proposed Quaternary TPW rate (5°/Myr), moving the pole 30.3° would take ~6 million years. Anatomically modern humans have existed for ~300,000 years. The pole has moved <0.3° in that time.",
    "references": [
        "Steinberger & Torsvik 2008 (TPW rates)",
        "Kirschvink 1997 (IITPW in Precambrian)",
        "Besse & Courtillot 2002 (Quaternary polar wander)",
    ],
}

# ============================================================
# 3. Obliquity Variation
# ============================================================
# Earth's axial tilt oscillates between 22.1° and 24.5° with ~41,000-year period.
# Current: ~23.44°. At 10,000 BP: ~24.14°
# This changes the tropics and arctic circles, NOT the equator itself.
# But let's check if the GC inclination matches any obliquity-related angle.

current_obliquity = 23.44
min_obliquity = 22.1
max_obliquity = 24.5
gc_inclination = ang_dist  # 30.318° — angle of GC relative to equator

# Check various obliquity-related angles
obliquity_angles = {
    "current_obliquity": current_obliquity,
    "max_obliquity": max_obliquity,
    "min_obliquity": min_obliquity,
    "2x_current_obliquity": 2 * current_obliquity,
    "2x_max_obliquity": 2 * max_obliquity,
    "obliquity_plus_complement": current_obliquity + (90 - current_obliquity),  # = 90
    "co-obliquity (90 - obliquity)": 90 - current_obliquity,
}

closest_match = None
closest_diff = 999
for name, angle in obliquity_angles.items():
    diff = abs(gc_inclination - angle)
    if diff < closest_diff:
        closest_diff = diff
        closest_match = name

results["obliquity_variation"] = {
    "gc_inclination_deg": round(gc_inclination, 3),
    "obliquity_related_angles": {k: round(v, 2) for k, v in obliquity_angles.items()},
    "closest_match": closest_match,
    "closest_match_offset_deg": round(closest_diff, 2),
    "conclusion": f"The Great Circle's inclination (30.3°) does not closely match any simple obliquity-related angle. Closest: {closest_match} (off by {closest_diff:.1f}°). No astronomical significance via obliquity.",
    "note": "Obliquity changes the positions of the tropics and arctic circles, not the equator. Even if the inclination matched an obliquity angle, there is no physical mechanism linking a surface great circle to obliquity.",
}

# ============================================================
# Summary
# ============================================================
results["overall_conclusion"] = (
    "The Great Circle was NOT Earth's equator at any point in human history. "
    "Axial precession does not move the geographic pole. True polar wander operates "
    "at ~1°/Myr, requiring ~30 million years to move the pole the necessary 30.3°. "
    "The circle's inclination does not match any obliquity-related angle. "
    "The 'ancient equator' claim is definitively refuted."
)

out_path = os.path.join(os.path.dirname(__file__), "frozen_equator_analysis.json")
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)
print("Frozen equator analysis complete.")
print(f"  GC pole to North Pole: {ang_dist:.3f}°")
print(f"  TPW time needed: ~{tpw_time_myr:.0f} Myr")
print(f"  Closest obliquity match: {closest_match} (off by {closest_diff:.1f}°)")
print(f"  Conclusion: DEFINITIVELY NULL")
print(f"  Output: {out_path}")
