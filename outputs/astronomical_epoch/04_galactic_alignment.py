#!/usr/bin/env python3
"""
Analysis 4: Galactic Plane Intersection
========================================
Two great circles on a sphere always intersect. This analysis computes the
mutual inclination between the Great Circle (projected onto the sky) and
the Milky Way (galactic plane) at each epoch, and tests significance via
Monte Carlo.
"""

import json
import csv
import math
import os
import random
import numpy as np
from astropy.coordinates import SkyCoord, FK5
from astropy.time import Time
import astropy.units as u

random.seed(42)
np.random.seed(42)

OUT_DIR = os.path.dirname(__file__)

# Great Circle pole geographic coordinates
GC_POLE_LAT = 59.682122
GC_POLE_LON = -138.646087

# North Galactic Pole (J2000 ICRS)
NGP_RA = 192.859  # degrees
NGP_DEC = 27.128  # degrees


def angular_sep(ra1, dec1, ra2, dec2):
    """Angular separation in degrees."""
    ra1, dec1, ra2, dec2 = map(math.radians, [ra1, dec1, ra2, dec2])
    cos_d = (math.sin(dec1) * math.sin(dec2) +
             math.cos(dec1) * math.cos(dec2) * math.cos(ra2 - ra1))
    return math.degrees(math.acos(max(-1, min(1, cos_d))))


def gc_pole_celestial(epoch_year):
    """Get the GC pole's celestial coordinates at epoch (from Analysis 3 method)."""
    t = Time(epoch_year, format='jyear')
    gc_in_fk5 = SkyCoord(ra=(360 + GC_POLE_LON) % 360 * u.deg,
                          dec=GC_POLE_LAT * u.deg,
                          frame=FK5(equinox=t))
    gc_icrs = gc_in_fk5.icrs
    return gc_icrs.ra.deg, gc_icrs.dec.deg


def mutual_inclination(pole1_ra, pole1_dec, pole2_ra, pole2_dec):
    """Compute mutual inclination between two great circles given their poles.

    Two great circles have mutual inclination = angular separation between their poles.
    If the poles are close, the great circles are nearly aligned.
    If the poles are 90° apart, the great circles are perpendicular.
    """
    sep = angular_sep(pole1_ra, pole1_dec, pole2_ra, pole2_dec)
    # Mutual inclination is the angle between the two planes.
    # If poles are separated by angle α, the planes are inclined by α.
    # But if α > 90°, the inclination is 180° - α (the supplementary angle).
    if sep > 90:
        sep = 180 - sep
    return sep


# ============================================================
# Compute mutual inclination across epochs
# ============================================================
print("=" * 70)
print("ANALYSIS 4: GALACTIC PLANE INTERSECTION")
print("=" * 70)

epochs = list(range(-23000, 3001, 100))
alignment_data = []

print(f"Computing mutual inclination with galactic plane over {len(epochs)} epochs...")

for epoch in epochs:
    try:
        gc_ra, gc_dec = gc_pole_celestial(epoch)
        incl = mutual_inclination(gc_ra, gc_dec, NGP_RA, NGP_DEC)
        alignment_data.append({
            "epoch": epoch,
            "gc_pole_ra": round(gc_ra, 3),
            "gc_pole_dec": round(gc_dec, 3),
            "mutual_inclination_deg": round(incl, 2),
        })
    except Exception:
        pass

# Find minimum inclination (maximum alignment)
min_entry = min(alignment_data, key=lambda x: x["mutual_inclination_deg"])
max_entry = max(alignment_data, key=lambda x: x["mutual_inclination_deg"])

print(f"\nMinimum mutual inclination: {min_entry['mutual_inclination_deg']:.2f}° at epoch {min_entry['epoch']}")
print(f"Maximum mutual inclination: {max_entry['mutual_inclination_deg']:.2f}° at epoch {max_entry['epoch']}")

# Current epoch
current = next(d for d in alignment_data if d["epoch"] == 2000)
print(f"Current (2000 CE): {current['mutual_inclination_deg']:.2f}°")

# ============================================================
# Monte Carlo: 10,000 random great circles on Earth
# ============================================================
print(f"\nMonte Carlo: 10,000 random great circles...")

N_MC = 10000
mc_min_inclinations = []

for i in range(N_MC):
    # Random pole on sphere
    z = random.uniform(-1, 1)
    rand_lat = math.degrees(math.asin(z))
    rand_lon = random.uniform(-180, 180)

    # Compute min inclination across all epochs for this random circle
    min_incl = 999
    # Sample fewer epochs for MC (every 500 years instead of 100)
    for epoch in range(-23000, 3001, 500):
        try:
            t = Time(epoch, format='jyear')
            rc_in_fk5 = SkyCoord(ra=(360 + rand_lon) % 360 * u.deg,
                                  dec=rand_lat * u.deg,
                                  frame=FK5(equinox=t))
            rc_icrs = rc_in_fk5.icrs
            incl = mutual_inclination(rc_icrs.ra.deg, rc_icrs.dec.deg, NGP_RA, NGP_DEC)
            if incl < min_incl:
                min_incl = incl
        except Exception:
            pass

    mc_min_inclinations.append(min_incl)

    if (i + 1) % 1000 == 0:
        print(f"  {i+1}/{N_MC} done...")

gc_min_incl = min_entry["mutual_inclination_deg"]
mc_min_inclinations.sort()
count_better = sum(1 for x in mc_min_inclinations if x <= gc_min_incl)
percentile = count_better / N_MC * 100

mc_mean = np.mean(mc_min_inclinations)
mc_std = np.std(mc_min_inclinations)
mc_median = np.median(mc_min_inclinations)

print(f"\nGC minimum inclination: {gc_min_incl:.2f}°")
print(f"MC mean minimum: {mc_mean:.2f}° ± {mc_std:.2f}°")
print(f"MC median: {mc_median:.2f}°")
print(f"GC percentile: {percentile:.1f}%")

# ============================================================
# Save results
# ============================================================
results = {
    "galactic_north_pole_j2000": {"ra": NGP_RA, "dec": NGP_DEC},
    "gc_pole_geographic": {"lat": GC_POLE_LAT, "lon": GC_POLE_LON},
    "alignment_sweep": {
        "n_epochs": len(alignment_data),
        "epoch_range": [alignment_data[0]["epoch"], alignment_data[-1]["epoch"]],
        "step_years": 100,
    },
    "minimum_inclination": {
        "epoch": min_entry["epoch"],
        "inclination_deg": min_entry["mutual_inclination_deg"],
        "gc_pole_ra": min_entry["gc_pole_ra"],
        "gc_pole_dec": min_entry["gc_pole_dec"],
    },
    "maximum_inclination": {
        "epoch": max_entry["epoch"],
        "inclination_deg": max_entry["mutual_inclination_deg"],
    },
    "current_epoch": {
        "inclination_deg": current["mutual_inclination_deg"],
    },
    "monte_carlo": {
        "n_trials": N_MC,
        "gc_min_inclination_deg": round(gc_min_incl, 2),
        "mc_mean_min_inclination_deg": round(mc_mean, 2),
        "mc_std_deg": round(mc_std, 2),
        "mc_median_deg": round(mc_median, 2),
        "gc_percentile": round(percentile, 1),
    },
    "conclusion": (
        f"The Great Circle's minimum mutual inclination with the galactic plane is "
        f"{gc_min_incl:.1f}° (at epoch {min_entry['epoch']}). "
        f"Monte Carlo testing ({N_MC} random circles) shows a mean minimum inclination of "
        f"{mc_mean:.1f}° ± {mc_std:.1f}°. The Great Circle is at the {percentile:.1f}th percentile. "
        f"{'This is within normal range — not significantly aligned with the galactic plane.' if percentile > 5 else 'This represents an unusually close alignment with the galactic plane.'}"
    ),
}

# Save alignment time series as CSV
csv_path = os.path.join(OUT_DIR, "galactic_alignment_timeseries.csv")
with open(csv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["epoch", "gc_pole_ra", "gc_pole_dec", "mutual_inclination_deg"])
    writer.writeheader()
    writer.writerows(alignment_data)

json_path = os.path.join(OUT_DIR, "galactic_alignment.json")
with open(json_path, "w") as f:
    json.dump(results, f, indent=2)

print(f"\nOutputs:")
print(f"  {csv_path}")
print(f"  {json_path}")
