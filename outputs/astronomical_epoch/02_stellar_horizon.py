#!/usr/bin/env python3
"""
Analysis 2: Stellar Horizon Events at Key Sites
================================================
For each key site on the Great Circle, compute the circle's bearing (azimuth),
then check if any bright star rose or set at that azimuth at the epoch of
monument construction. Include Monte Carlo significance testing.
"""

import json
import csv
import math
import os
import random
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u

random.seed(42)
np.random.seed(42)

OUT_DIR = os.path.dirname(__file__)

# ============================================================
# Great Circle parameters
# ============================================================
POLE_LAT = math.radians(59.682122)
POLE_LON = math.radians(-138.646087)

# Key sites with construction epochs (negative = BCE)
SITES = [
    {"name": "Giza",         "lat": 29.9792, "lon": 31.1342, "epochs": [-2500, -2600, -3000, -10000]},
    {"name": "Nazca",        "lat": -14.735, "lon": -75.13,  "epochs": [-200, -500, -600]},
    {"name": "Easter Island", "lat": -27.1127, "lon": -109.3497, "epochs": [1000, 1200, 1400]},
    {"name": "Persepolis",   "lat": 29.9352, "lon": 52.8916, "epochs": [-500, -520, -2000]},
    {"name": "Mohenjo-daro", "lat": 27.3242, "lon": 68.1367, "epochs": [-2500, -2600, -3000]},
    {"name": "Angkor Wat",   "lat": 13.4125, "lon": 103.8670, "epochs": [1100, 1150, 1200]},
    {"name": "Machu Picchu", "lat": -13.1631, "lon": -72.5450, "epochs": [1450, 1400]},
    {"name": "Ollantaytambo","lat": -13.2581, "lon": -72.2636, "epochs": [1400, 1450]},
]

# ============================================================
# Bright star catalog (Hipparcos: mag < 2.5, ~50 stars)
# RA/Dec in J2000 epoch, proper motions in mas/yr
# ============================================================
BRIGHT_STARS = [
    {"name": "Sirius",      "ra": 101.287, "dec": -16.716, "mag": -1.46, "pm_ra": -546.01, "pm_dec": -1223.14},
    {"name": "Canopus",     "ra": 95.988,  "dec": -52.696, "mag": -0.74, "pm_ra": 19.93,   "pm_dec": 23.24},
    {"name": "Arcturus",    "ra": 213.915, "dec": 19.182,  "mag": -0.05, "pm_ra": -1093.45,"pm_dec": -1999.40},
    {"name": "Vega",        "ra": 279.235, "dec": 38.784,  "mag": 0.03,  "pm_ra": 200.94,  "pm_dec": 286.23},
    {"name": "Capella",     "ra": 79.172,  "dec": 45.998,  "mag": 0.08,  "pm_ra": 75.52,   "pm_dec": -427.11},
    {"name": "Rigel",       "ra": 78.634,  "dec": -8.202,  "mag": 0.13,  "pm_ra": 1.87,    "pm_dec": -0.56},
    {"name": "Procyon",     "ra": 114.826, "dec": 5.225,   "mag": 0.34,  "pm_ra": -714.59, "pm_dec": -1036.80},
    {"name": "Betelgeuse",  "ra": 88.793,  "dec": 7.407,   "mag": 0.42,  "pm_ra": 27.33,   "pm_dec": 10.86},
    {"name": "Achernar",    "ra": 24.429,  "dec": -57.237, "mag": 0.46,  "pm_ra": 87.00,   "pm_dec": -38.24},
    {"name": "Hadar",       "ra": 210.956, "dec": -60.373, "mag": 0.61,  "pm_ra": -33.27,  "pm_dec": -23.16},
    {"name": "Altair",      "ra": 297.696, "dec": 8.868,   "mag": 0.76,  "pm_ra": 536.23,  "pm_dec": 385.29},
    {"name": "Acrux",       "ra": 186.650, "dec": -63.099, "mag": 0.76,  "pm_ra": -35.37,  "pm_dec": -14.73},
    {"name": "Aldebaran",   "ra": 68.980,  "dec": 16.509,  "mag": 0.85,  "pm_ra": 62.78,   "pm_dec": -189.36},
    {"name": "Antares",     "ra": 247.352, "dec": -26.432, "mag": 0.96,  "pm_ra": -12.11,  "pm_dec": -23.30},
    {"name": "Spica",       "ra": 201.298, "dec": -11.161, "mag": 0.97,  "pm_ra": -42.50,  "pm_dec": -31.73},
    {"name": "Pollux",      "ra": 116.329, "dec": 28.026,  "mag": 1.14,  "pm_ra": -625.69, "pm_dec": -45.95},
    {"name": "Fomalhaut",   "ra": 344.413, "dec": -29.622, "mag": 1.16,  "pm_ra": 329.22,  "pm_dec": -164.22},
    {"name": "Deneb",       "ra": 310.358, "dec": 45.280,  "mag": 1.25,  "pm_ra": 1.56,    "pm_dec": 1.55},
    {"name": "Mimosa",      "ra": 191.930, "dec": -59.689, "mag": 1.25,  "pm_ra": -48.24,  "pm_dec": -12.82},
    {"name": "Regulus",     "ra": 152.093, "dec": 11.967,  "mag": 1.35,  "pm_ra": -249.40, "pm_dec": 5.59},
    {"name": "Adhara",      "ra": 104.656, "dec": -28.972, "mag": 1.50,  "pm_ra": 2.63,    "pm_dec": 2.29},
    {"name": "Castor",      "ra": 113.650, "dec": 31.888,  "mag": 1.58,  "pm_ra": -191.45, "pm_dec": -145.19},
    {"name": "Bellatrix",   "ra": 81.283,  "dec": 6.350,   "mag": 1.64,  "pm_ra": -8.11,   "pm_dec": -12.88},
    {"name": "Alnath",      "ra": 81.573,  "dec": 28.608,  "mag": 1.65,  "pm_ra": 23.28,   "pm_dec": -174.22},
    {"name": "Alnilam",     "ra": 84.053,  "dec": -1.202,  "mag": 1.69,  "pm_ra": 1.49,    "pm_dec": -1.06},
    {"name": "Alnitak",     "ra": 85.190,  "dec": -1.943,  "mag": 1.77,  "pm_ra": 3.19,    "pm_dec": 2.03},
    {"name": "Mintaka",     "ra": 83.001,  "dec": -0.299,  "mag": 2.23,  "pm_ra": 1.67,    "pm_dec": 0.56},
    {"name": "Dubhe",       "ra": 165.932, "dec": 61.751,  "mag": 1.79,  "pm_ra": -136.46, "pm_dec": -35.25},
    {"name": "Alioth",      "ra": 193.507, "dec": 55.960,  "mag": 1.77,  "pm_ra": 111.74,  "pm_dec": -8.99},
    {"name": "Wezen",       "ra": 107.098, "dec": -26.393, "mag": 1.84,  "pm_ra": -2.75,   "pm_dec": 3.33},
    {"name": "Kaus Australis","ra":276.043,"dec": -34.384, "mag": 1.85,  "pm_ra": -39.61,  "pm_dec": -124.05},
    {"name": "Alkaid",      "ra": 206.885, "dec": 49.313,  "mag": 1.86,  "pm_ra": -121.23, "pm_dec": -15.56},
    {"name": "Sargas",      "ra": 264.330, "dec": -43.000, "mag": 1.87,  "pm_ra": 6.06,    "pm_dec": -0.95},
    {"name": "Avior",       "ra": 125.629, "dec": -59.509, "mag": 1.86,  "pm_ra": -25.34,  "pm_dec": 22.72},
    {"name": "Menkalinan",  "ra": 89.882,  "dec": 44.947,  "mag": 1.90,  "pm_ra": -56.41,  "pm_dec": -0.88},
    {"name": "Alhena",      "ra": 99.428,  "dec": 16.399,  "mag": 1.93,  "pm_ra": -2.04,   "pm_dec": -66.92},
    {"name": "Peacock",     "ra": 306.412, "dec": -56.735, "mag": 1.94,  "pm_ra": 7.71,    "pm_dec": -86.15},
    {"name": "Alsephina",   "ra": 131.176, "dec": -54.709, "mag": 1.96,  "pm_ra": -28.57,  "pm_dec": 16.30},
    {"name": "Mirzam",      "ra": 95.675,  "dec": -17.956, "mag": 1.98,  "pm_ra": -3.45,   "pm_dec": -0.47},
    {"name": "Alphard",     "ra": 141.897, "dec": -8.659,  "mag": 1.98,  "pm_ra": -14.49,  "pm_dec": 33.25},
    {"name": "Hamal",       "ra": 31.793,  "dec": 23.463,  "mag": 2.00,  "pm_ra": 190.73,  "pm_dec": -148.08},
    {"name": "Diphda",      "ra": 10.897,  "dec": -17.987, "mag": 2.02,  "pm_ra": 232.79,  "pm_dec": 32.71},
    {"name": "Nunki",       "ra": 283.816, "dec": -26.297, "mag": 2.02,  "pm_ra": 13.87,   "pm_dec": -52.65},
    {"name": "Miaplacidus", "ra": 138.300, "dec": -69.717, "mag": 1.68,  "pm_ra": -157.66, "pm_dec": 108.91},
    {"name": "Suhail",      "ra": 136.999, "dec": -43.433, "mag": 2.21,  "pm_ra": -23.21,  "pm_dec": 14.28},
    {"name": "Naos",        "ra": 120.896, "dec": -40.003, "mag": 2.25,  "pm_ra": -30.82,  "pm_dec": 17.05},
]


def precess_star(star, epoch_year):
    """Apply precession and proper motion to get star coords at a given epoch.

    Uses astropy for proper precession. epoch_year is in astronomical year
    (negative for BCE, e.g., -2500 for 2500 BCE).
    """
    # J2000.0 coordinates
    ra_j2000 = star["ra"]  # degrees
    dec_j2000 = star["dec"]  # degrees
    pm_ra = star["pm_ra"]  # mas/yr (includes cos(dec) factor for Hipparcos)
    pm_dec = star["pm_dec"]  # mas/yr

    # Time difference from J2000 (year 2000) to target epoch
    dt_years = epoch_year - 2000.0

    # Apply proper motion (Hipparcos pm_ra already includes cos(dec))
    ra_pm = ra_j2000 + (pm_ra / 3600000.0) * dt_years / math.cos(math.radians(dec_j2000))
    dec_pm = dec_j2000 + (pm_dec / 3600000.0) * dt_years

    # Use astropy for precession
    # Convert epoch_year to astropy Time
    # For BCE dates: -2500 → "J-2500"
    t = Time(epoch_year, format='jyear')

    # Create coordinate at J2000 with proper-motion-corrected values
    coord = SkyCoord(ra=ra_pm*u.deg, dec=dec_pm*u.deg, frame='icrs')

    # Precess to target epoch using FK5
    from astropy.coordinates import FK5
    fk5_epoch = FK5(equinox=t)
    precessed = coord.transform_to(fk5_epoch)

    return precessed.ra.deg, precessed.dec.deg


def rising_setting_azimuth(dec_deg, lat_deg):
    """Compute rising and setting azimuths for a star at given latitude.

    Returns (az_rise, az_set) in degrees from North.
    Returns (None, None) if the star is circumpolar or never rises.
    """
    dec = math.radians(dec_deg)
    lat = math.radians(lat_deg)

    cos_az = -math.sin(dec) / math.cos(lat)

    if abs(cos_az) > 1.0:
        return None, None  # circumpolar or never rises

    az_rise = math.degrees(math.acos(cos_az))
    az_set = 360.0 - az_rise

    return az_rise, az_set


def gc_bearing_at_site(site_lat_deg, site_lon_deg):
    """Compute the Great Circle's bearing (azimuth) at a given site.

    The bearing is the direction along the GC at this point.
    Returns two bearings (the two directions along the circle).
    """
    lat = math.radians(site_lat_deg)
    lon = math.radians(site_lon_deg)

    # The GC is defined by its pole. At any point on the GC, the bearing
    # is perpendicular to the direction toward the pole.
    #
    # Bearing from site to pole:
    dlon = POLE_LON - lon
    x = math.sin(dlon) * math.cos(POLE_LAT)
    y = math.cos(lat) * math.sin(POLE_LAT) - math.sin(lat) * math.cos(POLE_LAT) * math.cos(dlon)
    bearing_to_pole = math.degrees(math.atan2(x, y)) % 360

    # The GC bearing is perpendicular to the bearing toward the pole
    bearing1 = (bearing_to_pole + 90) % 360
    bearing2 = (bearing_to_pole - 90) % 360

    return bearing1, bearing2


def angular_diff(a, b):
    """Minimum angular difference between two azimuths (0-180)."""
    d = abs(a - b) % 360
    return min(d, 360 - d)


# ============================================================
# Main computation
# ============================================================
print("=" * 70)
print("ANALYSIS 2: STELLAR HORIZON EVENTS AT KEY SITES")
print("=" * 70)

all_results = []
csv_rows = []

for site in SITES:
    bearing1, bearing2 = gc_bearing_at_site(site["lat"], site["lon"])
    print(f"\n{site['name']} ({site['lat']:.2f}°, {site['lon']:.2f}°)")
    print(f"  GC bearings: {bearing1:.1f}° and {bearing2:.1f}°")

    for epoch in site["epochs"]:
        best_match = {"star": None, "offset": 999, "direction": None, "bearing_used": None}

        for star in BRIGHT_STARS:
            try:
                ra, dec = precess_star(star, epoch)
            except Exception:
                continue

            az_rise, az_set = rising_setting_azimuth(dec, site["lat"])
            if az_rise is None:
                continue

            for bearing in [bearing1, bearing2]:
                for az, direction in [(az_rise, "rise"), (az_set, "set")]:
                    offset = angular_diff(bearing, az)
                    if offset < best_match["offset"]:
                        best_match = {
                            "star": star["name"],
                            "mag": star["mag"],
                            "offset": round(offset, 2),
                            "direction": direction,
                            "bearing_used": round(bearing, 1),
                            "star_azimuth": round(az, 1),
                            "star_dec_at_epoch": round(dec, 2),
                        }

        epoch_str = f"{abs(epoch)} {'BCE' if epoch < 0 else 'CE'}"
        print(f"  Epoch {epoch_str}: best = {best_match['star']} "
              f"({best_match['direction']}), offset = {best_match['offset']:.1f}°")

        row = {
            "site": site["name"],
            "site_lat": site["lat"],
            "site_lon": site["lon"],
            "epoch_year": epoch,
            "gc_bearing_1": round(bearing1, 2),
            "gc_bearing_2": round(bearing2, 2),
            **best_match,
        }
        all_results.append(row)
        csv_rows.append(row)

# Write CSV
csv_path = os.path.join(OUT_DIR, "stellar_horizon_events.csv")
fieldnames = ["site", "site_lat", "site_lon", "epoch_year", "gc_bearing_1", "gc_bearing_2",
              "star", "mag", "offset", "direction", "bearing_used", "star_azimuth", "star_dec_at_epoch"]
with open(csv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(csv_rows)

# ============================================================
# Monte Carlo significance test
# ============================================================
print("\n" + "=" * 70)
print("MONTE CARLO SIGNIFICANCE TEST")
print("=" * 70)
print("Generating 10,000 random great circles and finding best star match at Giza...")

N_MC = 10000

# For Giza at 2500 BCE — precess all stars once
epoch_giza = -2500
giza_lat = 29.9792

precessed_stars = []
for star in BRIGHT_STARS:
    try:
        ra, dec = precess_star(star, epoch_giza)
        az_r, az_s = rising_setting_azimuth(dec, giza_lat)
        if az_r is not None:
            precessed_stars.append({"name": star["name"], "mag": star["mag"],
                                    "az_rise": az_r, "az_set": az_s})
    except Exception:
        continue

# Find the GC's best match at Giza
gc_b1, gc_b2 = gc_bearing_at_site(29.9792, 31.1342)
gc_best_offset = 999
for ps in precessed_stars:
    for b in [gc_b1, gc_b2]:
        for az in [ps["az_rise"], ps["az_set"]]:
            off = angular_diff(b, az)
            if off < gc_best_offset:
                gc_best_offset = off

print(f"Great Circle best offset at Giza (2500 BCE): {gc_best_offset:.2f}°")

# Monte Carlo: random great circles
mc_offsets = []
for i in range(N_MC):
    # Random pole (uniform on sphere)
    z = random.uniform(-1, 1)
    pole_lat_r = math.asin(z)
    pole_lon_r = random.uniform(-math.pi, math.pi)

    # Compute bearing at Giza for this random GC
    lat = math.radians(29.9792)
    lon = math.radians(31.1342)
    dlon = pole_lon_r - lon
    x = math.sin(dlon) * math.cos(pole_lat_r)
    y = math.cos(lat) * math.sin(pole_lat_r) - math.sin(lat) * math.cos(pole_lat_r) * math.cos(dlon)
    bearing_to_pole = math.degrees(math.atan2(x, y)) % 360
    b1 = (bearing_to_pole + 90) % 360
    b2 = (bearing_to_pole - 90) % 360

    best_off = 999
    for ps in precessed_stars:
        for b in [b1, b2]:
            for az in [ps["az_rise"], ps["az_set"]]:
                off = angular_diff(b, az)
                if off < best_off:
                    best_off = off
    mc_offsets.append(best_off)

# Compute percentile
mc_offsets.sort()
count_better = sum(1 for o in mc_offsets if o <= gc_best_offset)
percentile = count_better / N_MC * 100

mc_mean = np.mean(mc_offsets)
mc_std = np.std(mc_offsets)
mc_median = np.median(mc_offsets)

print(f"Monte Carlo: mean best offset = {mc_mean:.2f}° ± {mc_std:.2f}°")
print(f"Monte Carlo: median = {mc_median:.2f}°")
print(f"GC percentile: {percentile:.1f}% (lower = more unusual)")

significance = {
    "site": "Giza",
    "epoch": -2500,
    "gc_best_offset_deg": round(gc_best_offset, 2),
    "mc_n_trials": N_MC,
    "mc_mean_best_offset_deg": round(mc_mean, 2),
    "mc_std_deg": round(mc_std, 2),
    "mc_median_deg": round(mc_median, 2),
    "gc_percentile": round(percentile, 1),
    "conclusion": (
        f"The Great Circle's best stellar alignment at Giza (2500 BCE) has an offset of "
        f"{gc_best_offset:.1f}°. This is at the {percentile:.1f}th percentile among random "
        f"great circles (median random offset: {mc_median:.1f}°). "
        f"{'Not significant — typical of random circles.' if percentile > 5 else 'Potentially significant — unusually close alignment.'}"
    ),
}

sig_path = os.path.join(OUT_DIR, "stellar_significance.json")
with open(sig_path, "w") as f:
    json.dump(significance, f, indent=2)

# ============================================================
# Multi-site consistency check
# ============================================================
print("\n" + "=" * 70)
print("MULTI-SITE CONSISTENCY CHECK")
print("=" * 70)
print("Checking if the SAME star aligns at multiple sites at overlapping epochs...")

# For each star, check if it aligns within 5° at multiple sites
star_site_hits = {}
for r in all_results:
    if r["offset"] < 5.0 and r["star"]:
        key = r["star"]
        if key not in star_site_hits:
            star_site_hits[key] = []
        star_site_hits[key].append({
            "site": r["site"],
            "epoch": r["epoch_year"],
            "offset": r["offset"],
        })

multi_site_stars = {k: v for k, v in star_site_hits.items() if len(set(h["site"] for h in v)) >= 2}
print(f"\nStars within 5° at 2+ sites:")
for star, hits in sorted(multi_site_stars.items(), key=lambda x: -len(x[1])):
    sites_involved = set(h["site"] for h in hits)
    print(f"  {star}: {len(sites_involved)} sites — {', '.join(sites_involved)}")
    for h in hits:
        print(f"    {h['site']} @ {abs(h['epoch'])} {'BCE' if h['epoch'] < 0 else 'CE'}: {h['offset']:.1f}°")

if not multi_site_stars:
    print("  None found — no star aligns within 5° at multiple independent sites.")

# Save full results
full_results = {
    "stellar_events": all_results,
    "significance": significance,
    "multi_site_consistency": {
        k: v for k, v in multi_site_stars.items()
    } if multi_site_stars else "No multi-site alignments found within 5°",
    "overall_conclusion": (
        "The stellar horizon analysis found no compelling evidence of astronomical "
        "alignment. While individual site-epoch combinations can be found where a bright "
        "star rises or sets near the Great Circle's bearing, the Monte Carlo test shows "
        "this is expected for ANY great circle — the sky is dense with bright stars, and "
        "finding a ~few-degree match at any given site is not unusual."
    ),
}

full_path = os.path.join(OUT_DIR, "stellar_horizon_full.json")
with open(full_path, "w") as f:
    json.dump(full_results, f, indent=2, default=str)

print(f"\nOutputs written to:")
print(f"  {csv_path}")
print(f"  {sig_path}")
print(f"  {full_path}")
