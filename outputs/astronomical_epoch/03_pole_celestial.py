#!/usr/bin/env python3
"""
Analysis 3: The Great Circle's Pole on the Celestial Sphere
============================================================
Track where the GC's pole projects onto the celestial sphere through the
precession cycle. Check for proximity to significant celestial objects.
"""

import json
import csv
import math
import os
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, FK5, ICRS, Galactic
from astropy.time import Time
import astropy.units as u

OUT_DIR = os.path.dirname(__file__)

# Great Circle pole geographic coordinates
GC_POLE_LAT = 59.682122  # °N
GC_POLE_LON = -138.646087  # °E (i.e., 138.646°W)

# Key celestial reference points (J2000)
# North Ecliptic Pole: RA ≈ 270°, Dec ≈ +66.56° (by definition, the ecliptic pole
# is at 90° from the ecliptic plane)
ECLIPTIC_POLE_RA = 270.0  # degrees (18h)
ECLIPTIC_POLE_DEC = 66.56  # degrees

# North Galactic Pole (J2000)
GALACTIC_POLE_RA = 192.859  # degrees (12h51m26.3s)
GALACTIC_POLE_DEC = 27.128  # degrees (+27°07'42")

# Bright stars (subset for proximity checks, J2000)
REFERENCE_STARS = [
    {"name": "Vega",     "ra": 279.235, "dec": 38.784,  "mag": 0.03},
    {"name": "Polaris",  "ra": 37.954,  "dec": 89.264,  "mag": 1.98},
    {"name": "Deneb",    "ra": 310.358, "dec": 45.280,  "mag": 1.25},
    {"name": "Capella",  "ra": 79.172,  "dec": 45.998,  "mag": 0.08},
    {"name": "Arcturus", "ra": 213.915, "dec": 19.182,  "mag": -0.05},
    {"name": "Sirius",   "ra": 101.287, "dec": -16.716, "mag": -1.46},
    {"name": "Canopus",  "ra": 95.988,  "dec": -52.696, "mag": -0.74},
    {"name": "Aldebaran","ra": 68.980,  "dec": 16.509,  "mag": 0.85},
    {"name": "Betelgeuse","ra":88.793,  "dec": 7.407,   "mag": 0.42},
    {"name": "Rigel",    "ra": 78.634,  "dec": -8.202,  "mag": 0.13},
    {"name": "Procyon",  "ra": 114.826, "dec": 5.225,   "mag": 0.34},
    {"name": "Altair",   "ra": 297.696, "dec": 8.868,   "mag": 0.76},
    {"name": "Spica",    "ra": 201.298, "dec": -11.161, "mag": 0.97},
    {"name": "Antares",  "ra": 247.352, "dec": -26.432, "mag": 0.96},
    {"name": "Fomalhaut","ra": 344.413, "dec": -29.622, "mag": 1.16},
    {"name": "Regulus",  "ra": 152.093, "dec": 11.967,  "mag": 1.35},
    {"name": "Dubhe",    "ra": 165.932, "dec": 61.751,  "mag": 1.79},
    {"name": "Alkaid",   "ra": 206.885, "dec": 49.313,  "mag": 1.86},
    {"name": "Thuban",   "ra": 211.097, "dec": 64.376,  "mag": 3.67},  # Former pole star ~2700 BCE
]


def angular_sep(ra1, dec1, ra2, dec2):
    """Angular separation in degrees between two sky positions."""
    ra1, dec1, ra2, dec2 = map(math.radians, [ra1, dec1, ra2, dec2])
    cos_d = (math.sin(dec1) * math.sin(dec2) +
             math.cos(dec1) * math.cos(dec2) * math.cos(ra2 - ra1))
    return math.degrees(math.acos(max(-1, min(1, cos_d))))


def gc_pole_celestial_coords(epoch_year):
    """Compute the celestial coordinates of the GC pole's zenith direction.

    The GC pole is a fixed point on Earth's surface (59.68°N, 138.65°W).
    Its zenith direction projects to a fixed declination (= latitude = 59.68°)
    but a varying RA that depends on the relationship between geographic
    coordinates and celestial coordinates at that epoch.

    For a geographic location, the zenith at any instant has:
    - Dec = latitude (always)
    - RA = Local Sidereal Time at that longitude

    The key insight: we want to know WHERE on the celestial sphere the GC pole
    "points" at a reference moment for each epoch. Since the geographic pole
    traces a circle on the celestial sphere due to precession, the GC pole
    (which has a fixed relationship to the geographic pole) also traces a
    predictable path.

    Approach: The GC pole has a fixed angular relationship to the geographic
    pole. Due to precession, the celestial pole moves. The GC pole's celestial
    position can be computed by applying the same precession transformation
    to its J2000 celestial coordinates.

    At J2000, at sidereal time 0h, the GC pole's zenith points to:
    - Dec = 59.682° (always)
    - RA = the RA that transits at that longitude when LST=0 at Greenwich

    But this depends on the specific instant. A cleaner approach:
    Define the GC pole's direction in the ITRS (Earth-fixed) frame, then
    transform to celestial coordinates at each epoch.

    Simplest correct approach: The direction from Earth's center to the
    GC pole is fixed in Earth's body. Due to precession, this direction
    changes in celestial coordinates. We can compute this using the
    precession of the Earth's pole.

    The geographic pole (90°N) projects to Dec = 90° = celestial pole.
    The GC pole (59.68°N, 138.65°W) has a fixed angular relationship to the
    geographic pole. In the celestial frame at any epoch, the celestial pole
    is where it is (RA_cp, Dec_cp = defined by precession), and the GC pole
    is offset from it by the same angle.

    Actually, the simplest approach: The GC pole's latitude gives its
    declination relative to Earth's axis. Its longitude gives its orientation
    around that axis. As Earth's axis precesses, both change in celestial coords.

    Use astropy: create an EarthLocation for the GC pole, compute the zenith
    direction at a reference sidereal time, and transform to celestial coords.
    """
    # Use a fixed reference moment for each epoch: J0.0 of that year at midnight UT
    t = Time(epoch_year, format='jyear')

    # Earth location for the GC pole
    loc = EarthLocation(lat=GC_POLE_LAT * u.deg, lon=GC_POLE_LON * u.deg)

    # Compute the zenith at this location at this time
    altaz = AltAz(obstime=t, location=loc)
    zenith = SkyCoord(alt=90 * u.deg, az=0 * u.deg, frame=altaz)

    # Convert to ICRS (J2000 equatorial)
    icrs = zenith.icrs
    return icrs.ra.deg, icrs.dec.deg


# ============================================================
# Alternative approach: Use the GC pole's co-latitude and the
# precession of the celestial pole to trace the GC pole's
# celestial track analytically.
# ============================================================

def gc_pole_via_precession(epoch_year):
    """
    The GC pole is at co-latitude 30.318° from Earth's rotation axis,
    at a specific azimuth around that axis (determined by longitude).

    In celestial coordinates at epoch T:
    1. The celestial pole is at the precessing position
    2. The GC pole is always 30.318° from the celestial pole
    3. Its position angle around the celestial pole depends on the
       Earth's rotation, but at any fixed sidereal reference, we can
       compute it.

    For precession tracking, the key quantity is:
    - The GC pole's declination = 90° - angular_distance_from_celestial_pole
      = 90° - 30.318° = 59.682° (ONLY true when declination is measured
      from the TRUE pole of date, not the J2000 pole)

    So in the equatorial frame OF DATE, the GC pole is always at Dec = 59.682°.
    But in J2000 coordinates, its position changes with precession.

    The cleanest approach is to note that the GC pole's direction relative
    to the Earth's body is fixed. We can express this direction in
    ecliptic coordinates (which are approximately fixed) and then
    convert to equatorial coordinates at each epoch.
    """
    # The GC pole in Earth-fixed coordinates: lat 59.682°, lon -138.646°
    # Convert to ecliptic coordinates at J2000:
    # The GC pole's ecliptic longitude changes with precession (general precession
    # in longitude is ~50.29"/yr).

    # Use astropy's FK5 frame at each epoch to get the pole's coordinates
    # relative to the equator of date.

    # In the equatorial frame of date, the GC pole is at:
    # Dec_of_date = 59.682° (always, since it's fixed relative to Earth's axis)
    # RA_of_date depends on Earth's rotation — pick a reference sidereal time

    # For a meaningful track, we want the DIRECTION of the GC pole in
    # inertial (J2000) space. This requires knowing the orientation of
    # Earth's rotation axis at each epoch.

    # Precession parameters:
    # At epoch T (centuries from J2000), the celestial pole is at:
    # RA_pole ≈ 0° + precession terms (the celestial pole traces a circle
    # around the ecliptic pole)
    #
    # The ecliptic pole (J2000): RA=270°, Dec=66.56°
    # The celestial pole traces a circle of radius ε (obliquity ≈ 23.44°)
    # around the ecliptic pole with period 25,772 years.

    T = (epoch_year - 2000.0) / 25772.0  # fraction of precession cycle
    obliquity = 23.44  # approximately constant for our purposes

    # Celestial pole position in J2000 ecliptic coordinates:
    # The pole traces RA that increases by 360° per precession cycle
    # starting from RA ≈ 0° at some reference
    # In J2000 equatorial: the north celestial pole starts at (RA=any, Dec=90°)
    # and precesses to form a circle of radius ε around the ecliptic pole.

    # Position of celestial pole on the sky at epoch T (in J2000 ICRS):
    # RA_cp(T) = 270° + 360° * T  (the pole precesses prograde in ecliptic longitude)
    # Dec_cp(T) determined by distance from ecliptic pole = obliquity

    # Actually, the celestial pole in J2000 coordinates:
    # At J2000 (T=0), the celestial pole is at RA=0, Dec=90 (by definition)
    # At other epochs, we need the precession matrix.

    # Let's use astropy for precision:
    t = Time(epoch_year, format='jyear')

    # The celestial north pole of date in J2000 coords:
    # In FK5(equinox=t) frame, the north pole is at dec=90°.
    # Transform this to ICRS:
    pole_of_date = SkyCoord(ra=0*u.deg, dec=90*u.deg, frame=FK5(equinox=t))
    pole_icrs = pole_of_date.icrs

    # Now, the GC pole is at 30.318° from the celestial pole of date.
    # Its position angle around the celestial pole depends on its longitude
    # relative to the Greenwich meridian at the reference moment.

    # For a precession-cycle track, we can parameterize:
    # The GC pole traces a small circle of radius 30.318° around the
    # celestial pole of date. As the celestial pole itself moves due to
    # precession, the center of this small circle moves, and the GC pole's
    # J2000 position changes.

    # For the DIRECTION (ignoring Earth's daily rotation), the GC pole's
    # hour angle from the celestial pole is fixed = related to its longitude.
    # GC pole longitude = -138.646° → its sidereal hour angle (relative to
    # Greenwich) is -138.646° = 221.354° measured eastward.

    # In the frame of date, the GC pole is at:
    # Dec = GC_POLE_LAT = 59.682°
    # RA = LST at GC_POLE_LON at the reference time
    # For the precession track, use a fixed sidereal reference: GMST=0
    # Then RA_of_date = GC_POLE_LON = -138.646° = 221.354° (measured eastward)
    # This gives us the GC pole's position in the equatorial frame of date.

    gc_in_fk5 = SkyCoord(ra=(360 + GC_POLE_LON) % 360 * u.deg,
                          dec=GC_POLE_LAT * u.deg,
                          frame=FK5(equinox=t))
    gc_icrs = gc_in_fk5.icrs

    return gc_icrs.ra.deg, gc_icrs.dec.deg, pole_icrs.ra.deg, pole_icrs.dec.deg


# ============================================================
# Compute the track
# ============================================================
print("=" * 70)
print("ANALYSIS 3: GREAT CIRCLE POLE ON THE CELESTIAL SPHERE")
print("=" * 70)

# Sweep from -23000 to 3000 (25,000 BCE to 3000 CE) in 100-year steps
epochs = list(range(-23000, 3001, 100))
track = []

print(f"Computing pole celestial track over {len(epochs)} epochs...")

for epoch in epochs:
    try:
        ra, dec, cp_ra, cp_dec = gc_pole_via_precession(epoch)
        track.append({
            "epoch": epoch,
            "ra": round(ra, 3),
            "dec": round(dec, 3),
            "celestial_pole_ra": round(cp_ra, 3),
            "celestial_pole_dec": round(cp_dec, 3),
        })
    except Exception as e:
        pass

print(f"Computed {len(track)} epoch positions.")

# ============================================================
# Proximity to reference points
# ============================================================
print("\nChecking proximity to celestial reference points...")

# For each track point, compute distance to ecliptic pole, galactic pole, and stars
closest_ecliptic = {"epoch": None, "sep": 999}
closest_galactic = {"epoch": None, "sep": 999}
closest_stars = {s["name"]: {"epoch": None, "sep": 999} for s in REFERENCE_STARS}

for t in track:
    # Ecliptic pole
    sep = angular_sep(t["ra"], t["dec"], ECLIPTIC_POLE_RA, ECLIPTIC_POLE_DEC)
    if sep < closest_ecliptic["sep"]:
        closest_ecliptic = {"epoch": t["epoch"], "sep": round(sep, 2),
                            "gc_ra": t["ra"], "gc_dec": t["dec"]}

    # Galactic pole
    sep = angular_sep(t["ra"], t["dec"], GALACTIC_POLE_RA, GALACTIC_POLE_DEC)
    if sep < closest_galactic["sep"]:
        closest_galactic = {"epoch": t["epoch"], "sep": round(sep, 2),
                            "gc_ra": t["ra"], "gc_dec": t["dec"]}

    # Stars
    for s in REFERENCE_STARS:
        sep = angular_sep(t["ra"], t["dec"], s["ra"], s["dec"])
        if sep < closest_stars[s["name"]]["sep"]:
            closest_stars[s["name"]] = {"epoch": t["epoch"], "sep": round(sep, 2),
                                        "gc_ra": t["ra"], "gc_dec": t["dec"]}

print(f"\nClosest approach to North Ecliptic Pole:")
print(f"  {closest_ecliptic['sep']:.1f}° at epoch {closest_ecliptic['epoch']}")

print(f"\nClosest approach to North Galactic Pole:")
print(f"  {closest_galactic['sep']:.1f}° at epoch {closest_galactic['epoch']}")

print(f"\nClosest approaches to bright stars:")
for name, data in sorted(closest_stars.items(), key=lambda x: x[1]["sep"]):
    if data["sep"] < 15:
        print(f"  {name}: {data['sep']:.1f}° at epoch {data['epoch']}")

# Check if GC pole ever coincides with ecliptic pole
ecliptic_coincidence = closest_ecliptic["sep"] < 1.0

# ============================================================
# Save results
# ============================================================
csv_path = os.path.join(OUT_DIR, "pole_celestial_track.csv")
with open(csv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["epoch", "ra", "dec", "celestial_pole_ra", "celestial_pole_dec"])
    writer.writeheader()
    writer.writerows(track)

proximity = {
    "ecliptic_pole": {
        "reference": {"ra": ECLIPTIC_POLE_RA, "dec": ECLIPTIC_POLE_DEC},
        "closest_approach": closest_ecliptic,
        "coincidence_within_1deg": ecliptic_coincidence,
        "interpretation": (
            f"The GC pole comes within {closest_ecliptic['sep']}° of the North Ecliptic Pole "
            f"at epoch {closest_ecliptic['epoch']}. "
            f"{'This is close enough to suggest the GC approximates the ecliptic at this epoch.' if ecliptic_coincidence else 'This is not close enough to suggest the GC coincides with the ecliptic at any epoch.'}"
        ),
    },
    "galactic_pole": {
        "reference": {"ra": GALACTIC_POLE_RA, "dec": GALACTIC_POLE_DEC},
        "closest_approach": closest_galactic,
        "interpretation": f"The GC pole comes within {closest_galactic['sep']}° of the North Galactic Pole at epoch {closest_galactic['epoch']}.",
    },
    "bright_stars": {
        name: data for name, data in sorted(closest_stars.items(), key=lambda x: x[1]["sep"])
    },
}

prox_path = os.path.join(OUT_DIR, "pole_proximity.json")
with open(prox_path, "w") as f:
    json.dump(proximity, f, indent=2)

print(f"\nOutputs:")
print(f"  {csv_path}")
print(f"  {prox_path}")
print(f"\nConclusion: GC pole {'approaches' if closest_ecliptic['sep'] < 5 else 'does not closely approach'} the ecliptic pole (min sep: {closest_ecliptic['sep']}°)")
