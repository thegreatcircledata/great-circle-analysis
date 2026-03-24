#!/usr/bin/env python3
"""
Polynesian Navigation & Sacred Site Alignment — Master Data Module
===================================================================
Directive 09, Script 00

Consolidates all Polynesian sacred site data, voyaging routes, and
colonization chronology into a single importable module.

Sources:
- Easter Island ahu: DiNapoli, Lipo & Hunt (2019, 2020); Van Tilburg EISP
- Other Polynesian sites: UNESCO, NPS, archaeological surveys, OSM
- Voyaging routes: PVS, Lewis (1972), Gladwin (1970), Fitzpatrick & Callaghan (2009)
- Colonization dates: Wilmshurst et al. (2011), Kirch (2000)
"""

import math
import sys
import os

# ─── Import Easter Island ahu from existing dataset ───────────────────
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from easter_island_ahu_data import AHU_DATA

# ─── Import compiled Polynesian datasets ──────────────────────────────
from polynesian_sacred_sites import POLYNESIAN_SACRED_SITES
from polynesian_voyaging_data import (
    VOYAGING_ROUTES,
    COLONIZATION_CHRONOLOGY,
    POLYNESIAN_TRIANGLE,
)

# ─── Great Circle constants (consistent with all prior scripts) ───────
POLE_LAT = 59.682122
POLE_LON = -138.646087
R_EARTH = 6371.0
QUARTER_CIRC = R_EARTH * math.pi / 2  # ~10,007.5 km
N_MC = 10000

# ─── Distance thresholds for Pacific analysis ────────────────────────
THRESHOLDS_KM = [50, 100, 200, 500]

# ─── Merge all sacred sites into a unified list ──────────────────────
ALL_SACRED_SITES = []

# Easter Island ahu
for ahu in AHU_DATA:
    ALL_SACRED_SITES.append({
        "name": ahu["name"],
        "lat": ahu["lat"],
        "lon": ahu["lon"],
        "island": "Easter Island (Rapa Nui)",
        "island_group": "Easter Island",
        "site_type": "ahu",
        "notes": f"Platform {ahu['platform_length_m']}m, {ahu['n_moai']} moai",
    })

# All other Polynesian sites
for site in POLYNESIAN_SACRED_SITES:
    ALL_SACRED_SITES.append(site)


# ─── Major Polynesian island groups with centroids ───────────────────
ISLAND_GROUPS = [
    {"name": "Easter Island (Rapa Nui)", "lat": -27.11, "lon": -109.35},
    {"name": "Pitcairn Islands", "lat": -25.07, "lon": -130.10},
    {"name": "Henderson Island", "lat": -24.37, "lon": -128.32},
    {"name": "Mangareva / Gambier", "lat": -23.13, "lon": -134.97},
    {"name": "Marquesas Islands", "lat": -9.43, "lon": -140.07},
    {"name": "Tuamotu Archipelago", "lat": -17.00, "lon": -144.00},
    {"name": "Society Islands", "lat": -17.32, "lon": -149.75},
    {"name": "Cook Islands", "lat": -19.00, "lon": -159.50},
    {"name": "Samoa", "lat": -13.76, "lon": -172.10},
    {"name": "Tonga", "lat": -20.42, "lon": -175.20},
    {"name": "Fiji", "lat": -17.71, "lon": 178.07},
    {"name": "Hawaiian Islands", "lat": 20.79, "lon": -156.33},
    {"name": "New Zealand (Aotearoa)", "lat": -41.27, "lon": 174.78},
    {"name": "Niue", "lat": -19.05, "lon": -169.87},
    {"name": "Tokelau", "lat": -9.17, "lon": -171.82},
    {"name": "Tuvalu", "lat": -7.47, "lon": 179.20},
    {"name": "Wallis & Futuna", "lat": -13.28, "lon": -176.18},
    {"name": "Rapa Iti", "lat": -27.62, "lon": -144.34},
    {"name": "Line Islands", "lat": 1.87, "lon": -157.47},
    {"name": "Phoenix Islands", "lat": -3.85, "lon": -171.73},
]


# ─── Polynesian star compass directions ──────────────────────────────
# Traditional Hawaiian/Tahitian star compass with approximate declinations
# These are the bearing sectors used in Polynesian wayfinding
STAR_COMPASS = {
    # Rising positions (eastern horizon) at tropical latitudes
    "Hoku-pa'a (Polaris)": {"bearing": 0, "notes": "True north reference"},
    "Na Hiku (Big Dipper)": {"bearing": 15, "notes": "NNE rising"},
    "'A'a (Sirius)": {"bearing": 107, "notes": "ESE rising at 20°N"},
    "Hokule'a (Arcturus)": {"bearing": 67, "notes": "ENE rising, zenith star of Hawaii"},
    "Ke Ali'i o Kona (Canopus)": {"bearing": 153, "notes": "SSE rising"},
    "Newe (Southern Cross)": {"bearing": 170, "notes": "S at transit"},
    "Hikianalia (Spica)": {"bearing": 97, "notes": "E rising"},
    "Ka Makau Nui (Scorpius tail)": {"bearing": 135, "notes": "SE rising"},
    "Hanaiakamalama (S. Cross)": {"bearing": 160, "notes": "SSE rising at 20°N"},
}


# ─── Geometry helpers (consistent with prior scripts) ─────────────────

def haversine_km(lat1, lon1, lat2, lon2):
    """Haversine distance in km between two points."""
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlat = la2 - la1
    dlon = lo2 - lo1
    a = math.sin(dlat/2)**2 + math.cos(la1)*math.cos(la2)*math.sin(dlon/2)**2
    return R_EARTH * 2 * math.asin(math.sqrt(min(1.0, a)))


def dist_from_great_circle(lat, lon):
    """Distance (km) from a point to the Great Circle defined by the pole."""
    dist_to_pole = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(dist_to_pole - QUARTER_CIRC)


def bearing_of_great_circle_at_point(point_lat, point_lon):
    """
    Compute the bearing of the Great Circle at a given point.
    Returns two bearings (the circle goes both ways), both in [0, 360).
    """
    plat = math.radians(POLE_LAT)
    plon = math.radians(POLE_LON)
    qlat = math.radians(point_lat)
    qlon = math.radians(point_lon)

    dlon = plon - qlon
    x = math.sin(dlon) * math.cos(plat)
    y = (math.cos(qlat) * math.sin(plat) -
         math.sin(qlat) * math.cos(plat) * math.cos(dlon))
    bearing_to_pole = math.degrees(math.atan2(x, y)) % 360

    bearing1 = (bearing_to_pole + 90) % 360
    bearing2 = (bearing_to_pole - 90) % 360
    return bearing1, bearing2


def initial_bearing(lat1, lon1, lat2, lon2):
    """Initial bearing (forward azimuth) from point 1 to point 2."""
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlon = lo2 - lo1
    x = math.sin(dlon) * math.cos(la2)
    y = math.cos(la1)*math.sin(la2) - math.sin(la1)*math.cos(la2)*math.cos(dlon)
    return math.degrees(math.atan2(x, y)) % 360


# ═══════════════════════════════════════════════════════════════════════
# Self-test / summary
# ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 70)
    print("POLYNESIAN NAVIGATION DATA MODULE — SUMMARY")
    print("=" * 70)

    print(f"\nGreat Circle pole: ({POLE_LAT}, {POLE_LON})")
    print(f"Quarter-circumference: {QUARTER_CIRC:.1f} km")

    print(f"\nSacred sites: {len(ALL_SACRED_SITES)} total")
    groups = {}
    for s in ALL_SACRED_SITES:
        g = s["island_group"]
        groups[g] = groups.get(g, 0) + 1
    for g, n in sorted(groups.items(), key=lambda x: -x[1]):
        print(f"  {g:30s}: {n:3d} sites")

    print(f"\nIsland groups: {len(ISLAND_GROUPS)}")
    print(f"Voyaging routes: {len(VOYAGING_ROUTES)}")
    print(f"Colonization entries: {len(COLONIZATION_CHRONOLOGY)}")

    # Quick distance check: Easter Island should be ~0 km from circle
    ei_dist = dist_from_great_circle(-27.11, -109.35)
    print(f"\nEaster Island distance to Great Circle: {ei_dist:.1f} km")

    print(f"\nData module ready for import.")
