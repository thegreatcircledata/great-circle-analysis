#!/usr/bin/env python3
"""
Exhaustive Triplet Scan — Resolution Test 1
============================================
Among ALL possible triplets of famous archaeological sites,
where does Giza–Nazca–Easter Island rank for combined
geometric fit + monument-settlement divergence?

Steps:
  1. Build curated list of ~200 famous sites (Hancock + Pleiades major)
  2. For each triplet, compute best-fit great circle pole via SVD
  3. Compute geometric fit (max & mean distance of anchors to their GC)
  4. Compute Poisson-based D (monument–settlement divergence) along each GC
  5. Rank Giza–Nazca–Easter Island on D, fit, and combined score
"""

import csv
import json
import os
import sys
import time
import numpy as np
from itertools import combinations
from collections import Counter

import openpyxl

# ── Paths ─────────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HANCOCK_XLSX = os.path.join(BASE, "analysis", "hancock_archaeological_sites.xlsx")
PLEIADES_CSV = os.path.join(BASE, "data", "pleiades", "pleiades-places-latest.csv")
OUT_DIR = os.path.join(BASE, "outputs", "resolution_tests")
os.makedirs(OUT_DIR, exist_ok=True)

R_EARTH = 6371.0  # km
BAND_KM = 50.0    # half-width of corridor in km
BAND_RAD = BAND_KM / R_EARTH

# ── Coordinate conversions ────────────────────────────────────────────
def latlon_to_cart(lat, lon):
    """Convert lat/lon (degrees) to unit 3-vector."""
    la, lo = np.radians(lat), np.radians(lon)
    return np.array([np.cos(la) * np.cos(lo),
                     np.cos(la) * np.sin(lo),
                     np.sin(la)])

def cart_to_latlon(v):
    """Convert unit 3-vector to (lat, lon) degrees."""
    v = v / np.linalg.norm(v)
    lat = np.degrees(np.arcsin(np.clip(v[2], -1, 1)))
    lon = np.degrees(np.arctan2(v[1], v[0]))
    return lat, lon

# ── Great circle distance from a point to a pole's equator ───────────
def gc_distance_to_circle(site_cart, pole_cart):
    """
    Angular distance (radians) from a site to the great circle defined by pole.
    The great circle is the set of points at 90° from the pole.
    Distance = |acos(dot(site, pole))| - π/2
    """
    dot = np.clip(np.dot(site_cart, pole_cart), -1, 1)
    angular_dist_to_pole = np.arccos(dot)
    return abs(angular_dist_to_pole - np.pi / 2) * R_EARTH

# ── Load Hancock sites ────────────────────────────────────────────────
def load_hancock_sites():
    wb = openpyxl.load_workbook(HANCOCK_XLSX, read_only=True)
    ws = wb["Coordinates"]
    sites = []
    for i, row in enumerate(ws.iter_rows(min_row=2, values_only=True)):
        name, country, lat, lon, source = row[:5]
        if lat is None or lon is None:
            continue
        try:
            lat, lon = float(lat), float(lon)
        except (ValueError, TypeError):
            continue
        sites.append({"name": name, "lat": lat, "lon": lon, "source": "hancock"})
    wb.close()
    return sites

# ── Load Pleiades major sites ─────────────────────────────────────────
def load_pleiades_major_sites(hancock_names):
    """
    Select Pleiades sites that are:
    - featureType contains temple, pyramid, amphitheatre, or sanctuary
    - High connectivity (many connections to other places)
    - Have precise coordinates
    """
    monument_types = {"temple", "temple-2", "pyramid", "amphitheatre",
                      "amphitheatre-roman", "sanctuary", "sanctuary-2"}

    sites = []
    all_monuments = []  # for D calculation
    all_settlements = []  # for D calculation

    with open(PLEIADES_CSV, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            lat_str = row.get("reprLat", "")
            lon_str = row.get("reprLong", "")
            if not lat_str or not lon_str:
                continue
            try:
                lat = float(lat_str)
                lon = float(lon_str)
            except ValueError:
                continue

            ft = row.get("featureTypes", "").lower().strip()
            title = row.get("title", "")
            connections = row.get("connectsWith", "")
            n_connections = len(connections.split(",")) if connections.strip() else 0
            precision = row.get("locationPrecision", "")

            # Classify for D calculation
            ft_tokens = set(ft.replace(",", " ").split())
            is_monument = bool(ft_tokens & {"temple", "temple-2", "pyramid",
                                             "amphitheatre", "sanctuary", "church",
                                             "mosque", "synagogue", "fort",
                                             "monument", "tomb", "tumulus",
                                             "bath", "aqueduct", "bridge",
                                             "theatre", "stadium", "hippodrome",
                                             "arch", "lighthouse", "tower",
                                             "wall", "gate", "dam", "mine",
                                             "amphitheatre-roman", "sanctuary-2",
                                             "temple-2", "villa", "palace",
                                             "castellum", "castrum"})
            is_settlement = bool(ft_tokens & {"settlement", "settlement-modern",
                                               "urban", "village", "town",
                                               "city", "polis", "oppidum",
                                               "vicus", "municipium", "colonia",
                                               "civitas", "station", "port",
                                               "harbour", "settlement-modern"})

            cart = latlon_to_cart(lat, lon)
            if is_monument:
                all_monuments.append(cart)
            if is_settlement:
                all_settlements.append(cart)

            # Select as "famous" site for triplet anchors
            is_major_type = bool(ft_tokens & monument_types)
            if is_major_type and n_connections >= 1 and precision == "precise":
                # Deduplicate against Hancock
                name_lower = title.lower()
                if not any(h.lower() in name_lower or name_lower in h.lower()
                           for h in hancock_names):
                    sites.append({
                        "name": title,
                        "lat": lat,
                        "lon": lon,
                        "source": "pleiades",
                        "connections": n_connections,
                        "featureType": ft
                    })

    # Sort by connectivity, take top sites
    sites.sort(key=lambda s: s["connections"], reverse=True)

    monuments_arr = np.array(all_monuments) if all_monuments else np.zeros((0, 3))
    settlements_arr = np.array(all_settlements) if all_settlements else np.zeros((0, 3))

    return sites, monuments_arr, settlements_arr

# ── Poisson D computation (fast vectorized) ───────────────────────────
def compute_poisson_D(pole_cart, monuments_cart, settlements_cart):
    """
    Fast Poisson-approximation D for a great circle pole.

    For each category (monuments, settlements):
    - Count how many fall within BAND_KM of the great circle
    - Compute expected count under uniform Poisson model
    - Z = (observed - expected) / sqrt(expected)
    - D = Z_monuments - Z_settlements
    """
    # Angular distance from each point to the great circle
    # = |arccos(dot(point, pole)) - π/2|

    if len(monuments_cart) == 0 and len(settlements_cart) == 0:
        return 0.0, 0.0, 0.0, 0, 0

    # Fraction of sphere within ±BAND_RAD of equator
    # Area of band = 2π * R² * 2 * sin(band_rad) (for small band_rad ≈ 2π * R² * 2 * band_rad)
    # Fraction = 2 * sin(band_rad) ≈ 2 * band_rad (for small angles)
    frac = 2 * np.sin(BAND_RAD)  # fraction of sphere surface in band

    # Monument count
    if len(monuments_cart) > 0:
        dots_m = np.clip(monuments_cart @ pole_cart, -1, 1)
        ang_m = np.abs(np.arccos(dots_m) - np.pi / 2)
        n_mon = np.sum(ang_m < BAND_RAD)
        expected_mon = len(monuments_cart) * frac
        z_mon = (n_mon - expected_mon) / np.sqrt(max(expected_mon, 1))
    else:
        n_mon, z_mon = 0, 0.0

    # Settlement count
    if len(settlements_cart) > 0:
        dots_s = np.clip(settlements_cart @ pole_cart, -1, 1)
        ang_s = np.abs(np.arccos(dots_s) - np.pi / 2)
        n_set = np.sum(ang_s < BAND_RAD)
        expected_set = len(settlements_cart) * frac
        z_set = (n_set - expected_set) / np.sqrt(max(expected_set, 1))
    else:
        n_set, z_set = 0, 0.0

    D = z_mon - z_set
    return D, z_mon, z_set, int(n_mon), int(n_set)

# ── Main ──────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("EXHAUSTIVE TRIPLET SCAN — Resolution Test 1")
    print("=" * 60)

    # Step 1: Load Hancock sites
    print("\n[1/5] Loading Hancock sites...")
    hancock_sites = load_hancock_sites()
    hancock_names = [s["name"] for s in hancock_sites]
    print(f"  Hancock sites with coordinates: {len(hancock_sites)}")

    # Step 2: Load Pleiades major sites
    print("[2/5] Loading Pleiades major sites + monument/settlement arrays...")
    pleiades_sites, monuments_cart, settlements_cart = load_pleiades_major_sites(hancock_names)
    print(f"  Pleiades candidate sites: {len(pleiades_sites)}")
    print(f"  Pleiades monuments for D: {len(monuments_cart)}")
    print(f"  Pleiades settlements for D: {len(settlements_cart)}")

    # Step 3: Add globally famous sites not in Hancock or Pleiades
    # (Pleiades is Mediterranean-focused, so we need Asian/African/American sites)
    global_famous = [
        {"name": "Borobudur", "lat": -7.608, "lon": 110.204, "source": "curated"},
        {"name": "Prambanan", "lat": -7.752, "lon": 110.492, "source": "curated"},
        {"name": "Mohenjo-daro", "lat": 27.329, "lon": 68.138, "source": "curated"},
        {"name": "Harappa", "lat": 30.631, "lon": 72.868, "source": "curated"},
        {"name": "Sanchi Stupa", "lat": 23.479, "lon": 77.740, "source": "curated"},
        {"name": "Ellora Caves", "lat": 20.027, "lon": 75.179, "source": "curated"},
        {"name": "Ajanta Caves", "lat": 20.552, "lon": 75.700, "source": "curated"},
        {"name": "Hampi", "lat": 15.335, "lon": 76.460, "source": "curated"},
        {"name": "Sigiriya", "lat": 7.957, "lon": 80.760, "source": "curated"},
        {"name": "Bagan Temples", "lat": 21.172, "lon": 94.859, "source": "curated"},
        {"name": "Sukhothai", "lat": 17.017, "lon": 99.703, "source": "curated"},
        {"name": "Great Wall Badaling", "lat": 40.360, "lon": 116.014, "source": "curated"},
        {"name": "Terracotta Army Xian", "lat": 34.384, "lon": 109.278, "source": "curated"},
        {"name": "Forbidden City Beijing", "lat": 39.917, "lon": 116.397, "source": "curated"},
        {"name": "Temple of Heaven Beijing", "lat": 39.882, "lon": 116.407, "source": "curated"},
        {"name": "Gyeongju Tumuli", "lat": 35.834, "lon": 129.219, "source": "curated"},
        {"name": "Kofun Daisen (Nintoku)", "lat": 34.564, "lon": 135.488, "source": "curated"},
        {"name": "Ise Grand Shrine", "lat": 34.455, "lon": 136.726, "source": "curated"},
        {"name": "Great Zimbabwe", "lat": -20.274, "lon": 30.934, "source": "curated"},
        {"name": "Meroe Pyramids", "lat": 16.938, "lon": 33.749, "source": "curated"},
        {"name": "Aksum Stelae Field", "lat": 14.131, "lon": 38.720, "source": "curated"},
        {"name": "Lalibela Rock Churches", "lat": 12.032, "lon": 39.044, "source": "curated"},
        {"name": "Leptis Magna", "lat": 32.638, "lon": 14.289, "source": "curated"},
        {"name": "Volubilis", "lat": 34.074, "lon": -5.554, "source": "curated"},
        {"name": "Cahokia Monks Mound", "lat": 38.660, "lon": -90.062, "source": "curated"},
        {"name": "Poverty Point", "lat": 32.634, "lon": -91.408, "source": "curated"},
        {"name": "Mesa Verde", "lat": 37.184, "lon": -108.489, "source": "curated"},
        {"name": "Chaco Canyon Pueblo Bonito", "lat": 36.061, "lon": -107.961, "source": "curated"},
        {"name": "Caral", "lat": -10.893, "lon": -77.520, "source": "curated"},
        {"name": "Chan Chan", "lat": -8.106, "lon": -79.075, "source": "curated"},
        {"name": "Copan", "lat": 14.840, "lon": -89.141, "source": "curated"},
        {"name": "Monte Alban", "lat": 17.044, "lon": -96.768, "source": "curated"},
        {"name": "Palenque", "lat": 17.484, "lon": -92.046, "source": "curated"},
        {"name": "Uxmal", "lat": 20.359, "lon": -89.771, "source": "curated"},
        {"name": "El Tajin", "lat": 20.449, "lon": -97.378, "source": "curated"},
        {"name": "Persepolis", "lat": 29.935, "lon": 52.891, "source": "curated"},
        {"name": "Pasargadae", "lat": 30.194, "lon": 53.167, "source": "curated"},
        {"name": "Hattusa", "lat": 40.018, "lon": 34.615, "source": "curated"},
        {"name": "Troy Hisarlik", "lat": 39.957, "lon": 26.239, "source": "curated"},
        {"name": "Knossos", "lat": 35.298, "lon": 25.163, "source": "curated"},
        {"name": "Mycenae", "lat": 37.731, "lon": 22.757, "source": "curated"},
        {"name": "Delphi", "lat": 38.482, "lon": 22.501, "source": "curated"},
        {"name": "Olympia", "lat": 37.638, "lon": 21.630, "source": "curated"},
        {"name": "Ephesus", "lat": 37.941, "lon": 27.342, "source": "curated"},
        {"name": "Palmyra", "lat": 34.551, "lon": 38.269, "source": "curated"},
        {"name": "Petra", "lat": 30.329, "lon": 35.444, "source": "curated"},
        {"name": "Jerash", "lat": 32.281, "lon": 35.891, "source": "curated"},
        {"name": "Pompeii", "lat": 40.751, "lon": 14.487, "source": "curated"},
        {"name": "Colosseum Rome", "lat": 41.890, "lon": 12.492, "source": "curated"},
        {"name": "Parthenon Athens", "lat": 37.971, "lon": 23.727, "source": "curated"},
        {"name": "Stonehenge", "lat": 51.179, "lon": -1.826, "source": "curated"},
        {"name": "Avebury", "lat": 51.429, "lon": -1.854, "source": "curated"},
        {"name": "Newgrange", "lat": 53.695, "lon": -6.475, "source": "curated"},
        {"name": "Carnac Stones", "lat": 47.585, "lon": -3.077, "source": "curated"},
        {"name": "Skara Brae", "lat": 59.049, "lon": -3.343, "source": "curated"},
        {"name": "Nan Madol", "lat": 6.844, "lon": 158.334, "source": "curated"},
        {"name": "Tonga Haamonga", "lat": -21.139, "lon": -175.068, "source": "curated"},
    ]

    # Step 4: Combine and deduplicate
    # Keep all Hancock sites, add curated global sites, add top Pleiades sites to reach ~200
    target = 210
    n_pleiades_needed = max(0, target - len(hancock_sites) - len(global_famous))
    pleiades_add = pleiades_sites[:n_pleiades_needed]

    all_sites = hancock_sites + global_famous + pleiades_add
    print(f"  Pre-dedup: {len(hancock_sites)} Hancock + {len(global_famous)} curated + {len(pleiades_add)} Pleiades = {len(all_sites)}")

    # Remove near-duplicates (within 10 km)
    deduped = []
    for s in all_sites:
        s["cart"] = latlon_to_cart(s["lat"], s["lon"])
        is_dup = False
        for d in deduped:
            dist = np.arccos(np.clip(np.dot(s["cart"], d["cart"]), -1, 1)) * R_EARTH
            if dist < 10:
                is_dup = True
                break
        if not is_dup:
            deduped.append(s)

    all_sites = deduped
    print(f"\n  Combined site list (after dedup): {len(all_sites)}")

    # Identify Giza, Nazca, Easter Island indices
    giza_idx = nazca_idx = easter_idx = None
    for i, s in enumerate(all_sites):
        name = s["name"].lower()
        # Must match "giza" specifically (not "Cholula Great Pyramid")
        if "giza" in name and ("pyramid" in name or "khufu" in name or "plateau" in name):
            if giza_idx is None:
                giza_idx = i
        if "nazca" in name and "line" in name:
            if nazca_idx is None:
                nazca_idx = i
        if ("easter island" in name or "rapa nui" in name):
            if easter_idx is None:
                easter_idx = i

    print(f"  Giza index: {giza_idx} ({all_sites[giza_idx]['name'] if giza_idx is not None else 'NOT FOUND'})")
    print(f"  Nazca index: {nazca_idx} ({all_sites[nazca_idx]['name'] if nazca_idx is not None else 'NOT FOUND'})")
    print(f"  Easter Island index: {easter_idx} ({all_sites[easter_idx]['name'] if easter_idx is not None else 'NOT FOUND'})")

    if giza_idx is None or nazca_idx is None or easter_idx is None:
        # Add missing sites manually
        manual_adds = []
        if giza_idx is None:
            manual_adds.append({"name": "Great Pyramid of Giza", "lat": 29.9792, "lon": 31.1342, "source": "manual"})
        if nazca_idx is None:
            manual_adds.append({"name": "Nazca Lines", "lat": -14.735, "lon": -75.13, "source": "manual"})
        if easter_idx is None:
            manual_adds.append({"name": "Easter Island (Ahu Tongariki)", "lat": -27.1246, "lon": -109.2771, "source": "manual"})

        for s in manual_adds:
            s["cart"] = latlon_to_cart(s["lat"], s["lon"])
            all_sites.append(s)

        # Re-find indices
        for i, s in enumerate(all_sites):
            name = s["name"].lower()
            if ("giza" in name or "great pyramid" in name) and giza_idx is None:
                giza_idx = i
            if "nazca" in name and nazca_idx is None:
                nazca_idx = i
            if ("easter" in name or "rapa nui" in name) and easter_idx is None:
                easter_idx = i

        print(f"  After manual adds: {len(all_sites)} sites")
        print(f"  Giza: {giza_idx}, Nazca: {nazca_idx}, Easter: {easter_idx}")

    alison_triplet = tuple(sorted([giza_idx, nazca_idx, easter_idx]))

    N = len(all_sites)
    n_triplets = N * (N - 1) * (N - 2) // 6
    print(f"\n  Total sites: {N}")
    print(f"  Total triplets: {n_triplets:,}")

    # Step 4: Precompute cartesian coordinates array
    print("\n[3/5] Precomputing site cartesian coordinates...")
    carts = np.array([s["cart"] for s in all_sites])  # (N, 3)

    # Step 5: Run exhaustive triplet scan in BATCHED mode for speed
    print(f"[4/5] Running exhaustive triplet scan over {n_triplets:,} triplets...")

    # Generate all triplet indices
    print("  Generating triplet index array...")
    triplet_indices = np.array(list(combinations(range(N), 3)), dtype=np.int32)
    actual_n = len(triplet_indices)
    print(f"  {actual_n:,} triplets generated")

    # Find Alison triplet index
    alison_row = None
    for row_idx in range(actual_n):
        if tuple(triplet_indices[row_idx]) == alison_triplet:
            alison_row = row_idx
            break
    print(f"  Alison triplet at row: {alison_row}")

    # Precompute monument/settlement angular distances from ALL possible poles
    # Strategy: batch SVD, batch D computation

    # Fraction of sphere in band
    frac = 2 * np.sin(BAND_RAD)
    expected_mon = len(monuments_cart) * frac
    expected_set = len(settlements_cart) * frac
    sqrt_exp_mon = np.sqrt(max(expected_mon, 1))
    sqrt_exp_set = np.sqrt(max(expected_set, 1))

    # Pre-allocate
    all_D = np.zeros(actual_n, dtype=np.float32)
    all_fit_max = np.zeros(actual_n, dtype=np.float32)
    all_fit_mean = np.zeros(actual_n, dtype=np.float32)
    all_z_mon = np.zeros(actual_n, dtype=np.float32)
    all_z_set = np.zeros(actual_n, dtype=np.float32)
    all_n_mon = np.zeros(actual_n, dtype=np.int32)
    all_n_set = np.zeros(actual_n, dtype=np.int32)

    # Store top 200 by combined score
    top_combined = []

    # FAST D: instead of arccos, use the fact that |dot(point, pole)| < sin(BAND_RAD)
    # means the point is within BAND_RAD of the great circle.
    # This is because angular distance to GC = |arccos(dot) - π/2| < BAND_RAD
    # ⟺ cos(π/2 + BAND_RAD) < dot < cos(π/2 - BAND_RAD)
    # ⟺ -sin(BAND_RAD) < dot < sin(BAND_RAD)
    # ⟺ |dot| < sin(BAND_RAD)
    sin_band = np.sin(BAND_RAD)

    # Transpose monument/settlement arrays for fast matmul: (3, M) and (3, S)
    mon_T = monuments_cart.T.copy()  # (3, n_mon_total)
    set_T = settlements_cart.T.copy()  # (3, n_set_total)

    BATCH = 2000
    t0 = time.time()

    for batch_start in range(0, actual_n, BATCH):
        batch_end = min(batch_start + BATCH, actual_n)
        batch_idx = triplet_indices[batch_start:batch_end]  # (B, 3)
        B = len(batch_idx)

        # Compute poles for entire batch
        # For great circles: we need the normal to the best-fit plane THROUGH THE ORIGIN.
        # Cross product v1×v2 gives normal to plane through origin containing v1,v2.
        # For 3 points, the best-fit great circle pole is the right singular vector
        # with smallest singular value. But SVD is slow in a loop.
        #
        # Shortcut: for 3 unit vectors, the normal to the plane through origin
        # that best fits them is the eigenvector of A^T A with smallest eigenvalue,
        # where A = [v1; v2; v3] (3×3).
        # A^T A = sum of outer products. The smallest eigenvector of this 3×3 matrix
        # is our pole.
        #
        # Even faster: use cross product v1 × v2, then project v3 and measure residual.
        # The pole of the great circle through v1 and v2 is cross(v1, v2).
        # The distance of v3 from this great circle is the fit metric.
        # This is actually what we want: how close are 3 points to a great circle?
        # Use SVD batch via loop (still fast enough with numpy).

        v1 = carts[batch_idx[:, 0]]  # (B, 3)
        v2 = carts[batch_idx[:, 1]]  # (B, 3)
        v3 = carts[batch_idx[:, 2]]  # (B, 3)

        # For each triplet, compute pole via eigendecomposition of A^T A (3×3)
        # where A = [v1; v2; v3]. The eigenvector with smallest eigenvalue = pole.
        # A^T A = v1⊗v1 + v2⊗v2 + v3⊗v3
        # Batch: compute all M matrices at once using einsum
        # M[b] = v1[b]⊗v1[b] + v2[b]⊗v2[b] + v3[b]⊗v3[b]  shape (B, 3, 3)
        M = np.einsum('bi,bj->bij', v1, v1) + np.einsum('bi,bj->bij', v2, v2) + np.einsum('bi,bj->bij', v3, v3)
        # Batch eigendecomposition
        eigenvalues, eigenvectors = np.linalg.eigh(M)  # eigenvalues sorted ascending
        # Pole = eigenvector with smallest eigenvalue (index 0)
        poles = eigenvectors[:, :, 0]  # (B, 3)
        degenerate = eigenvalues[:, 0] < 1e-20

        # Convention: positive z
        flip = poles[:, 2] < 0
        poles[flip] *= -1

        # Geometric fit for each triplet: |dot(vi, pole)| * R_EARTH
        # Distance from point to great circle = arcsin(|dot(point, pole)|) * R
        # For small angles, ≈ |dot| * R
        dot1 = np.abs(np.sum(v1 * poles, axis=1))  # (B,)
        dot2 = np.abs(np.sum(v2 * poles, axis=1))
        dot3 = np.abs(np.sum(v3 * poles, axis=1))
        # Use arcsin for accuracy
        d1_km = np.arcsin(np.clip(dot1, 0, 1)) * R_EARTH
        d2_km = np.arcsin(np.clip(dot2, 0, 1)) * R_EARTH
        d3_km = np.arcsin(np.clip(dot3, 0, 1)) * R_EARTH

        fit_max_batch = np.maximum(np.maximum(d1_km, d2_km), d3_km)
        fit_mean_batch = (d1_km + d2_km + d3_km) / 3.0

        # D computation: poles @ mon_T gives (B, n_mon_total)
        # Count |dot| < sin_band for each pole
        dots_mon = poles @ mon_T  # (B, n_mon_total)
        n_mon_batch = np.sum(np.abs(dots_mon) < sin_band, axis=1)  # (B,)
        z_mon_batch = (n_mon_batch - expected_mon) / sqrt_exp_mon

        dots_set = poles @ set_T  # (B, n_set_total)
        n_set_batch = np.sum(np.abs(dots_set) < sin_band, axis=1)  # (B,)
        z_set_batch = (n_set_batch - expected_set) / sqrt_exp_set

        D_batch = z_mon_batch - z_set_batch

        # Handle degenerate
        D_batch[degenerate] = 0.0
        fit_max_batch[degenerate] = 9999.0
        fit_mean_batch[degenerate] = 9999.0

        # Store
        sl = slice(batch_start, batch_end)
        all_D[sl] = D_batch
        all_fit_max[sl] = fit_max_batch
        all_fit_mean[sl] = fit_mean_batch
        all_z_mon[sl] = z_mon_batch
        all_z_set[sl] = z_set_batch
        all_n_mon[sl] = n_mon_batch
        all_n_set[sl] = n_set_batch

        # Combined scores
        combined_batch = D_batch / (fit_mean_batch + 1.0)

        # Track top 200
        for b in range(B):
            gi = batch_start + b
            combined = float(combined_batch[b])
            if len(top_combined) < 200 or combined > top_combined[-1][0]:
                i, j, k = int(batch_idx[b, 0]), int(batch_idx[b, 1]), int(batch_idx[b, 2])
                pole_lat, pole_lon = cart_to_latlon(poles[b])
                detail = {
                    "sites": (all_sites[i]["name"], all_sites[j]["name"], all_sites[k]["name"]),
                    "indices": (i, j, k),
                    "pole_lat": round(pole_lat, 4),
                    "pole_lon": round(pole_lon, 4),
                    "fit_max_km": round(float(fit_max_batch[b]), 2),
                    "fit_mean_km": round(float(fit_mean_batch[b]), 2),
                    "d1_km": round(float(d1_km[b]), 2),
                    "d2_km": round(float(d2_km[b]), 2),
                    "d3_km": round(float(d3_km[b]), 2),
                    "D_approx": round(float(D_batch[b]), 4),
                    "z_mon": round(float(z_mon_batch[b]), 4),
                    "z_set": round(float(z_set_batch[b]), 4),
                    "n_mon": int(n_mon_batch[b]),
                    "n_set": int(n_set_batch[b]),
                    "combined_score": round(combined, 4)
                }
                if len(top_combined) < 200:
                    top_combined.append((combined, gi, detail))
                    if len(top_combined) == 200:
                        top_combined.sort(key=lambda x: -x[0])
                else:
                    top_combined[-1] = (combined, gi, detail)
                    top_combined.sort(key=lambda x: -x[0])

        if (batch_start // BATCH) % 50 == 0:
            elapsed = time.time() - t0
            done = batch_end
            rate = done / max(elapsed, 0.01)
            eta = (actual_n - done) / max(rate, 1)
            print(f"  {done:>10,} / {actual_n:,} ({100*done/actual_n:.1f}%) "
                  f"— {elapsed:.0f}s elapsed, ~{eta:.0f}s remaining")

    elapsed = time.time() - t0
    print(f"\n  Scan complete: {actual_n:,} triplets in {elapsed:.1f}s ({actual_n/elapsed:.0f}/s)")

    # Print Alison triplet info
    alison_idx_in_results = alison_row
    if alison_row is not None:
        ai, aj, ak = triplet_indices[alison_row]
        print(f"\n  *** ALISON TRIPLET (row {alison_row}) ***")
        print(f"      Sites: {all_sites[ai]['name']}, {all_sites[aj]['name']}, {all_sites[ak]['name']}")
        print(f"      Fit: max={all_fit_max[alison_row]:.2f} km, mean={all_fit_mean[alison_row]:.2f} km")
        print(f"      D={all_D[alison_row]:.4f}, Z_mon={all_z_mon[alison_row]:.4f}, Z_set={all_z_set[alison_row]:.4f}")
        print(f"      Mon={all_n_mon[alison_row]}, Set={all_n_set[alison_row]}")

    idx = actual_n

    elapsed = time.time() - t0
    print(f"\n  Scan complete: {idx:,} triplets in {elapsed:.1f}s ({idx/elapsed:.0f}/s)")

    # Step 6: Compute ranks
    print("\n[5/5] Computing ranks and producing output...")

    if alison_idx_in_results is None:
        print("  WARNING: Alison triplet not found in scan! Searching manually...")
        # Find it by name
        for entry in top_combined:
            names = set(n.lower() for n in entry[2]["sites"])
            if any("giza" in n or "pyramid" in n for n in names) and \
               any("nazca" in n for n in names) and \
               any("easter" in n or "rapa" in n for n in names):
                alison_idx_in_results = entry[1]
                print(f"  Found in top combined list at idx {entry[1]}")
                break

    # Get Alison values
    if alison_idx_in_results is not None:
        alison_D = float(all_D[alison_idx_in_results])
        alison_fit_max = float(all_fit_max[alison_idx_in_results])
        alison_fit_mean = float(all_fit_mean[alison_idx_in_results])
        alison_combined = alison_D / (alison_fit_mean + 1.0)
    else:
        print("  FATAL: Cannot find Alison triplet!")
        alison_D = alison_fit_max = alison_fit_mean = alison_combined = None

    # Rank on D
    rank_D = int(np.sum(all_D[:idx] > alison_D)) + 1 if alison_D is not None else None

    # Rank on fit (lower is better)
    rank_fit = int(np.sum(all_fit_mean[:idx] < alison_fit_mean)) + 1 if alison_fit_mean is not None else None

    # Rank on combined score
    all_combined = all_D[:idx] / (all_fit_mean[:idx] + 1.0)
    rank_combined = int(np.sum(all_combined > alison_combined)) + 1 if alison_combined is not None else None

    # Count triplets with D > 8 AND fit < 10 km
    high_D_low_fit = int(np.sum((all_D[:idx] > 8) & (all_fit_max[:idx] < 10)))

    # Percentiles
    pct_D = round(100 * (1 - rank_D / idx), 4) if rank_D else None
    pct_fit = round(100 * (1 - rank_fit / idx), 4) if rank_fit else None
    pct_combined = round(100 * (1 - rank_combined / idx), 4) if rank_combined else None

    print(f"\n  RESULTS:")
    print(f"  Total triplets scanned: {idx:,}")
    print(f"  Alison D = {alison_D:.4f}, rank = {rank_D:,} ({pct_D}th percentile)")
    print(f"  Alison fit_mean = {alison_fit_mean:.2f} km, rank = {rank_fit:,} ({pct_fit}th percentile)")
    print(f"  Alison combined = {alison_combined:.4f}, rank = {rank_combined:,} ({pct_combined}th percentile)")
    print(f"  Triplets with D > 8 AND fit_max < 10 km: {high_D_low_fit}")

    # ── Output files ──────────────────────────────────────────────────

    # Top 100 by combined score
    top100 = [entry[2] for entry in top_combined[:100]]
    with open(os.path.join(OUT_DIR, "triplet_scan.json"), "w") as f:
        json.dump({
            "description": "Top 100 triplets by combined score (D / (mean_fit + 1))",
            "total_sites": N,
            "total_triplets": idx,
            "site_list": [{"name": s["name"], "lat": round(s["lat"], 4), "lon": round(s["lon"], 4),
                          "source": s["source"]} for s in all_sites],
            "top_100_combined": top100,
            "stats": {
                "D_mean": round(float(np.mean(all_D[:idx])), 4),
                "D_std": round(float(np.std(all_D[:idx])), 4),
                "D_max": round(float(np.max(all_D[:idx])), 4),
                "D_min": round(float(np.min(all_D[:idx])), 4),
                "fit_mean_median": round(float(np.median(all_fit_mean[:idx])), 2),
                "fit_mean_mean": round(float(np.mean(all_fit_mean[:idx])), 2),
                "combined_max": round(float(np.max(all_combined)), 4),
                "combined_mean": round(float(np.mean(all_combined)), 4),
                "high_D_low_fit_count": high_D_low_fit,
            }
        }, f, indent=2)

    # Alison triplet ranks
    ranks = {
        "description": "Giza-Nazca-Easter Island triplet percentile rankings",
        "alison_triplet": {
            "sites": ["Great Pyramid of Giza", "Nazca Lines", "Easter Island"],
            "D_approx": alison_D,
            "fit_mean_km": alison_fit_mean,
            "fit_max_km": alison_fit_max,
            "combined_score": alison_combined,
        },
        "ranks": {
            "D_rank": rank_D,
            "D_percentile": pct_D,
            "fit_rank": rank_fit,
            "fit_percentile": pct_fit,
            "combined_rank": rank_combined,
            "combined_percentile": pct_combined,
        },
        "context": {
            "total_triplets": idx,
            "total_sites": N,
            "high_D_and_low_fit_count": high_D_low_fit,
        }
    }
    with open(os.path.join(OUT_DIR, "triplet_ranks.json"), "w") as f:
        json.dump(ranks, f, indent=2)

    # Verdict
    if rank_combined is not None:
        # Normalize to equivalent rank at 1.3M scale
        equiv_rank_1_3M = int(rank_combined * (1_300_000 / idx))
        if equiv_rank_1_3M <= 10:
            interpretation = "Genuinely exceptional — collinearity is not a Bible code artifact"
        elif equiv_rank_1_3M <= 100:
            interpretation = "Notable but not unique — several other triplets are comparable"
        elif equiv_rank_1_3M <= 1000:
            interpretation = "Unremarkable — selection bias likely explains it"
        elif equiv_rank_1_3M <= 10000:
            interpretation = "Weak — many triplets score comparably or better"
        else:
            interpretation = "Fully explained by selection — the 'Bible code' explanation wins"
    else:
        interpretation = "UNABLE TO DETERMINE"
        equiv_rank_1_3M = None

    verdict_md = f"""# Triplet Scan Verdict — Resolution Test 1

## Summary
Scanned **{idx:,}** triplets from **{N}** famous archaeological sites
(Hancock list + Pleiades major temples/pyramids/sanctuaries/amphitheatres).

## Alison Triplet (Giza–Nazca–Easter Island)

| Metric | Value | Rank | Percentile |
|--------|-------|------|------------|
| D (monument–settlement divergence) | {alison_D:.4f} | {rank_D:,} of {idx:,} | {pct_D}% |
| Geometric fit (mean, km) | {alison_fit_mean:.2f} | {rank_fit:,} of {idx:,} | {pct_fit}% |
| Combined score D/(fit+1) | {alison_combined:.4f} | **{rank_combined:,}** of {idx:,} | **{pct_combined}%** |

## Key Statistics
- Triplets with D > 8 AND fit_max < 10 km: **{high_D_low_fit}**
- Global max combined score: {float(np.max(all_combined)):.4f}
- Mean combined score: {float(np.mean(all_combined)):.4f}

## Top 5 Triplets by Combined Score

| Rank | Sites | D | Fit (mean km) | Combined |
|------|-------|---|---------------|----------|
"""
    for i, entry in enumerate(top100[:5]):
        verdict_md += f"| {i+1} | {', '.join(entry['sites'])} | {entry['D_approx']:.2f} | {entry['fit_mean_km']:.2f} | {entry['combined_score']:.2f} |\n"

    verdict_md += f"""
## Verdict

| Alison rank | Interpretation |
|-------------|---------------|
| Top 10 (equiv ~1.3M) | Genuinely exceptional — collinearity is not a Bible code artifact |
| Top 100 | Notable but not unique — several other triplets are comparable |
| Top 1,000 | Unremarkable — selection bias likely explains it |
| Top 1,000–10,000 | Weak — many triplets score comparably or better |
| Top 10,000+ | Fully explained by selection — the "Bible code" explanation wins |

**Alison triplet rank: {rank_combined:,} of {idx:,} (equivalent to ~{equiv_rank_1_3M:,} at 1.3M scale) → {interpretation}**
"""

    with open(os.path.join(OUT_DIR, "TRIPLET_VERDICT.md"), "w") as f:
        f.write(verdict_md)

    print(f"\n  Output files written to {OUT_DIR}/")
    print(f"  - triplet_scan.json")
    print(f"  - triplet_ranks.json")
    print(f"  - TRIPLET_VERDICT.md")
    print(f"\n  VERDICT: Rank {rank_combined:,} → {interpretation}")


if __name__ == "__main__":
    main()
