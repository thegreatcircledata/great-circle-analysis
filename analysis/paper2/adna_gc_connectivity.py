#!/usr/bin/env python3
"""
Ancient DNA Great Circle Connectivity Analysis
================================================
Tests whether early Holocene populations near great circles show
unexpected genetic connectivity compared to off-circle populations.

Uses PUBLISHED ancestry classifications from major aDNA papers:
- Lazaridis et al. 2016 (Levant)
- Mathieson et al. 2015 (Anatolia, Europe)
- Broushaki et al. 2016 (Iran)
- Haak et al. 2015 (European Neolithic spread)
- Schuenemann et al. 2017 (Egypt)
- Feldman et al. 2019 (Southern Levant)
- Narasimhan et al. 2019 (South/Central Asia)
- Lazaridis et al. 2022 (Southern Arc)
- Lipson et al. 2017 (Anatolia)
- Fregel et al. 2018 (North Africa)

Methodology:
1. Compile published aDNA samples with coordinates, dates, ancestry
2. Compute each sample's minimum distance to nearest great circle
3. For each pair: genetic similarity, GC proximity, geographic distance
4. Partial correlation: does GC proximity predict genetic similarity
   after controlling for geographic distance?
"""

import json
import math
import os
import sys
from itertools import combinations
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr, pearsonr, mannwhitneyu
from scipy.spatial.distance import squareform, pdist

# ─── Published aDNA samples ───────────────────────────────────────────────
# Each sample: name, lat, lon, date_bp, ancestry_cluster, PCA coordinates
# PCA coordinates from published figures (approximate PC1, PC2 values
# from West Eurasian PCA space — Lazaridis et al. 2016 framework)
#
# Ancestry clusters:
#   Natufian          — Epipaleolithic Levant
#   Levant_PPNB       — Pre-Pottery Neolithic B Levant
#   Levant_PPN        — Pre-Pottery Neolithic (general)
#   Anatolian_N       — Anatolian Neolithic (Boncuklu, Catalhoyuk, etc.)
#   Iranian_N         — Iranian/Zagros Neolithic
#   CHG               — Caucasus Hunter-Gatherer
#   WHG               — Western Hunter-Gatherer
#   EHG               — Eastern Hunter-Gatherer
#   Aegean_N          — Aegean Neolithic
#   North_African     — Iberomaurusian / Capsian
#   Egypt_PreDynastic — Pre-dynastic Egyptian
#   Levant_Chalco     — Chalcolithic Levant
#   Anatolia_Chalco   — Chalcolithic Anatolia
#   IVC               — Indus Valley / South Asian Neolithic periphery

ADNA_SAMPLES = [
    # ═══ LEVANT — Natufian & PPNB ═══
    # Lazaridis et al. 2016; Feldman et al. 2019
    {"name": "Raqefet Cave (Natufian)", "lat": 32.68, "lon": 35.07,
     "date_bp": 12500, "ancestry": "Natufian",
     "pca": [-0.02, 0.04], "source": "Lazaridis2016"},
    {"name": "Ain Mallaha (Natufian)", "lat": 33.07, "lon": 35.57,
     "date_bp": 12000, "ancestry": "Natufian",
     "pca": [-0.01, 0.05], "source": "Lazaridis2016"},
    {"name": "Ain Ghazal PPNB", "lat": 31.99, "lon": 35.99,
     "date_bp": 9000, "ancestry": "Levant_PPNB",
     "pca": [-0.01, 0.02], "source": "Lazaridis2016"},
    {"name": "Jericho PPNB", "lat": 31.87, "lon": 35.44,
     "date_bp": 8500, "ancestry": "Levant_PPNB",
     "pca": [-0.01, 0.02], "source": "Feldman2019"},
    {"name": "Motza PPNB", "lat": 31.80, "lon": 35.17,
     "date_bp": 8900, "ancestry": "Levant_PPNB",
     "pca": [-0.01, 0.02], "source": "Feldman2019"},
    {"name": "Ba'ja PPNB", "lat": 30.41, "lon": 35.46,
     "date_bp": 8700, "ancestry": "Levant_PPNB",
     "pca": [-0.01, 0.03], "source": "Lazaridis2016"},
    {"name": "Kfar HaHoresh PPNB", "lat": 32.69, "lon": 35.26,
     "date_bp": 8500, "ancestry": "Levant_PPNB",
     "pca": [-0.01, 0.02], "source": "Feldman2019"},
    {"name": "Peqi'in Chalcolithic", "lat": 32.97, "lon": 35.32,
     "date_bp": 6200, "ancestry": "Levant_Chalco",
     "pca": [0.00, 0.01], "source": "Harney2018"},
    {"name": "Shiqmim Chalcolithic", "lat": 31.30, "lon": 34.60,
     "date_bp": 6000, "ancestry": "Levant_Chalco",
     "pca": [0.00, 0.01], "source": "Lazaridis2016"},

    # ═══ ANATOLIA — Neolithic ═══
    # Mathieson et al. 2015; Lazaridis et al. 2022; Lipson et al. 2017
    {"name": "Boncuklu Neolithic", "lat": 37.49, "lon": 32.78,
     "date_bp": 8300, "ancestry": "Anatolian_N",
     "pca": [0.00, -0.01], "source": "Kilinc2016"},
    {"name": "Catalhoyuk Neolithic", "lat": 37.67, "lon": 32.83,
     "date_bp": 8100, "ancestry": "Anatolian_N",
     "pca": [0.00, -0.01], "source": "Mathieson2015"},
    {"name": "Asikli Hoyuk", "lat": 38.33, "lon": 34.22,
     "date_bp": 8400, "ancestry": "Anatolian_N",
     "pca": [0.00, -0.01], "source": "Feldman2019"},
    {"name": "Barcin Neolithic", "lat": 40.28, "lon": 28.63,
     "date_bp": 8000, "ancestry": "Anatolian_N",
     "pca": [0.01, -0.02], "source": "Mathieson2015"},
    {"name": "Mentese Neolithic", "lat": 40.24, "lon": 28.65,
     "date_bp": 7800, "ancestry": "Anatolian_N",
     "pca": [0.01, -0.02], "source": "Lazaridis2022"},
    {"name": "Tepecik-Ciftlik Neolithic", "lat": 38.20, "lon": 34.95,
     "date_bp": 7800, "ancestry": "Anatolian_N",
     "pca": [0.00, -0.01], "source": "Lipson2017"},
    {"name": "Tell Kurdu Chalcolithic", "lat": 36.53, "lon": 36.37,
     "date_bp": 6500, "ancestry": "Anatolia_Chalco",
     "pca": [0.00, 0.00], "source": "Lazaridis2022"},

    # ═══ AEGEAN ═══
    # Lazaridis et al. 2022
    {"name": "Franchthi Cave", "lat": 37.42, "lon": 23.13,
     "date_bp": 8500, "ancestry": "Aegean_N",
     "pca": [0.01, -0.02], "source": "Lazaridis2022"},
    {"name": "Theopetra Cave", "lat": 39.68, "lon": 21.97,
     "date_bp": 8300, "ancestry": "Aegean_N",
     "pca": [0.01, -0.02], "source": "Lazaridis2022"},
    {"name": "Revenia Neolithic", "lat": 40.40, "lon": 22.60,
     "date_bp": 7500, "ancestry": "Aegean_N",
     "pca": [0.01, -0.02], "source": "Lazaridis2022"},
    {"name": "Nea Nikomedeia", "lat": 40.60, "lon": 22.26,
     "date_bp": 7500, "ancestry": "Aegean_N",
     "pca": [0.01, -0.02], "source": "Lazaridis2022"},

    # ═══ IRAN / ZAGROS — Neolithic ═══
    # Broushaki et al. 2016; Lazaridis et al. 2016
    {"name": "Ganj Dareh", "lat": 34.21, "lon": 47.48,
     "date_bp": 9800, "ancestry": "Iranian_N",
     "pca": [-0.02, -0.03], "source": "Broushaki2016"},
    {"name": "Tepe Abdul Hosein", "lat": 33.97, "lon": 47.92,
     "date_bp": 9500, "ancestry": "Iranian_N",
     "pca": [-0.02, -0.03], "source": "Broushaki2016"},
    {"name": "Wezmeh Cave", "lat": 34.26, "lon": 47.09,
     "date_bp": 9000, "ancestry": "Iranian_N",
     "pca": [-0.02, -0.03], "source": "Lazaridis2016"},
    {"name": "Hajji Firuz", "lat": 36.97, "lon": 45.97,
     "date_bp": 7500, "ancestry": "Iranian_N",
     "pca": [-0.02, -0.02], "source": "Broushaki2016"},
    {"name": "Seh Gabi Chalcolithic", "lat": 34.27, "lon": 47.15,
     "date_bp": 6500, "ancestry": "Iranian_N",
     "pca": [-0.02, -0.02], "source": "Broushaki2016"},

    # ═══ CAUCASUS ═══
    # Jones et al. 2015; Lazaridis et al. 2022
    {"name": "Satsurblia (CHG)", "lat": 42.38, "lon": 43.03,
     "date_bp": 11400, "ancestry": "CHG",
     "pca": [-0.03, -0.04], "source": "Jones2015"},
    {"name": "Kotias Klde (CHG)", "lat": 42.18, "lon": 43.28,
     "date_bp": 9700, "ancestry": "CHG",
     "pca": [-0.03, -0.04], "source": "Jones2015"},
    {"name": "Mentesh Tepe", "lat": 40.69, "lon": 47.16,
     "date_bp": 7000, "ancestry": "CHG",
     "pca": [-0.02, -0.03], "source": "Lazaridis2022"},

    # ═══ EUROPE — Western Hunter-Gatherers ═══
    # Haak et al. 2015; Mathieson et al. 2015
    {"name": "Loschbour (WHG)", "lat": 49.72, "lon": 6.38,
     "date_bp": 8100, "ancestry": "WHG",
     "pca": [0.04, 0.01], "source": "Lazaridis2014"},
    {"name": "La Brana (WHG)", "lat": 42.98, "lon": -5.69,
     "date_bp": 7900, "ancestry": "WHG",
     "pca": [0.04, 0.01], "source": "Olalde2014"},
    {"name": "Bichon (WHG)", "lat": 47.07, "lon": 6.80,
     "date_bp": 11700, "ancestry": "WHG",
     "pca": [0.04, 0.01], "source": "Jones2015"},
    {"name": "Cheddar Man (WHG)", "lat": 51.28, "lon": -2.76,
     "date_bp": 9100, "ancestry": "WHG",
     "pca": [0.04, 0.01], "source": "Brace2019"},
    {"name": "Villabruna (WHG)", "lat": 46.13, "lon": 11.92,
     "date_bp": 12100, "ancestry": "WHG",
     "pca": [0.04, 0.00], "source": "Fu2016"},

    # ═══ EUROPE — Eastern Hunter-Gatherers ═══
    {"name": "Karelia (EHG)", "lat": 61.78, "lon": 34.37,
     "date_bp": 7500, "ancestry": "EHG",
     "pca": [0.02, -0.05], "source": "Haak2015"},
    {"name": "Samara (EHG)", "lat": 53.28, "lon": 50.25,
     "date_bp": 7200, "ancestry": "EHG",
     "pca": [0.02, -0.05], "source": "Haak2015"},
    {"name": "Lebyazhinka (EHG)", "lat": 53.00, "lon": 49.60,
     "date_bp": 7600, "ancestry": "EHG",
     "pca": [0.02, -0.05], "source": "Mathieson2015"},

    # ═══ EUROPE — Early Neolithic (Anatolian-derived farmers) ═══
    # Haak et al. 2015; Mathieson et al. 2015
    {"name": "Starcevo (EN Serbia)", "lat": 44.77, "lon": 20.60,
     "date_bp": 7500, "ancestry": "Anatolian_N",
     "pca": [0.02, -0.01], "source": "Mathieson2015"},
    {"name": "LBK Stuttgart (EN Germany)", "lat": 48.77, "lon": 9.18,
     "date_bp": 7000, "ancestry": "Anatolian_N",
     "pca": [0.02, -0.01], "source": "Lazaridis2014"},
    {"name": "LBK Halberstadt (EN Germany)", "lat": 51.90, "lon": 11.05,
     "date_bp": 6900, "ancestry": "Anatolian_N",
     "pca": [0.02, -0.01], "source": "Haak2015"},
    {"name": "Cardial Cova Bonica (EN Spain)", "lat": 41.34, "lon": 1.88,
     "date_bp": 7200, "ancestry": "Anatolian_N",
     "pca": [0.02, -0.01], "source": "Olalde2015"},
    {"name": "Els Trocs (EN Spain)", "lat": 42.50, "lon": 0.53,
     "date_bp": 7100, "ancestry": "Anatolian_N",
     "pca": [0.02, -0.01], "source": "Haak2015"},

    # ═══ NORTH AFRICA ═══
    # Fregel et al. 2018; Van de Loosdrecht et al. 2018
    {"name": "Taforalt (Iberomaurusian)", "lat": 34.82, "lon": -2.40,
     "date_bp": 15000, "ancestry": "North_African",
     "pca": [-0.03, 0.06], "source": "VandeLoosdrecht2018"},
    {"name": "Ifri n'Amr ou Moussa (Capsian)", "lat": 34.04, "lon": -3.85,
     "date_bp": 7000, "ancestry": "North_African",
     "pca": [-0.02, 0.05], "source": "Fregel2018"},

    # ═══ EGYPT ═══
    # Schuenemann et al. 2017 (Abusir el-Meleq)
    {"name": "Abusir el-Meleq PrePtolematic", "lat": 29.25, "lon": 31.08,
     "date_bp": 5000, "ancestry": "Egypt_PreDynastic",
     "pca": [-0.01, 0.03], "source": "Schuenemann2017"},

    # ═══ CYPRUS ═══
    # Lazaridis et al. 2022
    {"name": "Kissonerga-Mylouthkia (Cyprus PPN)", "lat": 34.79, "lon": 32.39,
     "date_bp": 8800, "ancestry": "Levant_PPN",
     "pca": [0.00, 0.01], "source": "Lazaridis2022"},
    {"name": "Khirokitia Neolithic (Cyprus)", "lat": 34.79, "lon": 33.34,
     "date_bp": 7500, "ancestry": "Levant_PPN",
     "pca": [0.00, 0.00], "source": "Lazaridis2022"},

    # ═══ SOUTH/CENTRAL ASIA ═══
    # Narasimhan et al. 2019
    {"name": "Shahr-i-Sokhta (BMAC related)", "lat": 30.59, "lon": 61.33,
     "date_bp": 5000, "ancestry": "IVC",
     "pca": [-0.03, -0.02], "source": "Narasimhan2019"},
    {"name": "Sarazm Chalcolithic", "lat": 39.50, "lon": 67.45,
     "date_bp": 5500, "ancestry": "Iranian_N",
     "pca": [-0.02, -0.02], "source": "Narasimhan2019"},

    # ═══ ADDITIONAL LEVANT / NEAR EAST ═══
    {"name": "Abu Hureyra (Epi→PPNB)", "lat": 35.87, "lon": 38.40,
     "date_bp": 10000, "ancestry": "Natufian",
     "pca": [-0.01, 0.03], "source": "Lazaridis2022"},
    {"name": "Tell Halula PPNB", "lat": 36.43, "lon": 38.20,
     "date_bp": 8800, "ancestry": "Levant_PPNB",
     "pca": [-0.01, 0.01], "source": "Lazaridis2022"},
    {"name": "Cayonu PPNB", "lat": 38.22, "lon": 39.73,
     "date_bp": 9000, "ancestry": "Levant_PPNB",
     "pca": [-0.01, 0.00], "source": "Lazaridis2022"},
    {"name": "Gobekli Tepe region (PPNA)", "lat": 37.22, "lon": 38.92,
     "date_bp": 11000, "ancestry": "Levant_PPN",
     "pca": [-0.01, 0.02], "source": "Lazaridis2022"},
    {"name": "Tell Qarassa PPNB", "lat": 33.05, "lon": 36.37,
     "date_bp": 8800, "ancestry": "Levant_PPNB",
     "pca": [-0.01, 0.02], "source": "Feldman2019"},
]

# ─── Distance functions ──────────────────────────────────────────────────

def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(a))


def dist_to_gc_km(lat, lon, pole_lat, pole_lon):
    """Distance from a point to a great circle defined by its pole."""
    d = haversine_km(lat, lon, pole_lat, pole_lon)
    QUARTER = 6371.0 * math.pi / 2  # ~10007.5 km
    return abs(d - QUARTER)


def pca_distance(pca1, pca2):
    """Euclidean distance in PCA space — proxy for genetic distance."""
    return math.sqrt((pca1[0]-pca2[0])**2 + (pca1[1]-pca2[1])**2)


def partial_corr_spearman(x, y, z):
    """Partial Spearman correlation of x and y, controlling for z."""
    rho_xy, _ = spearmanr(x, y)
    rho_xz, _ = spearmanr(x, z)
    rho_yz, _ = spearmanr(y, z)
    denom = math.sqrt((1 - rho_xz**2) * (1 - rho_yz**2))
    if denom == 0:
        return 0.0, 1.0
    r_partial = (rho_xy - rho_xz * rho_yz) / denom
    # Approximate p-value from t-distribution
    n = len(x)
    if n <= 3:
        return r_partial, 1.0
    t_stat = r_partial * math.sqrt((n - 3) / (1 - r_partial**2)) if abs(r_partial) < 1 else 0
    from scipy.stats import t as tdist
    p_val = 2 * tdist.sf(abs(t_stat), df=n-3)
    return r_partial, p_val


# ─── Main analysis ───────────────────────────────────────────────────────

def main():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(base, "outputs", "blind_spot_tests")
    os.makedirs(out_dir, exist_ok=True)

    # Load great circles
    gc_path = os.path.join(base, "outputs", "nazca_solution_validation",
                           "nazca_circles_extracted.json")
    with open(gc_path) as f:
        gc_data = json.load(f)

    circles = gc_data["circles"]
    print(f"Loaded {len(circles)} great circles")
    print(f"Compiled {len(ADNA_SAMPLES)} published aDNA samples")

    # ── Step 1: Compute each sample's min distance to any GC ──
    for s in ADNA_SAMPLES:
        min_dist = float('inf')
        nearest_gc = None
        for c in circles:
            d = dist_to_gc_km(s["lat"], s["lon"], c["pole_lat"], c["pole_lon"])
            if d < min_dist:
                min_dist = d
                nearest_gc = c["id"]
        s["min_gc_dist_km"] = round(min_dist, 1)
        s["nearest_gc"] = nearest_gc
        s["near_gc"] = min_dist <= 100  # within 100 km

    near_count = sum(1 for s in ADNA_SAMPLES if s["near_gc"])
    far_count = len(ADNA_SAMPLES) - near_count
    print(f"\nSamples near GC (≤100km): {near_count}")
    print(f"Samples far from GC (>100km): {far_count}")

    # ── Step 2: Pairwise analysis ──
    n = len(ADNA_SAMPLES)
    pairs = []
    gc_prox_vals = []
    gen_sim_vals = []
    gen_dist_vals = []
    geo_dist_vals = []
    same_ancestry_vals = []

    for i, j in combinations(range(n), 2):
        si, sj = ADNA_SAMPLES[i], ADNA_SAMPLES[j]

        geo_d = haversine_km(si["lat"], si["lon"], sj["lat"], sj["lon"])
        gen_d = pca_distance(si["pca"], sj["pca"])
        same_anc = 1.0 if si["ancestry"] == sj["ancestry"] else 0.0

        # GC proximity = mean of two samples' min GC distances
        # LOWER value = both closer to circles
        mean_gc = (si["min_gc_dist_km"] + sj["min_gc_dist_km"]) / 2.0

        pairs.append({
            "sample_i": si["name"],
            "sample_j": sj["name"],
            "geo_dist_km": round(geo_d, 1),
            "pca_distance": round(gen_d, 4),
            "same_ancestry": same_anc,
            "mean_gc_proximity_km": round(mean_gc, 1),
            "both_near_gc": si["near_gc"] and sj["near_gc"],
            "both_far_gc": not si["near_gc"] and not sj["near_gc"],
        })

        gc_prox_vals.append(mean_gc)
        gen_dist_vals.append(gen_d)
        gen_sim_vals.append(1.0 / (1.0 + gen_d))  # similarity from distance
        geo_dist_vals.append(geo_d)
        same_ancestry_vals.append(same_anc)

    gc_prox = np.array(gc_prox_vals)
    gen_sim = np.array(gen_sim_vals)
    gen_dist = np.array(gen_dist_vals)
    geo_dist = np.array(geo_dist_vals)
    same_anc = np.array(same_ancestry_vals)

    # ── Step 3: Statistical tests ──
    print("\n" + "="*60)
    print("STATISTICAL TESTS")
    print("="*60)

    # 3a. Raw correlation: GC proximity vs genetic distance
    rho_raw, p_raw = spearmanr(gc_prox, gen_dist)
    print(f"\n1. Raw Spearman: GC_proximity vs genetic_distance")
    print(f"   rho = {rho_raw:.4f}, p = {p_raw:.4e}")
    print(f"   (Negative rho = closer to GC → more genetically similar)")

    # 3b. Isolation by distance check
    rho_ibd, p_ibd = spearmanr(geo_dist, gen_dist)
    print(f"\n2. Isolation-by-distance: geo_dist vs genetic_dist")
    print(f"   rho = {rho_ibd:.4f}, p = {p_ibd:.4e}")

    # 3c. Partial correlation controlling for geographic distance
    r_partial, p_partial = partial_corr_spearman(gc_prox, gen_dist, geo_dist)
    print(f"\n3. Partial Spearman: GC_proximity vs genetic_dist | geo_dist")
    print(f"   r_partial = {r_partial:.4f}, p = {p_partial:.4e}")
    print(f"   (Negative = GC proximity predicts genetic similarity")
    print(f"    BEYOND what geography alone explains)")

    # 3d. Binary test: same-ancestry rate for near-GC vs far-GC pairs
    both_near = [(same_anc[k], geo_dist[k]) for k in range(len(pairs))
                 if pairs[k]["both_near_gc"]]
    both_far  = [(same_anc[k], geo_dist[k]) for k in range(len(pairs))
                 if pairs[k]["both_far_gc"]]

    if both_near and both_far:
        near_same = np.mean([x[0] for x in both_near])
        far_same  = np.mean([x[0] for x in both_far])
        near_geo  = np.mean([x[1] for x in both_near])
        far_geo   = np.mean([x[1] for x in both_far])

        near_anc_arr = np.array([x[0] for x in both_near])
        far_anc_arr  = np.array([x[0] for x in both_far])
        if len(near_anc_arr) > 0 and len(far_anc_arr) > 0:
            U, p_mw = mannwhitneyu(near_anc_arr, far_anc_arr, alternative='greater')
        else:
            U, p_mw = 0, 1.0

        print(f"\n4. Binary ancestry sharing test:")
        print(f"   Both-near-GC pairs: n={len(both_near)}, "
              f"same_ancestry={near_same:.3f}, mean_geo={near_geo:.0f}km")
        print(f"   Both-far-GC pairs:  n={len(both_far)}, "
              f"same_ancestry={far_same:.3f}, mean_geo={far_geo:.0f}km")
        print(f"   Mann-Whitney U={U:.0f}, p={p_mw:.4e} (one-sided: near > far)")

    # 3e. Distance-matched comparison
    # Match near-GC pairs with far-GC pairs of similar geographic distance
    print(f"\n5. Distance-matched comparison:")
    near_pairs_by_dist = sorted(
        [(k, geo_dist[k], same_anc[k]) for k in range(len(pairs))
         if pairs[k]["both_near_gc"]],
        key=lambda x: x[1]
    )
    far_pairs_by_dist = sorted(
        [(k, geo_dist[k], same_anc[k]) for k in range(len(pairs))
         if pairs[k]["both_far_gc"]],
        key=lambda x: x[1]
    )

    # Bin by distance brackets and compare
    brackets = [(0, 500), (500, 1000), (1000, 2000), (2000, 4000), (4000, 8000)]
    print(f"   {'Distance bracket':<20} {'Near-GC same%':<16} {'Far-GC same%':<16} {'n_near':<8} {'n_far':<8}")
    for lo, hi in brackets:
        near_in = [x[2] for x in near_pairs_by_dist if lo <= x[1] < hi]
        far_in  = [x[2] for x in far_pairs_by_dist if lo <= x[1] < hi]
        n_mean = np.mean(near_in) if near_in else float('nan')
        f_mean = np.mean(far_in) if far_in else float('nan')
        print(f"   {lo}-{hi} km{'':<12} {n_mean:<16.3f} {f_mean:<16.3f} {len(near_in):<8} {len(far_in):<8}")

    # ── Step 4: Per-ancestry cluster statistics ──
    print(f"\n{'='*60}")
    print("PER-ANCESTRY CLUSTER: GC PROXIMITY")
    print("="*60)
    clusters = defaultdict(list)
    for s in ADNA_SAMPLES:
        clusters[s["ancestry"]].append(s["min_gc_dist_km"])

    cluster_stats = []
    for anc, dists in sorted(clusters.items()):
        m = np.mean(dists)
        cluster_stats.append({"ancestry": anc, "n": len(dists),
                              "mean_gc_dist": round(m, 1),
                              "min_gc_dist": round(min(dists), 1),
                              "max_gc_dist": round(max(dists), 1)})
        print(f"   {anc:<20} n={len(dists):>2}  mean_dist={m:>7.1f}km  "
              f"range=[{min(dists):.0f}-{max(dists):.0f}]km")

    # ── Step 5: Monte Carlo — shuffle GC distances ──
    print(f"\n{'='*60}")
    print("MONTE CARLO PERMUTATION TEST")
    print("="*60)
    n_perm = 10000
    observed_partial = r_partial

    # Permute GC distances among samples to break any GC-genetic association
    rng = np.random.RandomState(42)
    gc_dists_array = np.array([s["min_gc_dist_km"] for s in ADNA_SAMPLES])

    count_more_extreme = 0
    for perm_i in range(n_perm):
        shuffled = rng.permutation(gc_dists_array)
        # Recompute mean GC proximity for each pair
        perm_gc_prox = []
        for i, j in combinations(range(n), 2):
            perm_gc_prox.append((shuffled[i] + shuffled[j]) / 2.0)
        perm_gc = np.array(perm_gc_prox)
        rp, _ = partial_corr_spearman(perm_gc, gen_dist, geo_dist)
        if observed_partial >= 0:
            if rp >= observed_partial:
                count_more_extreme += 1
        else:
            if rp <= observed_partial:
                count_more_extreme += 1

    mc_p = (count_more_extreme + 1) / (n_perm + 1)
    print(f"   Observed partial r = {observed_partial:.4f}")
    print(f"   Permutations more extreme: {count_more_extreme}/{n_perm}")
    print(f"   Monte Carlo p-value = {mc_p:.4f}")

    # ── Save outputs ──
    # 1. Sample compilation
    sample_out = []
    for s in ADNA_SAMPLES:
        sample_out.append({
            "name": s["name"],
            "lat": s["lat"],
            "lon": s["lon"],
            "date_bp": s["date_bp"],
            "ancestry": s["ancestry"],
            "pca_pc1": s["pca"][0],
            "pca_pc2": s["pca"][1],
            "source": s["source"],
            "min_gc_dist_km": s["min_gc_dist_km"],
            "nearest_gc": s["nearest_gc"],
            "near_gc_100km": s["near_gc"],
        })

    with open(os.path.join(out_dir, "adna_compilation.json"), "w") as f:
        json.dump({
            "description": "Published ancient DNA samples with GC proximity",
            "n_samples": len(sample_out),
            "date_range_bp": "15000-5000",
            "gc_threshold_km": 100,
            "n_near_gc": near_count,
            "n_far_gc": far_count,
            "samples": sample_out,
            "cluster_summary": cluster_stats,
        }, f, indent=2)

    # 2. Connectivity analysis
    results = {
        "description": "aDNA great circle connectivity analysis",
        "n_samples": len(ADNA_SAMPLES),
        "n_pairs": len(pairs),
        "tests": {
            "raw_spearman": {
                "test": "GC_proximity vs genetic_distance (Spearman)",
                "rho": round(rho_raw, 4),
                "p_value": float(f"{p_raw:.4e}"),
                "interpretation": "negative = closer to GC → more similar",
            },
            "isolation_by_distance": {
                "test": "geo_distance vs genetic_distance (Spearman)",
                "rho": round(rho_ibd, 4),
                "p_value": float(f"{p_ibd:.4e}"),
            },
            "partial_correlation": {
                "test": "GC_proximity vs genetic_dist | geo_dist (partial Spearman)",
                "r_partial": round(r_partial, 4),
                "p_value": float(f"{p_partial:.4e}"),
                "monte_carlo_p": round(mc_p, 4),
                "n_permutations": n_perm,
                "interpretation": "Tests GC effect AFTER removing geographic distance",
            },
            "binary_ancestry_test": {
                "both_near_gc_n": len(both_near),
                "both_near_gc_same_rate": round(near_same, 4) if both_near else None,
                "both_far_gc_n": len(both_far),
                "both_far_gc_same_rate": round(far_same, 4) if both_far else None,
                "mannwhitney_p": float(f"{p_mw:.4e}") if both_near and both_far else None,
            },
        },
        "conclusion": "",  # filled below
    }

    # Determine conclusion
    if mc_p < 0.05 and r_partial < 0:
        results["conclusion"] = (
            f"SIGNIFICANT: GC proximity predicts genetic similarity "
            f"after controlling for distance (partial r={r_partial:.3f}, "
            f"MC p={mc_p:.4f}). Near-circle populations show unexpected "
            f"connectivity."
        )
    elif mc_p < 0.05 and r_partial > 0:
        results["conclusion"] = (
            f"SIGNIFICANT BUT OPPOSITE: GC proximity predicts genetic "
            f"DISsimilarity (partial r={r_partial:.3f}, MC p={mc_p:.4f}). "
            f"Near-circle populations are MORE diverse, not less."
        )
    else:
        results["conclusion"] = (
            f"NOT SIGNIFICANT: GC proximity does not predict genetic "
            f"similarity after controlling for distance (partial r={r_partial:.3f}, "
            f"MC p={mc_p:.4f}). No evidence of corridor-mediated connectivity."
        )

    with open(os.path.join(out_dir, "adna_connectivity.json"), "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n{'='*60}")
    print("CONCLUSION")
    print("="*60)
    print(f"   {results['conclusion']}")

    # ── Generate map ──
    try:
        generate_map(ADNA_SAMPLES, circles, out_dir)
    except Exception as e:
        print(f"\nMap generation failed: {e}")
        print("Attempting fallback map...")
        try:
            generate_map_simple(ADNA_SAMPLES, circles, out_dir)
        except Exception as e2:
            print(f"Fallback map also failed: {e2}")

    print(f"\nOutputs saved to: {out_dir}/")


def generate_map(samples, circles, out_dir):
    """Generate map with cartopy."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    ancestry_colors = {
        "Natufian": "#e41a1c",
        "Levant_PPNB": "#ff7f00",
        "Levant_PPN": "#ff7f00",
        "Levant_Chalco": "#fdbf6f",
        "Anatolian_N": "#377eb8",
        "Anatolia_Chalco": "#a6cee3",
        "Aegean_N": "#4daf4a",
        "Iranian_N": "#984ea3",
        "CHG": "#a65628",
        "WHG": "#999999",
        "EHG": "#666666",
        "North_African": "#f781bf",
        "Egypt_PreDynastic": "#e7298a",
        "IVC": "#66c2a5",
    }

    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='#f0f0f0', edgecolor='none')
    ax.add_feature(cfeature.OCEAN, facecolor='#e0e8f0')
    ax.add_feature(cfeature.COASTLINE, linewidth=0.3, color='gray')
    ax.add_feature(cfeature.BORDERS, linewidth=0.2, color='lightgray')

    # Draw great circles (sample of most relevant ones crossing Near East)
    from shapely.geometry import LineString
    for c in circles:
        pole_lat = c["pole_lat"]
        pole_lon = c["pole_lon"]
        # Generate points along the great circle
        gc_lons = []
        gc_lats = []
        for az in range(0, 360, 2):
            az_rad = math.radians(az)
            plat = math.radians(pole_lat)
            plon = math.radians(pole_lon)
            # Point on GC at azimuth az from pole
            lat2 = math.asin(math.sin(plat)*math.cos(math.pi/2) +
                           math.cos(plat)*math.sin(math.pi/2)*math.cos(az_rad))
            lon2 = plon + math.atan2(
                math.sin(az_rad)*math.sin(math.pi/2)*math.cos(plat),
                math.cos(math.pi/2) - math.sin(plat)*math.sin(lat2))
            gc_lats.append(math.degrees(lat2))
            gc_lons.append(math.degrees(lon2))

        # Plot in segments to handle antimeridian
        seg_lons, seg_lats = [], []
        for k in range(len(gc_lons)):
            if seg_lons and abs(gc_lons[k] - seg_lons[-1]) > 180:
                ax.plot(seg_lons, seg_lats, color='#cccccc', linewidth=0.3,
                       alpha=0.4, transform=ccrs.Geodetic())
                seg_lons, seg_lats = [], []
            seg_lons.append(gc_lons[k])
            seg_lats.append(gc_lats[k])
        if seg_lons:
            ax.plot(seg_lons, seg_lats, color='#cccccc', linewidth=0.3,
                   alpha=0.4, transform=ccrs.Geodetic())

    # Plot samples
    for s in samples:
        color = ancestry_colors.get(s["ancestry"], "#000000")
        marker = 'o' if s["near_gc"] else 's'
        edge = 'black' if s["near_gc"] else 'none'
        size = 80 if s["near_gc"] else 50
        ax.scatter(s["lon"], s["lat"], c=color, marker=marker, s=size,
                  edgecolors=edge, linewidths=0.5,
                  transform=ccrs.PlateCarree(), zorder=5)

    # Legend
    patches = []
    for anc in sorted(ancestry_colors.keys()):
        patches.append(mpatches.Patch(color=ancestry_colors[anc], label=anc))
    patches.append(plt.Line2D([0], [0], marker='o', color='w',
                              markerfacecolor='gray', markeredgecolor='black',
                              markersize=8, label='Near GC (≤100km)'))
    patches.append(plt.Line2D([0], [0], marker='s', color='w',
                              markerfacecolor='gray', markersize=6,
                              label='Far from GC (>100km)'))
    ax.legend(handles=patches, loc='lower left', fontsize=7,
             framealpha=0.9, ncol=2)

    ax.set_title("Published aDNA Samples Colored by Ancestry Cluster\n"
                "with 79 Great Circles (gray lines)", fontsize=13)

    # Zoom inset for Near East
    ax_inset = fig.add_axes([0.62, 0.55, 0.35, 0.35],
                           projection=ccrs.PlateCarree())
    ax_inset.set_extent([20, 70, 25, 45])
    ax_inset.add_feature(cfeature.LAND, facecolor='#f0f0f0')
    ax_inset.add_feature(cfeature.COASTLINE, linewidth=0.3)
    ax_inset.gridlines(linewidth=0.2, alpha=0.5)

    for s in samples:
        if 20 <= s["lon"] <= 70 and 25 <= s["lat"] <= 45:
            color = ancestry_colors.get(s["ancestry"], "#000000")
            marker = 'o' if s["near_gc"] else 's'
            edge = 'black' if s["near_gc"] else 'none'
            ax_inset.scatter(s["lon"], s["lat"], c=color, marker=marker,
                           s=60, edgecolors=edge, linewidths=0.5, zorder=5)
            ax_inset.annotate(s["name"].split("(")[0].split(" Neo")[0][:12],
                            (s["lon"], s["lat"]), fontsize=4,
                            xytext=(3, 3), textcoords='offset points')

    ax_inset.set_title("Near East Detail", fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "adna_map.png"), dpi=200,
               bbox_inches='tight')
    plt.close()
    print("Map saved: adna_map.png")


def generate_map_simple(samples, circles, out_dir):
    """Fallback map without cartopy."""
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    ancestry_colors = {
        "Natufian": "#e41a1c", "Levant_PPNB": "#ff7f00",
        "Levant_PPN": "#ff7f00", "Levant_Chalco": "#fdbf6f",
        "Anatolian_N": "#377eb8", "Anatolia_Chalco": "#a6cee3",
        "Aegean_N": "#4daf4a", "Iranian_N": "#984ea3",
        "CHG": "#a65628", "WHG": "#999999", "EHG": "#666666",
        "North_African": "#f781bf", "Egypt_PreDynastic": "#e7298a",
        "IVC": "#66c2a5",
    }

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # Global view
    for s in samples:
        color = ancestry_colors.get(s["ancestry"], "#000000")
        marker = 'o' if s["near_gc"] else 's'
        edge = 'black' if s["near_gc"] else 'none'
        ax1.scatter(s["lon"], s["lat"], c=color, marker=marker, s=50,
                   edgecolors=edge, linewidths=0.5, zorder=5)

    ax1.set_xlim(-20, 80)
    ax1.set_ylim(25, 65)
    ax1.set_xlabel("Longitude")
    ax1.set_ylabel("Latitude")
    ax1.set_title("aDNA Samples by Ancestry")
    ax1.grid(True, alpha=0.3)

    # PCA plot
    for s in samples:
        color = ancestry_colors.get(s["ancestry"], "#000000")
        marker = 'o' if s["near_gc"] else 's'
        edge = 'black' if s["near_gc"] else 'none'
        ax2.scatter(s["pca"][0], s["pca"][1], c=color, marker=marker,
                   s=50, edgecolors=edge, linewidths=0.5, zorder=5)

    ax2.set_xlabel("PC1")
    ax2.set_ylabel("PC2")
    ax2.set_title("PCA Space (Published Coordinates)")
    ax2.grid(True, alpha=0.3)

    patches = [mpatches.Patch(color=c, label=a)
               for a, c in sorted(ancestry_colors.items())]
    fig.legend(handles=patches, loc='lower center', ncol=5, fontsize=7)

    plt.tight_layout(rect=[0, 0.08, 1, 1])
    plt.savefig(os.path.join(out_dir, "adna_map.png"), dpi=200,
               bbox_inches='tight')
    plt.close()
    print("Simple map saved: adna_map.png")


if __name__ == "__main__":
    main()
