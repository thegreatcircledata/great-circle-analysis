#!/usr/bin/env python3
"""
Analysis 2: Deep-Time Archaeological Site Enrichment
=====================================================
Directive 11 — Out-of-Africa Migration Overlay

Compiles all known archaeological sites >40,000 BP along the Out-of-Africa
southern dispersal route, computes distance to the Great Circle, and tests
enrichment via Monte Carlo.  Also performs the temporal layering test —
binning sites by epoch and testing for continuous occupation.

Sites compiled from:
  - Armitage et al. (2011) — Jebel Faya
  - Petraglia et al. (2007) — Jwalapuram
  - O'Connell et al. (2018) — SE Asia & Sahul dates
  - Groucutt et al. (2015) — Arabian/Indian sites
  - Rose et al. (2011) — Dhofar sites
"""

import csv, json, math, os, random, sys
import numpy as np
from collections import defaultdict

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
N_MC = 10000
CORRIDOR_KM = 200  # Wider band for deep-time sparse data

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "The Great Circle", "outputs", "out_of_africa_overlay")
P3K14C = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")
os.makedirs(OUT_DIR, exist_ok=True)

random.seed(42)
np.random.seed(42)

# ============================================================
# GEOMETRY
# ============================================================
def haversine(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))


def gc_distance(lat, lon):
    d = haversine(POLE_LAT, POLE_LON, lat, lon)
    return abs(d - QUARTER_CIRC)


def random_pole():
    z = random.uniform(-1, 1)
    theta = random.uniform(0, 2 * math.pi)
    return math.degrees(math.asin(z)), math.degrees(theta) - 180


def distance_to_gc(lat, lon, pole_lat, pole_lon):
    d = haversine(pole_lat, pole_lon, lat, lon)
    return abs(d - EARTH_R_KM * math.pi / 2)


# ============================================================
# DEEP-TIME ARCHAEOLOGICAL SITE DATABASE
# ============================================================
# Compiled from published literature. Each entry:
#   (lat, lon, name, age_bp, reference)
# age_bp = approximate median calibrated date in years Before Present

DEEP_TIME_SITES = [
    # ===== AFRICA — Departure Zone =====
    (-34.41, 21.22, "Blombos Cave", 100000, "Henshilwood et al. 2011"),
    (-31.58, 29.17, "Klasies River Mouth", 120000, "Deacon & Deacon 1999"),
    (-27.58, 32.00, "Border Cave", 170000, "Backwell et al. 2018"),
    (-29.10, 26.04, "Florisbad", 259000, "Grün et al. 1996"),
    ( 9.58, 41.83, "Dire Dawa / Porc-Epic Cave", 77000, "Pleurdeau 2006"),
    (-4.32, 35.35, "Mumba Rock Shelter", 130000, "Mehlman 1989"),
    (-1.83, 35.33, "Enkapune Ya Muto", 50000, "Ambrose 1998"),

    # ===== ARABIA =====
    (25.21, 56.02, "Jebel Faya", 125000, "Armitage et al. 2011"),
    (17.26, 54.10, "Aybut Al Auwal (Dhofar)", 106000, "Rose et al. 2011"),
    (26.00, 56.10, "Al Wusta", 85000, "Groucutt et al. 2018"),
    (17.10, 53.50, "Jebel Kawr", 45000, "Rose et al. 2011"),
    (21.00, 55.00, "Mudayy", 75000, "Delagnes et al. 2012"),
    (26.27, 56.34, "Jebel Barakah", 130000, "Wahida et al. 2009"),

    # ===== INDIA =====
    (15.37, 78.13, "Jwalapuram", 74000, "Petraglia et al. 2007"),
    (13.22, 77.63, "Attirampakkam", 385000, "Akhilesh et al. 2018"),  # deep time
    (22.94, 77.61, "Bhimbetka", 100000, "Misra 1985"),
    (24.80, 73.72, "16R Dune (Thar Desert)", 96000, "Blinkhorn et al. 2013"),
    (15.50, 76.50, "Jurreru Valley (above Toba)", 77000, "Petraglia et al. 2012"),
    (11.90, 79.80, "Attirampakkam / Pallavaram", 250000, "Pappu et al. 2011"),

    # ===== SE ASIA =====
    (20.22, 103.40, "Tam Pa Ling (Laos)", 63000, "Demeter et al. 2012"),
    ( 3.82, 113.77, "Niah Cave (Borneo)", 40000, "Barker et al. 2007"),
    (-8.53, 120.45, "Liang Bua (Flores)", 50000, "Sutikna et al. 2016"),
    (-4.95, 119.62, "Leang Burung 2 (Sulawesi)", 35000, "Bulbeck et al. 2004"),
    (-5.00, 119.50, "Maros karst (Sulawesi)", 45000, "Aubert et al. 2014"),
    (18.31, 103.75, "That Nang Ing (Laos)", 40000, "Demeter et al. 2015"),
    (14.50, 100.00, "Khorat Plateau (Thailand)", 40000, "Marwick 2009"),

    # ===== SAHUL / AUSTRALIA / PNG =====
    (-12.45, 132.89, "Madjedbebe", 65000, "Clarkson et al. 2017"),
    (-33.72, 142.98, "Lake Mungo", 42000, "Bowler et al. 2003"),
    (-34.67, 150.01, "Jerimalai", 42000, "O'Connor 2007"),
    (-8.53, 125.42, "Laili rock shelter (Timor)", 44000, "O'Connor et al. 2011"),
    (-5.78, 141.10, "Mololo Cave (PNG)", 55000, "Summerhayes et al. 2010"),
    (-3.60, 135.20, "Kebar Valley (PNG)", 26000, "Pasveer 2003"),
    (-2.47, 140.70, "Kosipe (PNG)", 49000, "Summerhayes et al. 2010"),
    (-6.30, 147.00, "Huon Peninsula (PNG)", 42000, "Groube et al. 1986"),
    (-9.42, 147.15, "Caution Bay (PNG)", 5000, "McNiven et al. 2011"),

    # ===== ADDITIONAL DEEP TIME SITES =====
    # Levant (for temporal completeness)
    (32.67, 35.57, "Skhul Cave", 100000, "Grün et al. 2005"),
    (32.67, 35.57, "Qafzeh Cave", 92000, "Valladas et al. 1988"),
    (31.40, 35.47, "Boker Tachtit", 47000, "Marks 1983"),

    # East Africa additional
    (-1.28, 36.25, "Olorgesailie", 320000, "Potts et al. 2018"),
    ( 0.80, 36.42, "Baringo (Kenya)", 240000, "McBrearty 1999"),
]


# ============================================================
# TEMPORAL EPOCH BINS
# ============================================================
EPOCHS = [
    (">50,000 BP",  50000, 999999),
    ("50,000–25,000 BP", 25000, 50000),
    ("25,000–12,000 BP", 12000, 25000),
    ("12,000–8,000 BP", 8000, 12000),
    ("8,000–5,000 BP", 5000, 8000),
    ("5,000–2,000 BP", 2000, 5000),
]


# ============================================================
# ANALYSIS
# ============================================================
def run_analysis():
    print("=" * 70)
    print("ANALYSIS 2: DEEP-TIME ARCHAEOLOGICAL SITE ENRICHMENT")
    print("=" * 70)

    # --- Compute distances ---
    print(f"\n--- Site Distances to Great Circle (corridor = {CORRIDOR_KM} km) ---")
    sites_with_dist = []
    for lat, lon, name, age_bp, ref in DEEP_TIME_SITES:
        d = gc_distance(lat, lon)
        on_corridor = d <= CORRIDOR_KM
        sites_with_dist.append((lat, lon, name, age_bp, ref, d, on_corridor))
        marker = "***" if on_corridor else "   "
        print(f"  {marker} {name:40s}  {age_bp:>7,d} BP  {d:8.1f} km  {ref}")

    n_total = len(sites_with_dist)
    n_on = sum(1 for *_, on in sites_with_dist if on)
    print(f"\n  Total sites: {n_total}")
    print(f"  On corridor (≤{CORRIDOR_KM} km): {n_on}  ({100*n_on/n_total:.1f}%)")
    print(f"  Off corridor: {n_total - n_on}")

    # --- Save CSV ---
    with open(os.path.join(OUT_DIR, "deep_time_sites.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lat", "lon", "name", "age_bp", "reference",
                     "distance_to_gc_km", "on_corridor"])
        for lat, lon, name, age_bp, ref, d, on in sites_with_dist:
            w.writerow([lat, lon, name, age_bp, ref, round(d, 2), on])

    # --- Monte Carlo enrichment test ---
    print(f"\n--- Monte Carlo Enrichment ({N_MC} random great circles) ---")
    mc_counts = []
    mc_fractions = []
    for i in range(N_MC):
        p_lat, p_lon = random_pole()
        count = sum(1 for lat, lon, *_ in DEEP_TIME_SITES
                    if distance_to_gc(lat, lon, p_lat, p_lon) <= CORRIDOR_KM)
        mc_counts.append(count)
        mc_fractions.append(count / n_total)
        if (i + 1) % 2000 == 0:
            print(f"  ... {i+1}/{N_MC} trials complete")

    mc_counts = np.array(mc_counts)
    mc_mean = np.mean(mc_counts)
    mc_std = np.std(mc_counts)
    z_score = (n_on - mc_mean) / mc_std if mc_std > 0 else 0
    p_value = np.mean(mc_counts >= n_on)

    print(f"\n  Observed on-corridor: {n_on}")
    print(f"  MC baseline mean:    {mc_mean:.2f}")
    print(f"  MC baseline std:     {mc_std:.2f}")
    print(f"  Z-score:             {z_score:.2f}")
    print(f"  p-value (≥obs):      {p_value:.4f}")

    # --- Enrichment at multiple thresholds ---
    thresholds = [100, 200, 300, 500]
    threshold_results = {}
    for thresh in thresholds:
        obs = sum(1 for *_, d, _ in sites_with_dist if d <= thresh)
        mc_t = []
        for i in range(N_MC):
            p_lat, p_lon = random_pole()
            c = sum(1 for lat, lon, *_ in DEEP_TIME_SITES
                    if distance_to_gc(lat, lon, p_lat, p_lon) <= thresh)
            mc_t.append(c)
        mc_t = np.array(mc_t)
        t_z = (obs - np.mean(mc_t)) / np.std(mc_t) if np.std(mc_t) > 0 else 0
        t_p = np.mean(mc_t >= obs)
        threshold_results[thresh] = {
            "observed": int(obs),
            "mc_mean": round(float(np.mean(mc_t)), 2),
            "mc_std": round(float(np.std(mc_t)), 2),
            "z_score": round(float(t_z), 2),
            "p_value": round(float(t_p), 4),
        }
        print(f"  {thresh:4d} km:  obs={obs:2d}  mc_mean={np.mean(mc_t):.1f}  "
              f"Z={t_z:.2f}  p={t_p:.4f}")

    # --- Temporal Layering Test ---
    print("\n--- Temporal Layering Test ---")
    epoch_results = {}
    for epoch_name, lo, hi in EPOCHS:
        epoch_sites = [(lat, lon, name, age, ref, d, on)
                       for lat, lon, name, age, ref, d, on in sites_with_dist
                       if lo <= age < hi]
        n_epoch = len(epoch_sites)
        n_epoch_on = sum(1 for *_, on in epoch_sites if on)

        if n_epoch == 0:
            epoch_results[epoch_name] = {"n_sites": 0, "n_on_corridor": 0,
                                          "fraction": 0, "note": "no sites in epoch"}
            print(f"  {epoch_name:25s}  n={n_epoch:3d}  on_corridor=0  (no data)")
            continue

        frac = n_epoch_on / n_epoch if n_epoch > 0 else 0

        # MC for this epoch's sites
        mc_e = []
        epoch_coords = [(lat, lon) for lat, lon, *_ in epoch_sites]
        for _ in range(N_MC):
            p_lat, p_lon = random_pole()
            c = sum(1 for lat, lon in epoch_coords
                    if distance_to_gc(lat, lon, p_lat, p_lon) <= CORRIDOR_KM)
            mc_e.append(c)
        mc_e = np.array(mc_e)
        e_z = (n_epoch_on - np.mean(mc_e)) / np.std(mc_e) if np.std(mc_e) > 0 else 0
        e_p = np.mean(mc_e >= n_epoch_on)

        epoch_results[epoch_name] = {
            "n_sites": n_epoch,
            "n_on_corridor": n_epoch_on,
            "fraction_on_corridor": round(frac, 3),
            "mc_mean": round(float(np.mean(mc_e)), 2),
            "mc_std": round(float(np.std(mc_e)), 2),
            "z_score": round(float(e_z), 2),
            "p_value": round(float(e_p), 4),
            "sites": [name for _, _, name, *_ in epoch_sites],
        }
        print(f"  {epoch_name:25s}  n={n_epoch:3d}  on={n_epoch_on:2d}  "
              f"frac={frac:.2f}  Z={e_z:.2f}  p={e_p:.4f}")

    # --- Load p3k14c for Holocene epochs ---
    print("\n--- Loading p3k14c for Holocene epoch enrichment ---")
    p3k_holocene = {name: {"n_sites": 0, "n_on_corridor": 0} for name, lo, hi in EPOCHS if hi <= 12000}
    if os.path.exists(P3K14C):
        p3k_count = 0
        p3k_on = 0
        with open(P3K14C, "r", encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f)
            seen = set()
            for row in reader:
                try:
                    lat = float(row.get("Lat") or row.get("lat", ""))
                    lon = float(row.get("Long") or row.get("lon", ""))
                    age = float(row.get("Age") or row.get("age", ""))
                except (ValueError, TypeError):
                    continue
                if abs(lat) < 0.01 and abs(lon) < 0.01:
                    continue
                site_key = (round(lat, 2), round(lon, 2))
                if site_key in seen:
                    continue
                seen.add(site_key)

                # Only Holocene epochs
                for epoch_name, lo, hi in EPOCHS:
                    if lo <= age < hi and hi <= 12000:
                        if epoch_name not in p3k_holocene:
                            p3k_holocene[epoch_name] = {"n_sites": 0, "n_on_corridor": 0}
                        d = gc_distance(lat, lon)
                        p3k_holocene[epoch_name]["n_sites"] += 1
                        if d <= CORRIDOR_KM:
                            p3k_holocene[epoch_name]["n_on_corridor"] += 1
                        p3k_count += 1
                        break

        print(f"  Loaded {p3k_count} unique p3k14c sites for Holocene epochs")
        for epoch_name in p3k_holocene:
            info = p3k_holocene[epoch_name]
            if info["n_sites"] > 0:
                frac = info["n_on_corridor"] / info["n_sites"]
                print(f"  {epoch_name:25s}  n={info['n_sites']:5d}  "
                      f"on={info['n_on_corridor']:4d}  frac={frac:.3f}")
                # Update epoch results with p3k14c data
                if epoch_name in epoch_results:
                    epoch_results[epoch_name]["p3k14c_n_sites"] = info["n_sites"]
                    epoch_results[epoch_name]["p3k14c_n_on_corridor"] = info["n_on_corridor"]
                    epoch_results[epoch_name]["p3k14c_fraction"] = round(frac, 4)
                else:
                    epoch_results[epoch_name] = {
                        "n_sites": info["n_sites"],
                        "n_on_corridor": info["n_on_corridor"],
                        "fraction_on_corridor": round(frac, 4),
                        "source": "p3k14c",
                    }
    else:
        print(f"  WARNING: p3k14c not found at {P3K14C}")

    # --- Continuity test ---
    print("\n--- Continuity Test ---")
    occupied_epochs = [name for name, _, _ in EPOCHS
                       if epoch_results.get(name, {}).get("n_on_corridor", 0) > 0
                       or epoch_results.get(name, {}).get("p3k14c_n_on_corridor", 0) > 0]
    total_epochs = len(EPOCHS)
    continuity = len(occupied_epochs) / total_epochs if total_epochs > 0 else 0
    print(f"  Epochs with on-corridor sites: {len(occupied_epochs)}/{total_epochs}")
    print(f"  Occupied epochs: {occupied_epochs}")
    print(f"  Continuity score: {continuity:.2f}")

    # --- Save results ---
    results = {
        "corridor_km": CORRIDOR_KM,
        "n_total_sites": n_total,
        "n_on_corridor": n_on,
        "fraction_on_corridor": round(n_on / n_total, 4),
        "mc_mean": round(float(mc_mean), 2),
        "mc_std": round(float(mc_std), 2),
        "z_score": round(float(z_score), 2),
        "p_value": round(float(p_value), 4),
        "n_monte_carlo": N_MC,
        "threshold_analysis": threshold_results,
        "temporal_layering": epoch_results,
        "continuity": {
            "total_epochs": total_epochs,
            "occupied_epochs": len(occupied_epochs),
            "occupied_epoch_names": occupied_epochs,
            "continuity_score": round(continuity, 3),
        },
    }

    with open(os.path.join(OUT_DIR, "temporal_layering.json"), "w") as f:
        json.dump({"temporal_layering": epoch_results}, f, indent=2)

    with open(os.path.join(OUT_DIR, "continuity_test.json"), "w") as f:
        json.dump(results["continuity"], f, indent=2)

    with open(os.path.join(OUT_DIR, "deep_time_enrichment.json"), "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to {OUT_DIR}/")

    # --- Timeline visualization ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(16, 8))

        # Plot all sites by age and distance
        on_ages = [age for _, _, _, age, _, d, on in sites_with_dist if on]
        on_dists = [d for _, _, _, age, _, d, on in sites_with_dist if on]
        off_ages = [age for _, _, _, age, _, d, on in sites_with_dist if not on]
        off_dists = [d for _, _, _, age, _, d, on in sites_with_dist if not on]

        ax.scatter(on_ages, on_dists, c='red', s=60, zorder=5,
                   label=f'On corridor (≤{CORRIDOR_KM} km)', edgecolors='darkred')
        ax.scatter(off_ages, off_dists, c='gray', s=40, zorder=4, alpha=0.5,
                   label=f'Off corridor (>{CORRIDOR_KM} km)', edgecolors='dimgray')

        # Corridor band
        ax.axhline(CORRIDOR_KM, color='red', linestyle='--', alpha=0.5, linewidth=1)

        # Label notable sites
        for lat, lon, name, age, ref, d, on in sites_with_dist:
            if on and age > 40000:
                ax.annotate(name, (age, d), fontsize=7, ha='left',
                            xytext=(5, 5), textcoords='offset points')

        # Epoch shading
        colors = ['#fee0d2', '#fc9272', '#de2d26', '#bdd7e7', '#6baed6', '#2171b5']
        for idx, (name, lo, hi) in enumerate(EPOCHS):
            ax.axvspan(lo, min(hi, 200000), alpha=0.08, color=colors[idx % len(colors)])
            y_pos = CORRIDOR_KM * 2.5 - idx * 30
            ax.text((lo + min(hi, 150000)) / 2, max(on_dists + off_dists) * 0.95 - idx * 50,
                    name, fontsize=7, ha='center', alpha=0.7)

        ax.set_xlabel("Age (years Before Present)", fontsize=12)
        ax.set_ylabel("Distance to Great Circle (km)", fontsize=12)
        ax.set_title("Deep-Time Archaeological Sites: Distance to Great Circle by Age", fontsize=14)
        ax.legend(loc='upper right')
        ax.set_xlim(0, max(age for _, _, _, age, *_ in sites_with_dist) * 1.05)
        ax.invert_xaxis()  # Older on right
        ax.grid(True, alpha=0.3)

        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, "deep_time_timeline.png"), dpi=150)
        plt.close(fig)
        print(f"Timeline saved to {OUT_DIR}/deep_time_timeline.png")
    except ImportError:
        print("matplotlib not available — skipping timeline visualization")

    return results


if __name__ == "__main__":
    run_analysis()
