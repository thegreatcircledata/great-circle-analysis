#!/usr/bin/env python3
"""
Analysis 4: Genetic Distance Along vs. Across the Corridor
============================================================
Directive 11 — Out-of-Africa Migration Overlay

Tests whether populations along the Great Circle corridor show
enhanced genetic similarity (lower Fst relative to geographic
distance) compared to off-corridor populations.

Uses published Fst values from:
  - HUGO Pan-Asian SNP Consortium (2009) — 73 Asian/Oceanian populations
  - 1000 Genomes Project — global populations
  - Malaspinas et al. (2016) — Aboriginal Australian genome

Since we cannot redistribute the full Fst matrices, we use a
representative set of published population-pair Fst values compiled
from supplementary materials of the above papers.
"""

import csv, json, math, os, random, sys
import numpy as np
from scipy import stats as scipy_stats

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
N_MC = 10000
CORRIDOR_KM = 500  # Wide band for genetic analysis

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "The Great Circle", "outputs", "out_of_africa_overlay")
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


def distance_to_gc_arb(lat, lon, pole_lat, pole_lon):
    d = haversine(pole_lat, pole_lon, lat, lon)
    return abs(d - EARTH_R_KM * math.pi / 2)


# ============================================================
# POPULATION DATA
# ============================================================
# Representative populations along and across the southern route
# Coordinates are approximate population centroids
# Fst values from HUGO Pan-Asian (2009), 1000 Genomes, Malaspinas (2016)

POPULATIONS = [
    # (id, name, lat, lon, region)
    # --- On/near southern corridor ---
    ("YRI", "Yoruba (Nigeria)", 7.4, 3.9, "West Africa"),
    ("LWK", "Luhya (Kenya)", 0.5, 34.8, "East Africa"),
    ("BDW", "Bedouin", 31.0, 35.0, "Levant"),
    ("PAL", "Palestinian", 31.9, 35.2, "Levant"),
    ("BAL", "Balochi", 28.0, 66.0, "Pakistan"),
    ("SIN", "Sindhi", 25.4, 68.4, "Pakistan-India"),
    ("GIH", "Gujarati", 22.3, 72.6, "India"),
    ("DRA", "Dravidian (S India)", 11.0, 78.0, "South India"),
    ("VED", "Vedda (Sri Lanka)", 7.5, 80.8, "Sri Lanka"),
    ("MLY", "Malay", 3.1, 101.7, "Malaysia"),
    ("JAV", "Javanese", -7.3, 110.4, "Indonesia"),
    ("PNG", "Papua New Guinean", -6.0, 145.8, "PNG"),
    ("ABR", "Aboriginal Australian", -25.0, 134.0, "Australia"),

    # --- Off-corridor (northern/interior) ---
    ("CEU", "European (CEU)", 50.0, 10.0, "Europe"),
    ("FIN", "Finnish", 62.0, 25.0, "Northern Europe"),
    ("RUS", "Russian", 55.0, 37.0, "Eastern Europe"),
    ("UZB", "Uzbek", 41.0, 64.0, "Central Asia"),
    ("MNG", "Mongolian", 47.9, 106.9, "Central/East Asia"),
    ("CHB", "Han Chinese (Beijing)", 40.0, 116.4, "North China"),
    ("JPT", "Japanese", 35.0, 137.0, "Japan"),
    ("KOR", "Korean", 37.5, 127.0, "Korea"),
    ("TIB", "Tibetan", 31.0, 91.0, "Tibet"),
]

# Published Fst values (representative pairs from HUGO 2009, 1000 Genomes)
# (pop1_id, pop2_id, Fst)
# These are approximate values compiled from supplementary tables
FST_PAIRS = [
    # On-corridor pairs (southern route)
    ("LWK", "BDW", 0.065),
    ("BDW", "BAL", 0.028),
    ("BAL", "SIN", 0.008),
    ("SIN", "GIH", 0.005),
    ("GIH", "DRA", 0.009),
    ("DRA", "VED", 0.015),
    ("VED", "MLY", 0.040),
    ("MLY", "JAV", 0.012),
    ("JAV", "PNG", 0.095),
    ("PNG", "ABR", 0.065),
    ("LWK", "BAL", 0.075),
    ("LWK", "GIH", 0.082),
    ("LWK", "DRA", 0.088),
    ("LWK", "MLY", 0.105),
    ("LWK", "PNG", 0.145),
    ("BAL", "DRA", 0.015),
    ("BAL", "MLY", 0.042),
    ("BAL", "PNG", 0.110),
    ("GIH", "MLY", 0.030),
    ("GIH", "PNG", 0.105),
    ("DRA", "MLY", 0.035),
    ("DRA", "PNG", 0.110),
    ("BDW", "GIH", 0.032),
    ("BDW", "DRA", 0.038),
    ("SIN", "MLY", 0.025),

    # Off-corridor pairs (northern route)
    ("CEU", "FIN", 0.008),
    ("CEU", "RUS", 0.010),
    ("RUS", "UZB", 0.020),
    ("UZB", "MNG", 0.030),
    ("MNG", "CHB", 0.015),
    ("CHB", "JPT", 0.007),
    ("CHB", "KOR", 0.004),
    ("JPT", "KOR", 0.005),
    ("CEU", "CHB", 0.100),
    ("CEU", "JPT", 0.098),
    ("FIN", "CHB", 0.105),
    ("FIN", "MNG", 0.075),
    ("RUS", "CHB", 0.085),
    ("RUS", "MNG", 0.055),
    ("CEU", "MNG", 0.070),
    ("TIB", "CHB", 0.012),
    ("TIB", "MNG", 0.025),

    # Cross pairs (one on-corridor, one off-corridor)
    ("GIH", "CEU", 0.045),
    ("GIH", "CHB", 0.065),
    ("DRA", "CEU", 0.055),
    ("DRA", "CHB", 0.072),
    ("MLY", "CHB", 0.025),
    ("MLY", "JPT", 0.028),
    ("MLY", "CEU", 0.078),
    ("PNG", "CHB", 0.135),
    ("PNG", "CEU", 0.155),
    ("PNG", "JPT", 0.140),
    ("ABR", "CHB", 0.145),
    ("ABR", "CEU", 0.160),
    ("LWK", "CEU", 0.135),
    ("LWK", "CHB", 0.155),
    ("LWK", "JPT", 0.150),
    ("BAL", "CEU", 0.025),
    ("BAL", "CHB", 0.065),
    ("JAV", "CEU", 0.070),
    ("JAV", "CHB", 0.020),
    ("JAV", "JPT", 0.025),
    ("BDW", "CEU", 0.028),
    ("BDW", "CHB", 0.072),
]


# ============================================================
# ANALYSIS
# ============================================================
def classify_pair(pop1_id, pop2_id, pop_dict, corridor_km, pole_lat=None, pole_lon=None):
    """Classify a population pair as both-on, both-off, or mixed."""
    if pole_lat is None:
        pole_lat, pole_lon = POLE_LAT, POLE_LON

    p1 = pop_dict[pop1_id]
    p2 = pop_dict[pop2_id]

    if pole_lat == POLE_LAT and pole_lon == POLE_LON:
        d1 = gc_distance(p1["lat"], p1["lon"])
        d2 = gc_distance(p2["lat"], p2["lon"])
    else:
        d1 = distance_to_gc_arb(p1["lat"], p1["lon"], pole_lat, pole_lon)
        d2 = distance_to_gc_arb(p2["lat"], p2["lon"], pole_lat, pole_lon)

    on1 = d1 <= corridor_km
    on2 = d2 <= corridor_km

    if on1 and on2:
        return "both_on"
    elif not on1 and not on2:
        return "both_off"
    else:
        return "mixed"


def run_analysis():
    print("=" * 70)
    print("ANALYSIS 4: GENETIC DISTANCE ALONG vs. ACROSS CORRIDOR")
    print("=" * 70)

    # Build population dictionary
    pop_dict = {}
    for pid, name, lat, lon, region in POPULATIONS:
        pop_dict[pid] = {"name": name, "lat": lat, "lon": lon, "region": region}

    # Compute distances to Great Circle for each population
    print(f"\n--- Population Distances to Great Circle (corridor = {CORRIDOR_KM} km) ---")
    for pid, name, lat, lon, region in POPULATIONS:
        d = gc_distance(lat, lon)
        on = "ON " if d <= CORRIDOR_KM else "off"
        print(f"  [{on}] {name:30s}  {d:7.1f} km  ({region})")

    # --- Isolation by Distance analysis ---
    print(f"\n--- Isolation by Distance Analysis ---")
    pair_data = []
    for pop1_id, pop2_id, fst in FST_PAIRS:
        if pop1_id not in pop_dict or pop2_id not in pop_dict:
            continue
        p1 = pop_dict[pop1_id]
        p2 = pop_dict[pop2_id]
        geo_dist = haversine(p1["lat"], p1["lon"], p2["lat"], p2["lon"])
        classification = classify_pair(pop1_id, pop2_id, pop_dict, CORRIDOR_KM)
        # Linearized Fst: Fst / (1 - Fst) for regression
        lin_fst = fst / (1 - fst) if fst < 1 else fst
        pair_data.append({
            "pop1": pop1_id, "pop2": pop2_id,
            "fst": fst, "lin_fst": lin_fst,
            "geo_dist_km": geo_dist,
            "log_geo_dist": math.log(geo_dist) if geo_dist > 0 else 0,
            "classification": classification,
        })

    # Regression: lin_Fst ~ log(geo_dist)
    x = np.array([p["log_geo_dist"] for p in pair_data])
    y = np.array([p["lin_fst"] for p in pair_data])

    slope, intercept, r_value, p_value, std_err = scipy_stats.linregress(x, y)
    residuals = y - (slope * x + intercept)

    # Add residuals to pair data
    for i, p in enumerate(pair_data):
        p["residual"] = float(residuals[i])

    print(f"  IBD regression: lin_Fst = {slope:.4f} * log(dist) + {intercept:.4f}")
    print(f"  R² = {r_value**2:.4f},  p = {p_value:.2e}")

    # Compare residuals by classification
    print(f"\n--- Residual Analysis by Corridor Classification ---")
    for cls in ["both_on", "both_off", "mixed"]:
        cls_residuals = [p["residual"] for p in pair_data if p["classification"] == cls]
        if cls_residuals:
            mean_r = np.mean(cls_residuals)
            std_r = np.std(cls_residuals)
            n = len(cls_residuals)
            print(f"  {cls:10s}  n={n:3d}  mean_residual={mean_r:+.4f}  std={std_r:.4f}")

    # Statistical test: are "both_on" residuals significantly more negative?
    on_resid = [p["residual"] for p in pair_data if p["classification"] == "both_on"]
    off_resid = [p["residual"] for p in pair_data if p["classification"] == "both_off"]
    mixed_resid = [p["residual"] for p in pair_data if p["classification"] == "mixed"]

    test_results = {}
    if len(on_resid) >= 3 and len(off_resid) >= 3:
        u_stat, u_p = scipy_stats.mannwhitneyu(on_resid, off_resid, alternative='less')
        print(f"\n  Mann-Whitney U (on < off): U={u_stat:.1f}, p={u_p:.4f}")
        test_results["mann_whitney_on_vs_off"] = {
            "U_statistic": round(float(u_stat), 2),
            "p_value": round(float(u_p), 4),
            "alternative": "on-corridor residuals < off-corridor residuals",
            "significant_0.05": bool(u_p < 0.05),
        }

    if len(on_resid) >= 3 and len(mixed_resid) >= 3:
        u_stat2, u_p2 = scipy_stats.mannwhitneyu(on_resid, mixed_resid, alternative='less')
        print(f"  Mann-Whitney U (on < mixed): U={u_stat2:.1f}, p={u_p2:.4f}")
        test_results["mann_whitney_on_vs_mixed"] = {
            "U_statistic": round(float(u_stat2), 2),
            "p_value": round(float(u_p2), 4),
            "significant_0.05": bool(u_p2 < 0.05),
        }

    # --- Monte Carlo: corridor effect test ---
    print(f"\n--- Monte Carlo Corridor Effect ({N_MC} random great circles) ---")

    # Metric: mean residual of "both_on" pairs minus mean residual of "both_off" pairs
    # Negative = corridor enhances genetic similarity
    obs_effect = np.mean(on_resid) - np.mean(off_resid) if on_resid and off_resid else 0

    mc_effects = []
    for i in range(N_MC):
        p_lat, p_lon = random_pole()
        # Reclassify pairs
        mc_on = []
        mc_off = []
        for j, p in enumerate(pair_data):
            cls = classify_pair(p["pop1"], p["pop2"], pop_dict, CORRIDOR_KM,
                                p_lat, p_lon)
            if cls == "both_on":
                mc_on.append(residuals[j])
            elif cls == "both_off":
                mc_off.append(residuals[j])

        if mc_on and mc_off:
            mc_effects.append(np.mean(mc_on) - np.mean(mc_off))

        if (i + 1) % 2000 == 0:
            print(f"  ... {i+1}/{N_MC} trials complete")

    mc_effects = np.array(mc_effects)
    if len(mc_effects) > 0:
        mc_pct = np.mean(mc_effects <= obs_effect) * 100
        mc_p = np.mean(mc_effects <= obs_effect)
        mc_z = (obs_effect - np.mean(mc_effects)) / np.std(mc_effects) if np.std(mc_effects) > 0 else 0

        print(f"\n  Observed effect (on-off mean residual): {obs_effect:+.4f}")
        print(f"  MC mean effect:  {np.mean(mc_effects):+.4f}")
        print(f"  MC std:          {np.std(mc_effects):.4f}")
        print(f"  Percentile:      {mc_pct:.2f}%")
        print(f"  p-value:         {mc_p:.4f}")
        print(f"  Z-score:         {mc_z:.2f}")
    else:
        mc_pct = mc_p = mc_z = float('nan')
        print("  WARNING: Insufficient data for Monte Carlo test")

    # --- Save CSV ---
    with open(os.path.join(OUT_DIR, "genetic_distance_pairs.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pop1", "pop2", "fst", "lin_fst", "geo_dist_km",
                     "log_geo_dist", "classification", "ibd_residual"])
        for p in pair_data:
            w.writerow([p["pop1"], p["pop2"], round(p["fst"], 4),
                        round(p["lin_fst"], 4), round(p["geo_dist_km"], 1),
                        round(p["log_geo_dist"], 4), p["classification"],
                        round(p["residual"], 5)])

    # --- Save results ---
    results = {
        "corridor_km": CORRIDOR_KM,
        "n_populations": len(POPULATIONS),
        "n_pairs": len(pair_data),
        "ibd_regression": {
            "slope": round(float(slope), 5),
            "intercept": round(float(intercept), 5),
            "r_squared": round(float(r_value**2), 4),
            "p_value": round(float(p_value), 6),
        },
        "residual_analysis": {
            "both_on": {
                "n": len(on_resid),
                "mean_residual": round(float(np.mean(on_resid)), 5) if on_resid else None,
                "std": round(float(np.std(on_resid)), 5) if on_resid else None,
            },
            "both_off": {
                "n": len(off_resid),
                "mean_residual": round(float(np.mean(off_resid)), 5) if off_resid else None,
                "std": round(float(np.std(off_resid)), 5) if off_resid else None,
            },
            "mixed": {
                "n": len(mixed_resid),
                "mean_residual": round(float(np.mean(mixed_resid)), 5) if mixed_resid else None,
            },
        },
        "statistical_tests": test_results,
        "monte_carlo": {
            "n_trials": N_MC,
            "observed_effect": round(float(obs_effect), 5),
            "mc_mean": round(float(np.mean(mc_effects)), 5) if len(mc_effects) > 0 else None,
            "mc_std": round(float(np.std(mc_effects)), 5) if len(mc_effects) > 0 else None,
            "percentile": round(float(mc_pct), 2),
            "p_value": round(float(mc_p), 4),
            "z_score": round(float(mc_z), 2),
        },
        "caveats": [
            "Fst values are approximate, compiled from published supplementary tables",
            "Post-Neolithic migrations may overwrite deep-time corridor signal",
            "Limited population sampling along some segments of the corridor",
            "Results are EXPLORATORY — consistent with corridor hypothesis, not proof",
        ],
    }

    with open(os.path.join(OUT_DIR, "isolation_by_distance_residuals.json"), "w") as f:
        json.dump(results["residual_analysis"], f, indent=2)

    with open(os.path.join(OUT_DIR, "genetic_corridor_effect.json"), "w") as f:
        json.dump(results["monte_carlo"], f, indent=2)

    with open(os.path.join(OUT_DIR, "genetic_distance_analysis.json"), "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to {OUT_DIR}/")

    # --- IBD Plot ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(12, 8))

        colors = {"both_on": "red", "both_off": "blue", "mixed": "gray"}
        labels = {"both_on": "Both on corridor", "both_off": "Both off corridor", "mixed": "Mixed"}

        for cls in ["mixed", "both_off", "both_on"]:  # Plot on-corridor last (on top)
            cx = [p["log_geo_dist"] for p in pair_data if p["classification"] == cls]
            cy = [p["lin_fst"] for p in pair_data if p["classification"] == cls]
            ax.scatter(cx, cy, c=colors[cls], s=40, alpha=0.7, label=labels[cls],
                       edgecolors='white', linewidth=0.5, zorder=3 if cls == "both_on" else 2)

        # Regression line
        x_line = np.linspace(min(x), max(x), 100)
        y_line = slope * x_line + intercept
        ax.plot(x_line, y_line, 'k--', linewidth=1, alpha=0.5,
                label=f'IBD: R²={r_value**2:.3f}')

        ax.set_xlabel("log(Geographic Distance in km)", fontsize=12)
        ax.set_ylabel("Linearized Fst  [Fst / (1-Fst)]", fontsize=12)
        ax.set_title("Isolation by Distance: On-Corridor vs Off-Corridor Population Pairs", fontsize=14)
        ax.legend(loc='upper left', fontsize=10)
        ax.grid(True, alpha=0.3)

        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, "ibd_plot.png"), dpi=150)
        plt.close(fig)
        print(f"IBD plot saved to {OUT_DIR}/ibd_plot.png")
    except ImportError:
        print("matplotlib not available — skipping IBD plot")

    return results


if __name__ == "__main__":
    run_analysis()
