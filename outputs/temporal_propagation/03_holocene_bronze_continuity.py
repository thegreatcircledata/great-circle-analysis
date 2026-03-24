#!/usr/bin/env python3
"""
Analysis 3: Early Holocene → Bronze Age Transition
====================================================
Directive 04 — Temporal Propagation & Wavefront Analysis

Tests whether Bronze Age monumental sites preferentially sit near
early Holocene occupation sites on the corridor, suggesting continuity,
or whether the two epochs are spatially uncorrelated (rediscovery).
"""

import csv, json, math, os, sys
import numpy as np

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
THRESHOLD_KM = 50.0
N_MC = 10000

# Epoch definitions
EARLY_HOLOCENE = (-10500, -6500)  # 10500–6500 BCE
BRONZE_AGE = (-3000, -2000)        # 3000–2000 BCE

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "The Great Circle", "outputs", "temporal_propagation")
os.makedirs(OUT_DIR, exist_ok=True)


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


def bearing(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
    return (math.degrees(math.atan2(x, y)) + 360) % 360


def position_along_circle(lat, lon):
    """Position along circle in degrees (Giza = 0)."""
    b = bearing(POLE_LAT, POLE_LON, lat, lon)
    b_giza = bearing(POLE_LAT, POLE_LON, 29.9792, 31.1342)
    return (b - b_giza + 360) % 360


def assign_cluster(lat, lon):
    if haversine(-27.1, -109.3, lat, lon) < 200:
        return "Easter Island"
    if 20 <= lat <= 36 and 25 <= lon <= 40:
        return "Egypt / Levant"
    if -25 <= lat <= 0 and -80 <= lon <= -60:
        return "Peru / Andes"
    if 25 <= lat <= 38 and 44 <= lon <= 60:
        return "Iran / Persia"
    if 20 <= lat <= 35 and 60 <= lon <= 80:
        return "Indus Valley"
    if 0 <= lat <= 25 and 95 <= lon <= 115:
        return "Southeast Asia"
    return "Other"


# ============================================================
# DATA
# ============================================================
def bp_to_bce(age_bp):
    if age_bp <= 2000:
        return 1950 - age_bp
    elif age_bp <= 5000:
        cal_bp = age_bp * 1.05 - 100
    elif age_bp <= 10000:
        cal_bp = age_bp * 1.08 - 250
    else:
        cal_bp = age_bp * 1.10 - 450
    return 1950 - cal_bp


MONUMENTAL_KW = ['temple', 'pyramid', 'tomb', 'mound', 'monument', 'cairn', 'megalith',
                  'henge', 'huaca', 'ceremonial', 'ritual', 'geoglyph', 'ahu', 'moai',
                  'stela', 'ziggurat', 'palace', 'citadel', 'fortress', 'dolmen',
                  'stone circle', 'standing stone', 'barrow', 'tumulus', 'sanctuary',
                  'necropolis', 'acropolis', 'mastaba']


def load_p3k14c_by_epoch():
    """Load p3k14c sites within 50km, split into early Holocene and Bronze Age."""
    early_holocene = []
    bronze_age = []

    p3k14c_path = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")
    with open(p3k14c_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["Lat"])
                lon = float(row["Long"])
                age_bp = float(row["Age"])
            except (ValueError, KeyError):
                continue
            if lat == 0 and lon == 0:
                continue
            dist = gc_distance(lat, lon)
            if dist > THRESHOLD_KM:
                continue
            cal_bce = bp_to_bce(age_bp)
            name = row.get("SiteName", "")
            site = {
                "name": name or "Unnamed",
                "lat": lat, "lon": lon,
                "date_bce": round(cal_bce),
                "circle_pos": position_along_circle(lat, lon),
                "cluster": assign_cluster(lat, lon),
                "source": "p3k14c",
            }

            if EARLY_HOLOCENE[0] <= cal_bce <= EARLY_HOLOCENE[1]:
                early_holocene.append(site)
            elif BRONZE_AGE[0] <= cal_bce <= BRONZE_AGE[1]:
                site["is_monumental"] = any(kw in (name or "").lower() for kw in MONUMENTAL_KW)
                bronze_age.append(site)

    print(f"Early Holocene (10500–6500 BCE): {len(early_holocene)} sites within {THRESHOLD_KM}km")
    print(f"Bronze Age (3000–2000 BCE): {len(bronze_age)} sites within {THRESHOLD_KM}km")
    return early_holocene, bronze_age


# ============================================================
# ANALYSIS
# ============================================================
def compute_nearest_distances(target_sites, reference_sites):
    """For each target site, find distance to nearest reference site."""
    if not reference_sites:
        return np.array([])
    ref_lats = np.array([s["lat"] for s in reference_sites])
    ref_lons = np.array([s["lon"] for s in reference_sites])

    distances = []
    for t in target_sites:
        min_d = float("inf")
        for r in reference_sites:
            d = haversine(t["lat"], t["lon"], r["lat"], r["lon"])
            if d < min_d:
                min_d = d
        distances.append(min_d)
    return np.array(distances)


def run_continuity_analysis(early_holocene, bronze_age):
    """Test spatial continuity between early Holocene and Bronze Age sites."""
    if len(early_holocene) < 3 or len(bronze_age) < 3:
        print("ERROR: Too few sites for analysis")
        return None

    # Observed mean distance from Bronze Age to nearest early Holocene
    obs_distances = compute_nearest_distances(bronze_age, early_holocene)
    obs_mean = float(np.mean(obs_distances))
    obs_median = float(np.median(obs_distances))
    print(f"\nObserved mean distance (Bronze→Holocene): {obs_mean:.1f} km")
    print(f"Observed median distance: {obs_median:.1f} km")

    # Monte Carlo: shuffle early Holocene positions along circle
    # Keep them on the circle but randomize their angular positions
    eh_positions = np.array([s["circle_pos"] for s in early_holocene])

    mc_means = np.zeros(N_MC)
    for i in range(N_MC):
        # Randomly shift all early Holocene positions by a uniform offset
        # This preserves internal structure but randomizes position along circle
        offset = np.random.uniform(0, 360)
        shuffled_positions = (eh_positions + offset) % 360

        # Convert shuffled positions back to lat/lon on the circle
        shuffled_eh = []
        for j, pos in enumerate(shuffled_positions):
            # Compute point on circle at this bearing from pole
            giza_bearing = bearing(POLE_LAT, POLE_LON, 29.9792, 31.1342)
            actual_bearing = (pos + giza_bearing) % 360
            # Destination point at quarter-circumference distance from pole
            lat1 = math.radians(POLE_LAT)
            lon1 = math.radians(POLE_LON)
            d = QUARTER_CIRC / EARTH_R_KM  # angular distance in radians
            brng = math.radians(actual_bearing)
            lat2 = math.asin(math.sin(lat1)*math.cos(d) +
                            math.cos(lat1)*math.sin(d)*math.cos(brng))
            lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d)*math.cos(lat1),
                                      math.cos(d) - math.sin(lat1)*math.sin(lat2))
            shuffled_eh.append({
                "lat": math.degrees(lat2),
                "lon": math.degrees(lon2),
            })

        dists = compute_nearest_distances(bronze_age, shuffled_eh)
        mc_means[i] = np.mean(dists)

    # Percentile: what fraction of MC trials have mean distance ≤ observed?
    percentile = float(np.mean(mc_means <= obs_mean))
    p_value_closer = percentile  # Fraction where random is AS CLOSE or closer

    if percentile < 0.05:
        interp = "Bronze Age sites are significantly closer to early Holocene sites than expected — evidence of CONTINUITY"
    elif percentile > 0.95:
        interp = "Bronze Age sites are significantly FARTHER from early Holocene sites than expected — evidence of AVOIDANCE"
    else:
        interp = "No significant spatial correlation between epochs — corridor may have been REDISCOVERED"

    print(f"Monte Carlo percentile: {percentile:.4f}")
    print(f"  {interp}")

    result = {
        "early_holocene_n": len(early_holocene),
        "bronze_age_n": len(bronze_age),
        "observed_mean_distance_km": round(obs_mean, 1),
        "observed_median_distance_km": round(obs_median, 1),
        "monte_carlo": {
            "n_trials": N_MC,
            "null_mean_distance_km": round(float(np.mean(mc_means)), 1),
            "null_std_km": round(float(np.std(mc_means)), 1),
            "percentile": round(percentile, 5),
            "p_value_continuity": round(percentile, 5),
        },
        "interpretation": interp,
        "cluster_breakdown": {},
    }

    # Per-cluster breakdown
    clusters = set(s["cluster"] for s in bronze_age)
    for c in sorted(clusters):
        ba_c = [s for s in bronze_age if s["cluster"] == c]
        eh_c = [s for s in early_holocene if s["cluster"] == c]
        if ba_c and eh_c:
            d = compute_nearest_distances(ba_c, eh_c)
            result["cluster_breakdown"][c] = {
                "bronze_age_n": len(ba_c),
                "early_holocene_n": len(eh_c),
                "mean_distance_km": round(float(np.mean(d)), 1),
                "median_distance_km": round(float(np.median(d)), 1),
            }
        elif ba_c:
            result["cluster_breakdown"][c] = {
                "bronze_age_n": len(ba_c),
                "early_holocene_n": 0,
                "note": "No early Holocene presence — corridor not occupied in this region",
            }

    # --- Plot ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

        # Left: map along circle
        eh_pos = [s["circle_pos"] for s in early_holocene]
        ba_pos = [s["circle_pos"] for s in bronze_age]

        ax1.scatter(eh_pos, [0]*len(eh_pos), c="#1E88E5", s=15, alpha=0.5,
                   label=f"Early Holocene (n={len(early_holocene)})")
        ax1.scatter(ba_pos, [1]*len(ba_pos), c="#E53935", s=15, alpha=0.5,
                   label=f"Bronze Age (n={len(bronze_age)})")
        ax1.set_xlabel("Position along Great Circle (degrees from Giza)")
        ax1.set_yticks([0, 1])
        ax1.set_yticklabels(["Early Holocene\n10500–6500 BCE", "Bronze Age\n3000–2000 BCE"])
        ax1.set_title("Site Positions Along Circle by Epoch")
        ax1.set_xlim(0, 360)
        ax1.legend(loc="upper right")
        ax1.grid(True, alpha=0.3, axis='x')

        # Right: MC null distribution vs observed
        ax2.hist(mc_means, bins=50, alpha=0.7, color="#999999", edgecolor="black", linewidth=0.5)
        ax2.axvline(obs_mean, color="#E53935", linewidth=2,
                    label=f"Observed = {obs_mean:.0f} km")
        ax2.set_xlabel("Mean distance to nearest Holocene site (km)")
        ax2.set_ylabel("Count")
        ax2.set_title(f"Monte Carlo: Bronze Age → Holocene Proximity\n"
                      f"Percentile = {percentile:.3f}")
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, "continuity_map.png"), dpi=150)
        plt.close()
        print("\nSaved continuity_map.png")
    except ImportError:
        print("\nMatplotlib not available — skipping plot")

    # Save
    with open(os.path.join(OUT_DIR, "holocene_bronze_age_continuity.json"), "w") as f:
        json.dump(result, f, indent=2)
    print("Saved holocene_bronze_age_continuity.json")

    return result


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("EARLY HOLOCENE → BRONZE AGE CONTINUITY TEST")
    print("=" * 60)

    early_holocene, bronze_age = load_p3k14c_by_epoch()
    result = run_continuity_analysis(early_holocene, bronze_age)

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)
