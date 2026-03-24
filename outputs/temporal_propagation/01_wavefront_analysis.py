#!/usr/bin/env python3
"""
Analysis 1: Space-Time Wavefront Plot & Propagation Tests
==========================================================
Directive 04 — Temporal Propagation & Wavefront Analysis

Compiles all dated monumental sites within 50km of the Great Circle,
computes their position along the circle (degrees from Giza),
and tests for directional propagation via Spearman rank correlation
with 10,000-trial Monte Carlo. Also computes pairwise cluster
synchrony vs. diffusion feasibility.
"""

import csv, json, math, os, sys
import numpy as np
from collections import defaultdict

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
THRESHOLD_KM = 50.0
N_MC = 10000
GIZA_LAT, GIZA_LON = 29.9792, 31.1342  # Reference point: 0 degrees

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "The Great Circle", "outputs", "temporal_propagation")
os.makedirs(OUT_DIR, exist_ok=True)

# Cluster definitions (from merge_all_sites.py)
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
    """Initial bearing from point 1 to point 2, in degrees [0, 360)."""
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
    return (math.degrees(math.atan2(x, y)) + 360) % 360


def position_along_circle(lat, lon):
    """
    Compute position along the Great Circle in degrees (0-360),
    with Giza = 0 degrees. The position is the bearing from the pole
    to the site's closest point on the circle.

    For a great circle defined by its pole, any point on the circle
    can be parameterized by the bearing from the pole.
    """
    # Bearing from pole to the site
    b = bearing(POLE_LAT, POLE_LON, lat, lon)
    # Bearing from pole to Giza
    b_giza = bearing(POLE_LAT, POLE_LON, GIZA_LAT, GIZA_LON)
    # Relative position (Giza = 0)
    return (b - b_giza + 360) % 360


# ============================================================
# DATA LOADING
# ============================================================
def bp_to_bce(age_bp):
    """Convert 14C years BP to approximate calendar BCE (simplified IntCal20)."""
    if age_bp <= 2000:
        return 1950 - age_bp
    elif age_bp <= 5000:
        cal_bp = age_bp * 1.05 - 100
    elif age_bp <= 10000:
        cal_bp = age_bp * 1.08 - 250
    else:
        cal_bp = age_bp * 1.10 - 450
    return 1950 - cal_bp


# Megalithic Portal type -> rough estimated date BCE
TYPE_DATES = {
    "Ancient Temple": (-3500, -500),
    "Ancient Village": (-6000, -1000),
    "Artificial Mound": (-4000, -1000),
    "Chambered Cairn": (-4500, -2500),
    "Dolmen": (-5000, -2000),
    "Geoglyph": (-1500, 500),
    "Henge": (-3500, -2000),
    "Long Barrow": (-4500, -3000),
    "Mastaba": (-3100, -2100),
    "Menhir": (-5000, -1500),
    "Passage Grave": (-4500, -2500),
    "Pyramid": (-2700, -2100),
    "Rock Art": (-10000, -500),
    "Round Barrow(s)": (-2500, -1000),
    "Sculptured Stone": (-800, 800),
    "Standing Stone": (-4000, -1500),
    "Standing Stones": (-4000, -1500),
    "Stone Circle": (-3500, -1500),
    "Stone Row / Alignment": (-3500, -1500),
    "Stone Setting": (-4000, -1500),
    "Temple": (-3500, -300),
    "Tomb": (-5000, -1000),
    "Tower": (-2000, 0),
    "Tumulus": (-3500, -1000),
}

# Monumental keywords (for p3k14c site names)
MONUMENTAL_KW = ['temple', 'pyramid', 'tomb', 'mound', 'monument', 'cairn', 'megalith',
                  'henge', 'huaca', 'ceremonial', 'ritual', 'geoglyph', 'ahu', 'moai',
                  'stela', 'ziggurat', 'palace', 'citadel', 'fortress', 'dolmen',
                  'stone circle', 'standing stone', 'barrow', 'tumulus', 'sanctuary',
                  'necropolis', 'acropolis']


def is_monumental_p3k14c(site_name):
    sn = (site_name or "").lower()
    return any(kw in sn for kw in MONUMENTAL_KW)


def is_monumental_type(type_str):
    t = (type_str or "").lower()
    non_monumental = ['village', 'settlement', 'farmstead', 'villa', 'camp', 'shelter',
                       'findspot', 'mine', 'quarry', 'field system', 'enclosure']
    if any(kw in t for kw in non_monumental):
        return False
    monumental = ['pyramid', 'temple', 'tomb', 'dolmen', 'henge', 'stone circle',
                   'standing stone', 'barrow', 'cairn', 'passage grave', 'menhir',
                   'geoglyph', 'rock art', 'mastaba', 'tumulus', 'ahu', 'moai',
                   'ancient temple', 'tower']
    return any(kw in t for kw in monumental)


def load_all_dated_sites():
    """Load dated sites from p3k14c, Pleiades, Megalithic Portal, and globe_sites_all.json."""
    sites = []

    # 1. p3k14c radiocarbon dates (highest quality)
    p3k14c_path = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")
    p3k_count = 0
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
            # Only include 10500 BCE to 500 CE
            if cal_bce < -10500 or cal_bce > 500:
                continue
            site_name = row.get("SiteName", "")
            sites.append({
                "name": site_name or "Unnamed",
                "lat": lat,
                "lon": lon,
                "date_bce": round(cal_bce),
                "cluster": assign_cluster(lat, lon),
                "source": "p3k14c",
                "is_monumental": is_monumental_p3k14c(site_name),
                "type": "radiocarbon",
            })
            p3k_count += 1
    print(f"p3k14c: {p3k_count} dated sites within {THRESHOLD_KM}km")

    # 2. Pleiades (with minDate)
    pleiades_path = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")
    pleiades_count = 0
    with open(pleiades_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["reprLat"])
                lon = float(row["reprLong"])
                min_date = float(row.get("minDate", ""))
            except (ValueError, KeyError):
                continue
            dist = gc_distance(lat, lon)
            if dist > THRESHOLD_KM:
                continue
            # min_date is in calendar years (negative = BCE)
            date_bce = min_date  # Already in the right format
            if date_bce < -10500 or date_bce > 500:
                continue
            feat_type = row.get("featureTypes", "")
            sites.append({
                "name": row.get("title", "Unnamed"),
                "lat": lat,
                "lon": lon,
                "date_bce": round(date_bce),
                "cluster": assign_cluster(lat, lon),
                "source": "pleiades",
                "is_monumental": is_monumental_type(feat_type),
                "type": feat_type,
            })
            pleiades_count += 1
    print(f"Pleiades: {pleiades_count} dated sites within {THRESHOLD_KM}km")

    # 3. Megalithic Portal sites (type-based dating from globe_sites_all.json)
    portal_path = os.path.join(BASE_DIR, "data", "supplementary", "globe_sites_all.json")
    with open(portal_path) as f:
        all_globe = json.load(f)
    portal_count = 0
    for s in all_globe:
        if s.get("source") != "portal":
            continue
        lat, lon = s["lat"], s["lon"]
        dist = gc_distance(lat, lon)
        if dist > THRESHOLD_KM:
            continue
        site_type = s.get("type", "")
        date_range = TYPE_DATES.get(site_type)
        if not date_range:
            # Try partial matching
            for key, val in TYPE_DATES.items():
                if key.lower() in site_type.lower():
                    date_range = val
                    break
        if not date_range:
            continue
        midpoint = (date_range[0] + date_range[1]) / 2
        sites.append({
            "name": s.get("name", "Unnamed"),
            "lat": lat,
            "lon": lon,
            "date_bce": round(midpoint),
            "cluster": assign_cluster(lat, lon),
            "source": "portal",
            "is_monumental": is_monumental_type(site_type),
            "type": site_type,
        })
        portal_count += 1
    print(f"Megalithic Portal: {portal_count} typed+dated sites within {THRESHOLD_KM}km")

    print(f"\nTotal dated sites: {len(sites)}")
    return sites


# ============================================================
# ANALYSIS
# ============================================================
def spearman_rho(x, y):
    """Compute Spearman rank correlation coefficient."""
    n = len(x)
    if n < 3:
        return 0.0
    rank_x = np.argsort(np.argsort(x)).astype(float)
    rank_y = np.argsort(np.argsort(y)).astype(float)
    d = rank_x - rank_y
    return 1 - 6 * np.sum(d**2) / (n * (n**2 - 1))


def run_wavefront_analysis(sites):
    """Main wavefront analysis with Spearman test and cluster synchrony."""
    # Compute position along circle for each site
    for s in sites:
        s["circle_pos"] = position_along_circle(s["lat"], s["lon"])

    # Filter to monumental sites only for the wavefront test
    monumental = [s for s in sites if s["is_monumental"]]
    print(f"\nMonumental sites for wavefront: {len(monumental)}")

    # --- 1. Spearman rank correlation ---
    positions = np.array([s["circle_pos"] for s in monumental])
    dates = np.array([s["date_bce"] for s in monumental])

    rho_obs = spearman_rho(positions, dates)

    # Monte Carlo: shuffle dates 10,000 times
    rho_null = np.zeros(N_MC)
    for i in range(N_MC):
        shuffled = np.random.permutation(dates)
        rho_null[i] = spearman_rho(positions, shuffled)

    p_value = np.mean(np.abs(rho_null) >= abs(rho_obs))

    if p_value < 0.01:
        interp = f"Significant directional propagation (rho={rho_obs:.3f}, p={p_value:.4f})"
    elif p_value < 0.05:
        interp = f"Marginal evidence of propagation (rho={rho_obs:.3f}, p={p_value:.4f})"
    else:
        interp = f"No significant propagation signal (rho={rho_obs:.3f}, p={p_value:.4f})"

    spearman_result = {
        "rho_observed": round(rho_obs, 4),
        "p_value": round(p_value, 5),
        "n_monumental_sites": len(monumental),
        "n_mc_trials": N_MC,
        "null_rho_mean": round(float(np.mean(rho_null)), 4),
        "null_rho_std": round(float(np.std(rho_null)), 4),
        "interpretation": interp,
    }
    print(f"\nSpearman: rho={rho_obs:.4f}, p={p_value:.5f}")
    print(f"  {interp}")

    # --- 2. Cluster synchrony matrix ---
    clusters = ["Egypt / Levant", "Iran / Persia", "Indus Valley",
                 "Peru / Andes", "Easter Island", "Southeast Asia"]
    cluster_sites = {c: [s for s in monumental if s["cluster"] == c] for c in clusters}

    # Known diffusion rates
    DIFFUSION_NEOLITHIC = 1.0  # km/year
    DIFFUSION_BRONZE_AGE = 5.0  # km/year

    sync_matrix = []
    for i, c1 in enumerate(clusters):
        for j, c2 in enumerate(clusters):
            if j <= i:
                continue
            s1 = cluster_sites[c1]
            s2 = cluster_sites[c2]
            if len(s1) < 2 or len(s2) < 2:
                continue

            dates1 = [s["date_bce"] for s in s1]
            dates2 = [s["date_bce"] for s in s2]
            mean1, mean2 = np.mean(dates1), np.mean(dates2)
            min1, max1 = min(dates1), max(dates1)
            min2, max2 = min(dates2), max(dates2)

            # Temporal overlap
            overlap_start = max(min1, min2)
            overlap_end = min(max1, max2)
            overlap_years = max(0, overlap_end - overlap_start)

            # Geographic distance along circle (use mean positions)
            pos1 = np.mean([s["circle_pos"] for s in s1])
            pos2 = np.mean([s["circle_pos"] for s in s2])
            angular_sep = min(abs(pos2 - pos1), 360 - abs(pos2 - pos1))
            km_sep = angular_sep / 360 * 2 * math.pi * EARTH_R_KM

            # Expected delay under diffusion models
            delay_neolithic = km_sep / DIFFUSION_NEOLITHIC
            delay_bronze = km_sep / DIFFUSION_BRONZE_AGE
            observed_delay = abs(mean1 - mean2)

            sync_matrix.append({
                "cluster_1": c1,
                "cluster_2": c2,
                "n_sites_1": len(s1),
                "n_sites_2": len(s2),
                "mean_date_1": round(mean1),
                "mean_date_2": round(mean2),
                "active_range_1": f"{round(min1)} to {round(max1)} BCE",
                "active_range_2": f"{round(min2)} to {round(max2)} BCE",
                "temporal_overlap_years": round(overlap_years),
                "observed_delay_years": round(observed_delay),
                "km_separation_along_circle": round(km_sep),
                "expected_delay_neolithic_yrs": round(delay_neolithic),
                "expected_delay_bronze_age_yrs": round(delay_bronze),
                "diffusion_feasible_neolithic": bool(observed_delay >= delay_neolithic * 0.5),
                "diffusion_feasible_bronze_age": bool(observed_delay >= delay_bronze * 0.5),
                "interpretation": (
                    "Independent emergence (overlap too large for diffusion)"
                    if observed_delay < delay_bronze * 0.5
                    else "Diffusion plausible at Bronze Age rate"
                    if observed_delay < delay_neolithic * 0.5
                    else "Diffusion plausible at Neolithic rate"
                ),
            })

    # --- 3. Create wavefront plot ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(14, 8))

        cluster_colors = {
            "Egypt / Levant": "#E63946",
            "Iran / Persia": "#F4A261",
            "Indus Valley": "#2A9D8F",
            "Peru / Andes": "#264653",
            "Easter Island": "#E76F51",
            "Southeast Asia": "#606C38",
            "Other": "#999999",
        }

        for cluster, color in cluster_colors.items():
            cs = [s for s in monumental if s["cluster"] == cluster]
            if not cs:
                continue
            x = [s["circle_pos"] for s in cs]
            y = [s["date_bce"] for s in cs]
            ax.scatter(x, y, c=color, s=20, alpha=0.6, label=f"{cluster} (n={len(cs)})",
                      edgecolors='none')

        ax.set_xlabel("Position along Great Circle (degrees from Giza)", fontsize=12)
        ax.set_ylabel("Date (BCE / CE)", fontsize=12)
        ax.set_title(f"Space-Time Wavefront — Monumental Sites on the Great Circle\n"
                     f"Spearman ρ = {rho_obs:.3f}, p = {p_value:.4f} (n={len(monumental)})",
                     fontsize=13)
        ax.legend(loc="upper right", fontsize=9)
        ax.invert_yaxis()  # Oldest at top
        ax.set_xlim(0, 360)
        ax.grid(True, alpha=0.3)

        # Add reference lines for major clusters
        for cluster in clusters:
            cs = cluster_sites.get(cluster, [])
            if len(cs) >= 3:
                mean_pos = np.mean([s["circle_pos"] for s in cs])
                ax.axvline(mean_pos, color=cluster_colors.get(cluster, "gray"),
                          alpha=0.2, linestyle='--')

        plt.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, "wavefront_plot.png"), dpi=150)
        plt.close()
        print("\nSaved wavefront_plot.png")
    except ImportError:
        print("\nMatplotlib not available — skipping plot")

    # --- 4. Save outputs ---
    # CSV of all dated sites
    csv_path = os.path.join(OUT_DIR, "dated_sites_along_circle.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "name", "lat", "lon", "date_bce", "circle_pos", "cluster",
            "source", "is_monumental", "type"
        ])
        writer.writeheader()
        for s in sorted(sites, key=lambda x: x["circle_pos"]):
            writer.writerow({k: s[k] for k in writer.fieldnames})
    print(f"Saved dated_sites_along_circle.csv ({len(sites)} sites)")

    # Spearman result
    with open(os.path.join(OUT_DIR, "spearman_correlation.json"), "w") as f:
        json.dump(spearman_result, f, indent=2)
    print("Saved spearman_correlation.json")

    # Cluster synchrony
    with open(os.path.join(OUT_DIR, "cluster_synchrony_matrix.json"), "w") as f:
        json.dump(sync_matrix, f, indent=2)
    print(f"Saved cluster_synchrony_matrix.json ({len(sync_matrix)} pairs)")

    return spearman_result, sync_matrix


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("WAVEFRONT ANALYSIS — Temporal Propagation Directive")
    print("=" * 60)

    sites = load_all_dated_sites()
    spearman, sync = run_wavefront_analysis(sites)

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)
