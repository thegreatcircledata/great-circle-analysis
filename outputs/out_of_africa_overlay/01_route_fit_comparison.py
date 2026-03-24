#!/usr/bin/env python3
"""
Analysis 1: Route Fit Comparison
=================================
Directive 11 — Out-of-Africa Migration Overlay

Digitises the consensus southern coastal dispersal route from published
waypoints, computes fit to the Great Circle, and runs Monte Carlo
significance testing against 10,000 random great circles.  Also computes
fit to the northern dispersal route for comparison.

Sources:
  - Reyes-Centeno et al. (2014) — southern route waypoints
  - Field & Lahr (2006) — GIS least-cost path OIS-4
  - Tassi et al. (2017) — multiple dispersal waves
  - Groucutt et al. (2015) — mapped sites & routes
"""

import csv, json, math, os, random, sys
import numpy as np

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
N_MC = 10000

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
    """Distance from a point to the Great Circle (km)."""
    d = haversine(POLE_LAT, POLE_LON, lat, lon)
    return abs(d - QUARTER_CIRC)


def random_pole():
    """Generate a uniformly random pole on the sphere."""
    z = random.uniform(-1, 1)
    theta = random.uniform(0, 2 * math.pi)
    lat = math.degrees(math.asin(z))
    lon = math.degrees(theta) - 180
    return lat, lon


def distance_to_gc(lat, lon, pole_lat, pole_lon):
    """Distance from point to arbitrary great circle defined by pole."""
    d = haversine(pole_lat, pole_lon, lat, lon)
    qc = EARTH_R_KM * math.pi / 2
    return abs(d - qc)

# ============================================================
# SOUTHERN DISPERSAL ROUTE — CONSENSUS WAYPOINTS
# ============================================================
# Compiled from Reyes-Centeno et al. (2014), Field & Lahr (2006),
# Groucutt et al. (2015), Petraglia et al. (2010), Bulbeck (2007).
#
# The southern coastal route runs from East Africa through the
# Bab el-Mandeb strait, along the Arabian coast, through coastal
# Iran (Makran), along coastal India, through SE Asia to Sahul.
# Waypoints are approximate centroid locations based on published
# maps — spaced at roughly equal intervals along the route.

SOUTHERN_ROUTE = [
    # (lat, lon, label)
    # East Africa — departure zone
    (-3.4,  37.8,  "East Africa — Rift Valley (departure)"),
    ( 7.0,  39.5,  "Horn of Africa — Ethiopian Rift"),
    (11.5,  43.0,  "Djibouti — Bab el-Mandeb approach"),
    (12.8,  45.0,  "Bab el-Mandeb strait crossing"),

    # Arabia — southern coastal route
    (14.5,  49.0,  "Yemen coast — Hadhramaut"),
    (17.0,  54.5,  "Oman — Dhofar (Nubian Complex sites)"),
    (22.5,  59.5,  "Oman — northeastern coast"),
    (25.3,  56.5,  "UAE — Jebel Faya region"),
    (26.5,  54.0,  "Strait of Hormuz"),

    # Iran — Makran coast
    (25.4,  57.8,  "Iran — western Makran coast"),
    (25.2,  60.5,  "Iran — central Makran coast"),
    (25.0,  63.0,  "Pakistan — Balochistan coast"),

    # India — southern coastal route
    (24.5,  67.0,  "Pakistan — Sindh coast / Indus delta"),
    (22.0,  69.0,  "India — Gujarat coast"),
    (19.0,  72.8,  "India — Mumbai / Konkan coast"),
    (15.0,  74.0,  "India — Goa / Karnataka coast"),
    (11.0,  76.0,  "India — Kerala / Western Ghats"),
    ( 8.0,  77.5,  "India — southern tip (Kanyakumari)"),

    # SE Asia — coastal route to Sunda
    ( 7.0,  80.0,  "Sri Lanka"),
    ( 6.0,  95.0,  "Andaman Sea / Nicobar Islands"),
    ( 4.0,  98.5,  "Sumatra — northern coast"),
    ( 2.0, 103.0,  "Malay Peninsula — Singapore Strait"),
    ( 0.5, 109.0,  "Borneo — western coast"),
    (-2.0, 112.0,  "Java Sea / Borneo south"),
    (-5.0, 119.0,  "Sulawesi / Flores"),

    # Sahul — Australia / PNG
    (-8.0, 131.0,  "Timor — Sahul shelf crossing"),
    (-12.4, 131.0, "Australia — Madjedbebe / Arnhem Land"),
    (-6.0, 143.0,  "PNG — southern coast"),
    (-5.5, 145.0,  "PNG — Huon Peninsula / highlands"),
]

# ============================================================
# NORTHERN DISPERSAL ROUTE — COMPARISON
# ============================================================
# The northern route through the Levant, Anatolia, Central Asia.
# From: Tassi et al. (2017), Reyes-Centeno et al. (2014).

NORTHERN_ROUTE = [
    # (lat, lon, label)
    # East Africa — same departure
    (-3.4,  37.8,  "East Africa — Rift Valley (departure)"),
    ( 7.0,  39.5,  "Horn of Africa — Ethiopian Rift"),

    # Nile corridor northward
    (15.6,  32.5,  "Sudan — Nile Valley"),
    (25.0,  32.8,  "Egypt — Upper Nile"),
    (30.0,  31.2,  "Egypt — Nile Delta"),

    # Levant
    (31.5,  34.5,  "Levant — Negev / Sinai"),
    (33.0,  35.5,  "Levant — Lebanon coast"),
    (36.0,  36.0,  "Syria — Orontes Valley"),

    # Anatolia & Caucasus
    (37.5,  38.0,  "SE Anatolia — Gobekli Tepe region"),
    (39.0,  43.0,  "Eastern Anatolia / Lake Van"),
    (41.5,  44.5,  "Caucasus — Georgia"),

    # Central Asia
    (40.0,  52.0,  "Turkmenistan — Caspian coast"),
    (39.0,  63.0,  "Uzbekistan — Amu Darya"),
    (38.0,  68.0,  "Tajikistan"),
    (37.0,  72.0,  "Afghanistan — Hindu Kush"),

    # East & Northeast Asia
    (35.0,  75.0,  "Kashmir / northern Pakistan"),
    (35.0,  82.0,  "Tibet — northern plateau edge"),
    (38.0,  90.0,  "Tarim Basin"),
    (40.0, 100.0,  "Gansu Corridor — China"),
    (38.0, 110.0,  "North China — Yellow River"),
    (35.0, 120.0,  "East China coast"),
    (25.0, 120.0,  "South China — Taiwan Strait"),
]


# ============================================================
# ANALYSIS
# ============================================================
def compute_route_fit(route, pole_lat, pole_lon):
    """Compute mean, median, rms distance of route waypoints to a great circle."""
    dists = [distance_to_gc(lat, lon, pole_lat, pole_lon) for lat, lon, _ in route]
    mean_d = np.mean(dists)
    median_d = np.median(dists)
    rms_d = np.sqrt(np.mean(np.array(dists)**2))
    max_d = max(dists)
    return dists, mean_d, median_d, rms_d, max_d


def run_analysis():
    print("=" * 70)
    print("ANALYSIS 1: ROUTE FIT COMPARISON")
    print("=" * 70)

    # --- Observed fits ---
    print("\n--- Southern Route Fit to Great Circle ---")
    s_dists, s_mean, s_median, s_rms, s_max = compute_route_fit(
        SOUTHERN_ROUTE, POLE_LAT, POLE_LON
    )
    for (lat, lon, label), d in zip(SOUTHERN_ROUTE, s_dists):
        print(f"  {label:50s}  {d:8.1f} km")
    print(f"\n  Mean distance:   {s_mean:.1f} km")
    print(f"  Median distance: {s_median:.1f} km")
    print(f"  RMS distance:    {s_rms:.1f} km")
    print(f"  Max distance:    {s_max:.1f} km")

    print("\n--- Northern Route Fit to Great Circle ---")
    n_dists, n_mean, n_median, n_rms, n_max = compute_route_fit(
        NORTHERN_ROUTE, POLE_LAT, POLE_LON
    )
    for (lat, lon, label), d in zip(NORTHERN_ROUTE, n_dists):
        print(f"  {label:50s}  {d:8.1f} km")
    print(f"\n  Mean distance:   {n_mean:.1f} km")
    print(f"  Median distance: {n_median:.1f} km")
    print(f"  RMS distance:    {n_rms:.1f} km")
    print(f"  Max distance:    {n_max:.1f} km")

    # --- Save waypoints CSV ---
    with open(os.path.join(OUT_DIR, "southern_route_waypoints.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lat", "lon", "label", "distance_to_gc_km"])
        for (lat, lon, label), d in zip(SOUTHERN_ROUTE, s_dists):
            w.writerow([lat, lon, label, round(d, 2)])

    with open(os.path.join(OUT_DIR, "northern_route_waypoints.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lat", "lon", "label", "distance_to_gc_km"])
        for (lat, lon, label), d in zip(NORTHERN_ROUTE, n_dists):
            w.writerow([lat, lon, label, round(d, 2)])

    # --- Monte Carlo significance test ---
    print(f"\n--- Monte Carlo ({N_MC} random great circles) ---")
    mc_south_means = []
    mc_south_rms = []
    mc_north_means = []
    mc_north_rms = []

    for i in range(N_MC):
        p_lat, p_lon = random_pole()
        _, sm, _, sr, _ = compute_route_fit(SOUTHERN_ROUTE, p_lat, p_lon)
        _, nm, _, nr, _ = compute_route_fit(NORTHERN_ROUTE, p_lat, p_lon)
        mc_south_means.append(sm)
        mc_south_rms.append(sr)
        mc_north_means.append(nm)
        mc_north_rms.append(nr)
        if (i + 1) % 2000 == 0:
            print(f"  ... {i+1}/{N_MC} trials complete")

    mc_south_means = np.array(mc_south_means)
    mc_south_rms = np.array(mc_south_rms)
    mc_north_means = np.array(mc_north_means)
    mc_north_rms = np.array(mc_north_rms)

    # Percentile ranks (lower = better fit)
    south_mean_pct = np.mean(mc_south_means <= s_mean) * 100
    south_rms_pct = np.mean(mc_south_rms <= s_rms) * 100
    north_mean_pct = np.mean(mc_north_means <= n_mean) * 100
    north_rms_pct = np.mean(mc_north_rms <= n_rms) * 100

    # p-values
    south_mean_p = np.mean(mc_south_means <= s_mean)
    south_rms_p = np.mean(mc_south_rms <= s_rms)
    north_mean_p = np.mean(mc_north_means <= n_mean)
    north_rms_p = np.mean(mc_north_rms <= n_rms)

    # Z-scores
    south_mean_z = (s_mean - np.mean(mc_south_means)) / np.std(mc_south_means)
    south_rms_z = (s_rms - np.mean(mc_south_rms)) / np.std(mc_south_rms)
    north_mean_z = (n_mean - np.mean(mc_north_means)) / np.std(mc_north_means)
    north_rms_z = (n_rms - np.mean(mc_north_rms)) / np.std(mc_north_rms)

    print(f"\n  SOUTHERN ROUTE:")
    print(f"    Mean distance: {s_mean:.1f} km  (MC mean: {np.mean(mc_south_means):.1f}, "
          f"MC std: {np.std(mc_south_means):.1f})")
    print(f"    Percentile (mean): {south_mean_pct:.2f}%  |  p = {south_mean_p:.4f}  |  Z = {south_mean_z:.2f}")
    print(f"    RMS distance:  {s_rms:.1f} km  (MC mean: {np.mean(mc_south_rms):.1f})")
    print(f"    Percentile (RMS):  {south_rms_pct:.2f}%  |  p = {south_rms_p:.4f}  |  Z = {south_rms_z:.2f}")

    print(f"\n  NORTHERN ROUTE:")
    print(f"    Mean distance: {n_mean:.1f} km  (MC mean: {np.mean(mc_north_means):.1f})")
    print(f"    Percentile (mean): {north_mean_pct:.2f}%  |  p = {north_mean_p:.4f}  |  Z = {north_mean_z:.2f}")
    print(f"    RMS distance:  {n_rms:.1f} km  (MC mean: {np.mean(mc_north_rms):.1f})")
    print(f"    Percentile (RMS):  {north_rms_pct:.2f}%  |  p = {north_rms_p:.4f}  |  Z = {north_rms_z:.2f}")

    # --- Differential test: southern vs northern ---
    # For each random circle, compute the difference (south_mean - north_mean)
    # Negative = fits southern better; positive = fits northern better
    mc_diff = mc_south_means - mc_north_means
    obs_diff = s_mean - n_mean
    diff_pct = np.mean(mc_diff <= obs_diff) * 100
    diff_p = np.mean(mc_diff <= obs_diff)

    print(f"\n  DIFFERENTIAL (South - North mean distance):")
    print(f"    Observed: {obs_diff:.1f} km  (negative = GC fits south better)")
    print(f"    MC mean diff: {np.mean(mc_diff):.1f} km")
    print(f"    Percentile: {diff_pct:.2f}%  |  p = {diff_p:.4f}")

    # --- Per-segment analysis ---
    print("\n--- Per-Segment Analysis (Southern Route) ---")
    segments = {
        "East Africa": [(lat, lon, lab) for lat, lon, lab in SOUTHERN_ROUTE if "Africa" in lab or "Horn" in lab or "Djibouti" in lab],
        "Arabia": [(lat, lon, lab) for lat, lon, lab in SOUTHERN_ROUTE if any(x in lab for x in ["Yemen", "Oman", "UAE", "Hormuz", "Mandeb"])],
        "Iran-Pakistan": [(lat, lon, lab) for lat, lon, lab in SOUTHERN_ROUTE if any(x in lab for x in ["Iran", "Makran", "Balochistan"])],
        "India": [(lat, lon, lab) for lat, lon, lab in SOUTHERN_ROUTE if "India" in lab or "Gujarat" in lab or "Mumbai" in lab],
        "SE Asia": [(lat, lon, lab) for lat, lon, lab in SOUTHERN_ROUTE if any(x in lab for x in ["Sri Lanka", "Andaman", "Sumatra", "Malay", "Borneo", "Java", "Sulawesi"])],
        "Sahul/PNG": [(lat, lon, lab) for lat, lon, lab in SOUTHERN_ROUTE if any(x in lab for x in ["Timor", "Australia", "PNG"])],
    }

    seg_results = {}
    for seg_name, seg_points in segments.items():
        if not seg_points:
            continue
        seg_dists = [gc_distance(lat, lon) for lat, lon, _ in seg_points]
        seg_mean = np.mean(seg_dists)
        seg_results[seg_name] = {
            "n_waypoints": len(seg_points),
            "mean_distance_km": round(float(seg_mean), 1),
            "min_distance_km": round(float(min(seg_dists)), 1),
            "max_distance_km": round(float(max(seg_dists)), 1),
        }
        print(f"  {seg_name:20s}  n={len(seg_points):2d}  mean={seg_mean:7.1f} km  "
              f"min={min(seg_dists):7.1f}  max={max(seg_dists):7.1f}")

    # --- Save results ---
    results = {
        "southern_route": {
            "n_waypoints": len(SOUTHERN_ROUTE),
            "mean_distance_km": round(float(s_mean), 2),
            "median_distance_km": round(float(s_median), 2),
            "rms_distance_km": round(float(s_rms), 2),
            "max_distance_km": round(float(s_max), 2),
            "mc_mean_baseline_km": round(float(np.mean(mc_south_means)), 2),
            "mc_std_km": round(float(np.std(mc_south_means)), 2),
            "percentile_mean": round(float(south_mean_pct), 2),
            "p_value_mean": round(float(south_mean_p), 5),
            "z_score_mean": round(float(south_mean_z), 2),
            "percentile_rms": round(float(south_rms_pct), 2),
            "p_value_rms": round(float(south_rms_p), 5),
            "z_score_rms": round(float(south_rms_z), 2),
            "per_waypoint_distances_km": [round(float(d), 2) for d in s_dists],
        },
        "northern_route": {
            "n_waypoints": len(NORTHERN_ROUTE),
            "mean_distance_km": round(float(n_mean), 2),
            "median_distance_km": round(float(n_median), 2),
            "rms_distance_km": round(float(n_rms), 2),
            "max_distance_km": round(float(n_max), 2),
            "mc_mean_baseline_km": round(float(np.mean(mc_north_means)), 2),
            "mc_std_km": round(float(np.std(mc_north_means)), 2),
            "percentile_mean": round(float(north_mean_pct), 2),
            "p_value_mean": round(float(north_mean_p), 5),
            "z_score_mean": round(float(north_mean_z), 2),
            "percentile_rms": round(float(north_rms_pct), 2),
            "p_value_rms": round(float(north_rms_p), 5),
            "z_score_rms": round(float(north_rms_z), 2),
        },
        "differential": {
            "observed_diff_km": round(float(obs_diff), 2),
            "mc_mean_diff_km": round(float(np.mean(mc_diff)), 2),
            "percentile": round(float(diff_pct), 2),
            "p_value": round(float(diff_p), 5),
        },
        "per_segment": seg_results,
        "n_monte_carlo": N_MC,
    }

    with open(os.path.join(OUT_DIR, "route_fit_comparison.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {OUT_DIR}/route_fit_comparison.json")

    # --- Generate map ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(18, 9))

        # Draw the Great Circle
        gc_lats = []
        gc_lons = []
        pole_lat_r = math.radians(POLE_LAT)
        pole_lon_r = math.radians(POLE_LON)
        for deg in range(361):
            theta = math.radians(deg)
            # Point on great circle at angular distance 90° from pole
            lat = math.asin(math.sin(pole_lat_r) * math.cos(math.pi/2) +
                           math.cos(pole_lat_r) * math.sin(math.pi/2) * math.cos(theta))
            lon = pole_lon_r + math.atan2(
                math.sin(theta) * math.sin(math.pi/2) * math.cos(pole_lat_r),
                math.cos(math.pi/2) - math.sin(pole_lat_r) * math.sin(lat)
            )
            gc_lats.append(math.degrees(lat))
            gc_lons.append(math.degrees(lon))

        # Split at discontinuities for plotting
        for i in range(len(gc_lons)):
            if gc_lons[i] > 180:
                gc_lons[i] -= 360
            if gc_lons[i] < -180:
                gc_lons[i] += 360

        # Plot in segments to avoid wrap-around lines
        seg_lats, seg_lons = [gc_lats[0]], [gc_lons[0]]
        for i in range(1, len(gc_lons)):
            if abs(gc_lons[i] - gc_lons[i-1]) > 180:
                ax.plot(seg_lons, seg_lats, 'b-', linewidth=2, alpha=0.6)
                seg_lats, seg_lons = [], []
            seg_lats.append(gc_lats[i])
            seg_lons.append(gc_lons[i])
        ax.plot(seg_lons, seg_lats, 'b-', linewidth=2, alpha=0.6, label="Great Circle")

        # Southern route
        s_lats = [lat for lat, lon, _ in SOUTHERN_ROUTE]
        s_lons = [lon for lat, lon, _ in SOUTHERN_ROUTE]
        ax.plot(s_lons, s_lats, 'r-o', markersize=5, linewidth=1.5, alpha=0.8, label="Southern Route")

        # Northern route
        n_lats = [lat for lat, lon, _ in NORTHERN_ROUTE]
        n_lons = [lon for lat, lon, _ in NORTHERN_ROUTE]
        ax.plot(n_lons, n_lats, 'g-s', markersize=4, linewidth=1.5, alpha=0.8, label="Northern Route")

        # Simple coastline from Natural Earth if available
        ne_path = os.path.join(BASE_DIR, "data", "natural_earth")
        if os.path.isdir(ne_path):
            # Try to load coastline
            pass  # Skip for now — matplotlib basemap not guaranteed

        ax.set_xlim(-180, 180)
        ax.set_ylim(-60, 75)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title("Great Circle vs Out-of-Africa Dispersal Routes")
        ax.legend(loc="lower left", fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_aspect("equal")

        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, "route_overlay_map.png"), dpi=150)
        plt.close(fig)
        print(f"Map saved to {OUT_DIR}/route_overlay_map.png")
    except ImportError:
        print("matplotlib not available — skipping map generation")

    return results


if __name__ == "__main__":
    run_analysis()
