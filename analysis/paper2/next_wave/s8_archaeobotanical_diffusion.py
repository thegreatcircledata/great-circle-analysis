#!/usr/bin/env python3
"""
Study 8: Archaeobotanical Diffusion Rate Test
===============================================
If the corridor functioned as a real diffusion channel, crop packages should
have spread faster along the GC bearing than perpendicular to it.

Tests:
  - Old World wheat-barley-lentil package diffusion waypoints
  - New World maize diffusion waypoints

Phases:
  1. Compute spread bearings & angular offsets from GC
  2. Circular statistics (V-test, Rayleigh test)
  3. Speed analysis (along-GC vs perpendicular components)
  4. Null model (1000 random great circles)
  5. New World equivalent (maize)
"""

import sys, os, math, json, warnings
import numpy as np
from scipy import stats as sp_stats
from pathlib import Path

sys.stdout = os.fdopen(sys.stdout.fileno(), "w", buffering=1)
warnings.filterwarnings("ignore", category=FutureWarning)

np.random.seed(42)

# ── Shared Constants ─────────────────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
MC_TRIALS = 1000

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path("/Users/elliotallan/megalith_site_research")
OUT_DIR = BASE / "outputs/next_wave/archaeobotanical_diffusion"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Geometry ─────────────────────────────────────────────────────────────────

def haversine_km(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2)**2
    return 2 * EARTH_R_KM * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def dist_from_gc(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    d = haversine_km(lat, lon, pole_lat, pole_lon)
    qc = EARTH_R_KM * math.pi / 2
    return abs(d - qc)


def forward_azimuth(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dlam = math.radians(lon2 - lon1)
    x = math.sin(dlam) * math.cos(phi2)
    y = (math.cos(phi1) * math.sin(phi2)
         - math.sin(phi1) * math.cos(phi2) * math.cos(dlam))
    return math.degrees(math.atan2(x, y)) % 360


def gc_bearing_at(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    bearing_to_pole = forward_azimuth(lat, lon, pole_lat, pole_lon)
    return (bearing_to_pole + 90) % 360


def angular_offset_0_90(spread_bearing, gc_bearing):
    """Compute angular offset in [0, 90]. 0 = along GC, 90 = perpendicular."""
    diff = abs(spread_bearing - gc_bearing) % 360
    if diff > 180:
        diff = 360 - diff
    if diff > 90:
        diff = 180 - diff
    return diff


# ── Hardcoded Data ───────────────────────────────────────────────────────────

CROP_DIFFUSION_WAYPOINTS = [
    ("Tell Abu Hureyra", 35.87, 38.40, 11000),
    ("Ohalo II", 32.72, 35.57, 21000),
    ("Çayönü", 38.22, 39.72, 8500),
    ("Jericho", 31.87, 35.44, 9000),
    ("Jarmo", 35.55, 44.95, 7000),
    ("Ali Kosh", 32.40, 47.55, 7500),
    ("Çatalhöyük", 37.67, 32.83, 7400),
    ("Knossos", 35.30, 25.16, 7000),
    ("Franchthi Cave", 37.42, 23.13, 7000),
    ("Karanovo", 42.52, 25.55, 6200),
    ("Lepenski Vir", 44.56, 22.02, 6000),
    ("Linearbandkeramik (avg)", 50.0, 10.0, 5500),
    ("Mehrgarh", 29.37, 67.62, 7000),
    ("Sarazm", 39.50, 67.45, 3500),
    ("Fayum A", 29.48, 30.63, 5200),
    ("Merimde", 30.37, 30.85, 4800),
    ("Nabta Playa", 22.51, 30.72, 7000),
]

MAIZE_DIFFUSION = [
    ("Tehuacan Valley", 18.45, -97.39, 7000),
    ("Guilá Naquitz", 16.95, -96.42, 6250),
    ("Oaxaca", 17.06, -96.72, 6000),
    ("Panama (Aguadulce)", 8.24, -80.55, 5000),
    ("Coastal Ecuador", -2.0, -80.7, 4000),
    ("Coastal Peru", -12.0, -76.8, 3000),
    ("SW United States", 33.5, -108.5, 2100),
    ("Eastern Woodlands", 37.0, -86.0, 1800),
]


# ── Phase 1: Compute Spread Bearings ────────────────────────────────────────

def compute_spread_data(waypoints, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    """Sort by date (oldest first), compute pairwise bearings and offsets."""
    sorted_wp = sorted(waypoints, key=lambda x: -x[3])  # oldest = largest BP
    pairs = []
    for i in range(len(sorted_wp) - 1):
        name1, lat1, lon1, date1 = sorted_wp[i]
        name2, lat2, lon2, date2 = sorted_wp[i + 1]
        spread_bearing = forward_azimuth(lat1, lon1, lat2, lon2)
        gc_bear = gc_bearing_at(lat1, lon1, pole_lat, pole_lon)
        offset = angular_offset_0_90(spread_bearing, gc_bear)
        distance = haversine_km(lat1, lon1, lat2, lon2)
        time_diff = date1 - date2  # years (BP difference)
        pairs.append({
            "from": name1, "to": name2,
            "from_lat": lat1, "from_lon": lon1,
            "to_lat": lat2, "to_lon": lon2,
            "date_from_bp": date1, "date_to_bp": date2,
            "spread_bearing": round(spread_bearing, 2),
            "gc_bearing": round(gc_bear, 2),
            "angular_offset_deg": round(offset, 2),
            "distance_km": round(distance, 2),
            "time_years": time_diff,
        })
    return pairs


# ── Phase 2: Circular Statistics ─────────────────────────────────────────────

def v_test(angles_deg, target_deg=0.0):
    """
    V-test for mean direction against a specified target.
    V = n * R_bar * cos(mean_angle - target)
    Under H0 (uniform), 2*V is approx chi-squared(1).
    """
    angles_rad = np.radians(angles_deg)
    target_rad = np.radians(target_deg)
    n = len(angles_rad)
    C = np.sum(np.cos(angles_rad))
    S = np.sum(np.sin(angles_rad))
    R = np.sqrt(C**2 + S**2)
    R_bar = R / n
    mean_angle = np.arctan2(S, C)
    V = n * R_bar * np.cos(mean_angle - target_rad)
    # Under H0, 2V ~ chi2(1) for one-sided test
    p_value = 1 - sp_stats.chi2.cdf(2 * V, df=1) if V > 0 else 1.0
    return {
        "V_statistic": round(float(V), 4),
        "p_value": round(float(p_value), 6),
        "mean_angle_deg": round(float(np.degrees(mean_angle)) % 360, 2),
        "R_bar": round(float(R_bar), 4),
        "n": n,
    }


def rayleigh_test(angles_deg):
    """Rayleigh test for non-uniformity of circular data."""
    angles_rad = np.radians(angles_deg)
    n = len(angles_rad)
    C = np.sum(np.cos(angles_rad))
    S = np.sum(np.sin(angles_rad))
    R = np.sqrt(C**2 + S**2)
    R_bar = R / n
    Z = n * R_bar**2
    # Approximation for p-value
    p_value = np.exp(-Z) * (1 + (2 * Z - Z**2) / (4 * n)
                            - (24 * Z - 132 * Z**2 + 76 * Z**3 - 9 * Z**4) / (288 * n**2))
    p_value = max(0.0, min(1.0, p_value))
    return {
        "Z_statistic": round(float(Z), 4),
        "R_bar": round(float(R_bar), 4),
        "p_value": round(float(p_value), 6),
        "n": n,
    }


def circular_stats(pairs):
    """Run circular statistics on angular offsets from a set of pairs."""
    offsets = np.array([p["angular_offset_deg"] for p in pairs])
    mean_offset = float(np.mean(offsets))
    median_offset = float(np.median(offsets))
    # V-test against 0 (aligned with GC)
    vt_0 = v_test(offsets, target_deg=0.0)
    # V-test against 90 (perpendicular to GC)
    vt_90 = v_test(offsets, target_deg=90.0)
    # Rayleigh test
    ray = rayleigh_test(offsets)
    return {
        "mean_angular_offset_deg": round(mean_offset, 2),
        "median_angular_offset_deg": round(median_offset, 2),
        "std_angular_offset_deg": round(float(np.std(offsets, ddof=1)), 2),
        "v_test_vs_0": vt_0,
        "v_test_vs_90": vt_90,
        "rayleigh_test": ray,
        "offsets": [round(float(x), 2) for x in offsets],
    }


# ── Phase 3: Speed Analysis ─────────────────────────────────────────────────

def speed_analysis(pairs):
    """Decompose spread into along-GC and perpendicular components."""
    results = []
    speeds_along = []
    speeds_perp = []
    for p in pairs:
        offset_rad = math.radians(p["angular_offset_deg"])
        dist = p["distance_km"]
        time_yrs = p["time_years"]
        if time_yrs <= 0:
            continue
        along = dist * math.cos(offset_rad)
        perp = dist * math.sin(offset_rad)
        speed_along = along / time_yrs
        speed_perp = perp / time_yrs
        speeds_along.append(speed_along)
        speeds_perp.append(speed_perp)
        results.append({
            "from": p["from"], "to": p["to"],
            "along_gc_km": round(along, 2),
            "perp_gc_km": round(perp, 2),
            "speed_along_km_yr": round(speed_along, 4),
            "speed_perp_km_yr": round(speed_perp, 4),
            "time_years": time_yrs,
        })

    speeds_along = np.array(speeds_along)
    speeds_perp = np.array(speeds_perp)

    # Wilcoxon signed-rank: speed_along > speed_perp
    if len(speeds_along) >= 6:
        diffs = speeds_along - speeds_perp
        try:
            stat, p_val = sp_stats.wilcoxon(diffs, alternative="greater")
        except ValueError:
            stat, p_val = float("nan"), 1.0
    else:
        stat, p_val = float("nan"), float("nan")

    return {
        "pairs": results,
        "mean_speed_along_km_yr": round(float(np.mean(speeds_along)), 4) if len(speeds_along) else None,
        "mean_speed_perp_km_yr": round(float(np.mean(speeds_perp)), 4) if len(speeds_perp) else None,
        "median_speed_along_km_yr": round(float(np.median(speeds_along)), 4) if len(speeds_along) else None,
        "median_speed_perp_km_yr": round(float(np.median(speeds_perp)), 4) if len(speeds_perp) else None,
        "wilcoxon_stat": round(float(stat), 4) if not np.isnan(stat) else None,
        "wilcoxon_p": round(float(p_val), 6) if not np.isnan(p_val) else None,
        "n_pairs": len(results),
    }


# ── Phase 4: Null Model ─────────────────────────────────────────────────────

def random_pole():
    """Generate a random pole position, uniform on sphere."""
    z = np.random.uniform(-1, 1)
    lat = np.degrees(np.arcsin(z))
    lon = np.random.uniform(-180, 180)
    return lat, lon


def null_model(waypoints, n_trials=MC_TRIALS):
    """Compare Alison GC mean offset against random GCs."""
    # Compute real mean offset
    real_pairs = compute_spread_data(waypoints)
    real_offsets = [p["angular_offset_deg"] for p in real_pairs]
    real_mean = np.mean(real_offsets)

    random_means = []
    for _ in range(n_trials):
        pole_lat, pole_lon = random_pole()
        sorted_wp = sorted(waypoints, key=lambda x: -x[3])
        offsets = []
        for i in range(len(sorted_wp) - 1):
            name1, lat1, lon1, date1 = sorted_wp[i]
            name2, lat2, lon2, date2 = sorted_wp[i + 1]
            spread_bearing = forward_azimuth(lat1, lon1, lat2, lon2)
            gc_bear = gc_bearing_at(lat1, lon1, pole_lat, pole_lon)
            offset = angular_offset_0_90(spread_bearing, gc_bear)
            offsets.append(offset)
        random_means.append(np.mean(offsets))

    random_means = np.array(random_means)
    percentile = float(np.sum(random_means <= real_mean) / len(random_means) * 100)

    return {
        "real_mean_offset_deg": round(float(real_mean), 2),
        "random_mean_offsets_mean": round(float(np.mean(random_means)), 2),
        "random_mean_offsets_std": round(float(np.std(random_means)), 2),
        "random_mean_offsets_5th": round(float(np.percentile(random_means, 5)), 2),
        "random_mean_offsets_95th": round(float(np.percentile(random_means, 95)), 2),
        "percentile_rank": round(percentile, 2),
        "n_trials": n_trials,
        "random_means": [round(float(x), 2) for x in random_means],
    }


# ── Plotting ─────────────────────────────────────────────────────────────────

def save_offset_histogram(offsets, null_means, real_mean, label, filename):
    """Histogram of angular offsets + null distribution comparison."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Left: histogram of individual offsets
        ax = axes[0]
        ax.hist(offsets, bins=np.arange(0, 95, 5), color="#4A90D9", edgecolor="white",
                alpha=0.85)
        ax.axvline(np.mean(offsets), color="#E74C3C", linewidth=2, linestyle="--",
                   label=f"Mean = {np.mean(offsets):.1f}°")
        ax.axvline(45, color="gray", linewidth=1, linestyle=":", label="Random expectation (45°)")
        ax.set_xlabel("Angular Offset from GC (degrees)")
        ax.set_ylabel("Count")
        ax.set_title(f"{label}: Spread Bearing Offsets")
        ax.legend(fontsize=9)
        ax.set_xlim(0, 90)

        # Right: null distribution
        ax = axes[1]
        ax.hist(null_means, bins=40, color="#95A5A6", edgecolor="white", alpha=0.7,
                label="Random GCs")
        ax.axvline(real_mean, color="#E74C3C", linewidth=2.5, linestyle="-",
                   label=f"Alison GC = {real_mean:.1f}°")
        pct = np.sum(np.array(null_means) <= real_mean) / len(null_means) * 100
        ax.set_xlabel("Mean Angular Offset (degrees)")
        ax.set_ylabel("Count")
        ax.set_title(f"Null Distribution (percentile: {pct:.1f}%)")
        ax.legend(fontsize=9)

        plt.tight_layout()
        plt.savefig(OUT_DIR / filename, dpi=150)
        plt.close()
        print(f"  -> saved {filename}")
    except ImportError:
        print("  [matplotlib not available, skipping plot]")


def save_speed_plot(speed_data, label, filename):
    """Paired dot plot of along-GC vs perp-GC speeds."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        pairs = speed_data["pairs"]
        if not pairs:
            return

        names = [f"{p['from'][:12]} -> {p['to'][:12]}" for p in pairs]
        along = [p["speed_along_km_yr"] for p in pairs]
        perp = [p["speed_perp_km_yr"] for p in pairs]
        y = np.arange(len(names))

        fig, ax = plt.subplots(figsize=(10, max(5, len(names) * 0.45)))
        ax.barh(y - 0.15, along, height=0.3, color="#2ECC71", alpha=0.85,
                label="Along GC")
        ax.barh(y + 0.15, perp, height=0.3, color="#E67E22", alpha=0.85,
                label="Perpendicular to GC")
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=8)
        ax.set_xlabel("Speed (km/yr)")
        ax.set_title(f"{label}: Diffusion Speed Components")
        ax.legend(fontsize=9)
        ax.axvline(0, color="gray", linewidth=0.5)
        plt.tight_layout()
        plt.savefig(OUT_DIR / filename, dpi=150)
        plt.close()
        print(f"  -> saved {filename}")
    except ImportError:
        print("  [matplotlib not available, skipping plot]")


def save_map_plot(waypoints, pairs, label, filename):
    """Simple map of waypoints with spread arrows and GC bearing indicators."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        sorted_wp = sorted(waypoints, key=lambda x: -x[3])
        lats = [w[1] for w in sorted_wp]
        lons = [w[2] for w in sorted_wp]
        names = [w[0] for w in sorted_wp]
        dates = [w[3] for w in sorted_wp]

        fig, ax = plt.subplots(figsize=(12, 8))

        # Draw arrows for each consecutive pair
        for p in pairs:
            ax.annotate("", xy=(p["to_lon"], p["to_lat"]),
                        xytext=(p["from_lon"], p["from_lat"]),
                        arrowprops=dict(arrowstyle="->", color="#E74C3C",
                                        linewidth=1.5, alpha=0.6))

        # Plot sites
        scatter = ax.scatter(lons, lats, c=dates, cmap="viridis_r", s=80,
                             edgecolors="black", linewidth=0.5, zorder=5)
        for i, name in enumerate(names):
            ax.annotate(name, (lons[i], lats[i]), fontsize=7,
                        xytext=(5, 5), textcoords="offset points")

        cbar = plt.colorbar(scatter, ax=ax, shrink=0.6)
        cbar.set_label("Date (BP)")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(f"{label}: Crop Diffusion Waypoints")
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(OUT_DIR / filename, dpi=150)
        plt.close()
        print(f"  -> saved {filename}")
    except ImportError:
        print("  [matplotlib not available, skipping plot]")


# ===========================================================================
#  MAIN EXECUTION
# ===========================================================================

print("=" * 70)
print("Study 8: Archaeobotanical Diffusion Rate Test")
print("=" * 70)

# ── Old World Wheat-Barley-Lentil Package ────────────────────────────────────

print("\n" + "=" * 70)
print("OLD WORLD: Wheat-Barley-Lentil Package")
print("=" * 70)

# Phase 1
print("\n--- Phase 1: Compute Spread Bearings ---")
ow_pairs = compute_spread_data(CROP_DIFFUSION_WAYPOINTS)
print(f"  Computed {len(ow_pairs)} consecutive pairs (sorted by date)")
print(f"  {'From':<25} {'To':<25} {'Offset':>8} {'Dist_km':>10} {'Yrs':>6}")
for p in ow_pairs:
    print(f"  {p['from']:<25} {p['to']:<25} {p['angular_offset_deg']:>8.1f}° "
          f"{p['distance_km']:>10.1f} {p['time_years']:>6}")

# Phase 2
print("\n--- Phase 2: Circular Statistics ---")
ow_circ = circular_stats(ow_pairs)
print(f"  Mean angular offset:   {ow_circ['mean_angular_offset_deg']:.2f}°")
print(f"  Median angular offset: {ow_circ['median_angular_offset_deg']:.2f}°")
print(f"  Std:                   {ow_circ['std_angular_offset_deg']:.2f}°")
print(f"  V-test vs 0° (along GC):  V={ow_circ['v_test_vs_0']['V_statistic']:.4f}, "
      f"p={ow_circ['v_test_vs_0']['p_value']:.6f}")
print(f"  V-test vs 90° (perp GC):  V={ow_circ['v_test_vs_90']['V_statistic']:.4f}, "
      f"p={ow_circ['v_test_vs_90']['p_value']:.6f}")
print(f"  Rayleigh test:             Z={ow_circ['rayleigh_test']['Z_statistic']:.4f}, "
      f"p={ow_circ['rayleigh_test']['p_value']:.6f}")

# Phase 3
print("\n--- Phase 3: Speed Analysis ---")
ow_speed = speed_analysis(ow_pairs)
print(f"  Mean speed along GC:   {ow_speed['mean_speed_along_km_yr']} km/yr")
print(f"  Mean speed perp to GC: {ow_speed['mean_speed_perp_km_yr']} km/yr")
print(f"  Median speed along GC: {ow_speed['median_speed_along_km_yr']} km/yr")
print(f"  Median speed perp GC:  {ow_speed['median_speed_perp_km_yr']} km/yr")
wstat = ow_speed['wilcoxon_stat']
wpval = ow_speed['wilcoxon_p']
print(f"  Wilcoxon signed-rank:  stat={wstat}, p={wpval}")

# Phase 4
print("\n--- Phase 4: Null Model (1000 random GCs) ---")
ow_null = null_model(CROP_DIFFUSION_WAYPOINTS)
print(f"  Alison GC mean offset:    {ow_null['real_mean_offset_deg']:.2f}°")
print(f"  Random GCs mean(mean):    {ow_null['random_mean_offsets_mean']:.2f}°")
print(f"  Random GCs std:           {ow_null['random_mean_offsets_std']:.2f}°")
print(f"  Random 5th-95th pctile:   {ow_null['random_mean_offsets_5th']:.2f}° - "
      f"{ow_null['random_mean_offsets_95th']:.2f}°")
print(f"  Alison percentile rank:   {ow_null['percentile_rank']:.2f}%")

# Plots
print("\n--- Plots ---")
save_offset_histogram(ow_circ["offsets"], ow_null["random_means"],
                      ow_null["real_mean_offset_deg"],
                      "Old World", "ow_offset_histogram.png")
save_speed_plot(ow_speed, "Old World", "ow_speed_comparison.png")
save_map_plot(CROP_DIFFUSION_WAYPOINTS, ow_pairs, "Old World", "ow_diffusion_map.png")


# ── New World: Maize Diffusion ───────────────────────────────────────────────

print("\n" + "=" * 70)
print("NEW WORLD: Maize Diffusion")
print("=" * 70)

# Phase 1
print("\n--- Phase 1: Compute Spread Bearings ---")
nw_pairs = compute_spread_data(MAIZE_DIFFUSION)
print(f"  Computed {len(nw_pairs)} consecutive pairs (sorted by date)")
print(f"  {'From':<25} {'To':<25} {'Offset':>8} {'Dist_km':>10} {'Yrs':>6}")
for p in nw_pairs:
    print(f"  {p['from']:<25} {p['to']:<25} {p['angular_offset_deg']:>8.1f}° "
          f"{p['distance_km']:>10.1f} {p['time_years']:>6}")

# Phase 2
print("\n--- Phase 2: Circular Statistics ---")
nw_circ = circular_stats(nw_pairs)
print(f"  Mean angular offset:   {nw_circ['mean_angular_offset_deg']:.2f}°")
print(f"  Median angular offset: {nw_circ['median_angular_offset_deg']:.2f}°")
print(f"  Std:                   {nw_circ['std_angular_offset_deg']:.2f}°")
print(f"  V-test vs 0° (along GC):  V={nw_circ['v_test_vs_0']['V_statistic']:.4f}, "
      f"p={nw_circ['v_test_vs_0']['p_value']:.6f}")
print(f"  V-test vs 90° (perp GC):  V={nw_circ['v_test_vs_90']['V_statistic']:.4f}, "
      f"p={nw_circ['v_test_vs_90']['p_value']:.6f}")
print(f"  Rayleigh test:             Z={nw_circ['rayleigh_test']['Z_statistic']:.4f}, "
      f"p={nw_circ['rayleigh_test']['p_value']:.6f}")

# Phase 3
print("\n--- Phase 3: Speed Analysis ---")
nw_speed = speed_analysis(nw_pairs)
print(f"  Mean speed along GC:   {nw_speed['mean_speed_along_km_yr']} km/yr")
print(f"  Mean speed perp to GC: {nw_speed['mean_speed_perp_km_yr']} km/yr")
print(f"  Median speed along GC: {nw_speed['median_speed_along_km_yr']} km/yr")
print(f"  Median speed perp GC:  {nw_speed['median_speed_perp_km_yr']} km/yr")
wstat_nw = nw_speed['wilcoxon_stat']
wpval_nw = nw_speed['wilcoxon_p']
print(f"  Wilcoxon signed-rank:  stat={wstat_nw}, p={wpval_nw}")

# Phase 4
print("\n--- Phase 4: Null Model (1000 random GCs) ---")
nw_null = null_model(MAIZE_DIFFUSION)
print(f"  Alison GC mean offset:    {nw_null['real_mean_offset_deg']:.2f}°")
print(f"  Random GCs mean(mean):    {nw_null['random_mean_offsets_mean']:.2f}°")
print(f"  Random GCs std:           {nw_null['random_mean_offsets_std']:.2f}°")
print(f"  Random 5th-95th pctile:   {nw_null['random_mean_offsets_5th']:.2f}° - "
      f"{nw_null['random_mean_offsets_95th']:.2f}°")
print(f"  Alison percentile rank:   {nw_null['percentile_rank']:.2f}%")

# Plots
print("\n--- Plots ---")
save_offset_histogram(nw_circ["offsets"], nw_null["random_means"],
                      nw_null["real_mean_offset_deg"],
                      "New World (Maize)", "nw_offset_histogram.png")
save_speed_plot(nw_speed, "New World (Maize)", "nw_speed_comparison.png")
save_map_plot(MAIZE_DIFFUSION, nw_pairs, "New World (Maize)", "nw_diffusion_map.png")


# ===========================================================================
#  VERDICT
# ===========================================================================
print("\n" + "=" * 70)
print("VERDICT")
print("=" * 70)

# Evaluate Old World
ow_mean = ow_circ["mean_angular_offset_deg"]
ow_vp = ow_circ["v_test_vs_0"]["p_value"]
ow_pctile = ow_null["percentile_rank"]

if ow_mean < 30 and ow_vp < 0.05 and ow_pctile < 10:
    ow_verdict = "SUPPORTED"
    ow_explanation = (
        f"Old World mean offset={ow_mean:.1f}° (<30°), V-test p={ow_vp:.4f} (<0.05), "
        f"percentile={ow_pctile:.1f}% (<10%). Diffusion aligns with GC corridor."
    )
elif ow_mean > 60:
    ow_verdict = "COUNTER-EVIDENCE"
    ow_explanation = (
        f"Old World mean offset={ow_mean:.1f}° (>60°). Diffusion is predominantly "
        "perpendicular to GC, counter to corridor hypothesis."
    )
else:
    ow_verdict = "NULL"
    ow_explanation = (
        f"Old World mean offset={ow_mean:.1f}°, V-test p={ow_vp:.4f}, "
        f"percentile={ow_pctile:.1f}%. No significant alignment with GC corridor."
    )

# Evaluate New World
nw_mean = nw_circ["mean_angular_offset_deg"]
nw_vp = nw_circ["v_test_vs_0"]["p_value"]
nw_pctile = nw_null["percentile_rank"]

if nw_mean < 30 and nw_vp < 0.05 and nw_pctile < 10:
    nw_verdict = "SUPPORTED"
    nw_explanation = (
        f"New World mean offset={nw_mean:.1f}° (<30°), V-test p={nw_vp:.4f} (<0.05), "
        f"percentile={nw_pctile:.1f}% (<10%). Maize diffusion aligns with GC corridor."
    )
elif nw_mean > 60:
    nw_verdict = "COUNTER-EVIDENCE"
    nw_explanation = (
        f"New World mean offset={nw_mean:.1f}° (>60°). Maize diffusion is predominantly "
        "perpendicular to GC, counter to corridor hypothesis."
    )
else:
    nw_verdict = "NULL"
    nw_explanation = (
        f"New World mean offset={nw_mean:.1f}°, V-test p={nw_vp:.4f}, "
        f"percentile={nw_pctile:.1f}%. No significant alignment with GC corridor."
    )

# Overall
if ow_verdict == "SUPPORTED" and nw_verdict == "SUPPORTED":
    overall_verdict = "SUPPORTED"
    overall_explanation = "Both Old and New World crop diffusion show GC-aligned spread."
elif ow_verdict == "COUNTER-EVIDENCE" or nw_verdict == "COUNTER-EVIDENCE":
    overall_verdict = "COUNTER-EVIDENCE"
    overall_explanation = "At least one dataset shows perpendicular spread to GC."
else:
    overall_verdict = "NULL"
    overall_explanation = "Neither dataset shows statistically significant GC-aligned diffusion."

print(f"\n  Old World verdict:  {ow_verdict}")
print(f"    {ow_explanation}")
print(f"\n  New World verdict:  {nw_verdict}")
print(f"    {nw_explanation}")
print(f"\n  OVERALL:            {overall_verdict}")
print(f"    {overall_explanation}")


# ===========================================================================
#  SAVE OUTPUTS
# ===========================================================================
print("\n" + "=" * 70)
print("SAVING OUTPUTS")
print("=" * 70)

# Strip random_means from null results for JSON (too large)
ow_null_out = {k: v for k, v in ow_null.items() if k != "random_means"}
nw_null_out = {k: v for k, v in nw_null.items() if k != "random_means"}

# Strip offsets list from circ results for JSON readability
ow_circ_out = {k: v for k, v in ow_circ.items() if k != "offsets"}
nw_circ_out = {k: v for k, v in nw_circ.items() if k != "offsets"}

results = {
    "study": "S8: Archaeobotanical Diffusion Rate Test",
    "overall_verdict": overall_verdict,
    "overall_explanation": overall_explanation,
    "old_world": {
        "verdict": ow_verdict,
        "explanation": ow_explanation,
        "pairs": ow_pairs,
        "circular_stats": ow_circ_out,
        "speed_analysis": {k: v for k, v in ow_speed.items() if k != "pairs"},
        "speed_pairs": ow_speed["pairs"],
        "null_model": ow_null_out,
    },
    "new_world": {
        "verdict": nw_verdict,
        "explanation": nw_explanation,
        "pairs": nw_pairs,
        "circular_stats": nw_circ_out,
        "speed_analysis": {k: v for k, v in nw_speed.items() if k != "pairs"},
        "speed_pairs": nw_speed["pairs"],
        "null_model": nw_null_out,
    },
    "parameters": {
        "pole_lat": POLE_LAT, "pole_lon": POLE_LON,
        "mc_trials": MC_TRIALS,
        "random_seed": 42,
    },
}

with open(OUT_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("  -> saved results.json")

# ── RESULTS.md ───────────────────────────────────────────────────────────────
md = [
    "# Study 8: Archaeobotanical Diffusion Rate Test",
    "",
    f"**Overall Verdict: {overall_verdict}**",
    "",
    f"{overall_explanation}",
    "",
    "---",
    "",
    "## Old World: Wheat-Barley-Lentil Package",
    "",
    f"**Verdict: {ow_verdict}**",
    "",
    f"{ow_explanation}",
    "",
    "### Phase 1: Spread Bearings",
    "",
    f"- {len(ow_pairs)} consecutive pairs sorted chronologically",
    "",
    "| From | To | Offset (deg) | Distance (km) | Time (yrs) |",
    "|------|----|-------------|---------------|------------|",
]
for p in ow_pairs:
    md.append(f"| {p['from']} | {p['to']} | {p['angular_offset_deg']:.1f} | "
              f"{p['distance_km']:.1f} | {p['time_years']} |")

md += [
    "",
    "### Phase 2: Circular Statistics",
    "",
    f"- Mean angular offset: {ow_circ['mean_angular_offset_deg']:.2f} deg",
    f"- Median angular offset: {ow_circ['median_angular_offset_deg']:.2f} deg",
    f"- Std: {ow_circ['std_angular_offset_deg']:.2f} deg",
    f"- V-test vs 0 deg (along GC): V={ow_circ['v_test_vs_0']['V_statistic']:.4f}, "
    f"p={ow_circ['v_test_vs_0']['p_value']:.6f}",
    f"- V-test vs 90 deg (perp GC): V={ow_circ['v_test_vs_90']['V_statistic']:.4f}, "
    f"p={ow_circ['v_test_vs_90']['p_value']:.6f}",
    f"- Rayleigh test: Z={ow_circ['rayleigh_test']['Z_statistic']:.4f}, "
    f"p={ow_circ['rayleigh_test']['p_value']:.6f}",
    "",
    "### Phase 3: Speed Analysis",
    "",
    f"- Mean speed along GC: {ow_speed['mean_speed_along_km_yr']} km/yr",
    f"- Mean speed perp to GC: {ow_speed['mean_speed_perp_km_yr']} km/yr",
    f"- Median speed along GC: {ow_speed['median_speed_along_km_yr']} km/yr",
    f"- Median speed perp to GC: {ow_speed['median_speed_perp_km_yr']} km/yr",
    f"- Wilcoxon signed-rank (along > perp): stat={ow_speed['wilcoxon_stat']}, "
    f"p={ow_speed['wilcoxon_p']}",
    "",
    "### Phase 4: Null Model",
    "",
    f"- Alison GC mean offset: {ow_null['real_mean_offset_deg']:.2f} deg",
    f"- Random GCs mean: {ow_null['random_mean_offsets_mean']:.2f} deg "
    f"(std={ow_null['random_mean_offsets_std']:.2f})",
    f"- 5th-95th percentile: {ow_null['random_mean_offsets_5th']:.2f} - "
    f"{ow_null['random_mean_offsets_95th']:.2f} deg",
    f"- Alison percentile rank: {ow_null['percentile_rank']:.2f}%",
    "",
    "---",
    "",
    "## New World: Maize Diffusion",
    "",
    f"**Verdict: {nw_verdict}**",
    "",
    f"{nw_explanation}",
    "",
    "### Phase 1: Spread Bearings",
    "",
    f"- {len(nw_pairs)} consecutive pairs sorted chronologically",
    "",
    "| From | To | Offset (deg) | Distance (km) | Time (yrs) |",
    "|------|----|-------------|---------------|------------|",
]
for p in nw_pairs:
    md.append(f"| {p['from']} | {p['to']} | {p['angular_offset_deg']:.1f} | "
              f"{p['distance_km']:.1f} | {p['time_years']} |")

md += [
    "",
    "### Phase 2: Circular Statistics",
    "",
    f"- Mean angular offset: {nw_circ['mean_angular_offset_deg']:.2f} deg",
    f"- Median angular offset: {nw_circ['median_angular_offset_deg']:.2f} deg",
    f"- Std: {nw_circ['std_angular_offset_deg']:.2f} deg",
    f"- V-test vs 0 deg (along GC): V={nw_circ['v_test_vs_0']['V_statistic']:.4f}, "
    f"p={nw_circ['v_test_vs_0']['p_value']:.6f}",
    f"- V-test vs 90 deg (perp GC): V={nw_circ['v_test_vs_90']['V_statistic']:.4f}, "
    f"p={nw_circ['v_test_vs_90']['p_value']:.6f}",
    f"- Rayleigh test: Z={nw_circ['rayleigh_test']['Z_statistic']:.4f}, "
    f"p={nw_circ['rayleigh_test']['p_value']:.6f}",
    "",
    "### Phase 3: Speed Analysis",
    "",
    f"- Mean speed along GC: {nw_speed['mean_speed_along_km_yr']} km/yr",
    f"- Mean speed perp to GC: {nw_speed['mean_speed_perp_km_yr']} km/yr",
    f"- Median speed along GC: {nw_speed['median_speed_along_km_yr']} km/yr",
    f"- Median speed perp to GC: {nw_speed['median_speed_perp_km_yr']} km/yr",
    f"- Wilcoxon signed-rank (along > perp): stat={nw_speed['wilcoxon_stat']}, "
    f"p={nw_speed['wilcoxon_p']}",
    "",
    "### Phase 4: Null Model",
    "",
    f"- Alison GC mean offset: {nw_null['real_mean_offset_deg']:.2f} deg",
    f"- Random GCs mean: {nw_null['random_mean_offsets_mean']:.2f} deg "
    f"(std={nw_null['random_mean_offsets_std']:.2f})",
    f"- 5th-95th percentile: {nw_null['random_mean_offsets_5th']:.2f} - "
    f"{nw_null['random_mean_offsets_95th']:.2f} deg",
    f"- Alison percentile rank: {nw_null['percentile_rank']:.2f}%",
    "",
    "---",
    "",
    "## Parameters",
    "",
    f"- Pole: ({POLE_LAT}, {POLE_LON})",
    f"- Monte Carlo trials: {MC_TRIALS}",
    f"- Random seed: 42",
    "",
    "## Criteria",
    "",
    "- SUPPORTED: Mean angular offset < 30 deg AND V-test p < 0.05 AND Alison ranks in bottom 10% of random GCs",
    "- NULL: Mean angular offset ~45 deg AND V-test not significant",
    "- COUNTER-EVIDENCE: Mean angular offset > 60 deg",
    "",
]

with open(OUT_DIR / "RESULTS.md", "w") as f:
    f.write("\n".join(md))
print("  -> saved RESULTS.md")

print("\nDone.")
