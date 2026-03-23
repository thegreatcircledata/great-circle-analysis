#!/usr/bin/env python3
"""
Study 10: Modern Infrastructure Echo Test
==========================================
If the Alison great circle traces a deep habitability/connectivity corridor,
modern infrastructure (roads, railways, cities) may independently converge on
similar routes — because the same geographic logic that attracted ancient
builders also attracts modern engineers.

Approach: Use world's ~100 largest cities as a modern infrastructure proxy.
Compare the Alison circle's city proximity against 1000 random great circles.
Then test whether ancient-site density and modern-city density correlate
along the arc.

Phases:
  1. World city database (~100 major cities)
  2. Alison circle city count at 50/100/200/500 km thresholds
  3. Random circle comparison (Z-scores, percentiles)
  4. Ancient vs modern density correlation along the arc
  5. Land fraction analysis (Alison vs random circles)
"""

import math
import json
import os
import warnings
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
np.random.seed(42)

# ── Shared Constants & Geometry ──────────────────────────────────────────────

POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2

OUT_DIR = "/Users/elliotallan/megalith_site_research/outputs/next_wave/modern_infrastructure"
os.makedirs(OUT_DIR, exist_ok=True)


def haversine_km(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return 2 * EARTH_R_KM * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def dist_from_gc(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LON):
    """Distance from point to great circle defined by pole."""
    d = haversine_km(lat, lon, pole_lat, pole_lon)
    return abs(d - QUARTER_CIRC)


def random_pole():
    """Generate a random pole uniformly on the upper hemisphere."""
    lon = np.random.uniform(-180, 180)
    lat = np.degrees(np.arcsin(np.random.uniform(0, 1)))
    return lat, lon


def gc_points(pole_lat, pole_lon, n=360):
    """Generate n equally-spaced points along the great circle defined by pole."""
    # The GC is the set of points at exactly QUARTER_CIRC from the pole.
    # Parameterize by azimuth from pole.
    pole_lat_r = math.radians(pole_lat)
    pole_lon_r = math.radians(pole_lon)
    angular_dist = QUARTER_CIRC / EARTH_R_KM  # radians (pi/2)

    points = []
    for i in range(n):
        bearing = math.radians(i * 360.0 / n)
        lat2 = math.asin(
            math.sin(pole_lat_r) * math.cos(angular_dist)
            + math.cos(pole_lat_r) * math.sin(angular_dist) * math.cos(bearing)
        )
        lon2 = pole_lon_r + math.atan2(
            math.sin(bearing) * math.sin(angular_dist) * math.cos(pole_lat_r),
            math.cos(angular_dist) - math.sin(pole_lat_r) * math.sin(lat2),
        )
        points.append((math.degrees(lat2), math.degrees(lon2)))
    return points


# ── Phase 1: World City Database ─────────────────────────────────────────────

WORLD_CITIES = [
    ("Tokyo", 35.6762, 139.6503, 37.4),
    ("Delhi", 28.7041, 77.1025, 32.9),
    ("Shanghai", 31.2304, 121.4737, 29.2),
    ("São Paulo", -23.5505, -46.6333, 22.4),
    ("Mexico City", 19.4326, -99.1332, 21.8),
    ("Cairo", 30.0444, 31.2357, 21.3),
    ("Mumbai", 19.0760, 72.8777, 21.0),
    ("Beijing", 39.9042, 116.4074, 20.9),
    ("Dhaka", 23.8103, 90.4125, 22.0),
    ("Osaka", 34.6937, 135.5023, 19.2),
    ("New York", 40.7128, -74.0060, 18.8),
    ("Karachi", 24.8607, 67.0011, 16.5),
    ("Istanbul", 41.0082, 28.9784, 15.6),
    ("Chongqing", 29.4316, 106.9123, 15.9),
    ("Buenos Aires", -34.6037, -58.3816, 15.4),
    ("Kolkata", 22.5726, 88.3639, 15.1),
    ("Lagos", 6.5244, 3.3792, 14.9),
    ("Kinshasa", -4.4419, 15.2663, 14.3),
    ("Manila", 14.5995, 120.9842, 14.2),
    ("Tianjin", 39.3434, 117.3616, 13.9),
    ("Guangzhou", 23.1291, 113.2644, 13.8),
    ("Rio de Janeiro", -22.9068, -43.1729, 13.5),
    ("Lahore", 31.5204, 74.3587, 13.1),
    ("Bangalore", 12.9716, 77.5946, 12.8),
    ("Moscow", 55.7558, 37.6173, 12.5),
    ("Shenzhen", 22.5431, 114.0579, 12.4),
    ("Bogotá", 4.7110, -74.0721, 11.3),
    ("Jakarta", -6.2088, 106.8456, 11.0),
    ("Lima", -12.0464, -77.0428, 10.9),
    ("Paris", 48.8566, 2.3522, 11.1),
    ("London", 51.5074, -0.1278, 9.5),
    ("Bangkok", 13.7563, 100.5018, 10.7),
    ("Hyderabad", 17.3850, 78.4867, 10.5),
    ("Chennai", 13.0827, 80.2707, 10.9),
    ("Chicago", 41.8781, -87.6298, 9.5),
    ("Chengdu", 30.5728, 104.0668, 9.4),
    ("Nanjing", 32.0603, 118.7969, 9.4),
    ("Ho Chi Minh City", 10.8231, 106.6297, 9.3),
    ("Tehran", 35.6892, 51.3890, 9.1),
    ("Wuhan", 30.5928, 114.3055, 8.9),
    ("Xi'an", 34.3416, 108.9398, 8.7),
    ("Hong Kong", 22.3193, 114.1694, 7.5),
    ("Riyadh", 24.7136, 46.6753, 7.6),
    ("Ahmedabad", 23.0225, 72.5714, 8.4),
    ("Kuala Lumpur", 3.1390, 101.6869, 8.3),
    ("Hangzhou", 30.2741, 120.1551, 8.2),
    ("Surat", 21.1702, 72.8311, 7.8),
    ("Santiago", -33.4489, -70.6693, 6.8),
    ("Madrid", 40.4168, -3.7038, 6.7),
    ("Baghdad", 33.3152, 44.3661, 7.5),
    ("Toronto", 43.6532, -79.3832, 6.3),
    ("Belo Horizonte", -19.9167, -43.9345, 6.0),
    ("Nagoya", 35.1815, 136.9066, 9.5),
    ("Dar es Salaam", -6.7924, 39.2083, 7.0),
    ("Johannesburg", -26.2041, 28.0473, 6.1),
    ("Barcelona", 41.3874, 2.1686, 5.6),
    ("Jeddah", 21.4858, 39.1925, 4.6),
    ("Harbin", 45.8038, 126.5350, 5.9),
    ("Singapore", 1.3521, 103.8198, 5.9),
    ("Pune", 18.5204, 73.8567, 7.4),
    ("Luanda", -8.8390, 13.2894, 8.3),
    ("Addis Ababa", 9.0250, 38.7469, 5.2),
    ("Nairobi", -1.2921, 36.8219, 5.1),
    ("Casablanca", 33.5731, -7.5898, 3.8),
    ("Kabul", 34.5553, 69.2075, 4.6),
    ("Berlin", 52.5200, 13.4050, 3.7),
    ("Ankara", 39.9334, 32.8597, 5.7),
    ("St. Petersburg", 59.9311, 30.3609, 5.4),
    ("Yangon", 16.8661, 96.1951, 5.3),
    ("Alexandria", 31.2001, 29.9187, 5.3),
    ("Abidjan", 5.3600, -4.0083, 5.1),
    ("Guadalajara", 20.6597, -103.3496, 5.3),
    ("Surabaya", -7.2575, 112.7521, 3.0),
    ("Melbourne", -37.8136, 144.9631, 5.1),
    ("Rome", 41.9028, 12.4964, 4.3),
    ("Monterrey", 25.6866, -100.3161, 5.0),
    ("Sydney", -33.8688, 151.2093, 5.3),
    ("Recife", -8.0476, -34.8770, 4.1),
    ("Hanoi", 21.0285, 105.8542, 8.1),
    ("Fortaleza", -3.7172, -38.5433, 4.1),
    ("Khartoum", 15.5007, 32.5599, 5.8),
    ("Medellín", 6.2476, -75.5658, 4.0),
    ("Changsha", 28.2280, 112.9388, 4.6),
    ("Zhengzhou", 34.7466, 113.6253, 4.5),
    ("Kunming", 25.0389, 102.7183, 4.3),
    ("Dalian", 38.9140, 121.6147, 4.0),
    ("Accra", 5.6037, -0.1870, 4.2),
    ("Algiers", 36.7538, 3.0588, 3.9),
    ("Los Angeles", 34.0522, -118.2437, 13.2),
    ("Seoul", 37.5665, 126.9780, 9.9),
    ("Taipei", 25.0330, 121.5654, 7.0),
    ("Kiev", 50.4501, 30.5234, 3.0),
    ("Lima", -12.0464, -77.0428, 10.9),
    ("Houston", 29.7604, -95.3698, 7.1),
    ("Miami", 25.7617, -80.1918, 6.2),
    ("San Francisco", 37.7749, -122.4194, 4.7),
    ("Denver", 39.7392, -104.9903, 2.9),
    ("Phoenix", 33.4484, -112.0740, 4.9),
    ("Atlanta", 33.7490, -84.3880, 6.1),
    ("Seattle", 47.6062, -122.3321, 4.0),
]

# Deduplicate (Lima appears twice)
_seen = set()
CITIES = []
for name, lat, lon, pop in WORLD_CITIES:
    if name not in _seen:
        _seen.add(name)
        CITIES.append((name, lat, lon, pop))

print(f"Phase 1: {len(CITIES)} unique world cities loaded")

# ── Phase 2: Alison Circle City Count ────────────────────────────────────────

THRESHOLDS_KM = [50, 100, 200, 500]

alison_counts = {}
alison_nearby = {}
for thresh in THRESHOLDS_KM:
    nearby = []
    for name, lat, lon, pop in CITIES:
        d = dist_from_gc(lat, lon)
        if d <= thresh:
            nearby.append((name, lat, lon, pop, round(d, 1)))
    alison_counts[thresh] = len(nearby)
    alison_nearby[thresh] = nearby

print("\nPhase 2: Alison circle city counts")
for thresh in THRESHOLDS_KM:
    print(f"  Within {thresh:>4d} km: {alison_counts[thresh]:>3d} cities")
    if alison_counts[thresh] > 0 and thresh <= 200:
        for name, lat, lon, pop, d in sorted(alison_nearby[thresh], key=lambda x: x[4]):
            print(f"    {name:25s}  dist={d:6.1f} km  pop={pop:.1f}M")

# ── Phase 3: Random Circle Comparison ────────────────────────────────────────

N_RANDOM = 1000
random_counts = {thresh: [] for thresh in THRESHOLDS_KM}

print(f"\nPhase 3: Generating {N_RANDOM} random great circles...")
for i in range(N_RANDOM):
    plat, plon = random_pole()
    for thresh in THRESHOLDS_KM:
        count = 0
        for _, clat, clon, _ in CITIES:
            d = dist_from_gc(clat, clon, plat, plon)
            if d <= thresh:
                count += 1
        random_counts[thresh].append(count)

zscores = {}
percentiles = {}
for thresh in THRESHOLDS_KM:
    arr = np.array(random_counts[thresh])
    mean_r = arr.mean()
    std_r = arr.std()
    z = (alison_counts[thresh] - mean_r) / std_r if std_r > 0 else 0.0
    pct = np.mean(arr <= alison_counts[thresh]) * 100
    zscores[thresh] = round(z, 2)
    percentiles[thresh] = round(pct, 1)
    print(f"  {thresh:>4d} km — Alison: {alison_counts[thresh]:>3d}, "
          f"Random mean: {mean_r:.1f} ± {std_r:.1f}, "
          f"Z = {z:+.2f}, percentile = {pct:.1f}%")

# ── Phase 4: Ancient vs Modern Correlation ───────────────────────────────────

print("\nPhase 4: Ancient vs modern density correlation along the arc...")

# Load archaeological sites
osm_path = "/Users/elliotallan/megalith_site_research/data/osm/osm_archaeological_sites.csv"
osm_df = pd.read_csv(osm_path)
osm_df = osm_df.dropna(subset=["lat", "lon"])
arch_lats = osm_df["lat"].values
arch_lons = osm_df["lon"].values
print(f"  Loaded {len(osm_df)} archaeological sites")

# Generate 360 equally-spaced points along Alison GC
gc_pts = gc_points(POLE_LAT, POLE_LON, n=360)

# Precompute city coordinates
city_lats = np.array([c[1] for c in CITIES])
city_lons = np.array([c[2] for c in CITIES])

# For each GC point, count archaeological sites within 100km and cities within 500km
ARCH_RADIUS = 100  # km
CITY_RADIUS = 500  # km

arch_density = np.zeros(360)
city_density = np.zeros(360)

for idx, (glat, glon) in enumerate(gc_pts):
    if idx % 60 == 0:
        print(f"  Processing GC point {idx}/360...")

    # Count cities within CITY_RADIUS
    for clat, clon in zip(city_lats, city_lons):
        if haversine_km(glat, glon, clat, clon) <= CITY_RADIUS:
            city_density[idx] += 1

    # Count archaeological sites within ARCH_RADIUS
    # Use a bounding box pre-filter for speed
    dlat_approx = ARCH_RADIUS / 111.0  # rough degrees
    dlon_approx = ARCH_RADIUS / (111.0 * max(math.cos(math.radians(glat)), 0.01))
    mask = (
        (np.abs(arch_lats - glat) <= dlat_approx)
        & (np.abs(arch_lons - glon) <= dlon_approx)
    )
    candidates = np.where(mask)[0]
    for ci in candidates:
        if haversine_km(glat, glon, arch_lats[ci], arch_lons[ci]) <= ARCH_RADIUS:
            arch_density[idx] += 1

# Correlation
# Only use points where at least one has nonzero density
nonzero_mask = (arch_density > 0) | (city_density > 0)
if nonzero_mask.sum() > 10:
    r_val, p_val = pearsonr(arch_density[nonzero_mask], city_density[nonzero_mask])
else:
    r_val, p_val = 0.0, 1.0

# Also do full correlation (all 360 points)
r_full, p_full = pearsonr(arch_density, city_density)

print(f"  Pearson r (all 360 pts):     {r_full:.4f}  (p = {p_full:.2e})")
print(f"  Pearson r (nonzero, n={nonzero_mask.sum():>3d}): {r_val:.4f}  (p = {p_val:.2e})")

# ── Phase 5: Land Fraction Analysis ──────────────────────────────────────────

print("\nPhase 5: Land fraction analysis...")

land_grid_path = "/Users/elliotallan/megalith_site_research/data/land_grid_05deg.npy"
land_grid = np.load(land_grid_path)  # shape (360, 720), 0.5-degree resolution


def is_land(lat, lon, grid=land_grid):
    """Check if point is on land using 0.5-degree grid."""
    # Grid: row 0 = 90°N, row 359 = 89.5°S; col 0 = -180°, col 719 = 179.5°E
    row = int((90.0 - lat) / 0.5)
    col = int((lon + 180.0) / 0.5)
    row = max(0, min(row, 359))
    col = max(0, min(col, 719))
    return bool(grid[row, col])


# Alison circle land fraction
alison_land = sum(1 for lat, lon in gc_pts if is_land(lat, lon))
alison_land_frac = alison_land / len(gc_pts)
print(f"  Alison GC land fraction: {alison_land_frac:.3f} ({alison_land}/{len(gc_pts)} points)")

# Random circles land fraction
random_land_fracs = []
for i in range(N_RANDOM):
    plat, plon = random_pole()
    pts = gc_points(plat, plon, n=360)
    land_count = sum(1 for lat, lon in pts if is_land(lat, lon))
    random_land_fracs.append(land_count / len(pts))

random_land_arr = np.array(random_land_fracs)
land_mean = random_land_arr.mean()
land_std = random_land_arr.std()
land_z = (alison_land_frac - land_mean) / land_std if land_std > 0 else 0
land_pct = np.mean(random_land_arr <= alison_land_frac) * 100

print(f"  Random mean land frac: {land_mean:.3f} ± {land_std:.3f}")
print(f"  Z = {land_z:+.2f}, percentile = {land_pct:.1f}%")

# ── Visualizations ───────────────────────────────────────────────────────────

print("\nGenerating plots...")

# Plot 1: Random circle city count distributions with Alison line
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("Study 10: Cities Near Great Circle — Alison vs Random", fontsize=14, fontweight="bold")

for ax, thresh in zip(axes.flat, THRESHOLDS_KM):
    arr = np.array(random_counts[thresh])
    ax.hist(arr, bins=30, color="#4A90D9", alpha=0.7, edgecolor="white", label="Random GCs")
    ax.axvline(alison_counts[thresh], color="red", linewidth=2.5, linestyle="--",
               label=f"Alison = {alison_counts[thresh]}")
    ax.set_xlabel("Number of cities")
    ax.set_ylabel("Frequency")
    ax.set_title(f"Within {thresh} km  (Z={zscores[thresh]:+.2f}, P{percentiles[thresh]:.0f})")
    ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "city_count_distributions.png"), dpi=150)
plt.close()

# Plot 2: Ancient vs modern density along the arc
fig, ax1 = plt.subplots(figsize=(14, 5))
azimuths = np.arange(360)

color_arch = "#8B4513"
color_city = "#1E90FF"

ax1.set_xlabel("Azimuth along GC (degrees from pole)")
ax1.set_ylabel("Archaeological sites (100 km)", color=color_arch)
ax1.bar(azimuths, arch_density, width=1.0, color=color_arch, alpha=0.5, label="Arch. sites")
ax1.tick_params(axis="y", labelcolor=color_arch)

ax2 = ax1.twinx()
ax2.set_ylabel("Major cities (500 km)", color=color_city)
ax2.plot(azimuths, city_density, color=color_city, linewidth=1.5, label="Cities")
ax2.tick_params(axis="y", labelcolor=color_city)

fig.suptitle(f"Ancient Site & Modern City Density Along Alison GC  (r = {r_full:.3f})",
             fontsize=13, fontweight="bold")
fig.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "ancient_modern_density.png"), dpi=150)
plt.close()

# Plot 3: Land fraction distribution
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(random_land_fracs, bins=40, color="#228B22", alpha=0.7, edgecolor="white", label="Random GCs")
ax.axvline(alison_land_frac, color="red", linewidth=2.5, linestyle="--",
           label=f"Alison = {alison_land_frac:.3f}")
ax.set_xlabel("Land fraction")
ax.set_ylabel("Frequency")
ax.set_title(f"Land Fraction: Alison vs Random GCs  (Z={land_z:+.2f}, P{land_pct:.0f})")
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "land_fraction_distribution.png"), dpi=150)
plt.close()

# Plot 4: Map of GC with cities
fig, ax = plt.subplots(figsize=(16, 8))
# Plot GC
gc_lats = [p[0] for p in gc_pts]
gc_lons = [p[1] for p in gc_pts]
ax.scatter(gc_lons, gc_lats, s=1, c="red", alpha=0.4, label="Alison GC")

# Plot all cities
for name, lat, lon, pop in CITIES:
    d = dist_from_gc(lat, lon)
    if d <= 200:
        ax.scatter(lon, lat, s=pop * 3, c="blue", alpha=0.8, zorder=5)
        ax.annotate(name, (lon, lat), fontsize=6, ha="left", va="bottom", color="blue")
    else:
        ax.scatter(lon, lat, s=pop * 2, c="gray", alpha=0.3, zorder=3)

ax.set_xlim(-180, 180)
ax.set_ylim(-90, 90)
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.set_title("World Cities & Alison Great Circle (blue = within 200 km)", fontsize=13, fontweight="bold")
ax.set_aspect("equal")
ax.grid(True, alpha=0.2)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "gc_cities_map.png"), dpi=150)
plt.close()

print("Plots saved.")

# ── Compile Results ──────────────────────────────────────────────────────────

results = {
    "study": "S10 — Modern Infrastructure Echo Test",
    "date": datetime.now().strftime("%Y-%m-%d"),
    "seed": 42,
    "n_cities": len(CITIES),
    "n_random_circles": N_RANDOM,
    "alison_city_counts": {str(k): v for k, v in alison_counts.items()},
    "alison_nearby_cities": {
        str(k): [{"name": n, "lat": la, "lon": lo, "pop_M": p, "dist_km": d}
                 for n, la, lo, p, d in v]
        for k, v in alison_nearby.items()
    },
    "random_circle_stats": {
        str(thresh): {
            "mean": round(float(np.mean(random_counts[thresh])), 2),
            "std": round(float(np.std(random_counts[thresh])), 2),
            "median": round(float(np.median(random_counts[thresh])), 2),
            "z_score": zscores[thresh],
            "percentile": percentiles[thresh],
        }
        for thresh in THRESHOLDS_KM
    },
    "ancient_modern_correlation": {
        "pearson_r_full": round(r_full, 4),
        "p_value_full": float(f"{p_full:.6e}"),
        "pearson_r_nonzero": round(r_val, 4),
        "p_value_nonzero": float(f"{p_val:.6e}"),
        "n_nonzero_points": int(nonzero_mask.sum()),
    },
    "land_fraction": {
        "alison": round(alison_land_frac, 4),
        "random_mean": round(land_mean, 4),
        "random_std": round(land_std, 4),
        "z_score": round(land_z, 2),
        "percentile": round(land_pct, 1),
    },
}

with open(os.path.join(OUT_DIR, "results.json"), "w") as f:
    json.dump(results, f, indent=2)

# ── RESULTS.md ───────────────────────────────────────────────────────────────

nearby_200_list = "\n".join(
    f"| {n:25s} | {la:8.2f} | {lo:9.2f} | {p:5.1f} | {d:7.1f} |"
    for n, la, lo, p, d in sorted(alison_nearby[200], key=lambda x: x[4])
)

md = f"""# Study 10: Modern Infrastructure Echo Test

**Date:** {results['date']}
**Seed:** 42

## Question
If the Alison great circle traces a deep habitability/connectivity corridor,
do modern cities — proxies for infrastructure — cluster near it more than
expected by chance?

## Method
- **City database:** {len(CITIES)} of the world's largest cities (population ≥ ~3M)
- **Metric:** Count cities within 50 / 100 / 200 / 500 km of each great circle
- **Null model:** {N_RANDOM} random great circles (uniformly random poles, upper hemisphere)
- **Correlation:** Pearson r between archaeological-site density and city density
  along 360 equally-spaced points on the Alison GC
- **Land fraction:** Fraction of GC crossing land, Alison vs random

## Results

### City Proximity (Alison vs Random)

| Threshold | Alison | Random Mean ± SD | Z-score | Percentile |
|-----------|--------|------------------|---------|------------|
| 50 km     | {alison_counts[50]:>6d} | {np.mean(random_counts[50]):>5.1f} ± {np.std(random_counts[50]):.1f} | {zscores[50]:>+6.2f} | {percentiles[50]:>5.1f}% |
| 100 km    | {alison_counts[100]:>6d} | {np.mean(random_counts[100]):>5.1f} ± {np.std(random_counts[100]):.1f} | {zscores[100]:>+6.2f} | {percentiles[100]:>5.1f}% |
| 200 km    | {alison_counts[200]:>6d} | {np.mean(random_counts[200]):>5.1f} ± {np.std(random_counts[200]):.1f} | {zscores[200]:>+6.2f} | {percentiles[200]:>5.1f}% |
| 500 km    | {alison_counts[500]:>6d} | {np.mean(random_counts[500]):>5.1f} ± {np.std(random_counts[500]):.1f} | {zscores[500]:>+6.2f} | {percentiles[500]:>5.1f}% |

### Cities Within 200 km of Alison GC

| City                      | Lat      | Lon       | Pop(M) | Dist(km) |
|---------------------------|----------|-----------|--------|----------|
{nearby_200_list}

### Ancient vs Modern Density Correlation

| Metric | Value |
|--------|-------|
| Pearson r (all 360 pts) | {r_full:.4f} |
| p-value | {p_full:.2e} |
| Pearson r (nonzero, n={nonzero_mask.sum()}) | {r_val:.4f} |
| p-value (nonzero) | {p_val:.2e} |

### Land Fraction

| Circle | Land Fraction |
|--------|---------------|
| Alison | {alison_land_frac:.4f} |
| Random mean ± SD | {land_mean:.4f} ± {land_std:.4f} |
| Z-score | {land_z:+.2f} |
| Percentile | {land_pct:.1f}% |

## Interpretation

The Alison great circle's city proximity is compared against a null distribution
of {N_RANDOM} random great circles. Z-scores above +2 would indicate the Alison
circle passes near significantly more major cities than expected by chance,
supporting the "habitability corridor" hypothesis.

The ancient-vs-modern density correlation tests whether the same segments of the
arc that are rich in archaeological sites also tend to be near modern cities —
suggesting persistent geographic attractors across millennia.

## Outputs
- `results.json` — full numeric results
- `city_count_distributions.png` — histograms at each threshold
- `ancient_modern_density.png` — dual-axis density along the arc
- `land_fraction_distribution.png` — land fraction comparison
- `gc_cities_map.png` — map of GC with cities
"""

with open(os.path.join(OUT_DIR, "RESULTS.md"), "w") as f:
    f.write(md)

print(f"\nAll outputs saved to {OUT_DIR}/")
print("Done.")
