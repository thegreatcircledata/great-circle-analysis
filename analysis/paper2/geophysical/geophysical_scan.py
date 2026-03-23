#!/usr/bin/env python3
"""
Geophysical Property Scan
==========================
Test whether the Alison Great Circle correlates with measurable geophysical
properties: geoid anomaly, magnetic field, topographic consistency, and
ocean current alignment.

Data sources:
- EGM96 5-arcminute geoid grid (downloaded PGM file)
- ETOPO elevation via Open-Elevation API (SRTM/ASTER)
- IGRF magnetic field model via NOAA API
- Ocean currents: analytical approach using known major current systems
"""

import json, math, os, sys, time, struct
import numpy as np
import requests
from geographiclib.geodesic import Geodesic

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE = os.path.expanduser("~/megalith_site_research")
POLE_LAT, POLE_LNG = 59.682122, -138.646087
EARTH_R = 6371.0
N_RANDOM = 100
N_POINTS = 360  # 1 per degree of arc

np.random.seed(42)
geod = Geodesic.WGS84


# ============================================================
# GENERATE GREAT CIRCLE POINTS
# ============================================================
def great_circle_points(pole_lat, pole_lon, n_points=360):
    """Generate n_points evenly spaced along a great circle defined by its pole."""
    points = []
    for i in range(n_points):
        bearing = i * (360.0 / n_points)
        # Go quarter-circumference from pole
        result = geod.Direct(pole_lat, pole_lon, bearing, EARTH_R * 1000 * math.pi / 2)
        lat, lon = result['lat2'], result['lon2']
        points.append((lat, lon))
    return np.array(points)


def random_pole():
    """Generate a random pole (uniformly distributed on sphere)."""
    z = np.random.uniform(-1, 1)
    theta = np.random.uniform(0, 2 * np.pi)
    lat = np.degrees(np.arcsin(z))
    lon = np.degrees(theta) - 180
    return lat, lon


print("Generating circle points...")
alison_pts = great_circle_points(POLE_LAT, POLE_LNG, N_POINTS)
print(f"Alison circle: {len(alison_pts)} points, lat range [{alison_pts[:,0].min():.1f}, {alison_pts[:,0].max():.1f}]")

random_poles = [random_pole() for _ in range(N_RANDOM)]
print(f"Generated {N_RANDOM} random poles")


# ============================================================
# TEST 1: GEOID ANOMALY (EGM96)
# ============================================================
print(f"\n{'='*70}")
print("TEST 1: GEOID ANOMALY (EGM96)")
print("="*70)

# Load EGM96 grid
pgm_path = os.path.join(BASE, "data/geophysical/geoids/egm96-5.pgm")

with open(pgm_path, 'rb') as f:
    offset, scale = -108.0, 0.003
    # Skip header
    while True:
        line = f.readline()
        if line.startswith(b'#'):
            text = line.decode().strip()
            if 'Offset' in text:
                offset = float(text.split()[-1])
            if 'Scale' in text:
                scale = float(text.split()[-1])
        elif not line.startswith(b'P'):
            # dimensions line
            parts = line.strip().split()
            width, height = int(parts[0]), int(parts[1])
            maxval = int(f.readline().strip())
            break

    raw = f.read(width * height * 2)
    geoid_grid = offset + scale * np.frombuffer(raw, dtype='>u2').reshape(height, width).astype(float)

print(f"EGM96 grid loaded: {width}x{height}, range [{geoid_grid.min():.1f}, {geoid_grid.max():.1f}] m")

# Resolution: 5 arcminutes = 1/12 degree
GEOID_RES = 12  # cells per degree


def sample_geoid(lat, lon):
    """Bilinear interpolation of geoid height at (lat, lon)."""
    # Grid origin: 90N, 0E; going south and east
    lon = lon % 360  # ensure 0-360
    y = (90 - lat) * GEOID_RES
    x = lon * GEOID_RES
    y0, x0 = int(y), int(x)
    y1, x1 = min(y0 + 1, height - 1), (x0 + 1) % width
    fy, fx = y - y0, x - x0
    val = (geoid_grid[y0, x0] * (1-fx)*(1-fy) +
           geoid_grid[y0, x1] * fx*(1-fy) +
           geoid_grid[y1, x0] * (1-fx)*fy +
           geoid_grid[y1, x1] * fx*fy)
    return float(val)


def geoid_stats(pts):
    vals = np.array([sample_geoid(lat, lon) for lat, lon in pts])
    return {
        "mean": float(np.mean(vals)),
        "std": float(np.std(vals)),
        "min": float(np.min(vals)),
        "max": float(np.max(vals)),
        "range": float(np.ptp(vals)),
        "gradient_mean": float(np.mean(np.abs(np.diff(vals)))),
        "values": vals
    }


# Alison circle geoid
t0 = time.time()
alison_geoid = geoid_stats(alison_pts)
print(f"Alison geoid: mean={alison_geoid['mean']:.1f}m, std={alison_geoid['std']:.1f}m, "
      f"range={alison_geoid['range']:.1f}m, gradient={alison_geoid['gradient_mean']:.2f}m/deg")

# Random circles
random_geoid_means = []
random_geoid_stds = []
random_geoid_ranges = []
random_geoid_gradients = []

for i, (plat, plon) in enumerate(random_poles):
    pts = great_circle_points(plat, plon, N_POINTS)
    stats = geoid_stats(pts)
    random_geoid_means.append(stats["mean"])
    random_geoid_stds.append(stats["std"])
    random_geoid_ranges.append(stats["range"])
    random_geoid_gradients.append(stats["gradient_mean"])
    if (i + 1) % 25 == 0:
        print(f"  {i+1}/{N_RANDOM} random circles done [{time.time()-t0:.1f}s]")

# Z-scores
def z_of(val, arr):
    return (val - np.mean(arr)) / np.std(arr) if np.std(arr) > 0 else 0

geoid_z_mean = z_of(alison_geoid["mean"], random_geoid_means)
geoid_z_std = z_of(alison_geoid["std"], random_geoid_stds)
geoid_z_range = z_of(alison_geoid["range"], random_geoid_ranges)
geoid_z_gradient = z_of(alison_geoid["gradient_mean"], random_geoid_gradients)

print(f"\nGeoid Z-scores (Alison vs {N_RANDOM} random circles):")
print(f"  Mean geoid height: Z = {geoid_z_mean:+.2f} (Alison={alison_geoid['mean']:.1f}, "
      f"random={np.mean(random_geoid_means):.1f}±{np.std(random_geoid_means):.1f})")
print(f"  Std of geoid:      Z = {geoid_z_std:+.2f} (Alison={alison_geoid['std']:.1f}, "
      f"random={np.mean(random_geoid_stds):.1f}±{np.std(random_geoid_stds):.1f})")
print(f"  Range of geoid:    Z = {geoid_z_range:+.2f}")
print(f"  Mean gradient:     Z = {geoid_z_gradient:+.2f}")

geoid_percentile_mean = float(np.mean(np.array(random_geoid_means) < alison_geoid["mean"]) * 100)
print(f"  Alison mean geoid percentile: {geoid_percentile_mean:.1f}%")

test1_result = {
    "alison": {
        "mean": round(alison_geoid["mean"], 2),
        "std": round(alison_geoid["std"], 2),
        "range": round(alison_geoid["range"], 2),
        "gradient_mean": round(alison_geoid["gradient_mean"], 3)
    },
    "random_100": {
        "mean_of_means": round(float(np.mean(random_geoid_means)), 2),
        "std_of_means": round(float(np.std(random_geoid_means)), 2),
        "mean_of_stds": round(float(np.mean(random_geoid_stds)), 2),
        "mean_of_gradients": round(float(np.mean(random_geoid_gradients)), 3)
    },
    "z_scores": {
        "mean_height": round(geoid_z_mean, 2),
        "std_height": round(geoid_z_std, 2),
        "range": round(geoid_z_range, 2),
        "gradient": round(geoid_z_gradient, 2)
    },
    "percentile_mean": round(geoid_percentile_mean, 1),
    "interpretation": ""
}

if abs(geoid_z_mean) > 2:
    test1_result["interpretation"] = f"Alison circle follows a geoid {'ridge' if geoid_z_mean > 0 else 'trough'} (Z={geoid_z_mean:+.2f})"
elif abs(geoid_z_gradient) > 2:
    test1_result["interpretation"] = f"Alison circle follows an unusual geoid gradient (Z={geoid_z_gradient:+.2f})"
else:
    test1_result["interpretation"] = "No significant geoid anomaly correlation"


# ============================================================
# TEST 2: MAGNETIC FIELD (IGRF via NOAA API)
# ============================================================
print(f"\n{'='*70}")
print("TEST 2: CRUSTAL MAGNETIC FIELD (IGRF MODEL)")
print("="*70)
print("Note: Using IGRF model via NOAA API. This gives the main field, not")
print("crustal anomalies (EMAG2v3 grid would be needed for that).")
print("IGRF total intensity still tests whether the circle follows a")
print("magnetic field strength contour.")

# Sample every 10 degrees to keep API calls manageable (36 points per circle)
MAGNETIC_STEP = 10
mag_indices = list(range(0, N_POINTS, MAGNETIC_STEP))


def sample_magnetic_batch(pts, indices):
    """Sample IGRF total intensity at selected points via NOAA API."""
    values = []
    for idx in indices:
        lat, lon = pts[idx]
        url = (f"https://www.ngdc.noaa.gov/geomag-web/calculators/calculateIgrfwmm"
               f"?lat1={lat:.4f}&lon1={lon:.4f}&key=zNEw7&resultFormat=json&model=IGRF")
        try:
            r = requests.get(url, timeout=15)
            data = r.json()
            if 'result' in data and len(data['result']) > 0:
                ti = data['result'][0].get('totalintensity', None)
                if ti is not None:
                    values.append(float(ti))
                    continue
        except Exception:
            pass
        values.append(np.nan)
    return np.array(values)


print(f"\nSampling {len(mag_indices)} points along Alison circle...")
t0 = time.time()
alison_mag = sample_magnetic_batch(alison_pts, mag_indices)
valid = ~np.isnan(alison_mag)
print(f"  Got {valid.sum()}/{len(mag_indices)} valid readings [{time.time()-t0:.1f}s]")
if valid.sum() > 0:
    print(f"  Total intensity: mean={np.nanmean(alison_mag):.0f} nT, "
          f"std={np.nanstd(alison_mag):.0f} nT, "
          f"range=[{np.nanmin(alison_mag):.0f}, {np.nanmax(alison_mag):.0f}] nT")

# Sample 20 random circles (API rate limiting prevents 100)
N_MAG_RANDOM = 20
random_mag_means = []
random_mag_stds = []

print(f"\nSampling {N_MAG_RANDOM} random circles ({len(mag_indices)} points each)...")
for i in range(N_MAG_RANDOM):
    plat, plon = random_poles[i]
    pts = great_circle_points(plat, plon, N_POINTS)
    vals = sample_magnetic_batch(pts, mag_indices)
    v = vals[~np.isnan(vals)]
    if len(v) > 5:
        random_mag_means.append(float(np.mean(v)))
        random_mag_stds.append(float(np.std(v)))
    if (i + 1) % 5 == 0:
        print(f"  {i+1}/{N_MAG_RANDOM} done [{time.time()-t0:.1f}s]")

if random_mag_means and valid.sum() > 5:
    mag_z_mean = z_of(float(np.nanmean(alison_mag)), random_mag_means)
    mag_z_std = z_of(float(np.nanstd(alison_mag)), random_mag_stds)
    print(f"\nMagnetic Z-scores:")
    print(f"  Mean intensity: Z = {mag_z_mean:+.2f} (Alison={np.nanmean(alison_mag):.0f}, "
          f"random={np.mean(random_mag_means):.0f}±{np.std(random_mag_means):.0f})")
    print(f"  Std intensity:  Z = {mag_z_std:+.2f}")
else:
    mag_z_mean = mag_z_std = 0
    print("Insufficient magnetic data for comparison")

test2_result = {
    "note": "IGRF main field model (not crustal anomaly EMAG2v3). Tests field strength contour tracking.",
    "n_points_per_circle": len(mag_indices),
    "n_random_circles": N_MAG_RANDOM,
    "alison": {
        "mean_nT": round(float(np.nanmean(alison_mag)), 1) if valid.sum() > 0 else None,
        "std_nT": round(float(np.nanstd(alison_mag)), 1) if valid.sum() > 0 else None,
        "n_valid": int(valid.sum())
    },
    "z_scores": {
        "mean_intensity": round(mag_z_mean, 2),
        "std_intensity": round(mag_z_std, 2)
    }
}


# ============================================================
# TEST 3: TOPOGRAPHIC CONSISTENCY (Open-Elevation API + EGM96 as proxy)
# ============================================================
print(f"\n{'='*70}")
print("TEST 3: TOPOGRAPHIC CONSISTENCY")
print("="*70)
print("Using Open-Elevation API (SRTM/ASTER-based) for elevation sampling.")

# Batch query the API (supports POST with multiple locations)
def batch_elevation(pts, batch_size=100):
    """Get elevations for points via Open-Elevation API."""
    elevations = np.full(len(pts), np.nan)
    for start in range(0, len(pts), batch_size):
        batch = pts[start:start+batch_size]
        locations = [{"latitude": float(lat), "longitude": float(lon)} for lat, lon in batch]
        try:
            r = requests.post(
                "https://api.open-elevation.com/api/v1/lookup",
                json={"locations": locations},
                timeout=30
            )
            if r.status_code == 200:
                data = r.json()
                for j, result in enumerate(data.get("results", [])):
                    elevations[start + j] = result.get("elevation", np.nan)
        except Exception as e:
            print(f"  API error at batch {start}: {e}")
        time.sleep(0.5)  # rate limit
    return elevations


# Sample every 4 degrees (90 points) for API efficiency
TOPO_STEP = 4
topo_indices = list(range(0, N_POINTS, TOPO_STEP))
topo_pts = alison_pts[topo_indices]

print(f"Sampling {len(topo_pts)} points along Alison circle...")
t0 = time.time()
alison_elev = batch_elevation(topo_pts)
valid_e = ~np.isnan(alison_elev)
print(f"  Got {valid_e.sum()}/{len(topo_pts)} valid elevations [{time.time()-t0:.1f}s]")

if valid_e.sum() > 10:
    elev_valid = alison_elev[valid_e]
    ocean_frac = float(np.mean(elev_valid <= 0))
    low_frac = float(np.mean(elev_valid <= 200))
    print(f"  Elevation: mean={np.mean(elev_valid):.0f}m, std={np.std(elev_valid):.0f}m")
    print(f"  Ocean fraction: {ocean_frac*100:.1f}%")
    print(f"  Low elevation (≤200m): {low_frac*100:.1f}%")

    # Sample random circles (smaller set due to API limits)
    N_TOPO_RANDOM = 30
    random_elev_means = []
    random_elev_stds = []
    random_ocean_fracs = []
    random_low_fracs = []

    print(f"\nSampling {N_TOPO_RANDOM} random circles...")
    for i in range(N_TOPO_RANDOM):
        plat, plon = random_poles[i]
        rpts = great_circle_points(plat, plon, N_POINTS)
        rpts_sub = rpts[topo_indices]
        relev = batch_elevation(rpts_sub)
        rv = relev[~np.isnan(relev)]
        if len(rv) > 10:
            random_elev_means.append(float(np.mean(rv)))
            random_elev_stds.append(float(np.std(rv)))
            random_ocean_fracs.append(float(np.mean(rv <= 0)))
            random_low_fracs.append(float(np.mean(rv <= 200)))
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{N_TOPO_RANDOM} done [{time.time()-t0:.1f}s]")

    if random_elev_means:
        topo_z_mean = z_of(float(np.mean(elev_valid)), random_elev_means)
        topo_z_std = z_of(float(np.std(elev_valid)), random_elev_stds)
        topo_z_ocean = z_of(ocean_frac, random_ocean_fracs)
        topo_z_low = z_of(low_frac, random_low_fracs)

        print(f"\nTopographic Z-scores:")
        print(f"  Mean elevation: Z = {topo_z_mean:+.2f} (Alison={np.mean(elev_valid):.0f}m, "
              f"random={np.mean(random_elev_means):.0f}±{np.std(random_elev_means):.0f}m)")
        print(f"  Std elevation:  Z = {topo_z_std:+.2f} (Alison={np.std(elev_valid):.0f}m, "
              f"random={np.mean(random_elev_stds):.0f}m)")
        print(f"  Ocean fraction: Z = {topo_z_ocean:+.2f} (Alison={ocean_frac*100:.1f}%, "
              f"random={np.mean(random_ocean_fracs)*100:.1f}%)")
        print(f"  Low elev frac:  Z = {topo_z_low:+.2f} (Alison={low_frac*100:.1f}%, "
              f"random={np.mean(random_low_fracs)*100:.1f}%)")
    else:
        topo_z_mean = topo_z_std = topo_z_ocean = topo_z_low = 0
else:
    topo_z_mean = topo_z_std = topo_z_ocean = topo_z_low = 0
    ocean_frac = low_frac = 0
    random_elev_means = random_elev_stds = random_ocean_fracs = random_low_fracs = []
    print("Insufficient elevation data")

test3_result = {
    "n_points_per_circle": len(topo_indices),
    "n_random_circles": N_TOPO_RANDOM if random_elev_means else 0,
    "alison": {
        "mean_elevation_m": round(float(np.mean(alison_elev[valid_e])), 1) if valid_e.sum() > 0 else None,
        "std_elevation_m": round(float(np.std(alison_elev[valid_e])), 1) if valid_e.sum() > 0 else None,
        "ocean_fraction": round(ocean_frac, 3),
        "low_elevation_fraction": round(low_frac, 3),
        "n_valid": int(valid_e.sum())
    },
    "z_scores": {
        "mean_elevation": round(topo_z_mean, 2),
        "std_elevation": round(topo_z_std, 2),
        "ocean_fraction": round(topo_z_ocean, 2),
        "low_elevation_fraction": round(topo_z_low, 2)
    }
}


# ============================================================
# TEST 4: OCEAN CURRENT ALIGNMENT
# ============================================================
print(f"\n{'='*70}")
print("TEST 4: OCEAN CURRENT ALIGNMENT")
print("="*70)
print("Using known major ocean current systems with approximate bearings.")

# Major ocean currents with approximate bounding boxes and directions
MAJOR_CURRENTS = [
    # (name, lat_min, lat_max, lon_min, lon_max, bearing_deg, speed_note)
    ("Gulf Stream", 25, 45, -80, -40, 45, "strong NE"),
    ("North Atlantic Current", 40, 60, -40, -10, 60, "NE"),
    ("Canary Current", 15, 35, -25, -10, 200, "SSW"),
    ("North Equatorial Current", 5, 20, -80, -15, 270, "W"),
    ("South Equatorial Current", -10, 0, -40, 10, 270, "W"),
    ("Brazil Current", -35, -10, -55, -35, 200, "SSW"),
    ("Benguela Current", -35, -15, 5, 20, 350, "N"),
    ("Agulhas Current", -40, -25, 25, 45, 240, "WSW"),
    ("Somali Current", -5, 15, 45, 55, 30, "NE (summer)"),
    ("Humboldt Current", -45, -5, -85, -70, 350, "N"),
    ("South Pacific Current", -55, -40, -120, -70, 90, "E"),
    ("Kuroshio Current", 20, 40, 125, 150, 45, "NE"),
    ("North Pacific Current", 35, 50, 150, -140, 90, "E"),
    ("California Current", 25, 45, -130, -115, 170, "S"),
    ("Antarctic Circumpolar", -65, -50, -180, 180, 90, "E"),
    ("East Australian Current", -40, -20, 150, 165, 200, "SSW"),
]


def circle_bearing_at(pts, idx):
    """Approximate bearing of the circle at a given point."""
    i0 = (idx - 1) % len(pts)
    i1 = (idx + 1) % len(pts)
    lat0, lon0 = pts[i0]
    lat1, lon1 = pts[i1]
    result = geod.Inverse(lat0, lon0, lat1, lon1)
    return result['azi1'] % 360


def angle_diff(a, b):
    """Minimum angle between two bearings (0-180)."""
    d = abs(a - b) % 360
    return min(d, 360 - d)


# Find ocean-crossing segments of Alison circle (elevation ≤ 0)
ocean_segments = []
for i, idx in enumerate(topo_indices):
    if valid_e[i] and alison_elev[i] <= 0:
        lat, lon = alison_pts[idx]
        bearing = circle_bearing_at(alison_pts, idx)

        # Check which currents this point is in
        for cname, clat_min, clat_max, clon_min, clon_max, cbearing, cnote in MAJOR_CURRENTS:
            # Handle wrap-around for Pacific
            if clon_min > clon_max:
                in_box = (clat_min <= lat <= clat_max and (lon >= clon_min or lon <= clon_max))
            else:
                in_box = (clat_min <= lat <= clat_max and clon_min <= lon <= clon_max)

            if in_box:
                alignment = angle_diff(bearing, cbearing)
                ocean_segments.append({
                    "lat": round(lat, 2), "lon": round(lon, 2),
                    "circle_bearing": round(bearing, 1),
                    "current": cname,
                    "current_bearing": cbearing,
                    "alignment_angle": round(alignment, 1),
                    "aligned": alignment < 30,
                    "note": cnote
                })

n_aligned = sum(1 for s in ocean_segments if s["aligned"])
n_total = len(ocean_segments)
print(f"\nOcean-current intersections: {n_total}")
print(f"Aligned (< 30°): {n_aligned} ({n_aligned/n_total*100:.0f}%)" if n_total > 0 else "No ocean crossings found")

print(f"\n{'Current':<25} {'Lat':>6} {'Lon':>7} {'Circle°':>8} {'Current°':>9} {'Angle':>6} {'Aligned'}")
print("-" * 75)
for s in ocean_segments:
    mark = "YES" if s["aligned"] else "no"
    print(f"  {s['current']:<23} {s['lat']:>6.1f} {s['lon']:>7.1f} "
          f"{s['circle_bearing']:>8.1f} {s['current_bearing']:>9} "
          f"{s['alignment_angle']:>6.1f} {mark}")

# Compare to random circles
if n_total > 0:
    random_alignment_rates = []
    for i in range(min(N_RANDOM, 50)):
        plat, plon = random_poles[i]
        rpts = great_circle_points(plat, plon, N_POINTS)
        # Use geoid grid to approximate ocean (geoid doesn't directly give this,
        # but we can check if elevation would be ocean from the sampled elevations)
        # Instead, use a simple lat/lon ocean mask based on known geography
        r_aligned = 0
        r_total = 0
        for idx in topo_indices:
            rlat, rlon = rpts[idx]
            rbearing = circle_bearing_at(rpts, idx)
            for cname, clat_min, clat_max, clon_min, clon_max, cbearing, cnote in MAJOR_CURRENTS:
                if clon_min > clon_max:
                    in_box = (clat_min <= rlat <= clat_max and (rlon >= clon_min or rlon <= clon_max))
                else:
                    in_box = (clat_min <= rlat <= clat_max and clon_min <= rlon <= clon_max)
                if in_box:
                    r_total += 1
                    if angle_diff(rbearing, cbearing) < 30:
                        r_aligned += 1
        if r_total > 0:
            random_alignment_rates.append(r_aligned / r_total)

    if random_alignment_rates:
        alison_rate = n_aligned / n_total
        ocean_z = z_of(alison_rate, random_alignment_rates)
        print(f"\nAlignment rate Z-score: {ocean_z:+.2f} (Alison={alison_rate*100:.0f}%, "
              f"random={np.mean(random_alignment_rates)*100:.0f}±{np.std(random_alignment_rates)*100:.0f}%)")
    else:
        ocean_z = 0
else:
    ocean_z = 0
    random_alignment_rates = []

test4_result = {
    "n_ocean_current_intersections": n_total,
    "n_aligned_lt30deg": n_aligned,
    "alignment_rate": round(n_aligned / n_total, 3) if n_total > 0 else 0,
    "z_score_alignment": round(ocean_z, 2),
    "segments": ocean_segments,
    "random_mean_rate": round(float(np.mean(random_alignment_rates)), 3) if random_alignment_rates else None,
    "note": "Uses known major current bounding boxes and directions; not a full ocean current model"
}


# ============================================================
# OVERALL SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("OVERALL SUMMARY")
print("="*70)

summary_rows = [
    ("Geoid mean height", geoid_z_mean),
    ("Geoid std", geoid_z_std),
    ("Geoid gradient", geoid_z_gradient),
    ("Magnetic mean intensity", mag_z_mean),
    ("Magnetic std intensity", mag_z_std),
    ("Topographic mean elev", topo_z_mean),
    ("Topographic std elev", topo_z_std),
    ("Ocean fraction", topo_z_ocean),
    ("Low elev fraction", topo_z_low),
    ("Ocean current alignment", ocean_z),
]

print(f"\n{'Property':<30} {'Z-score':>10} {'Significant'}")
print("-" * 55)
for name, z in summary_rows:
    sig = "YES" if abs(z) > 2 else "no"
    print(f"  {name:<28} {z:>+10.2f}  {sig}")

sig_count = sum(1 for _, z in summary_rows if abs(z) > 2)
print(f"\nSignificant correlations (|Z| > 2): {sig_count} out of {len(summary_rows)}")

if sig_count == 0:
    overall = "NO GEOPHYSICAL CORRELATION — the circle does not track any measured Earth property"
elif sig_count <= 2:
    overall = "WEAK/MARGINAL correlation with some properties"
else:
    overall = "SIGNIFICANT geophysical correlation detected"
print(f"Overall verdict: {overall}")


# ============================================================
# SAVE
# ============================================================
output = {
    "meta": {
        "name": "Geophysical Property Scan",
        "description": "Tests whether the Alison Great Circle correlates with geoid, magnetic, topographic, or ocean current properties",
        "n_random_circles": N_RANDOM,
        "n_points_per_circle": N_POINTS,
        "pole": {"lat": POLE_LAT, "lon": POLE_LNG}
    },
    "test1_geoid_anomaly": test1_result,
    "test2_magnetic_field": test2_result,
    "test3_topographic_consistency": test3_result,
    "test4_ocean_current_alignment": test4_result,
    "summary": {
        "z_scores": {name: round(z, 2) for name, z in summary_rows},
        "n_significant": sig_count,
        "verdict": overall
    }
}

outpath = os.path.join(BASE, "results/geophysical_property_scan.json")
with open(outpath, "w") as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nResults saved to {outpath}")
