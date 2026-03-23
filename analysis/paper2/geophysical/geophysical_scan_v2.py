#!/usr/bin/env python3
"""
Geophysical Property Scan v2
==============================
Test whether the Alison Great Circle correlates with measurable geophysical
properties of the Earth.

Data sources (all publicly available):
- EGM96 5-arcminute geoid grid (NGA) — adequate proxy for EGM2008 at circle scale
- EMAG2v3 crustal magnetic anomaly grid, 2-arcminute (NOAA NCEI)
- ETOPO1 global relief model, 1-arcminute (NOAA NCEI)
- Ocean currents: major current systems with directions (analytical model)

Tests:
1. Geoid anomaly — does the circle follow a geoid ridge, trough, or gradient?
2. Crustal magnetic anomaly — is there a magnetic signature?
3. Topographic consistency — unusually consistent elevation or coastal zones?
4. Ocean current alignment — does the circle align with navigable currents?
"""

import json, math, os, sys, time, struct, zipfile
import numpy as np
import netCDF4 as nc
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
# UTILITY: Generate great circle points
# ============================================================
def great_circle_points(pole_lat, pole_lon, n_points=360):
    """Generate n_points evenly spaced along a great circle defined by its pole."""
    points = []
    for i in range(n_points):
        bearing = i * (360.0 / n_points)
        result = geod.Direct(pole_lat, pole_lon, bearing, EARTH_R * 1000 * math.pi / 2)
        points.append((result['lat2'], result['lon2']))
    return np.array(points)


def random_pole():
    """Generate a random pole (uniformly distributed on sphere)."""
    z = np.random.uniform(-1, 1)
    theta = np.random.uniform(0, 2 * np.pi)
    return np.degrees(np.arcsin(z)), np.degrees(theta) - 180


def z_of(val, arr):
    s = np.std(arr)
    return (val - np.mean(arr)) / s if s > 0 else 0.0


def percentile_of(val, arr):
    return float(np.mean(np.array(arr) < val) * 100)


# Generate all circles upfront
print("Generating circle points...")
alison_pts = great_circle_points(POLE_LAT, POLE_LNG, N_POINTS)
print(f"Alison circle: {len(alison_pts)} points, lat [{alison_pts[:,0].min():.1f}, {alison_pts[:,0].max():.1f}]")

random_poles_list = [random_pole() for _ in range(N_RANDOM)]
print(f"Generated {N_RANDOM} random poles")

# Pre-generate random circle points for reuse across tests
print("Pre-generating random circle points...")
t0 = time.time()
random_circle_pts = [great_circle_points(plat, plon, N_POINTS) for plat, plon in random_poles_list]
print(f"  Done [{time.time()-t0:.1f}s]")


# ============================================================
# TEST 1: GEOID ANOMALY (EGM96)
# ============================================================
print(f"\n{'='*70}")
print("TEST 1: GEOID ANOMALY (EGM96 5-arcminute)")
print("="*70)
print("Note: EGM96 used (already available). At great-circle scale,")
print("EGM96 and EGM2008 differ by <1m RMS — negligible for this test.")

pgm_path = os.path.join(BASE, "data/geophysical/geoids/egm96-5.pgm")

with open(pgm_path, 'rb') as f:
    offset, scale = -108.0, 0.003
    while True:
        line = f.readline()
        if line.startswith(b'#'):
            text = line.decode().strip()
            if 'Offset' in text:
                offset = float(text.split()[-1])
            if 'Scale' in text:
                scale = float(text.split()[-1])
        elif not line.startswith(b'P'):
            parts = line.strip().split()
            g_width, g_height = int(parts[0]), int(parts[1])
            maxval = int(f.readline().strip())
            break
    raw = f.read(g_width * g_height * 2)
    geoid_grid = offset + scale * np.frombuffer(raw, dtype='>u2').reshape(g_height, g_width).astype(float)

print(f"EGM96 grid: {g_width}x{g_height}, range [{geoid_grid.min():.1f}, {geoid_grid.max():.1f}] m")
GEOID_RES = 12  # cells per degree (5 arcminutes)


def sample_geoid(lat, lon):
    lon = lon % 360
    y = (90 - lat) * GEOID_RES
    x = lon * GEOID_RES
    y0, x0 = int(y), int(x)
    y1 = min(y0 + 1, g_height - 1)
    x1 = (x0 + 1) % g_width
    fy, fx = y - y0, x - x0
    return float(
        geoid_grid[y0, x0] * (1-fx)*(1-fy) +
        geoid_grid[y0, x1] * fx*(1-fy) +
        geoid_grid[y1, x0] * (1-fx)*fy +
        geoid_grid[y1, x1] * fx*fy
    )


def geoid_stats(pts):
    vals = np.array([sample_geoid(lat, lon) for lat, lon in pts])
    return {
        "mean": float(np.mean(vals)),
        "std": float(np.std(vals)),
        "range": float(np.ptp(vals)),
        "gradient_mean": float(np.mean(np.abs(np.diff(vals)))),
    }


t0 = time.time()
alison_geoid = geoid_stats(alison_pts)
print(f"Alison geoid: mean={alison_geoid['mean']:.2f}m, std={alison_geoid['std']:.2f}m, "
      f"range={alison_geoid['range']:.2f}m, gradient={alison_geoid['gradient_mean']:.3f}m/deg")

random_geoid = {"means": [], "stds": [], "ranges": [], "gradients": []}
for i, rpts in enumerate(random_circle_pts):
    stats = geoid_stats(rpts)
    random_geoid["means"].append(stats["mean"])
    random_geoid["stds"].append(stats["std"])
    random_geoid["ranges"].append(stats["range"])
    random_geoid["gradients"].append(stats["gradient_mean"])
    if (i + 1) % 25 == 0:
        print(f"  {i+1}/{N_RANDOM} random circles [{time.time()-t0:.1f}s]")

geoid_z = {
    "mean_height": z_of(alison_geoid["mean"], random_geoid["means"]),
    "std_height": z_of(alison_geoid["std"], random_geoid["stds"]),
    "range": z_of(alison_geoid["range"], random_geoid["ranges"]),
    "gradient": z_of(alison_geoid["gradient_mean"], random_geoid["gradients"]),
}

for k, v in geoid_z.items():
    print(f"  Z({k}) = {v:+.2f}")

geoid_pct = percentile_of(alison_geoid["mean"], random_geoid["means"])
print(f"  Percentile (mean height): {geoid_pct:.1f}%")

if abs(geoid_z["mean_height"]) > 2:
    geoid_interp = f"Circle follows geoid {'ridge' if geoid_z['mean_height'] > 0 else 'trough'} (Z={geoid_z['mean_height']:+.2f})"
elif abs(geoid_z["gradient"]) > 2:
    geoid_interp = f"Circle follows unusual geoid gradient (Z={geoid_z['gradient']:+.2f})"
else:
    geoid_interp = "No significant geoid correlation"

test1_result = {
    "data_source": "EGM96 5-arcminute geoid grid (NGA)",
    "note": "EGM96 and EGM2008 differ by <1m RMS at great-circle wavelengths",
    "alison": {k: round(v, 3) for k, v in alison_geoid.items()},
    "random_100": {
        "mean_of_means": round(float(np.mean(random_geoid["means"])), 3),
        "std_of_means": round(float(np.std(random_geoid["means"])), 3),
        "mean_of_stds": round(float(np.mean(random_geoid["stds"])), 3),
        "mean_of_gradients": round(float(np.mean(random_geoid["gradients"])), 3),
    },
    "z_scores": {k: round(v, 2) for k, v in geoid_z.items()},
    "percentile_mean": round(geoid_pct, 1),
    "interpretation": geoid_interp,
}
print(f"  => {geoid_interp}")


# ============================================================
# TEST 2: CRUSTAL MAGNETIC ANOMALY (EMAG2v3)
# ============================================================
print(f"\n{'='*70}")
print("TEST 2: CRUSTAL MAGNETIC ANOMALY (EMAG2v3)")
print("="*70)

emag_zip = os.path.join(BASE, "data/geophysical/emag2/EMAG2_V3_20170530.zip")
emag_bin = os.path.join(BASE, "data/geophysical/emag2/emag2v3_grid.npy")

# EMAG2v3 grid parameters (from data documentation)
EMAG_LON_MIN, EMAG_LON_MAX = 0.0333333, 360.0  # roughly 0 to 360
EMAG_LAT_MIN, EMAG_LAT_MAX = -89.9833, 89.9833
EMAG_RES = 30  # cells per degree (2 arcminutes)
EMAG_NCOLS = 10800  # 360 * 30
EMAG_NROWS = 5400   # 180 * 30 (but actual coverage is -72 to 72 for valid data)

# Convert CSV to numpy binary grid (one-time operation)
if not os.path.exists(emag_bin):
    print("Converting EMAG2v3 CSV to binary grid (one-time, ~5 min)...")
    # The CSV columns are: id, ?, lon, lat, ?, anomaly_nT, ?, ?
    # Grid: 10800 cols (lon) x 5400 rows (lat), from -90 to +90, 0 to 360
    grid = np.full((EMAG_NROWS, EMAG_NCOLS), np.nan, dtype=np.float32)

    t0 = time.time()
    row_count = 0
    with zipfile.ZipFile(emag_zip) as z:
        with z.open('EMAG2_V3_20170530.csv') as f:
            for line in f:
                parts = line.decode().strip().split(',')
                if len(parts) >= 6:
                    try:
                        lon = float(parts[2])
                        lat = float(parts[3])
                        anomaly = float(parts[5])

                        # Convert to grid indices
                        col = int(round((lon % 360) * EMAG_RES)) % EMAG_NCOLS
                        row = int(round((lat + 90) * EMAG_RES)) % EMAG_NROWS

                        if -90000 < anomaly < 90000:  # skip sentinel values (99999, -99999)
                            grid[row, col] = anomaly
                    except (ValueError, IndexError):
                        pass

                row_count += 1
                if row_count % 5_000_000 == 0:
                    elapsed = time.time() - t0
                    pct_valid = np.sum(~np.isnan(grid)) / grid.size * 100
                    print(f"  {row_count/1e6:.0f}M rows, {pct_valid:.1f}% valid [{elapsed:.0f}s]")

    np.save(emag_bin, grid)
    n_valid = np.sum(~np.isnan(grid))
    print(f"  Grid saved: {n_valid}/{grid.size} valid values ({n_valid/grid.size*100:.1f}%)")
    print(f"  Total time: {time.time()-t0:.0f}s")
else:
    print("Loading pre-computed EMAG2v3 binary grid...")

emag_grid = np.load(emag_bin)
n_valid = np.sum(~np.isnan(emag_grid))
print(f"EMAG2v3 grid: {emag_grid.shape}, {n_valid} valid values ({n_valid/emag_grid.size*100:.1f}%)")


def sample_emag(lat, lon):
    """Sample EMAG2v3 crustal magnetic anomaly at (lat, lon) with bilinear interpolation."""
    lon = lon % 360
    y = (lat + 90) * EMAG_RES
    x = lon * EMAG_RES

    y0, x0 = int(y), int(x)
    y1 = min(y0 + 1, EMAG_NROWS - 1)
    x1 = (x0 + 1) % EMAG_NCOLS
    fy, fx = y - y0, x - x0

    vals = [emag_grid[y0, x0], emag_grid[y0, x1],
            emag_grid[y1, x0], emag_grid[y1, x1]]

    # If any neighbor is NaN, use nearest valid
    if any(np.isnan(v) for v in vals):
        valid = [v for v in vals if not np.isnan(v)]
        return float(np.mean(valid)) if valid else np.nan

    return float(vals[0]*(1-fx)*(1-fy) + vals[1]*fx*(1-fy) +
                 vals[2]*(1-fx)*fy + vals[3]*fx*fy)


def magnetic_stats(pts):
    vals = np.array([sample_emag(lat, lon) for lat, lon in pts])
    valid = vals[~np.isnan(vals)]
    if len(valid) < 10:
        return None
    return {
        "mean": float(np.mean(valid)),
        "std": float(np.std(valid)),
        "abs_mean": float(np.mean(np.abs(valid))),
        "range": float(np.ptp(valid)),
        "n_valid": int(len(valid)),
    }


t0 = time.time()
alison_mag = magnetic_stats(alison_pts)
if alison_mag:
    print(f"Alison magnetic: mean={alison_mag['mean']:.1f}nT, std={alison_mag['std']:.1f}nT, "
          f"|mean|={alison_mag['abs_mean']:.1f}nT, valid={alison_mag['n_valid']}/{N_POINTS}")

    random_mag = {"means": [], "stds": [], "abs_means": [], "ranges": []}
    for i, rpts in enumerate(random_circle_pts):
        stats = magnetic_stats(rpts)
        if stats:
            random_mag["means"].append(stats["mean"])
            random_mag["stds"].append(stats["std"])
            random_mag["abs_means"].append(stats["abs_mean"])
            random_mag["ranges"].append(stats["range"])
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{N_RANDOM} random circles [{time.time()-t0:.1f}s]")

    mag_z = {
        "mean_anomaly": z_of(alison_mag["mean"], random_mag["means"]),
        "std_anomaly": z_of(alison_mag["std"], random_mag["stds"]),
        "abs_mean_anomaly": z_of(alison_mag["abs_mean"], random_mag["abs_means"]),
    }

    for k, v in mag_z.items():
        print(f"  Z({k}) = {v:+.2f}")

    if abs(mag_z["mean_anomaly"]) > 2:
        mag_interp = f"Circle follows {'positive' if mag_z['mean_anomaly'] > 0 else 'negative'} magnetic anomaly (Z={mag_z['mean_anomaly']:+.2f})"
    elif abs(mag_z["std_anomaly"]) > 2:
        mag_interp = f"Circle has unusually {'high' if mag_z['std_anomaly'] > 0 else 'low'} magnetic variability (Z={mag_z['std_anomaly']:+.2f})"
    else:
        mag_interp = "No significant crustal magnetic anomaly correlation"

    test2_result = {
        "data_source": "EMAG2v3 (2-arcminute, NOAA NCEI, doi:10.7289/V5H70CVX)",
        "n_random_circles_valid": len(random_mag["means"]),
        "alison": {k: round(v, 2) for k, v in alison_mag.items()},
        "random": {
            "mean_of_means": round(float(np.mean(random_mag["means"])), 2),
            "std_of_means": round(float(np.std(random_mag["means"])), 2),
            "mean_of_stds": round(float(np.mean(random_mag["stds"])), 2),
            "mean_of_abs_means": round(float(np.mean(random_mag["abs_means"])), 2),
        },
        "z_scores": {k: round(v, 2) for k, v in mag_z.items()},
        "interpretation": mag_interp,
    }
else:
    mag_z = {"mean_anomaly": 0, "std_anomaly": 0, "abs_mean_anomaly": 0}
    mag_interp = "Insufficient magnetic data"
    test2_result = {
        "data_source": "EMAG2v3 (2-arcminute, NOAA NCEI)",
        "error": "Insufficient valid data points",
    }

print(f"  => {mag_interp}")


# ============================================================
# TEST 3: TOPOGRAPHIC CONSISTENCY (ETOPO1)
# ============================================================
print(f"\n{'='*70}")
print("TEST 3: TOPOGRAPHIC CONSISTENCY (ETOPO1)")
print("="*70)

etopo_path = os.path.join(BASE, "data/geophysical/etopo/ETOPO1_Ice_g_gmt4.grd")
etopo_ds = nc.Dataset(etopo_path)
etopo_x = etopo_ds.variables['x'][:]  # -180 to 180
etopo_y = etopo_ds.variables['y'][:]  # -90 to 90
# Don't load entire z array (~900MB) — sample on demand
print(f"ETOPO1 grid: {len(etopo_x)}x{len(etopo_y)}, "
      f"lon [{etopo_x[0]:.1f}, {etopo_x[-1]:.1f}], lat [{etopo_y[0]:.1f}, {etopo_y[-1]:.1f}]")

ETOPO_RES = 60  # cells per degree (1 arcminute)
ETOPO_LON_OFFSET = 180.0  # etopo_x starts at -180


def sample_etopo(lat, lon):
    """Sample ETOPO1 elevation at (lat, lon) via nearest neighbor (1-arcminute is fine)."""
    lon_idx = int(round((lon + 180.0) * ETOPO_RES))
    lat_idx = int(round((lat + 90.0) * ETOPO_RES))
    lon_idx = max(0, min(lon_idx, len(etopo_x) - 1))
    lat_idx = max(0, min(lat_idx, len(etopo_y) - 1))
    return float(etopo_ds.variables['z'][lat_idx, lon_idx])


def batch_sample_etopo(pts):
    """Sample elevations for an array of points efficiently."""
    elevs = np.zeros(len(pts))
    for i, (lat, lon) in enumerate(pts):
        elevs[i] = sample_etopo(lat, lon)
    return elevs


def topo_stats(pts):
    elevs = batch_sample_etopo(pts)
    return {
        "mean_elevation": float(np.mean(elevs)),
        "std_elevation": float(np.std(elevs)),
        "median_elevation": float(np.median(elevs)),
        "ocean_fraction": float(np.mean(elevs <= 0)),
        "low_elevation_fraction": float(np.mean(elevs <= 200)),
        "coastal_fraction": float(np.mean((elevs > -200) & (elevs <= 200))),
        "n_sign_changes": int(np.sum(np.diff(np.sign(elevs)) != 0)),
    }


t0 = time.time()
alison_topo = topo_stats(alison_pts)
print(f"Alison topo: mean={alison_topo['mean_elevation']:.0f}m, std={alison_topo['std_elevation']:.0f}m, "
      f"ocean={alison_topo['ocean_fraction']*100:.1f}%, "
      f"low_elev={alison_topo['low_elevation_fraction']*100:.1f}%, "
      f"coastal={alison_topo['coastal_fraction']*100:.1f}%, "
      f"sign_changes={alison_topo['n_sign_changes']}")

random_topo = {k: [] for k in alison_topo.keys()}
for i, rpts in enumerate(random_circle_pts):
    stats = topo_stats(rpts)
    for k in random_topo:
        random_topo[k].append(stats[k])
    if (i + 1) % 25 == 0:
        print(f"  {i+1}/{N_RANDOM} random circles [{time.time()-t0:.1f}s]")

topo_z = {}
for k in ["mean_elevation", "std_elevation", "ocean_fraction",
          "low_elevation_fraction", "coastal_fraction", "n_sign_changes"]:
    topo_z[k] = z_of(alison_topo[k], random_topo[k])

for k, v in topo_z.items():
    print(f"  Z({k}) = {v:+.2f}")

# Interpretation
topo_findings = []
if abs(topo_z["ocean_fraction"]) > 2:
    topo_findings.append(f"{'more' if topo_z['ocean_fraction'] > 0 else 'less'} ocean than random (Z={topo_z['ocean_fraction']:+.2f})")
if abs(topo_z["std_elevation"]) > 2:
    topo_findings.append(f"{'more' if topo_z['std_elevation'] > 0 else 'less'} elevation variance (Z={topo_z['std_elevation']:+.2f})")
if abs(topo_z["low_elevation_fraction"]) > 2:
    topo_findings.append(f"{'more' if topo_z['low_elevation_fraction'] > 0 else 'fewer'} low-elevation zones (Z={topo_z['low_elevation_fraction']:+.2f})")
if abs(topo_z["coastal_fraction"]) > 2:
    topo_findings.append(f"{'more' if topo_z['coastal_fraction'] > 0 else 'fewer'} coastal zones (Z={topo_z['coastal_fraction']:+.2f})")
if abs(topo_z["n_sign_changes"]) > 2:
    topo_findings.append(f"{'more' if topo_z['n_sign_changes'] > 0 else 'fewer'} land-ocean crossings (Z={topo_z['n_sign_changes']:+.2f})")

topo_interp = "; ".join(topo_findings) if topo_findings else "No significant topographic anomaly"

test3_result = {
    "data_source": "ETOPO1 Ice Surface, 1-arcminute (NOAA NCEI)",
    "n_points_per_circle": N_POINTS,
    "n_random_circles": N_RANDOM,
    "alison": {k: round(v, 3) for k, v in alison_topo.items()},
    "random_100": {
        "mean_of_" + k: round(float(np.mean(random_topo[k])), 3)
        for k in alison_topo.keys()
    },
    "z_scores": {k: round(v, 2) for k, v in topo_z.items()},
    "percentile_ocean_fraction": round(percentile_of(alison_topo["ocean_fraction"], random_topo["ocean_fraction"]), 1),
    "interpretation": topo_interp,
}
print(f"  => {topo_interp}")


# ============================================================
# TEST 4: OCEAN CURRENT ALIGNMENT
# ============================================================
print(f"\n{'='*70}")
print("TEST 4: OCEAN CURRENT ALIGNMENT")
print("="*70)
print("Using known major ocean current systems with bearings.")
print("For each ocean-crossing point, checks alignment with nearest current.")

# Major ocean currents: (name, lat_min, lat_max, lon_min, lon_max, bearing_deg)
# lon_min > lon_max means the box wraps around the dateline
MAJOR_CURRENTS = [
    ("Gulf Stream", 25, 45, -80, -40, 45),
    ("North Atlantic Current", 40, 60, -40, -10, 60),
    ("Canary Current", 15, 35, -25, -10, 200),
    ("North Equatorial Current (Atl)", 5, 20, -80, -15, 270),
    ("South Equatorial Current (Atl)", -10, 0, -40, 10, 270),
    ("Brazil Current", -35, -10, -55, -35, 200),
    ("Benguela Current", -35, -15, 5, 20, 350),
    ("Agulhas Current", -40, -25, 25, 45, 240),
    ("Somali Current", -5, 15, 45, 55, 30),
    ("Humboldt Current", -45, -5, -85, -70, 350),
    ("South Pacific Current", -55, -40, -120, -70, 90),
    ("Kuroshio Current", 20, 40, 125, 150, 45),
    ("North Pacific Current", 35, 50, 150, 220, 90),  # wraps dateline
    ("California Current", 25, 45, -130, -115, 170),
    ("Antarctic Circumpolar", -65, -50, -180, 180, 90),
    ("East Australian Current", -40, -20, 150, 165, 200),
    ("Mozambique Current", -25, -10, 35, 45, 200),
    ("North Equatorial Current (Pac)", 5, 20, 120, 280, 270),  # wraps dateline
    ("South Equatorial Current (Pac)", -10, 0, 150, 280, 270),  # wraps dateline
    ("Indian Ocean Gyre (N)", 0, 15, 55, 100, 90),
    ("Indian Ocean Gyre (S)", -30, -10, 40, 110, 270),
]


def circle_bearing_at(pts, idx):
    i0 = (idx - 1) % len(pts)
    i1 = (idx + 1) % len(pts)
    result = geod.Inverse(pts[i0][0], pts[i0][1], pts[i1][0], pts[i1][1])
    return result['azi1'] % 360


def angle_diff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


def in_current_box(lat, lon, clat_min, clat_max, clon_min, clon_max):
    if not (clat_min <= lat <= clat_max):
        return False
    # Normalize lon to [-180, 180] and also check [0, 360] representation
    lon180 = ((lon + 180) % 360) - 180
    if clon_min <= clon_max:
        return clon_min <= lon180 <= clon_max or clon_min <= (lon180 + 360) <= clon_max
    else:
        # Wraps around — e.g., 150 to 220 (= -140)
        lon360 = lon % 360
        clon_max_360 = clon_max % 360
        return lon360 >= clon_min or lon360 <= clon_max_360


def ocean_current_alignment(pts, elevations):
    """Compute ocean current alignment for a circle given its elevations."""
    segments = []
    for i in range(len(pts)):
        if elevations[i] > 0:
            continue  # Skip land
        lat, lon = pts[i]
        bearing = circle_bearing_at(pts, i)

        for cname, clat_min, clat_max, clon_min, clon_max, cbearing in MAJOR_CURRENTS:
            if in_current_box(lat, lon, clat_min, clat_max, clon_min, clon_max):
                alignment = angle_diff(bearing, cbearing)
                segments.append({
                    "lat": round(lat, 2), "lon": round(lon, 2),
                    "circle_bearing": round(bearing, 1),
                    "current": cname,
                    "current_bearing": cbearing,
                    "alignment_angle": round(alignment, 1),
                    "aligned": alignment < 30,
                })

    if not segments:
        return 0, 0, segments
    n_aligned = sum(1 for s in segments if s["aligned"])
    return n_aligned, len(segments), segments


# Get Alison circle elevations (already computed above or sample now)
print("Computing ocean current alignment...")
t0 = time.time()
alison_elevs = batch_sample_etopo(alison_pts)
a_aligned, a_total, a_segments = ocean_current_alignment(alison_pts, alison_elevs)
a_rate = a_aligned / a_total if a_total > 0 else 0
print(f"Alison: {a_aligned}/{a_total} aligned (<30°), rate={a_rate*100:.1f}%")

# Deduplicate segments by current name for display
seen = set()
unique_segments = []
for s in a_segments:
    key = s["current"]
    if key not in seen:
        seen.add(key)
        unique_segments.append(s)

print(f"\n{'Current':<30} {'Lat':>6} {'Lon':>7} {'Circle°':>8} {'Current°':>9} {'Angle':>6} {'Aligned'}")
print("-" * 80)
for s in unique_segments:
    mark = "YES" if s["aligned"] else "no"
    print(f"  {s['current']:<28} {s['lat']:>6.1f} {s['lon']:>7.1f} "
          f"{s['circle_bearing']:>8.1f} {s['current_bearing']:>9} "
          f"{s['alignment_angle']:>6.1f} {mark}")

# Random circles
random_alignment_rates = []
for i, rpts in enumerate(random_circle_pts):
    r_elevs = batch_sample_etopo(rpts)
    r_aligned, r_total, _ = ocean_current_alignment(rpts, r_elevs)
    if r_total > 0:
        random_alignment_rates.append(r_aligned / r_total)
    if (i + 1) % 25 == 0:
        print(f"  {i+1}/{N_RANDOM} random circles [{time.time()-t0:.1f}s]")

if random_alignment_rates and a_total > 0:
    ocean_z = z_of(a_rate, random_alignment_rates)
    print(f"\nAlignment rate Z-score: {ocean_z:+.2f} (Alison={a_rate*100:.0f}%, "
          f"random={np.mean(random_alignment_rates)*100:.0f}±{np.std(random_alignment_rates)*100:.0f}%)")
else:
    ocean_z = 0.0

if abs(ocean_z) > 2:
    ocean_interp = f"Circle is unusually {'well' if ocean_z > 0 else 'poorly'}-aligned with ocean currents (Z={ocean_z:+.2f})"
else:
    ocean_interp = "No significant ocean current alignment"

test4_result = {
    "data_source": "Analytical model of 21 major ocean current systems with ETOPO1 ocean mask",
    "n_ocean_current_intersections": a_total,
    "n_aligned_lt30deg": a_aligned,
    "alignment_rate": round(a_rate, 3) if a_total > 0 else 0,
    "z_score_alignment": round(ocean_z, 2),
    "unique_currents_crossed": [s["current"] for s in unique_segments],
    "segments_summary": unique_segments,
    "random_100": {
        "n_valid": len(random_alignment_rates),
        "mean_rate": round(float(np.mean(random_alignment_rates)), 3) if random_alignment_rates else None,
        "std_rate": round(float(np.std(random_alignment_rates)), 3) if random_alignment_rates else None,
    },
    "interpretation": ocean_interp,
}
print(f"  => {ocean_interp}")


# Close ETOPO dataset
etopo_ds.close()


# ============================================================
# OVERALL SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("OVERALL SUMMARY")
print("="*70)

summary_rows = [
    ("Geoid mean height", geoid_z["mean_height"]),
    ("Geoid std", geoid_z["std_height"]),
    ("Geoid gradient", geoid_z["gradient"]),
    ("Magnetic mean anomaly", mag_z.get("mean_anomaly", 0)),
    ("Magnetic std anomaly", mag_z.get("std_anomaly", 0)),
    ("Magnetic |mean| anomaly", mag_z.get("abs_mean_anomaly", 0)),
    ("Topographic mean elev", topo_z["mean_elevation"]),
    ("Topographic std elev", topo_z["std_elevation"]),
    ("Ocean fraction", topo_z["ocean_fraction"]),
    ("Low elev fraction", topo_z["low_elevation_fraction"]),
    ("Coastal fraction", topo_z["coastal_fraction"]),
    ("Land-ocean crossings", topo_z["n_sign_changes"]),
    ("Ocean current alignment", ocean_z),
]

print(f"\n{'Property':<30} {'Z-score':>10} {'Significant'}")
print("-" * 55)
for name, z in summary_rows:
    sig = "YES" if abs(z) > 2 else "no"
    print(f"  {name:<28} {z:>+10.2f}  {sig}")

sig_count = sum(1 for _, z in summary_rows if abs(z) > 2)
total_tests = len(summary_rows)
print(f"\nSignificant (|Z| > 2): {sig_count}/{total_tests}")

# Bonferroni correction note
bonf_threshold = 0.05 / total_tests
bonf_z = float(abs(np.percentile(np.random.standard_normal(100000), bonf_threshold/2 * 100)))
# Approximate: for 13 tests, Bonferroni threshold is |Z| > ~2.94
bonf_sig = sum(1 for _, z in summary_rows if abs(z) > 2.94)

if sig_count == 0:
    verdict = "NO geophysical correlation — circle does not track any Earth property"
elif bonf_sig > 0:
    verdict = f"SIGNIFICANT after Bonferroni correction ({bonf_sig}/{total_tests} survive)"
elif sig_count >= 3:
    verdict = f"MODERATE correlation ({sig_count}/{total_tests} nominally significant, none survive Bonferroni)"
else:
    verdict = f"WEAK/MARGINAL ({sig_count}/{total_tests} nominally significant, likely chance)"

print(f"Bonferroni-corrected (|Z| > 2.94 for {total_tests} tests): {bonf_sig}")
print(f"Verdict: {verdict}")


# ============================================================
# SAVE
# ============================================================
output = {
    "meta": {
        "name": "Geophysical Property Scan v2",
        "description": "Tests whether the Alison Great Circle correlates with geoid, crustal magnetic, topographic, or ocean current properties",
        "version": 2,
        "n_random_circles": N_RANDOM,
        "n_points_per_circle": N_POINTS,
        "pole": {"lat": POLE_LAT, "lon": POLE_LNG},
        "data_sources": {
            "geoid": "EGM96 5-arcminute (NGA)",
            "magnetic": "EMAG2v3 2-arcminute crustal anomaly (NOAA NCEI, doi:10.7289/V5H70CVX)",
            "topographic": "ETOPO1 Ice Surface 1-arcminute (NOAA NCEI)",
            "ocean_currents": "Analytical model of 21 major current systems",
        },
    },
    "test1_geoid_anomaly": test1_result,
    "test2_crustal_magnetic_anomaly": test2_result,
    "test3_topographic_consistency": test3_result,
    "test4_ocean_current_alignment": test4_result,
    "summary": {
        "z_scores": {name: round(z, 2) for name, z in summary_rows},
        "n_nominally_significant": sig_count,
        "n_bonferroni_significant": bonf_sig,
        "bonferroni_threshold_z": 2.94,
        "verdict": verdict,
    }
}

outpath = os.path.join(BASE, "results/geophysical_property_scan.json")
with open(outpath, "w") as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nResults saved to {outpath}")
