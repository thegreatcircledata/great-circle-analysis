#!/usr/bin/env python3
"""
Directive 2: Anti-Divergence Circles
=====================================
Find the great circle that maximizes SETTLEMENT clustering while minimizing
MONUMENT clustering — the inverse of the Alison pattern.

If the Alison circle is monument-enriched / settlement-depleted,
the anti-divergence circle should be settlement-enriched / monument-depleted.

Steps:
  1. Coarse scan at 5° resolution → find anti-divergence pole
  2. Refine at 1° resolution
  3. Trace the circle and characterize terrain/geography
  4. Map comparison with Alison circle
  5. Angular relationship between poles

Outputs:
  outputs/extended_controls/anti_divergence.json
  outputs/extended_controls/anti_divergence_map.png
  outputs/extended_controls/terrain_comparison.json
"""

import csv, json, math, os, sys, time
import numpy as np

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
np.random.seed(42)

# ============================================================
# CONFIGURATION
# ============================================================
ALISON_POLE_LAT = 59.682122
ALISON_POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * np.pi / 2
THRESHOLD_KM = 50
STRIP_FRACTION = np.sin(THRESHOLD_KM / EARTH_R_KM)

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "extended_controls")
os.makedirs(OUT_DIR, exist_ok=True)

# Pleiades type classifications (matching existing codebase)
MONUMENTAL_TYPES = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus", "shrine"
}
SETTLEMENT_TYPES = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement", "townhouse", "production"
}

# ============================================================
# VECTORIZED HELPER FUNCTIONS
# ============================================================
def dist_from_gc_vec(site_lats, site_lons, pole_lat, pole_lon):
    """Vectorized distance from each site to the great circle defined by pole."""
    lat1r = np.radians(pole_lat)
    lon1r = np.radians(pole_lon)
    lat2r = np.radians(site_lats)
    lon2r = np.radians(site_lons)
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat/2)**2 + np.cos(lat1r)*np.cos(lat2r)*np.sin(dlon/2)**2
    d = EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QUARTER_CIRC)

def count_within_vec(site_lats, site_lons, pole_lat, pole_lon):
    dists = dist_from_gc_vec(site_lats, site_lons, pole_lat, pole_lon)
    return int(np.sum(dists <= THRESHOLD_KM))

def poisson_z(observed, expected):
    """Poisson-approximation Z-score."""
    if expected <= 0:
        return 0.0
    return (observed - expected) / np.sqrt(expected)

def haversine_deg(lat1, lon1, lat2, lon2):
    """Angular distance in degrees between two points."""
    lat1r, lon1r = np.radians(lat1), np.radians(lon1)
    lat2r, lon2r = np.radians(lat2), np.radians(lon2)
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat/2)**2 + np.cos(lat1r)*np.cos(lat2r)*np.sin(dlon/2)**2
    return np.degrees(2 * np.arcsin(np.sqrt(min(1.0, a))))

def to_xyz(lat, lon):
    """Convert lat/lon to unit 3D Cartesian."""
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    return np.array([
        np.cos(lat_r) * np.cos(lon_r),
        np.cos(lat_r) * np.sin(lon_r),
        np.sin(lat_r)
    ])

def trace_great_circle(pole_lat, pole_lon, n_points=360):
    """Trace points along the great circle defined by its pole.
    Returns list of (lat, lon) tuples."""
    pole = to_xyz(pole_lat, pole_lon)
    # Find two orthogonal vectors in the plane perpendicular to the pole
    # Pick an arbitrary vector not parallel to pole
    if abs(pole[2]) < 0.9:
        arbitrary = np.array([0, 0, 1])
    else:
        arbitrary = np.array([1, 0, 0])
    u = np.cross(pole, arbitrary)
    u = u / np.linalg.norm(u)
    v = np.cross(pole, u)
    v = v / np.linalg.norm(v)

    points = []
    for i in range(n_points):
        theta = 2 * np.pi * i / n_points
        p = np.cos(theta) * u + np.sin(theta) * v
        lat = np.degrees(np.arcsin(np.clip(p[2], -1, 1)))
        lon = np.degrees(np.arctan2(p[1], p[0]))
        points.append((lat, lon))
    return points

# ============================================================
# LOAD PLEIADES DATA
# ============================================================
print("=" * 70)
print("DIRECTIVE 2: ANTI-DIVERGENCE CIRCLES")
print("=" * 70)
print("\nLoading Pleiades data...")

pleiades_path = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")
mon_list = []
set_list = []

with open(pleiades_path, encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            lat = float(row["reprLat"])
            lon = float(row["reprLong"])
            if lat == 0 and lon == 0:
                continue
            feature_types = {t.strip() for t in row.get("featureTypes", "").split(",")}
            is_mon = bool(feature_types & MONUMENTAL_TYPES)
            is_set = bool(feature_types & SETTLEMENT_TYPES)
            if is_mon:
                mon_list.append((lat, lon))
            elif is_set:
                set_list.append((lat, lon))
        except (ValueError, KeyError, TypeError):
            pass

mon_arr = np.array(mon_list)
set_arr = np.array(set_list)
mon_lats, mon_lons = mon_arr[:, 0], mon_arr[:, 1]
set_lats, set_lons = set_arr[:, 0], set_arr[:, 1]

N_MON = len(mon_arr)
N_SET = len(set_arr)
EXPECTED_MON = N_MON * STRIP_FRACTION
EXPECTED_SET = N_SET * STRIP_FRACTION

print(f"Monumental sites: {N_MON}")
print(f"Settlement sites: {N_SET}")
print(f"Strip fraction (50km): {STRIP_FRACTION:.6f}")
print(f"Expected monuments per circle: {EXPECTED_MON:.1f}")
print(f"Expected settlements per circle: {EXPECTED_SET:.1f}")

# ============================================================
# STEP 1: COARSE SCAN — 5° resolution
# ============================================================
print(f"\n{'='*70}")
print("STEP 1: COARSE SCAN (5° resolution)")
print(f"{'='*70}")

COARSE_STEP = 5
coarse_results = []
best_anti_D = -999
best_anti_pole = None

# Use antipodal symmetry: only scan one hemisphere
poles_coarse = []
for lat in np.arange(-90, 90 + COARSE_STEP, COARSE_STEP):
    for lon in np.arange(-180, 180, COARSE_STEP):
        if lat > 0:
            poles_coarse.append((lat, lon))
        elif lat == 0 and lon >= 0:
            poles_coarse.append((lat, lon))
# Deduplicate polar
poles_dedup = []
seen_polar = False
for lat, lon in poles_coarse:
    if abs(lat) == 90:
        if not seen_polar:
            poles_dedup.append((lat, lon))
            seen_polar = True
    else:
        poles_dedup.append((lat, lon))
poles_coarse = poles_dedup

print(f"Total coarse poles: {len(poles_coarse)}")
t0 = time.time()

for i, (plat, plon) in enumerate(poles_coarse):
    mon_near = count_within_vec(mon_lats, mon_lons, plat, plon)
    set_near = count_within_vec(set_lats, set_lons, plat, plon)

    z_mon = poisson_z(mon_near, EXPECTED_MON)
    z_set = poisson_z(set_near, EXPECTED_SET)

    # Anti-divergence: settlement enrichment minus monument enrichment
    # High anti_D = lots of settlements near, few monuments near
    anti_D = z_set - z_mon

    # Also compute the original Alison-style D for reference
    alison_D = z_mon - z_set

    entry = {
        'lat': float(plat), 'lon': float(plon),
        'mon_near': mon_near, 'set_near': set_near,
        'z_mon': round(z_mon, 3), 'z_set': round(z_set, 3),
        'anti_D': round(anti_D, 3), 'alison_D': round(alison_D, 3)
    }
    coarse_results.append(entry)

    if anti_D > best_anti_D:
        best_anti_D = anti_D
        best_anti_pole = (plat, plon)

    if (i + 1) % 500 == 0:
        print(f"  Scanned {i+1}/{len(poles_coarse)} poles...")

elapsed = time.time() - t0
print(f"\nCoarse scan complete in {elapsed:.1f}s")
print(f"Best anti-divergence pole (coarse): ({best_anti_pole[0]:.0f}°, {best_anti_pole[1]:.0f}°)")
print(f"  anti_D = {best_anti_D:.3f}")

# Also compute Alison pole anti_D for comparison
alison_mon_near = count_within_vec(mon_lats, mon_lons, ALISON_POLE_LAT, ALISON_POLE_LON)
alison_set_near = count_within_vec(set_lats, set_lons, ALISON_POLE_LAT, ALISON_POLE_LON)
alison_z_mon = poisson_z(alison_mon_near, EXPECTED_MON)
alison_z_set = poisson_z(alison_set_near, EXPECTED_SET)
alison_anti_D = alison_z_set - alison_z_mon
alison_D = alison_z_mon - alison_z_set
print(f"\nAlison pole anti_D = {alison_anti_D:.3f}  (alison_D = {alison_D:.3f})")
print(f"  Monuments near: {alison_mon_near}, Settlements near: {alison_set_near}")

# Top 10 anti-divergence poles from coarse scan
coarse_sorted = sorted(coarse_results, key=lambda x: x['anti_D'], reverse=True)
print("\nTop 10 anti-divergence poles (coarse):")
for rank, entry in enumerate(coarse_sorted[:10], 1):
    print(f"  #{rank}: ({entry['lat']:6.1f}°, {entry['lon']:7.1f}°)  "
          f"anti_D={entry['anti_D']:+.3f}  mon={entry['mon_near']}  set={entry['set_near']}  "
          f"z_mon={entry['z_mon']:+.3f}  z_set={entry['z_set']:+.3f}")

# ============================================================
# STEP 2: REFINE at 1° resolution
# ============================================================
print(f"\n{'='*70}")
print("STEP 2: REFINE BEST POLE (1° resolution)")
print(f"{'='*70}")

# Refine around top 5 coarse peaks (they may be in different regions)
# First, identify distinct peaks (>10° apart)
distinct_peaks = []
for entry in coarse_sorted:
    is_distinct = True
    for pk in distinct_peaks:
        if haversine_deg(entry['lat'], entry['lon'], pk['lat'], pk['lon']) < 15:
            is_distinct = False
            break
    if is_distinct:
        distinct_peaks.append(entry)
    if len(distinct_peaks) >= 5:
        break

print(f"Refining around {len(distinct_peaks)} distinct peaks...")

FINE_STEP = 1
fine_results = []
best_fine_anti_D = -999
best_fine_pole = None

for pk_idx, peak in enumerate(distinct_peaks):
    lat0, lon0 = peak['lat'], peak['lon']
    print(f"\n  Peak {pk_idx+1}: ({lat0:.0f}°, {lon0:.0f}°) coarse anti_D={peak['anti_D']:.3f}")

    for lat in np.arange(lat0 - 5, lat0 + 5.5, FINE_STEP):
        for lon in np.arange(lon0 - 5, lon0 + 5.5, FINE_STEP):
            # Wrap coordinates
            if lat > 90:
                lat = 90
            if lat < -90:
                lat = -90
            lon_w = ((lon + 180) % 360) - 180

            mon_near = count_within_vec(mon_lats, mon_lons, lat, lon_w)
            set_near = count_within_vec(set_lats, set_lons, lat, lon_w)
            z_mon = poisson_z(mon_near, EXPECTED_MON)
            z_set = poisson_z(set_near, EXPECTED_SET)
            anti_D = z_set - z_mon

            entry = {
                'lat': round(float(lat), 1), 'lon': round(float(lon_w), 1),
                'mon_near': mon_near, 'set_near': set_near,
                'z_mon': round(z_mon, 3), 'z_set': round(z_set, 3),
                'anti_D': round(anti_D, 3),
                'peak_parent': pk_idx
            }
            fine_results.append(entry)

            if anti_D > best_fine_anti_D:
                best_fine_anti_D = anti_D
                best_fine_pole = (float(lat), float(lon_w))

fine_sorted = sorted(fine_results, key=lambda x: x['anti_D'], reverse=True)
print(f"\nBest anti-divergence pole (refined): ({best_fine_pole[0]:.1f}°, {best_fine_pole[1]:.1f}°)")
print(f"  anti_D = {best_fine_anti_D:.3f}")

# Top 10 refined results
print("\nTop 10 refined anti-divergence poles:")
fine_distinct = []
for entry in fine_sorted:
    is_distinct = True
    for pk in fine_distinct:
        if haversine_deg(entry['lat'], entry['lon'], pk['lat'], pk['lon']) < 5:
            is_distinct = False
            break
    if is_distinct:
        fine_distinct.append(entry)
    if len(fine_distinct) >= 10:
        break

for rank, entry in enumerate(fine_distinct, 1):
    print(f"  #{rank}: ({entry['lat']:6.1f}°, {entry['lon']:7.1f}°)  "
          f"anti_D={entry['anti_D']:+.3f}  mon={entry['mon_near']}  set={entry['set_near']}  "
          f"z_mon={entry['z_mon']:+.3f}  z_set={entry['z_set']:+.3f}")

# ============================================================
# STEP 3: CHARACTERIZE THE ANTI-DIVERGENCE CIRCLE
# ============================================================
print(f"\n{'='*70}")
print("STEP 3: TRACE AND CHARACTERIZE CIRCLES")
print(f"{'='*70}")

# Trace both circles
anti_pole = best_fine_pole
anti_circle = trace_great_circle(anti_pole[0], anti_pole[1], n_points=720)
alison_circle = trace_great_circle(ALISON_POLE_LAT, ALISON_POLE_LON, n_points=720)

# Geographic region classification
def classify_region(lat, lon):
    """Simple geographic region classifier."""
    if lat < -60:
        return "Antarctica"
    elif lat > 66:
        return "Arctic"
    elif -35 <= lat <= 37 and -20 <= lon <= 55:
        return "Africa"
    elif 25 <= lat <= 72 and -12 <= lon <= 45:
        return "Europe"
    elif 10 <= lat <= 55 and 55 <= lon <= 145:
        return "Asia (mainland)"
    elif -10 <= lat <= 10 and 95 <= lon <= 150:
        return "SE Asia / Indonesia"
    elif -50 <= lat <= -10 and 110 <= lon <= 180:
        return "Australia/Oceania"
    elif 7 <= lat <= 72 and -170 <= lon <= -50:
        return "North America"
    elif -57 <= lat <= 15 and -85 <= lon <= -30:
        return "South America"
    elif 25 <= lat <= 45 and 25 <= lon <= 65:
        return "Middle East"
    else:
        return "Ocean/Other"

# Classify regions for both circles
def region_breakdown(circle_pts):
    """Count points per region."""
    from collections import Counter
    regions = Counter()
    for lat, lon in circle_pts:
        regions[classify_region(lat, lon)] += 1
    total = len(circle_pts)
    return {k: {'count': v, 'fraction': round(v/total, 3)} for k, v in regions.most_common()}

anti_regions = region_breakdown(anti_circle)
alison_regions = region_breakdown(alison_circle)

print("\nAnti-divergence circle regions:")
for region, info in anti_regions.items():
    print(f"  {region}: {info['count']} pts ({info['fraction']*100:.1f}%)")

print("\nAlison circle regions:")
for region, info in alison_regions.items():
    print(f"  {region}: {info['count']} pts ({info['fraction']*100:.1f}%)")

# Latitude distribution
anti_lats = [p[0] for p in anti_circle]
alison_lats = [p[0] for p in alison_circle]

print(f"\nLatitude range — Anti-divergence: [{min(anti_lats):.1f}°, {max(anti_lats):.1f}°]")
print(f"Latitude range — Alison:          [{min(alison_lats):.1f}°, {max(alison_lats):.1f}°]")

# ============================================================
# STEP 4: LIST SITES ON / NEAR EACH CIRCLE
# ============================================================
print(f"\n{'='*70}")
print("STEP 4: SITES NEAR THE ANTI-DIVERGENCE CIRCLE")
print(f"{'='*70}")

# Get all monuments and settlements near anti-divergence circle
anti_mon_dists = dist_from_gc_vec(mon_lats, mon_lons, anti_pole[0], anti_pole[1])
anti_set_dists = dist_from_gc_vec(set_lats, set_lons, anti_pole[0], anti_pole[1])

anti_mon_count = int(np.sum(anti_mon_dists <= THRESHOLD_KM))
anti_set_count = int(np.sum(anti_set_dists <= THRESHOLD_KM))
print(f"Monuments within 50km: {anti_mon_count}  (expected: {EXPECTED_MON:.1f})")
print(f"Settlements within 50km: {anti_set_count}  (expected: {EXPECTED_SET:.1f})")
print(f"Settlement/Monument ratio: {anti_set_count/max(anti_mon_count,1):.2f}")

# Compare to Alison
print(f"\nAlison circle comparison:")
print(f"Monuments within 50km: {alison_mon_near}  (expected: {EXPECTED_MON:.1f})")
print(f"Settlements within 50km: {alison_set_near}  (expected: {EXPECTED_SET:.1f})")
print(f"Settlement/Monument ratio: {alison_set_near/max(alison_mon_near,1):.2f}")

# ============================================================
# STEP 5: ANGULAR RELATIONSHIP
# ============================================================
print(f"\n{'='*70}")
print("STEP 5: ANGULAR RELATIONSHIP BETWEEN POLES")
print(f"{'='*70}")

angular_sep = haversine_deg(anti_pole[0], anti_pole[1], ALISON_POLE_LAT, ALISON_POLE_LON)
print(f"Anti-divergence pole: ({anti_pole[0]:.1f}°, {anti_pole[1]:.1f}°)")
print(f"Alison pole:          ({ALISON_POLE_LAT:.1f}°, {ALISON_POLE_LON:.1f}°)")
print(f"Angular separation:   {angular_sep:.1f}°")

if angular_sep < 20:
    interp = "NEAR-IDENTICAL — both patterns occupy similar geography (unexpected)"
elif angular_sep < 70:
    interp = "MODERATE — circles intersect but have distinct paths"
elif 70 <= angular_sep <= 110:
    interp = "NEAR-PERPENDICULAR — monuments and settlements occupy orthogonal terrain types"
elif angular_sep > 160:
    interp = "NEAR-ANTIPODAL — same circle, opposite enrichment (D sign flip)"
else:
    interp = "OBLIQUE — intermediate relationship"

print(f"Interpretation:       {interp}")

# Also check antipodal separation (since a pole and its antipode define the same GC)
anti_antipodal_sep = haversine_deg(-anti_pole[0], anti_pole[1] + 180 if anti_pole[1] < 0 else anti_pole[1] - 180,
                                    ALISON_POLE_LAT, ALISON_POLE_LON)
print(f"Antipodal separation: {anti_antipodal_sep:.1f}° (checking if they define the same GC)")
effective_sep = min(angular_sep, anti_antipodal_sep)
print(f"Effective separation: {effective_sep:.1f}° (min of direct and antipodal)")

# ============================================================
# STEP 6: STATISTICAL CONTEXT
# ============================================================
print(f"\n{'='*70}")
print("STEP 6: STATISTICAL CONTEXT")
print(f"{'='*70}")

# Where does the anti-divergence pole rank among all poles for anti_D?
anti_d_values = [r['anti_D'] for r in coarse_results]
anti_d_arr = np.array(anti_d_values)
anti_percentile = 100 * np.sum(anti_d_arr <= best_anti_D) / len(anti_d_arr)

# Also compute the D distribution stats
print(f"Anti-D distribution across {len(coarse_results)} coarse poles:")
print(f"  Mean:   {np.mean(anti_d_arr):.3f}")
print(f"  Std:    {np.std(anti_d_arr):.3f}")
print(f"  Min:    {np.min(anti_d_arr):.3f}")
print(f"  Max:    {np.max(anti_d_arr):.3f}")
print(f"  Best anti-D percentile: {anti_percentile:.1f}%")

# Z-score of the best anti_D relative to the distribution
anti_D_zscore = (best_anti_D - np.mean(anti_d_arr)) / np.std(anti_d_arr)
print(f"  Z-score of best anti_D: {anti_D_zscore:.2f}")

# Compare: what's the best ALISON-style D?
alison_d_values = [r['alison_D'] for r in coarse_results]
alison_d_arr = np.array(alison_d_values)
best_alison_D = np.max(alison_d_arr)
alison_D_zscore = (best_alison_D - np.mean(alison_d_arr)) / np.std(alison_d_arr)
print(f"\nAlison-D (monument enrichment) best: {best_alison_D:.3f}  Z={alison_D_zscore:.2f}")
print(f"Anti-D (settlement enrichment) best:  {best_anti_D:.3f}  Z={anti_D_zscore:.2f}")

symmetry = best_anti_D / best_alison_D if best_alison_D != 0 else 0
print(f"Symmetry ratio (anti_D / alison_D):   {symmetry:.2f}")
print(f"  (1.0 = perfectly symmetric; >1 = settlement pattern stronger; <1 = monument pattern stronger)")

# ============================================================
# SAVE RESULTS
# ============================================================
print(f"\n{'='*70}")
print("SAVING RESULTS")
print(f"{'='*70}")

# Main results JSON
results = {
    'anti_divergence_pole': {
        'lat': round(anti_pole[0], 1),
        'lon': round(anti_pole[1], 1),
        'anti_D_coarse': round(best_anti_D, 3),
        'anti_D_refined': round(best_fine_anti_D, 3),
        'mon_near': anti_mon_count,
        'set_near': anti_set_count,
        'z_mon': round(poisson_z(anti_mon_count, EXPECTED_MON), 3),
        'z_set': round(poisson_z(anti_set_count, EXPECTED_SET), 3),
    },
    'alison_pole': {
        'lat': ALISON_POLE_LAT,
        'lon': ALISON_POLE_LON,
        'alison_D': round(alison_D, 3),
        'anti_D': round(alison_anti_D, 3),
        'mon_near': alison_mon_near,
        'set_near': alison_set_near,
        'z_mon': round(alison_z_mon, 3),
        'z_set': round(alison_z_set, 3),
    },
    'angular_relationship': {
        'direct_separation_deg': round(angular_sep, 1),
        'antipodal_separation_deg': round(anti_antipodal_sep, 1),
        'effective_separation_deg': round(effective_sep, 1),
        'interpretation': interp,
    },
    'statistics': {
        'n_coarse_poles': len(coarse_results),
        'anti_D_mean': round(float(np.mean(anti_d_arr)), 3),
        'anti_D_std': round(float(np.std(anti_d_arr)), 3),
        'anti_D_percentile': round(anti_percentile, 1),
        'anti_D_zscore': round(anti_D_zscore, 2),
        'best_alison_D': round(float(best_alison_D), 3),
        'symmetry_ratio': round(symmetry, 2),
    },
    'top10_anti_divergence_poles': [
        {
            'rank': i+1,
            'lat': e['lat'], 'lon': e['lon'],
            'anti_D': e['anti_D'],
            'mon_near': e['mon_near'], 'set_near': e['set_near'],
            'z_mon': e['z_mon'], 'z_set': e['z_set']
        }
        for i, e in enumerate(fine_distinct[:10])
    ],
    'coarse_scan': coarse_results,
    'data_summary': {
        'n_monuments': N_MON,
        'n_settlements': N_SET,
        'threshold_km': THRESHOLD_KM,
        'strip_fraction': round(STRIP_FRACTION, 6),
        'expected_mon': round(EXPECTED_MON, 1),
        'expected_set': round(EXPECTED_SET, 1),
    }
}

with open(os.path.join(OUT_DIR, "anti_divergence.json"), 'w') as f:
    json.dump(results, f, indent=2)
print("Saved: anti_divergence.json")

# Terrain/region comparison JSON
terrain_comparison = {
    'anti_divergence_circle': {
        'pole': {'lat': round(anti_pole[0], 1), 'lon': round(anti_pole[1], 1)},
        'regions': anti_regions,
        'latitude_range': [round(min(anti_lats), 1), round(max(anti_lats), 1)],
        'settlement_monument_ratio': round(anti_set_count / max(anti_mon_count, 1), 2),
    },
    'alison_circle': {
        'pole': {'lat': ALISON_POLE_LAT, 'lon': ALISON_POLE_LON},
        'regions': alison_regions,
        'latitude_range': [round(min(alison_lats), 1), round(max(alison_lats), 1)],
        'settlement_monument_ratio': round(alison_set_near / max(alison_mon_near, 1), 2),
    },
    'angular_separation': round(effective_sep, 1),
}

with open(os.path.join(OUT_DIR, "terrain_comparison.json"), 'w') as f:
    json.dump(terrain_comparison, f, indent=2)
print("Saved: terrain_comparison.json")

# ============================================================
# STEP 7: MAP VISUALIZATION
# ============================================================
print(f"\n{'='*70}")
print("STEP 7: GENERATING MAP")
print(f"{'='*70}")

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(16, 9))

    # World map background (simple coastline approximation)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_facecolor('#1a1a2e')
    ax.set_aspect('equal')

    # Plot site distributions
    ax.scatter(mon_lons, mon_lats, s=2, c='#FFD700', alpha=0.15, label=f'Monuments ({N_MON})', zorder=2)
    ax.scatter(set_lons, set_lats, s=2, c='#4169E1', alpha=0.15, label=f'Settlements ({N_SET})', zorder=2)

    # Plot Alison circle
    al_lats = [p[0] for p in alison_circle]
    al_lons = [p[1] for p in alison_circle]
    # Sort by longitude to avoid wrap-around lines
    for i in range(len(al_lons)):
        ax.plot(al_lons[i], al_lats[i], '.', color='#FFD700', markersize=1.5, alpha=0.8, zorder=3)

    # Plot anti-divergence circle
    ad_lats = [p[0] for p in anti_circle]
    ad_lons = [p[1] for p in anti_circle]
    for i in range(len(ad_lons)):
        ax.plot(ad_lons[i], ad_lats[i], '.', color='#4169E1', markersize=1.5, alpha=0.8, zorder=3)

    # Plot poles
    ax.plot(ALISON_POLE_LON, ALISON_POLE_LAT, '*', color='#FFD700', markersize=15,
            markeredgecolor='white', markeredgewidth=1, zorder=5, label='Alison pole')
    ax.plot(anti_pole[1], anti_pole[0], '*', color='#4169E1', markersize=15,
            markeredgecolor='white', markeredgewidth=1, zorder=5, label='Anti-divergence pole')

    # Highlight sites on anti-divergence circle
    anti_mon_mask = anti_mon_dists <= THRESHOLD_KM
    anti_set_mask = anti_set_dists <= THRESHOLD_KM
    if np.any(anti_mon_mask):
        ax.scatter(mon_lons[anti_mon_mask], mon_lats[anti_mon_mask], s=20, c='#FFD700',
                   edgecolors='white', linewidth=0.5, alpha=0.9, zorder=4)
    if np.any(anti_set_mask):
        ax.scatter(set_lons[anti_set_mask], set_lats[anti_set_mask], s=20, c='#4169E1',
                   edgecolors='white', linewidth=0.5, alpha=0.9, zorder=4)

    # Grid
    ax.grid(True, alpha=0.2, color='gray')
    ax.set_xlabel('Longitude', fontsize=12)
    ax.set_ylabel('Latitude', fontsize=12)
    ax.set_title(
        f'Anti-Divergence Circle vs Alison Circle\n'
        f'Gold = Alison (monument-enriched, D={alison_D:+.2f})  |  '
        f'Blue = Anti-divergence (settlement-enriched, anti_D={best_fine_anti_D:+.2f})\n'
        f'Angular separation: {effective_sep:.1f}°',
        fontsize=13, color='white', pad=15
    )
    ax.tick_params(colors='gray')

    # Legend
    leg = ax.legend(loc='lower left', fontsize=9, facecolor='#1a1a2e', edgecolor='gray',
                    labelcolor='white')

    # Annotation box
    stats_text = (
        f"Anti-div circle: {anti_set_count} settlements, {anti_mon_count} monuments\n"
        f"Alison circle: {alison_set_near} settlements, {alison_mon_near} monuments\n"
        f"Expected: {EXPECTED_SET:.0f} set, {EXPECTED_MON:.0f} mon per circle"
    )
    ax.text(175, -80, stats_text, fontsize=8, color='white', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#1a1a2e', edgecolor='gray', alpha=0.8))

    fig.patch.set_facecolor('#0f0f1a')
    plt.tight_layout()
    map_path = os.path.join(OUT_DIR, "anti_divergence_map.png")
    plt.savefig(map_path, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close()
    print(f"Saved: anti_divergence_map.png")

    # ============================================================
    # BONUS: Anti-D heatmap
    # ============================================================
    print("Generating anti_D heatmap...")

    # Build 2D grid from coarse results
    lats_unique = sorted(set(r['lat'] for r in coarse_results))
    lons_unique = sorted(set(r['lon'] for r in coarse_results))

    lat_to_idx = {lat: i for i, lat in enumerate(lats_unique)}
    lon_to_idx = {lon: i for i, lon in enumerate(lons_unique)}

    grid = np.full((len(lats_unique), len(lons_unique)), np.nan)
    for r in coarse_results:
        if r['lat'] in lat_to_idx and r['lon'] in lon_to_idx:
            grid[lat_to_idx[r['lat']], lon_to_idx[r['lon']]] = r['anti_D']

    fig2, ax2 = plt.subplots(1, 1, figsize=(16, 8))
    lon_arr = np.array(lons_unique)
    lat_arr = np.array(lats_unique)
    im = ax2.pcolormesh(lon_arr, lat_arr, grid, cmap='RdYlBu', shading='auto')
    plt.colorbar(im, ax=ax2, label='Anti-D (settlement enrichment - monument enrichment)')

    # Mark the best anti-divergence pole
    ax2.plot(anti_pole[1], anti_pole[0], '*', color='blue', markersize=20,
             markeredgecolor='white', markeredgewidth=2, zorder=5, label='Anti-div pole')
    ax2.plot(ALISON_POLE_LON, ALISON_POLE_LAT, '*', color='gold', markersize=20,
             markeredgecolor='white', markeredgewidth=2, zorder=5, label='Alison pole')

    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')
    ax2.set_title('Anti-Divergence Score Heatmap (settlement enrichment − monument enrichment)')
    ax2.legend(loc='lower left')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    heatmap_path = os.path.join(OUT_DIR, "anti_divergence_heatmap.png")
    plt.savefig(heatmap_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: anti_divergence_heatmap.png")

except ImportError as e:
    print(f"WARNING: Could not generate maps (matplotlib not available): {e}")

# ============================================================
# FINAL SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")
print(f"""
Anti-Divergence Pole: ({anti_pole[0]:.1f}°, {anti_pole[1]:.1f}°)
  anti_D = {best_fine_anti_D:.3f}  (z_set - z_mon)
  Settlements within 50km: {anti_set_count}  (expected {EXPECTED_SET:.0f})
  Monuments within 50km:   {anti_mon_count}  (expected {EXPECTED_MON:.0f})

Alison Pole: ({ALISON_POLE_LAT:.1f}°, {ALISON_POLE_LON:.1f}°)
  alison_D = {alison_D:.3f}  (z_mon - z_set)
  Monuments within 50km:   {alison_mon_near}  (expected {EXPECTED_MON:.0f})
  Settlements within 50km: {alison_set_near}  (expected {EXPECTED_SET:.0f})

Angular Separation: {effective_sep:.1f}° ({interp})

Statistical Context:
  Anti-D percentile: {anti_percentile:.1f}%
  Anti-D Z-score: {anti_D_zscore:.2f}
  Symmetry ratio: {symmetry:.2f}
""")

print("Outputs saved to: outputs/extended_controls/")
print("  - anti_divergence.json")
print("  - terrain_comparison.json")
print("  - anti_divergence_map.png")
print("  - anti_divergence_heatmap.png")
print("\nDone.")
