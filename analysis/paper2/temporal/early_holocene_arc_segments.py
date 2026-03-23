#!/usr/bin/env python3
"""
Early Holocene Arc Segment Decomposition
==========================================
For p3k14c sites in the 10,500–8,500 BP range:
  1. Identify sites within 50km of Alison's Great Circle
  2. Project each near-circle site onto the circle to get its azimuthal position
  3. Bin into 12 × 30-degree arc segments
  4. Compare observed vs expected counts per segment (MC baseline)

The azimuthal position is computed as the bearing from the great circle pole
to the site — this gives the "longitude" in the pole-centered frame, which
maps directly to position along the great circle.
"""

import csv, math, json, os, sys
import numpy as np
from collections import defaultdict, Counter

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
THRESHOLD_KM = 50
N_TRIALS = 1000
N_SEGMENTS = 12
SEGMENT_WIDTH = 360.0 / N_SEGMENTS  # 30 degrees

# Early Holocene in 14C years BP
AGE_MIN_BP = 8500
AGE_MAX_BP = 10500

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "early_holocene_arcs")
os.makedirs(OUT_DIR, exist_ok=True)


# ============================================================
# HELPERS
# ============================================================
def gc_dist(lat1, lon1, lat2, lon2):
    """Haversine distance in km."""
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))


def dist_from_circle(lat, lon):
    """Distance from site to the great circle (equator of the pole)."""
    d = gc_dist(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QUARTER_CIRC)


def bearing_from_pole(lat, lon):
    """
    Compute initial bearing from the pole to the site.
    This gives the azimuthal position along the great circle.
    Returns degrees [0, 360).
    """
    lat1 = math.radians(POLE_LAT)
    lon1 = math.radians(POLE_LON)
    lat2 = math.radians(lat)
    lon2 = math.radians(lon)
    dlon = lon2 - lon1
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(dlon)
    bearing = math.degrees(math.atan2(x, y))
    return bearing % 360


def segment_index(bearing):
    """Map bearing [0,360) to segment index [0, N_SEGMENTS)."""
    return int(bearing / SEGMENT_WIDTH) % N_SEGMENTS


def segment_label(idx):
    """Human-readable label for a segment."""
    start = idx * SEGMENT_WIDTH
    end = start + SEGMENT_WIDTH
    return f"{start:.0f}°–{end:.0f}°"


# ============================================================
# LOAD p3k14c DATA
# ============================================================
print("=" * 70)
print("EARLY HOLOCENE ARC SEGMENT DECOMPOSITION")
print(f"  Age range: {AGE_MIN_BP}–{AGE_MAX_BP} BP")
print(f"  Great circle threshold: {THRESHOLD_KM} km")
print(f"  Segments: {N_SEGMENTS} × {SEGMENT_WIDTH:.0f}°")
print("=" * 70)

p3k_file = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")

# Deduplicate by SiteID, keeping oldest age
site_groups = defaultdict(list)
with open(p3k_file, encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            lat = float(row['Lat'])
            lon = float(row['Long'])
            age_bp = float(row['Age']) if row['Age'] else None
        except (ValueError, TypeError):
            continue
        if lat == 0 and lon == 0:
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue
        if age_bp is None:
            continue
        site_id = row.get('SiteID', '') or f"anon_{lat:.4f}_{lon:.4f}"
        site_groups[site_id].append({
            'lat': lat, 'lon': lon, 'age_bp': age_bp,
            'site_name': row.get('SiteName', ''),
            'country': row.get('Country', ''),
            'continent': row.get('Continent', ''),
        })

# Aggregate: mean location, oldest age per site
all_sites = []
for sid, rows in site_groups.items():
    ages = [r['age_bp'] for r in rows]
    all_sites.append({
        'lat': np.mean([r['lat'] for r in rows]),
        'lon': np.mean([r['lon'] for r in rows]),
        'age_bp': max(ages),
        'site_name': rows[0]['site_name'],
        'country': rows[0]['country'],
        'continent': rows[0]['continent'],
        'site_id': sid,
    })

print(f"\nTotal unique p3k14c sites: {len(all_sites)}")

# Filter to early Holocene
early_holo = [s for s in all_sites if AGE_MIN_BP <= s['age_bp'] <= AGE_MAX_BP]
print(f"Sites in {AGE_MIN_BP}–{AGE_MAX_BP} BP: {len(early_holo)}")

if len(early_holo) < 10:
    print("ERROR: Too few sites in this time range. Aborting.")
    sys.exit(1)

# ============================================================
# IDENTIFY NEAR-CIRCLE SITES & COMPUTE ARC POSITIONS
# ============================================================
near_circle = []
for s in early_holo:
    d = dist_from_circle(s['lat'], s['lon'])
    if d <= THRESHOLD_KM:
        b = bearing_from_pole(s['lat'], s['lon'])
        seg = segment_index(b)
        near_circle.append({
            **s,
            'dist_km': round(d, 2),
            'bearing': round(b, 2),
            'segment': seg,
            'segment_label': segment_label(seg),
        })

print(f"Near-circle sites (≤{THRESHOLD_KM} km): {len(near_circle)}")
print(f"Fraction: {len(near_circle)/len(early_holo):.4f}")

# Count per segment
observed_counts = Counter(s['segment'] for s in near_circle)
observed_per_seg = [observed_counts.get(i, 0) for i in range(N_SEGMENTS)]

print(f"\nObserved per 30° segment:")
for i in range(N_SEGMENTS):
    bar = "█" * observed_per_seg[i]
    print(f"  {segment_label(i):>10s}  {observed_per_seg[i]:>3d}  {bar}")

# ============================================================
# ALSO DECOMPOSE ALL EARLY HOLOCENE SITES (not just near-circle)
# to show the spatial coverage
# ============================================================
all_bearings = [bearing_from_pole(s['lat'], s['lon']) for s in early_holo]
all_seg_counts = Counter(segment_index(b) for b in all_bearings)
all_per_seg = [all_seg_counts.get(i, 0) for i in range(N_SEGMENTS)]

print(f"\nAll early Holocene sites per 30° segment (spatial coverage):")
for i in range(N_SEGMENTS):
    print(f"  {segment_label(i):>10s}  {all_per_seg[i]:>4d}")

# ============================================================
# MONTE CARLO: EXPECTED SEGMENT DISTRIBUTION
# ============================================================
print(f"\nRunning {N_TRIALS}-trial Monte Carlo...")
np.random.seed(42)

lats_all_eh = np.array([s['lat'] for s in early_holo])
lons_all_eh = np.array([s['lon'] for s in early_holo])
n_eh = len(early_holo)

# For each trial: shuffle sites, count how many fall near-circle per segment
# Use distribution-matched random: pick from the same pool with replacement + jitter
mc_segment_counts = np.zeros((N_TRIALS, N_SEGMENTS))
mc_total_near = []

for t in range(N_TRIALS):
    # Draw n_eh random points matched to the spatial distribution
    idx = np.random.randint(0, n_eh, size=n_eh)
    rlat = lats_all_eh[idx] + np.random.normal(0, 2, n_eh)
    rlon = lons_all_eh[idx] + np.random.normal(0, 2, n_eh)
    rlat = np.clip(rlat, -90, 90)
    rlon = np.clip(rlon, -180, 180)

    # Compute distances from circle (vectorized)
    pole_lat_r = np.radians(POLE_LAT)
    pole_lon_r = np.radians(POLE_LON)
    lat2_r = np.radians(rlat)
    lon2_r = np.radians(rlon)
    dlat = lat2_r - pole_lat_r
    dlon = lon2_r - pole_lon_r
    a = np.sin(dlat/2)**2 + np.cos(pole_lat_r)*np.cos(lat2_r)*np.sin(dlon/2)**2
    d_pole = EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    d_circle = np.abs(d_pole - QUARTER_CIRC)

    near_mask = d_circle <= THRESHOLD_KM
    mc_total_near.append(int(np.sum(near_mask)))

    # For near-circle random points, compute bearings and segments
    if np.sum(near_mask) > 0:
        near_lats = rlat[near_mask]
        near_lons = rlon[near_mask]
        for nl, nln in zip(near_lats, near_lons):
            b = bearing_from_pole(nl, nln)
            seg = segment_index(b)
            mc_segment_counts[t, seg] += 1

# Expected counts per segment
expected_per_seg = np.mean(mc_segment_counts, axis=0)
std_per_seg = np.std(mc_segment_counts, axis=0)
expected_total = np.mean(mc_total_near)
std_total = np.std(mc_total_near)

total_observed = len(near_circle)
z_total = (total_observed - expected_total) / std_total if std_total > 0 else 0

print(f"\nOverall: observed={total_observed}, expected={expected_total:.1f} ± {std_total:.1f}, Z={z_total:.2f}")

# Per-segment Z-scores
print(f"\n{'Segment':>10s}  {'Obs':>4s}  {'Exp':>6s}  {'Std':>5s}  {'Z':>6s}  {'Enrich':>7s}")
print("-" * 50)
segment_results = []
for i in range(N_SEGMENTS):
    obs = observed_per_seg[i]
    exp = expected_per_seg[i]
    std = std_per_seg[i]
    z = (obs - exp) / std if std > 0 else 0
    enrich = obs / exp if exp > 0 else float('inf') if obs > 0 else 0
    seg_res = {
        "segment": i,
        "label": segment_label(i),
        "bearing_start": i * SEGMENT_WIDTH,
        "bearing_end": (i + 1) * SEGMENT_WIDTH,
        "observed": obs,
        "expected": round(float(exp), 2),
        "std": round(float(std), 2),
        "z_score": round(float(z), 2),
        "enrichment": round(float(enrich), 3),
        "all_sites_in_segment": all_per_seg[i],
    }
    segment_results.append(seg_res)
    print(f"  {segment_label(i):>8s}  {obs:>4d}  {exp:>6.1f}  {std:>5.1f}  {z:>6.2f}  {enrich:>7.2f}x")


# ============================================================
# IDENTIFY WHAT'S IN THE HOT SEGMENTS
# ============================================================
# Sort segments by Z-score to find drivers
hot_segments = sorted(segment_results, key=lambda s: s['z_score'], reverse=True)

print("\n" + "=" * 70)
print("TOP ENRICHED SEGMENTS")
print("=" * 70)
for seg in hot_segments[:5]:
    print(f"\n  {seg['label']} — Z={seg['z_score']:.2f}, observed={seg['observed']}, expected={seg['expected']:.1f}")
    # List the sites in this segment
    sites_in_seg = [s for s in near_circle if s['segment'] == seg['segment']]
    for s in sites_in_seg:
        print(f"    {s['site_name'][:60]:<60s}  {s['country']:<15s}  {s['continent']:<15s}  "
              f"({s['lat']:.2f}, {s['lon']:.2f})  {s['age_bp']:.0f} BP  d={s['dist_km']:.1f}km  "
              f"bearing={s['bearing']:.1f}°")


# ============================================================
# GEOGRAPHIC CONTEXT: Map segments to approximate regions
# ============================================================
print("\n" + "=" * 70)
print("GEOGRAPHIC CONTEXT PER SEGMENT")
print("=" * 70)

# For each segment, show where the great circle passes through
# by generating points on the circle at each segment center
segment_geo = []
for i in range(N_SEGMENTS):
    bearing_center = (i + 0.5) * SEGMENT_WIDTH
    # Point on the great circle at this bearing from the pole
    # Using forward geodesic from pole at distance = quarter circumference
    br = math.radians(bearing_center)
    d_r = QUARTER_CIRC / EARTH_R_KM  # angular distance in radians
    pole_lat_r = math.radians(POLE_LAT)
    pole_lon_r = math.radians(POLE_LON)

    lat2 = math.asin(math.sin(pole_lat_r) * math.cos(d_r) +
                      math.cos(pole_lat_r) * math.sin(d_r) * math.cos(br))
    lon2 = pole_lon_r + math.atan2(
        math.sin(br) * math.sin(d_r) * math.cos(pole_lat_r),
        math.cos(d_r) - math.sin(pole_lat_r) * math.sin(lat2)
    )
    lat2_deg = math.degrees(lat2)
    lon2_deg = math.degrees(lon2)
    # Normalize longitude
    if lon2_deg > 180:
        lon2_deg -= 360
    elif lon2_deg < -180:
        lon2_deg += 360

    segment_geo.append({
        "segment": i,
        "label": segment_label(i),
        "circle_lat": round(lat2_deg, 2),
        "circle_lon": round(lon2_deg, 2),
    })
    print(f"  {segment_label(i):>10s}  circle passes through ({lat2_deg:>7.2f}, {lon2_deg:>8.2f})")


# ============================================================
# CONTINENT BREAKDOWN OF NEAR-CIRCLE SITES
# ============================================================
print("\n" + "=" * 70)
print("CONTINENT BREAKDOWN OF NEAR-CIRCLE SITES")
print("=" * 70)
cont_counts = Counter(s['continent'] for s in near_circle)
for cont, n in cont_counts.most_common():
    print(f"  {cont:<20s}  {n}")


# ============================================================
# SAVE RESULTS
# ============================================================
output = {
    "meta": {
        "analysis": "Early Holocene Arc Segment Decomposition",
        "date": "2026-03-20",
        "age_range_bp": [AGE_MIN_BP, AGE_MAX_BP],
        "threshold_km": THRESHOLD_KM,
        "n_segments": N_SEGMENTS,
        "segment_width_deg": SEGMENT_WIDTH,
        "monte_carlo_trials": N_TRIALS,
        "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
    },
    "summary": {
        "total_early_holocene_sites": len(early_holo),
        "near_circle_sites": total_observed,
        "expected_near_circle": round(float(expected_total), 2),
        "std_near_circle": round(float(std_total), 2),
        "z_score_total": round(float(z_total), 2),
        "fraction_near_circle": round(total_observed / len(early_holo), 5),
    },
    "segments": segment_results,
    "segment_geographic_context": segment_geo,
    "near_circle_sites": sorted(near_circle, key=lambda s: s['bearing']),
    "continent_breakdown": dict(cont_counts),
    "hot_segments_ranked": [
        {"label": s['label'], "z_score": s['z_score'], "observed": s['observed'],
         "expected": s['expected']}
        for s in hot_segments
    ],
}

out_path = os.path.join(OUT_DIR, "early_holocene_arc_segments.json")
with open(out_path, "w") as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nSaved: {out_path}")

results_path = os.path.join(BASE_DIR, "results", "early_holocene_arc_segments.json")
with open(results_path, "w") as f:
    json.dump(output, f, indent=2, default=str)
print(f"Saved: {results_path}")

print("\nDONE.")
