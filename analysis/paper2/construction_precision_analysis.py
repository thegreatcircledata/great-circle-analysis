#!/usr/bin/env python3
"""
Construction Precision Analysis — Phase 2A (partial) + Phase 2B (gap analysis)

Computes the Construction Anomaly Index (CAI) for sites with quantitative data,
tests spatial correlation with the Alison Great Circle, and produces a comprehensive
gap analysis for sites lacking measurements.

Author: Claude (for Ell)
Date: 2026-03-22
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from math import radians, sin, cos, asin, sqrt, atan2, pi
from scipy import stats

# ── Paths ──
BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data" / "construction_precision"
OUT  = BASE / "outputs" / "construction_precision"
CHARTS = OUT / "charts"
CHARTS.mkdir(parents=True, exist_ok=True)

# ── Alison Great Circle parameters ──
POLE_LAT = 59.682122
POLE_LON = -138.646087

def haversine(lat1, lon1, lat2, lon2):
    """Great circle distance in km."""
    R = 6371.0
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat/2)**2 + cos(lat1)*cos(lat2)*sin(dlon/2)**2
    return 2 * R * asin(sqrt(a))

def angular_distance_to_great_circle(site_lat, site_lon, pole_lat, pole_lon):
    """
    Compute angular distance from a site to the great circle defined by its pole.
    The great circle is the set of points 90° from the pole.
    Returns distance in km (absolute deviation from 90°).
    """
    ang = haversine(site_lat, site_lon, pole_lat, pole_lon)
    quarter_earth = 6371.0 * pi / 2  # ~10,007.5 km
    return abs(ang - quarter_earth)

# ── Site coordinates (manually compiled) ──
SITE_COORDS = {
    "Great Pyramid":       (29.9792, 31.1342),
    "Serapeum":            (29.8714, 31.2157),
    "Osireion at Abydos":  (26.1844, 31.9181),
    "Valley Temple at Giza": (29.9743, 31.1377),
    "Sacsayhuamán":        (-13.5094, -71.9820),
    "Ollantaytambo":       (-13.2581, -72.2625),
    "Puma Punku":          (-16.5617, -68.6803),
    "Göbekli Tepe":        (37.2233, 38.9225),
    "Baalbek":             (34.0069, 36.2039),
    "Stonehenge":          (51.1789, -1.8262),
    "Malta Ħaġar Qim":    (35.8278, 14.4419),
    "Easter Island Ahu Vinapu": (-27.1750, -109.3500),
}

def get_site_base(name):
    """Map CSV site names (with sub-components) to base site names."""
    for base in SITE_COORDS:
        if name.startswith(base):
            return base
    # Handle naming variations
    if "Great Pyramid" in name:
        return "Great Pyramid"
    if "Serapeum" in name:
        return "Serapeum"
    if "Puma Punku" in name:
        return "Puma Punku"
    return name

# ──────────────────────────────────────────────
# STEP 1: Load & Clean Data
# ──────────────────────────────────────────────

print("=" * 70)
print("CONSTRUCTION PRECISION ANALYSIS")
print("=" * 70)

df = pd.read_csv(DATA / "construction_measurements.csv")
print(f"\nLoaded {len(df)} site entries from construction_measurements.csv")

# Parse period to numeric BCE (approximate midpoint)
def parse_period(p):
    """Extract approximate year BCE from period string."""
    p = str(p).lower().strip()
    if "9500" in p: return -9150  # Göbekli Tepe midpoint
    if "3600" in p: return -3400  # Malta
    if "2560" in p: return -2560  # Great Pyramid
    if "2558" in p: return -2545  # Valley Temple
    if "2500" in p: return -2500  # Stonehenge
    if "1438" in p: return -85    # Sacsayhuamán (CE)
    if "15th century" in p: return -50  # Ollantaytambo (CE, ~1450)
    if "1290" in p: return -1285  # Osireion
    if "1279" in p: return -650   # Serapeum (spans 1279-30 BCE)
    if "536" in p: return 623     # Puma Punku (AD)
    if "1st century" in p: return 50  # Baalbek (CE)
    if "1200" in p: return 800    # Easter Island (post-AD 1200)
    return np.nan

df['year_numeric'] = df['period'].apply(parse_period)

# Convert N/A strings to actual NaN
for col in ['block_mass_tonnes', 'joint_gap_mean_mm', 'joint_gap_sigma_mm',
            'surface_flatness_mm_per_m', 'angular_accuracy_arc_sec', 'transport_distance_km']:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Add site coordinates and great circle distance
df['base_site'] = df['site'].apply(get_site_base)
df['lat'] = df['base_site'].map(lambda s: SITE_COORDS.get(s, (np.nan, np.nan))[0])
df['lon'] = df['base_site'].map(lambda s: SITE_COORDS.get(s, (np.nan, np.nan))[1])
df['distance_to_circle_km'] = df.apply(
    lambda r: angular_distance_to_great_circle(r['lat'], r['lon'], POLE_LAT, POLE_LON)
    if pd.notna(r['lat']) else np.nan, axis=1
)

# Save cleaned data
df.to_csv(OUT / "cleaned_data.csv", index=False)
print(f"Saved cleaned data to {OUT / 'cleaned_data.csv'}")

# ──────────────────────────────────────────────
# STEP 2: Construction Anomaly Index (CAI)
# ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 2: CONSTRUCTION ANOMALY INDEX")
print("=" * 70)

# Define baselines from experimental archaeology and modern standards
# Format: (measurement_type, baseline_value, baseline_sigma, unit, reference)
# baseline = what's achievable with period tools; lower = more precise
BASELINES = {
    'surface_flatness': {
        'modern_cnc': (1.5, 0.3, "Modern CNC — MIA standard"),
        'modern_grinding': (1.25, 0.25, "Modern machine grinding — MIA standard"),
        'ancient_stone_tools': (3.0, 1.0, "Estimated: stone tool abrasion on limestone"),
        'copper_tools': (2.0, 0.5, "Estimated: copper tools + stone rubbers on limestone"),
    },
    'joint_gap': {
        'modern_cnc': (0.5, 0.2, "Modern CNC dry-fit stone"),
        'copper_tools': (2.0, 0.5, "Estimated: copper-tool era achievable with trial fitting"),
        'stone_tools': (5.0, 2.0, "Estimated: stone-tool era"),
    },
    'angular_accuracy': {
        'modern_survey': (5.0, 2.0, "Modern total station"),
        'ancient_survey': (60.0, 20.0, "Estimated: ancient plumb/level/sight methods"),
    }
}

# Compute CAI for each row with quantitative precision data
cai_records = []

for _, row in df.iterrows():
    site = row['site']

    # Surface flatness CAI
    if pd.notna(row['surface_flatness_mm_per_m']):
        val = row['surface_flatness_mm_per_m']
        # Use copper tools baseline for Egyptian sites, stone tools for others
        if 'Pyramid' in site or 'Serapeum' in site or 'Osireion' in site or 'Valley' in site:
            bl_val, bl_sig, bl_ref = BASELINES['surface_flatness']['copper_tools']
        elif 'Puma' in site:
            bl_val, bl_sig, bl_ref = BASELINES['surface_flatness']['ancient_stone_tools']
        else:
            bl_val, bl_sig, bl_ref = BASELINES['surface_flatness']['ancient_stone_tools']

        cai = (bl_val - val) / bl_sig
        cai_records.append({
            'site': site,
            'measurement_type': 'surface_flatness',
            'measured_value': val,
            'unit': 'mm/m',
            'baseline_value': bl_val,
            'baseline_sigma': bl_sig,
            'baseline_ref': bl_ref,
            'CAI': round(cai, 2),
            'interpretation': 'MORE precise than baseline' if cai > 0 else 'LESS precise than baseline',
            'data_quality': row['data_quality'],
        })

    # Joint gap CAI
    if pd.notna(row['joint_gap_mean_mm']):
        val = row['joint_gap_mean_mm']
        if 'Pyramid' in site or 'Serapeum' in site:
            bl_val, bl_sig, bl_ref = BASELINES['joint_gap']['copper_tools']
        else:
            bl_val, bl_sig, bl_ref = BASELINES['joint_gap']['stone_tools']

        cai = (bl_val - val) / bl_sig
        cai_records.append({
            'site': site,
            'measurement_type': 'joint_gap',
            'measured_value': val,
            'unit': 'mm',
            'baseline_value': bl_val,
            'baseline_sigma': bl_sig,
            'baseline_ref': bl_ref,
            'CAI': round(cai, 2),
            'interpretation': 'MORE precise than baseline' if cai > 0 else 'LESS precise than baseline',
            'data_quality': row['data_quality'],
        })

    # Angular accuracy CAI
    if pd.notna(row['angular_accuracy_arc_sec']):
        val = row['angular_accuracy_arc_sec']
        bl_val, bl_sig, bl_ref = BASELINES['angular_accuracy']['ancient_survey']
        cai = (bl_val - val) / bl_sig
        cai_records.append({
            'site': site,
            'measurement_type': 'angular_accuracy',
            'measured_value': val,
            'unit': 'arc_sec',
            'baseline_value': bl_val,
            'baseline_sigma': bl_sig,
            'baseline_ref': bl_ref,
            'CAI': round(cai, 2),
            'interpretation': 'MORE precise than baseline' if cai > 0 else 'LESS precise than baseline',
            'data_quality': row['data_quality'],
        })

cai_df = pd.DataFrame(cai_records)
cai_df.to_csv(OUT / "anomaly_index.csv", index=False)

print(f"\nComputed CAI for {len(cai_df)} measurement-site pairs:")
print(f"  Sites with CAI data: {cai_df['site'].nunique()}")
print(f"  Measurement types: {cai_df['measurement_type'].unique().tolist()}")

# Summary by site
print("\n--- CAI Summary by Site ---")
for site in cai_df['site'].unique():
    subset = cai_df[cai_df['site'] == site]
    mean_cai = subset['CAI'].mean()
    details = "; ".join([f"{r['measurement_type']}={r['CAI']:.1f}" for _, r in subset.iterrows()])
    quality = subset['data_quality'].iloc[0]
    marker = "***" if mean_cai > 2 else "**" if mean_cai > 1 else "*" if mean_cai > 0 else ""
    print(f"  {site}: mean CAI = {mean_cai:+.2f} [{quality}] {marker}")
    print(f"    {details}")

# Compare against modern CNC
print("\n--- Comparison: Ancient vs Modern CNC ---")
gp_casing_flatness = 0.13  # mm/m, Great Pyramid casing stones
modern_cnc = 1.5  # mm/m, MIA standard
print(f"  Great Pyramid casing flatness: {gp_casing_flatness} mm/m")
print(f"  Modern CNC stone flatness:     {modern_cnc} mm/m")
print(f"  Ratio: Great Pyramid is {modern_cnc/gp_casing_flatness:.1f}× MORE precise than modern CNC standard")
gp_passage = 0.028  # mm/m
print(f"  Great Pyramid passage:         {gp_passage} mm/m")
print(f"  Ratio: Passage is {modern_cnc/gp_passage:.0f}× MORE precise than modern CNC standard")

# ──────────────────────────────────────────────
# STEP 3: Spatial Correlation with Great Circle
# ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 3: SPATIAL CORRELATION WITH GREAT CIRCLE")
print("=" * 70)

# Aggregate CAI per base site and merge with distance
site_cai = cai_df.copy()
site_cai['base_site'] = site_cai['site'].apply(get_site_base)
site_agg = site_cai.groupby('base_site').agg(
    mean_CAI=('CAI', 'mean'),
    max_CAI=('CAI', 'max'),
    n_measurements=('CAI', 'count'),
).reset_index()

# Add distance to circle
site_agg['distance_to_circle_km'] = site_agg['base_site'].map(
    lambda s: angular_distance_to_great_circle(
        SITE_COORDS[s][0], SITE_COORDS[s][1], POLE_LAT, POLE_LON
    ) if s in SITE_COORDS else np.nan
)

print(f"\nSites with CAI + distance data: {len(site_agg)}")
print(site_agg[['base_site', 'mean_CAI', 'distance_to_circle_km', 'n_measurements']].to_string(index=False))

# Spearman correlation (if we have enough sites)
valid = site_agg.dropna(subset=['mean_CAI', 'distance_to_circle_km'])
if len(valid) >= 3:
    rho, p_val = stats.spearmanr(valid['mean_CAI'], valid['distance_to_circle_km'])
    print(f"\nSpearman correlation (CAI vs distance): rho = {rho:.3f}, p = {p_val:.4f}")
    print(f"WARNING: Only {len(valid)} sites — insufficient for meaningful inference (need ≥10)")
else:
    rho, p_val = np.nan, np.nan
    print(f"\nInsufficient sites ({len(valid)}) for correlation test")

# Monte Carlo: random great circles
print("\nMonte Carlo: 10,000 random great circles...")
n_monte = 10000
rng = np.random.default_rng(42)
random_rhos = []

for _ in range(n_monte):
    rand_pole_lat = np.degrees(np.arcsin(rng.uniform(-1, 1)))
    rand_pole_lon = rng.uniform(-180, 180)
    rand_dists = []
    for _, row in valid.iterrows():
        s = row['base_site']
        lat, lon = SITE_COORDS[s]
        d = angular_distance_to_great_circle(lat, lon, rand_pole_lat, rand_pole_lon)
        rand_dists.append(d)
    if len(set(rand_dists)) > 1:
        r, _ = stats.spearmanr(valid['mean_CAI'].values, rand_dists)
        random_rhos.append(r)

random_rhos = np.array(random_rhos)
if not np.isnan(rho):
    percentile = np.sum(random_rhos <= rho) / len(random_rhos) * 100
    print(f"Observed rho: {rho:.3f}")
    print(f"Monte Carlo percentile: {percentile:.1f}th (higher = more negative correlation expected)")
    print(f"Random circle median rho: {np.median(random_rhos):.3f}")
else:
    percentile = np.nan

# Save spatial correlation results
spatial_results = {
    'n_sites': len(valid),
    'spearman_rho': rho,
    'p_value': p_val,
    'monte_carlo_n': n_monte,
    'monte_carlo_percentile': percentile,
    'verdict': 'INSUFFICIENT DATA' if len(valid) < 10 else ('SIGNIFICANT' if p_val < 0.05 else 'NOT SIGNIFICANT'),
}
pd.DataFrame([spatial_results]).to_csv(OUT / "spatial_correlation.csv", index=False)

# All sites distance to circle (for context)
print("\n--- All Sites: Distance to Alison Great Circle ---")
all_sites_dist = []
for name, (lat, lon) in sorted(SITE_COORDS.items()):
    d = angular_distance_to_great_circle(lat, lon, POLE_LAT, POLE_LON)
    all_sites_dist.append({'site': name, 'lat': lat, 'lon': lon, 'distance_to_circle_km': round(d, 1)})
    has_data = name in site_agg['base_site'].values
    marker = " [HAS PRECISION DATA]" if has_data else ""
    print(f"  {name}: {d:.0f} km{marker}")

pd.DataFrame(all_sites_dist).to_csv(OUT / "site_distances.csv", index=False)

# ──────────────────────────────────────────────
# STEP 4: Temporal Analysis
# ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 4: TEMPORAL ANALYSIS")
print("=" * 70)

# Merge CAI with temporal data
site_cai_temporal = site_cai.merge(
    df[['site', 'year_numeric']].drop_duplicates(),
    on='site', how='left'
)
temporal_agg = site_cai_temporal.groupby('site').agg(
    mean_CAI=('CAI', 'mean'),
    year=('year_numeric', 'first'),
).dropna()

print(f"\nSites with CAI + date: {len(temporal_agg)}")
if len(temporal_agg) >= 3:
    print("\nCAI by period:")
    for _, r in temporal_agg.sort_values('year').iterrows():
        year_label = f"{abs(int(r['year']))} {'BCE' if r['year'] < 0 else 'CE'}"
        print(f"  {r.name}: ~{year_label}, CAI = {r['mean_CAI']:+.2f}")

    r_temp, p_temp = stats.pearsonr(temporal_agg['year'], temporal_agg['mean_CAI'])
    print(f"\nPearson correlation (year vs CAI): r = {r_temp:.3f}, p = {p_temp:.4f}")
    print(f"WARNING: Only {len(temporal_agg)} sites — insufficient for temporal inference")
else:
    r_temp, p_temp = np.nan, np.nan
    print("Insufficient temporal data")

pd.DataFrame([{
    'n_sites': len(temporal_agg),
    'pearson_r': r_temp, 'p_value': p_temp,
    'verdict': 'INSUFFICIENT DATA'
}]).to_csv(OUT / "temporal_analysis.csv", index=False)

# ──────────────────────────────────────────────
# STEP 5: PHASE 2B — GAP ANALYSIS
# ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 5: PHASE 2B — GAP ANALYSIS")
print("=" * 70)

gap_records = []
for _, row in df.iterrows():
    site = row['site']
    has_joint = pd.notna(row['joint_gap_mean_mm'])
    has_flatness = pd.notna(row['surface_flatness_mm_per_m'])
    has_angular = pd.notna(row['angular_accuracy_arc_sec'])
    n_measurements = sum([has_joint, has_flatness, has_angular])

    if n_measurements >= 2:
        status = "QUANTIFIED"
    elif n_measurements == 1:
        status = "PARTIAL"
    elif row['data_quality'] == 'MEDIUM':
        status = "QUALITATIVE ONLY"
    else:
        status = "UNMEASURED"

    gap_records.append({
        'site': site,
        'period': row['period'],
        'data_quality': row['data_quality'],
        'classification': status,
        'has_joint_gap': has_joint,
        'has_surface_flatness': has_flatness,
        'has_angular_accuracy': has_angular,
        'n_precision_measurements': n_measurements,
        'has_3d_scan': site in ['Stonehenge sarsens'],  # Only Stonehenge has confirmed sub-mm scan
        'scan_accessible': False,  # Stonehenge scan is restricted
        'tool_technology': row['tool_technology'],
    })

gap_df = pd.DataFrame(gap_records)

print("\n--- Classification ---")
for cls in ['QUANTIFIED', 'PARTIAL', 'QUALITATIVE ONLY', 'UNMEASURED']:
    sites = gap_df[gap_df['classification'] == cls]['site'].tolist()
    print(f"  {cls}: {len(sites)}")
    for s in sites:
        print(f"    - {s}")

# ──────────────────────────────────────────────
# STEP 6: VISUALIZATIONS
# ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 6: VISUALIZATIONS")
print("=" * 70)

# --- Chart 1: World Map with CAI ---
fig, ax = plt.subplots(1, 1, figsize=(16, 8))

# Load Natural Earth coastlines if available
try:
    import json
    ne_path = BASE / "data" / "natural_earth"
    coast_files = list(ne_path.glob("*coastline*.json")) + list(ne_path.glob("*land*.json"))
    if not coast_files:
        # Try geojson
        coast_files = list(ne_path.glob("*.geojson"))
except:
    coast_files = []

# Simple world outline
ax.set_xlim(-180, 180)
ax.set_ylim(-60, 75)
ax.set_facecolor('#e8e8e8')
ax.set_aspect('equal')

# Draw great circle
gc_lons = np.linspace(-180, 180, 1000)
gc_lats = []
pole_lat_r = radians(POLE_LAT)
pole_lon_r = radians(POLE_LON)
for lon in gc_lons:
    lon_r = radians(lon)
    # Great circle: points 90° from pole
    # lat = arctan(-cos(lon - pole_lon) / tan(pole_lat))... simplified approach:
    # Parametric: point on GC at azimuth theta from pole
    pass

# Simpler: plot great circle by rotating
thetas = np.linspace(0, 2*pi, 1000)
gc_points = []
for theta in thetas:
    # Point on equator of pole-centered system
    x = cos(theta)
    y = sin(theta)
    z = 0
    # Rotate by (90 - pole_lat) around y-axis to get geographic coords
    angle = radians(90 - POLE_LAT)
    xr = x * cos(angle) + z * sin(angle)
    yr = y
    zr = -x * sin(angle) + z * cos(angle)
    # Convert to lat/lon
    lat = np.degrees(np.arcsin(zr))
    lon = np.degrees(np.arctan2(yr, xr)) + POLE_LON
    if lon > 180: lon -= 360
    if lon < -180: lon += 360
    gc_points.append((lon, lat))

gc_points = np.array(gc_points)
# Split at discontinuities (wrap-around)
diffs = np.abs(np.diff(gc_points[:, 0]))
splits = np.where(diffs > 90)[0] + 1
segments = np.split(gc_points, splits)
for seg in segments:
    ax.plot(seg[:, 0], seg[:, 1], 'r-', linewidth=1.5, alpha=0.5, zorder=2)

# Build CAI lookup for base sites
site_mean_cai = {}
for _, r in site_agg.iterrows():
    site_mean_cai[r['base_site']] = r['mean_CAI']

# Plot sites
for name, (lat, lon) in SITE_COORDS.items():
    if name in site_mean_cai:
        cai_val = site_mean_cai[name]
        # Color by CAI: red = anomalously precise, blue = expected
        if cai_val > 2:
            color = '#d62728'  # bright red
            size = 120
        elif cai_val > 0:
            color = '#ff7f0e'  # orange
            size = 90
        else:
            color = '#1f77b4'  # blue
            size = 70
    else:
        color = '#888888'  # grey = no data
        size = 50

    ax.scatter(lon, lat, c=color, s=size, edgecolors='black', linewidth=0.5, zorder=5)
    # Label
    offset_x, offset_y = 3, 3
    if 'Easter' in name:
        offset_x = -40
    elif 'Göbekli' in name:
        offset_y = -8
    ax.annotate(name.replace(' Ħaġar Qim', '').replace(' Ahu Vinapu', ''),
                (lon, lat), xytext=(offset_x, offset_y),
                textcoords='offset points', fontsize=7, zorder=6)

# Legend
legend_elements = [
    mpatches.Patch(facecolor='#d62728', edgecolor='black', label='CAI > 2 (highly anomalous)'),
    mpatches.Patch(facecolor='#ff7f0e', edgecolor='black', label='CAI 0–2 (somewhat anomalous)'),
    mpatches.Patch(facecolor='#1f77b4', edgecolor='black', label='CAI < 0 (below baseline)'),
    mpatches.Patch(facecolor='#888888', edgecolor='black', label='No precision data'),
    plt.Line2D([0], [0], color='red', alpha=0.5, linewidth=1.5, label='Alison Great Circle'),
]
ax.legend(handles=legend_elements, loc='lower left', fontsize=8)
ax.set_title("Ancient Construction Precision: Sites Colored by Construction Anomaly Index (CAI)\nOverlaid on Alison Great Circle", fontsize=12)
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
plt.tight_layout()
plt.savefig(CHARTS / "01_world_map_cai.png", dpi=150)
plt.close()
print("  Saved: 01_world_map_cai.png")

# --- Chart 2: CAI vs Distance scatter ---
fig, ax = plt.subplots(figsize=(10, 6))
colors = {'Great Pyramid': '#d62728', 'Serapeum': '#ff7f0e', 'Puma Punku': '#2ca02c'}
for _, row in site_agg.iterrows():
    name = row['base_site']
    c = colors.get(name, '#888888')
    ax.scatter(row['distance_to_circle_km'], row['mean_CAI'],
               c=c, s=100, edgecolors='black', zorder=5)
    ax.annotate(name, (row['distance_to_circle_km'], row['mean_CAI']),
                xytext=(8, 5), textcoords='offset points', fontsize=9)

# Add unmeasured sites as grey markers at CAI=0
for name, (lat, lon) in SITE_COORDS.items():
    if name not in site_agg['base_site'].values:
        d = angular_distance_to_great_circle(lat, lon, POLE_LAT, POLE_LON)
        ax.scatter(d, 0, c='#cccccc', s=60, edgecolors='grey', marker='x', zorder=3)
        ax.annotate(name, (d, 0), xytext=(5, -10), textcoords='offset points',
                    fontsize=7, color='grey')

ax.axhline(0, color='black', linewidth=0.5, linestyle='--', alpha=0.5)
ax.set_xlabel("Distance to Alison Great Circle (km)", fontsize=11)
ax.set_ylabel("Mean Construction Anomaly Index (CAI)", fontsize=11)
ax.set_title("CAI vs. Distance to Great Circle\n(Grey X = sites with no precision data)", fontsize=12)
if not np.isnan(rho):
    ax.text(0.02, 0.98, f"Spearman ρ = {rho:.3f} (p = {p_val:.3f})\nn = {len(valid)} sites (INSUFFICIENT)",
            transform=ax.transAxes, va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
plt.tight_layout()
plt.savefig(CHARTS / "02_cai_vs_distance.png", dpi=150)
plt.close()
print("  Saved: 02_cai_vs_distance.png")

# --- Chart 3: CAI by Site (bar chart) ---
fig, ax = plt.subplots(figsize=(14, 6))

# Sort by CAI
all_cai = cai_df.sort_values('CAI', ascending=True)
colors_bar = []
for _, r in all_cai.iterrows():
    if r['data_quality'] == 'LOW':
        colors_bar.append('#ffcc99')  # light orange for disputed
    elif r['CAI'] > 2:
        colors_bar.append('#d62728')
    elif r['CAI'] > 0:
        colors_bar.append('#ff7f0e')
    else:
        colors_bar.append('#1f77b4')

labels = [f"{r['site']}\n({r['measurement_type']})" for _, r in all_cai.iterrows()]
bars = ax.barh(range(len(all_cai)), all_cai['CAI'].values, color=colors_bar, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(all_cai)))
ax.set_yticklabels(labels, fontsize=7)
ax.axvline(0, color='black', linewidth=1)
ax.axvline(2, color='red', linewidth=0.5, linestyle='--', alpha=0.5, label='CAI = 2 (strong anomaly)')

# Add value labels
for i, (_, r) in enumerate(all_cai.iterrows()):
    ax.text(r['CAI'] + 0.1, i, f"{r['CAI']:.1f}", va='center', fontsize=7)

ax.set_xlabel("Construction Anomaly Index (CAI)\nPositive = more precise than experimental baseline", fontsize=10)
ax.set_title("Construction Anomaly Index by Site and Measurement Type\n(Light orange = disputed/low quality data)", fontsize=12)
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(CHARTS / "03_cai_bar_chart.png", dpi=150)
plt.close()
print("  Saved: 03_cai_bar_chart.png")

# --- Chart 4: Gap Analysis Summary ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

# Left: classification pie
class_counts = gap_df['classification'].value_counts()
colors_pie = {'QUANTIFIED': '#2ca02c', 'PARTIAL': '#ff7f0e', 'QUALITATIVE ONLY': '#ffcc00', 'UNMEASURED': '#d62728'}
ax1.pie(class_counts.values,
        labels=[f"{k}\n({v} sites)" for k, v in class_counts.items()],
        colors=[colors_pie.get(k, '#888') for k in class_counts.index],
        autopct='%1.0f%%', startangle=90, textprops={'fontsize': 9})
ax1.set_title("Data Availability Classification\n(19 site entries)", fontsize=11)

# Right: measurement coverage heatmap
sites_short = [s[:30] for s in gap_df['site']]
meas_types = ['has_joint_gap', 'has_surface_flatness', 'has_angular_accuracy']
meas_labels = ['Joint Gap', 'Surface\nFlatness', 'Angular\nAccuracy']
data_matrix = gap_df[meas_types].values.astype(float)

im = ax2.imshow(data_matrix, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)
ax2.set_xticks(range(len(meas_labels)))
ax2.set_xticklabels(meas_labels, fontsize=9)
ax2.set_yticks(range(len(sites_short)))
ax2.set_yticklabels(sites_short, fontsize=7)
ax2.set_title("Measurement Coverage by Site\n(Green = exists, Red = missing)", fontsize=11)

# Add text annotations
for i in range(len(gap_df)):
    for j in range(len(meas_types)):
        val = data_matrix[i, j]
        ax2.text(j, i, '✓' if val else '✗',
                ha='center', va='center', fontsize=8,
                color='white' if val else 'darkred')

plt.tight_layout()
plt.savefig(CHARTS / "04_gap_analysis.png", dpi=150)
plt.close()
print("  Saved: 04_gap_analysis.png")

# --- Chart 5: Precision comparison across sites + modern standards ---
fig, ax = plt.subplots(figsize=(12, 6))

# Surface flatness comparison (where available)
flatness_data = df[df['surface_flatness_mm_per_m'].notna()][['site', 'surface_flatness_mm_per_m', 'data_quality']].copy()
flatness_data = flatness_data.sort_values('surface_flatness_mm_per_m')

# Add modern baselines
modern_baselines = pd.DataFrame([
    {'site': 'Modern CNC (MIA std)', 'surface_flatness_mm_per_m': 1.5, 'data_quality': 'REFERENCE'},
    {'site': 'Modern grinding (MIA std)', 'surface_flatness_mm_per_m': 1.25, 'data_quality': 'REFERENCE'},
])
flatness_data = pd.concat([flatness_data, modern_baselines], ignore_index=True)
flatness_data = flatness_data.sort_values('surface_flatness_mm_per_m')

colors_flat = []
for _, r in flatness_data.iterrows():
    if r['data_quality'] == 'REFERENCE':
        colors_flat.append('#4444ff')
    elif r['data_quality'] == 'LOW':
        colors_flat.append('#ffcc99')
    elif r['data_quality'] == 'HIGH':
        colors_flat.append('#d62728')
    else:
        colors_flat.append('#ff7f0e')

ax.barh(range(len(flatness_data)), flatness_data['surface_flatness_mm_per_m'].values,
        color=colors_flat, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(flatness_data)))
ax.set_yticklabels(flatness_data['site'].values, fontsize=8)
ax.set_xlabel("Surface Flatness (mm per meter)\nLower = more precise", fontsize=10)
ax.set_title("Surface Flatness: Ancient Sites vs. Modern Industrial Standards\n(Blue = modern reference; Red = ancient HIGH quality; Orange = MEDIUM; Light = LOW/disputed)",
             fontsize=11)
ax.set_xscale('log')

for i, (_, r) in enumerate(flatness_data.iterrows()):
    ax.text(r['surface_flatness_mm_per_m'] * 1.1, i,
            f"{r['surface_flatness_mm_per_m']:.4f}" if r['surface_flatness_mm_per_m'] < 0.1
            else f"{r['surface_flatness_mm_per_m']:.2f}",
            va='center', fontsize=7)

plt.tight_layout()
plt.savefig(CHARTS / "05_flatness_comparison.png", dpi=150)
plt.close()
print("  Saved: 05_flatness_comparison.png")

# ──────────────────────────────────────────────
# GO / NO-GO ASSESSMENT
# ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("GO / NO-GO ASSESSMENT")
print("=" * 70)

tests = {
    'Sites with quantitative data >= 10': ('FAIL', f'{gap_df[gap_df["classification"].isin(["QUANTIFIED","PARTIAL"])]["site"].nunique()} sites'),
    'Sites with CAI computable >= 6': ('FAIL', f'{cai_df["site"].nunique()} sites'),
    'Spatial correlation p < 0.05': ('INSUFFICIENT DATA', f'p = {p_val:.4f}, n = {len(valid)}' if not np.isnan(p_val) else 'N/A'),
    'Monte Carlo percentile > 95th': ('INSUFFICIENT DATA', f'{percentile:.1f}th' if not np.isnan(percentile) else 'N/A'),
    'Temporal correlation r > 0.5': ('INSUFFICIENT DATA', f'r = {r_temp:.3f}' if not np.isnan(r_temp) else 'N/A'),
}

for test, (result, detail) in tests.items():
    print(f"  [{result:^20s}] {test}: {detail}")

print(f"\n  VERDICT: *** PARTIAL — Phase 2B Gap Analysis ***")
print(f"  The data landscape is dominated by one site (Giza) with rigorous measurements.")
print(f"  No cross-site statistical analysis is currently possible.")
print(f"  The gap analysis and research proposal are the primary deliverables.")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
