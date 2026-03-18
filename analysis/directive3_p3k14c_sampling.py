#!/usr/bin/env python3
"""
Directive 3: p3k14c Sampling Intensity Control
================================================
Tests whether the p3k14c radiocarbon signal (Z=7.21, driven by South America
Z=9.17) reflects genuine site clustering or Peruvian sampling density.

Sub-tests:
  3A: On-line vs off-line sampling density within South America
  3B: Grid-cell sampling intensity normalization
  3C: Within-Peru control
  3D: Africa and Asia without South America
"""

import csv, math, random, json, os, time
import numpy as np
from scipy import stats as sp_stats
from collections import defaultdict

# ============================================================
# CONFIGURATION
# ============================================================
POLE_LAT = 59.682122
POLE_LNG = -138.646087
EARTH_R_KM = 6371
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
THRESHOLD_KM = 50
N_MC_TRIALS = 200
GRID_DEG = 1  # 1° × 1° grid cells

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "outputs", "p3k14c_sampling_control")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def haversine_km(lat1, lon1, lat2, lon2):
    R = EARTH_R_KM
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(min(1.0, a)))

def dist_from_gc(lat, lon, pole_lat=POLE_LAT, pole_lon=POLE_LNG):
    d = haversine_km(lat, lon, pole_lat, pole_lon)
    return abs(d - QUARTER_CIRC)

def rand_matched(site_lats, site_lons):
    lat = random.choice(site_lats) + random.gauss(0, 2)
    lon = random.choice(site_lons) + random.gauss(0, 2)
    return max(-90, min(90, lat)), max(-180, min(180, lon))

def mc_zscore(sites, threshold=THRESHOLD_KM, n_trials=N_MC_TRIALS):
    """Run distribution-matched Monte Carlo and return Z, obs, mu, sigma."""
    site_lats = [s[0] for s in sites]
    site_lons = [s[1] for s in sites]
    n = len(sites)
    observed = sum(1 for lat, lon in sites if dist_from_gc(lat, lon) <= threshold)

    rand_counts = []
    for _ in range(n_trials):
        pts = [rand_matched(site_lats, site_lons) for _ in range(n)]
        rand_counts.append(sum(1 for lat, lon in pts if dist_from_gc(lat, lon) <= threshold))

    mu = np.mean(rand_counts)
    sigma = np.std(rand_counts)
    z = (observed - mu) / sigma if sigma > 0 else 0
    p_val = 1 - sp_stats.norm.cdf(z) if z > 0 else 1.0
    return {
        "z_score": round(float(z), 2),
        "observed": int(observed),
        "expected": round(float(mu), 2),
        "std": round(float(sigma), 2),
        "p_value": float(p_val),
        "n_sites": n
    }

def grid_key(lat, lon, deg=GRID_DEG):
    return (int(math.floor(lat / deg)) * deg, int(math.floor(lon / deg)) * deg)

# ============================================================
# LOAD P3K14C DATA
# ============================================================
print("=" * 70)
print("LOADING P3K14C DATA")
print("=" * 70)

all_rows = []  # (lat, lon, siteid, continent, country)
csv_path = os.path.join(BASE_DIR, "p3k14c_data.csv")

with open(csv_path, encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            lat = float(row["Lat"])
            lon = float(row["Long"])
            if lat == 0 and lon == 0: continue
            siteid = row.get("SiteID", "").strip()
            continent = row.get("Continent", "").strip()
            country = row.get("Country", "").strip()
            all_rows.append((lat, lon, siteid, continent, country))
        except (ValueError, KeyError, TypeError):
            pass

print(f"Total rows: {len(all_rows)}")

# Deduplicate by SiteID for unique sites
unique_by_siteid = {}
for lat, lon, siteid, continent, country in all_rows:
    if siteid and siteid not in unique_by_siteid:
        unique_by_siteid[siteid] = (lat, lon, continent, country)

# Also deduplicate by coordinate for sites without SiteID
unique_by_coord = {}
for lat, lon, siteid, continent, country in all_rows:
    key = (round(lat, 4), round(lon, 4))
    if key not in unique_by_coord:
        unique_by_coord[key] = (lat, lon, continent, country)

print(f"Unique sites by SiteID: {len(unique_by_siteid)}")
print(f"Unique sites by coordinate: {len(unique_by_coord)}")

# Build unique site lists
all_unique_sites = [(lat, lon, continent, country) for lat, lon, continent, country in unique_by_coord.values()]

# Continent breakdown
continent_counts = defaultdict(int)
for _, _, cont, _ in all_unique_sites:
    continent_counts[cont] += 1
print(f"Continent breakdown: {dict(continent_counts)}")

# Identify South American and Peruvian sites
sa_rows = [(lat, lon, siteid, cont, country) for lat, lon, siteid, cont, country in all_rows if cont == "South America"]
sa_unique = [(lat, lon) for lat, lon, cont, country in all_unique_sites if cont == "South America"]
peru_rows = [(lat, lon, siteid, cont, country) for lat, lon, siteid, cont, country in all_rows if country == "Peru"]
peru_unique_coords = set()
peru_unique = []
for lat, lon, siteid, cont, country in peru_rows:
    key = (round(lat, 4), round(lon, 4))
    if key not in peru_unique_coords:
        peru_unique_coords.add(key)
        peru_unique.append((lat, lon))

print(f"South America: {len(sa_rows)} rows, {len(sa_unique)} unique sites")
print(f"Peru: {len(peru_rows)} rows, {len(peru_unique)} unique sites")

# ============================================================
# TEST 3A: ON-LINE vs OFF-LINE SAMPLING DENSITY (South America)
# ============================================================
print("\n" + "=" * 70)
print("TEST 3A: ON-LINE vs OFF-LINE SAMPLING DENSITY (South America)")
print("=" * 70)

# All SA rows: classify as on-line or off-line
sa_online_rows = []
sa_offline_rows = []
sa_online_siteids = set()
sa_offline_siteids = set()

for lat, lon, siteid, cont, country in sa_rows:
    d = dist_from_gc(lat, lon)
    if d <= THRESHOLD_KM:
        sa_online_rows.append((lat, lon, siteid))
        if siteid: sa_online_siteids.add(siteid)
    else:
        sa_offline_rows.append((lat, lon, siteid))
        if siteid: sa_offline_siteids.add(siteid)

# Unique sites on-line/off-line
sa_online_unique_coords = set()
sa_offline_unique_coords = set()
for lat, lon, siteid in sa_online_rows:
    sa_online_unique_coords.add((round(lat, 4), round(lon, 4)))
for lat, lon, siteid in sa_offline_rows:
    sa_offline_unique_coords.add((round(lat, 4), round(lon, 4)))

online_sites = len(sa_online_unique_coords)
offline_sites = len(sa_offline_unique_coords)
online_dates = len(sa_online_rows)
offline_dates = len(sa_offline_rows)

dates_per_site_online = online_dates / online_sites if online_sites > 0 else 0
dates_per_site_offline = offline_dates / offline_sites if offline_sites > 0 else 0

print(f"On-line (within {THRESHOLD_KM}km):")
print(f"  Unique sites: {online_sites}")
print(f"  Total C14 dates: {online_dates}")
print(f"  Dates per site: {dates_per_site_online:.2f}")
print(f"Off-line (>{THRESHOLD_KM}km):")
print(f"  Unique sites: {offline_sites}")
print(f"  Total C14 dates: {offline_dates}")
print(f"  Dates per site: {dates_per_site_offline:.2f}")

# Z-score using all rows vs unique sites only
all_sa_coords = [(lat, lon) for lat, lon, _, _, _ in sa_rows]
all_rows_z = mc_zscore(all_sa_coords)
unique_sites_z = mc_zscore(sa_unique)

print(f"\nZ-score (all rows): {all_rows_z['z_score']}")
print(f"Z-score (unique sites): {unique_sites_z['z_score']}")

sampling_density = {
    "online": {
        "unique_sites": online_sites,
        "total_c14_dates": online_dates,
        "dates_per_site": round(dates_per_site_online, 2)
    },
    "offline": {
        "unique_sites": offline_sites,
        "total_c14_dates": offline_dates,
        "dates_per_site": round(dates_per_site_offline, 2)
    },
    "ratio_online_offline": round(dates_per_site_online / dates_per_site_offline, 2) if dates_per_site_offline > 0 else None,
    "z_all_rows": all_rows_z,
    "z_unique_sites": unique_sites_z,
    "z_difference": round(all_rows_z["z_score"] - unique_sites_z["z_score"], 2),
    "diagnosis": ""
}

z_diff = abs(all_rows_z["z_score"] - unique_sites_z["z_score"])
if z_diff < 1:
    sampling_density["diagnosis"] = "CLEAR: All-rows and unique-site Z-scores are similar — sampling density is NOT inflating the signal"
elif z_diff < 3:
    sampling_density["diagnosis"] = "MODERATE: Some Z-score difference — sampling density may partially inflate the signal"
else:
    sampling_density["diagnosis"] = "CONCERN: Large Z-score difference — sampling density IS inflating the signal"

print(f"Diagnosis: {sampling_density['diagnosis']}")

with open(os.path.join(OUT_DIR, "sampling_density.json"), "w") as f:
    json.dump(sampling_density, f, indent=2)

# ============================================================
# TEST 3B: GRID-CELL SAMPLING INTENSITY
# ============================================================
print("\n" + "=" * 70)
print("TEST 3B: 1° GRID-CELL SAMPLING INTENSITY (South America)")
print("=" * 70)

# For each 1° cell in South America, count C14 dates
sa_grid_counts = defaultdict(int)
for lat, lon, siteid, cont, country in sa_rows:
    gk = grid_key(lat, lon)
    sa_grid_counts[gk] += 1

# Classify cells as near-circle or far-from-circle
near_circle_cells = {}
far_circle_cells = {}

for (glat, glon), count in sa_grid_counts.items():
    # Cell center
    clat = glat + GRID_DEG / 2
    clon = glon + GRID_DEG / 2
    d = dist_from_gc(clat, clon)
    if d <= THRESHOLD_KM:
        near_circle_cells[(glat, glon)] = count
    else:
        far_circle_cells[(glat, glon)] = count

# For latitude-matched comparison, group far cells by latitude band
near_lats = [k[0] for k in near_circle_cells]
lat_range = (min(near_lats) - 5, max(near_lats) + 5) if near_lats else (-60, 15)

# Filter far cells to similar latitudes
far_matched = {k: v for k, v in far_circle_cells.items()
               if lat_range[0] <= k[0] <= lat_range[1]}

near_counts = list(near_circle_cells.values())
far_counts = list(far_matched.values())

grid_analysis = {
    "grid_resolution_deg": GRID_DEG,
    "near_circle_cells": {
        "n_cells": len(near_circle_cells),
        "total_dates": sum(near_counts),
        "mean_dates_per_cell": round(float(np.mean(near_counts)), 2) if near_counts else 0,
        "median_dates_per_cell": round(float(np.median(near_counts)), 2) if near_counts else 0,
        "max_dates_in_cell": max(near_counts) if near_counts else 0,
    },
    "far_circle_cells_latmatched": {
        "n_cells": len(far_matched),
        "total_dates": sum(far_counts),
        "mean_dates_per_cell": round(float(np.mean(far_counts)), 2) if far_counts else 0,
        "median_dates_per_cell": round(float(np.median(far_counts)), 2) if far_counts else 0,
        "max_dates_in_cell": max(far_counts) if far_counts else 0,
        "lat_range": list(lat_range),
    },
    "ratio_near_far_mean": None,
    "mannwhitney_p": None,
    "diagnosis": ""
}

if near_counts and far_counts:
    grid_analysis["ratio_near_far_mean"] = round(
        np.mean(near_counts) / np.mean(far_counts), 2) if np.mean(far_counts) > 0 else None
    try:
        stat, p = sp_stats.mannwhitneyu(near_counts, far_counts, alternative='greater')
        grid_analysis["mannwhitney_p"] = float(p)
        if p < 0.05:
            grid_analysis["diagnosis"] = f"CONCERN: Near-circle cells have significantly higher C14 density (Mann-Whitney p={p:.4f})"
        else:
            grid_analysis["diagnosis"] = f"CLEAR: No significant difference in C14 density (Mann-Whitney p={p:.4f})"
    except:
        grid_analysis["diagnosis"] = "Could not compute Mann-Whitney test"

print(f"Near-circle cells: {grid_analysis['near_circle_cells']['n_cells']} "
      f"(mean {grid_analysis['near_circle_cells']['mean_dates_per_cell']} dates/cell)")
print(f"Far-circle cells (lat-matched): {grid_analysis['far_circle_cells_latmatched']['n_cells']} "
      f"(mean {grid_analysis['far_circle_cells_latmatched']['mean_dates_per_cell']} dates/cell)")
print(f"Diagnosis: {grid_analysis['diagnosis']}")

with open(os.path.join(OUT_DIR, "grid_analysis.json"), "w") as f:
    json.dump(grid_analysis, f, indent=2)

# ============================================================
# TEST 3C: WITHIN-PERU CONTROL
# ============================================================
print("\n" + "=" * 70)
print("TEST 3C: WITHIN-PERU CONTROL")
print("=" * 70)

if len(peru_unique) >= 10:
    peru_z = mc_zscore(peru_unique)
    print(f"Peru-only Z-score: {peru_z['z_score']} "
          f"(obs={peru_z['observed']}, exp={peru_z['expected']}, n={peru_z['n_sites']})")

    # On-line vs off-line within Peru
    peru_online = sum(1 for lat, lon in peru_unique if dist_from_gc(lat, lon) <= THRESHOLD_KM)
    peru_offline = len(peru_unique) - peru_online

    peru_control = {
        "peru_only_z": peru_z,
        "peru_total_unique_sites": len(peru_unique),
        "peru_online": peru_online,
        "peru_offline": peru_offline,
        "peru_online_fraction": round(peru_online / len(peru_unique), 4) if peru_unique else 0,
        "diagnosis": ""
    }

    if peru_z["z_score"] > 3:
        peru_control["diagnosis"] = (f"Signal SURVIVES within Peru (Z={peru_z['z_score']}). "
                                     f"The Great Circle enrichment exists WITHIN the most-sampled region.")
    elif peru_z["z_score"] > 2:
        peru_control["diagnosis"] = (f"Marginal signal within Peru (Z={peru_z['z_score']}). "
                                     f"Some enrichment but not conclusive.")
    else:
        peru_control["diagnosis"] = (f"No significant signal within Peru alone (Z={peru_z['z_score']}). "
                                     f"Signal is a Peru-vs-rest-of-South-America artifact.")
else:
    peru_control = {"diagnosis": f"Insufficient Peru sites ({len(peru_unique)}) for meaningful test"}

print(f"Diagnosis: {peru_control['diagnosis']}")

with open(os.path.join(OUT_DIR, "peru_control.json"), "w") as f:
    json.dump(peru_control, f, indent=2)

# ============================================================
# TEST 3D: SIGNAL WITHOUT SOUTH AMERICA
# ============================================================
print("\n" + "=" * 70)
print("TEST 3D: SIGNAL WITHOUT SOUTH AMERICA")
print("=" * 70)

non_sa_sites = [(lat, lon) for lat, lon, cont, country in all_unique_sites if cont != "South America"]
africa_sites = [(lat, lon) for lat, lon, cont, country in all_unique_sites if cont == "Africa"]
asia_sites = [(lat, lon) for lat, lon, cont, country in all_unique_sites if cont == "Asia"]

print(f"Non-SA sites: {len(non_sa_sites)}")
print(f"Africa sites: {len(africa_sites)}")
print(f"Asia sites: {len(asia_sites)}")

non_sa_z = mc_zscore(non_sa_sites) if len(non_sa_sites) >= 10 else None
africa_z = mc_zscore(africa_sites) if len(africa_sites) >= 10 else None
asia_z = mc_zscore(asia_sites) if len(asia_sites) >= 10 else None

non_sa_signal = {
    "non_sa": non_sa_z,
    "africa": africa_z,
    "asia": asia_z,
    "n_non_sa": len(non_sa_sites),
    "n_africa": len(africa_sites),
    "n_asia": len(asia_sites),
    "diagnosis": ""
}

if non_sa_z:
    print(f"Non-SA combined Z = {non_sa_z['z_score']}")
    if non_sa_z["z_score"] > 3:
        non_sa_signal["diagnosis"] = (f"Signal SURVIVES without South America (Z={non_sa_z['z_score']}). "
                                      f"Replication holds on independent geographic grounds.")
    else:
        non_sa_signal["diagnosis"] = (f"Signal does NOT survive without South America (Z={non_sa_z['z_score']}). "
                                      f"SA drives the replication.")

if africa_z:
    print(f"Africa Z = {africa_z['z_score']}")
if asia_z:
    print(f"Asia Z = {asia_z['z_score']}")
print(f"Diagnosis: {non_sa_signal['diagnosis']}")

with open(os.path.join(OUT_DIR, "non_sa_signal.json"), "w") as f:
    json.dump(non_sa_signal, f, indent=2)

# ============================================================
# COMBINED RESULTS
# ============================================================
print("\n" + "=" * 70)
print("COMBINED RESULTS")
print("=" * 70)

overall = {
    "test_3a_sampling_density": sampling_density,
    "test_3b_grid_analysis": grid_analysis,
    "test_3c_peru_control": peru_control,
    "test_3d_non_sa_signal": non_sa_signal,
}

# Overall diagnosis
diagnoses = []
if "CLEAR" in sampling_density.get("diagnosis", ""):
    diagnoses.append("3A: Sampling density is NOT inflating the signal")
else:
    diagnoses.append("3A: " + sampling_density.get("diagnosis", ""))

diagnoses.append("3B: " + grid_analysis.get("diagnosis", ""))
diagnoses.append("3C: " + peru_control.get("diagnosis", ""))
diagnoses.append("3D: " + non_sa_signal.get("diagnosis", ""))

overall["summary_diagnoses"] = diagnoses

# Write RESULTS.md
sa_z_all = all_rows_z["z_score"]
sa_z_uniq = unique_sites_z["z_score"]

results_md = f"""# Directive 3: p3k14c Sampling Intensity Control — Results

## Test 3A: On-line vs Off-line Sampling Density (South America)

| Metric | On-line (<{THRESHOLD_KM}km) | Off-line (>{THRESHOLD_KM}km) |
|--------|-----------|----------|
| Unique sites | {sampling_density['online']['unique_sites']} | {sampling_density['offline']['unique_sites']} |
| Total C14 dates | {sampling_density['online']['total_c14_dates']} | {sampling_density['offline']['total_c14_dates']} |
| Dates per site | {sampling_density['online']['dates_per_site']} | {sampling_density['offline']['dates_per_site']} |

On-line/off-line ratio: {sampling_density['ratio_online_offline']}

Z-score (all rows): {sa_z_all}
Z-score (unique sites only): {sa_z_uniq}
Difference: {sampling_density['z_difference']}

**Diagnosis:** {sampling_density['diagnosis']}

## Test 3B: Grid-Cell Sampling Intensity

Near-circle cells: {grid_analysis['near_circle_cells']['n_cells']} (mean {grid_analysis['near_circle_cells']['mean_dates_per_cell']} dates/cell)
Far-circle cells (lat-matched): {grid_analysis['far_circle_cells_latmatched']['n_cells']} (mean {grid_analysis['far_circle_cells_latmatched']['mean_dates_per_cell']} dates/cell)

**Diagnosis:** {grid_analysis['diagnosis']}

## Test 3C: Within-Peru Control

{peru_control.get('diagnosis', 'N/A')}

## Test 3D: Signal Without South America

Non-SA sites: {len(non_sa_sites)}
Non-SA Z = {non_sa_z['z_score'] if non_sa_z else 'N/A'}
Africa Z = {africa_z['z_score'] if africa_z else 'N/A'}
Asia Z = {asia_z['z_score'] if asia_z else 'N/A'}

**Diagnosis:** {non_sa_signal['diagnosis']}

## Summary
"""

for d in diagnoses:
    results_md += f"- {d}\n"

with open(os.path.join(OUT_DIR, "RESULTS.md"), "w") as f:
    f.write(results_md)

print(f"\nAll results saved to {OUT_DIR}/")
print("Done.")
