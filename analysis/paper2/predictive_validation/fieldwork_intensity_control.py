#!/usr/bin/env python3
"""
Fieldwork Intensity Control for Predictive Validation
======================================================
Test whether the p3k14c post-2001 great circle signal (Z=+16.74 at 50km)
survives after controlling for regional research intensity.

1. Density of post-2001 dates by 5°×5° grid cells
2. Near-circle enrichment ratio per grid cell
3. Remove top 5 most research-intensive cells, recompute Z
4. Continental decomposition of Z-scores
"""

import csv, json, math, os, sys, re, time
import numpy as np

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE = os.path.expanduser("~/megalith_site_research")
P3K14C_CSV = os.path.join(BASE, "data/p3k14c/p3k14c_data.csv")

POLE_LAT, POLE_LNG = 59.682122, -138.646087
EARTH_R = 6371.0
QUARTER_CIRC = EARTH_R * np.pi / 2
THRESHOLD_KM = 50
GRID_SIZE = 5  # degrees
N_MC = 1000

np.random.seed(42)


# ============================================================
# HELPERS
# ============================================================
def dist_from_gc_vec(lats, lons):
    lat1r = np.radians(POLE_LAT)
    lon1r = np.radians(POLE_LNG)
    lat2r = np.radians(lats)
    lon2r = np.radians(lons)
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat/2)**2 + np.cos(lat1r)*np.cos(lat2r)*np.sin(dlon/2)**2
    d = EARTH_R * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QUARTER_CIRC)


def mc_zscore(lats, lons, n_trials=N_MC):
    n = len(lats)
    if n < 5:
        return {"z_score": 0, "observed": 0, "random_mean": 0, "random_std": 0,
                "enrichment": 0, "n_sites": n, "note": "too few sites"}

    dists = dist_from_gc_vec(lats, lons)
    observed = int(np.sum(dists <= THRESHOLD_KM))

    rand_counts = np.empty(n_trials)
    for t in range(n_trials):
        idx = np.random.randint(0, n, size=n)
        rlats = np.clip(lats[idx] + np.random.normal(0, 2, n), -90, 90)
        rlons = np.clip(lons[idx] + np.random.normal(0, 2, n), -180, 180)
        rdists = dist_from_gc_vec(rlats, rlons)
        rand_counts[t] = int(np.sum(rdists <= THRESHOLD_KM))

    mu = float(np.mean(rand_counts))
    sigma = float(np.std(rand_counts))
    z = (observed - mu) / sigma if sigma > 0 else 0
    enrich = observed / mu if mu > 0 else 0

    return {
        "z_score": round(z, 2),
        "observed": observed,
        "random_mean": round(mu, 2),
        "random_std": round(sigma, 2),
        "enrichment": round(enrich, 2),
        "n_sites": n
    }


def grid_cell(lat, lon):
    """Return (lat_floor, lon_floor) for a 5°×5° grid cell."""
    return (int(np.floor(lat / GRID_SIZE)) * GRID_SIZE,
            int(np.floor(lon / GRID_SIZE)) * GRID_SIZE)


def cell_label(lat, lon):
    ns = "N" if lat >= 0 else "S"
    ew = "E" if lon >= 0 else "W"
    return f"{abs(lat)}-{abs(lat)+GRID_SIZE}{ns}, {abs(lon)}-{abs(lon)+GRID_SIZE}{ew}"


def classify_continent(lat, lon):
    if 35 <= lat <= 72 and -12 <= lon <= 45:
        return "Europe"
    if 10 <= lat <= 45 and 25 <= lon <= 75:
        return "Middle East/Central Asia"
    if -40 <= lat <= 37 and -25 <= lon <= 55:
        return "Africa"
    if -60 <= lat <= 15 and -90 <= lon <= -30:
        return "South America"
    if 15 <= lat <= 72 and -170 <= lon <= -50:
        return "North America"
    if -10 <= lat <= 55 and 60 <= lon <= 145:
        return "South/East Asia"
    if -50 <= lat <= -10 and 110 <= lon <= 180:
        return "Oceania"
    return "Other"


def extract_pub_year(ref, source):
    """Extract publication year from reference/source fields."""
    m = re.findall(r'\((\d{4})\)', ref)
    if m:
        return int(m[-1])
    m = re.findall(r'(?:^|[.\s:])(\d{4})(?:[.\s,;:\-]|$)', ref)
    if m:
        return int(m[0])
    m = re.findall(r'(\d{4})', source)
    if m:
        return int(m[-1])
    return None


# ============================================================
# LOAD POST-2001 p3k14c DATA
# ============================================================
print("=" * 70)
print("LOADING POST-2001 p3k14c DATA")
print("=" * 70)

all_lats, all_lons = [], []

with open(P3K14C_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            lat = float(row["Lat"])
            lon = float(row["Long"])
            if lat == 0 and lon == 0:
                continue
        except (ValueError, TypeError):
            continue

        ref = row.get("Reference", "").strip()
        source = row.get("Source", "").strip()
        pub_year = extract_pub_year(ref, source)
        if not pub_year or pub_year < 1950 or pub_year > 2025:
            continue
        if pub_year <= 2001:
            continue

        all_lats.append(lat)
        all_lons.append(lon)

lats = np.array(all_lats)
lons = np.array(all_lons)
dists = dist_from_gc_vec(lats, lons)
near_mask = dists <= THRESHOLD_KM

print(f"Post-2001 p3k14c dates: {len(lats)}")
print(f"Within {THRESHOLD_KM}km of circle: {int(np.sum(near_mask))}")


# ============================================================
# 1. GRID CELL DENSITY
# ============================================================
print(f"\n{'=' * 70}")
print("STEP 1: DENSITY BY 5°×5° GRID CELLS")
print("=" * 70)

cell_total = {}
cell_near = {}

for i in range(len(lats)):
    cell = grid_cell(lats[i], lons[i])
    cell_total[cell] = cell_total.get(cell, 0) + 1
    if near_mask[i]:
        cell_near[cell] = cell_near.get(cell, 0) + 1

# Sort by total count
top_cells = sorted(cell_total.items(), key=lambda x: -x[1])

print(f"\nTotal grid cells with data: {len(cell_total)}")
print(f"Grid cells with near-circle dates: {len(cell_near)}")

print(f"\n{'Cell':<30} {'Total':>8} {'Near':>6} {'Ratio':>8} {'Continent'}")
print("-" * 75)
grid_results = []
for cell, total in top_cells[:30]:
    near = cell_near.get(cell, 0)
    ratio = near / total * 100 if total > 0 else 0
    cont = classify_continent(cell[0] + GRID_SIZE/2, cell[1] + GRID_SIZE/2)
    label = cell_label(cell[0], cell[1])
    print(f"  {label:<28} {total:>8} {near:>6} {ratio:>7.1f}% {cont}")
    grid_results.append({
        "cell_lat": cell[0], "cell_lon": cell[1], "label": label,
        "total": total, "near_circle": near, "pct_near": round(ratio, 2),
        "continent": cont
    })

# All cells with near-circle dates
print(f"\n--- All cells with near-circle dates ---")
print(f"{'Cell':<30} {'Total':>8} {'Near':>6} {'Ratio':>8} {'Continent'}")
print("-" * 75)
near_cells = sorted([(c, cell_total[c], cell_near[c]) for c in cell_near],
                    key=lambda x: -x[2])
for cell, total, near in near_cells:
    ratio = near / total * 100
    cont = classify_continent(cell[0] + GRID_SIZE/2, cell[1] + GRID_SIZE/2)
    label = cell_label(cell[0], cell[1])
    print(f"  {label:<28} {total:>8} {near:>6} {ratio:>7.1f}% {cont}")


# ============================================================
# 2. ENRICHMENT UNIFORMITY
# ============================================================
print(f"\n{'=' * 70}")
print("STEP 2: ENRICHMENT UNIFORMITY ACROSS GRID CELLS")
print("=" * 70)

cells_with_near = [(c, cell_total[c], cell_near[c]) for c in cell_near]
cells_with_near.sort(key=lambda x: -x[2])

total_near_dates = int(np.sum(near_mask))
top3_near = sum(x[2] for x in cells_with_near[:3])
top5_near = sum(x[2] for x in cells_with_near[:5])

print(f"\nTotal near-circle dates: {total_near_dates}")
print(f"Top 3 cells contribute: {top3_near} ({top3_near/total_near_dates*100:.1f}%)")
print(f"Top 5 cells contribute: {top5_near} ({top5_near/total_near_dates*100:.1f}%)")
print(f"Remaining cells contribute: {total_near_dates - top5_near} ({(total_near_dates-top5_near)/total_near_dates*100:.1f}%)")

concentration = "CONCENTRATED" if top3_near / total_near_dates > 0.6 else \
                "MODERATELY SPREAD" if top3_near / total_near_dates > 0.4 else \
                "WELL DISTRIBUTED"
print(f"\nEnrichment pattern: {concentration}")


# ============================================================
# 3. REMOVE TOP 5 CELLS, RECOMPUTE Z
# ============================================================
print(f"\n{'=' * 70}")
print("STEP 3: Z-SCORE AFTER REMOVING TOP 5 RESEARCH-INTENSIVE CELLS")
print("=" * 70)

top5_cells = set(x[0] for x in cells_with_near[:5])
print(f"\nRemoving cells:")
for cell, total, near in cells_with_near[:5]:
    label = cell_label(cell[0], cell[1])
    cont = classify_continent(cell[0] + GRID_SIZE/2, cell[1] + GRID_SIZE/2)
    print(f"  {label} ({cont}): {total} total, {near} near")

# Filter
keep_mask = np.array([grid_cell(lats[i], lons[i]) not in top5_cells for i in range(len(lats))])
filt_lats = lats[keep_mask]
filt_lons = lons[keep_mask]

print(f"\nSites remaining after removal: {len(filt_lats)} (removed {len(lats) - len(filt_lats)})")
print(f"Computing Z-score ({N_MC} MC trials)...")
t0 = time.time()
z_minus5 = mc_zscore(filt_lats, filt_lons)
print(f"  Z = {z_minus5['z_score']:+.2f} (obs={z_minus5['observed']}, "
      f"mean={z_minus5['random_mean']:.1f}, {z_minus5['enrichment']:.2f}x) [{time.time()-t0:.1f}s]")

# Also remove top 10
top10_cells = set(x[0] for x in top_cells[:10])
keep10 = np.array([grid_cell(lats[i], lons[i]) not in top10_cells for i in range(len(lats))])
filt10_lats = lats[keep10]
filt10_lons = lons[keep10]
print(f"\nAfter removing top 10 TOTAL cells: {len(filt10_lats)} sites remain")
print(f"Computing Z-score...")
t0 = time.time()
z_minus10 = mc_zscore(filt10_lats, filt10_lons)
print(f"  Z = {z_minus10['z_score']:+.2f} (obs={z_minus10['observed']}, "
      f"mean={z_minus10['random_mean']:.1f}, {z_minus10['enrichment']:.2f}x) [{time.time()-t0:.1f}s]")

# Full dataset Z for comparison
print(f"\nFull post-2001 dataset Z-score:")
t0 = time.time()
z_full = mc_zscore(lats, lons)
print(f"  Z = {z_full['z_score']:+.2f} (obs={z_full['observed']}, "
      f"mean={z_full['random_mean']:.1f}, {z_full['enrichment']:.2f}x) [{time.time()-t0:.1f}s]")


# ============================================================
# 4. CONTINENTAL DECOMPOSITION
# ============================================================
print(f"\n{'=' * 70}")
print("STEP 4: CONTINENTAL Z-SCORE DECOMPOSITION")
print("=" * 70)

continents = {}
for i in range(len(lats)):
    cont = classify_continent(lats[i], lons[i])
    if cont not in continents:
        continents[cont] = {"lats": [], "lons": []}
    continents[cont]["lats"].append(lats[i])
    continents[cont]["lons"].append(lons[i])

continent_results = {}
print(f"\n{'Continent':<25} {'N':>8} {'Obs':>6} {'Mean':>7} {'Z':>8} {'Enrich':>8}")
print("-" * 70)

for cont in sorted(continents.keys(), key=lambda c: -len(continents[c]["lats"])):
    clats = np.array(continents[cont]["lats"])
    clons = np.array(continents[cont]["lons"])
    t0 = time.time()
    result = mc_zscore(clats, clons, n_trials=N_MC)
    elapsed = time.time() - t0
    continent_results[cont] = result
    print(f"  {cont:<23} {result['n_sites']:>8} {result['observed']:>6} "
          f"{result['random_mean']:>7.1f} {result['z_score']:>+8.2f} "
          f"{result['enrichment']:>7.2f}x  [{elapsed:.1f}s]")


# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'=' * 70}")
print("SUMMARY")
print("=" * 70)

print(f"\n  Full post-2001 Z:            {z_full['z_score']:+.2f}")
print(f"  After removing top 5 cells:  {z_minus5['z_score']:+.2f}")
print(f"  After removing top 10 cells: {z_minus10['z_score']:+.2f}")

sig_continents = [(c, r["z_score"]) for c, r in continent_results.items() if r["z_score"] > 2]
sig_continents.sort(key=lambda x: -x[1])
if sig_continents:
    print(f"\n  Continents with Z > 2:")
    for c, z in sig_continents:
        print(f"    {c}: Z = {z:+.2f}")

neg_continents = [(c, r["z_score"]) for c, r in continent_results.items() if r["z_score"] < -2]
if neg_continents:
    print(f"\n  Continents with Z < -2:")
    for c, z in sorted(neg_continents, key=lambda x: x[1]):
        print(f"    {c}: Z = {z:+.2f}")

# Verdict
if z_minus5["z_score"] > 2:
    verdict = "ROBUST — signal survives removal of top 5 research-intensive cells"
elif z_minus5["z_score"] > 0:
    verdict = "WEAKENED — signal positive but not significant after cell removal"
else:
    verdict = "FRAGILE — signal disappears when top cells removed; driven by research intensity"

print(f"\n  Verdict: {verdict}")


# ============================================================
# SAVE
# ============================================================
output = {
    "meta": {
        "name": "Fieldwork Intensity Control — Post-2001 p3k14c",
        "description": "Tests whether the post-2001 great circle signal survives controlling for regional research intensity",
        "methodology": {
            "grid_size_deg": GRID_SIZE,
            "threshold_km": THRESHOLD_KM,
            "n_mc_trials": N_MC,
            "baseline": "Distribution-matched MC (±2° Gaussian jitter)",
            "data_split": "p3k14c publication year > 2001"
        }
    },
    "full_dataset": {
        "n_sites": z_full["n_sites"],
        "z_score_50km": z_full
    },
    "grid_cell_analysis": {
        "n_cells_total": len(cell_total),
        "n_cells_with_near_circle": len(cell_near),
        "top_30_cells": grid_results,
        "all_cells_with_near_circle": [
            {"cell_lat": c[0][0], "cell_lon": c[0][1],
             "label": cell_label(c[0][0], c[0][1]),
             "total": c[1], "near_circle": c[2],
             "pct_near": round(c[2]/c[1]*100, 2),
             "continent": classify_continent(c[0][0]+GRID_SIZE/2, c[0][1]+GRID_SIZE/2)}
            for c in near_cells
        ],
        "concentration": {
            "top3_near": top3_near,
            "top5_near": top5_near,
            "total_near": total_near_dates,
            "top3_pct": round(top3_near/total_near_dates*100, 1),
            "top5_pct": round(top5_near/total_near_dates*100, 1),
            "pattern": concentration
        }
    },
    "robustness_tests": {
        "remove_top5_near_cells": {
            "cells_removed": [
                {"cell": cell_label(c[0][0], c[0][1]),
                 "total": c[1], "near": c[2],
                 "continent": classify_continent(c[0][0]+GRID_SIZE/2, c[0][1]+GRID_SIZE/2)}
                for c in cells_with_near[:5]
            ],
            "n_remaining": len(filt_lats),
            "result": z_minus5
        },
        "remove_top10_total_cells": {
            "n_remaining": len(filt10_lats),
            "result": z_minus10
        }
    },
    "continental_decomposition": continent_results,
    "verdict": verdict
}

outpath = os.path.join(BASE, "results/fieldwork_intensity_control.json")
with open(outpath, "w") as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved to {outpath}")
