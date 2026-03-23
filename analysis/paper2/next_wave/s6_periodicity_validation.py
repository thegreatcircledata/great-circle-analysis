#!/usr/bin/env python3
"""
Study 6: 36-Degree Periodicity Validation on Data-Driven Clusters
=================================================================
Previous analysis found 36-degree periodicity in inter-site arc distances
among 6 hand-picked sites. This study validates whether the periodicity
also appears in data-driven clusters (not cherry-picked).

Phases:
  1. Identify data-driven clusters via KDE on arc positions
  2. Pairwise arc distance Rayleigh test for 36-deg periodicity
  3. Multi-period scan with BH FDR correction
  4. Monte Carlo null distribution
"""

import math
import json
import os
import warnings
from datetime import datetime
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
from scipy.special import i0  # for Rayleigh
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# ── Shared Constants & Geometry ──────────────────────────────────────────────

POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2

GIZA_LAT = 29.9792
GIZA_LON = 31.1342

OUT_DIR = "/Users/elliotallan/megalith_site_research/outputs/next_wave/periodicity_36deg"
os.makedirs(OUT_DIR, exist_ok=True)


def haversine_km(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return 2 * EARTH_R_KM * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def dist_from_gc(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QUARTER_CIRC)


def forward_azimuth(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dlam = math.radians(lon2 - lon1)
    x = math.sin(dlam) * math.cos(phi2)
    y = math.cos(phi1) * math.sin(phi2) - math.sin(phi1) * math.cos(phi2) * math.cos(dlam)
    return math.degrees(math.atan2(x, y)) % 360


def arc_position(lat, lon):
    """Angular distance along GC from Giza (eastward), in degrees 0-360."""
    az_site = forward_azimuth(POLE_LAT, POLE_LON, lat, lon)
    az_giza = forward_azimuth(POLE_LAT, POLE_LON, GIZA_LAT, GIZA_LON)
    return (az_site - az_giza) % 360


def rayleigh_test(angles_deg):
    """Rayleigh test for uniformity on circular data.
    Returns (Z, p_value, mean_resultant_length, mean_direction_deg).
    """
    angles_rad = np.radians(angles_deg)
    n = len(angles_rad)
    if n == 0:
        return 0.0, 1.0, 0.0, 0.0
    C = np.sum(np.cos(angles_rad))
    S = np.sum(np.sin(angles_rad))
    R_bar = math.sqrt(C ** 2 + S ** 2) / n
    Z = n * R_bar ** 2
    # Approximate p-value (Mardia & Jupp)
    p = math.exp(-Z) * (1 + (2 * Z - Z ** 2) / (4 * n) - (24 * Z - 132 * Z ** 2 + 76 * Z ** 3 - 9 * Z ** 4) / (288 * n ** 2))
    p = max(0.0, min(1.0, p))
    mean_dir = math.degrees(math.atan2(S, C)) % 360
    return Z, p, R_bar, mean_dir


# ── Phase 0: Load & Merge Sites ─────────────────────────────────────────────

def load_sites():
    print("Loading sites from 3 databases...")

    # Pleiades
    df_pl = pd.read_csv(
        "/Users/elliotallan/megalith_site_research/data/pleiades/pleiades-places-latest.csv",
        usecols=["reprLat", "reprLong", "title"],
        low_memory=False,
    )
    df_pl = df_pl.dropna(subset=["reprLat", "reprLong"])
    df_pl = df_pl.rename(columns={"reprLat": "lat", "reprLong": "lon", "title": "name"})
    df_pl["source"] = "pleiades"

    # p3k14c — group by SiteID first
    df_p3 = pd.read_csv(
        "/Users/elliotallan/megalith_site_research/data/p3k14c/p3k14c_data.csv",
        usecols=["Lat", "Long", "SiteName", "SiteID"],
        low_memory=False,
    )
    df_p3 = df_p3.dropna(subset=["Lat", "Long"])
    df_p3 = df_p3.groupby("SiteID", as_index=False).agg(
        {"Lat": "mean", "Long": "mean", "SiteName": "first"}
    )
    df_p3 = df_p3.rename(columns={"Lat": "lat", "Long": "lon", "SiteName": "name"})
    df_p3["source"] = "p3k14c"

    # OSM
    df_osm = pd.read_csv(
        "/Users/elliotallan/megalith_site_research/data/osm/osm_archaeological_sites.csv",
        usecols=["lat", "lon", "name"],
        low_memory=False,
    )
    df_osm = df_osm.dropna(subset=["lat", "lon"])
    df_osm["source"] = "osm"

    # Combine
    all_sites = pd.concat(
        [df_pl[["lat", "lon", "name", "source"]],
         df_p3[["lat", "lon", "name", "source"]],
         df_osm[["lat", "lon", "name", "source"]]],
        ignore_index=True,
    )
    print(f"  Total raw sites: {len(all_sites):,}")

    # Filter within 50km of GC
    all_sites["gc_dist"] = all_sites.apply(
        lambda r: dist_from_gc(r["lat"], r["lon"]), axis=1
    )
    near_gc = all_sites[all_sites["gc_dist"] <= 50.0].copy()
    print(f"  Sites within 50km of GC: {len(near_gc):,}")

    # Deduplicate within 10km
    near_gc = near_gc.sort_values("gc_dist").reset_index(drop=True)
    keep = np.ones(len(near_gc), dtype=bool)
    lats = near_gc["lat"].values
    lons = near_gc["lon"].values
    for i in range(len(near_gc)):
        if not keep[i]:
            continue
        for j in range(i + 1, len(near_gc)):
            if not keep[j]:
                continue
            if haversine_km(lats[i], lons[i], lats[j], lons[j]) < 10.0:
                keep[j] = False
    deduped = near_gc[keep].copy().reset_index(drop=True)
    print(f"  After 10km dedup: {len(deduped):,}")

    # Compute arc positions
    deduped["arc_pos"] = deduped.apply(
        lambda r: arc_position(r["lat"], r["lon"]), axis=1
    )
    return deduped


# ── Phase 1: KDE Cluster Detection ──────────────────────────────────────────

def find_clusters(sites):
    print("\nPhase 1: KDE cluster detection...")
    arc_vals = sites["arc_pos"].values

    # KDE with bandwidth ~2 degrees
    # We work on a circular domain [0, 360); mirror data for boundary handling
    extended = np.concatenate([arc_vals - 360, arc_vals, arc_vals + 360])
    kde = gaussian_kde(extended, bw_method=2.0 / np.std(extended))

    x_grid = np.linspace(0, 360, 3601)
    density = kde(x_grid)

    # Normalize so it integrates to ~1 over [0, 360]
    density = density / (np.sum(density) * (x_grid[1] - x_grid[0]))

    mean_density = np.mean(density)
    threshold = 2.0 * mean_density

    # Find local maxima
    maxima_idx = argrelextrema(density, np.greater, order=10)[0]
    peak_positions = []
    peak_densities = []
    for idx in maxima_idx:
        if density[idx] > threshold:
            peak_positions.append(x_grid[idx])
            peak_densities.append(density[idx])

    # Merge peaks within 3 degrees
    merged_peaks = []
    merged_densities = []
    used = set()
    for i, pos in enumerate(peak_positions):
        if i in used:
            continue
        group_pos = [pos]
        group_den = [peak_densities[i]]
        for j in range(i + 1, len(peak_positions)):
            if j in used:
                continue
            arc_diff = min(abs(peak_positions[j] - pos), 360 - abs(peak_positions[j] - pos))
            if arc_diff < 3.0:
                group_pos.append(peak_positions[j])
                group_den.append(peak_densities[j])
                used.add(j)
        best = np.argmax(group_den)
        merged_peaks.append(group_pos[best])
        merged_densities.append(group_den[best])

    print(f"  Found {len(merged_peaks)} cluster centers (density > 2x mean)")
    for i, (pos, den) in enumerate(zip(merged_peaks, merged_densities)):
        print(f"    Cluster {i+1}: arc_pos = {pos:.1f} deg, density = {den:.4f}")

    # Plot KDE
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(x_grid, density, "b-", lw=1.2, label="KDE density")
    ax.axhline(threshold, color="red", ls="--", lw=0.8, label=f"2x mean ({threshold:.4f})")
    for pos in merged_peaks:
        ax.axvline(pos, color="green", ls=":", alpha=0.7, lw=0.8)
    ax.scatter(merged_peaks, merged_densities, c="green", zorder=5, s=40, label="Cluster centers")
    ax.set_xlabel("Arc position from Giza (degrees, eastward)")
    ax.set_ylabel("Density")
    ax.set_title(f"Study 6 Phase 1: KDE of Site Arc Positions ({len(merged_peaks)} clusters)")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 360)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "kde_arc_positions.png"), dpi=150)
    plt.close(fig)

    return np.array(merged_peaks), x_grid, density, threshold


# ── Phase 2: 36-Degree Rayleigh Test ────────────────────────────────────────

def test_period(cluster_centers, period):
    """Rayleigh test for a given period on pairwise distances."""
    pairs = list(combinations(cluster_centers, 2))
    dists = []
    for a, b in pairs:
        d = min(abs(a - b), 360 - abs(a - b))
        dists.append(d)
    dists = np.array(dists)

    # Map residuals to circular variable
    residuals = dists % period
    circular_angles = residuals * (360.0 / period)
    Z, p, R_bar, mean_dir = rayleigh_test(circular_angles)
    return Z, p, R_bar, mean_dir, dists


def phase2(cluster_centers):
    print("\nPhase 2: 36-degree Rayleigh test on pairwise arc distances...")
    n_clusters = len(cluster_centers)
    n_pairs = n_clusters * (n_clusters - 1) // 2
    print(f"  {n_clusters} clusters -> {n_pairs} pairwise distances")

    Z, p, R_bar, mean_dir, dists = test_period(cluster_centers, 36.0)
    print(f"  Rayleigh Z = {Z:.4f}, p = {p:.6f}, R_bar = {R_bar:.4f}, mean_dir = {mean_dir:.1f} deg")
    return Z, p, R_bar, mean_dir, dists


# ── Phase 3: Multi-Period Scan ───────────────────────────────────────────────

def phase3(cluster_centers):
    print("\nPhase 3: Multi-period scan with BH FDR correction...")
    periods = [10, 12, 15, 18, 20, 22.5, 24, 30, 36, 40, 45, 60, 72, 90, 120]
    results = []
    pvals = []
    for period in periods:
        Z, p, R_bar, mean_dir, _ = test_period(cluster_centers, period)
        results.append({
            "period": period,
            "rayleigh_Z": round(Z, 4),
            "p_value": p,
            "R_bar": round(R_bar, 4),
            "mean_dir": round(mean_dir, 1),
        })
        pvals.append(p)
        sig = "*" if p < 0.05 else ""
        print(f"    Period {period:6.1f} deg:  Z={Z:7.3f}  p={p:.6f}  R={R_bar:.4f}  {sig}")

    # BH FDR correction
    reject, pvals_corrected, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
    for i, r in enumerate(results):
        r["p_corrected"] = round(float(pvals_corrected[i]), 6)
        r["survives_bh"] = bool(reject[i])
        r["p_value"] = round(r["p_value"], 6)

    surviving = [r for r in results if r["survives_bh"]]
    print(f"\n  Periods surviving BH correction (alpha=0.05): {len(surviving)}")
    for r in surviving:
        print(f"    {r['period']} deg: p_raw={r['p_value']:.6f}, p_corrected={r['p_corrected']:.6f}")

    # Periodogram plot
    fig, ax = plt.subplots(figsize=(10, 5))
    period_vals = [r["period"] for r in results]
    neg_log_p = [-math.log10(max(r["p_value"], 1e-20)) for r in results]
    neg_log_p_corr = [-math.log10(max(r["p_corrected"], 1e-20)) for r in results]

    ax.bar([x - 0.8 for x in range(len(periods))], neg_log_p, width=0.8,
           label="Raw p-value", alpha=0.7, color="steelblue")
    ax.bar([x + 0.0 for x in range(len(periods))], neg_log_p_corr, width=0.8,
           label="BH-corrected p-value", alpha=0.7, color="darkorange")
    ax.axhline(-math.log10(0.05), color="red", ls="--", lw=1, label="p = 0.05")
    ax.axhline(-math.log10(0.01), color="darkred", ls="--", lw=1, label="p = 0.01")
    ax.set_xticks(range(len(periods)))
    ax.set_xticklabels([str(p) for p in period_vals], rotation=45, ha="right")
    ax.set_xlabel("Period (degrees)")
    ax.set_ylabel("-log10(p)")
    ax.set_title("Study 6 Phase 3: Multi-Period Scan (Rayleigh Test)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "periodogram.png"), dpi=150)
    plt.close(fig)

    return results


# ── Phase 4: Monte Carlo Null ────────────────────────────────────────────────

def phase4(cluster_centers, observed_Z):
    print("\nPhase 4: Monte Carlo null distribution (1000 iterations)...")
    n_clusters = len(cluster_centers)
    np.random.seed(42)
    null_Z = []
    for _ in range(1000):
        fake_centers = np.random.uniform(0, 360, size=n_clusters)
        Z, _, _, _, _ = test_period(fake_centers, 36.0)
        null_Z.append(Z)
    null_Z = np.array(null_Z)
    percentile = np.mean(null_Z <= observed_Z) * 100

    print(f"  Observed Z = {observed_Z:.4f}")
    print(f"  Null Z: mean = {null_Z.mean():.4f}, std = {null_Z.std():.4f}, "
          f"max = {null_Z.max():.4f}")
    print(f"  Percentile rank: {percentile:.1f}%")

    mc_pvalue = np.mean(null_Z >= observed_Z)
    print(f"  Monte Carlo p-value: {mc_pvalue:.4f}")

    return percentile, mc_pvalue, null_Z


# ── Verdict ──────────────────────────────────────────────────────────────────

def compute_verdict(p36_raw, p36_corrected, survives_bh_36, period_results, mc_pvalue):
    sig_periods = [r for r in period_results if r["survives_bh"]]

    if p36_raw > 0.05:
        return "FALSIFIED"
    if survives_bh_36 and mc_pvalue < 0.05:
        if len(sig_periods) >= 4:
            return "WEAKENED"
        return "CONFIRMED"
    if p36_raw < 0.01 and not survives_bh_36:
        return "WEAKENED"
    if p36_raw < 0.05:
        return "WEAKENED"
    return "FALSIFIED"


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("Study 6: 36-Degree Periodicity Validation on Data-Driven Clusters")
    print("=" * 70)

    sites = load_sites()

    # Phase 1
    cluster_centers, x_grid, density, threshold = find_clusters(sites)
    n_clusters = len(cluster_centers)

    if n_clusters < 3:
        print("\nERROR: Too few clusters for meaningful pairwise analysis.")
        return

    # Phase 2
    Z36, p36, R36, md36, pairwise_dists = phase2(cluster_centers)

    # Phase 3
    period_results = phase3(cluster_centers)

    # Find 36-deg entry
    entry_36 = next(r for r in period_results if r["period"] == 36)
    survives_36 = entry_36["survives_bh"]

    # Phase 4
    percentile, mc_pvalue, null_Z = phase4(cluster_centers, Z36)

    # Verdict
    verdict = compute_verdict(p36, entry_36["p_corrected"], survives_36, period_results, mc_pvalue)

    print(f"\n{'='*70}")
    print(f"VERDICT: {verdict}")
    print(f"{'='*70}")

    # ── Save results.json ────────────────────────────────────────────────
    results = {
        "study": "S6: 36-Degree Periodicity Validation on Data-Driven Clusters",
        "date": datetime.now().strftime("%Y-%m-%d"),
        "verdict": verdict,
        "sites_total": int(len(sites)),
        "n_clusters": int(n_clusters),
        "cluster_centers_deg": [round(c, 2) for c in sorted(cluster_centers)],
        "n_pairwise_distances": int(len(pairwise_dists)),
        "phase2_36deg": {
            "rayleigh_Z": round(Z36, 4),
            "p_value": round(p36, 6),
            "R_bar": round(R36, 4),
            "mean_direction_deg": round(md36, 1),
        },
        "phase3_multi_period": period_results,
        "phase3_surviving_bh": [r["period"] for r in period_results if r["survives_bh"]],
        "phase4_monte_carlo": {
            "n_iterations": 1000,
            "observed_Z": round(Z36, 4),
            "null_mean_Z": round(float(null_Z.mean()), 4),
            "null_std_Z": round(float(null_Z.std()), 4),
            "percentile_rank": round(percentile, 1),
            "mc_p_value": round(float(mc_pvalue), 4),
        },
    }

    with open(os.path.join(OUT_DIR, "results.json"), "w") as f:
        json.dump(results, f, indent=2)

    # ── Save RESULTS.md ──────────────────────────────────────────────────
    surviving = results["phase3_surviving_bh"]
    surv_str = ", ".join(str(s) for s in surviving) if surviving else "None"

    md = f"""# Study 6: 36-Degree Periodicity Validation on Data-Driven Clusters

**Date:** {results['date']}
**Verdict:** {verdict}

## Overview

Previous analysis found 36-degree periodicity in inter-site arc distances
among 6 hand-picked sites along the proposed great circle (GC). This study
tests whether that periodicity persists in data-driven clusters identified
from {results['sites_total']:,} archaeological sites within 50 km of the GC.

## Phase 1: Cluster Detection

- Sites within 50 km of GC: {results['sites_total']:,}
- KDE bandwidth: ~2 degrees
- Clusters identified (density > 2x mean): **{n_clusters}**
- Cluster arc positions (degrees from Giza, eastward):
  {', '.join(f'{c:.1f}' for c in sorted(cluster_centers))}

## Phase 2: 36-Degree Rayleigh Test

Pairwise arc distances between {n_clusters} cluster centers ({results['n_pairwise_distances']} pairs)
tested for 36-degree periodicity.

| Metric | Value |
|--------|-------|
| Rayleigh Z | {results['phase2_36deg']['rayleigh_Z']:.4f} |
| p-value | {results['phase2_36deg']['p_value']:.6f} |
| Mean resultant length (R) | {results['phase2_36deg']['R_bar']:.4f} |
| Mean direction | {results['phase2_36deg']['mean_direction_deg']:.1f} deg |

## Phase 3: Multi-Period Scan

Tested periods: 10, 12, 15, 18, 20, 22.5, 24, 30, 36, 40, 45, 60, 72, 90, 120 degrees.
BH FDR correction applied.

| Period | Rayleigh Z | p (raw) | p (BH) | Survives? |
|--------|-----------|---------|--------|-----------|
"""
    for r in period_results:
        flag = "Yes" if r["survives_bh"] else ""
        md += f"| {r['period']} | {r['rayleigh_Z']:.4f} | {r['p_value']:.6f} | {r['p_corrected']:.6f} | {flag} |\n"

    md += f"""
**Periods surviving BH correction:** {surv_str}

## Phase 4: Monte Carlo Null Distribution

- 1000 random cluster sets (uniform on 0-360 deg arc, same N={n_clusters})
- Observed 36-deg Rayleigh Z: {results['phase4_monte_carlo']['observed_Z']:.4f}
- Null distribution: mean = {results['phase4_monte_carlo']['null_mean_Z']:.4f}, std = {results['phase4_monte_carlo']['null_std_Z']:.4f}
- Percentile rank of observed Z: {results['phase4_monte_carlo']['percentile_rank']:.1f}%
- Monte Carlo p-value: {results['phase4_monte_carlo']['mc_p_value']:.4f}

## Verdict: {verdict}

"""
    if verdict == "CONFIRMED":
        md += ("The 36-degree periodicity is statistically significant on data-driven clusters, "
               "survives multiple-comparison correction, and exceeds Monte Carlo null expectations.\n")
    elif verdict == "WEAKENED":
        md += ("The 36-degree signal shows some statistical significance but is weakened by "
               "one or more factors: multiple other periods are also significant, the signal "
               "does not survive BH correction, or Monte Carlo comparison is inconclusive.\n")
    else:
        md += ("The 36-degree periodicity is NOT statistically significant on data-driven clusters. "
               "The original finding appears to depend on hand-picked site selection.\n")

    md += """
## Outputs

- `results.json` — Full numerical results
- `kde_arc_positions.png` — KDE density plot of site arc positions with cluster centers
- `periodogram.png` — Multi-period scan bar chart
"""

    with open(os.path.join(OUT_DIR, "RESULTS.md"), "w") as f:
        f.write(md)

    print(f"\nOutputs written to {OUT_DIR}/")
    print("  results.json")
    print("  RESULTS.md")
    print("  kde_arc_positions.png")
    print("  periodogram.png")


if __name__ == "__main__":
    main()
