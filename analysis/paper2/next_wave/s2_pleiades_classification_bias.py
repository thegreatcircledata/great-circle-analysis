#!/usr/bin/env python3
"""
Study 2: Pleiades Classification Bias Test
===========================================
Tests whether Pleiades (built primarily by classicists studying the Mediterranean)
has systematically different monument-to-settlement classification ratios near the
great circle vs. away from it — a potential cataloguing artifact.

Phases:
  1. Regional Classification Profile — Spearman correlation of monument_ratio
     with mean distance-to-GC per region.
  2. Within-Region Divergence — on-corridor vs off-corridor monument_ratio
     within each region.
  3. Unclassified Site Analysis — compare fraction unclassified on/off corridor.
"""

import csv
import json
import math
import os
import sys
import textwrap

import numpy as np
from scipy import stats

# ── buffered stdout ──────────────────────────────────────────────────────────
sys.stdout = os.fdopen(sys.stdout.fileno(), "w", buffering=1)

np.random.seed(42)

# ── paths ────────────────────────────────────────────────────────────────────
DATA_CSV = "/Users/elliotallan/megalith_site_research/data/pleiades/pleiades-places-latest.csv"
OUT_DIR  = "/Users/elliotallan/megalith_site_research/outputs/next_wave/classification_bias"
os.makedirs(OUT_DIR, exist_ok=True)

# ── shared constants & geometry ──────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
THRESHOLD_KM = 50


def haversine_km(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return 2 * EARTH_R_KM * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def dist_from_gc(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QUARTER_CIRC)


# ── type sets ────────────────────────────────────────────────────────────────
MONUMENTAL_TYPES = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus", "shrine",
}
SETTLEMENT_TYPES = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement", "townhouse", "production",
}

# ── region definitions ───────────────────────────────────────────────────────
REGIONS = {
    "Egypt/Nubia":          {"lat": (22, 31), "lon": (24, 36)},
    "Levant":               {"lat": (29, 37), "lon": (34, 37)},
    "Anatolia":             {"lat": (36, 42), "lon": (26, 45)},
    "Greece/Balkans":       {"lat": (35, 42), "lon": (19, 30)},
    "Italy":                {"lat": (36, 47), "lon": (7, 19)},
    "North Africa (non-E)": {"lat": (25, 37), "lon": (-10, 24)},
    "Iberia":               {"lat": (36, 44), "lon": (-10, 4)},
    "Mesopotamia":          {"lat": (30, 37), "lon": (38, 48)},
    "Iran/Central Asia":    {"lat": (25, 42), "lon": (48, 75)},
    "South Asia":           {"lat": (5, 35),  "lon": (65, 90)},
}


def in_region(lat, lon, bounds):
    return (bounds["lat"][0] <= lat <= bounds["lat"][1] and
            bounds["lon"][0] <= lon <= bounds["lon"][1])


# ── classify a featureTypes string ───────────────────────────────────────────
def classify_site(feature_str):
    """Return ('monument', 'settlement', 'unclassified', or 'both')."""
    if not feature_str or not feature_str.strip():
        return "unclassified"
    tokens = {t.strip().lower() for t in feature_str.split(",")}
    is_mon = bool(tokens & MONUMENTAL_TYPES)
    is_set = bool(tokens & SETTLEMENT_TYPES)
    if is_mon and is_set:
        return "both"
    if is_mon:
        return "monument"
    if is_set:
        return "settlement"
    return "unclassified"


# ══════════════════════════════════════════════════════════════════════════════
#  LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("Study 2: Pleiades Classification Bias Test")
print("=" * 70)
print()

print("Loading Pleiades data …")
sites = []
skipped = 0
with open(DATA_CSV, newline="", encoding="utf-8") as fh:
    reader = csv.DictReader(fh)
    for row in reader:
        try:
            lat = float(row["reprLat"])
            lon = float(row["reprLong"])
        except (ValueError, TypeError):
            skipped += 1
            continue
        feat = row.get("featureTypes", "")
        cat = classify_site(feat)
        dgc = dist_from_gc(lat, lon)
        sites.append({
            "lat": lat,
            "lon": lon,
            "featureTypes": feat,
            "category": cat,
            "dist_gc_km": dgc,
        })

print(f"  Loaded {len(sites):,} sites with coordinates  (skipped {skipped} without coords)")
cats = {}
for s in sites:
    cats[s["category"]] = cats.get(s["category"], 0) + 1
for c in sorted(cats):
    print(f"    {c:15s}: {cats[c]:>6,}")
print()

# ══════════════════════════════════════════════════════════════════════════════
#  PHASE 1 — Regional Classification Profile
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("PHASE 1: Regional Classification Profile")
print("=" * 70)
print()

region_stats = {}
for rname, bounds in REGIONS.items():
    n_mon = 0
    n_set = 0
    n_both = 0
    n_unc = 0
    dists = []
    for s in sites:
        if not in_region(s["lat"], s["lon"], bounds):
            continue
        dists.append(s["dist_gc_km"])
        if s["category"] == "monument":
            n_mon += 1
        elif s["category"] == "settlement":
            n_set += 1
        elif s["category"] == "both":
            n_both += 1
        else:
            n_unc += 1
    total = n_mon + n_set + n_both + n_unc
    denom = n_mon + n_set
    mon_ratio = n_mon / denom if denom > 0 else float("nan")
    mean_dist = float(np.mean(dists)) if dists else float("nan")
    region_stats[rname] = {
        "n_monument": n_mon,
        "n_settlement": n_set,
        "n_both": n_both,
        "n_unclassified": n_unc,
        "total": total,
        "monument_ratio": mon_ratio,
        "mean_dist_gc_km": mean_dist,
    }

# Print table
print(f"{'Region':<22s} {'Mon':>5s} {'Set':>5s} {'Both':>5s} {'Unc':>5s} {'Total':>6s} {'MonRat':>7s} {'MeanDist':>9s}")
print("-" * 70)
for rname in REGIONS:
    rs = region_stats[rname]
    print(f"{rname:<22s} {rs['n_monument']:>5d} {rs['n_settlement']:>5d} "
          f"{rs['n_both']:>5d} {rs['n_unclassified']:>5d} {rs['total']:>6d} "
          f"{rs['monument_ratio']:>7.3f} {rs['mean_dist_gc_km']:>9.1f}")
print()

# Spearman correlation: monument_ratio vs mean_dist_gc
valid = [(r, region_stats[r]) for r in REGIONS
         if not math.isnan(region_stats[r]["monument_ratio"]) and region_stats[r]["total"] >= 10]
ratios = [v[1]["monument_ratio"] for v in valid]
dists_r = [v[1]["mean_dist_gc_km"] for v in valid]

if len(valid) >= 4:
    spearman_r, spearman_p = stats.spearmanr(ratios, dists_r)
else:
    spearman_r, spearman_p = float("nan"), float("nan")

print(f"Spearman correlation (monument_ratio vs mean_dist_to_GC):")
print(f"  rho = {spearman_r:+.4f},  p = {spearman_p:.4f}  (n = {len(valid)} regions)")
if spearman_p > 0.1:
    phase1_verdict = "NO significant correlation — classification ratio not driven by GC proximity"
elif spearman_p <= 0.05 and spearman_r < 0:
    phase1_verdict = "SIGNIFICANT negative correlation — regions closer to GC have higher monument ratio"
elif spearman_p <= 0.05 and spearman_r > 0:
    phase1_verdict = "SIGNIFICANT positive correlation — regions farther from GC have higher monument ratio"
else:
    phase1_verdict = "MARGINAL correlation (0.05 < p < 0.1)"
print(f"  → {phase1_verdict}")
print()

# ══════════════════════════════════════════════════════════════════════════════
#  PHASE 2 — Within-Region Divergence
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("PHASE 2: Within-Region Divergence (on-corridor vs off-corridor)")
print("=" * 70)
print()

ON_CORRIDOR_KM = 50
OFF_CORRIDOR_KM = 200

within_region = {}
for rname, bounds in REGIONS.items():
    on_mon = on_set = off_mon = off_set = 0
    for s in sites:
        if not in_region(s["lat"], s["lon"], bounds):
            continue
        if s["category"] not in ("monument", "settlement", "both"):
            continue
        is_mon = s["category"] in ("monument", "both")
        is_set = s["category"] in ("settlement", "both")
        if s["dist_gc_km"] < ON_CORRIDOR_KM:
            if is_mon:
                on_mon += 1
            if is_set:
                on_set += 1
        elif s["dist_gc_km"] > OFF_CORRIDOR_KM:
            if is_mon:
                off_mon += 1
            if is_set:
                off_set += 1

    on_total = on_mon + on_set
    off_total = off_mon + off_set
    on_ratio = on_mon / on_total if on_total > 0 else float("nan")
    off_ratio = off_mon / off_total if off_total > 0 else float("nan")

    # Fisher exact test if both have data
    if on_total > 0 and off_total > 0:
        table = [[on_mon, on_set], [off_mon, off_set]]
        odds_ratio, fisher_p = stats.fisher_exact(table)
    else:
        odds_ratio, fisher_p = float("nan"), float("nan")

    within_region[rname] = {
        "on_mon": on_mon, "on_set": on_set, "on_total": on_total, "on_ratio": on_ratio,
        "off_mon": off_mon, "off_set": off_set, "off_total": off_total, "off_ratio": off_ratio,
        "fisher_odds": odds_ratio, "fisher_p": fisher_p,
    }

print(f"{'Region':<22s} {'On_M':>5s} {'On_S':>5s} {'OnRat':>7s} "
      f"{'Off_M':>5s} {'Off_S':>5s} {'OffRat':>7s} {'Fisher_p':>9s}")
print("-" * 78)
for rname in REGIONS:
    wr = within_region[rname]
    on_r = f"{wr['on_ratio']:.3f}" if not math.isnan(wr["on_ratio"]) else "   n/a"
    off_r = f"{wr['off_ratio']:.3f}" if not math.isnan(wr["off_ratio"]) else "   n/a"
    fp = f"{wr['fisher_p']:.4f}" if not math.isnan(wr["fisher_p"]) else "      n/a"
    print(f"{rname:<22s} {wr['on_mon']:>5d} {wr['on_set']:>5d} {on_r:>7s} "
          f"{wr['off_mon']:>5d} {wr['off_set']:>5d} {off_r:>7s} {fp:>9s}")
print()

# Summarise within-region findings
sig_regions = [r for r in REGIONS if within_region[r]["fisher_p"] < 0.05
               and within_region[r]["on_total"] >= 5 and within_region[r]["off_total"] >= 5]
testable_regions = [r for r in REGIONS
                    if within_region[r]["on_total"] >= 5 and within_region[r]["off_total"] >= 5]

print(f"Testable regions (on >= 5 AND off >= 5): {len(testable_regions)}")
print(f"Regions with Fisher p < 0.05: {len(sig_regions)}")
for r in sig_regions:
    wr = within_region[r]
    print(f"  {r}: on_ratio={wr['on_ratio']:.3f} vs off_ratio={wr['off_ratio']:.3f} (p={wr['fisher_p']:.4f})")
print()

# ══════════════════════════════════════════════════════════════════════════════
#  PHASE 3 — Unclassified Site Analysis
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("PHASE 3: Unclassified Site Analysis")
print("=" * 70)
print()

on_unc = on_total_all = off_unc = off_total_all = 0
for s in sites:
    if s["dist_gc_km"] < ON_CORRIDOR_KM:
        on_total_all += 1
        if s["category"] == "unclassified":
            on_unc += 1
    elif s["dist_gc_km"] > OFF_CORRIDOR_KM:
        off_total_all += 1
        if s["category"] == "unclassified":
            off_unc += 1

on_unc_frac = on_unc / on_total_all if on_total_all > 0 else float("nan")
off_unc_frac = off_unc / off_total_all if off_total_all > 0 else float("nan")

print(f"On-corridor  (<{ON_CORRIDOR_KM} km):  {on_unc:>5d} unclassified / {on_total_all:>6d} total = {on_unc_frac:.4f}")
print(f"Off-corridor (>{OFF_CORRIDOR_KM} km): {off_unc:>5d} unclassified / {off_total_all:>6d} total = {off_unc_frac:.4f}")

if on_total_all > 0 and off_total_all > 0:
    table3 = [[on_unc, on_total_all - on_unc], [off_unc, off_total_all - off_unc]]
    odds3, p3 = stats.fisher_exact(table3)
    print(f"Fisher exact: OR = {odds3:.3f}, p = {p3:.4f}")
else:
    odds3, p3 = float("nan"), float("nan")
    print("Insufficient data for Fisher test")
print()

# ══════════════════════════════════════════════════════════════════════════════
#  OVERALL VERDICT
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("VERDICT")
print("=" * 70)
print()

# Decision logic
fatal = (spearman_p <= 0.05 and spearman_r < 0 and len(sig_regions) == 0)
concern = (spearman_p <= 0.1 and len(testable_regions) < 3)
if fatal:
    verdict = "FATAL"
    verdict_text = ("Monument ratio strongly correlates with GC proximity AND "
                    "within-region divergence disappears. Classification bias is "
                    "a plausible explanation for apparent monument clustering.")
elif concern:
    verdict_text = ("Regional classification varies AND within-region sample sizes "
                    "too small to confirm or deny bias.")
    verdict = "CONCERN"
else:
    verdict = "CLEAR"
    verdict_text = ("Monument ratio is NOT significantly correlated with distance-to-GC "
                    "at the regional level (Spearman p > 0.1), OR within-region divergence "
                    "persists, indicating the pattern is not a classification artifact.")

print(f"  {verdict}: {verdict_text}")
print()

# ══════════════════════════════════════════════════════════════════════════════
#  PLOTS
# ══════════════════════════════════════════════════════════════════════════════
print("Generating plots …")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

# ── Plot 1: monument_ratio vs mean_dist_gc scatter ──
ax = axes[0]
for rname in REGIONS:
    rs = region_stats[rname]
    if math.isnan(rs["monument_ratio"]):
        continue
    ax.scatter(rs["mean_dist_gc_km"], rs["monument_ratio"],
               s=max(20, rs["total"] / 5), zorder=3)
    ax.annotate(rname, (rs["mean_dist_gc_km"], rs["monument_ratio"]),
                fontsize=7, ha="left", va="bottom")
ax.set_xlabel("Mean Distance to Great Circle (km)")
ax.set_ylabel("Monument Ratio")
ax.set_title(f"Phase 1: Region Monument Ratio vs GC Distance\n"
             f"Spearman ρ={spearman_r:+.3f}, p={spearman_p:.3f}")
ax.grid(alpha=0.3)

# ── Plot 2: within-region on vs off corridor ──
ax = axes[1]
rnames_plot = [r for r in REGIONS
               if not math.isnan(within_region[r]["on_ratio"])
               and not math.isnan(within_region[r]["off_ratio"])]
x_pos = np.arange(len(rnames_plot))
bar_w = 0.35
on_vals = [within_region[r]["on_ratio"] for r in rnames_plot]
off_vals = [within_region[r]["off_ratio"] for r in rnames_plot]
ax.bar(x_pos - bar_w / 2, on_vals, bar_w, label=f"On-corridor (<{ON_CORRIDOR_KM} km)", color="steelblue")
ax.bar(x_pos + bar_w / 2, off_vals, bar_w, label=f"Off-corridor (>{OFF_CORRIDOR_KM} km)", color="coral")
ax.set_xticks(x_pos)
ax.set_xticklabels([r[:12] for r in rnames_plot], rotation=45, ha="right", fontsize=7)
ax.set_ylabel("Monument Ratio")
ax.set_title("Phase 2: On- vs Off-Corridor Monument Ratio")
ax.legend(fontsize=7)
ax.grid(axis="y", alpha=0.3)

# ── Plot 3: unclassified fraction ──
ax = axes[2]
labels = [f"On (<{ON_CORRIDOR_KM} km)", f"Off (>{OFF_CORRIDOR_KM} km)"]
fracs = [on_unc_frac, off_unc_frac]
totals = [on_total_all, off_total_all]
bars = ax.bar(labels, fracs, color=["steelblue", "coral"], width=0.5)
for bar, frac, tot in zip(bars, fracs, totals):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
            f"{frac:.3f}\n(n={tot:,})", ha="center", va="bottom", fontsize=9)
ax.set_ylabel("Fraction Unclassified")
ax.set_title(f"Phase 3: Unclassified Fraction\nFisher p={p3:.4f}" if not math.isnan(p3) else "Phase 3: Unclassified Fraction")
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plot_path = os.path.join(OUT_DIR, "classification_bias_results.png")
plt.savefig(plot_path, dpi=150)
plt.close()
print(f"  Saved → {plot_path}")

# ══════════════════════════════════════════════════════════════════════════════
#  SAVE results.json
# ══════════════════════════════════════════════════════════════════════════════
results = {
    "study": "S2: Pleiades Classification Bias Test",
    "verdict": verdict,
    "verdict_text": verdict_text,
    "phase1": {
        "description": "Spearman correlation of regional monument_ratio with mean distance-to-GC",
        "spearman_rho": round(spearman_r, 4) if not math.isnan(spearman_r) else None,
        "spearman_p": round(spearman_p, 4) if not math.isnan(spearman_p) else None,
        "n_regions": len(valid),
        "region_stats": {r: {k: (round(v, 4) if isinstance(v, float) and not math.isnan(v) else v)
                             for k, v in region_stats[r].items()}
                         for r in REGIONS},
    },
    "phase2": {
        "description": "Within-region on-corridor vs off-corridor monument ratio (Fisher exact test)",
        "on_corridor_threshold_km": ON_CORRIDOR_KM,
        "off_corridor_threshold_km": OFF_CORRIDOR_KM,
        "testable_regions": len(testable_regions),
        "significant_regions": len(sig_regions),
        "significant_region_names": sig_regions,
        "per_region": {r: {k: (round(v, 4) if isinstance(v, float) and not math.isnan(v) else
                              (None if isinstance(v, float) and math.isnan(v) else v))
                           for k, v in within_region[r].items()}
                       for r in REGIONS},
    },
    "phase3": {
        "description": "Fraction of sites unclassified: on-corridor vs off-corridor",
        "on_unclassified": on_unc,
        "on_total": on_total_all,
        "on_fraction": round(on_unc_frac, 4) if not math.isnan(on_unc_frac) else None,
        "off_unclassified": off_unc,
        "off_total": off_total_all,
        "off_fraction": round(off_unc_frac, 4) if not math.isnan(off_unc_frac) else None,
        "fisher_odds": round(odds3, 4) if not math.isnan(odds3) else None,
        "fisher_p": round(p3, 4) if not math.isnan(p3) else None,
    },
}

json_path = os.path.join(OUT_DIR, "results.json")
with open(json_path, "w") as fh:
    json.dump(results, fh, indent=2)
print(f"  Saved → {json_path}")

# ══════════════════════════════════════════════════════════════════════════════
#  SAVE RESULTS.md
# ══════════════════════════════════════════════════════════════════════════════
md_lines = []
md_lines.append("# Study 2: Pleiades Classification Bias Test")
md_lines.append("")
md_lines.append(f"**Verdict: {verdict}**")
md_lines.append("")
md_lines.append(verdict_text)
md_lines.append("")

md_lines.append("## Phase 1: Regional Classification Profile")
md_lines.append("")
md_lines.append(f"Spearman rank correlation of monument_ratio with mean distance-to-GC:")
md_lines.append(f"- rho = {spearman_r:+.4f}, p = {spearman_p:.4f} (n = {len(valid)} regions)")
md_lines.append(f"- {phase1_verdict}")
md_lines.append("")
md_lines.append("| Region | Mon | Set | Both | Unc | Total | MonRat | MeanDist (km) |")
md_lines.append("|--------|----:|----:|-----:|----:|------:|-------:|--------------:|")
for rname in REGIONS:
    rs = region_stats[rname]
    md_lines.append(f"| {rname} | {rs['n_monument']} | {rs['n_settlement']} | "
                    f"{rs['n_both']} | {rs['n_unclassified']} | {rs['total']} | "
                    f"{rs['monument_ratio']:.3f} | {rs['mean_dist_gc_km']:.1f} |")
md_lines.append("")

md_lines.append("## Phase 2: Within-Region Divergence")
md_lines.append("")
md_lines.append(f"On-corridor: <{ON_CORRIDOR_KM} km from GC | Off-corridor: >{OFF_CORRIDOR_KM} km from GC")
md_lines.append("")
md_lines.append("| Region | On_Mon | On_Set | On_Ratio | Off_Mon | Off_Set | Off_Ratio | Fisher p |")
md_lines.append("|--------|-------:|-------:|---------:|--------:|--------:|----------:|---------:|")
for rname in REGIONS:
    wr = within_region[rname]
    on_r = f"{wr['on_ratio']:.3f}" if not math.isnan(wr["on_ratio"]) else "n/a"
    off_r = f"{wr['off_ratio']:.3f}" if not math.isnan(wr["off_ratio"]) else "n/a"
    fp = f"{wr['fisher_p']:.4f}" if not math.isnan(wr["fisher_p"]) else "n/a"
    md_lines.append(f"| {rname} | {wr['on_mon']} | {wr['on_set']} | {on_r} | "
                    f"{wr['off_mon']} | {wr['off_set']} | {off_r} | {fp} |")
md_lines.append("")
md_lines.append(f"Testable regions (both groups >= 5): {len(testable_regions)}")
md_lines.append(f"Significant divergence (p < 0.05): {len(sig_regions)}")
md_lines.append("")

md_lines.append("## Phase 3: Unclassified Site Analysis")
md_lines.append("")
md_lines.append(f"- On-corridor:  {on_unc} / {on_total_all} = {on_unc_frac:.4f}")
md_lines.append(f"- Off-corridor: {off_unc} / {off_total_all} = {off_unc_frac:.4f}")
if not math.isnan(p3):
    md_lines.append(f"- Fisher exact: OR = {odds3:.3f}, p = {p3:.4f}")
md_lines.append("")

md_lines.append("## Decision Criteria")
md_lines.append("")
md_lines.append("- **CLEAR**: monument_ratio NOT correlated with distance-to-GC (Spearman p > 0.1), OR divergence persists within-region")
md_lines.append("- **CONCERN**: Regional classification varies AND within-region sample sizes too small")
md_lines.append("- **FATAL**: monument_ratio strongly correlates with GC proximity AND within-region divergence disappears")
md_lines.append("")

md_path = os.path.join(OUT_DIR, "RESULTS.md")
with open(md_path, "w") as fh:
    fh.write("\n".join(md_lines))
print(f"  Saved → {md_path}")

print()
print("=" * 70)
print("DONE")
print("=" * 70)
