#!/usr/bin/env python3
"""
Study 4: Monument Type Decomposition (Mortuary vs. Ceremonial vs. Administrative)
==================================================================================
Decomposes the monument-corridor divergence signal by functional sub-type.
Tests whether the signal is a general monument phenomenon or driven by
specific categories -- particularly mortuary sites in the Egypt segment.

Phases:
  1. Sub-type classification of all Pleiades monumental sites
  2. Per-type divergence (Z-scores, MC with 1000 trials)
  3. Temporal cross-check (pre-2000 BCE, 2000 BCE-0 CE, 0 CE-present)
  4. Egypt/Levant segment isolation
"""

import sys
import os
import math
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout = os.fdopen(sys.stdout.fileno(), "w", buffering=1)

# ============================================================
# CONSTANTS
# ============================================================
POLE_LAT = 59.682122
POLE_LON = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
THRESHOLD_KM = 50
MC_TRIALS = 1000

PLEIADES_CSV = "/Users/elliotallan/megalith_site_research/data/pleiades/pleiades-places-latest.csv"
OUTPUT_DIR = "/Users/elliotallan/megalith_site_research/outputs/next_wave/monument_types"

np.random.seed(42)

# ============================================================
# TYPE SETS
# ============================================================
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

MORTUARY_TYPES = {
    "tomb", "tumulus", "necropolis", "mausoleum", "pyramid",
    "mastaba", "catacomb", "cemetery", "burial",
}
CEREMONIAL_TYPES = {
    "temple", "sanctuary", "shrine", "altar", "church",
    "monastery", "synagogue", "sacred-grove",
    # Also match Pleiades variants
    "temple-2", "church-2",
}
ADMINISTRATIVE_TYPES = {
    "palace", "fortress", "fort", "castle", "wall",
    "aqueduct", "bridge", "bath", "theater", "amphitheater",
    "arch", "hippodrome", "stadium", "forum",
    # Pleiades variants
    "amphitheatre", "theatre",
}

# ============================================================
# GEOMETRY
# ============================================================
def haversine_km(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return 2 * EARTH_R_KM * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def dist_from_gc(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QUARTER_CIRC)


# ============================================================
# DATA LOADING
# ============================================================
def load_pleiades():
    """Load Pleiades CSV and return list of dicts with parsed fields."""
    sites = []
    with open(PLEIADES_CSV, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["reprLat"])
                lon = float(row["reprLong"])
            except (ValueError, TypeError):
                continue
            ft_raw = row.get("featureTypes", "")
            feature_types = {t.strip() for t in ft_raw.split(",") if t.strip()}
            try:
                min_date = float(row["minDate"]) if row["minDate"] else None
            except ValueError:
                min_date = None
            try:
                max_date = float(row["maxDate"]) if row["maxDate"] else None
            except ValueError:
                max_date = None
            sites.append({
                "lat": lat,
                "lon": lon,
                "feature_types": feature_types,
                "min_date": min_date,
                "max_date": max_date,
                "title": row.get("title", ""),
                "dist_gc": dist_from_gc(lat, lon),
            })
    return sites


def classify_site(feature_types):
    """Classify a site into sub-type categories. Returns set of categories."""
    cats = set()
    if feature_types & MORTUARY_TYPES:
        cats.add("mortuary")
    if feature_types & CEREMONIAL_TYPES:
        cats.add("ceremonial")
    if feature_types & ADMINISTRATIVE_TYPES:
        cats.add("administrative")
    return cats


# ============================================================
# MONTE CARLO
# ============================================================
def mc_on_corridor(pool_dists, observed_n, observed_on, n_trials=MC_TRIALS):
    """
    Monte Carlo: randomly sample observed_n sites from pool_dists,
    count how many fall within THRESHOLD_KM. Return Z-score and p-value.
    """
    pool = np.array(pool_dists)
    counts = np.zeros(n_trials)
    for i in range(n_trials):
        sample = np.random.choice(pool, size=observed_n, replace=False)
        counts[i] = np.sum(sample < THRESHOLD_KM)
    mean_mc = np.mean(counts)
    std_mc = np.std(counts)
    z = (observed_on - mean_mc) / std_mc if std_mc > 0 else 0.0
    p = np.mean(counts >= observed_on)
    return z, p, mean_mc, std_mc


# ============================================================
# PHASE 1: Sub-Type Classification
# ============================================================
def phase1(sites):
    print("=" * 60)
    print("PHASE 1: Sub-Type Classification")
    print("=" * 60)

    monuments = [s for s in sites if s["feature_types"] & MONUMENTAL_TYPES]
    settlements = [s for s in sites if s["feature_types"] & SETTLEMENT_TYPES]

    # Classify monuments into sub-types
    mortuary = [s for s in monuments if "mortuary" in classify_site(s["feature_types"])]
    ceremonial = [s for s in monuments if "ceremonial" in classify_site(s["feature_types"])]
    administrative = [s for s in monuments if "administrative" in classify_site(s["feature_types"])]
    # Other monumental: in MONUMENTAL_TYPES but not in any sub-type
    other_mon = [s for s in monuments if not classify_site(s["feature_types"])]

    print(f"  Total Pleiades sites loaded: {len(sites)}")
    print(f"  Monumental sites:    {len(monuments)}")
    print(f"    Mortuary:          {len(mortuary)}")
    print(f"    Ceremonial:        {len(ceremonial)}")
    print(f"    Administrative:    {len(administrative)}")
    print(f"    Other monumental:  {len(other_mon)}")
    print(f"  Settlement sites:    {len(settlements)}")
    print()

    return {
        "monuments": monuments,
        "settlements": settlements,
        "mortuary": mortuary,
        "ceremonial": ceremonial,
        "administrative": administrative,
        "other_monumental": other_mon,
    }


# ============================================================
# PHASE 2: Per-Type Divergence
# ============================================================
def phase2(groups, all_monument_dists):
    print("=" * 60)
    print("PHASE 2: Per-Type Divergence (MC null = full monument pool)")
    print("=" * 60)

    results = {}
    for label in ["monuments", "settlements", "mortuary", "ceremonial", "administrative", "other_monumental"]:
        grp = groups[label]
        n_total = len(grp)
        if n_total == 0:
            print(f"  {label}: N=0, skipping")
            results[label] = {"n": 0, "on_corridor": 0, "z": None, "p": None}
            continue
        on_corr = sum(1 for s in grp if s["dist_gc"] < THRESHOLD_KM)
        frac = on_corr / n_total

        # For settlements, use settlement pool as null; for monument sub-types, use monument pool
        if label == "settlements":
            pool_dists = [s["dist_gc"] for s in grp]
            # MC: sample from settlement pool
            z, p, mc_mean, mc_std = mc_on_corridor(pool_dists, n_total, on_corr)
        else:
            # Use FULL monument pool as null
            z, p, mc_mean, mc_std = mc_on_corridor(all_monument_dists, min(n_total, len(all_monument_dists)), on_corr)

        print(f"  {label:20s}  N={n_total:5d}  on-corridor={on_corr:4d} ({frac:.3f})  "
              f"Z={z:+.2f}  p={p:.4f}  MC_mean={mc_mean:.1f}")
        results[label] = {
            "n": n_total,
            "on_corridor": on_corr,
            "fraction": round(frac, 4),
            "z_score": round(z, 3),
            "p_value": round(p, 4),
            "mc_mean": round(mc_mean, 2),
            "mc_std": round(mc_std, 3),
        }

    # Compute divergence relative to settlement Z
    sett_z = results["settlements"].get("z_score")
    if sett_z is not None:
        print(f"\n  Divergence (sub_type_Z - settlement_Z, settlement_Z = {sett_z:+.2f}):")
        for label in ["monuments", "mortuary", "ceremonial", "administrative", "other_monumental"]:
            z = results[label].get("z_score")
            if z is not None:
                div = round(z - sett_z, 3)
                results[label]["divergence_vs_settlement"] = div
                print(f"    {label:20s}  divergence = {div:+.3f}")
    print()
    return results


# ============================================================
# PHASE 3: Temporal Cross-Check
# ============================================================
def phase3(groups, all_monument_dists):
    print("=" * 60)
    print("PHASE 3: Temporal Cross-Check")
    print("=" * 60)

    periods = {
        "pre_2000BCE": lambda d: d is not None and d < -2000,
        "2000BCE_0CE": lambda d: d is not None and -2000 <= d <= 0,
        "0CE_present": lambda d: d is not None and d > 0,
    }

    temporal_results = {}
    for label in ["mortuary", "ceremonial", "administrative"]:
        grp = groups[label]
        temporal_results[label] = {}
        print(f"\n  --- {label} ---")
        for period_name, date_filter in periods.items():
            subset = [s for s in grp if date_filter(s["min_date"])]
            n = len(subset)
            if n < 5:
                print(f"    {period_name:15s}  N={n} (too few, skipping)")
                temporal_results[label][period_name] = {"n": n, "skipped": True}
                continue
            on_corr = sum(1 for s in subset if s["dist_gc"] < THRESHOLD_KM)
            frac = on_corr / n
            z, p, mc_mean, mc_std = mc_on_corridor(all_monument_dists, min(n, len(all_monument_dists)), on_corr)
            print(f"    {period_name:15s}  N={n:4d}  on-corridor={on_corr:3d} ({frac:.3f})  Z={z:+.2f}  p={p:.4f}")
            temporal_results[label][period_name] = {
                "n": n,
                "on_corridor": on_corr,
                "fraction": round(frac, 4),
                "z_score": round(z, 3),
                "p_value": round(p, 4),
            }
    print()
    return temporal_results


# ============================================================
# PHASE 4: Egypt/Levant Segment Isolation
# ============================================================
def phase4(groups):
    print("=" * 60)
    print("PHASE 4: Egypt/Levant Segment Isolation")
    print("=" * 60)

    # Egypt/Levant box: lat 25-35, lon 25-40, within 200km of GC
    def in_egypt_box(s):
        return (25 <= s["lat"] <= 35 and 25 <= s["lon"] <= 40 and s["dist_gc"] < 200)

    egypt_mon = [s for s in groups["monuments"] if in_egypt_box(s)]
    egypt_on_corr = [s for s in egypt_mon if s["dist_gc"] < THRESHOLD_KM]

    # Classify the on-corridor Egypt monuments
    egypt_mortuary_on = [s for s in egypt_on_corr if "mortuary" in classify_site(s["feature_types"])]
    egypt_ceremonial_on = [s for s in egypt_on_corr if "ceremonial" in classify_site(s["feature_types"])]
    egypt_admin_on = [s for s in egypt_on_corr if "administrative" in classify_site(s["feature_types"])]
    egypt_other_on = [s for s in egypt_on_corr if not classify_site(s["feature_types"])]

    n_on = len(egypt_on_corr)
    print(f"  Egypt/Levant box (lat 25-35, lon 25-40, <200km from GC):")
    print(f"    Monumental sites in box:       {len(egypt_mon)}")
    print(f"    On-corridor (<50km):           {n_on}")
    if n_on > 0:
        print(f"      Mortuary:    {len(egypt_mortuary_on):3d}  ({len(egypt_mortuary_on)/n_on:.1%})")
        print(f"      Ceremonial:  {len(egypt_ceremonial_on):3d}  ({len(egypt_ceremonial_on)/n_on:.1%})")
        print(f"      Administrative: {len(egypt_admin_on):3d}  ({len(egypt_admin_on)/n_on:.1%})")
        print(f"      Other:       {len(egypt_other_on):3d}  ({len(egypt_other_on)/n_on:.1%})")

    # List some example mortuary sites
    print(f"\n  Example on-corridor mortuary sites in Egypt segment:")
    for s in sorted(egypt_mortuary_on, key=lambda x: x["dist_gc"])[:10]:
        print(f"    {s['title'][:50]:50s}  dist={s['dist_gc']:.1f}km  types={','.join(s['feature_types'])}")

    egypt_results = {
        "monuments_in_box": len(egypt_mon),
        "on_corridor": n_on,
        "mortuary_on_corridor": len(egypt_mortuary_on),
        "ceremonial_on_corridor": len(egypt_ceremonial_on),
        "administrative_on_corridor": len(egypt_admin_on),
        "other_on_corridor": len(egypt_other_on),
        "mortuary_fraction": round(len(egypt_mortuary_on) / n_on, 4) if n_on > 0 else None,
    }
    print()
    return egypt_results


# ============================================================
# PLOTTING
# ============================================================
def plot_divergence_bar(phase2_results):
    """Bar chart of Z-scores by sub-type."""
    labels = ["monuments", "mortuary", "ceremonial", "administrative", "other_monumental", "settlements"]
    display = ["All\nMonuments", "Mortuary", "Ceremonial", "Admin.", "Other\nMonum.", "Settlements"]
    colors = ["#2196F3", "#E53935", "#FF9800", "#4CAF50", "#9E9E9E", "#78909C"]

    zs = [phase2_results[l].get("z_score", 0) or 0 for l in labels]

    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.bar(display, zs, color=colors, edgecolor="black", linewidth=0.5)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.axhline(1.96, color="red", linewidth=0.7, linestyle="--", alpha=0.5, label="Z=1.96 (p<0.05)")
    ax.axhline(-1.96, color="red", linewidth=0.7, linestyle="--", alpha=0.5)

    for bar, z in zip(bars, zs):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.15 * (1 if z >= 0 else -1),
                f"Z={z:+.1f}", ha="center", va="bottom" if z >= 0 else "top", fontsize=9, fontweight="bold")

    ax.set_ylabel("Z-score (on-corridor enrichment vs. monument null)", fontsize=10)
    ax.set_title("Study 4: Monument Type Decomposition\nOn-Corridor Enrichment by Sub-Type", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "divergence_by_type.png"), dpi=150)
    plt.close()
    print("  Saved divergence_by_type.png")


def plot_temporal(temporal_results):
    """Grouped bar chart: Z-score by sub-type and period."""
    periods = ["pre_2000BCE", "2000BCE_0CE", "0CE_present"]
    period_labels = ["Pre-2000 BCE", "2000 BCE - 0 CE", "0 CE - Present"]
    sub_types = ["mortuary", "ceremonial", "administrative"]
    colors = ["#E53935", "#FF9800", "#4CAF50"]

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(periods))
    width = 0.25

    for i, (st, color) in enumerate(zip(sub_types, colors)):
        zs = []
        for p in periods:
            entry = temporal_results.get(st, {}).get(p, {})
            z = entry.get("z_score", None)
            zs.append(z if z is not None else 0)
        bars = ax.bar(x + i * width, zs, width, label=st.capitalize(), color=color, edgecolor="black", linewidth=0.4)
        for bar, z_val, p_name in zip(bars, zs, periods):
            entry = temporal_results.get(st, {}).get(p_name, {})
            n = entry.get("n", 0)
            if entry.get("skipped"):
                ax.text(bar.get_x() + bar.get_width() / 2, 0.1, f"N={n}", ha="center", fontsize=7, style="italic")
            else:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.1,
                        f"Z={z_val:+.1f}\nN={n}", ha="center", va="bottom", fontsize=7)

    ax.axhline(0, color="black", linewidth=0.8)
    ax.axhline(1.96, color="red", linewidth=0.7, linestyle="--", alpha=0.5)
    ax.set_xticks(x + width)
    ax.set_xticklabels(period_labels)
    ax.set_ylabel("Z-score")
    ax.set_title("Study 4: Temporal Cross-Check by Sub-Type", fontsize=12, fontweight="bold")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "temporal_by_type.png"), dpi=150)
    plt.close()
    print("  Saved temporal_by_type.png")


def plot_egypt_pie(egypt_results):
    """Pie chart of on-corridor Egypt monument sub-types."""
    labels = ["Mortuary", "Ceremonial", "Administrative", "Other"]
    sizes = [
        egypt_results["mortuary_on_corridor"],
        egypt_results["ceremonial_on_corridor"],
        egypt_results["administrative_on_corridor"],
        egypt_results["other_on_corridor"],
    ]
    colors = ["#E53935", "#FF9800", "#4CAF50", "#9E9E9E"]

    # Filter out zero slices
    filtered = [(l, s, c) for l, s, c in zip(labels, sizes, colors) if s > 0]
    if not filtered:
        print("  No Egypt on-corridor monuments to plot.")
        return
    labels_f, sizes_f, colors_f = zip(*filtered)

    fig, ax = plt.subplots(figsize=(6, 6))
    wedges, texts, autotexts = ax.pie(sizes_f, labels=labels_f, colors=colors_f, autopct="%1.0f%%",
                                       startangle=90, textprops={"fontsize": 11})
    ax.set_title("Egypt/Levant Segment: On-Corridor Monument Types\n(lat 25-35, lon 25-40, <50km from GC)",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "egypt_segment_pie.png"), dpi=150)
    plt.close()
    print("  Saved egypt_segment_pie.png")


# ============================================================
# RESULTS MARKDOWN
# ============================================================
def write_results_md(phase2_results, temporal_results, egypt_results):
    sett_z = phase2_results["settlements"].get("z_score", 0) or 0
    mon_z = phase2_results["monuments"].get("z_score", 0) or 0
    mort_z = phase2_results["mortuary"].get("z_score", 0) or 0
    cere_z = phase2_results["ceremonial"].get("z_score", 0) or 0
    admin_z = phase2_results["administrative"].get("z_score", 0) or 0

    # Determine which sub-type drives the signal
    sub_zs = {"mortuary": mort_z, "ceremonial": cere_z, "administrative": admin_z}
    driver = max(sub_zs, key=sub_zs.get)
    driver_z = sub_zs[driver]

    mort_frac = egypt_results.get("mortuary_fraction")
    mort_frac_str = f"{mort_frac:.0%}" if mort_frac is not None else "N/A"

    md = f"""# Study 4: Monument Type Decomposition
## Mortuary vs. Ceremonial vs. Administrative

### Key Question
Is the monument-corridor divergence signal a general monument phenomenon,
or is it driven by specific sub-types (particularly mortuary sites)?

### Phase 2: Per-Type Z-Scores

| Sub-Type | N | On-Corridor | Z-Score | Divergence vs. Settlement |
|----------|---|-------------|---------|---------------------------|
| All Monuments | {phase2_results['monuments']['n']} | {phase2_results['monuments']['on_corridor']} | {mon_z:+.2f} | {phase2_results['monuments'].get('divergence_vs_settlement', 'N/A')} |
| Mortuary | {phase2_results['mortuary']['n']} | {phase2_results['mortuary']['on_corridor']} | {mort_z:+.2f} | {phase2_results['mortuary'].get('divergence_vs_settlement', 'N/A')} |
| Ceremonial | {phase2_results['ceremonial']['n']} | {phase2_results['ceremonial']['on_corridor']} | {cere_z:+.2f} | {phase2_results['ceremonial'].get('divergence_vs_settlement', 'N/A')} |
| Administrative | {phase2_results['administrative']['n']} | {phase2_results['administrative']['on_corridor']} | {admin_z:+.2f} | {phase2_results['administrative'].get('divergence_vs_settlement', 'N/A')} |
| Settlements (baseline) | {phase2_results['settlements']['n']} | {phase2_results['settlements']['on_corridor']} | {sett_z:+.2f} | -- |

### Phase 3: Temporal Cross-Check

"""
    for label in ["mortuary", "ceremonial", "administrative"]:
        md += f"**{label.capitalize()}:**\n\n"
        md += "| Period | N | On-Corridor | Z-Score |\n|--------|---|-------------|--------|\n"
        for period in ["pre_2000BCE", "2000BCE_0CE", "0CE_present"]:
            entry = temporal_results.get(label, {}).get(period, {})
            if entry.get("skipped"):
                md += f"| {period} | {entry.get('n', 0)} | -- | (too few) |\n"
            else:
                md += f"| {period} | {entry.get('n', 0)} | {entry.get('on_corridor', 0)} | {entry.get('z_score', 0):+.2f} |\n"
        md += "\n"

    md += f"""### Phase 4: Egypt/Levant Segment

- Monuments in Egypt box (lat 25-35, lon 25-40, <200km): {egypt_results['monuments_in_box']}
- On-corridor (<50km): {egypt_results['on_corridor']}
- Mortuary fraction of on-corridor: **{mort_frac_str}**

### Conclusion

The strongest on-corridor enrichment comes from **{driver}** sites (Z={driver_z:+.2f}),
compared to the overall monument Z={mon_z:+.2f} and settlement baseline Z={sett_z:+.2f}.

"""
    if driver == "mortuary" and driver_z > mon_z:
        md += ("The signal is disproportionately driven by mortuary monuments. "
               "In the Egypt/Levant segment, mortuary sites account for "
               f"{mort_frac_str} of on-corridor monuments, confirming that "
               "the divergence signal is not a general monument phenomenon "
               "but is concentrated in funerary architecture.\n")
    elif abs(mort_z - cere_z) < 0.5 and abs(mort_z - admin_z) < 0.5:
        md += ("The signal is broadly distributed across all monument sub-types, "
               "suggesting it is a general monument phenomenon rather than being "
               "driven by any single functional category.\n")
    else:
        md += (f"The {driver} sub-type shows the strongest signal, but the "
               "pattern varies across categories. See temporal cross-check "
               "for period-specific details.\n")

    path = os.path.join(OUTPUT_DIR, "RESULTS.md")
    with open(path, "w") as f:
        f.write(md)
    print(f"  Saved RESULTS.md")


# ============================================================
# MAIN
# ============================================================
def main():
    print("Study 4: Monument Type Decomposition")
    print("=" * 60)

    # Load data
    print("Loading Pleiades data...")
    sites = load_pleiades()
    print(f"  Loaded {len(sites)} sites with valid coordinates.\n")

    # Phase 1
    groups = phase1(sites)

    # Precompute GC distances for all monuments (used as MC null pool)
    all_monument_dists = [s["dist_gc"] for s in groups["monuments"]]
    print(f"  Monument pool for MC null: {len(all_monument_dists)} sites")
    print(f"  Monument pool on-corridor rate: {sum(1 for d in all_monument_dists if d < THRESHOLD_KM) / len(all_monument_dists):.4f}\n")

    # Phase 2
    phase2_results = phase2(groups, all_monument_dists)

    # Phase 3
    temporal_results = phase3(groups, all_monument_dists)

    # Phase 4
    egypt_results = phase4(groups)

    # Save JSON
    all_results = {
        "study": "S4_monument_type_decomposition",
        "mc_trials": MC_TRIALS,
        "threshold_km": THRESHOLD_KM,
        "phase2_per_type": phase2_results,
        "phase3_temporal": temporal_results,
        "phase4_egypt_segment": egypt_results,
    }
    json_path = os.path.join(OUTPUT_DIR, "results.json")
    with open(json_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"Saved results.json")

    # Plots
    print("\nGenerating plots...")
    plot_divergence_bar(phase2_results)
    plot_temporal(temporal_results)
    plot_egypt_pie(egypt_results)

    # RESULTS.md
    print("\nWriting RESULTS.md...")
    write_results_md(phase2_results, temporal_results, egypt_results)

    print("\nStudy 4 complete.")


if __name__ == "__main__":
    main()
