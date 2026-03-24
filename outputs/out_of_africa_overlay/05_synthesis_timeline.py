#!/usr/bin/env python3
"""
Analysis 5: Full 60,000-Year Synthesis Timeline
=================================================
Directive 11 — Out-of-Africa Migration Overlay

Combines results from all analyses in this directive PLUS Directives 01
and 04 into a single master timeline figure and data file covering
60,000+ years of corridor activity.
"""

import csv, json, math, os, sys
import numpy as np

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "The Great Circle", "outputs", "out_of_africa_overlay")
os.makedirs(OUT_DIR, exist_ok=True)

EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
CORRIDOR_KM = 200

def haversine(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1; dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))

def gc_distance(lat, lon):
    return abs(haversine(POLE_LAT, POLE_LON, lat, lon) - QUARTER_CIRC)

# ============================================================
# LOAD ALL DATA SOURCES
# ============================================================

# 1. Deep-time sites from Analysis 2
DEEP_TIME_SITES = []
dt_path = os.path.join(OUT_DIR, "deep_time_sites.csv")
if os.path.exists(dt_path):
    with open(dt_path) as f:
        for row in csv.DictReader(f):
            DEEP_TIME_SITES.append({
                "name": row["name"],
                "age_bp": int(row["age_bp"]),
                "distance_km": float(row["distance_to_gc_km"]),
                "on_corridor": row["on_corridor"] == "True",
                "source": "deep_time_compilation",
            })

# 2. p3k14c radiocarbon dates (Holocene)
P3K14C = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")

# 3. Key events timeline
KEY_EVENTS = [
    (125000, "Jebel Faya — first AMH tools in Arabia"),
    (74000,  "Toba super-eruption"),
    (65000,  "Madjedbebe — earliest humans in Australia"),
    (55000,  "Mololo Cave — earliest humans in PNG (ON CORRIDOR)"),
    (50000,  "Southern coastal migration wave"),
    (26000,  "Last Glacial Maximum begins"),
    (20000,  "LGM peak — sea level -125m"),
    (12000,  "Younger Dryas — climate shock"),
    (10500,  "Early Holocene campsites (Directive 01)"),
    (9000,   "Agriculture begins — multiple independent origins along corridor"),
    (8000,   "Persian Gulf floods — Arabian segment disrupted"),
    (6500,   "Domestication centers cluster on corridor (p=0.002, Directive 01)"),
    (5000,   "Bronze Age monuments begin (Directive 04)"),
    (4500,   "Pyramid construction — Egypt"),
    (3000,   "Monument-settlement continuity confirmed (p=0.004, Directive 04)"),
]

# 4. Corridor viability from Analysis 3
VIABILITY = {}
viab_path = os.path.join(OUT_DIR, "corridor_viability_by_epoch.json")
if os.path.exists(viab_path):
    with open(viab_path) as f:
        VIABILITY = json.load(f)

# 5. Green Arabia phases
GREEN_ARABIA = [
    (130000, 71000, "MIS 5"),
    (60000,  50000, "MIS 3 phase 1"),
    (44000,  35000, "MIS 3 phase 2"),
    (11000,   5500, "Holocene humid"),
]

# ============================================================
# BUILD TIMELINE DATA
# ============================================================
def run_analysis():
    print("=" * 70)
    print("ANALYSIS 5: FULL 60,000-YEAR SYNTHESIS TIMELINE")
    print("=" * 70)

    # --- Bin all on-corridor sites by millennium ---
    print("\n--- Building temporal bins (1,000-year intervals) ---")
    bins = {}  # {millennium_start: count}
    max_age = 130000

    # Deep-time sites
    for site in DEEP_TIME_SITES:
        if site["distance_km"] <= 500:  # Broader for timeline
            millennium = (site["age_bp"] // 1000) * 1000
            bins[millennium] = bins.get(millennium, 0) + 1

    # p3k14c sites (Holocene portion)
    p3k_count = 0
    if os.path.exists(P3K14C):
        print("  Loading p3k14c for Holocene bins...")
        seen = set()
        with open(P3K14C, "r", encoding="utf-8", errors="replace") as f:
            for row in csv.DictReader(f):
                try:
                    lat = float(row.get("Lat") or row.get("lat", ""))
                    lon = float(row.get("Long") or row.get("lon", ""))
                    age = float(row.get("Age") or row.get("age", ""))
                except (ValueError, TypeError):
                    continue
                if abs(lat) < 0.01 and abs(lon) < 0.01:
                    continue
                if age > 12000:
                    continue  # Only Holocene
                site_key = (round(lat, 2), round(lon, 2))
                if site_key in seen:
                    continue
                seen.add(site_key)
                if gc_distance(lat, lon) <= CORRIDOR_KM:
                    millennium = (int(age) // 1000) * 1000
                    bins[millennium] = bins.get(millennium, 0) + 1
                    p3k_count += 1
        print(f"  Added {p3k_count} on-corridor p3k14c sites")

    # Print bin summary
    print(f"\n  Total bins with sites: {len(bins)}")
    for age in sorted(bins.keys(), reverse=True):
        if bins[age] > 0:
            bar = "#" * min(bins[age], 50)
            print(f"  {age:>7,d} BP: {bins[age]:4d} {bar}")

    # --- Corridor viability interpolation ---
    viability_ages = sorted([int(k) for k in VIABILITY.keys()])
    viability_scores = {int(k): v["mean_viability"] for k, v in VIABILITY.items()}

    # --- Save timeline data ---
    timeline_data = {
        "bins": {str(k): v for k, v in sorted(bins.items())},
        "key_events": [{"age_bp": a, "label": l} for a, l in KEY_EVENTS],
        "corridor_viability": {str(k): v for k, v in viability_scores.items()},
        "green_arabia_phases": [{"start_bp": s, "end_bp": e, "label": l} for s, e, l in GREEN_ARABIA],
        "summary": {
            "total_on_corridor_deep_time": sum(1 for s in DEEP_TIME_SITES if s["distance_km"] <= 500),
            "total_on_corridor_p3k14c": p3k_count,
            "earliest_on_corridor_bp": max((s["age_bp"] for s in DEEP_TIME_SITES if s["distance_km"] <= 200), default=0),
            "latest_on_corridor_bp": min((s["age_bp"] for s in DEEP_TIME_SITES if s["distance_km"] <= 200), default=0),
            "occupied_millennia": len([k for k, v in bins.items() if v > 0]),
        },
    }

    with open(os.path.join(OUT_DIR, "timeline_data.json"), "w") as f:
        json.dump(timeline_data, f, indent=2)
    print(f"\nTimeline data saved.")

    # --- Master timeline figure ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyBboxPatch
        import matplotlib.patches as mpatches

        fig, axes = plt.subplots(2, 1, figsize=(20, 12), sharex=True,
                                  gridspec_kw={'height_ratios': [3, 1]})

        # === TOP PANEL: Site density ===
        ax1 = axes[0]

        # Sort bins and plot
        all_ages = sorted(bins.keys())
        counts = [bins.get(a, 0) for a in all_ages]

        # Split into deep-time and Holocene for different colors
        dt_ages = [a for a in all_ages if a > 12000]
        dt_counts = [bins.get(a, 0) for a in dt_ages]
        hol_ages = [a for a in all_ages if a <= 12000]
        hol_counts = [bins.get(a, 0) for a in hol_ages]

        ax1.bar(dt_ages, dt_counts, width=800, color='#d62728', alpha=0.7,
                label='Deep-time sites (compiled)')
        ax1.bar(hol_ages, hol_counts, width=800, color='#1f77b4', alpha=0.7,
                label='Holocene sites (p3k14c)')

        # Key events as annotations
        y_max = max(counts) if counts else 1
        for age, label in KEY_EVENTS:
            if age <= max_age:
                ax1.axvline(age, color='gray', linestyle=':', alpha=0.3, linewidth=0.5)
                # Stagger labels
                y_pos = y_max * 0.85 - (KEY_EVENTS.index((age, label)) % 4) * y_max * 0.12
                ax1.annotate(label, (age, y_pos), fontsize=6.5, ha='left',
                             rotation=0, alpha=0.8,
                             arrowprops=dict(arrowstyle='-', color='gray', alpha=0.3),
                             xytext=(age + 500, y_pos))

        # Green Arabia shading
        for s, e, l in GREEN_ARABIA:
            ax1.axvspan(e, s, alpha=0.06, color='green')

        ax1.set_ylabel("On-Corridor Sites per 1,000 Years", fontsize=12)
        ax1.set_title("THE GREAT CIRCLE: 130,000 YEARS OF CORRIDOR ACTIVITY", fontsize=16, fontweight='bold')
        ax1.legend(loc='upper left', fontsize=10)
        ax1.grid(True, alpha=0.2)
        ax1.set_xlim(0, 130000)

        # === BOTTOM PANEL: Corridor viability ===
        ax2 = axes[1]

        if viability_ages:
            v_ages = sorted(viability_scores.keys())
            v_scores = [viability_scores[a] for a in v_ages]
            ax2.plot(v_ages, v_scores, 'go-', linewidth=2, markersize=8)
            ax2.fill_between(v_ages, v_scores, alpha=0.2, color='green')
            ax2.axhline(1, color='orange', linestyle='--', alpha=0.5, linewidth=1)
            ax2.set_ylim(0, 3.5)
        else:
            ax2.text(65000, 1.5, "(Viability data from Analysis 3)", ha='center', fontsize=12)

        # Green Arabia shading
        for s, e, l in GREEN_ARABIA:
            ax2.axvspan(e, s, alpha=0.06, color='green')
            ax2.text((s + e) / 2, 3.2, l, fontsize=7, ha='center', color='green', alpha=0.6)

        ax2.set_ylabel("Corridor Viability\n(0=blocked, 3=excellent)", fontsize=10)
        ax2.set_xlabel("Years Before Present", fontsize=12)
        ax2.grid(True, alpha=0.2)
        ax2.invert_xaxis()

        # Add epoch labels at top
        epoch_labels = [
            (130000, 71000, "MIS 5\n(interglacial)"),
            (71000, 27000, "MIS 3-4\n(glacial)"),
            (27000, 12000, "LGM →\ndeglaciation"),
            (12000, 5000, "Early\nHolocene"),
            (5000, 0, "Bronze\nAge+"),
        ]
        for s, e, l in epoch_labels:
            ax1.text((s + e) / 2, y_max * 1.05, l, fontsize=8, ha='center',
                     style='italic', color='dimgray')

        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, "master_timeline_60k.png"),
                    dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"Master timeline figure saved to {OUT_DIR}/master_timeline_60k.png")

    except ImportError:
        print("matplotlib not available — skipping figure")

    return timeline_data


if __name__ == "__main__":
    run_analysis()
