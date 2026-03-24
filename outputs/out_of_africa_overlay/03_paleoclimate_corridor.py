#!/usr/bin/env python3
"""
Analysis 3: Paleoclimate Corridor Persistence
===============================================
Directive 11 — Out-of-Africa Migration Overlay

Assesses whether the Great Circle corridor was continuously habitable
across major climate epochs from OIS-4 (~74,000 BP) through the Holocene.
Samples the overland portion of the circle and evaluates viability at
each epoch based on sea level, aridity, and known "Green Arabia" phases.

Uses:
  - Lambeck et al. (2014) sea level curve
  - Green Arabia phases from Groucutt & Petraglia (2012)
  - Modern elevation as baseline (SRTM where available)
"""

import csv, json, math, os, sys
import numpy as np

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "The Great Circle", "outputs", "out_of_africa_overlay")
os.makedirs(OUT_DIR, exist_ok=True)

np.random.seed(42)

# ============================================================
# GEOMETRY
# ============================================================
def haversine(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))


def gc_distance(lat, lon):
    d = haversine(POLE_LAT, POLE_LON, lat, lon)
    return abs(d - QUARTER_CIRC)


def points_on_great_circle(n=360):
    """Generate n evenly spaced points on the Great Circle."""
    pole_lat_r = math.radians(POLE_LAT)
    pole_lon_r = math.radians(POLE_LON)
    points = []
    for i in range(n):
        theta = math.radians(i * 360.0 / n)
        lat = math.asin(math.sin(pole_lat_r) * math.cos(math.pi/2) +
                        math.cos(pole_lat_r) * math.sin(math.pi/2) * math.cos(theta))
        lon = pole_lon_r + math.atan2(
            math.sin(theta) * math.sin(math.pi/2) * math.cos(pole_lat_r),
            math.cos(math.pi/2) - math.sin(pole_lat_r) * math.sin(lat)
        )
        points.append((math.degrees(lat), math.degrees(lon)))
    return points


# ============================================================
# SEA LEVEL CURVE (Lambeck et al. 2014, simplified)
# ============================================================
# Approximate global mean sea level relative to present (meters)
# Key data points from Lambeck et al. 2014, Table S1
SEA_LEVEL_CURVE = {
    # age_bp: sea_level_meters_relative_to_present
    0:      0,
    5000:  -2,
    8000:  -15,
    10000: -40,
    12000: -60,
    15000: -90,
    18000: -120,
    20000: -125,  # LGM minimum
    22000: -120,
    25000: -100,
    30000: -80,
    35000: -70,
    40000: -70,
    45000: -65,
    50000: -70,
    55000: -80,
    60000: -80,
    65000: -60,  # OIS-4 start
    70000: -70,
    74000: -80,  # Toba eruption / OIS-4 peak
    80000: -50,
    85000: -20,  # MIS 5a warm phase
    90000: -30,
    100000: -15,
    110000: -20,
    120000: +5,  # MIS 5e interglacial peak
    125000: +6,
    130000: 0,
}

def interpolate_sea_level(age_bp):
    """Linearly interpolate sea level from the curve."""
    ages = sorted(SEA_LEVEL_CURVE.keys())
    if age_bp <= ages[0]:
        return SEA_LEVEL_CURVE[ages[0]]
    if age_bp >= ages[-1]:
        return SEA_LEVEL_CURVE[ages[-1]]
    for i in range(len(ages) - 1):
        if ages[i] <= age_bp <= ages[i+1]:
            frac = (age_bp - ages[i]) / (ages[i+1] - ages[i])
            return SEA_LEVEL_CURVE[ages[i]] * (1 - frac) + SEA_LEVEL_CURVE[ages[i+1]] * frac
    return 0


# ============================================================
# GREEN ARABIA PHASES (Groucutt & Petraglia 2012)
# ============================================================
GREEN_ARABIA_PHASES = [
    # (start_bp, end_bp, label)
    (130000, 71000, "MIS 5 — Full interglacial (multiple humid pulses)"),
    (60000,  50000, "MIS 3 humid phase 1"),
    (44000,  35000, "MIS 3 humid phase 2"),
    (11000,   5500, "Early-Mid Holocene humid phase (Green Arabia)"),
]

def is_green_arabia(age_bp):
    for start, end, _ in GREEN_ARABIA_PHASES:
        if end <= age_bp <= start:
            return True
    return False


# ============================================================
# PERSIAN GULF DRY PERIODS
# ============================================================
# The Gulf was dry when sea level was below approximately -80m
# (Gulf maximum depth ~90m, but much is shallower)
GULF_THRESHOLD = -60  # Conservative: most of Gulf exposed below -60m

def is_gulf_dry(age_bp):
    sl = interpolate_sea_level(age_bp)
    return sl < GULF_THRESHOLD


# ============================================================
# CORRIDOR SEGMENTS & VIABILITY
# ============================================================
CORRIDOR_SEGMENTS = [
    {
        "name": "East Africa (departure)",
        "lon_range": (30, 45),
        "lat_range": (-10, 15),
        "viability_notes": "Always habitable — human origin zone",
        "bottleneck": False,
    },
    {
        "name": "Arabian crossing (Bab el-Mandeb → Dhofar)",
        "lon_range": (43, 58),
        "lat_range": (10, 28),
        "viability_notes": "Requires Green Arabia; Bab el-Mandeb strait narrower at low sea level",
        "bottleneck": True,
    },
    {
        "name": "Persian Gulf / Strait of Hormuz",
        "lon_range": (48, 57),
        "lat_range": (24, 30),
        "viability_notes": "Gulf was dry land during LGM and OIS-4; passable during low sea level",
        "bottleneck": True,
    },
    {
        "name": "Makran coast (Iran → Pakistan)",
        "lon_range": (57, 67),
        "lat_range": (23, 27),
        "viability_notes": "Coastal route; wider continental shelf at low sea level",
        "bottleneck": False,
    },
    {
        "name": "Indian subcontinent",
        "lon_range": (67, 80),
        "lat_range": (8, 25),
        "viability_notes": "Always habitable; monsoon variations affect productivity",
        "bottleneck": False,
    },
    {
        "name": "SE Asia (Sunda Shelf)",
        "lon_range": (95, 120),
        "lat_range": (-10, 20),
        "viability_notes": "Sunda Shelf was dry land at low sea level — more habitable, not less",
        "bottleneck": False,
    },
    {
        "name": "Sahul crossing (Timor → Australia/PNG)",
        "lon_range": (120, 150),
        "lat_range": (-15, 0),
        "viability_notes": "Water crossings always required (~70-90 km); narrower at low sea level",
        "bottleneck": True,
    },
]


# ============================================================
# ANALYSIS
# ============================================================
def run_analysis():
    print("=" * 70)
    print("ANALYSIS 3: PALEOCLIMATE CORRIDOR PERSISTENCE")
    print("=" * 70)

    # --- Key epochs ---
    KEY_EPOCHS = [
        (0,     "Present"),
        (5000,  "Mid Holocene"),
        (8000,  "Persian Gulf flooding / early agriculture"),
        (12000, "Younger Dryas / early Holocene"),
        (20000, "Last Glacial Maximum (LGM)"),
        (35000, "MIS 3 — warm interstadial"),
        (50000, "MIS 3 — early occupation"),
        (65000, "OIS-4 — initial southern dispersal"),
        (74000, "Toba eruption epoch"),
        (85000, "MIS 5a — warm phase"),
        (125000, "MIS 5e — last interglacial peak"),
    ]

    # --- Sea level at each epoch ---
    print("\n--- Sea Level at Key Epochs ---")
    sl_results = {}
    for age, label in KEY_EPOCHS:
        sl = interpolate_sea_level(age)
        green = is_green_arabia(age)
        gulf_dry = is_gulf_dry(age)
        sl_results[age] = {
            "label": label,
            "sea_level_m": round(sl, 1),
            "green_arabia": green,
            "persian_gulf_dry": gulf_dry,
        }
        print(f"  {age:>7,d} BP  ({label:40s})  SL: {sl:+6.1f} m  "
              f"Green Arabia: {'YES' if green else 'no '}  "
              f"Gulf dry: {'YES' if gulf_dry else 'no '}")

    # --- Segment viability at each epoch ---
    print("\n--- Corridor Segment Viability by Epoch ---")
    segment_viability = {}
    for seg in CORRIDOR_SEGMENTS:
        seg_name = seg["name"]
        segment_viability[seg_name] = {}
        for age, label in KEY_EPOCHS:
            sl = interpolate_sea_level(age)
            green = is_green_arabia(age)
            gulf_dry = is_gulf_dry(age)

            # Viability score: 0=impassable, 1=marginal, 2=good, 3=excellent
            if "Arabia" in seg_name:
                if green:
                    score = 3  # Green Arabia = excellent
                elif gulf_dry:
                    score = 1  # Can traverse dry regions with coast route
                else:
                    score = 0  # Hyper-arid without green phase
            elif "Persian Gulf" in seg_name:
                if gulf_dry:
                    score = 3  # Dry land — excellent
                elif sl < -30:
                    score = 2  # Partially exposed
                else:
                    score = 1  # Underwater but Hormuz strait passable
            elif "Sahul" in seg_name:
                if sl < -60:
                    score = 2  # Shorter water crossings
                elif sl < -30:
                    score = 1  # Moderate crossings
                else:
                    score = 1  # Always requires water crossing
            elif "Sunda" in seg_name or "SE Asia" in seg_name:
                if sl < -40:
                    score = 3  # Sunda shelf = huge landmass
                else:
                    score = 2  # Island hopping viable
            else:
                score = 2  # Generally passable (Africa, India, Makran)

            segment_viability[seg_name][age] = {
                "label": label,
                "viability_score": score,
            }

    # Print viability matrix
    print(f"\n  {'Segment':45s}", end="")
    for age, _ in KEY_EPOCHS:
        print(f" {age//1000:4d}k", end="")
    print()
    print("  " + "-" * (45 + 6 * len(KEY_EPOCHS)))
    viability_chars = {0: "X", 1: ".", 2: "o", 3: "O"}
    for seg in CORRIDOR_SEGMENTS:
        seg_name = seg["name"]
        print(f"  {seg_name:45s}", end="")
        for age, _ in KEY_EPOCHS:
            score = segment_viability[seg_name][age]["viability_score"]
            print(f"    {viability_chars[score]} ", end="")
        print()

    print(f"\n  Legend: O=excellent  o=good  .=marginal  X=impassable")

    # --- Overall corridor viability score per epoch ---
    print("\n--- Overall Corridor Viability Score per Epoch ---")
    corridor_scores = {}
    for age, label in KEY_EPOCHS:
        scores = [segment_viability[seg["name"]][age]["viability_score"]
                  for seg in CORRIDOR_SEGMENTS]
        min_score = min(scores)
        mean_score = np.mean(scores)
        # Bottleneck score = minimum score across bottleneck segments
        bottleneck_scores = [segment_viability[seg["name"]][age]["viability_score"]
                             for seg in CORRIDOR_SEGMENTS if seg.get("bottleneck")]
        bottleneck_min = min(bottleneck_scores) if bottleneck_scores else min_score

        viable = min_score > 0
        corridor_scores[age] = {
            "label": label,
            "mean_viability": round(float(mean_score), 2),
            "min_viability": int(min_score),
            "bottleneck_min": int(bottleneck_min),
            "traversable": viable,
            "bottleneck_segment": [seg["name"] for seg in CORRIDOR_SEGMENTS
                                   if segment_viability[seg["name"]][age]["viability_score"] == min_score][0],
        }
        status = "TRAVERSABLE" if viable else "BLOCKED"
        print(f"  {age:>7,d} BP  mean={mean_score:.1f}  min={min_score}  "
              f"bottleneck={bottleneck_min}  [{status}]  "
              f"bottleneck: {corridor_scores[age]['bottleneck_segment']}")

    # --- Persian Gulf analysis ---
    print("\n--- Persian Gulf Dry Periods ---")
    gulf_dry_periods = []
    current_dry = None
    for age in range(0, 131000, 500):
        dry = is_gulf_dry(age)
        if dry and current_dry is None:
            current_dry = age
        elif not dry and current_dry is not None:
            gulf_dry_periods.append((current_dry, age))
            current_dry = None
    if current_dry is not None:
        gulf_dry_periods.append((current_dry, 130000))

    for start, end in gulf_dry_periods:
        sl_start = interpolate_sea_level(start)
        sl_end = interpolate_sea_level(end)
        print(f"  {start:>7,d} – {end:>7,d} BP  "
              f"(SL: {sl_start:+.0f}m to {sl_end:+.0f}m)")

    gulf_results = {
        "dry_periods": [{"start_bp": s, "end_bp": e} for s, e in gulf_dry_periods],
        "threshold_m": GULF_THRESHOLD,
        "note": "Gulf exposed as dry land when sea level below threshold",
    }

    # --- Bottleneck analysis ---
    print("\n--- Bottleneck Analysis ---")
    bottleneck_results = {}
    for seg in CORRIDOR_SEGMENTS:
        if not seg.get("bottleneck"):
            continue
        seg_name = seg["name"]
        # Find epochs where this segment is impassable or marginal
        difficult_epochs = [(age, label) for age, label in KEY_EPOCHS
                           if segment_viability[seg_name][age]["viability_score"] <= 1]
        bottleneck_results[seg_name] = {
            "viability_notes": seg["viability_notes"],
            "difficult_epochs": [{"age_bp": a, "label": l,
                                  "score": segment_viability[seg_name][a]["viability_score"]}
                                 for a, l in difficult_epochs],
            "n_difficult": len(difficult_epochs),
            "n_total": len(KEY_EPOCHS),
        }
        print(f"\n  {seg_name}:")
        print(f"    {seg['viability_notes']}")
        print(f"    Difficult/marginal in {len(difficult_epochs)}/{len(KEY_EPOCHS)} epochs:")
        for a, l in difficult_epochs:
            score = segment_viability[seg_name][a]["viability_score"]
            print(f"      {a:>7,d} BP ({l}) — score {score}")

    # --- Save results ---
    results = {
        "sea_level_by_epoch": sl_results,
        "segment_viability": {
            seg["name"]: {
                str(age): segment_viability[seg["name"]][age]
                for age, _ in KEY_EPOCHS
            }
            for seg in CORRIDOR_SEGMENTS
        },
        "corridor_viability_by_epoch": {
            str(age): corridor_scores[age] for age in corridor_scores
        },
        "persian_gulf_dry_periods": gulf_results,
        "bottleneck_analysis": bottleneck_results,
        "green_arabia_phases": [
            {"start_bp": s, "end_bp": e, "label": l}
            for s, e, l in GREEN_ARABIA_PHASES
        ],
    }

    with open(os.path.join(OUT_DIR, "corridor_viability_by_epoch.json"), "w") as f:
        json.dump(results["corridor_viability_by_epoch"], f, indent=2)

    with open(os.path.join(OUT_DIR, "bottleneck_analysis.json"), "w") as f:
        json.dump(results["bottleneck_analysis"], f, indent=2)

    with open(os.path.join(OUT_DIR, "persian_gulf_dry_periods.json"), "w") as f:
        json.dump(results["persian_gulf_dry_periods"], f, indent=2)

    with open(os.path.join(OUT_DIR, "paleoclimate_corridor.json"), "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to {OUT_DIR}/")

    # --- Multi-panel figure ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(3, 1, figsize=(16, 14), sharex=True)

        ages = np.arange(0, 131000, 500)

        # Panel 1: Sea level curve
        ax1 = axes[0]
        sls = [interpolate_sea_level(a) for a in ages]
        ax1.plot(ages, sls, 'b-', linewidth=1.5)
        ax1.fill_between(ages, sls, 0, where=[s < 0 for s in sls], alpha=0.1, color='blue')
        ax1.axhline(0, color='gray', linestyle='-', alpha=0.3)
        ax1.axhline(GULF_THRESHOLD, color='red', linestyle='--', alpha=0.5,
                     label=f'Gulf exposure threshold ({GULF_THRESHOLD}m)')
        ax1.set_ylabel("Sea Level (m rel. to present)")
        ax1.set_title("Global Sea Level (after Lambeck et al. 2014)")
        ax1.legend(loc='lower right')
        ax1.grid(True, alpha=0.3)

        # Shade Green Arabia
        for s, e, l in GREEN_ARABIA_PHASES:
            ax1.axvspan(e, s, alpha=0.1, color='green')

        # Panel 2: Corridor viability
        ax2 = axes[1]
        epoch_ages = sorted(corridor_scores.keys())
        means = [corridor_scores[a]["mean_viability"] for a in epoch_ages]
        mins = [corridor_scores[a]["min_viability"] for a in epoch_ages]
        bn = [corridor_scores[a]["bottleneck_min"] for a in epoch_ages]

        ax2.plot(epoch_ages, means, 'go-', linewidth=2, markersize=8, label='Mean viability')
        ax2.plot(epoch_ages, mins, 'rs-', linewidth=1.5, markersize=6, label='Min (bottleneck)')
        ax2.axhline(1, color='orange', linestyle='--', alpha=0.5, label='Marginal threshold')
        ax2.set_ylabel("Viability Score (0-3)")
        ax2.set_title("Corridor Viability Over Time")
        ax2.legend(loc='lower right')
        ax2.set_ylim(-0.2, 3.5)
        ax2.grid(True, alpha=0.3)

        # Shade Green Arabia
        for s, e, l in GREEN_ARABIA_PHASES:
            ax2.axvspan(e, s, alpha=0.1, color='green')

        # Panel 3: Segment-level heatmap
        ax3 = axes[2]
        seg_names = [seg["name"] for seg in CORRIDOR_SEGMENTS]
        heatmap_data = np.zeros((len(seg_names), len(KEY_EPOCHS)))
        for i, seg in enumerate(CORRIDOR_SEGMENTS):
            for j, (age, _) in enumerate(KEY_EPOCHS):
                heatmap_data[i, j] = segment_viability[seg["name"]][age]["viability_score"]

        im = ax3.imshow(heatmap_data, aspect='auto', cmap='RdYlGn', vmin=0, vmax=3)
        ax3.set_yticks(range(len(seg_names)))
        ax3.set_yticklabels([s[:35] for s in seg_names], fontsize=9)
        ax3.set_xticks(range(len(KEY_EPOCHS)))
        ax3.set_xticklabels([f"{a//1000}k" for a, _ in KEY_EPOCHS], fontsize=9, rotation=45)
        ax3.set_title("Segment Viability Matrix (green=good, red=blocked)")
        fig.colorbar(im, ax=ax3, shrink=0.5, label="Viability (0-3)")

        fig.suptitle("Paleoclimate Corridor Persistence: 130,000 BP to Present", fontsize=16, y=1.01)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, "paleoclimate_corridor_timelapse.png"),
                    dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Figure saved to {OUT_DIR}/paleoclimate_corridor_timelapse.png")
    except ImportError:
        print("matplotlib not available — skipping figure generation")

    return results


if __name__ == "__main__":
    run_analysis()
