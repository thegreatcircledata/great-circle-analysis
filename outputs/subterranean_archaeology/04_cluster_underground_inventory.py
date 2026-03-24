#!/usr/bin/env python3
"""
Directive 12, Script 04: Underground Structures Beneath Known Circle Sites
==========================================================================
Catalogs known underground features at each of the 6 major Great Circle
clusters. Compares underground feature density at circle sites vs
comparable control sites.

Output:
  - cluster_underground_inventory.json
  - underground_density_comparison.json
  - underground_features_map.png
"""

import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))

# ── Cluster Underground Inventories ───────────────────────────────────
# Each feature classified as:
#   natural, minimally_modified, architecturally_significant, monumental_underground

CLUSTER_INVENTORY = {
    "Egypt (Memphis/Giza)": {
        "center": {"lat": 29.5, "lon": 31.2},
        "features": [
            {"name": "Osiris Shaft beneath Giza causeway",
             "type": "monumental_underground", "depth_m": 30, "levels": 3,
             "notes": "Partially flooded, 3 levels, sarcophagus at bottom"},
            {"name": "Subterranean Chamber of the Great Pyramid",
             "type": "architecturally_significant", "depth_m": 30,
             "notes": "Carved from bedrock 30m below base, unfinished"},
            {"name": "ScanPyramids void (2017 muon tomography)",
             "type": "architecturally_significant", "depth_m": None,
             "notes": "~30m long void above Grand Gallery, discovered by muon tomography"},
            {"name": "Saqqara Serapeum",
             "type": "monumental_underground", "depth_m": 12,
             "notes": "~400m of tunnels, 24 massive granite sarcophagi (60-80 tons each)"},
            {"name": "Saqqara ibis catacombs",
             "type": "monumental_underground", "depth_m": 10,
             "notes": "~4 million mummified ibis, km of tunnels"},
            {"name": "Saqqara baboon catacombs",
             "type": "architecturally_significant", "depth_m": 8,
             "notes": "Adjacent to ibis catacombs"},
            {"name": "Djoser Step Pyramid substructure",
             "type": "monumental_underground", "depth_m": 28,
             "notes": "28m deep shaft, 5.7km of tunnels beneath pyramid"},
            {"name": "Valley of the Kings tomb complex",
             "type": "monumental_underground", "depth_m": 50,
             "notes": "63 tombs, KV5 alone has 121+ chambers (Weeks 1995)"},
            {"name": "Valley of the Queens",
             "type": "monumental_underground", "depth_m": 15,
             "notes": "~80 tombs including Nefertari's tomb"},
            {"name": "Osireion at Abydos",
             "type": "monumental_underground", "depth_m": 15,
             "notes": "Megalithic subterranean temple, granite blocks, water table"},
            {"name": "Giza Plateau trial passages / tunnels",
             "type": "architecturally_significant", "depth_m": 10,
             "notes": "Exploratory rock-cut passages near Great Pyramid"},
            {"name": "Abu Rawash pyramid shaft",
             "type": "architecturally_significant", "depth_m": 21,
             "notes": "Deep shaft beneath Djedefre's pyramid"},
            {"name": "Medium pyramid substructure",
             "type": "architecturally_significant", "depth_m": 15,
             "notes": "Descending passage and burial chamber"},
        ],
    },
    "Peru (Nazca/Cahuachi)": {
        "center": {"lat": -14.8, "lon": -75.1},
        "features": [
            {"name": "Chavín de Huántar underground galleries",
             "type": "monumental_underground", "depth_m": 15,
             "notes": "3+ levels, acoustic properties, Lanzón monolith chamber"},
            {"name": "Nazca puquios (36+ underground aqueducts)",
             "type": "architecturally_significant", "depth_m": 10,
             "notes": "Spiral access points, still functioning, UNESCO tentative list"},
            {"name": "Cahuachi buried construction phases",
             "type": "architecturally_significant", "depth_m": 5,
             "notes": "Multiple building phases buried beneath mound construction"},
            {"name": "Nazca Líneas subterranean water management",
             "type": "architecturally_significant", "depth_m": None,
             "notes": "Underground filtration galleries connected to geoglyphs"},
            {"name": "Wari underground chambers at Pikillacta",
             "type": "architecturally_significant", "depth_m": 5,
             "notes": "Subterranean storage/ritual chambers"},
        ],
    },
    "Iran (Persepolis)": {
        "center": {"lat": 29.9, "lon": 52.9},
        "features": [
            {"name": "Naqsh-e Rostam rock-cut tombs",
             "type": "monumental_underground", "depth_m": 15,
             "notes": "4 Achaemenid royal tombs carved into cliff face, Ka'ba-ye Zartosht"},
            {"name": "Persepolis drainage tunnel system",
             "type": "architecturally_significant", "depth_m": 5,
             "notes": "Extensive underground water management beneath the platform"},
            {"name": "Persepolis Treasury foundations",
             "type": "architecturally_significant", "depth_m": 3,
             "notes": "Sub-floor storage and foundation chambers"},
            {"name": "Buried city beneath Persepolis (magnetometry 2020s)",
             "type": "architecturally_significant", "depth_m": 3,
             "notes": "Magnetometry revealed extensive buried structures beyond visible ruins"},
            {"name": "Naqsh-e Rajab rock reliefs",
             "type": "minimally_modified", "depth_m": 0,
             "notes": "Sasanian rock-cut reliefs in natural rock face"},
        ],
    },
    "Indus Valley (Mohenjo-daro)": {
        "center": {"lat": 27.3, "lon": 68.1},
        "features": [
            {"name": "Great Bath and waterproofing/drainage system",
             "type": "architecturally_significant", "depth_m": 3,
             "notes": "Bitumen-sealed, earliest known public water management structure"},
            {"name": "Underground drainage/sewer network",
             "type": "architecturally_significant", "depth_m": 2,
             "notes": "Covered brick drains beneath streets, one of earliest urban sewers"},
            {"name": "Unexcavated lower waterlogged layers",
             "type": "natural", "depth_m": 3,
             "notes": "Estimated 3+ meters of stratigraphy never excavated due to water table"},
            {"name": "Well network (~700 wells)",
             "type": "architecturally_significant", "depth_m": 10,
             "notes": "Brick-lined wells, highest density of any ancient city"},
        ],
    },
    "Easter Island": {
        "center": {"lat": -27.1, "lon": -109.4},
        "features": [
            {"name": "Ana Kai Tangata (painted cave)",
             "type": "minimally_modified", "depth_m": 5,
             "notes": "Sooty tern (manutara) paintings, ceremony of the Birdman cult"},
            {"name": "Ana Te Pahu (lava tube complex)",
             "type": "natural", "depth_m": 8,
             "notes": "Largest lava tube, agricultural use, banana cultivation inside"},
            {"name": "Ana O Keke (virgin cave)",
             "type": "minimally_modified", "depth_m": 10,
             "notes": "Ritual seclusion cave on Poike cliff face"},
            {"name": "Numerous habitation lava tubes",
             "type": "natural", "depth_m": 5,
             "notes": "Multiple lava tubes used as shelters across the island"},
            {"name": "Ana Heu Neru (meditation cave)",
             "type": "minimally_modified", "depth_m": 3,
             "notes": "Cave at Orongo with petroglyphs"},
        ],
    },
    "Southeast Asia": {
        "center": {"lat": 15.0, "lon": 105.0},
        "features": [
            {"name": "Plain of Jars (underground contexts)",
             "type": "minimally_modified", "depth_m": 2,
             "notes": "Some jars associated with cave burial practices"},
            {"name": "Tam Pa Ling cave (Laos)",
             "type": "natural", "depth_m": 5,
             "notes": "Oldest modern human fossils in mainland SE Asia (~46-63 ka)"},
            {"name": "Tham Lod rock shelter (Thailand)",
             "type": "natural", "depth_m": 3,
             "notes": "Paleolithic habitation and coffin cave burials"},
        ],
    },
}

# Control sites: major archaeological sites NOT on the circle
CONTROL_SITES = {
    "Göbekli Tepe": {
        "center": {"lat": 37.223, "lon": 38.923},
        "gc_distance_km": None,
        "features": [
            {"name": "Enclosures A-H (partially subterranean)", "type": "architecturally_significant"},
            {"name": "Buried enclosures (ground-penetrating radar)", "type": "architecturally_significant"},
        ],
    },
    "Stonehenge": {
        "center": {"lat": 51.179, "lon": -1.826},
        "gc_distance_km": None,
        "features": [
            {"name": "Aubrey holes (shallow pits)", "type": "minimally_modified"},
            {"name": "No significant underground features", "type": "natural"},
        ],
    },
    "Teotihuacan": {
        "center": {"lat": 19.692, "lon": -98.844},
        "gc_distance_km": None,
        "features": [
            {"name": "Tunnel beneath Temple of the Feathered Serpent", "type": "monumental_underground"},
            {"name": "Cave beneath Pyramid of the Sun", "type": "natural"},
            {"name": "Underground chambers with offerings", "type": "architecturally_significant"},
        ],
    },
    "Angkor Wat": {
        "center": {"lat": 13.412, "lon": 103.867},
        "gc_distance_km": None,
        "features": [
            {"name": "Moat and water management", "type": "architecturally_significant"},
            {"name": "No significant underground features", "type": "natural"},
        ],
    },
    "Machu Picchu": {
        "center": {"lat": -13.163, "lon": -72.545},
        "gc_distance_km": None,
        "features": [
            {"name": "Royal Tomb (cave beneath Temple of the Sun)", "type": "minimally_modified"},
            {"name": "Drainage system", "type": "architecturally_significant"},
        ],
    },
    "Tikal": {
        "center": {"lat": 17.222, "lon": -89.624},
        "gc_distance_km": None,
        "features": [
            {"name": "Burial tombs beneath temples", "type": "architecturally_significant"},
            {"name": "No major underground complexes", "type": "natural"},
        ],
    },
}

FEATURE_WEIGHT = {
    "natural": 1,
    "minimally_modified": 2,
    "architecturally_significant": 3,
    "monumental_underground": 5,
}


def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlat, dlon = la2 - la1, lo2 - lo1
    a = math.sin(dlat / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin(dlon / 2) ** 2
    return R * 2 * math.asin(math.sqrt(min(1.0, a)))


def gc_distance(lat, lon):
    POLE_LAT, POLE_LON = 59.682122, -138.646087
    R = 6371.0
    d = haversine_km(POLE_LAT, POLE_LON, lat, lon)
    return abs(d - R * math.pi / 2)


def main():
    print("=" * 70)
    print("UNDERGROUND STRUCTURES BENEATH KNOWN CIRCLE SITES")
    print("=" * 70)

    # ── Cluster inventory ──
    print(f"\n{'Cluster':<30s}  {'Features':>8s}  {'Monumental':>10s}  {'Score':>6s}")
    print(f"{'-' * 30}  {'-' * 8}  {'-' * 10}  {'-' * 6}")

    cluster_summary = {}
    for name, data in CLUSTER_INVENTORY.items():
        features = data["features"]
        n_features = len(features)
        n_monumental = sum(1 for f in features if f["type"] == "monumental_underground")
        score = sum(FEATURE_WEIGHT.get(f["type"], 1) for f in features)
        cluster_summary[name] = {
            "n_features": n_features,
            "n_monumental": n_monumental,
            "weighted_score": score,
            "features": features,
        }
        print(f"  {name:<28s}  {n_features:>8d}  {n_monumental:>10d}  {score:>6d}")

    # ── Control sites ──
    print(f"\n--- Control Sites (NOT on circle) ---")
    print(f"{'Site':<30s}  {'Features':>8s}  {'Monumental':>10s}  {'Score':>6s}  {'GC Dist':>8s}")
    print(f"{'-' * 30}  {'-' * 8}  {'-' * 10}  {'-' * 6}  {'-' * 8}")

    control_summary = {}
    for name, data in CONTROL_SITES.items():
        features = data["features"]
        n_features = len(features)
        n_monumental = sum(1 for f in features if f["type"] == "monumental_underground")
        score = sum(FEATURE_WEIGHT.get(f["type"], 1) for f in features)
        d = gc_distance(data["center"]["lat"], data["center"]["lon"])
        control_summary[name] = {
            "n_features": n_features,
            "n_monumental": n_monumental,
            "weighted_score": score,
            "gc_distance_km": round(d, 1),
        }
        print(f"  {name:<28s}  {n_features:>8d}  {n_monumental:>10d}  {score:>6d}  {d:>7.1f}km")

    # ── Comparison ──
    circle_scores = [v["weighted_score"] for v in cluster_summary.values()]
    control_scores = [v["weighted_score"] for v in control_summary.values()]
    circle_features = [v["n_features"] for v in cluster_summary.values()]
    control_features = [v["n_features"] for v in control_summary.values()]

    import numpy as np
    from scipy import stats

    print(f"\n--- Density Comparison ---")
    print(f"  Circle clusters:  mean features = {np.mean(circle_features):.1f}, "
          f"mean score = {np.mean(circle_scores):.1f}")
    print(f"  Control sites:    mean features = {np.mean(control_features):.1f}, "
          f"mean score = {np.mean(control_scores):.1f}")

    # Mann-Whitney U test
    u_feat, p_feat = stats.mannwhitneyu(circle_features, control_features, alternative="greater")
    u_score, p_score = stats.mannwhitneyu(circle_scores, control_scores, alternative="greater")
    print(f"\n  Mann-Whitney U (features): U={u_feat:.0f}, p={p_feat:.4f} "
          f"{'→ Significant' if p_feat < 0.05 else ''}")
    print(f"  Mann-Whitney U (weighted score): U={u_score:.0f}, p={p_score:.4f} "
          f"{'→ Significant' if p_score < 0.05 else ''}")

    # ── Special note: 5 of 6 clusters have documented subsurface anomalies ──
    print(f"\n--- GPR/Subsurface Anomaly Cross-Reference ---")
    gpr_status = {
        "Egypt (Memphis/Giza)": "Yes — El-khteeb 2025 GPR, ScanPyramids muon tomography",
        "Peru (Nazca/Cahuachi)": "Yes — Cahuachi buried phases, puquio surveys",
        "Iran (Persepolis)": "Yes — 2020s magnetometry beneath Persepolis",
        "Indus Valley (Mohenjo-daro)": "Yes — unexcavated waterlogged layers documented",
        "Easter Island": "Partial — lava tube surveys, no systematic GPR",
        "Southeast Asia": "Limited — Tam Pa Ling excavation context only",
    }
    for cluster, status in gpr_status.items():
        print(f"  {cluster:<30s}: {status}")

    # ── Save results ──
    inventory_output = {
        "analysis": "Underground structures beneath known Great Circle sites",
        "cluster_inventory": {},
        "control_sites": {},
    }
    for name, data in CLUSTER_INVENTORY.items():
        inventory_output["cluster_inventory"][name] = {
            "center": data["center"],
            "n_features": cluster_summary[name]["n_features"],
            "n_monumental": cluster_summary[name]["n_monumental"],
            "weighted_score": cluster_summary[name]["weighted_score"],
            "features": [
                {"name": f["name"], "type": f["type"],
                 "depth_m": f.get("depth_m"), "notes": f.get("notes", "")}
                for f in data["features"]
            ],
        }
    for name, data in CONTROL_SITES.items():
        inventory_output["control_sites"][name] = control_summary[name]

    with open(os.path.join(BASE, "cluster_underground_inventory.json"), "w") as f:
        json.dump(inventory_output, f, indent=2)

    density_output = {
        "analysis": "Underground feature density: circle sites vs control",
        "circle_clusters": {
            "mean_features": round(float(np.mean(circle_features)), 1),
            "mean_weighted_score": round(float(np.mean(circle_scores)), 1),
        },
        "control_sites": {
            "mean_features": round(float(np.mean(control_features)), 1),
            "mean_weighted_score": round(float(np.mean(control_scores)), 1),
        },
        "mann_whitney_features": {
            "U": float(u_feat), "p_value": round(float(p_feat), 4),
            "significant": bool(p_feat < 0.05),
        },
        "mann_whitney_score": {
            "U": float(u_score), "p_value": round(float(p_score), 4),
            "significant": bool(p_score < 0.05),
        },
        "caveat": "Circle clusters include Egypt which has exceptionally rich underground heritage. "
                   "Control selection is not strictly matched. Treat as descriptive, not definitive.",
    }
    with open(os.path.join(BASE, "underground_density_comparison.json"), "w") as f:
        json.dump(density_output, f, indent=2)

    print(f"\nSaved cluster_underground_inventory.json")
    print(f"Saved underground_density_comparison.json")

    # ── Visualization ──
    make_comparison_figure(cluster_summary, control_summary)


def make_comparison_figure(cluster_summary, control_summary):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Feature count comparison
    names_c = list(cluster_summary.keys())
    scores_c = [cluster_summary[n]["n_features"] for n in names_c]
    names_ctrl = list(control_summary.keys())
    scores_ctrl = [control_summary[n]["n_features"] for n in names_ctrl]

    short_names_c = [n.split("(")[0].strip()[:15] for n in names_c]
    short_names_ctrl = [n[:15] for n in names_ctrl]

    ax1.barh(range(len(names_c)), scores_c, color="#F44336", alpha=0.7, label="Circle clusters")
    ax1.barh(range(len(names_c), len(names_c) + len(names_ctrl)), scores_ctrl,
             color="#2196F3", alpha=0.7, label="Control sites")
    ax1.set_yticks(range(len(names_c) + len(names_ctrl)))
    ax1.set_yticklabels(short_names_c + short_names_ctrl, fontsize=7)
    ax1.set_xlabel("Number of underground features")
    ax1.set_title("Underground Feature Count")
    ax1.legend(fontsize=8)

    # Weighted score comparison
    scores_c_w = [cluster_summary[n]["weighted_score"] for n in names_c]
    scores_ctrl_w = [control_summary[n]["weighted_score"] for n in names_ctrl]

    ax2.barh(range(len(names_c)), scores_c_w, color="#F44336", alpha=0.7, label="Circle clusters")
    ax2.barh(range(len(names_c), len(names_c) + len(names_ctrl)), scores_ctrl_w,
             color="#2196F3", alpha=0.7, label="Control sites")
    ax2.set_yticks(range(len(names_c) + len(names_ctrl)))
    ax2.set_yticklabels(short_names_c + short_names_ctrl, fontsize=7)
    ax2.set_xlabel("Weighted underground score")
    ax2.set_title("Weighted Underground Complexity Score")
    ax2.legend(fontsize=8)

    fig.suptitle("Underground Features: Circle Sites vs Control Sites", fontsize=13)
    fig.tight_layout()
    path = os.path.join(BASE, "underground_features_map.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
