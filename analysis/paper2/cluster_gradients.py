#!/usr/bin/env python3
"""
Test 5: Cluster Gradients — "earlier = closer" across monument clusters.

For each cluster region along the great circle, tests whether older
radiocarbon-dated sites tend to be closer to the GC (negative correlation
between age and distance from GC).

Uses p3k14c radiocarbon database, deduplicated to unique sites.
"""

import csv
import json
import math
import os
import numpy as np
from scipy import stats

# --- Great Circle parameters ---
POLE_LAT = math.radians(59.682122)
POLE_LON = math.radians(-138.646087)
R_EARTH = 6371.0
QUARTER_CIRC = R_EARTH * math.pi / 2  # ~10,007.5 km


def haversine(lat1, lon1, lat2, lon2):
    """Great circle distance in km between two points (radians)."""
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 2 * R_EARTH * math.asin(math.sqrt(a))


def gc_distance(lat_deg, lon_deg):
    """Signed distance from GC in km (negative = south/inside)."""
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    d_to_pole = haversine(lat, lon, POLE_LAT, POLE_LON)
    return d_to_pole - QUARTER_CIRC


# --- Cluster definitions ---
# Each cluster: name, lat_range, lon_range, description
# Generous bounding boxes to capture sites near the GC in each region
CLUSTERS = {
    "Memphis_Egypt": {
        "lat_min": 25.0, "lat_max": 32.0,
        "lon_min": 29.0, "lon_max": 35.0,
        "description": "Nile Valley pyramids and temples (Egypt)"
    },
    "Nazca_Peru": {
        "lat_min": -17.0, "lat_max": -13.0,
        "lon_min": -77.0, "lon_max": -73.0,
        "description": "Nazca/Ica region geoglyphs and sites (Peru)"
    },
    "Easter_Island": {
        "lat_min": -28.0, "lat_max": -26.5,
        "lon_min": -110.0, "lon_max": -108.5,
        "description": "Rapa Nui ahu and settlements"
    },
    "Persepolis_Iran": {
        "lat_min": 27.0, "lat_max": 34.0,
        "lon_min": 49.0, "lon_max": 56.0,
        "description": "Fars/Isfahan region sites (Iran)"
    },
    "Levant_PPNB": {
        "lat_min": 29.0, "lat_max": 37.0,
        "lon_min": 34.0, "lon_max": 42.0,
        "description": "Levant Pre-Pottery Neolithic and later (Israel/Jordan/Syria/Lebanon)"
    },
    "Indus_Valley": {
        "lat_min": 24.0, "lat_max": 32.0,
        "lon_min": 66.0, "lon_max": 74.0,
        "description": "Indus Valley / Sindh-Punjab sites (Pakistan/India)"
    },
    "Anatolia": {
        "lat_min": 36.0, "lat_max": 42.0,
        "lon_min": 26.0, "lon_max": 44.0,
        "description": "Anatolian Neolithic and Bronze Age sites (Turkey)"
    },
}

# Maximum distance from GC to include (km) — generous to get enough sites
MAX_GC_DIST = 500  # km


def load_p3k14c():
    """Load p3k14c, deduplicate to unique sites (oldest age, mean coords).

    Sites WITH a SiteID are grouped by SiteID.
    Sites WITHOUT a SiteID are grouped by rounded coordinates (2 decimal places).
    """
    path = os.path.join(os.path.dirname(__file__), "..", "data", "p3k14c", "p3k14c_data.csv")
    sites = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["Lat"])
                lon = float(row["Long"])
                age = int(row["Age"])  # radiocarbon years BP
                sid = row.get("SiteID", "").strip()
                name = row.get("SiteName", "")
                country = row.get("Country", "")
            except (ValueError, KeyError):
                continue

            # Use SiteID if available, otherwise round coordinates as key
            if sid:
                key = f"sid_{sid}"
            else:
                key = f"coord_{round(lat, 2)}_{round(lon, 2)}"

            if key not in sites:
                sites[key] = {
                    "lats": [lat], "lons": [lon],
                    "oldest_age": age, "name": name, "country": country
                }
            else:
                sites[key]["lats"].append(lat)
                sites[key]["lons"].append(lon)
                if age > sites[key]["oldest_age"]:
                    sites[key]["oldest_age"] = age

    # Deduplicate: mean coords, oldest age
    result = []
    for key, s in sites.items():
        result.append({
            "site_id": key,
            "name": s["name"],
            "country": s["country"],
            "lat": np.mean(s["lats"]),
            "lon": np.mean(s["lons"]),
            "age_bp": s["oldest_age"],
        })
    return result


def filter_cluster(sites, cluster):
    """Filter sites within cluster bounding box AND within MAX_GC_DIST of GC."""
    c = CLUSTERS[cluster]
    filtered = []
    for s in sites:
        if c["lat_min"] <= s["lat"] <= c["lat_max"] and c["lon_min"] <= s["lon"] <= c["lon_max"]:
            d = gc_distance(s["lat"], s["lon"])
            if abs(d) <= MAX_GC_DIST:
                filtered.append({**s, "gc_distance_km": abs(d)})
    return filtered


def test_cluster(sites, cluster_name, n_permutations=10000):
    """
    Test correlation between age (older = larger BP) and gc_distance.
    If earlier = closer, we expect NEGATIVE r: older sites (high BP) have
    smaller gc_distance.
    """
    filtered = filter_cluster(sites, cluster_name)
    n = len(filtered)

    if n < 5:
        return {
            "cluster": cluster_name,
            "description": CLUSTERS[cluster_name]["description"],
            "n_sites": n,
            "status": "SKIPPED — fewer than 5 dated sites",
        }

    ages = np.array([s["age_bp"] for s in filtered])
    dists = np.array([s["gc_distance_km"] for s in filtered])

    # Pearson correlation
    r, p_parametric = stats.pearsonr(ages, dists)

    # Also Spearman (robust to outliers)
    rho, p_spearman = stats.spearmanr(ages, dists)

    # Permutation test: shuffle ages, compute r
    rng = np.random.default_rng(42)
    perm_rs = np.empty(n_permutations)
    for i in range(n_permutations):
        shuffled = rng.permutation(ages)
        perm_rs[i] = np.corrcoef(shuffled, dists)[0, 1]

    # One-tailed: proportion of permutation r <= observed r (testing r < 0)
    p_perm = np.mean(perm_rs <= r)

    # Age and distance statistics
    age_range = [int(ages.min()), int(ages.max())]
    dist_range = [round(float(dists.min()), 1), round(float(dists.max()), 1)]

    # Example sites (5 oldest)
    sorted_sites = sorted(filtered, key=lambda s: -s["age_bp"])
    examples = [
        {"name": s["name"], "age_bp": s["age_bp"],
         "gc_distance_km": round(s["gc_distance_km"], 1)}
        for s in sorted_sites[:5]
    ]

    return {
        "cluster": cluster_name,
        "description": CLUSTERS[cluster_name]["description"],
        "n_sites": n,
        "pearson_r": round(r, 4),
        "pearson_p": round(p_parametric, 6),
        "spearman_rho": round(rho, 4),
        "spearman_p": round(p_spearman, 6),
        "permutation_p": round(p_perm, 4),
        "n_permutations": n_permutations,
        "age_range_bp": age_range,
        "distance_range_km": dist_range,
        "mean_age_bp": int(np.mean(ages)),
        "mean_distance_km": round(float(np.mean(dists)), 1),
        "oldest_5_sites": examples,
        "interpretation": (
            f"r={r:.4f}: {'earlier monuments ARE closer to GC' if r < 0 else 'earlier monuments are NOT closer to GC'} "
            f"(permutation p={p_perm:.4f})"
        ),
    }


def main():
    print("Loading p3k14c database...")
    sites = load_p3k14c()
    print(f"  {len(sites)} unique sites loaded")

    results = {
        "test": "Test 5: Cluster Gradients — earlier = closer?",
        "method": (
            "For each geographic cluster near the GC, compute Pearson r and Spearman rho "
            "between radiocarbon age (BP) and absolute distance from GC (km). "
            "If earlier = closer, expect r < 0 (older sites closer to GC). "
            "Permutation test (10,000 shuffles) for significance."
        ),
        "data_source": "p3k14c radiocarbon database (deduplicated to unique sites)",
        "gc_pole": {"lat": 59.682122, "lon": -138.646087},
        "max_gc_distance_km": MAX_GC_DIST,
        "reference_memphis_result": {
            "note": "Previous Memphis-only pyramid test found r=-0.428, p=0.029 (N=26)",
        },
        "clusters": [],
        "pattern_summary": {},
    }

    negative_r_clusters = []
    significant_clusters = []

    for cluster_name in CLUSTERS:
        print(f"\nTesting {cluster_name}...")
        result = test_cluster(sites, cluster_name)
        results["clusters"].append(result)

        if result.get("pearson_r") is not None:
            print(f"  N={result['n_sites']}, r={result.get('pearson_r', 'N/A')}, "
                  f"perm_p={result.get('permutation_p', 'N/A')}")
            if result["pearson_r"] < 0:
                negative_r_clusters.append(cluster_name)
            if result.get("permutation_p", 1) < 0.05:
                significant_clusters.append(cluster_name)
        else:
            print(f"  {result.get('status', 'N/A')}")

    # Pattern summary
    tested = [c for c in results["clusters"] if "pearson_r" in c]
    results["pattern_summary"] = {
        "clusters_tested": len(tested),
        "clusters_skipped": len(results["clusters"]) - len(tested),
        "clusters_with_negative_r": negative_r_clusters,
        "clusters_significant_p05": significant_clusters,
        "is_pattern": len(negative_r_clusters) >= 2,
        "is_memphis_specific": len(negative_r_clusters) <= 1 and "Memphis_Egypt" in negative_r_clusters,
        "verdict": (
            f"{len(negative_r_clusters)} of {len(tested)} tested clusters show r < 0 (earlier = closer). "
            f"{len(significant_clusters)} significant at p < 0.05. "
            + ("PATTERN detected across multiple clusters."
               if len(negative_r_clusters) >= 2
               else "Effect appears Memphis-specific.")
        ),
    }

    # Write output
    out_dir = os.path.join(os.path.dirname(__file__), "..", "outputs", "resolution_tests")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "cluster_gradients.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n{'='*60}")
    print(f"Results written to {out_path}")
    print(f"\nVERDICT: {results['pattern_summary']['verdict']}")


if __name__ == "__main__":
    main()
