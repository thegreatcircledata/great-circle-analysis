#!/usr/bin/env python3
"""
Great Circle Proximity Test — Standalone CLI Tool
===================================================
Tests whether archaeological sites cluster near a proposed great circle
more than expected by chance, using distribution-matched Monte Carlo.

Usage:
    python great_circle_test.py --input sites.csv --output results.json
    python great_circle_test.py --input sites.csv --pole-lat 59.68 --pole-lon -138.65 --trials 200
    python great_circle_test.py --input sites.csv --thresholds 25,50,100,200

Input CSV must have columns: lat, lon (and optionally: name, type)
Output is a JSON file with Z-scores at each threshold.

Requires only Python 3.7+ stdlib.
"""

import math, random, csv, json, time, argparse, os, sys

# ============================================================
# CONSTANTS
# ============================================================
EARTH_R_KM = 6371.0

# ============================================================
# MATH
# ============================================================
def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))

def gc_distance(lat, lon, pole_lat, pole_lon, quarter_circ):
    return abs(haversine_km(lat, lon, pole_lat, pole_lon) - quarter_circ)

def rand_matched(site_lats, site_lons):
    lat = random.choice(site_lats) + random.gauss(0, 2)
    lon = random.choice(site_lons) + random.gauss(0, 2)
    return max(-90, min(90, lat)), max(-180, min(180, lon))

# ============================================================
# MAIN
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description="Test archaeological site proximity to a great circle")
    parser.add_argument("--input", required=True, help="Input CSV (must have lat, lon columns)")
    parser.add_argument("--output", default="great_circle_results.json", help="Output JSON path")
    parser.add_argument("--pole-lat", type=float, default=59.682122,
                        help="Great circle pole latitude (default: 59.682122)")
    parser.add_argument("--pole-lon", type=float, default=-138.646087,
                        help="Great circle pole longitude (default: -138.646087)")
    parser.add_argument("--trials", type=int, default=200, help="Monte Carlo trials (default: 200)")
    parser.add_argument("--thresholds", default="25,50,100,200",
                        help="Comma-separated distance thresholds in km")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    args = parser.parse_args()

    random.seed(args.seed)
    thresholds = [int(t) for t in args.thresholds.split(",")]
    quarter_circ = EARTH_R_KM * math.pi / 2

    # Load sites
    sites = []
    with open(args.input, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["lat"])
                lon = float(row["lon"])
                if -90 <= lat <= 90 and -180 <= lon <= 180:
                    sites.append({
                        "lat": lat, "lon": lon,
                        "name": row.get("name", ""),
                        "type": row.get("type", "")
                    })
            except (ValueError, KeyError):
                continue

    n = len(sites)
    print(f"Loaded {n} sites from {args.input}")
    print(f"Pole: ({args.pole_lat}, {args.pole_lon})")
    print(f"Quarter circumference: {quarter_circ:.1f} km")
    print(f"Thresholds: {thresholds} km")
    print(f"Trials: {args.trials}")
    print()

    # Compute distances
    for s in sites:
        s["gc_dist"] = gc_distance(s["lat"], s["lon"], args.pole_lat, args.pole_lon, quarter_circ)

    site_lats = [s["lat"] for s in sites]
    site_lons = [s["lon"] for s in sites]

    # Observed counts
    observed = {}
    for t in thresholds:
        obs = sum(1 for s in sites if s["gc_dist"] <= t)
        observed[t] = obs
        print(f"  Within {t:>4}km: {obs} sites")

    # Monte Carlo
    print(f"\nRunning {args.trials} distribution-matched Monte Carlo trials...")
    baseline = {t: [] for t in thresholds}
    t_start = time.time()

    for trial in range(args.trials):
        pts = [rand_matched(site_lats, site_lons) for _ in range(n)]
        for t in thresholds:
            baseline[t].append(
                sum(1 for lat, lon in pts
                    if gc_distance(lat, lon, args.pole_lat, args.pole_lon, quarter_circ) <= t))
        if (trial + 1) % 50 == 0:
            elapsed = time.time() - t_start
            rate = (trial + 1) / elapsed
            eta = (args.trials - trial - 1) / rate
            print(f"  {trial+1}/{args.trials} ({elapsed:.0f}s, ETA {eta:.0f}s)")

    # Results
    print(f"\n{'Threshold':>10} | {'Observed':>8} | {'Mean':>8} | {'Std':>8} | {'Enrich':>7} | {'Z-score':>8} | {'p-value':>8}")
    print("-" * 75)

    results = {"pole": {"lat": args.pole_lat, "lon": args.pole_lon},
               "n_sites": n, "trials": args.trials, "seed": args.seed,
               "thresholds": {}}

    for t in thresholds:
        obs = observed[t]
        mu = sum(baseline[t]) / args.trials
        sigma = (sum((x - mu)**2 for x in baseline[t]) / args.trials) ** 0.5
        z = (obs - mu) / sigma if sigma > 0 else 0
        enrich = obs / mu if mu > 0 else float('inf')
        p = sum(1 for b in baseline[t] if b >= obs) / args.trials

        results["thresholds"][str(t)] = {
            "observed": obs, "baseline_mean": round(mu, 1),
            "baseline_std": round(sigma, 1), "z_score": round(z, 2),
            "enrichment": round(enrich, 2), "p_value": p
        }
        print(f"  {t:>6}km | {obs:>8} | {mu:>8.1f} | {sigma:>8.1f} | {enrich:>6.1f}x | {z:>8.2f} | {p:>8.4f}")

    # Top sites near circle
    nearby = sorted([s for s in sites if s["gc_dist"] <= 50], key=lambda x: x["gc_dist"])
    results["nearest_sites"] = [
        {"name": s["name"], "type": s["type"], "lat": s["lat"], "lon": s["lon"],
         "distance_km": round(s["gc_dist"], 1)}
        for s in nearby[:50]
    ]

    elapsed = time.time() - t_start
    results["elapsed_seconds"] = round(elapsed, 1)

    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {args.output} ({elapsed:.0f}s)")

if __name__ == "__main__":
    main()
