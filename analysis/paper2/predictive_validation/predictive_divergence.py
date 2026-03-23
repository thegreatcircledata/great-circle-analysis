#!/usr/bin/env python3
"""
Predictive Validation: Monument-Settlement Divergence on Post-2001 Data
========================================================================
Using ONLY Pleiades sites added after 2011 (post-Alison's 2001 proposal),
test whether the monument-settlement divergence D = Z_monument - Z_settlement
is positive. A positive D in post-2001 data means the divergence pattern
itself has predictive power, not just the overall alignment.

Split criterion: Pleiades 'created' field.
  - 2009-2011 = Barrington Atlas batch (pre-2001 knowledge)
  - 2012+     = post-Atlas additions (post-2001 proxy)

Also run the full dataset as reference to compare.
"""

import csv, math, os, sys, json, time
import numpy as np

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

# ============================================================
# CONFIGURATION
# ============================================================
POLE_LAT = 59.682122
POLE_LNG = -138.646087
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * np.pi / 2
THRESHOLD_KM = 50
N_MC = 1000

BASE_DIR = os.path.expanduser("~/megalith_site_research")

MONUMENTAL_TYPES = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus", "shrine"
}
SETTLEMENT_TYPES = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement", "townhouse", "production"
}

np.random.seed(42)

# ============================================================
# VECTORIZED HELPERS
# ============================================================
def dist_from_gc_vec(lats, lons):
    lat1r = np.radians(POLE_LAT)
    lon1r = np.radians(POLE_LNG)
    lat2r = np.radians(lats)
    lon2r = np.radians(lons)
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1r) * np.cos(lat2r) * np.sin(dlon / 2) ** 2
    d = EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))
    return np.abs(d - QUARTER_CIRC)


def mc_zscore(lats, lons, n_trials=N_MC):
    """Compute Z-score with distribution-matched MC baseline."""
    n = len(lats)
    if n == 0:
        return {"z_score": 0, "observed": 0, "random_mean": 0, "random_std": 0,
                "enrichment": 0, "n_sites": 0}

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


def list_near_circle(lats, lons, titles, types_list, max_km=50):
    """Return list of sites within threshold, sorted by distance."""
    dists = dist_from_gc_vec(lats, lons)
    mask = dists <= max_km
    near = []
    for i in np.where(mask)[0]:
        near.append({
            "title": titles[i],
            "types": types_list[i],
            "lat": round(float(lats[i]), 4),
            "lon": round(float(lons[i]), 4),
            "distance_km": round(float(dists[i]), 1)
        })
    near.sort(key=lambda x: x["distance_km"])
    return near


# ============================================================
# LOAD AND SPLIT PLEIADES DATA
# ============================================================
print("=" * 70)
print("LOADING PLEIADES DATA — SPLIT BY CREATION DATE")
print("=" * 70)

pleiades_path = os.path.join(BASE_DIR, "data/pleiades/pleiades-places-latest.csv")

# Storage: {group: {category: [(lat, lon, title, types_str), ...]}}
groups = {
    "pre_2001": {"monument_all": [], "settlement_all": [],
                 "monument_ancient": [], "settlement_ancient": [],
                 "unclassified": []},
    "post_2001": {"monument_all": [], "settlement_all": [],
                  "monument_ancient": [], "settlement_ancient": [],
                  "unclassified": []},
    "full": {"monument_all": [], "settlement_all": [],
             "monument_ancient": [], "settlement_ancient": [],
             "unclassified": []}
}

with open(pleiades_path, encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            lat = float(row["reprLat"])
            lon = float(row["reprLong"])
            if lat == 0 and lon == 0:
                continue
        except (ValueError, TypeError):
            continue

        created = row.get("created", "")
        created_year = int(created[:4]) if created and len(created) >= 4 else 9999
        group = "pre_2001" if created_year <= 2011 else "post_2001"

        feature_types = {t.strip() for t in row.get("featureTypes", "").split(",")}
        is_monumental = bool(feature_types & MONUMENTAL_TYPES)
        is_settlement = bool(feature_types & SETTLEMENT_TYPES)

        try:
            min_date = int(row.get("minDate", "0"))
        except (ValueError, TypeError):
            min_date = 0

        title = row.get("title", "")
        types_str = row.get("featureTypes", "")
        entry = (lat, lon, title, types_str)

        if is_monumental:
            groups[group]["monument_all"].append(entry)
            groups["full"]["monument_all"].append(entry)
            if min_date < -2000:
                groups[group]["monument_ancient"].append(entry)
                groups["full"]["monument_ancient"].append(entry)
        elif is_settlement:
            groups[group]["settlement_all"].append(entry)
            groups["full"]["settlement_all"].append(entry)
            if min_date < -2000:
                groups[group]["settlement_ancient"].append(entry)
                groups["full"]["settlement_ancient"].append(entry)
        else:
            groups[group]["unclassified"].append(entry)
            groups["full"]["unclassified"].append(entry)

# Print counts
for gname in ["pre_2001", "post_2001", "full"]:
    g = groups[gname]
    total = sum(len(v) for v in g.values())
    print(f"\n{gname.upper()} ({total} total):")
    for cat in ["monument_all", "settlement_all", "monument_ancient", "settlement_ancient", "unclassified"]:
        print(f"  {cat:>25}: {len(g[cat]):>6}")


# ============================================================
# COMPUTE DIVERGENCE FOR EACH GROUP
# ============================================================
def compute_divergence(group_data, group_name, n_mc=N_MC):
    """Compute monument-settlement divergence D for a group."""
    print(f"\n{'=' * 70}")
    print(f"DIVERGENCE TEST: {group_name}")
    print(f"{'=' * 70}")

    results = {}
    for period_label, mon_key, set_key in [
        ("all_periods", "monument_all", "settlement_all"),
        ("ancient_pre2000BCE", "monument_ancient", "settlement_ancient")
    ]:
        mon_data = group_data[mon_key]
        set_data = group_data[set_key]

        if len(mon_data) == 0 or len(set_data) == 0:
            print(f"\n  {period_label}: SKIPPED (monument={len(mon_data)}, settlement={len(set_data)})")
            results[period_label] = {"skipped": True, "reason": "insufficient data"}
            continue

        mon_lats = np.array([x[0] for x in mon_data])
        mon_lons = np.array([x[1] for x in mon_data])
        mon_titles = [x[2] for x in mon_data]
        mon_types = [x[3] for x in mon_data]

        set_lats = np.array([x[0] for x in set_data])
        set_lons = np.array([x[1] for x in set_data])
        set_titles = [x[2] for x in set_data]
        set_types = [x[3] for x in set_data]

        print(f"\n  {period_label}:")
        print(f"    Monuments: {len(mon_data)}")
        print(f"    Settlements: {len(set_data)}")

        t0 = time.time()
        mon_z = mc_zscore(mon_lats, mon_lons, n_mc)
        set_z = mc_zscore(set_lats, set_lons, n_mc)
        elapsed = time.time() - t0

        D = round(mon_z["z_score"] - set_z["z_score"], 2)

        print(f"    Monument Z:   {mon_z['z_score']:+.2f} (obs={mon_z['observed']}, "
              f"mean={mon_z['random_mean']:.1f}, {mon_z['enrichment']:.2f}x)")
        print(f"    Settlement Z: {set_z['z_score']:+.2f} (obs={set_z['observed']}, "
              f"mean={set_z['random_mean']:.1f}, {set_z['enrichment']:.2f}x)")
        print(f"    D = {D:+.2f}  [{elapsed:.1f}s]")

        if D > 0:
            print(f"    ** POSITIVE DIVERGENCE: monuments cluster more than settlements **")
        else:
            print(f"    ** NEGATIVE/ZERO DIVERGENCE: no monument preference **")

        # List near-circle monuments for post-2001
        mon_near = list_near_circle(mon_lats, mon_lons, mon_titles, mon_types)
        set_near = list_near_circle(set_lats, set_lons, set_titles, set_types)

        results[period_label] = {
            "monument": mon_z,
            "settlement": set_z,
            "D": D,
            "monuments_within_50km": mon_near[:20],
            "settlements_within_50km": set_near[:20],
            "n_monuments_near": len(mon_near),
            "n_settlements_near": len(set_near)
        }

    return results


# Run for all three groups
all_results = {}
for gname in ["post_2001", "pre_2001", "full"]:
    label = {
        "post_2001": "POST-2001 (Pleiades created 2012+)",
        "pre_2001": "PRE-2001 (Pleiades created 2009-2011, Barrington Atlas)",
        "full": "FULL DATASET (all Pleiades)"
    }[gname]
    all_results[gname] = compute_divergence(groups[gname], label)


# ============================================================
# SUMMARY TABLE
# ============================================================
print(f"\n{'=' * 70}")
print("SUMMARY: MONUMENT-SETTLEMENT DIVERGENCE")
print(f"{'=' * 70}")

print(f"\n{'Group':<35} | {'Period':<20} | {'Mon Z':>7} | {'Set Z':>7} | {'D':>7} | {'Verdict'}")
print("-" * 100)
for gname in ["pre_2001", "post_2001", "full"]:
    for period in ["all_periods", "ancient_pre2000BCE"]:
        r = all_results[gname].get(period, {})
        if r.get("skipped"):
            print(f"  {gname:<33} | {period:<20} | {'—':>7} | {'—':>7} | {'—':>7} | SKIPPED")
            continue
        mon_z = r["monument"]["z_score"]
        set_z = r["settlement"]["z_score"]
        D = r["D"]
        verdict = "PREDICTIVE" if D > 2 else "WEAK" if D > 0 else "NO DIVERGENCE"
        print(f"  {gname:<33} | {period:<20} | {mon_z:>+7.2f} | {set_z:>+7.2f} | {D:>+7.2f} | {verdict}")

# Reference values from original analysis
print(f"\n  Reference (original Alison analysis):")
print(f"  {'all Pleiades (original)':<33} | {'ancient_pre2000BCE':<20} | {'—':>7} | {'—':>7} | {'+12.78':>7} | PREDICTIVE")
print(f"  {'all Pleiades (original)':<33} | {'all_periods':<20} | {'—':>7} | {'—':>7} | {'+7.36':>7} | PREDICTIVE")


# ============================================================
# INTERPRETATION
# ============================================================
print(f"\n{'=' * 70}")
print("INTERPRETATION")
print(f"{'=' * 70}")

post_all = all_results["post_2001"].get("all_periods", {})
post_anc = all_results["post_2001"].get("ancient_pre2000BCE", {})

if not post_all.get("skipped") and post_all["D"] > 0:
    print(f"\nAll-periods post-2001: D = {post_all['D']:+.2f}")
    print(f"  Monuments added after Alison's proposal STILL cluster more than settlements.")
    if post_all["D"] > 2:
        print(f"  The divergence is substantial (D > 2). PREDICTIVE POWER CONFIRMED.")
    else:
        print(f"  The divergence is positive but modest (D ≤ 2).")

if not post_anc.get("skipped"):
    if post_anc["D"] > 0:
        print(f"\nAncient post-2001: D = {post_anc['D']:+.2f}")
        print(f"  Even for ancient monuments discovered/digitized after 2001, the pattern holds.")
    else:
        print(f"\nAncient post-2001: D = {post_anc['D']:+.2f}")
        print(f"  Ancient subset shows no divergence in post-2001 data.")
elif post_anc.get("skipped"):
    print(f"\nAncient post-2001: SKIPPED due to insufficient data in one category.")

print(f"\nCaveat: Pleiades 'created' dates are database entry dates (all ≥2009),")
print(f"not discovery dates. The 2009-2011 batch = Barrington Atlas (pub. 2000).")
print(f"Post-2011 = later additions, many representing newly studied sites.")


# ============================================================
# SAVE RESULTS
# ============================================================
output = {
    "meta": {
        "name": "Predictive Validation — Monument-Settlement Divergence",
        "description": "Tests whether D = Z_monument - Z_settlement > 0 in post-2001 Pleiades data",
        "methodology": {
            "divergence_D": "Z_monument - Z_settlement at 50km threshold",
            "baseline": "Distribution-matched MC (±2° Gaussian jitter, independent lat/lon shuffle)",
            "n_mc_trials": N_MC,
            "threshold_km": THRESHOLD_KM,
            "pleiades_split": "created ≤2011 (Barrington Atlas batch) vs created ≥2012",
            "monument_types": sorted(MONUMENTAL_TYPES),
            "settlement_types": sorted(SETTLEMENT_TYPES),
            "ancient_cutoff": "minDate < -2000 BCE"
        },
        "reference_values": {
            "original_D_ancient": 12.78,
            "original_D_all": 7.36
        }
    },
    "pre_2001": {
        "n_total": sum(len(v) for v in groups["pre_2001"].values()),
        "results": all_results["pre_2001"]
    },
    "post_2001": {
        "n_total": sum(len(v) for v in groups["post_2001"].values()),
        "results": all_results["post_2001"]
    },
    "full": {
        "n_total": sum(len(v) for v in groups["full"].values()),
        "results": all_results["full"]
    }
}

# Clean numpy types for JSON serialization
def clean_for_json(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, dict):
        return {k: clean_for_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [clean_for_json(v) for v in obj]
    return obj

outpath = os.path.join(BASE_DIR, "results/predictive_validation_divergence.json")
with open(outpath, "w") as f:
    json.dump(clean_for_json(output), f, indent=2)
print(f"\nResults saved to {outpath}")
