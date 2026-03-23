#!/usr/bin/env python3
"""
Predictive Validation Test
===========================
Split databases into pre-2001 and post-2001 sites, compute great circle
Z-scores for each group separately. If the signal holds for post-2001 sites,
the circle has genuine predictive power beyond selection bias.

Databases:
1. Pleiades — split by 'created' field. 2009-2011 batch = Barrington Atlas
   digitization (sites known well before 2001). 2012+ = later additions.
2. p3k14c — split by publication year extracted from Reference field.
3. Megalithic Portal — KML files (no creation dates, used as whole-dataset reference).

Jim Alison proposed the great circle c. 2001.
"""

import csv
import json
import math
import os
import random
import re
import xml.etree.ElementTree as ET
import glob
from collections import Counter

# ============================================================
# CONFIGURATION
# ============================================================
POLE_LAT = 59.682122
POLE_LNG = -138.646087
EARTH_R_KM = 6371
QUARTER_CIRC = EARTH_R_KM * math.pi / 2  # ~10,008 km
THRESHOLDS = [25, 50, 100]
N_TRIALS = 1000
BASE = "/Users/elliotallan/megalith_site_research"

random.seed(42)

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def haversine_km(lat1, lon1, lat2, lon2):
    R = EARTH_R_KM
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return R * 2 * math.asin(math.sqrt(min(1.0, a)))


def dist_from_great_circle(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LNG)
    return abs(d - QUARTER_CIRC)


def compute_z_scores(sites, n_trials=N_TRIALS, thresholds=THRESHOLDS):
    """Compute Z-scores using distribution-matched random baseline."""
    n = len(sites)
    if n == 0:
        return {t: {"observed": 0, "random_mean": 0, "random_std": 0,
                     "z_score": 0, "enrichment": 0, "n_sites": 0} for t in thresholds}

    # Count observed
    distances = [dist_from_great_circle(s["lat"], s["lon"]) for s in sites]
    observed = {t: sum(1 for d in distances if d <= t) for t in thresholds}

    # Distribution-matched baseline
    lats = [s["lat"] for s in sites]
    lons = [s["lon"] for s in sites]

    rand_counts = {t: [] for t in thresholds}
    for trial in range(n_trials):
        # Shuffle lats and lons independently + Gaussian jitter
        trial_lats = [random.choice(lats) + random.gauss(0, 2) for _ in range(n)]
        trial_lons = [random.choice(lons) + random.gauss(0, 2) for _ in range(n)]
        trial_lats = [max(-90, min(90, la)) for la in trial_lats]
        trial_lons = [max(-180, min(180, lo)) for lo in trial_lons]
        for t in thresholds:
            cnt = sum(1 for la, lo in zip(trial_lats, trial_lons)
                      if dist_from_great_circle(la, lo) <= t)
            rand_counts[t].append(cnt)

    results = {}
    for t in thresholds:
        obs = observed[t]
        mu = sum(rand_counts[t]) / n_trials
        sigma = (sum((x - mu) ** 2 for x in rand_counts[t]) / n_trials) ** 0.5
        z = (obs - mu) / sigma if sigma > 0 else 0
        enrich = obs / mu if mu > 0 else 0
        pct_near = obs / n * 100 if n > 0 else 0
        results[t] = {
            "observed": obs,
            "random_mean": round(mu, 2),
            "random_std": round(sigma, 2),
            "z_score": round(z, 2),
            "enrichment": round(enrich, 2),
            "pct_near": round(pct_near, 3),
            "n_sites": n
        }
    return results


# ============================================================
# 1. PLEIADES DATABASE
# ============================================================
print("=" * 70)
print("DATABASE 1: PLEIADES")
print("=" * 70)

pleiades_pre = []  # created 2009-2011 (Barrington Atlas batch, all pre-2001 knowledge)
pleiades_post = []  # created 2012+ (later additions, many post-2001 discoveries)

with open(os.path.join(BASE, "data/pleiades/pleiades-places-latest.csv")) as f:
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

        site = {"lat": lat, "lon": lon, "title": row.get("title", ""),
                "types": row.get("featureTypes", ""), "created": created[:10]}

        if created_year <= 2011:
            pleiades_pre.append(site)
        else:
            pleiades_post.append(site)

print(f"Pre-2001 knowledge (created 2009-2011, Barrington Atlas batch): {len(pleiades_pre)}")
print(f"Post-2001 additions (created 2012+): {len(pleiades_post)}")

print("\nComputing Z-scores for PRE-2001 Pleiades sites...")
pleiades_pre_z = compute_z_scores(pleiades_pre)
print("Computing Z-scores for POST-2001 Pleiades sites...")
pleiades_post_z = compute_z_scores(pleiades_post)

print(f"\n{'':>8} | {'PRE-2001 (n=' + str(len(pleiades_pre)) + ')':>30} | {'POST-2001 (n=' + str(len(pleiades_post)) + ')':>30}")
print(f"{'Thresh':>8} | {'Obs':>5} {'Mean':>7} {'Z':>7} {'Enrich':>7} | {'Obs':>5} {'Mean':>7} {'Z':>7} {'Enrich':>7}")
print("-" * 80)
for t in THRESHOLDS:
    pre = pleiades_pre_z[t]
    post = pleiades_post_z[t]
    print(f"  {t:>4}km | {pre['observed']:>5} {pre['random_mean']:>7.1f} {pre['z_score']:>7.2f} {pre['enrichment']:>6.2f}x"
          f" | {post['observed']:>5} {post['random_mean']:>7.1f} {post['z_score']:>7.2f} {post['enrichment']:>6.2f}x")


# ============================================================
# 2. p3k14c DATABASE
# ============================================================
print(f"\n{'=' * 70}")
print("DATABASE 2: p3k14c RADIOCARBON DATES")
print("=" * 70)

p3k_pre = []
p3k_post = []
p3k_no_year = 0

with open(os.path.join(BASE, "data/p3k14c/p3k14c_data.csv")) as f:
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

        # Extract publication year using multiple patterns
        pub_year = None
        # Pattern 1: (YYYY) in reference
        m = re.findall(r'\((\d{4})\)', ref)
        if m:
            pub_year = int(m[-1])
        # Pattern 2: standalone YYYY in reference (Author YYYY, Author.YYYY, YYYY:page)
        if not pub_year:
            m = re.findall(r'(?:^|[.\s:])(\d{4})(?:[.\s,;:\-]|$)', ref)
            if m:
                pub_year = int(m[0])
        # Pattern 3: year in Source field (e.g., 'Goldberg_2016')
        if not pub_year:
            m = re.findall(r'(\d{4})', source)
            if m:
                pub_year = int(m[-1])

        if not pub_year or pub_year < 1950 or pub_year > 2025:
            p3k_no_year += 1
            continue

        site = {"lat": lat, "lon": lon, "name": row.get("SiteName", ""),
                "pub_year": pub_year, "country": row.get("Country", "")}

        if pub_year <= 2001:
            p3k_pre.append(site)
        else:
            p3k_post.append(site)

print(f"Pre-2001 (published ≤2001): {len(p3k_pre)}")
print(f"Post-2001 (published >2001): {len(p3k_post)}")
print(f"No year extracted: {p3k_no_year}")

print("\nComputing Z-scores for PRE-2001 p3k14c dates...")
p3k_pre_z = compute_z_scores(p3k_pre)
print("Computing Z-scores for POST-2001 p3k14c dates...")
p3k_post_z = compute_z_scores(p3k_post)

print(f"\n{'':>8} | {'PRE-2001 (n=' + str(len(p3k_pre)) + ')':>30} | {'POST-2001 (n=' + str(len(p3k_post)) + ')':>30}")
print(f"{'Thresh':>8} | {'Obs':>5} {'Mean':>7} {'Z':>7} {'Enrich':>7} | {'Obs':>5} {'Mean':>7} {'Z':>7} {'Enrich':>7}")
print("-" * 80)
for t in THRESHOLDS:
    pre = p3k_pre_z[t]
    post = p3k_post_z[t]
    print(f"  {t:>4}km | {pre['observed']:>5} {pre['random_mean']:>7.1f} {pre['z_score']:>7.2f} {pre['enrichment']:>6.2f}x"
          f" | {post['observed']:>5} {post['random_mean']:>7.1f} {post['z_score']:>7.2f} {post['enrichment']:>6.2f}x")


# ============================================================
# 3. MEGALITHIC PORTAL (whole dataset, no date split possible)
# ============================================================
print(f"\n{'=' * 70}")
print("DATABASE 3: MEGALITHIC PORTAL (reference — no date split available)")
print("=" * 70)

portal_sites = []
kml_dir = os.path.join(BASE, "data/megalithic_portal")
kml_files = glob.glob(os.path.join(kml_dir, "MegP_*.kml"))

for filepath in sorted(kml_files):
    filename = os.path.basename(filepath)
    site_type = filename.replace("MegP_", "").replace(".kml", "").replace("_", " ").strip()
    try:
        tree = ET.parse(filepath)
        root = tree.getroot()
        for pm in root.iter():
            if 'Placemark' not in pm.tag:
                continue
            name_el = coords_el = None
            for child in pm.iter():
                if 'name' in child.tag and child.tag.endswith('name') and name_el is None:
                    name_el = child
                if 'coordinates' in child.tag:
                    coords_el = child
            if coords_el is not None and coords_el.text:
                try:
                    parts = coords_el.text.strip().split(',')
                    lon, lat = float(parts[0]), float(parts[1])
                    if -90 <= lat <= 90 and -180 <= lon <= 180 and (lat != 0 or lon != 0):
                        portal_sites.append({"lat": lat, "lon": lon,
                                             "name": name_el.text if name_el else "Unknown",
                                             "type": site_type})
                except (ValueError, IndexError):
                    continue
    except Exception:
        continue

# Dedup
seen = set()
portal_deduped = []
for s in portal_sites:
    key = (round(s["lat"], 3), round(s["lon"], 3))
    if key not in seen:
        seen.add(key)
        portal_deduped.append(s)

print(f"Total Portal sites (deduped): {len(portal_deduped)}")
print("Computing Z-scores for full Portal dataset...")
portal_z = compute_z_scores(portal_deduped)

print(f"\n{'Thresh':>8} | {'Obs':>5} {'Mean':>7} {'Z':>7} {'Enrich':>7}")
print("-" * 45)
for t in THRESHOLDS:
    r = portal_z[t]
    print(f"  {t:>4}km | {r['observed']:>5} {r['random_mean']:>7.1f} {r['z_score']:>7.2f} {r['enrichment']:>6.2f}x")


# ============================================================
# 4. GEOGRAPHIC OVERLAP CHECK
# ============================================================
print(f"\n{'=' * 70}")
print("GEOGRAPHIC COMPARISON: Pre vs Post site distributions")
print("=" * 70)

def region_breakdown(sites):
    regions = Counter()
    for s in sites:
        lat, lon = s["lat"], s["lon"]
        if 25 < lat < 45 and -15 < lon < 45:
            regions["Mediterranean/Europe"] += 1
        elif 10 < lat < 45 and 25 < lon < 75:
            regions["Middle East/Central Asia"] += 1
        elif -60 < lat < 15 and -90 < lon < -30:
            regions["Americas"] += 1
        elif -40 < lat < 15 and -25 < lon < 55:
            regions["Africa"] += 1
        elif -10 < lat < 70 and 60 < lon < 180:
            regions["East/SE Asia & Oceania"] += 1
        else:
            regions["Other"] += 1
    return dict(regions)

print("\nPleiades pre-2001:")
for r, c in sorted(region_breakdown(pleiades_pre).items(), key=lambda x: -x[1]):
    print(f"  {c:>6} ({c/len(pleiades_pre)*100:5.1f}%) — {r}")
print("\nPleiades post-2001:")
for r, c in sorted(region_breakdown(pleiades_post).items(), key=lambda x: -x[1]):
    print(f"  {c:>6} ({c/len(pleiades_post)*100:5.1f}%) — {r}")

print("\np3k14c pre-2001:")
for r, c in sorted(region_breakdown(p3k_pre).items(), key=lambda x: -x[1]):
    print(f"  {c:>6} ({c/len(p3k_pre)*100:5.1f}%) — {r}")
print("\np3k14c post-2001:")
for r, c in sorted(region_breakdown(p3k_post).items(), key=lambda x: -x[1]):
    print(f"  {c:>6} ({c/len(p3k_post)*100:5.1f}%) — {r}")


# ============================================================
# 5. INTERPRETATION
# ============================================================
print(f"\n{'=' * 70}")
print("INTERPRETATION")
print("=" * 70)

# Compare Z-scores at 50km threshold
pleiades_pre_50z = pleiades_pre_z[50]["z_score"]
pleiades_post_50z = pleiades_post_z[50]["z_score"]
p3k_pre_50z = p3k_pre_z[50]["z_score"]
p3k_post_50z = p3k_post_z[50]["z_score"]

print(f"\nZ-score comparison at 50km threshold:")
print(f"  Pleiades pre-2001:  Z = {pleiades_pre_50z:+.2f}")
print(f"  Pleiades post-2001: Z = {pleiades_post_50z:+.2f}")
print(f"  p3k14c pre-2001:    Z = {p3k_pre_50z:+.2f}")
print(f"  p3k14c post-2001:   Z = {p3k_post_50z:+.2f}")

# Assess predictive power
if pleiades_post_50z >= pleiades_pre_50z * 0.5 and pleiades_post_50z > 2:
    pleiades_verdict = "PREDICTIVE — post-2001 signal is comparably strong"
elif pleiades_post_50z > 2:
    pleiades_verdict = "PARTIALLY PREDICTIVE — post-2001 signal present but weaker"
elif pleiades_post_50z > 0:
    pleiades_verdict = "WEAK — post-2001 signal exists but not significant"
else:
    pleiades_verdict = "SELECTION BIAS — signal absent in post-2001 data"

if p3k_post_50z >= p3k_pre_50z * 0.5 and p3k_post_50z > 2:
    p3k_verdict = "PREDICTIVE — post-2001 signal is comparably strong"
elif p3k_post_50z > 2:
    p3k_verdict = "PARTIALLY PREDICTIVE — post-2001 signal present but weaker"
elif p3k_post_50z > 0:
    p3k_verdict = "WEAK — post-2001 signal exists but not significant"
else:
    p3k_verdict = "SELECTION BIAS — signal absent in post-2001 data"

print(f"\nPleiades verdict: {pleiades_verdict}")
print(f"p3k14c verdict:   {p3k_verdict}")

# Note on Pleiades date proxy
print(f"\n** IMPORTANT CAVEAT **")
print(f"Pleiades 'created' dates reflect when entries were digitized (all ≥2009),")
print(f"not when sites were discovered. The 2009-2011 batch represents the")
print(f"Barrington Atlas of the Greek & Roman World (pub. 2000) — sites known")
print(f"well before Alison's 2001 proposal. Post-2011 entries are later additions")
print(f"to the database, many representing newly studied or more obscure sites.")
print(f"This is an imperfect but meaningful proxy for pre/post-2001 knowledge.")


# ============================================================
# 6. SAVE RESULTS
# ============================================================
output = {
    "meta": {
        "name": "Predictive Validation Test",
        "description": "Split databases pre/post 2001 (Alison's circle proposal) to test predictive power vs selection bias",
        "methodology": {
            "baseline": "Distribution-matched random (independent shuffle of real lats/lons + 2° Gaussian jitter)",
            "n_trials": N_TRIALS,
            "thresholds_km": THRESHOLDS,
            "random_seed": 42,
            "pleiades_split": "created ≤2011 (Barrington Atlas batch) vs created ≥2012 (later additions)",
            "p3k14c_split": "Publication year ≤2001 vs >2001, extracted from Reference field"
        }
    },
    "pleiades": {
        "pre_2001": {
            "n_sites": len(pleiades_pre),
            "split_criterion": "Pleiades created date 2009-2011 (Barrington Atlas digitization)",
            "results": {str(t): pleiades_pre_z[t] for t in THRESHOLDS},
            "geographic_distribution": region_breakdown(pleiades_pre)
        },
        "post_2001": {
            "n_sites": len(pleiades_post),
            "split_criterion": "Pleiades created date 2012+ (post-Atlas additions)",
            "results": {str(t): pleiades_post_z[t] for t in THRESHOLDS},
            "geographic_distribution": region_breakdown(pleiades_post)
        },
        "verdict": pleiades_verdict
    },
    "p3k14c": {
        "pre_2001": {
            "n_sites": len(p3k_pre),
            "split_criterion": "Reference publication year ≤ 2001",
            "results": {str(t): p3k_pre_z[t] for t in THRESHOLDS},
            "geographic_distribution": region_breakdown(p3k_pre)
        },
        "post_2001": {
            "n_sites": len(p3k_post),
            "split_criterion": "Reference publication year > 2001",
            "results": {str(t): p3k_post_z[t] for t in THRESHOLDS},
            "geographic_distribution": region_breakdown(p3k_post)
        },
        "verdict": p3k_verdict
    },
    "portal_reference": {
        "n_sites": len(portal_deduped),
        "note": "No date metadata in KML files — used as whole-dataset reference only",
        "results": {str(t): portal_z[t] for t in THRESHOLDS}
    },
    "z_score_summary_50km": {
        "pleiades_pre": pleiades_pre_50z,
        "pleiades_post": pleiades_post_50z,
        "p3k14c_pre": p3k_pre_50z,
        "p3k14c_post": p3k_post_50z,
        "portal_full": portal_z[50]["z_score"]
    }
}

outpath = os.path.join(BASE, "results/predictive_validation.json")
with open(outpath, "w") as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved to {outpath}")
