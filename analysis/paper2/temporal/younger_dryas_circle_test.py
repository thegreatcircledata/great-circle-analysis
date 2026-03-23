#!/usr/bin/env python3
"""
Younger Dryas Great Circle Test
================================
Tests whether radiocarbon date density near the Great Circle shows
anomalous behavior during the Younger Dryas (12,800–11,600 BP).

Approach:
  1. Bin ALL p3k14c dates globally by 500-year windows (20,000–5,000 BP)
  2. Bin dates within 200km of the Great Circle
  3. Compute ratio: near-circle / global density per bin
  4. Test whether 12,800–11,600 BP bins show anomalous ratio
  5. Compare to 100 random great circles
  6. Check for any pre-5000 BCE signal in raw density

Output: results/younger_dryas_test.json
"""

import csv, math, json, os, sys, time
import numpy as np

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
THRESHOLD_KM = 200
N_RANDOM_CIRCLES = 100

BASE_DIR = os.path.expanduser("~/megalith_site_research")
os.makedirs(os.path.join(BASE_DIR, "results"), exist_ok=True)

# 500-year bins from 20,000 BP to 5,000 BP
BIN_EDGES_BP = list(range(5000, 20500, 500))  # 5000, 5500, ..., 20000
BINS_BP = [(BIN_EDGES_BP[i], BIN_EDGES_BP[i+1]) for i in range(len(BIN_EDGES_BP)-1)]
# Each bin: [low_bp, high_bp) — e.g., [5000, 5500) means ages 5000–5499 BP

# Younger Dryas window (14C years BP — the YD onset is ~10,800 14C BP,
# which calibrates to ~12,800 cal BP; the YD end ~10,000 14C BP calibrates
# to ~11,600 cal BP. We test both 14C BP bins and approximate cal BP bins.)
YD_14C_START = 10000   # 14C years BP
YD_14C_END = 11000     # 14C years BP (the bin containing the YD)
YD_CAL_START = 11500   # cal BP
YD_CAL_END = 13000     # cal BP (bracket the YD onset generously)


# ============================================================
# VECTORIZED HELPERS
# ============================================================
def haversine_vec(lat1, lon1, lats, lons):
    """Haversine distance from single point to array of points (km)."""
    lat1_r, lon1_r = np.radians(lat1), np.radians(lon1)
    lat2_r = np.radians(lats)
    lon2_r = np.radians(lons)
    dlat = lat2_r - lat1_r
    dlon = lon2_r - lon1_r
    a = np.sin(dlat/2)**2 + np.cos(lat1_r)*np.cos(lat2_r)*np.sin(dlon/2)**2
    return EARTH_R_KM * 2 * np.arcsin(np.sqrt(np.minimum(1.0, a)))


def gc_dist_from_pole(pole_lat, pole_lon, lats, lons):
    """Distance from each point to the great circle defined by pole (km)."""
    d = haversine_vec(pole_lat, pole_lon, lats, lons)
    return np.abs(d - QUARTER_CIRC)


def c14_to_cal_bp(age_bp):
    """
    Rough 14C BP -> cal BP conversion (simplified IntCal20 approximation).
    """
    if age_bp <= 0:
        return 0
    if age_bp < 2500:
        return age_bp * 1.0
    elif age_bp < 5000:
        return age_bp * 1.05 + 50
    elif age_bp < 8000:
        return age_bp * 1.08 + 100
    elif age_bp < 12000:
        return age_bp * 1.12 + 200
    else:
        return age_bp * 1.15 + 300


def random_pole():
    """Generate a uniformly random point on the sphere (pole for random circle)."""
    z = np.random.uniform(-1, 1)
    theta = np.random.uniform(0, 2 * np.pi)
    lat = np.degrees(np.arcsin(z))
    lon = np.degrees(theta) - 180
    return lat, lon


# ============================================================
# LOAD p3k14c DATA — ALL DATES (no deduplication by site)
# ============================================================
print("=" * 70)
print("YOUNGER DRYAS GREAT CIRCLE TEST")
print("=" * 70)

print("\nLoading p3k14c data (all 170k+ dates)...")
p3k_file = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")

all_lats = []
all_lons = []
all_ages_bp = []
all_ages_cal_bp = []
skipped = 0

with open(p3k_file, encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            lat = float(row['Lat'])
            lon = float(row['Long'])
            age_bp = float(row['Age']) if row['Age'] else None
        except (ValueError, TypeError):
            skipped += 1
            continue
        if age_bp is None or age_bp <= 0:
            skipped += 1
            continue
        if lat == 0 and lon == 0:
            skipped += 1
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            skipped += 1
            continue

        all_lats.append(lat)
        all_lons.append(lon)
        all_ages_bp.append(age_bp)
        all_ages_cal_bp.append(c14_to_cal_bp(age_bp))

all_lats = np.array(all_lats)
all_lons = np.array(all_lons)
all_ages_bp = np.array(all_ages_bp)
all_ages_cal_bp = np.array(all_ages_cal_bp)

print(f"  Loaded: {len(all_lats)} dates, skipped: {skipped}")

# ============================================================
# STEP 1 & 2: BIN GLOBALLY AND NEAR-CIRCLE (14C BP bins)
# ============================================================
print("\nComputing distances from Great Circle...")
gc_dists = gc_dist_from_pole(POLE_LAT, POLE_LON, all_lats, all_lons)
near_mask = gc_dists <= THRESHOLD_KM
print(f"  Dates within {THRESHOLD_KM}km of circle: {np.sum(near_mask)} / {len(all_lats)}")

# Bin in 14C BP
print("\nBinning dates by 500-year windows (14C BP)...")
global_counts = {}
near_counts = {}

for low_bp, high_bp in BINS_BP:
    age_mask = (all_ages_bp >= low_bp) & (all_ages_bp < high_bp)
    global_counts[(low_bp, high_bp)] = int(np.sum(age_mask))
    near_counts[(low_bp, high_bp)] = int(np.sum(age_mask & near_mask))

# Also bin in cal BP (approximate)
global_counts_cal = {}
near_counts_cal = {}
for low_bp, high_bp in BINS_BP:
    age_mask = (all_ages_cal_bp >= low_bp) & (all_ages_cal_bp < high_bp)
    global_counts_cal[(low_bp, high_bp)] = int(np.sum(age_mask))
    near_counts_cal[(low_bp, high_bp)] = int(np.sum(age_mask & near_mask))

# ============================================================
# STEP 3: COMPUTE RATIOS
# ============================================================
print("\nComputing near-circle / global ratios...")

# Overall fraction of dates near the circle (baseline expectation)
overall_near_frac = np.sum(near_mask) / len(all_lats)
print(f"  Overall near-circle fraction: {overall_near_frac:.6f} ({np.sum(near_mask)}/{len(all_lats)})")

bins_table = []
for low_bp, high_bp in BINS_BP:
    g = global_counts[(low_bp, high_bp)]
    n = near_counts[(low_bp, high_bp)]
    ratio = n / g if g > 0 else 0.0
    # Normalized ratio: divide by overall fraction to get enrichment
    enrichment = ratio / overall_near_frac if overall_near_frac > 0 and g > 0 else 0.0

    g_cal = global_counts_cal[(low_bp, high_bp)]
    n_cal = near_counts_cal[(low_bp, high_bp)]
    ratio_cal = n_cal / g_cal if g_cal > 0 else 0.0
    enrichment_cal = ratio_cal / overall_near_frac if overall_near_frac > 0 and g_cal > 0 else 0.0

    bins_table.append({
        "bin_low_bp": low_bp,
        "bin_high_bp": high_bp,
        "label": f"{high_bp}–{low_bp} BP",
        "global_count_14c": g,
        "near_count_14c": n,
        "ratio_14c": round(ratio, 6),
        "enrichment_14c": round(enrichment, 3),
        "global_count_cal": g_cal,
        "near_count_cal": n_cal,
        "ratio_cal": round(ratio_cal, 6),
        "enrichment_cal": round(enrichment_cal, 3),
    })

# Print table
print(f"\n{'Bin (BP)':<20} {'Global':>8} {'Near':>6} {'Ratio':>10} {'Enrich':>8}  |  {'G(cal)':>8} {'N(cal)':>6} {'E(cal)':>8}")
print("-" * 95)
for b in bins_table:
    yd_marker = " <-- YD" if YD_14C_START <= b['bin_low_bp'] < YD_14C_END else ""
    print(f"  {b['label']:<18} {b['global_count_14c']:>8} {b['near_count_14c']:>6} "
          f"{b['ratio_14c']:>10.6f} {b['enrichment_14c']:>8.3f}  |  "
          f"{b['global_count_cal']:>8} {b['near_count_cal']:>6} {b['enrichment_cal']:>8.3f}{yd_marker}")

# ============================================================
# STEP 4: YOUNGER DRYAS ANOMALY TEST
# ============================================================
print("\n" + "=" * 70)
print("YOUNGER DRYAS ANOMALY TEST")
print("=" * 70)

# Define YD bins and non-YD bins (using 14C BP)
yd_bins_14c = [b for b in bins_table if b['bin_low_bp'] >= YD_14C_START and b['bin_high_bp'] <= YD_14C_END + 500]
non_yd_bins_14c = [b for b in bins_table if not (b['bin_low_bp'] >= YD_14C_START and b['bin_high_bp'] <= YD_14C_END + 500)]

# Also check cal BP YD window
yd_bins_cal = [b for b in bins_table if b['bin_low_bp'] >= YD_CAL_START and b['bin_high_bp'] <= YD_CAL_END + 500]
non_yd_bins_cal = [b for b in bins_table if not (b['bin_low_bp'] >= YD_CAL_START and b['bin_high_bp'] <= YD_CAL_END + 500)]

def yd_anomaly_stats(yd_bins, non_yd_bins, key_enrich, key_global):
    """Compare YD enrichment to non-YD enrichment."""
    # Only use bins with sufficient data
    yd_e = [b[key_enrich] for b in yd_bins if b[key_global] >= 50]
    non_yd_e = [b[key_enrich] for b in non_yd_bins if b[key_global] >= 50]

    if not yd_e or not non_yd_e:
        return {"insufficient_data": True}

    yd_mean = np.mean(yd_e)
    non_yd_mean = np.mean(non_yd_e)
    non_yd_std = np.std(non_yd_e)
    z_score = (yd_mean - non_yd_mean) / non_yd_std if non_yd_std > 0 else 0

    return {
        "yd_mean_enrichment": round(float(yd_mean), 4),
        "non_yd_mean_enrichment": round(float(non_yd_mean), 4),
        "non_yd_std": round(float(non_yd_std), 4),
        "z_score": round(float(z_score), 3),
        "yd_bins_used": len(yd_e),
        "non_yd_bins_used": len(non_yd_e),
        "interpretation": (
            "significant dip" if z_score < -2 else
            "marginal dip" if z_score < -1 else
            "significant spike" if z_score > 2 else
            "marginal spike" if z_score > 1 else
            "no significant anomaly"
        ),
    }

yd_14c_result = yd_anomaly_stats(yd_bins_14c, non_yd_bins_14c, 'enrichment_14c', 'global_count_14c')
yd_cal_result = yd_anomaly_stats(yd_bins_cal, non_yd_bins_cal, 'enrichment_cal', 'global_count_cal')

print(f"\n  14C BP bins in YD window ({YD_14C_START}–{YD_14C_END} BP):")
for k, v in yd_14c_result.items():
    print(f"    {k}: {v}")

print(f"\n  Cal BP bins in YD window ({YD_CAL_START}–{YD_CAL_END} cal BP):")
for k, v in yd_cal_result.items():
    print(f"    {k}: {v}")

# ============================================================
# STEP 5: COMPARE TO 100 RANDOM GREAT CIRCLES
# ============================================================
print("\n" + "=" * 70)
print("RANDOM CIRCLE COMPARISON (100 circles)")
print("=" * 70)

np.random.seed(42)

random_circle_results = []
t0 = time.time()

for i in range(N_RANDOM_CIRCLES):
    pole_lat, pole_lon = random_pole()
    dists = gc_dist_from_pole(pole_lat, pole_lon, all_lats, all_lons)
    near = dists <= THRESHOLD_KM
    near_frac = np.sum(near) / len(all_lats)

    # Bin and compute enrichments (14C BP)
    enrichments_14c = []
    for low_bp, high_bp in BINS_BP:
        age_mask = (all_ages_bp >= low_bp) & (all_ages_bp < high_bp)
        g = int(np.sum(age_mask))
        n = int(np.sum(age_mask & near))
        ratio = n / g if g > 0 else 0
        e = ratio / near_frac if near_frac > 0 and g > 0 else 0
        enrichments_14c.append(e)

    # YD enrichment for this random circle
    yd_indices = [j for j, (lo, hi) in enumerate(BINS_BP)
                  if lo >= YD_14C_START and hi <= YD_14C_END + 500
                  and global_counts[(lo, hi)] >= 50]
    non_yd_indices = [j for j, (lo, hi) in enumerate(BINS_BP)
                      if not (lo >= YD_14C_START and hi <= YD_14C_END + 500)
                      and global_counts[(lo, hi)] >= 50]

    if yd_indices and non_yd_indices:
        yd_e = np.mean([enrichments_14c[j] for j in yd_indices])
        non_yd_e = np.mean([enrichments_14c[j] for j in non_yd_indices])
        non_yd_s = np.std([enrichments_14c[j] for j in non_yd_indices])
        yd_z = (yd_e - non_yd_e) / non_yd_s if non_yd_s > 0 else 0
    else:
        yd_e = yd_z = float('nan')

    random_circle_results.append({
        "pole_lat": round(pole_lat, 4),
        "pole_lon": round(pole_lon, 4),
        "n_near": int(np.sum(near)),
        "near_fraction": round(float(near_frac), 6),
        "yd_mean_enrichment": round(float(yd_e), 4),
        "yd_z_score": round(float(yd_z), 3),
        "all_enrichments_14c": [round(e, 4) for e in enrichments_14c],
    })

    if (i + 1) % 25 == 0:
        print(f"  Processed {i+1}/100 random circles ({time.time()-t0:.1f}s)")

# Alison circle's YD z-score vs random distribution
alison_yd_z = yd_14c_result.get('z_score', float('nan'))
random_yd_zs = [r['yd_z_score'] for r in random_circle_results if not np.isnan(r['yd_z_score'])]

if random_yd_zs:
    pct_more_extreme = sum(1 for z in random_yd_zs if abs(z) >= abs(alison_yd_z)) / len(random_yd_zs)
    random_yd_z_mean = np.mean(random_yd_zs)
    random_yd_z_std = np.std(random_yd_zs)
    # How many random circles show a "significant" YD anomaly?
    n_sig_dip = sum(1 for z in random_yd_zs if z < -2)
    n_sig_spike = sum(1 for z in random_yd_zs if z > 2)
else:
    pct_more_extreme = random_yd_z_mean = random_yd_z_std = float('nan')
    n_sig_dip = n_sig_spike = 0

random_comparison = {
    "n_random_circles": N_RANDOM_CIRCLES,
    "alison_circle_yd_z": alison_yd_z,
    "random_yd_z_mean": round(float(random_yd_z_mean), 3),
    "random_yd_z_std": round(float(random_yd_z_std), 3),
    "pct_random_more_extreme": round(float(pct_more_extreme), 3),
    "n_random_with_sig_yd_dip": n_sig_dip,
    "n_random_with_sig_yd_spike": n_sig_spike,
    "interpretation": (
        "Alison circle shows NO unique YD anomaly vs random circles"
        if pct_more_extreme > 0.05 else
        "Alison circle shows a RARE YD anomaly (p < 0.05 vs random circles)"
    ),
}

print(f"\n  Alison circle YD Z-score: {alison_yd_z}")
print(f"  Random circles — mean YD Z: {random_yd_z_mean:.3f} +/- {random_yd_z_std:.3f}")
print(f"  % random more extreme: {pct_more_extreme:.1%}")
print(f"  Random circles with sig YD dip (Z<-2): {n_sig_dip}/100")
print(f"  Random circles with sig YD spike (Z>2): {n_sig_spike}/100")
print(f"  → {random_comparison['interpretation']}")

# ============================================================
# STEP 6: ANY SIGNAL BEFORE 5000 BCE?
# ============================================================
print("\n" + "=" * 70)
print("PRE-5000 BCE SIGNAL CHECK")
print("=" * 70)

# 5000 BCE ≈ 6950 BP (14C) or ~7000 cal BP
# Check enrichment pattern across all ancient bins
pre_5000bce_threshold = 7000  # 14C BP, roughly 5000 BCE

ancient_bins = [b for b in bins_table if b['bin_low_bp'] >= pre_5000bce_threshold]
recent_bins = [b for b in bins_table if b['bin_low_bp'] < pre_5000bce_threshold]

ancient_enrichments = [b['enrichment_14c'] for b in ancient_bins if b['global_count_14c'] >= 50]
recent_enrichments = [b['enrichment_14c'] for b in recent_bins if b['global_count_14c'] >= 50]

if ancient_enrichments and recent_enrichments:
    ancient_mean = np.mean(ancient_enrichments)
    recent_mean = np.mean(recent_enrichments)
    all_enrichments = ancient_enrichments + recent_enrichments
    overall_std = np.std(all_enrichments)

    # Per-bin Z-scores relative to all bins
    overall_mean = np.mean(all_enrichments)

    print(f"\n  Enrichment — ancient (>{pre_5000bce_threshold} BP): mean={ancient_mean:.4f} (n={len(ancient_enrichments)} bins)")
    print(f"  Enrichment — recent (<{pre_5000bce_threshold} BP): mean={recent_mean:.4f} (n={len(recent_enrichments)} bins)")
    print(f"  Overall mean: {overall_mean:.4f}, std: {overall_std:.4f}")

    # Find any bins with notably high/low enrichment
    print(f"\n  Per-bin enrichments (ancient, >50 global dates):")
    for b in ancient_bins:
        if b['global_count_14c'] >= 50:
            e = b['enrichment_14c']
            z = (e - overall_mean) / overall_std if overall_std > 0 else 0
            marker = " ***" if abs(z) > 2 else " *" if abs(z) > 1.5 else ""
            print(f"    {b['label']:<20} global={b['global_count_14c']:>6}  near={b['near_count_14c']:>4}  "
                  f"enrichment={e:.3f}  z={z:+.2f}{marker}")

    pre_signal = {
        "ancient_mean_enrichment": round(float(ancient_mean), 4),
        "recent_mean_enrichment": round(float(recent_mean), 4),
        "overall_mean": round(float(overall_mean), 4),
        "overall_std": round(float(overall_std), 4),
        "ancient_vs_recent_diff": round(float(ancient_mean - recent_mean), 4),
        "notable_ancient_bins": [],
    }

    for b in ancient_bins:
        if b['global_count_14c'] >= 50:
            e = b['enrichment_14c']
            z = (e - overall_mean) / overall_std if overall_std > 0 else 0
            if abs(z) > 1.5:
                pre_signal['notable_ancient_bins'].append({
                    "bin": b['label'],
                    "enrichment": round(e, 4),
                    "z_from_mean": round(float(z), 3),
                    "global_count": b['global_count_14c'],
                    "near_count": b['near_count_14c'],
                })

    if pre_signal['notable_ancient_bins']:
        pre_signal['interpretation'] = (
            f"Found {len(pre_signal['notable_ancient_bins'])} ancient bins with notable "
            f"enrichment anomalies (|Z|>1.5). See notable_ancient_bins."
        )
    else:
        pre_signal['interpretation'] = (
            "No ancient bins show notable enrichment anomalies (all |Z| < 1.5). "
            "The Great Circle preference appears to be a Holocene phenomenon."
        )
else:
    pre_signal = {"insufficient_data": True}

# ============================================================
# SAVE RESULTS
# ============================================================
print("\n" + "=" * 70)
print("SAVING RESULTS")
print("=" * 70)

output = {
    "meta": {
        "analysis": "Younger Dryas Great Circle Density Test",
        "date": "2026-03-20",
        "description": (
            "Tests whether radiocarbon date density near the Alison Great Circle "
            "shows anomalous behavior during the Younger Dryas cold period "
            "(12,800–11,600 cal BP / ~10,800–10,000 14C BP)."
        ),
        "methodology": {
            "database": "p3k14c (170,150 dates)",
            "n_dates_used": int(len(all_lats)),
            "threshold_km": THRESHOLD_KM,
            "n_near_circle": int(np.sum(near_mask)),
            "overall_near_fraction": round(float(overall_near_frac), 6),
            "bins": "500-year windows from 20,000 BP to 5,000 BP (14C years BP)",
            "n_bins": len(BINS_BP),
            "yd_14c_window": f"{YD_14C_START}–{YD_14C_END} 14C BP",
            "yd_cal_window": f"{YD_CAL_START}–{YD_CAL_END} cal BP",
            "cal_conversion": "simplified IntCal20 polynomial approximation",
            "n_random_circles": N_RANDOM_CIRCLES,
            "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
        },
    },
    "bins_14c_bp": bins_table,
    "younger_dryas_test": {
        "14c_bp_bins": yd_14c_result,
        "cal_bp_bins": yd_cal_result,
    },
    "random_circle_comparison": random_comparison,
    "random_circle_details": [
        {k: v for k, v in r.items() if k != 'all_enrichments_14c'}
        for r in random_circle_results
    ],
    "pre_5000bce_signal": pre_signal,
    "summary": {
        "q1_yd_anomaly_on_circle": yd_14c_result.get('interpretation', 'unknown'),
        "q2_unique_vs_random": random_comparison['interpretation'],
        "q3_pre_5000bce_signal": pre_signal.get('interpretation', 'unknown'),
    },
}

out_path = os.path.join(BASE_DIR, "results", "younger_dryas_test.json")
with open(out_path, "w") as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nSaved: {out_path}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
for k, v in output['summary'].items():
    print(f"  {k}: {v}")
print("\nDONE.")
