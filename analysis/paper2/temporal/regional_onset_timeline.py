#!/usr/bin/env python3
"""
Regional Onset Timeline Along the Great Circle
=================================================
For each 30-degree arc segment of the Cyclopean Great Circle, compute the
earliest p3k14c radiocarbon date within 200 km of the circle.

This tests whether human activity onset progresses sequentially around the
circle (suggesting migration along it) or appears independently at different
times (suggesting independent settlement).

Output: results/regional_onset_timeline.json
"""

import csv, math, json, os
import numpy as np
from collections import defaultdict

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
THRESHOLD_KM = 200  # proximity to circle
BIN_WIDTH_YR = 500  # for density computation

BASE_DIR = os.path.expanduser("~/megalith_site_research")

SEGMENT_LABELS = {
    0: "0-30: North Africa / Sahara",
    30: "30-60: Egypt / Levant",
    60: "60-90: Egypt / Nile",
    90: "90-120: Iran / Mesopotamia",
    120: "120-150: South Asia / Indus",
    150: "150-180: SE Asia",
    180: "180-210: Pacific / PNG",
    210: "210-240: Pacific Ocean",
    240: "240-270: Pacific / Easter Island",
    270: "270-300: South America / Peru",
    300: "300-330: South America / Amazon / Atlantic",
    330: "330-360: Atlantic / North Africa",
}

# ============================================================
# HELPERS
# ============================================================
def gc_dist_from_circle(lat, lon):
    """Distance from a point to the Great Circle (km)."""
    lat1_r, lon1_r = math.radians(POLE_LAT), math.radians(POLE_LON)
    lat2_r, lon2_r = math.radians(lat), math.radians(lon)
    dlat = lat2_r - lat1_r
    dlon = lon2_r - lon1_r
    a = math.sin(dlat/2)**2 + math.cos(lat1_r)*math.cos(lat2_r)*math.sin(dlon/2)**2
    d = EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))
    return abs(d - QUARTER_CIRC)


def bearing_from_pole(lat, lon):
    """Bearing (0-360) from the pole to a point. This determines position along the circle."""
    lat1_r, lon1_r = math.radians(POLE_LAT), math.radians(POLE_LON)
    lat2_r, lon2_r = math.radians(lat), math.radians(lon)
    dlon = lon2_r - lon1_r
    x = math.sin(dlon) * math.cos(lat2_r)
    y = math.cos(lat1_r) * math.sin(lat2_r) - math.sin(lat1_r) * math.cos(lat2_r) * math.cos(dlon)
    bearing = math.degrees(math.atan2(x, y))
    return bearing % 360


def c14_to_cal(age_bp):
    """Rough conversion from radiocarbon years BP to calendar years BCE/CE."""
    if age_bp <= 0:
        return 1950
    if age_bp < 2500:
        cal_bp = age_bp * 1.0
    elif age_bp < 5000:
        cal_bp = age_bp * 1.05 + 50
    elif age_bp < 8000:
        cal_bp = age_bp * 1.08 + 100
    elif age_bp < 12000:
        cal_bp = age_bp * 1.12 + 200
    else:
        cal_bp = age_bp * 1.15 + 300
    return 1950 - cal_bp


# ============================================================
# LOAD p3k14c DATA (all dates, no classification filter)
# ============================================================
print("=" * 70)
print("REGIONAL ONSET TIMELINE ALONG THE GREAT CIRCLE")
print("=" * 70)

print("\nLoading p3k14c data...")
p3k_file = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")

# Load all individual dates (not deduplicated — we want all dates for density)
all_dates = []
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
        if lat == 0 and lon == 0:
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue
        if age_bp is None:
            continue
        all_dates.append({
            'lat': lat, 'lon': lon, 'age_bp': age_bp,
            'cal_year': c14_to_cal(age_bp),
            'site_name': row.get('SiteName', ''),
            'country': row.get('Country', ''),
        })

print(f"  Total dates loaded: {len(all_dates)}, skipped: {skipped}")

# ============================================================
# FILTER TO WITHIN 200 KM OF THE CIRCLE & ASSIGN SEGMENTS
# ============================================================
print(f"\nFiltering to within {THRESHOLD_KM} km of the Great Circle...")

segment_dates = defaultdict(list)  # segment_start -> list of cal_years
segment_sites = defaultdict(list)  # segment_start -> list of site details

for d in all_dates:
    dist = gc_dist_from_circle(d['lat'], d['lon'])
    if dist <= THRESHOLD_KM:
        bearing = bearing_from_pole(d['lat'], d['lon'])
        seg_start = int(bearing // 30) * 30
        segment_dates[seg_start].append(d['cal_year'])
        segment_sites[seg_start].append({
            'cal_year': d['cal_year'],
            'age_bp': d['age_bp'],
            'site_name': d['site_name'],
            'country': d['country'],
            'lat': d['lat'],
            'lon': d['lon'],
            'dist_km': round(dist, 1),
        })

total_near = sum(len(v) for v in segment_dates.values())
print(f"  Dates within {THRESHOLD_KM} km: {total_near}")

# ============================================================
# COMPUTE METRICS PER SEGMENT
# ============================================================
print(f"\n{'='*70}")
print(f"{'Segment':<45} {'N':>5} {'Earliest':>10} {'Median':>10} {'Onset':>10}")
print("-" * 85)

results_segments = []

for seg_start in range(0, 360, 30):
    label = SEGMENT_LABELS.get(seg_start, f"{seg_start}-{seg_start+30}")
    dates = segment_dates.get(seg_start, [])
    sites = segment_sites.get(seg_start, [])

    if not dates:
        print(f"  {label:<43} {'—':>5} {'—':>10} {'—':>10} {'—':>10}")
        results_segments.append({
            "segment_start": seg_start,
            "segment_end": seg_start + 30,
            "label": label,
            "n_dates": 0,
            "earliest_cal_year": None,
            "median_cal_year": None,
            "onset_cal_year": None,
            "earliest_site": None,
        })
        continue

    earliest = min(dates)
    median = float(np.median(dates))
    n = len(dates)

    # Compute density onset: first 500-year bin with >= 1 date per bin
    # Actually, find when density exceeds 1 date per 500-year bin
    # i.e., the start of the earliest 500-year window containing > 1 date
    # More precisely: bin all dates into 500-year bins, find the first bin
    # with count >= 2 (density > 1 per 500-year bin)
    bin_counts = defaultdict(int)
    for yr in dates:
        bin_start = int(yr // BIN_WIDTH_YR) * BIN_WIDTH_YR
        if yr < 0:
            bin_start = -int((-yr) // BIN_WIDTH_YR + 1) * BIN_WIDTH_YR
        bin_counts[bin_start] += 1

    # Find earliest bin with count >= 2 (onset of significant activity)
    onset = None
    for bin_start in sorted(bin_counts.keys()):
        if bin_counts[bin_start] >= 2:
            onset = bin_start
            break

    # Find the earliest site details
    earliest_site = min(sites, key=lambda s: s['cal_year'])

    def fmt_year(y):
        if y is None:
            return "—"
        if y < 0:
            return f"{int(abs(y))} BCE"
        return f"{int(y)} CE"

    print(f"  {label:<43} {n:>5} {fmt_year(earliest):>10} {fmt_year(median):>10} {fmt_year(onset):>10}")

    results_segments.append({
        "segment_start": seg_start,
        "segment_end": seg_start + 30,
        "label": label,
        "n_dates": n,
        "earliest_cal_year": round(earliest),
        "median_cal_year": round(median),
        "onset_cal_year": onset,
        "earliest_site": {
            "name": earliest_site['site_name'],
            "country": earliest_site['country'],
            "cal_year": round(earliest_site['cal_year']),
            "age_bp": earliest_site['age_bp'],
            "lat": round(earliest_site['lat'], 4),
            "lon": round(earliest_site['lon'], 4),
        },
        "bin_counts": {str(k): v for k, v in sorted(bin_counts.items())},
    })

# ============================================================
# SEQUENTIAL vs INDEPENDENT ANALYSIS
# ============================================================
print(f"\n{'='*70}")
print("SEQUENTIAL vs INDEPENDENT ONSET ANALYSIS")
print("=" * 70)

# Extract segments with onset dates
onset_data = [(r['segment_start'], r['onset_cal_year'], r['earliest_cal_year'])
              for r in results_segments if r['onset_cal_year'] is not None]

if len(onset_data) >= 3:
    # Sort by segment position (around the circle)
    onset_data.sort(key=lambda x: x[0])
    segments_with_onset = [x[0] for x in onset_data]
    onset_years = [x[1] for x in onset_data]
    earliest_years = [x[2] for x in onset_data]

    # Test 1: Correlation between segment position and onset year
    # If migration along the circle, we'd expect a correlation
    positions = np.array(segments_with_onset, dtype=float)
    onsets = np.array(onset_years, dtype=float)
    earliests = np.array(earliest_years, dtype=float)

    # Pearson correlation
    if np.std(onsets) > 0 and np.std(positions) > 0:
        onset_corr = float(np.corrcoef(positions, onsets)[0, 1])
        earliest_corr = float(np.corrcoef(positions, earliests)[0, 1])
    else:
        onset_corr = 0.0
        earliest_corr = 0.0

    # Test 2: Are adjacent segments more similar in onset than distant segments?
    adj_diffs = []
    non_adj_diffs = []
    for i in range(len(onset_data)):
        for j in range(i + 1, len(onset_data)):
            seg_diff = abs(onset_data[i][0] - onset_data[j][0])
            seg_diff = min(seg_diff, 360 - seg_diff)  # circular
            time_diff = abs(onset_data[i][1] - onset_data[j][1])
            if seg_diff <= 30:
                adj_diffs.append(time_diff)
            else:
                non_adj_diffs.append(time_diff)

    avg_adj = np.mean(adj_diffs) if adj_diffs else None
    avg_non_adj = np.mean(non_adj_diffs) if non_adj_diffs else None

    # Test 3: Range of onset years (narrow = simultaneous, wide = sequential)
    onset_range = max(onset_years) - min(onset_years)
    onset_std = float(np.std(onset_years))

    # Print onset sequence
    print("\nOnset sequence (sorted by position on circle):")
    for seg, onset_yr, earliest_yr in onset_data:
        label = SEGMENT_LABELS.get(seg, f"{seg}-{seg+30}")
        def fmt(y):
            return f"{int(abs(y))} BCE" if y < 0 else f"{int(y)} CE"
        print(f"  {label:<43} onset: {fmt(onset_yr):>10}  earliest: {fmt(earliest_yr):>10}")

    print(f"\nCorrelation (position vs onset year):    r = {onset_corr:.3f}")
    print(f"Correlation (position vs earliest year): r = {earliest_corr:.3f}")
    print(f"Onset year range: {int(onset_range)} years")
    print(f"Onset year std dev: {int(onset_std)} years")
    if avg_adj is not None and avg_non_adj is not None:
        print(f"Avg time diff (adjacent segments):     {int(avg_adj)} years")
        print(f"Avg time diff (non-adjacent segments): {int(avg_non_adj)} years")
        adjacency_ratio = avg_adj / avg_non_adj if avg_non_adj > 0 else None
        if adjacency_ratio is not None:
            print(f"Adjacency ratio: {adjacency_ratio:.2f} (< 1.0 = adjacent segments more similar)")

    # Interpretation
    interpretation = {}
    interpretation['onset_position_correlation'] = round(onset_corr, 3)
    interpretation['earliest_position_correlation'] = round(earliest_corr, 3)
    interpretation['onset_range_years'] = int(onset_range)
    interpretation['onset_std_years'] = int(onset_std)

    if avg_adj is not None and avg_non_adj is not None:
        interpretation['avg_adjacent_diff_years'] = int(avg_adj)
        interpretation['avg_non_adjacent_diff_years'] = int(avg_non_adj)
        interpretation['adjacency_ratio'] = round(avg_adj / avg_non_adj, 3) if avg_non_adj > 0 else None

    # Decision logic
    sequential_evidence = 0
    independent_evidence = 0

    if abs(onset_corr) > 0.5:
        sequential_evidence += 1
    else:
        independent_evidence += 1

    if avg_adj is not None and avg_non_adj is not None and avg_adj < avg_non_adj * 0.7:
        sequential_evidence += 1
    else:
        independent_evidence += 1

    if onset_range > 5000:
        # Wide range could support either — sequential if correlated, independent if not
        if abs(onset_corr) > 0.5:
            sequential_evidence += 1
        else:
            independent_evidence += 1
    elif onset_range < 2000:
        # Narrow range suggests roughly simultaneous (independent)
        independent_evidence += 1

    if sequential_evidence > independent_evidence:
        verdict = "SEQUENTIAL — onset progresses around the circle, suggesting migration or cultural diffusion along it"
    elif independent_evidence > sequential_evidence:
        verdict = "INDEPENDENT — onset appears at different times without clear directional progression, suggesting independent settlement"
    else:
        verdict = "AMBIGUOUS — evidence does not clearly favor sequential or independent onset"

    interpretation['verdict'] = verdict
    interpretation['sequential_evidence_score'] = sequential_evidence
    interpretation['independent_evidence_score'] = independent_evidence

    print(f"\nVERDICT: {verdict}")

else:
    interpretation = {"verdict": "INSUFFICIENT DATA — fewer than 3 segments have onset dates"}
    print(f"\n{interpretation['verdict']}")

# ============================================================
# SAVE RESULTS
# ============================================================
output = {
    "meta": {
        "analysis": "Regional Onset Timeline Along the Great Circle",
        "date": "2026-03-20",
        "methodology": {
            "data_source": "p3k14c radiocarbon database (all dates, no classification filter)",
            "proximity_threshold_km": THRESHOLD_KM,
            "segment_width_degrees": 30,
            "density_bin_width_years": BIN_WIDTH_YR,
            "onset_definition": "earliest 500-year bin with >= 2 dates",
            "c14_calibration": "simplified IntCal20 piecewise linear approximation",
            "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
        },
        "question": "Does onset of human activity progress sequentially around the circle (migration) or appear independently (independent settlement)?",
        "total_dates_loaded": len(all_dates),
        "dates_within_threshold": total_near,
    },
    "segments": results_segments,
    "interpretation": interpretation,
}

results_path = os.path.join(BASE_DIR, "results", "regional_onset_timeline.json")
with open(results_path, "w") as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nSaved: {results_path}")

print("\nDONE.")
