#!/usr/bin/env python3
"""
Test 3: Memphis Perpendicular Line Test

For each bearing 0°–175° (step 5°), draw a line through the necropolis centroid
and compute Pearson r between pyramid date_bce and perpendicular distance to that line.

If only the Great Circle bearing (~80°) produces a strong negative r, the
date-distance gradient is direction-specific. If most bearings do, it's
generic radial expansion.
"""

import json
import numpy as np
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Load pyramid data
with open('outputs/pyramid_date_distance/correlation.json') as f:
    data = json.load(f)

pyramids = data['pyramids']

# Extract arrays
lats = np.array([p['lat'] for p in pyramids])
lons = np.array([p['lon'] for p in pyramids])
dates = np.array([p['date_bce'] for p in pyramids])

# Centroid
centroid_lat = np.mean(lats)
centroid_lon = np.mean(lons)

print(f"Centroid: {centroid_lat:.4f}°N, {centroid_lon:.4f}°E")
print(f"N pyramids: {len(pyramids)}")


def perpendicular_distance_km(lat, lon, c_lat, c_lon, bearing_deg):
    """
    Signed perpendicular distance from point (lat, lon) to a great-circle line
    through (c_lat, c_lon) at the given bearing.

    Uses the cross-track distance formula:
      d_xt = asin(sin(d_13/R) * sin(θ_13 - θ_12)) * R

    where:
      d_13 = distance from center to point
      θ_13 = bearing from center to point
      θ_12 = bearing of the reference line

    Positive = right of line (looking along bearing), negative = left.
    """
    R = 6371.0  # Earth radius km

    lat1, lon1 = np.radians(c_lat), np.radians(c_lon)
    lat2, lon2 = np.radians(lat), np.radians(lon)

    # Angular distance from center to point
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    d_13 = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    # Bearing from center to point
    theta_13 = np.arctan2(
        np.sin(dlon) * np.cos(lat2),
        np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    )

    theta_12 = np.radians(bearing_deg)

    # Cross-track distance (signed)
    d_xt = np.arcsin(np.sin(d_13) * np.sin(theta_13 - theta_12)) * R

    return d_xt


# Scan bearings 0° to 175° in 5° steps
bearings = list(range(0, 180, 5))
results = []

for bearing in bearings:
    dists = [perpendicular_distance_km(p['lat'], p['lon'], centroid_lat, centroid_lon, bearing)
             for p in pyramids]
    r, p_val = pearsonr(dates, dists)
    results.append({
        'bearing': bearing,
        'r': round(r, 6),
        'p': round(p_val, 6),
        'abs_r': round(abs(r), 6)
    })

# Find best negative r
best_neg = min(results, key=lambda x: x['r'])
gc_result = next(r for r in results if r['bearing'] == 80)

# Summary stats
rs = [r['r'] for r in results]
neg_count = sum(1 for r in rs if r < -0.3)
strong_neg_bearings = [res['bearing'] for res in results if res['r'] < -0.3]

print(f"\n--- Bearing Scan Results ---")
print(f"GC bearing (80°): r = {gc_result['r']:.4f}, p = {gc_result['p']:.4f}")
print(f"Strongest negative r: bearing = {best_neg['bearing']}°, r = {best_neg['r']:.4f}, p = {best_neg['p']:.4f}")
print(f"Bearings with r < -0.3: {neg_count}/{len(bearings)} ({strong_neg_bearings})")
print(f"Mean r across all bearings: {np.mean(rs):.4f}")
print(f"Std r across all bearings: {np.std(rs):.4f}")
print(f"Min r: {min(rs):.4f}, Max r: {max(rs):.4f}")

# Determine verdict
if neg_count <= 6:  # ≤15° on each side of some peak
    if abs(best_neg['bearing'] - 80) <= 15 or abs(best_neg['bearing'] - 80 + 180) <= 15:
        verdict = "DIRECTION-SPECIFIC: Only GC bearing ±15° shows r < -0.3 — gradient is specific to GC orientation"
    else:
        verdict = "DIRECTION-SPECIFIC: Strong negative r is concentrated but NOT at GC bearing"
elif neg_count >= len(bearings) * 0.5:
    verdict = "GENERIC RADIAL: Most bearings show r < -0.3 — necropolis grew outward from center"
else:
    # Check if GC is near the peak
    if best_neg['bearing'] in range(65, 96):
        verdict = "MILD DIRECTIONAL: GC bearing is strongest but others are close — mild directional preference aligned with GC"
    else:
        verdict = "MILD DIRECTIONAL: Moderate concentration of negative r, peak not at GC bearing"

print(f"\nVERDICT: {verdict}")

# Save JSON
output = {
    'centroid_lat': round(centroid_lat, 6),
    'centroid_lon': round(centroid_lon, 6),
    'n_pyramids': len(pyramids),
    'gc_bearing_deg': 80,
    'bearing_scan': results,
    'gc_result': gc_result,
    'strongest_negative': best_neg,
    'bearings_with_r_below_neg03': strong_neg_bearings,
    'n_bearings_r_below_neg03': neg_count,
    'mean_r': round(float(np.mean(rs)), 6),
    'std_r': round(float(np.std(rs)), 6),
    'verdict': verdict
}

with open('outputs/resolution_tests/bearing_scan.json', 'w') as f:
    json.dump(output, f, indent=2)

print("\nSaved: outputs/resolution_tests/bearing_scan.json")

# --- Polar Plot ---
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': 'polar'})

# Convert bearings to radians (duplicate to full circle for visual symmetry)
theta = np.radians(bearings)
r_vals = [res['r'] for res in results]

# Mirror: bearing B and B+180 produce r and -r for signed distance,
# so show the full 360° by reflecting
theta_full = np.concatenate([theta, theta + np.pi])
r_full = np.concatenate([r_vals, [-r for r in r_vals]])

# Color by strength of negative correlation
colors = ['#d32f2f' if r < -0.3 else '#ff9800' if r < -0.15 else '#4caf50' for r in r_full]

# Plot as bars
bars = ax.bar(theta_full, [abs(r) for r in r_full], width=np.radians(4.5),
              bottom=0, alpha=0.7, color=colors, edgecolor='white', linewidth=0.5)

# Mark GC bearing
gc_theta = np.radians(80)
ax.annotate('GC\n80°', xy=(gc_theta, abs(gc_result['r'])),
            fontsize=11, fontweight='bold', ha='center', va='bottom',
            color='#1a237e')
ax.plot([gc_theta, gc_theta], [0, abs(gc_result['r'])], 'b-', linewidth=2, alpha=0.8)

# Mark GC bearing + 180
gc_theta2 = np.radians(260)
ax.plot([gc_theta2, gc_theta2], [0, abs(gc_result['r'])], 'b-', linewidth=2, alpha=0.8)

# Mark strongest
best_theta = np.radians(best_neg['bearing'])
ax.annotate(f"Peak\n{best_neg['bearing']}°\nr={best_neg['r']:.3f}",
            xy=(best_theta, abs(best_neg['r'])),
            fontsize=9, fontweight='bold', ha='center', va='bottom',
            color='#b71c1c')

ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_title('Pearson r (date vs perp. distance) by bearing through necropolis centroid\n'
             'Red = r < -0.3 | Orange = -0.3 to -0.15 | Green = > -0.15',
             fontsize=12, pad=20)
ax.set_rlabel_position(135)
ax.set_rmax(0.6)

plt.tight_layout()
plt.savefig('outputs/resolution_tests/bearing_scan_plot.png', dpi=150, bbox_inches='tight')
print("Saved: outputs/resolution_tests/bearing_scan_plot.png")
