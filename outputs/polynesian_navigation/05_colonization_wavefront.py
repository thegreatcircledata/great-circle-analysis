#!/usr/bin/env python3
"""
Analysis 4: Colonization Wavefront and the Great Circle
=========================================================
Directive 09, Script 05

Tests whether islands closer to the Great Circle were colonized earlier
than islands farther away.

Methods:
  A. Spearman correlation: colonization date vs distance to circle
  B. Monte Carlo: 10,000 random circles → distribution of Spearman rho
  C. Visualization: timeline colored by circle distance

Output:
  - colonization_vs_distance.json
  - colonization_timeline.png
"""

import json
import math
import os
import sys

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ─── Import data ─────────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from polynesian_voyaging_data import COLONIZATION_CHRONOLOGY

# ─── Constants ────────────────────────────────────────────────────────
POLE_LAT = 59.682122
POLE_LON = -138.646087
R_EARTH = 6371.0
QUARTER_CIRC = R_EARTH * math.pi / 2
N_MC = 10000


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1 = math.radians(lat1), math.radians(lon1)
    la2, lo2 = math.radians(lat2), math.radians(lon2)
    dlat = la2 - la1
    dlon = lo2 - lo1
    a = math.sin(dlat/2)**2 + math.cos(la1)*math.cos(la2)*math.sin(dlon/2)**2
    return R_EARTH * 2 * math.asin(math.sqrt(min(1.0, a)))


def dist_from_gc(lat, lon, pole_lat, pole_lon):
    d = haversine_km(lat, lon, pole_lat, pole_lon)
    qc = R_EARTH * math.pi / 2
    return abs(d - qc)


# ═══════════════════════════════════════════════════════════════════════
# Compile data
# ═══════════════════════════════════════════════════════════════════════

print("=" * 70)
print("COLONIZATION WAVEFRONT vs GREAT CIRCLE DISTANCE")
print("=" * 70)

island_data = []
for entry in COLONIZATION_CHRONOLOGY:
    lat, lon = entry["centroid"]
    d = dist_from_gc(lat, lon, POLE_LAT, POLE_LON)
    island_data.append({
        "name": entry["island_group"],
        "lat": lat,
        "lon": lon,
        "colonization_ce": entry["colonization_ce"],
        "date_range_ce": entry["date_range_ce"],
        "distance_to_circle_km": round(d, 1),
    })

# Sort by colonization date
island_data.sort(key=lambda x: x["colonization_ce"])

print(f"\n{'Island Group':35s}  {'Date (CE)':>10s}  {'Dist (km)':>10s}")
print(f"{'-'*35}  {'-'*10}  {'-'*10}")
for ig in island_data:
    print(f"  {ig['name']:33s}  {ig['colonization_ce']:>10d}  {ig['distance_to_circle_km']:>10.1f}")


# ═══════════════════════════════════════════════════════════════════════
# Spearman correlation: distance vs colonization date
# ═══════════════════════════════════════════════════════════════════════

dates = np.array([ig["colonization_ce"] for ig in island_data])
distances = np.array([ig["distance_to_circle_km"] for ig in island_data])

rho, spearman_p = stats.spearmanr(distances, dates)

print(f"\n--- Spearman correlation (distance → colonization date) ---")
print(f"  ρ = {rho:.4f}")
print(f"  p-value = {spearman_p:.6f}")
if rho < 0 and spearman_p < 0.05:
    print(f"  → Closer islands colonized LATER (negative correlation)")
elif rho > 0 and spearman_p < 0.05:
    print(f"  → Closer islands colonized EARLIER (positive — but wrong direction!)")
    print(f"  → Wait: positive rho means FARTHER islands colonized LATER")
    print(f"       i.e., distance and date increase together")
else:
    print(f"  → No significant correlation")

# Also test: Pearson correlation
pearson_r, pearson_p = stats.pearsonr(distances, dates)
print(f"\n  Pearson r = {pearson_r:.4f}, p = {pearson_p:.6f}")

# Interpretation note
print(f"\n  Interpretation:")
if rho < 0:
    print(f"    Negative ρ: islands CLOSER to the circle were colonized LATER.")
    print(f"    This is OPPOSITE to what we'd expect if the circle traces a colonization route.")
elif rho > 0:
    print(f"    Positive ρ: islands FARTHER from the circle were colonized LATER.")
    print(f"    This is consistent with the circle tracing an early colonization corridor,")
    print(f"    but note the confound that West Polynesia (colonized earliest) happens to be")
    print(f"    far from the circle, while East Polynesia (colonized later) includes Easter Island")
    print(f"    which is on the circle.")
else:
    print(f"    No clear pattern.")


# ═══════════════════════════════════════════════════════════════════════
# Separate analysis: East Polynesia only
# ═══════════════════════════════════════════════════════════════════════

print(f"\n--- East Polynesia only (post-AD 1000) ---")

east_poly = [ig for ig in island_data if ig["colonization_ce"] >= 1000]
if len(east_poly) >= 3:
    ep_dates = np.array([ig["colonization_ce"] for ig in east_poly])
    ep_dists = np.array([ig["distance_to_circle_km"] for ig in east_poly])
    ep_rho, ep_p = stats.spearmanr(ep_dists, ep_dates)
    print(f"  n = {len(east_poly)} island groups")
    print(f"  Spearman ρ = {ep_rho:.4f}, p = {ep_p:.6f}")
    if ep_p < 0.05:
        print(f"  → {'Closer colonized earlier' if ep_rho > 0 else 'Closer colonized later'}")
    else:
        print(f"  → No significant correlation within East Polynesia")
else:
    ep_rho, ep_p = None, None
    print(f"  Too few data points")


# ═══════════════════════════════════════════════════════════════════════
# Monte Carlo: random great circles
# ═══════════════════════════════════════════════════════════════════════

print(f"\n--- Monte Carlo: {N_MC} random great circles ---")

np.random.seed(42)
mc_rhos = np.zeros(N_MC)

for i in range(N_MC):
    # Random pole
    z = np.random.uniform(-1, 1)
    phi = np.random.uniform(0, 2 * math.pi)
    pole_lat_r = math.degrees(math.asin(z))
    pole_lon_r = math.degrees(phi) - 180

    mc_dists = np.array([dist_from_gc(ig["lat"], ig["lon"], pole_lat_r, pole_lon_r)
                         for ig in island_data])
    mc_rhos[i], _ = stats.spearmanr(mc_dists, dates)

# How does actual rho compare?
if rho < 0:
    mc_p_value = np.mean(mc_rhos <= rho)  # fraction with even more negative rho
else:
    mc_p_value = np.mean(mc_rhos >= rho)  # fraction with even more positive rho

print(f"  Actual Spearman ρ: {rho:.4f}")
print(f"  MC mean ρ: {np.mean(mc_rhos):.4f} ± {np.std(mc_rhos):.4f}")
print(f"  MC range: [{np.min(mc_rhos):.4f}, {np.max(mc_rhos):.4f}]")
print(f"  Percentile rank: {np.mean(mc_rhos <= rho)*100:.1f}%")
print(f"  MC p-value (actual more extreme): {mc_p_value:.4f}")


# ═══════════════════════════════════════════════════════════════════════
# Visualization
# ═══════════════════════════════════════════════════════════════════════

print(f"\nGenerating colonization timeline plot...")

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Left: Scatter plot — date vs distance
ax = axes[0]
scatter = ax.scatter(distances, dates, s=100, c=distances, cmap='RdYlBu_r',
                     edgecolors='black', linewidth=0.5, zorder=5)

# Labels
for ig in island_data:
    ax.annotate(ig["name"].replace(" (Rapa Nui)", "").replace(" (Aotearoa)", ""),
                (ig["distance_to_circle_km"], ig["colonization_ce"]),
                xytext=(5, 5), textcoords="offset points", fontsize=7)

# Error bars for date ranges
for ig in island_data:
    dr = ig["date_range_ce"]
    ax.plot([ig["distance_to_circle_km"]]*2, [dr[0], dr[1]],
            'k-', linewidth=0.5, alpha=0.5, zorder=3)

# Regression line
if len(distances) > 2:
    z = np.polyfit(distances, dates, 1)
    p = np.poly1d(z)
    x_range = np.linspace(distances.min(), distances.max(), 100)
    ax.plot(x_range, p(x_range), 'r--', alpha=0.5,
            label=f'ρ={rho:.3f}, p={spearman_p:.3f}')

ax.set_xlabel("Distance to Great Circle (km)", fontsize=11)
ax.set_ylabel("Colonization Date (CE)", fontsize=11)
ax.set_title("Polynesian Colonization Date\nvs Distance to Great Circle", fontsize=13)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.colorbar(scatter, ax=ax, label="Distance to GC (km)")

# Right: Monte Carlo distribution of rho
ax2 = axes[1]
ax2.hist(mc_rhos, bins=50, density=True, alpha=0.7, color='steelblue',
         edgecolor='black', linewidth=0.5)
ax2.axvline(rho, color='red', linewidth=2, linestyle='--',
            label=f'Actual ρ = {rho:.3f}')
ax2.axvline(0, color='gray', linewidth=1, linestyle=':')

ax2.set_xlabel("Spearman ρ (distance vs date)", fontsize=11)
ax2.set_ylabel("Density", fontsize=11)
ax2.set_title(f"Monte Carlo Distribution of ρ\n(n={N_MC} random great circles)", fontsize=13)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig("colonization_timeline.png", dpi=150)
print("Saved colonization_timeline.png")


# ═══════════════════════════════════════════════════════════════════════
# Save JSON
# ═══════════════════════════════════════════════════════════════════════

output = {
    "analysis": "Colonization wavefront vs Great Circle distance",
    "pole": {"lat": POLE_LAT, "lon": POLE_LON},
    "n_islands": len(island_data),
    "island_data": island_data,
    "spearman_test": {
        "rho": round(float(rho), 4),
        "p_value": round(float(spearman_p), 6),
        "significant": bool(spearman_p < 0.05),
        "direction": "positive" if rho > 0 else "negative",
        "interpretation": (
            "Positive ρ: farther islands colonized later (distance grows with date)"
            if rho > 0 else
            "Negative ρ: closer islands colonized later"
        ),
    },
    "pearson_test": {
        "r": round(float(pearson_r), 4),
        "p_value": round(float(pearson_p), 6),
    },
    "east_polynesia_only": {
        "n": len(east_poly),
        "rho": round(float(ep_rho), 4) if ep_rho is not None else None,
        "p_value": round(float(ep_p), 6) if ep_p is not None else None,
    },
    "monte_carlo": {
        "n_trials": N_MC,
        "actual_rho": round(float(rho), 4),
        "mc_mean_rho": round(float(np.mean(mc_rhos)), 4),
        "mc_std_rho": round(float(np.std(mc_rhos)), 4),
        "mc_p_value": round(float(mc_p_value), 4),
        "percentile": round(float(np.mean(mc_rhos <= rho) * 100), 1),
    },
    "caveat": (
        "The colonization chronology is debated, especially early dates. "
        "West Polynesia (Fiji/Tonga/Samoa, ~1000-800 BCE) was colonized millennia "
        "before East Polynesia (~AD 1000-1300), creating a bimodal distribution. "
        "The 'long pause' between West and East Polynesian settlement (~1800 years) "
        "dominates the variance. Easter Island's position on the circle is a single "
        "data point and its late colonization date (~AD 1230) works against the "
        "hypothesis that proximity to the circle predicts earlier settlement."
    ),
}

with open("colonization_vs_distance.json", "w") as f:
    json.dump(output, f, indent=2)

print("Saved colonization_vs_distance.json")
