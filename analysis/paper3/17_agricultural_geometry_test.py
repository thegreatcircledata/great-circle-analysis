#!/usr/bin/env python3
"""
TEST 5: Solar/Agricultural Geometry
Does the circle trace a path through latitudes with optimal growing season
characteristics for early agriculture?

Uses a latitude-based growing season model (Koppen approximation) since
WorldClim rasters require large downloads. The key insight is that early
agriculture requires BOTH sufficient growing season AND sufficient
seasonality to drive storage and surplus planning.

Agricultural optimality index = growing_season_months * seasonality_factor
where seasonality_factor = std(monthly_temp) / mean(monthly_temp)
"""
import numpy as np
import json
import os

ALISON_POLE = (59.682122, -138.646087)
EARTH_R = 6371.0
N_MC = 1000

def to_xyz(lat, lon):
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    return np.array([np.cos(lat_r)*np.cos(lon_r), np.cos(lat_r)*np.sin(lon_r), np.sin(lat_r)])

def gc_points(pole_lat, pole_lon, n=360):
    pole = to_xyz(pole_lat, pole_lon)
    if abs(pole[2]) < 0.9:
        ref = np.array([0, 0, 1])
    else:
        ref = np.array([1, 0, 0])
    u = np.cross(pole, ref)
    u /= np.linalg.norm(u)
    v = np.cross(pole, u)
    v /= np.linalg.norm(v)
    angles = np.linspace(0, 2*np.pi, n, endpoint=False)
    points = np.outer(np.cos(angles), u) + np.outer(np.sin(angles), v)
    lats = np.degrees(np.arcsin(np.clip(points[:, 2], -1, 1)))
    lons = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    return lats, lons

def growing_season_months(lat):
    """Estimate months with mean temperature > 5°C based on latitude.
    Uses a simple sinusoidal temperature model."""
    abs_lat = abs(lat)
    if abs_lat > 80:
        return 0
    elif abs_lat < 10:
        return 12  # tropical — always warm
    else:
        # Approximate: growing season decreases roughly linearly from
        # 12 months at equator to 0 at poles, with an inflection around 60°
        # Based on Koppen climate classification patterns
        if abs_lat < 25:
            return 12
        elif abs_lat < 35:
            return 11
        elif abs_lat < 40:
            return 9
        elif abs_lat < 45:
            return 8
        elif abs_lat < 50:
            return 7
        elif abs_lat < 55:
            return 6
        elif abs_lat < 60:
            return 5
        elif abs_lat < 65:
            return 3
        elif abs_lat < 70:
            return 2
        else:
            return 1

def seasonality_factor(lat):
    """Temperature seasonality increases with latitude.
    Returns a 0-1 factor: 0 at equator (no seasonality), ~1 at 60°."""
    abs_lat = abs(lat)
    if abs_lat < 5:
        return 0.05
    elif abs_lat > 60:
        return 1.0
    else:
        return abs_lat / 60.0

def agricultural_optimality(lat):
    """Combined index: growing season × seasonality.
    Peaks where there's enough warm season for crops AND enough cold season
    to drive storage/surplus behavior.

    The 'Goldilocks zone' for early agriculture: ~25-35° latitude."""
    gs = growing_season_months(lat)
    sf = seasonality_factor(lat)
    return gs * sf

def compute_circle_stats(pole_lat, pole_lon, n_points=360):
    """Compute agricultural statistics along a great circle."""
    lats, lons = gc_points(pole_lat, pole_lon, n_points)

    gs_values = np.array([growing_season_months(lat) for lat in lats])
    sf_values = np.array([seasonality_factor(lat) for lat in lats])
    ao_values = np.array([agricultural_optimality(lat) for lat in lats])

    # Also compute the fraction of the circle in the "agricultural sweet spot" (20-40° latitude)
    sweet_spot = np.sum((np.abs(lats) >= 20) & (np.abs(lats) <= 40)) / len(lats)

    return {
        "mean_growing_season": float(np.mean(gs_values)),
        "mean_seasonality": float(np.mean(sf_values)),
        "mean_ag_optimality": float(np.mean(ao_values)),
        "max_ag_optimality": float(np.max(ao_values)),
        "sweet_spot_fraction": float(sweet_spot),
        "mean_abs_lat": float(np.mean(np.abs(lats))),
        "std_abs_lat": float(np.std(np.abs(lats)))
    }

def main():
    print("=" * 60)
    print("TEST 5: SOLAR/AGRICULTURAL GEOMETRY")
    print("=" * 60)

    # Alison circle
    alison = compute_circle_stats(*ALISON_POLE)
    print(f"\nAlison circle:")
    print(f"  Mean growing season: {alison['mean_growing_season']:.1f} months")
    print(f"  Mean seasonality: {alison['mean_seasonality']:.3f}")
    print(f"  Mean ag optimality: {alison['mean_ag_optimality']:.2f}")
    print(f"  Sweet spot (20-40°) fraction: {alison['sweet_spot_fraction']:.3f}")
    print(f"  Mean |latitude|: {alison['mean_abs_lat']:.1f}°")

    # Monte Carlo
    print(f"\nRunning {N_MC} random circles...")
    rng = np.random.default_rng(42)

    mc_gs = []
    mc_ao = []
    mc_ss = []
    mc_sf = []

    for i in range(N_MC):
        z = rng.uniform(-1, 1)
        phi = rng.uniform(0, 2*np.pi)
        pole_lat = np.degrees(np.arcsin(z))
        pole_lon = np.degrees(phi) - 180

        stats = compute_circle_stats(pole_lat, pole_lon, 180)
        mc_gs.append(stats['mean_growing_season'])
        mc_ao.append(stats['mean_ag_optimality'])
        mc_ss.append(stats['sweet_spot_fraction'])
        mc_sf.append(stats['mean_seasonality'])

    mc_gs = np.array(mc_gs)
    mc_ao = np.array(mc_ao)
    mc_ss = np.array(mc_ss)
    mc_sf = np.array(mc_sf)

    # P-values
    p_gs = float(np.mean(mc_gs >= alison['mean_growing_season']))
    p_ao = float(np.mean(mc_ao >= alison['mean_ag_optimality']))
    p_ss = float(np.mean(mc_ss >= alison['sweet_spot_fraction']))

    print(f"\n=== GROWING SEASON ===")
    print(f"Alison: {alison['mean_growing_season']:.1f} months")
    print(f"Random: {np.mean(mc_gs):.1f} ± {np.std(mc_gs):.1f}")
    print(f"P(random >= Alison): {p_gs:.4f}")

    print(f"\n=== AGRICULTURAL OPTIMALITY (growing season × seasonality) ===")
    print(f"Alison: {alison['mean_ag_optimality']:.2f}")
    print(f"Random: {np.mean(mc_ao):.2f} ± {np.std(mc_ao):.2f}")
    print(f"P(random >= Alison): {p_ao:.4f}")

    print(f"\n=== SWEET SPOT FRACTION (20-40° latitude) ===")
    print(f"Alison: {alison['sweet_spot_fraction']:.3f}")
    print(f"Random: {np.mean(mc_ss):.3f} ± {np.std(mc_ss):.3f}")
    print(f"P(random >= Alison): {p_ss:.4f}")

    results = {
        "analysis": "Solar/Agricultural Geometry Test",
        "model": "Latitude-based growing season approximation (no external data required)",
        "note": "Growing season estimated from latitude; real WorldClim data would give more precise results",
        "alison": alison,
        "mc": {
            "n_trials": N_MC,
            "growing_season": {
                "alison": round(alison['mean_growing_season'], 1),
                "mc_mean": round(float(np.mean(mc_gs)), 1),
                "mc_std": round(float(np.std(mc_gs)), 1),
                "p_value": round(p_gs, 4),
                "significant": p_gs < 0.05
            },
            "ag_optimality": {
                "alison": round(alison['mean_ag_optimality'], 2),
                "mc_mean": round(float(np.mean(mc_ao)), 2),
                "mc_std": round(float(np.std(mc_ao)), 2),
                "p_value": round(p_ao, 4),
                "significant": p_ao < 0.05
            },
            "sweet_spot_fraction": {
                "alison": round(alison['sweet_spot_fraction'], 3),
                "mc_mean": round(float(np.mean(mc_ss)), 3),
                "mc_std": round(float(np.std(mc_ss)), 3),
                "p_value": round(p_ss, 4),
                "significant": p_ss < 0.05
            }
        },
        "conclusion": ""
    }

    any_sig = p_gs < 0.05 or p_ao < 0.05 or p_ss < 0.05
    if any_sig:
        results["conclusion"] = "SIGNIFICANT: "
        if p_ao < 0.05:
            results["conclusion"] += f"Alison circle passes through more agriculturally optimal latitudes than {(1-p_ao)*100:.1f}% of random circles. "
        if p_ss < 0.05:
            results["conclusion"] += f"Alison has {alison['sweet_spot_fraction']*100:.0f}% of its path in the 20-40° agricultural sweet spot (random: {np.mean(mc_ss)*100:.0f}%). "
    else:
        results["conclusion"] = f"NOT SIGNIFICANT: Alison circle has typical agricultural geometry (optimality p={p_ao:.4f}, sweet spot p={p_ss:.4f})"

    out_path = os.path.expanduser("~/megalith_site_research/paper3/results/causal_mechanisms/solar_agricultural_test.json")
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)

    md = f"""# Test 5: Solar/Agricultural Geometry

## Result: {'SIGNIFICANT' if any_sig else 'NOT SIGNIFICANT'}

### Growing Season
- Alison mean: {alison['mean_growing_season']:.1f} months
- Random mean: {np.mean(mc_gs):.1f} ± {np.std(mc_gs):.1f} months
- p = {p_gs:.4f}

### Agricultural Optimality (growing season × seasonality)
- Alison: {alison['mean_ag_optimality']:.2f}
- Random: {np.mean(mc_ao):.2f} ± {np.std(mc_ao):.2f}
- p = {p_ao:.4f}

### Sweet Spot (20-40° latitude) Fraction
- Alison: {alison['sweet_spot_fraction']*100:.1f}%
- Random: {np.mean(mc_ss)*100:.1f}% ± {np.std(mc_ss)*100:.1f}%
- p = {p_ss:.4f}

### Conclusion
{results['conclusion']}

### Caveat
This test uses a latitude-based growing season approximation, not actual climate data.
Continental interiors at the same latitude can have very different growing seasons.
A more precise test would use WorldClim monthly temperature rasters.
"""

    md_path = os.path.expanduser("~/megalith_site_research/paper3/results/causal_mechanisms/solar_agricultural_test.md")
    with open(md_path, 'w') as f:
        f.write(md)

    print(f"\nSaved to {out_path}")
    print(results["conclusion"])

if __name__ == "__main__":
    main()
