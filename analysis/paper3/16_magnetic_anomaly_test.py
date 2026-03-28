#!/usr/bin/env python3
"""
TEST 6: Magnetic Anomaly / Inclination
Does the Great Circle follow a magnetic anomaly contour or constant inclination band?
Uses pyIGRF for the International Geomagnetic Reference Field.
"""
import numpy as np
import json
import os

try:
    import pyigrf
    HAS_PYIGRF = True
except ImportError:
    HAS_PYIGRF = False

ALISON_POLE = (59.682122, -138.646087)
EARTH_R = 6371.0
N_MC = 1000

def to_xyz(lat, lon):
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    return np.array([np.cos(lat_r)*np.cos(lon_r), np.cos(lat_r)*np.sin(lon_r), np.sin(lat_r)])

def gc_points(pole_lat, pole_lon, n=360):
    """Generate n evenly spaced points along a great circle."""
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

def compute_magnetic_stats(pole_lat, pole_lon, year=2020.0, n_points=360):
    """Compute magnetic field statistics along a great circle."""
    lats, lons = gc_points(pole_lat, pole_lon, n_points)

    inclinations = []
    declinations = []
    intensities = []

    for lat, lon in zip(lats, lons):
        try:
            # pyigrf.igrf_value returns (D, I, H, X, Y, Z, F)
            # D=declination, I=inclination, H=horizontal, X,Y,Z components, F=total field
            result = pyigrf.igrf_value(lat, lon, 0, year)
            inclinations.append(result[1])  # I (inclination/dip)
            declinations.append(result[0])  # D (declination)
            intensities.append(result[6])   # F (total field intensity)
        except:
            pass

    if len(inclinations) < 10:
        return None

    inclinations = np.array(inclinations)
    intensities = np.array(intensities)

    return {
        "mean_inclination": float(np.mean(inclinations)),
        "std_inclination": float(np.std(inclinations)),
        "cv_inclination": float(np.std(inclinations) / (abs(np.mean(inclinations)) + 1e-10)),
        "range_inclination": float(np.max(inclinations) - np.min(inclinations)),
        "mean_intensity": float(np.mean(intensities)),
        "std_intensity": float(np.std(intensities)),
        "cv_intensity": float(np.std(intensities) / (np.mean(intensities) + 1e-10)),
        "range_intensity": float(np.max(intensities) - np.min(intensities)),
        "n_points": len(inclinations)
    }

def main():
    print("=" * 60)
    print("TEST 6: MAGNETIC ANOMALY / INCLINATION")
    print("=" * 60)

    if not HAS_PYIGRF:
        print("ERROR: pyigrf not installed")
        return

    # Alison circle
    print("\nComputing Alison circle magnetic field...")
    alison_stats = compute_magnetic_stats(*ALISON_POLE)
    if alison_stats is None:
        print("ERROR: Could not compute magnetic field for Alison circle")
        return

    print(f"  Mean inclination: {alison_stats['mean_inclination']:.1f}°")
    print(f"  Inclination std: {alison_stats['std_inclination']:.1f}°")
    print(f"  Inclination range: {alison_stats['range_inclination']:.1f}°")
    print(f"  Inclination CV: {alison_stats['cv_inclination']:.3f}")
    print(f"  Mean intensity: {alison_stats['mean_intensity']:.0f} nT")
    print(f"  Intensity std: {alison_stats['std_intensity']:.0f} nT")
    print(f"  Intensity CV: {alison_stats['cv_intensity']:.3f}")

    # Monte Carlo
    print(f"\nRunning {N_MC} random circles...")
    rng = np.random.default_rng(42)

    mc_incl_cv = []
    mc_incl_std = []
    mc_intensity_cv = []
    mc_incl_range = []

    for i in range(N_MC):
        z = rng.uniform(-1, 1)
        phi = rng.uniform(0, 2*np.pi)
        pole_lat = np.degrees(np.arcsin(z))
        pole_lon = np.degrees(phi) - 180

        stats = compute_magnetic_stats(pole_lat, pole_lon, n_points=180)  # fewer points for speed
        if stats is not None:
            mc_incl_cv.append(stats['cv_inclination'])
            mc_incl_std.append(stats['std_inclination'])
            mc_intensity_cv.append(stats['cv_intensity'])
            mc_incl_range.append(stats['range_inclination'])

        if (i+1) % 100 == 0:
            print(f"  {i+1}/{N_MC}")

    mc_incl_cv = np.array(mc_incl_cv)
    mc_incl_std = np.array(mc_incl_std)
    mc_intensity_cv = np.array(mc_intensity_cv)
    mc_incl_range = np.array(mc_incl_range)

    # Test 1: Does Alison follow a constant inclination band? (low CV = following contour)
    p_incl_cv = float(np.mean(mc_incl_cv <= alison_stats['cv_inclination']))
    p_incl_range = float(np.mean(mc_incl_range <= alison_stats['range_inclination']))
    p_intensity_cv = float(np.mean(mc_intensity_cv <= alison_stats['cv_intensity']))

    print(f"\n=== INCLINATION CONTOUR TEST ===")
    print(f"Low CV = circle follows a constant inclination = potential migration contour")
    print(f"Alison inclination CV: {alison_stats['cv_inclination']:.3f}")
    print(f"Random mean CV: {np.mean(mc_incl_cv):.3f} ± {np.std(mc_incl_cv):.3f}")
    print(f"P(random CV ≤ Alison CV): {p_incl_cv:.4f}")
    print(f"  {'SIGNIFICANT — Alison has LOWER inclination variance than most circles' if p_incl_cv < 0.05 else 'NOT SIGNIFICANT'}")

    print(f"\n=== INCLINATION RANGE TEST ===")
    print(f"Alison inclination range: {alison_stats['range_inclination']:.1f}°")
    print(f"Random mean range: {np.mean(mc_incl_range):.1f}° ± {np.std(mc_incl_range):.1f}°")
    print(f"P(random range ≤ Alison range): {p_incl_range:.4f}")

    print(f"\n=== INTENSITY VARIANCE TEST ===")
    print(f"Alison intensity CV: {alison_stats['cv_intensity']:.3f}")
    print(f"Random mean CV: {np.mean(mc_intensity_cv):.3f} ± {np.std(mc_intensity_cv):.3f}")
    print(f"P(random CV ≤ Alison CV): {p_intensity_cv:.4f}")

    # Results
    results = {
        "analysis": "Magnetic Anomaly / Inclination Test",
        "data": "IGRF-13 via pyigrf (epoch 2020.0)",
        "alison": alison_stats,
        "mc": {
            "n_trials": N_MC,
            "inclination_cv": {
                "alison": round(alison_stats['cv_inclination'], 4),
                "mc_mean": round(float(np.mean(mc_incl_cv)), 4),
                "mc_std": round(float(np.std(mc_incl_cv)), 4),
                "p_value": round(p_incl_cv, 4),
                "significant": p_incl_cv < 0.05,
                "interpretation": "Low CV means circle follows constant inclination (potential migration contour)"
            },
            "inclination_range": {
                "alison": round(alison_stats['range_inclination'], 1),
                "mc_mean": round(float(np.mean(mc_incl_range)), 1),
                "p_value": round(p_incl_range, 4),
                "significant": p_incl_range < 0.05
            },
            "intensity_cv": {
                "alison": round(alison_stats['cv_intensity'], 4),
                "mc_mean": round(float(np.mean(mc_intensity_cv)), 4),
                "p_value": round(p_intensity_cv, 4),
                "significant": p_intensity_cv < 0.05
            }
        },
        "conclusion": ""
    }

    any_sig = p_incl_cv < 0.05 or p_incl_range < 0.05 or p_intensity_cv < 0.05
    if any_sig:
        results["conclusion"] = "SIGNIFICANT: "
        if p_incl_cv < 0.05:
            results["conclusion"] += f"Alison circle follows a more constant magnetic inclination than {(1-p_incl_cv)*100:.1f}% of random circles. "
        if p_intensity_cv < 0.05:
            results["conclusion"] += f"Alison circle has more constant field intensity than {(1-p_intensity_cv)*100:.1f}% of random circles."
    else:
        results["conclusion"] = f"NOT SIGNIFICANT: Alison circle has typical magnetic field properties (inclination CV p={p_incl_cv:.4f}, intensity CV p={p_intensity_cv:.4f})"

    out_path = os.path.expanduser("~/megalith_site_research/paper3/results/causal_mechanisms/magnetic_anomaly_test.json")
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)

    # Markdown
    sig_str = "SIGNIFICANT" if any_sig else "NOT SIGNIFICANT"
    md = f"""# Test 6: Magnetic Anomaly / Inclination

## Result: {sig_str}

### Inclination Contour Test
Does the circle follow a constant magnetic inclination (dip angle)?
- Alison inclination CV: {alison_stats['cv_inclination']:.3f}
- Random mean CV: {np.mean(mc_incl_cv):.3f} ± {np.std(mc_incl_cv):.3f}
- p-value: {p_incl_cv:.4f}

### Inclination Range Test
- Alison range: {alison_stats['range_inclination']:.1f}°
- Random mean range: {np.mean(mc_incl_range):.1f}°
- p-value: {p_incl_range:.4f}

### Intensity Variance Test
- Alison intensity CV: {alison_stats['cv_intensity']:.3f}
- Random mean CV: {np.mean(mc_intensity_cv):.3f}
- p-value: {p_intensity_cv:.4f}

### Conclusion
{results['conclusion']}
"""

    md_path = os.path.expanduser("~/megalith_site_research/paper3/results/causal_mechanisms/magnetic_anomaly_test.md")
    with open(md_path, 'w') as f:
        f.write(md)

    print(f"\nSaved to {out_path}")

if __name__ == "__main__":
    main()
