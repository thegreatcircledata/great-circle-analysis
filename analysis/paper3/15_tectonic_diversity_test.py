#!/usr/bin/env python3
"""
TEST 4: Tectonic Diversity
Does the Great Circle cross more tectonic plate boundaries per unit distance
than random great circles?

Uses Bird (2003) plate boundary dataset from GitHub.
"""
import numpy as np
import json
import os
import time
import urllib.request

ALISON_POLE = (59.682122, -138.646087)
EARTH_R = 6371.0
N_MC = 1000

def to_xyz(lat, lon):
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    return np.array([np.cos(lat_r)*np.cos(lon_r), np.cos(lat_r)*np.sin(lon_r), np.sin(lat_r)])

def gc_points(pole_lat, pole_lon, n=3600):
    """Generate n evenly spaced points along a great circle defined by its pole."""
    pole = to_xyz(pole_lat, pole_lon)
    # Find two orthogonal vectors in the plane perpendicular to the pole
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

    lats = np.degrees(np.arcsin(points[:, 2]))
    lons = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    return lats, lons

def download_plate_boundaries():
    """Download Bird (2003) plate boundary data from GitHub."""
    cache_dir = os.path.expanduser("~/megalith_site_research/data/tectonic")
    os.makedirs(cache_dir, exist_ok=True)

    # Download the GeoJSON version
    url = "https://raw.githubusercontent.com/fraxen/tectonicplates/master/GeoJSON/PB2002_boundaries.json"
    filepath = os.path.join(cache_dir, "PB2002_boundaries.json")

    if not os.path.exists(filepath):
        print(f"Downloading plate boundaries from {url}...")
        urllib.request.urlretrieve(url, filepath)
        print(f"Saved to {filepath}")

    with open(filepath) as f:
        data = json.load(f)

    boundaries = []
    for feature in data['features']:
        btype = feature['properties'].get('Name', 'unknown')
        coords = feature['geometry']['coordinates']
        # LineString or MultiLineString
        if feature['geometry']['type'] == 'LineString':
            segments = [coords]
        else:
            segments = coords

        for seg in segments:
            lons = [c[0] for c in seg]
            lats = [c[1] for c in seg]
            boundaries.append({
                'type': btype,
                'lats': lats,
                'lons': lons
            })

    return boundaries

def count_boundary_crossings(gc_lats, gc_lons, boundaries, corridor_km=50):
    """Count how many plate boundary segments the great circle crosses (within corridor)."""
    crossings = set()
    crossing_types = []

    for bi, b in enumerate(boundaries):
        for j in range(len(b['lats'])):
            blat, blon = b['lats'][j], b['lons'][j]
            # Quick distance check to any GC point
            bxyz = to_xyz(blat, blon)

            # Check distance to nearest GC point (subsample for speed)
            for k in range(0, len(gc_lats), 10):  # check every 10th point
                gxyz = to_xyz(gc_lats[k], gc_lons[k])
                dot = np.clip(np.dot(bxyz, gxyz), -1, 1)
                dist = np.arccos(dot) * EARTH_R
                if dist < corridor_km:
                    if bi not in crossings:
                        crossings.add(bi)
                        # Classify boundary type
                        btype = b['type']
                        if '-' in btype:
                            # Extract boundary type from name like "EU-AF"
                            crossing_types.append('plate_boundary')
                        else:
                            crossing_types.append(btype)
                    break

    return len(crossings), crossing_types

def gc_distance_from_pole(pole_xyz, point_lat, point_lon):
    """Distance from a point to the great circle defined by a pole."""
    point_xyz = to_xyz(point_lat, point_lon)
    dot = np.clip(np.dot(pole_xyz, point_xyz), -1.0, 1.0)
    angular_dist = np.arccos(dot)
    return abs(angular_dist - np.pi/2) * EARTH_R

def count_crossings_fast(pole_lat, pole_lon, boundary_points, corridor_km=50):
    """Fast crossing count using pre-computed boundary point XYZ."""
    pole_xyz = to_xyz(pole_lat, pole_lon)

    # Vectorized distance computation
    dots = np.clip(boundary_points @ pole_xyz, -1, 1)
    angular_dists = np.arccos(dots)
    gc_dists = np.abs(angular_dists - np.pi/2) * EARTH_R

    return int(np.sum(gc_dists < corridor_km))

def main():
    print("=" * 60)
    print("TEST 4: TECTONIC DIVERSITY")
    print("=" * 60)

    # Download plate boundaries
    boundaries = download_plate_boundaries()
    print(f"Loaded {len(boundaries)} plate boundary segments")

    # Pre-compute all boundary point XYZ for fast MC
    all_points = []
    boundary_indices = []
    for bi, b in enumerate(boundaries):
        for j in range(len(b['lats'])):
            xyz = to_xyz(b['lats'][j], b['lons'][j])
            all_points.append(xyz)
            boundary_indices.append(bi)

    all_points = np.array(all_points)
    boundary_indices = np.array(boundary_indices)
    print(f"Total boundary vertices: {len(all_points)}")

    # Count unique boundaries within corridor
    def count_unique_boundaries(pole_lat, pole_lon, corridor_km=50):
        pole_xyz = to_xyz(pole_lat, pole_lon)
        dots = np.clip(all_points @ pole_xyz, -1, 1)
        angular_dists = np.arccos(dots)
        gc_dists = np.abs(angular_dists - np.pi/2) * EARTH_R
        near = gc_dists < corridor_km
        unique_boundaries = len(set(boundary_indices[near]))
        return unique_boundaries

    # Alison circle
    alison_crossings = count_unique_boundaries(*ALISON_POLE)
    print(f"\nAlison circle: {alison_crossings} plate boundaries within 50 km")

    # Also count at different corridor widths
    for ckm in [10, 25, 50, 100, 200]:
        n = count_unique_boundaries(*ALISON_POLE, corridor_km=ckm)
        print(f"  {ckm:4d} km corridor: {n} boundaries")

    # Compute boundary crossing density (per 1000 km of circle)
    circle_length = 2 * np.pi * EARTH_R  # ~40,030 km
    alison_density = alison_crossings / (circle_length / 1000)
    print(f"\nAlison crossing density: {alison_density:.2f} per 1000 km")

    # Monte Carlo
    print(f"\nRunning {N_MC} random circles...")
    rng = np.random.default_rng(42)
    mc_crossings = []

    for i in range(N_MC):
        z = rng.uniform(-1, 1)
        phi = rng.uniform(0, 2*np.pi)
        r = np.sqrt(1 - z**2)
        pole_lat = np.degrees(np.arcsin(z))
        pole_lon = np.degrees(phi) - 180

        n = count_unique_boundaries(pole_lat, pole_lon)
        mc_crossings.append(n)

        if (i+1) % 100 == 0:
            print(f"  {i+1}/{N_MC}")

    mc_crossings = np.array(mc_crossings)

    # Statistics
    mc_mean = np.mean(mc_crossings)
    mc_std = np.std(mc_crossings)
    mc_median = np.median(mc_crossings)
    percentile = np.mean(mc_crossings <= alison_crossings) * 100
    p_value = np.mean(mc_crossings >= alison_crossings)

    print(f"\nRandom circles: mean={mc_mean:.1f}, median={mc_median:.0f}, std={mc_std:.1f}")
    print(f"Alison: {alison_crossings} boundaries")
    print(f"Percentile: {percentile:.1f}%")
    print(f"P(random >= Alison): {p_value:.4f}")

    # Classify boundary types near Alison circle
    pole_xyz = to_xyz(*ALISON_POLE)
    dots = np.clip(all_points @ pole_xyz, -1, 1)
    angular_dists = np.arccos(dots)
    gc_dists = np.abs(angular_dists - np.pi/2) * EARTH_R
    near = gc_dists < 50
    near_boundaries = set(boundary_indices[near])

    boundary_names = []
    for bi in near_boundaries:
        b = boundaries[bi]
        boundary_names.append(b['type'])

    print(f"\nBoundary segments near Alison circle:")
    for name in sorted(boundary_names)[:20]:
        print(f"  {name}")

    # Results
    results = {
        "analysis": "Tectonic Diversity Test",
        "data": "Bird (2003) plate boundaries via GitHub",
        "corridor_km": 50,
        "alison": {
            "crossings": alison_crossings,
            "density_per_1000km": round(alison_density, 2),
            "boundary_names": sorted(boundary_names)[:30]
        },
        "mc": {
            "n_trials": N_MC,
            "mean": round(mc_mean, 1),
            "median": int(mc_median),
            "std": round(mc_std, 1),
            "p95": int(np.percentile(mc_crossings, 95)),
            "p99": int(np.percentile(mc_crossings, 99)),
            "max": int(np.max(mc_crossings))
        },
        "alison_percentile": round(percentile, 1),
        "p_value": round(p_value, 4),
        "significant": p_value < 0.05,
        "corridor_sensitivity": {},
        "conclusion": ""
    }

    # Corridor sensitivity
    for ckm in [10, 25, 50, 100, 200]:
        alison_n = count_unique_boundaries(*ALISON_POLE, corridor_km=ckm)
        rng2 = np.random.default_rng(42)
        mc_n = []
        for i in range(N_MC):
            z = rng2.uniform(-1, 1)
            phi = rng2.uniform(0, 2*np.pi)
            pole_lat = np.degrees(np.arcsin(z))
            pole_lon = np.degrees(phi) - 180
            mc_n.append(count_unique_boundaries(pole_lat, pole_lon, corridor_km=ckm))
        mc_n = np.array(mc_n)
        p = float(np.mean(mc_n >= alison_n))
        results["corridor_sensitivity"][str(ckm)] = {
            "alison": alison_n,
            "mc_mean": round(float(np.mean(mc_n)), 1),
            "p_value": round(p, 4)
        }

    if p_value < 0.05:
        results["conclusion"] = f"SIGNIFICANT: Alison circle crosses more plate boundaries ({alison_crossings}) than {100-percentile:.1f}% of random circles (p={p_value:.4f})"
    else:
        results["conclusion"] = f"NOT SIGNIFICANT: Alison circle crosses {alison_crossings} boundaries, which is typical (percentile={percentile:.1f}%, p={p_value:.4f})"

    out_path = os.path.expanduser("~/megalith_site_research/paper3/results/causal_mechanisms/tectonic_diversity_test.json")
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)

    # Markdown summary
    md = f"""# Test 4: Tectonic Diversity

## Result: {'SIGNIFICANT' if p_value < 0.05 else 'NOT SIGNIFICANT'}

The Alison Great Circle crosses **{alison_crossings} plate boundary segments** within 50 km
({alison_density:.2f} per 1,000 km of arc).

Random great circles cross a mean of {mc_mean:.1f} ± {mc_std:.1f} boundaries.

**Alison percentile: {percentile:.1f}%**
**p-value: {p_value:.4f}**

## Corridor Width Sensitivity
"""
    for ckm, vals in results["corridor_sensitivity"].items():
        md += f"- {ckm} km: Alison = {vals['alison']}, random mean = {vals['mc_mean']}, p = {vals['p_value']}\n"

    md_path = os.path.expanduser("~/megalith_site_research/paper3/results/causal_mechanisms/tectonic_diversity_test.md")
    with open(md_path, 'w') as f:
        f.write(md)

    print(f"\nSaved to {out_path}")
    print(results["conclusion"])

if __name__ == "__main__":
    main()
