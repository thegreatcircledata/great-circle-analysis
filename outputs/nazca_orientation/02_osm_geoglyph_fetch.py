#!/usr/bin/env python3
"""
Fetch Nazca geoglyph data from OpenStreetMap via Overpass API.
Extract polygon geometries and compute orientations from the long axis.
"""

import json
import math
import urllib.request
import urllib.parse
import time
import numpy as np
from collections import Counter

OUTPUT_DIR = "."

# ─── Overpass query ────────────────────────────────────────────────
# Search for archaeological sites, geoglyphs, and Nazca-related features
# within 60km of Nazca Pampa center

OVERPASS_QUERY = """
[out:json][timeout:120];
(
  // Archaeological sites near Nazca
  nwr["historic"="archaeological_site"](around:60000,-14.735,-75.130);
  // Any features tagged as geoglyphs
  nwr["site_type"="geoglyph"](around:60000,-14.735,-75.130);
  nwr["archaeological_site"="geoglyph"](around:60000,-14.735,-75.130);
  // Tourism attractions (some geoglyphs tagged this way)
  nwr["tourism"="attraction"]["name"~"[Nn]azca|[Nn]asca|[Gg]eoglyph|[Ll]ine"](around:60000,-14.735,-75.130);
);
out geom;
"""

print("=" * 65)
print("FETCHING NAZCA FEATURES FROM OPENSTREETMAP")
print("=" * 65)

url = "https://overpass-api.de/api/interpreter"
data = urllib.parse.urlencode({"data": OVERPASS_QUERY}).encode()

print("  Sending Overpass query...")
req = urllib.request.Request(url, data=data,
                             headers={"User-Agent": "NazcaResearch/1.0"})

try:
    with urllib.request.urlopen(req, timeout=180) as resp:
        result = json.loads(resp.read().decode())
except Exception as e:
    print(f"  ERROR: {e}")
    print("  Saving empty result.")
    with open("osm_nazca_raw.json", "w") as f:
        json.dump({"error": str(e), "elements": []}, f)
    raise SystemExit(1)

elements = result.get("elements", [])
print(f"  Retrieved {len(elements)} elements")

# Save raw data
with open("osm_nazca_raw.json", "w") as f:
    json.dump(result, f, indent=2)


# ─── Parse elements ───────────────────────────────────────────────

def extract_geometry(elem):
    """Extract coordinate list from an OSM element."""
    etype = elem["type"]
    if etype == "node":
        return [(elem["lat"], elem["lon"])]
    elif etype == "way":
        geom = elem.get("geometry", [])
        return [(p["lat"], p["lon"]) for p in geom]
    elif etype == "relation":
        # Extract all way members
        coords = []
        for member in elem.get("members", []):
            if member.get("geometry"):
                coords.extend([(p["lat"], p["lon"]) for p in member["geometry"]])
        return coords
    return []


def compute_orientation(coords):
    """
    Compute the orientation (long axis bearing) of a polygon/way using PCA.
    Returns bearing in degrees [0, 180) — axial data (no preferred direction).
    """
    if len(coords) < 3:
        return None

    # Convert to local Cartesian (km) centered on centroid
    lats = [c[0] for c in coords]
    lons = [c[1] for c in coords]
    clat = np.mean(lats)
    clon = np.mean(lons)

    # Approximate km offsets
    km_per_deg_lat = 111.32
    km_per_deg_lon = 111.32 * math.cos(math.radians(clat))

    x = [(lon - clon) * km_per_deg_lon for lon in lons]
    y = [(lat - clat) * km_per_deg_lat for lat in lats]

    points = np.column_stack([x, y])

    # PCA to find long axis
    cov = np.cov(points.T)
    if cov.shape != (2, 2):
        return None

    eigenvalues, eigenvectors = np.linalg.eigh(cov)

    # The eigenvector with the largest eigenvalue is the long axis
    long_axis = eigenvectors[:, np.argmax(eigenvalues)]

    # Convert to bearing (degrees from north, clockwise)
    # long_axis is (x, y) where x=east, y=north
    bearing = math.degrees(math.atan2(long_axis[0], long_axis[1])) % 180

    # Elongation ratio
    eig_sorted = sorted(eigenvalues, reverse=True)
    if eig_sorted[1] > 0:
        elongation = math.sqrt(eig_sorted[0] / eig_sorted[1])
    else:
        elongation = float('inf')

    return bearing, elongation


# ─── Process features ─────────────────────────────────────────────

print("\n" + "=" * 65)
print("PROCESSING FEATURES")
print("=" * 65)

features = []
tag_counts = Counter()

for elem in elements:
    tags = elem.get("tags", {})
    coords = extract_geometry(elem)

    name = tags.get("name", tags.get("name:en", "unnamed"))
    site_type = tags.get("site_type", tags.get("archaeological_site", "unknown"))

    tag_counts[site_type] += 1

    if len(coords) < 3:
        continue  # Need at least 3 points for orientation

    result = compute_orientation(coords)
    if result is None:
        continue

    bearing, elongation = result

    centroid_lat = np.mean([c[0] for c in coords])
    centroid_lon = np.mean([c[1] for c in coords])

    features.append({
        "osm_id": elem.get("id"),
        "osm_type": elem["type"],
        "name": name,
        "site_type": site_type,
        "tags": {k: v for k, v in tags.items() if k in
                 ["name", "name:en", "historic", "site_type",
                  "archaeological_site", "tourism", "description"]},
        "centroid_lat": round(centroid_lat, 6),
        "centroid_lon": round(centroid_lon, 6),
        "n_vertices": len(coords),
        "orientation_deg": round(bearing, 2),
        "elongation": round(elongation, 2),
    })

print(f"  Total elements: {len(elements)}")
print(f"  Features with computable orientation: {len(features)}")
print(f"\n  Site type breakdown:")
for st, count in tag_counts.most_common():
    print(f"    {st}: {count}")

# Filter for elongated features (elongation > 1.5 = clearly directional)
elongated = [f for f in features if f["elongation"] > 1.5]
print(f"\n  Elongated features (ratio > 1.5): {len(elongated)}")

# Print sample
print("\n  Sample features:")
for f in sorted(features, key=lambda x: -x["elongation"])[:20]:
    print(f"    {f['name'][:35]:35s}  type={f['site_type']:15s}  "
          f"bearing={f['orientation_deg']:6.1f}°  elong={f['elongation']:.1f}  "
          f"({f['centroid_lat']:.4f}, {f['centroid_lon']:.4f})")

# Save
with open("osm_geoglyph_orientations.json", "w") as f:
    json.dump({
        "source": "OpenStreetMap Overpass API",
        "query_center": {"lat": -14.735, "lon": -75.130},
        "query_radius_km": 60,
        "total_elements": len(elements),
        "features_with_orientation": len(features),
        "elongated_features": len(elongated),
        "features": features,
    }, f, indent=2)

print(f"\nSaved to osm_geoglyph_orientations.json")
