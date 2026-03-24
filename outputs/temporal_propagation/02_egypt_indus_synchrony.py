#!/usr/bin/env python3
"""
Analysis 2: Egypt–Indus Synchrony Deep Dive
=============================================
Directive 04 — Temporal Propagation & Wavefront Analysis

Tests whether the overlap between Memphis pyramid-building (2750–2500 BCE)
and Mohenjo-daro/Harappa urban phase (2600–2300 BCE) is statistically
significant, and whether diffusion could explain the synchrony.
"""

import csv, json, math, os, sys
import numpy as np

# ============================================================
# CONFIGURATION
# ============================================================
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087
THRESHOLD_KM = 100.0  # Wider threshold for deep dive (directive says 100km)
N_MC = 10000

BASE_DIR = os.path.expanduser("~/megalith_site_research")
OUT_DIR = os.path.join(BASE_DIR, "The Great Circle", "outputs", "temporal_propagation")
os.makedirs(OUT_DIR, exist_ok=True)

# Region bounding boxes
EGYPT_BOX = {"lat": (27, 32), "lon": (29, 33)}  # Memphis/Nile Delta focus
INDUS_BOX = {"lat": (24, 32), "lon": (66, 76)}   # Indus Valley focus

# Known diffusion rates (km/year)
DIFFUSION_NEOLITHIC = 1.0
DIFFUSION_BRONZE_AGE = 5.0
DIFFUSION_MARITIME = 10.0  # Coastal/riverine trade


# ============================================================
# GEOMETRY
# ============================================================
def haversine(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return EARTH_R_KM * 2 * math.asin(math.sqrt(min(1.0, a)))


def gc_distance(lat, lon):
    d = haversine(POLE_LAT, POLE_LON, lat, lon)
    return abs(d - QUARTER_CIRC)


def bearing(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
    return (math.degrees(math.atan2(x, y)) + 360) % 360


def in_box(lat, lon, box):
    return box["lat"][0] <= lat <= box["lat"][1] and box["lon"][0] <= lon <= box["lon"][1]


# ============================================================
# DATA
# ============================================================
def bp_to_bce(age_bp):
    if age_bp <= 2000:
        return 1950 - age_bp
    elif age_bp <= 5000:
        cal_bp = age_bp * 1.05 - 100
    elif age_bp <= 10000:
        cal_bp = age_bp * 1.08 - 250
    else:
        cal_bp = age_bp * 1.10 - 450
    return 1950 - cal_bp


# Monumental keywords
MONUMENTAL_KW = ['temple', 'pyramid', 'tomb', 'mound', 'monument', 'cairn', 'megalith',
                  'henge', 'huaca', 'ceremonial', 'ritual', 'geoglyph', 'ahu', 'moai',
                  'stela', 'ziggurat', 'palace', 'citadel', 'fortress', 'dolmen',
                  'stone circle', 'standing stone', 'barrow', 'tumulus', 'sanctuary',
                  'necropolis', 'acropolis', 'mastaba']


# Pleiades monumental feature types
PLEIADES_MONUMENTAL = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "amphitheatre",
    "hippodrome", "aqueduct", "lighthouse", "arch", "theatre",
    "basilica", "acropolis", "nuraghe", "tumulus", "shrine",
}


def is_monumental(name="", feat_type=""):
    name_lower = (name or "").lower()
    feat_lower = (feat_type or "").lower()
    if any(kw in name_lower for kw in MONUMENTAL_KW):
        return True
    if any(kw in feat_lower for kw in PLEIADES_MONUMENTAL):
        return True
    return False


def load_egypt_sites():
    """Load all dated sites in the Egypt cluster, within 100km of circle."""
    sites = []

    # p3k14c
    p3k14c_path = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")
    with open(p3k14c_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["Lat"])
                lon = float(row["Long"])
                age_bp = float(row["Age"])
            except (ValueError, KeyError):
                continue
            if lat == 0 and lon == 0:
                continue
            if not in_box(lat, lon, EGYPT_BOX):
                continue
            dist = gc_distance(lat, lon)
            if dist > THRESHOLD_KM:
                continue
            cal_bce = bp_to_bce(age_bp)
            if cal_bce < -5000 or cal_bce > 500:
                continue
            name = row.get("SiteName", "")
            sites.append({
                "name": name or "Unnamed",
                "lat": lat, "lon": lon,
                "date_bce": round(cal_bce),
                "source": "p3k14c",
                "is_monumental": is_monumental(name),
            })

    # Pleiades
    pleiades_path = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")
    with open(pleiades_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["reprLat"])
                lon = float(row["reprLong"])
                min_date = float(row.get("minDate", ""))
            except (ValueError, KeyError):
                continue
            if not in_box(lat, lon, EGYPT_BOX):
                continue
            dist = gc_distance(lat, lon)
            if dist > THRESHOLD_KM:
                continue
            if min_date < -5000 or min_date > 500:
                continue
            feat = row.get("featureTypes", "")
            sites.append({
                "name": row.get("title", "Unnamed"),
                "lat": lat, "lon": lon,
                "date_bce": round(min_date),
                "source": "pleiades",
                "is_monumental": is_monumental(feat_type=feat),
            })

    # Megalithic Portal (known pyramids)
    portal_path = os.path.join(BASE_DIR, "data", "supplementary", "globe_sites_all.json")
    with open(portal_path) as f:
        all_globe = json.load(f)
    for s in all_globe:
        if s.get("source") != "portal":
            continue
        lat, lon = s["lat"], s["lon"]
        if not in_box(lat, lon, EGYPT_BOX):
            continue
        dist = gc_distance(lat, lon)
        if dist > THRESHOLD_KM:
            continue
        site_type = s.get("type", "")
        type_lower = site_type.lower()
        # Assign known dates for Egyptian pyramids
        if "pyramid" in type_lower or "mastaba" in type_lower:
            date_bce = -2600  # Peak Old Kingdom
        elif "temple" in type_lower:
            date_bce = -2500
        elif "tomb" in type_lower:
            date_bce = -2800
        else:
            continue
        sites.append({
            "name": s.get("name", "Unnamed"),
            "lat": lat, "lon": lon,
            "date_bce": date_bce,
            "source": "portal",
            "is_monumental": True,
        })

    print(f"Egypt cluster: {len(sites)} dated sites ({sum(1 for s in sites if s['is_monumental'])} monumental)")
    return sites


def load_indus_sites():
    """Load all dated sites in the Indus cluster, within 100km of circle."""
    sites = []

    # p3k14c
    p3k14c_path = os.path.join(BASE_DIR, "data", "p3k14c", "p3k14c_data.csv")
    with open(p3k14c_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["Lat"])
                lon = float(row["Long"])
                age_bp = float(row["Age"])
            except (ValueError, KeyError):
                continue
            if lat == 0 and lon == 0:
                continue
            if not in_box(lat, lon, INDUS_BOX):
                continue
            dist = gc_distance(lat, lon)
            if dist > THRESHOLD_KM:
                continue
            cal_bce = bp_to_bce(age_bp)
            if cal_bce < -5000 or cal_bce > 500:
                continue
            name = row.get("SiteName", "")
            sites.append({
                "name": name or "Unnamed",
                "lat": lat, "lon": lon,
                "date_bce": round(cal_bce),
                "source": "p3k14c",
                "is_monumental": is_monumental(name),
            })

    # Pleiades
    pleiades_path = os.path.join(BASE_DIR, "data", "pleiades", "pleiades-places-latest.csv")
    with open(pleiades_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["reprLat"])
                lon = float(row["reprLong"])
                min_date = float(row.get("minDate", ""))
            except (ValueError, KeyError):
                continue
            if not in_box(lat, lon, INDUS_BOX):
                continue
            dist = gc_distance(lat, lon)
            if dist > THRESHOLD_KM:
                continue
            if min_date < -5000 or min_date > 500:
                continue
            feat = row.get("featureTypes", "")
            sites.append({
                "name": row.get("title", "Unnamed"),
                "lat": lat, "lon": lon,
                "date_bce": round(min_date),
                "source": "pleiades",
                "is_monumental": is_monumental(feat_type=feat),
            })

    # Known Indus sites from literature (Possehl 2002, Kenoyer 1998)
    known_indus = [
        {"name": "Mohenjo-daro", "lat": 27.324, "lon": 68.138, "date_bce": -2600, "is_monumental": True},
        {"name": "Harappa", "lat": 30.631, "lon": 72.869, "date_bce": -2600, "is_monumental": True},
        {"name": "Dholavira", "lat": 23.886, "lon": 70.212, "date_bce": -2650, "is_monumental": True},
        {"name": "Ganweriwala", "lat": 29.307, "lon": 71.183, "date_bce": -2500, "is_monumental": True},
        {"name": "Rakhigarhi", "lat": 29.281, "lon": 76.116, "date_bce": -2600, "is_monumental": True},
        {"name": "Kalibangan", "lat": 29.470, "lon": 74.130, "date_bce": -2700, "is_monumental": True},
        {"name": "Lothal", "lat": 22.522, "lon": 72.250, "date_bce": -2400, "is_monumental": True},
        {"name": "Chanhu-daro", "lat": 26.153, "lon": 68.291, "date_bce": -2500, "is_monumental": True},
        {"name": "Banawali", "lat": 29.595, "lon": 75.384, "date_bce": -2500, "is_monumental": True},
        {"name": "Surkotada", "lat": 23.616, "lon": 70.833, "date_bce": -2300, "is_monumental": True},
    ]
    for site in known_indus:
        dist = gc_distance(site["lat"], site["lon"])
        if dist <= THRESHOLD_KM:
            sites.append({
                "name": site["name"],
                "lat": site["lat"], "lon": site["lon"],
                "date_bce": site["date_bce"],
                "source": "literature",
                "is_monumental": site["is_monumental"],
            })

    # South Asia merged dataset
    sa_path = os.path.join(BASE_DIR, "data", "south_asia")
    for fname in ["indus_sites_bates2019.csv", "south_asia_merged.csv"]:
        fpath = os.path.join(sa_path, fname)
        if not os.path.exists(fpath):
            continue
        with open(fpath) as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    lat = float(row.get("lat", row.get("latitude", "")))
                    lon = float(row.get("lon", row.get("longitude", "")))
                except (ValueError, KeyError):
                    continue
                if not in_box(lat, lon, INDUS_BOX):
                    continue
                dist = gc_distance(lat, lon)
                if dist > THRESHOLD_KM:
                    continue
                # Estimate date from period if available
                period = (row.get("period", "") or row.get("classification", "")).lower()
                if "mature" in period or "urban" in period:
                    date_bce = -2500
                elif "early" in period:
                    date_bce = -3300
                elif "late" in period:
                    date_bce = -2000
                else:
                    continue  # Skip undated
                name = row.get("name", row.get("site_name", "Unnamed"))
                sites.append({
                    "name": name or "Unnamed",
                    "lat": lat, "lon": lon,
                    "date_bce": date_bce,
                    "source": "south_asia_db",
                    "is_monumental": True,  # Indus cities are monumental/urban
                })

    print(f"Indus cluster: {len(sites)} dated sites ({sum(1 for s in sites if s['is_monumental'])} monumental)")
    return sites


# ============================================================
# ANALYSIS
# ============================================================
def compute_construction_rate(dates, bin_width=100):
    """Construction starts per century, returned as (bin_centers, rates)."""
    if not dates:
        return [], []
    min_d, max_d = min(dates), max(dates)
    bins = np.arange(min_d - bin_width, max_d + bin_width * 2, bin_width)
    counts, edges = np.histogram(dates, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    return centers.tolist(), counts.tolist()


def run_synchrony_analysis(egypt_sites, indus_sites):
    """Test Egypt–Indus temporal synchrony."""
    # Filter to monumental only
    egypt_mon = [s for s in egypt_sites if s["is_monumental"]]
    indus_mon = [s for s in indus_sites if s["is_monumental"]]

    egypt_dates = np.array([s["date_bce"] for s in egypt_mon])
    indus_dates = np.array([s["date_bce"] for s in indus_mon])

    print(f"\nEgypt monumental: n={len(egypt_mon)}, range={egypt_dates.min():.0f} to {egypt_dates.max():.0f}")
    print(f"Indus monumental: n={len(indus_mon)}, range={indus_dates.min():.0f} to {indus_dates.max():.0f}")

    # Peak construction period
    egypt_peak_bins, egypt_peak_rates = compute_construction_rate(egypt_dates.tolist())
    indus_peak_bins, indus_peak_rates = compute_construction_rate(indus_dates.tolist())

    egypt_peak_idx = np.argmax(egypt_peak_rates) if egypt_peak_rates else 0
    indus_peak_idx = np.argmax(indus_peak_rates) if indus_peak_rates else 0
    egypt_peak = egypt_peak_bins[egypt_peak_idx] if egypt_peak_bins else 0
    indus_peak = indus_peak_bins[indus_peak_idx] if indus_peak_bins else 0

    observed_offset = abs(egypt_peak - indus_peak)
    print(f"Egypt peak: ~{egypt_peak:.0f} BCE")
    print(f"Indus peak: ~{indus_peak:.0f} BCE")
    print(f"Observed offset: {observed_offset:.0f} years")

    # Monte Carlo: draw both clusters from uniform 4000–1000 BCE
    mc_offsets = np.zeros(N_MC)
    for i in range(N_MC):
        # Resample Egypt dates uniformly
        fake_egypt = np.random.uniform(-4000, -1000, len(egypt_mon))
        fake_indus = np.random.uniform(-4000, -1000, len(indus_mon))
        # Compute peaks
        e_bins, e_rates = compute_construction_rate(fake_egypt.tolist())
        i_bins, i_rates = compute_construction_rate(fake_indus.tolist())
        if e_rates and i_rates:
            e_peak = e_bins[np.argmax(e_rates)]
            i_peak = i_bins[np.argmax(i_rates)]
            mc_offsets[i] = abs(e_peak - i_peak)

    p_value = np.mean(mc_offsets <= observed_offset)
    print(f"Monte Carlo p-value (offsets ≤ observed): {p_value:.4f}")

    # Diffusion feasibility
    # Memphis (29.87°N, 31.25°E) to Mohenjo-daro (27.32°N, 68.14°E)
    memphis_lat, memphis_lon = 29.87, 31.25
    mohenjo_lat, mohenjo_lon = 27.324, 68.138
    direct_km = haversine(memphis_lat, memphis_lon, mohenjo_lat, mohenjo_lon)

    # Along the circle (via Iran) — approximate as 1.5x direct
    along_circle_km = direct_km * 1.5

    # Account for geographic barriers (Arabian desert)
    # Realistic route: Nile → Sinai → Levant → Mesopotamia → Iran → Indus
    # Estimated realistic route distance
    realistic_km = direct_km * 2.0  # Roughly double for overland detour

    diffusion = {
        "direct_distance_km": round(direct_km),
        "estimated_circle_route_km": round(along_circle_km),
        "estimated_realistic_route_km": round(realistic_km),
        "observed_delay_years": round(observed_offset),
        "neolithic_rate_1km_yr": {
            "expected_delay_direct": round(direct_km / DIFFUSION_NEOLITHIC),
            "expected_delay_realistic": round(realistic_km / DIFFUSION_NEOLITHIC),
            "feasible": observed_offset >= direct_km / DIFFUSION_NEOLITHIC * 0.5,
        },
        "bronze_age_rate_5km_yr": {
            "expected_delay_direct": round(direct_km / DIFFUSION_BRONZE_AGE),
            "expected_delay_realistic": round(realistic_km / DIFFUSION_BRONZE_AGE),
            "feasible": observed_offset >= direct_km / DIFFUSION_BRONZE_AGE * 0.5,
        },
        "maritime_rate_10km_yr": {
            "expected_delay_direct": round(direct_km / DIFFUSION_MARITIME),
            "expected_delay_realistic": round(realistic_km / DIFFUSION_MARITIME),
            "feasible": observed_offset >= direct_km / DIFFUSION_MARITIME * 0.5,
        },
        "note": "Arabian desert between Egypt and Indus makes direct overland diffusion unlikely. Coastal or Mesopotamian relay routes more plausible.",
    }

    synchrony_result = {
        "egypt": {
            "n_monumental": len(egypt_mon),
            "earliest": round(float(egypt_dates.min())),
            "latest": round(float(egypt_dates.max())),
            "mean": round(float(egypt_dates.mean())),
            "peak_century": round(egypt_peak),
        },
        "indus": {
            "n_monumental": len(indus_mon),
            "earliest": round(float(indus_dates.min())),
            "latest": round(float(indus_dates.max())),
            "mean": round(float(indus_dates.mean())),
            "peak_century": round(indus_peak),
        },
        "observed_offset_years": round(observed_offset),
        "monte_carlo": {
            "n_trials": N_MC,
            "p_value": round(p_value, 5),
            "null_mean_offset": round(float(mc_offsets.mean())),
            "null_median_offset": round(float(np.median(mc_offsets))),
            "interpretation": (
                f"The observed offset of {observed_offset:.0f} years between Egypt and Indus peaks "
                f"has p={p_value:.4f} under a uniform-date null. "
                + ("This is statistically significant (p<0.05) — the synchrony is unlikely by chance."
                   if p_value < 0.05 else
                   "This is NOT statistically significant — such overlap is expected by chance.")
            ),
        },
    }

    # --- Plot ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

        # Top: construction rate curves
        if egypt_peak_bins and indus_peak_bins:
            ax1.bar([b for b in egypt_peak_bins], egypt_peak_rates,
                   width=80, alpha=0.6, color="#E63946", label="Egypt")
            ax1.bar([b for b in indus_peak_bins], indus_peak_rates,
                   width=80, alpha=0.6, color="#2A9D8F", label="Indus")
            ax1.set_xlabel("Date (BCE)")
            ax1.set_ylabel("Construction starts per century")
            ax1.set_title("Egypt vs. Indus Construction Rate Timeline")
            ax1.legend()
            ax1.invert_xaxis()
            ax1.grid(True, alpha=0.3)

        # Bottom: Monte Carlo null distribution
        ax2.hist(mc_offsets, bins=50, alpha=0.7, color="#999999", edgecolor="black", linewidth=0.5)
        ax2.axvline(observed_offset, color="#E63946", linewidth=2,
                    label=f"Observed offset = {observed_offset:.0f} yr")
        ax2.set_xlabel("Peak offset (years)")
        ax2.set_ylabel("Count")
        ax2.set_title(f"Monte Carlo Null Distribution (n={N_MC}, p={p_value:.4f})")
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, "dual_construction_timeline.png"), dpi=150)
        plt.close()
        print("\nSaved dual_construction_timeline.png")
    except ImportError:
        print("\nMatplotlib not available — skipping plot")

    # Save outputs
    with open(os.path.join(OUT_DIR, "egypt_indus_synchrony.json"), "w") as f:
        json.dump(synchrony_result, f, indent=2)
    print("Saved egypt_indus_synchrony.json")

    with open(os.path.join(OUT_DIR, "diffusion_feasibility.json"), "w") as f:
        json.dump(diffusion, f, indent=2)
    print("Saved diffusion_feasibility.json")

    return synchrony_result, diffusion


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("EGYPT–INDUS SYNCHRONY DEEP DIVE")
    print("=" * 60)

    egypt = load_egypt_sites()
    indus = load_indus_sites()
    sync, diff = run_synchrony_analysis(egypt, indus)

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)
