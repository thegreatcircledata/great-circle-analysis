#!/usr/bin/env python3
"""
Deep Dive Directive 3: High-Resolution Climate Reconstruction

Question: Does the monument-settlement divergence track climate more precisely
than the coarse Soreq Cave test (r = -0.77, p = 0.025 on 7-8 bins)?

Uses actual downloaded NOAA paleoclimate data:
  - Soreq Cave δ18O (Bar-Matthews 2003) — multi-centennial
  - Grant et al. 2012 updated Soreq — higher resolution, Bayesian chronology
  - LC21 Mediterranean marine core (Grant 2012) — sub-centennial Holocene
  - Jeita Cave Lebanon (Cheng 2016) — high-res, 20 ka

Tests:
  - Cross-correlation at 100yr, 250yr, 500yr, 1000yr bin resolutions
  - Lead-lag analysis (-500yr to +500yr)
  - Multi-proxy comparison
  - 4.2 ka event timing precision
"""

import json
import math
import os
import sys
import csv
import numpy as np
from scipy import stats as sp_stats
from scipy.interpolate import interp1d
from pathlib import Path

# ── Constants ───────────────────────────────────────────────────────────────
EARTH_R_KM = 6371.0
QUARTER_CIRC = EARTH_R_KM * math.pi / 2
POLE_LAT = 59.682122
POLE_LON = -138.646087

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
OUT_DIR = PROJECT_ROOT / "outputs" / "deep_dive_tests"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PALEO_DIR = PROJECT_ROOT / "data" / "paleoclimate"
PLEIADES_PATH = PROJECT_ROOT / "data" / "pleiades" / "pleiades-places-latest.csv"
DIVERGENCE_PATH = PROJECT_ROOT / "results" / "temporal_divergence_bins.json"

# ── Pleiades site classification ────────────────────────────────────────────
MONUMENTAL_TYPES = {
    "pyramid", "temple", "temple-2", "sanctuary", "monument",
    "architecturalcomplex", "tomb", "church", "church-2",
    "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
    "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus", "shrine"
}
SETTLEMENT_TYPES = {
    "settlement", "settlement-modern", "village", "villa", "farmstead",
    "station", "findspot", "port", "mine", "mine-2", "quarry",
    "bath", "cistern", "bridge", "fortified-settlement", "townhouse", "production"
}

# ── Geometry helpers ────────────────────────────────────────────────────────
def haversine_km(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi/2)**2 + math.cos(phi1)*math.cos(phi2)*math.sin(dlam/2)**2
    return 2 * EARTH_R_KM * math.atan2(math.sqrt(a), math.sqrt(1-a))

def dist_from_gc(lat, lon):
    d = haversine_km(lat, lon, POLE_LAT, POLE_LON)
    return abs(d - QUARTER_CIRC)

# ── Load paleoclimate data ──────────────────────────────────────────────────

def load_noaa_file(filepath, skip_comments=True):
    """Load NOAA tab-delimited paleoclimate file, skipping comment lines."""
    data = []
    headers = None
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            if headers is None:
                headers = line.split("\t")
                continue
            parts = line.split("\t")
            row = {}
            for i, h in enumerate(headers):
                if i < len(parts):
                    val = parts[i].strip()
                    try:
                        row[h] = float(val)
                    except ValueError:
                        row[h] = val if val else None
                else:
                    row[h] = None
            data.append(row)
    return headers, data


def load_soreq_2003():
    """Load Bar-Matthews 2003 Soreq Cave data."""
    fp = PALEO_DIR / "soreq_2003-noaa.txt"
    if not fp.exists():
        print(f"  WARNING: {fp} not found")
        return None
    headers, data = load_noaa_file(fp)
    print(f"  Soreq 2003: {len(data)} rows, columns: {headers}")

    # Extract age (ka BP) and d18O
    ages_ka = []
    d18o = []
    for row in data:
        age = row.get("age_calkaBP")
        val = row.get("d18O_vpdb")
        if age is not None and val is not None:
            ages_ka.append(float(age))
            d18o.append(float(val))

    return {
        "name": "Soreq Cave (Bar-Matthews 2003)",
        "ages_ka_bp": np.array(ages_ka),
        "d18O": np.array(d18o),
        "proxy_type": "speleothem_d18O",
        "location": "31.45°N, 35.03°E",
        "interpretation": "More negative = wetter conditions",
    }


def load_grant2012_soreq():
    """Load Grant 2012 updated Soreq record (rescaled chronology)."""
    fp = PALEO_DIR / "grant2012soreq-noaa.txt"
    if not fp.exists():
        print(f"  WARNING: {fp} not found")
        return None
    headers, data = load_noaa_file(fp)
    print(f"  Grant 2012 Soreq: {len(data)} rows, columns: {headers}")

    ages_ka = []
    d18o = []
    for row in data:
        age = row.get("age_calkaBP_rescale")
        val = row.get("d18OcarbVPDB_rescale")
        if age is not None and val is not None:
            try:
                a = float(age)
                v = float(val)
                ages_ka.append(a)
                d18o.append(v)
            except (ValueError, TypeError):
                continue

    return {
        "name": "Soreq Cave (Grant 2012 rescaled)",
        "ages_ka_bp": np.array(ages_ka),
        "d18O": np.array(d18o),
        "proxy_type": "speleothem_d18O",
        "location": "31.45°N, 35.03°E",
        "interpretation": "More negative = wetter conditions, Bayesian-modelled chronology",
    }


def load_lc21():
    """Load LC21 marine core G. ruber d18O (Grant 2012)."""
    fp = PALEO_DIR / "grant2012-ruber-noaa.txt"
    if not fp.exists():
        print(f"  WARNING: {fp} not found")
        return None
    headers, data = load_noaa_file(fp)
    print(f"  LC21 G.ruber: {len(data)} rows, columns: {headers}")

    ages_ka = []
    d18o = []
    for row in data:
        age = row.get("age_calkaBP")
        # Find d18O column - might have different exact name
        val = None
        for key in row:
            if "d18O" in key or "rub" in key.lower():
                val = row[key]
                break
        if age is not None and val is not None:
            try:
                a = float(age)
                v = float(val)
                ages_ka.append(a)
                d18o.append(v)
            except (ValueError, TypeError):
                continue

    return {
        "name": "LC21 Core G. ruber (Grant 2012)",
        "ages_ka_bp": np.array(ages_ka),
        "d18O": np.array(d18o),
        "proxy_type": "marine_foram_d18O",
        "location": "35.667°N, 26.583°E (Aegean Sea)",
        "interpretation": "More negative = warmer SST / lighter ice",
    }


def load_jeita():
    """Load Jeita Cave Lebanon speleothem data."""
    fp = PALEO_DIR / "jeita2015iso-j3-noaa.txt"
    if not fp.exists():
        print(f"  WARNING: {fp} not found")
        return None
    headers, data = load_noaa_file(fp)
    print(f"  Jeita Cave: {len(data)} rows, columns: {headers}")

    ages_ka = []
    d18o = []
    for row in data:
        age = row.get("age_calBP")
        val = row.get("d18OcarbVPDB")
        if age is not None and val is not None:
            try:
                a = float(age) / 1000  # convert yr BP to ka BP
                v = float(val)
                ages_ka.append(a)
                d18o.append(v)
            except (ValueError, TypeError):
                continue

    return {
        "name": "Jeita Cave J3 (Cheng 2016)",
        "ages_ka_bp": np.array(ages_ka),
        "d18O": np.array(d18o),
        "proxy_type": "speleothem_d18O",
        "location": "33.95°N, 35.65°E (Lebanon)",
        "interpretation": "More negative = wetter conditions",
    }


# ── Compute monument-settlement divergence at variable resolution ──────────

def load_pleiades_sites():
    """Load Pleiades sites with type and date classification."""
    sites = []
    with open(PLEIADES_PATH, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row.get("reprLat", ""))
                lon = float(row.get("reprLong", ""))
            except ValueError:
                continue

            ft_str = row.get("featureTypes", "")
            feature_types = {t.strip() for t in ft_str.split(",") if t.strip()}

            try:
                min_date = int(row.get("minDate", ""))
                max_date = int(row.get("maxDate", ""))
            except ValueError:
                continue

            is_monumental = bool(feature_types & MONUMENTAL_TYPES)
            is_settlement = bool(feature_types & SETTLEMENT_TYPES)

            if not (is_monumental or is_settlement):
                continue

            gc_dist = dist_from_gc(lat, lon)

            sites.append({
                "lat": lat, "lon": lon,
                "min_date": min_date, "max_date": max_date,
                "is_monumental": is_monumental,
                "is_settlement": is_settlement,
                "gc_dist_km": gc_dist,
            })

    print(f"  Loaded {len(sites)} classified Pleiades sites")
    return sites


def compute_divergence_at_resolution(sites, bin_width_years, gc_threshold_km=50,
                                       time_range=(-7000, 2000)):
    """Compute monument-settlement divergence D at specified temporal resolution.

    For each time bin:
      - Count monuments and settlements within gc_threshold_km of the GC
      - Compute enrichment ratio for each
      - D = (monument_fraction_on_GC) - (settlement_fraction_on_GC)
    """
    bin_start = time_range[0]
    bin_end = time_range[1]

    bins = []
    t = bin_start
    while t < bin_end:
        t_end = t + bin_width_years

        # Sites active in this bin: minDate <= t_end and maxDate >= t
        mon_total = 0
        mon_on_gc = 0
        set_total = 0
        set_on_gc = 0

        for s in sites:
            if s["max_date"] < t or s["min_date"] > t_end:
                continue  # not active in this bin

            on_gc = s["gc_dist_km"] < gc_threshold_km

            if s["is_monumental"]:
                mon_total += 1
                if on_gc:
                    mon_on_gc += 1
            if s["is_settlement"]:
                set_total += 1
                if on_gc:
                    set_on_gc += 1

        # Compute fractions
        mon_frac = mon_on_gc / mon_total if mon_total > 0 else 0
        set_frac = set_on_gc / set_total if set_total > 0 else 0
        D = mon_frac - set_frac if (mon_total >= 10 and set_total >= 10) else None

        bins.append({
            "bin_start": t,
            "bin_end": t_end,
            "bin_mid": (t + t_end) / 2,
            "bin_mid_ka_bp": -(t + t_end) / 2000 + 1.950,  # convert CE year to ka BP (ref 1950)
            "mon_total": mon_total,
            "mon_on_gc": mon_on_gc,
            "mon_frac": mon_frac,
            "set_total": set_total,
            "set_on_gc": set_on_gc,
            "set_frac": set_frac,
            "divergence_D": D,
        })

        t = t_end

    return bins


# ── Cross-correlation analysis ──────────────────────────────────────────────

def bin_climate_data(climate_record, bin_width_ka, time_range_ka=(0, 9)):
    """Bin climate proxy data into uniform bins."""
    if climate_record is None:
        return None

    ages = climate_record["ages_ka_bp"]
    values = climate_record["d18O"]

    # Sort by age
    order = np.argsort(ages)
    ages = ages[order]
    values = values[order]

    # Filter to time range
    mask = (ages >= time_range_ka[0]) & (ages <= time_range_ka[1])
    ages = ages[mask]
    values = values[mask]

    if len(ages) < 5:
        return None

    # Create interpolation function
    interp_func = interp1d(ages, values, kind='linear', fill_value='extrapolate')

    # Bin
    bins = []
    t = time_range_ka[0]
    while t < time_range_ka[1]:
        t_end = t + bin_width_ka
        # Average all data points in this bin
        in_bin = (ages >= t) & (ages < t_end)
        if np.sum(in_bin) > 0:
            bin_val = np.mean(values[in_bin])
        else:
            # Interpolate
            bin_val = float(interp_func((t + t_end) / 2))

        bins.append({
            "bin_start_ka": t,
            "bin_end_ka": t_end,
            "bin_mid_ka": (t + t_end) / 2,
            "value": bin_val,
            "n_points": int(np.sum(in_bin)),
        })
        t = t_end

    return bins


def correlate_divergence_climate(div_bins, climate_bins, lag_ka=0):
    """Correlate divergence D with climate proxy at a given lag.

    lag_ka > 0 means climate LEADS divergence (climate change happens first).
    """
    if climate_bins is None:
        return None

    # Match bins by midpoint time
    climate_dict = {round(b["bin_mid_ka"], 3): b["value"] for b in climate_bins}

    d_vals = []
    c_vals = []
    mids = []

    for db in div_bins:
        if db["divergence_D"] is None:
            continue
        # Convert divergence bin mid to ka BP
        div_ka = db["bin_mid_ka_bp"]
        # Look up climate at (div_ka + lag_ka) -- climate leads means we need earlier climate
        climate_ka = div_ka + lag_ka

        # Find closest climate bin
        best_match = None
        best_dist = float('inf')
        for ck, cv in climate_dict.items():
            dist = abs(ck - climate_ka)
            if dist < best_dist:
                best_dist = dist
                best_match = cv

        # Only match if within half a bin width
        bin_width_ka = climate_bins[0]["bin_end_ka"] - climate_bins[0]["bin_start_ka"]
        if best_dist <= bin_width_ka * 0.75 and best_match is not None:
            d_vals.append(db["divergence_D"])
            c_vals.append(best_match)
            mids.append(div_ka)

    if len(d_vals) < 4:
        return {
            "lag_ka": lag_ka,
            "n_matched": len(d_vals),
            "pearson_r": None,
            "pearson_p": None,
            "spearman_rho": None,
            "spearman_p": None,
        }

    r, p = sp_stats.pearsonr(d_vals, c_vals)
    rho, p_rho = sp_stats.spearmanr(d_vals, c_vals)

    return {
        "lag_ka": lag_ka,
        "n_matched": len(d_vals),
        "pearson_r": round(float(r), 4),
        "pearson_p": round(float(p), 6),
        "spearman_rho": round(float(rho), 4),
        "spearman_p": round(float(p_rho), 6),
        "matched_ka": [round(m, 3) for m in mids],
        "D_values": [round(d, 4) for d in d_vals],
        "climate_values": [round(c, 4) for c in c_vals],
    }


# ══════════════════════════════════════════════════════════════════════════
# MAIN ANALYSIS
# ══════════════════════════════════════════════════════════════════════════

def main():
    print(f"\n{'='*70}")
    print(f"DIRECTIVE 3: HIGH-RESOLUTION CLIMATE RECONSTRUCTION")
    print(f"{'='*70}")

    results = {
        "meta": {
            "analysis": "Deep Dive Directive 3: High-Resolution Paleoclimate Correlation",
            "date": "2026-03-21",
            "pole": {"lat": POLE_LAT, "lon": POLE_LON},
        },
    }

    # ── Step 1: Load climate data ──────────────────────────────────────
    print(f"\n{'─'*50}")
    print("Step 1: Loading paleoclimate data")
    print(f"{'─'*50}")

    climate_records = {}
    for loader, key in [
        (load_soreq_2003, "soreq_2003"),
        (load_grant2012_soreq, "grant2012_soreq"),
        (load_lc21, "lc21"),
        (load_jeita, "jeita"),
    ]:
        record = loader()
        if record is not None:
            climate_records[key] = record
            # Count Holocene points (0-12 ka BP)
            holocene_mask = record["ages_ka_bp"] <= 12
            n_hol = np.sum(holocene_mask)
            age_range = (np.min(record["ages_ka_bp"]), np.max(record["ages_ka_bp"]))
            print(f"    {record['name']}: {len(record['ages_ka_bp'])} total points, "
                  f"{n_hol} in Holocene, range {age_range[0]:.2f}-{age_range[1]:.2f} ka BP")

    results["climate_data_summary"] = {
        key: {
            "name": rec["name"],
            "n_total": len(rec["ages_ka_bp"]),
            "n_holocene": int(np.sum(rec["ages_ka_bp"] <= 12)),
            "age_range_ka": [round(float(np.min(rec["ages_ka_bp"])), 3),
                             round(float(np.max(rec["ages_ka_bp"])), 3)],
        }
        for key, rec in climate_records.items()
    }

    # ── Step 2: Compute divergence at multiple resolutions ─────────────
    print(f"\n{'─'*50}")
    print("Step 2: Computing monument-settlement divergence")
    print(f"{'─'*50}")

    print("  Loading Pleiades sites...")
    sites = load_pleiades_sites()

    resolutions = [100, 250, 500, 1000]
    divergence_series = {}

    for res in resolutions:
        bins = compute_divergence_at_resolution(sites, res, gc_threshold_km=50,
                                                  time_range=(-7000, 2000))
        valid_bins = [b for b in bins if b["divergence_D"] is not None]
        divergence_series[res] = bins
        print(f"  Resolution {res}yr: {len(bins)} total bins, {len(valid_bins)} with valid D")

    # ── Step 3: Multi-scale cross-correlation ──────────────────────────
    print(f"\n{'─'*50}")
    print("Step 3: Multi-scale cross-correlation")
    print(f"{'─'*50}")

    multiscale_results = {}

    for res in resolutions:
        bin_width_ka = res / 1000.0
        div_bins = divergence_series[res]

        multiscale_results[f"{res}yr"] = {}

        for clim_key, clim_rec in climate_records.items():
            # Bin climate at same resolution
            clim_bins = bin_climate_data(clim_rec, bin_width_ka, time_range_ka=(0, 9))
            if clim_bins is None:
                continue

            # Lag-0 correlation
            corr = correlate_divergence_climate(div_bins, clim_bins, lag_ka=0)
            if corr is None:
                continue

            multiscale_results[f"{res}yr"][clim_key] = {
                "lag_0": corr,
            }

            if corr["pearson_r"] is not None:
                sig = "***" if corr["pearson_p"] < 0.001 else "**" if corr["pearson_p"] < 0.01 else "*" if corr["pearson_p"] < 0.05 else "ns"
                print(f"  {res}yr x {clim_rec['name'][:25]:25s}: r={corr['pearson_r']:.4f} p={corr['pearson_p']:.6f} "
                      f"rho={corr['spearman_rho']:.4f} n={corr['n_matched']} {sig}")
            else:
                print(f"  {res}yr x {clim_rec['name'][:25]:25s}: insufficient data (n={corr['n_matched']})")

    results["multiscale_correlations"] = multiscale_results

    # ── Step 4: Lead-lag analysis at best resolution ──────────────────
    print(f"\n{'─'*50}")
    print("Step 4: Lead-lag analysis")
    print(f"{'─'*50}")

    # Use 500yr bins for lead-lag (good balance of resolution and sample size)
    lead_lag_results = {}
    lags_ka = np.arange(-0.5, 0.51, 0.1)  # -500yr to +500yr in 100yr steps

    for clim_key, clim_rec in climate_records.items():
        clim_bins = bin_climate_data(clim_rec, 0.5, time_range_ka=(0, 9))
        if clim_bins is None:
            continue

        lag_series = []
        for lag in lags_ka:
            corr = correlate_divergence_climate(divergence_series[500], clim_bins, lag_ka=lag)
            if corr and corr["pearson_r"] is not None:
                lag_series.append(corr)

        if lag_series:
            # Find optimal lag
            best = max(lag_series, key=lambda x: abs(x["pearson_r"]) if x["pearson_r"] is not None else 0)
            lead_lag_results[clim_key] = {
                "all_lags": lag_series,
                "optimal_lag_ka": best["lag_ka"],
                "optimal_r": best["pearson_r"],
                "optimal_p": best["pearson_p"],
                "optimal_n": best["n_matched"],
            }
            print(f"  {clim_rec['name'][:30]:30s}: optimal lag = {best['lag_ka']*1000:.0f}yr, "
                  f"r={best['pearson_r']:.4f}, p={best['pearson_p']:.6f}")

            # Print full lag series
            for ls in lag_series:
                marker = " <-- OPTIMAL" if ls["lag_ka"] == best["lag_ka"] else ""
                if ls["pearson_r"] is not None:
                    print(f"    lag={ls['lag_ka']*1000:+5.0f}yr: r={ls['pearson_r']:.4f} (p={ls['pearson_p']:.4f}) n={ls['n_matched']}{marker}")

    results["lead_lag"] = lead_lag_results

    # ── Step 5: 4.2 ka event analysis ─────────────────────────────────
    print(f"\n{'─'*50}")
    print("Step 5: 4.2 ka event timing")
    print(f"{'─'*50}")

    # At 100yr resolution, what happens around 4.2 ka BP (2250 BCE)?
    event_results = {}

    for clim_key, clim_rec in climate_records.items():
        ages = clim_rec["ages_ka_bp"]
        d18o = clim_rec["d18O"]

        # Focus on 3.5-5.0 ka BP window (1550-3050 BCE)
        mask = (ages >= 3.5) & (ages <= 5.5)
        if np.sum(mask) < 3:
            continue

        sub_ages = ages[mask]
        sub_d18o = d18o[mask]

        # Sort by age
        order = np.argsort(sub_ages)
        sub_ages = sub_ages[order]
        sub_d18o = sub_d18o[order]

        # Find the d18O shift — look for maximum rate of change
        if len(sub_ages) >= 5:
            # Compute rate of change (d(d18O)/d(age))
            rates = np.diff(sub_d18o) / np.diff(sub_ages)
            # The 4.2ka drying should show up as a positive rate (d18O becomes less negative)
            max_rate_idx = np.argmax(rates)
            event_age = (sub_ages[max_rate_idx] + sub_ages[max_rate_idx + 1]) / 2

            event_results[clim_key] = {
                "name": clim_rec["name"],
                "n_points_3.5_5.5ka": int(np.sum(mask)),
                "max_drying_rate_age_ka": round(float(event_age), 3),
                "max_drying_rate_age_bce": round(float(event_age * 1000 - 1950), 0),
                "d18o_before_event": round(float(np.mean(sub_d18o[sub_ages < event_age])), 3),
                "d18o_after_event": round(float(np.mean(sub_d18o[sub_ages > event_age])), 3) if np.sum(sub_ages > event_age) > 0 else None,
            }
            print(f"  {clim_rec['name'][:30]:30s}: max drying rate at {event_age:.3f} ka BP "
                  f"({event_age*1000-1950:.0f} BCE), "
                  f"n={int(np.sum(mask))} points in window")

    # Divergence collapse timing
    div_100 = divergence_series[100]
    event_div_bins = [b for b in div_100 if b["divergence_D"] is not None
                      and -3500 <= b["bin_mid"] <= -1500]
    if event_div_bins:
        # Find the bin with largest D drop
        d_values = [(b["bin_mid"], b["divergence_D"]) for b in event_div_bins]
        for i in range(1, len(d_values)):
            d_drop = d_values[i-1][1] - d_values[i][1]
            if d_drop > 0:
                print(f"  Divergence drop: {d_values[i-1][0]:.0f} to {d_values[i][0]:.0f} CE: "
                      f"D {d_values[i-1][1]:.3f} -> {d_values[i][1]:.3f} (Δ={d_drop:.3f})")

    results["event_4_2ka"] = event_results

    # ── Step 6: Multi-proxy comparison ────────────────────────────────
    print(f"\n{'─'*50}")
    print("Step 6: Multi-proxy consistency check")
    print(f"{'─'*50}")

    # At 500yr resolution, do all proxies show the same correlation pattern?
    consistency = {}
    for clim_key in climate_records:
        entry = multiscale_results.get("500yr", {}).get(clim_key, {}).get("lag_0", {})
        if entry and entry.get("pearson_r") is not None:
            consistency[clim_key] = {
                "r": entry["pearson_r"],
                "p": entry["pearson_p"],
                "n": entry["n_matched"],
            }
            print(f"  {climate_records[clim_key]['name'][:35]:35s}: r={entry['pearson_r']:.4f} (p={entry['pearson_p']:.4f})")

    results["multiproxy_consistency"] = consistency

    # Count how many proxies show significant (p<0.05) negative correlation
    sig_negative = sum(1 for v in consistency.values() if v["r"] < 0 and v["p"] < 0.05)
    total_proxies = len(consistency)
    print(f"\n  {sig_negative}/{total_proxies} proxies show significant negative correlation "
          f"(more negative d18O = wetter ↔ higher divergence)")

    # ── Interpretation ──────────────────────────────────────────────
    print(f"\n{'='*70}")
    print("INTERPRETATION")
    print(f"{'='*70}")

    interpretation = {
        "key_question": "Does the climate-divergence correlation survive at high resolution?",
        "findings": [],
        "caveats": [
            "Divergence D at fine resolution has small sample sizes per bin",
            "Pleiades dating precision (minDate/maxDate) limits temporal resolution",
            "Speleothem and marine d18O respond to different climate forcings",
            "The 4.2 ka event timing varies by proxy and location",
            "Correlation ≠ causation — climate and construction may share common drivers",
            "BH correction should be applied for multiple resolution tests",
        ],
    }

    # Summarize key findings
    for res in resolutions:
        res_data = multiscale_results.get(f"{res}yr", {})
        for clim_key, data in res_data.items():
            lag0 = data.get("lag_0", {})
            r = lag0.get("pearson_r")
            p = lag0.get("pearson_p")
            n = lag0.get("n_matched", 0)
            if r is not None:
                interpretation["findings"].append(
                    f"{res}yr bins, {climate_records[clim_key]['name'][:25]}: r={r:.4f}, p={p:.6f}, n={n}"
                )

    results["interpretation"] = interpretation

    for f in interpretation["findings"]:
        print(f"  • {f}")
    print("\n  Caveats:")
    for c in interpretation["caveats"]:
        print(f"    – {c}")

    # ── Generate plots ──────────────────────────────────────────────
    print(f"\n{'─'*50}")
    print("Generating plots...")
    try:
        make_plots(climate_records, divergence_series, multiscale_results, lead_lag_results, results)
    except Exception as e:
        print(f"  Plot error: {e}")
        import traceback
        traceback.print_exc()

    # Save results
    def jsonify(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(f"Not JSON serializable: {type(obj)}")

    with open(OUT_DIR / "climate_highres.json", "w") as f:
        json.dump(results, f, indent=2, default=jsonify)
    print(f"\n  Saved: {OUT_DIR / 'climate_highres.json'}")

    return results


def make_plots(climate_records, divergence_series, multiscale_results, lead_lag_results, results):
    """Generate all climate analysis plots."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # ── Plot 1: Multi-scale correlation summary ──
    fig, ax = plt.subplots(figsize=(10, 6))

    resolutions = [100, 250, 500, 1000]
    colors = {"soreq_2003": "green", "grant2012_soreq": "blue",
              "lc21": "red", "jeita": "purple"}

    for clim_key in climate_records:
        rs = []
        ps = []
        valid_res = []
        for res in resolutions:
            entry = multiscale_results.get(f"{res}yr", {}).get(clim_key, {}).get("lag_0", {})
            if entry and entry.get("pearson_r") is not None:
                rs.append(entry["pearson_r"])
                ps.append(entry["pearson_p"])
                valid_res.append(res)

        if valid_res:
            label = climate_records[clim_key]["name"][:25]
            ax.plot(valid_res, rs, 'o-', color=colors.get(clim_key, "gray"),
                   label=label, linewidth=2, markersize=8)

            # Mark significant points
            for r_val, p_val, res_val in zip(rs, ps, valid_res):
                if p_val < 0.05:
                    ax.plot(res_val, r_val, 'o', color=colors.get(clim_key, "gray"),
                           markersize=14, markeredgewidth=2, markeredgecolor='black',
                           markerfacecolor='none')

    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel("Bin resolution (years)")
    ax.set_ylabel("Pearson r (D vs δ18O)")
    ax.set_title("Climate-Divergence Correlation at Multiple Time Scales")
    ax.legend(fontsize=8, loc='best')
    ax.set_xscale('log')
    ax.set_xticks(resolutions)
    ax.set_xticklabels([str(r) for r in resolutions])
    ax.annotate("● = significant (p<0.05)", xy=(0.02, 0.02), xycoords='axes fraction',
                fontsize=8, style='italic')

    plt.tight_layout()
    plt.savefig(OUT_DIR / "climate_multiscale.png", dpi=150, bbox_inches="tight")
    print(f"  Saved: {OUT_DIR / 'climate_multiscale.png'}")
    plt.close()

    # ── Plot 2: Lead-lag analysis ──
    fig, axes = plt.subplots(1, len(lead_lag_results), figsize=(5*len(lead_lag_results), 5),
                              squeeze=False)
    axes = axes[0]

    for i, (clim_key, lag_data) in enumerate(lead_lag_results.items()):
        ax = axes[i]
        lags = [l["lag_ka"] * 1000 for l in lag_data["all_lags"]]  # convert to years
        rs = [l["pearson_r"] for l in lag_data["all_lags"]]

        ax.plot(lags, rs, 'o-', color=colors.get(clim_key, "gray"), linewidth=2, markersize=6)
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.axvline(x=0, color='gray', linestyle=':', alpha=0.5)

        # Mark optimal
        opt_lag = lag_data["optimal_lag_ka"] * 1000
        opt_r = lag_data["optimal_r"]
        ax.plot(opt_lag, opt_r, '*', color='red', markersize=15, zorder=5)

        ax.set_xlabel("Lag (years, + = climate leads)")
        ax.set_ylabel("Pearson r")
        ax.set_title(f"{climate_records[clim_key]['name'][:25]}\n"
                    f"Optimal: {opt_lag:.0f}yr, r={opt_r:.3f}")

    plt.suptitle("Lead-Lag Analysis: Climate vs Divergence", fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUT_DIR / "climate_leadlag.png", dpi=150, bbox_inches="tight")
    print(f"  Saved: {OUT_DIR / 'climate_leadlag.png'}")
    plt.close()

    # ── Plot 3: 4.2 ka event zoom ──
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    # Climate proxies around 4.2 ka
    ax = axes[0]
    for clim_key, clim_rec in climate_records.items():
        ages = clim_rec["ages_ka_bp"]
        d18o = clim_rec["d18O"]
        mask = (ages >= 3.0) & (ages <= 6.0)
        if np.sum(mask) < 3:
            continue
        sub_ages = ages[mask]
        sub_d18o = d18o[mask]
        order = np.argsort(sub_ages)
        ax.plot(sub_ages[order], sub_d18o[order], '.-',
               color=colors.get(clim_key, "gray"),
               label=clim_rec["name"][:25], alpha=0.8, markersize=3)

    ax.axvline(x=4.2, color='orange', linestyle='--', linewidth=2, label='4.2 ka event')
    ax.set_ylabel("δ18O (‰ VPDB)")
    ax.set_title("Climate Proxies Around the 4.2 ka Event")
    ax.legend(fontsize=7, loc='best')
    ax.invert_xaxis()

    # Divergence around 4.2 ka
    ax = axes[1]
    for res in [100, 250, 500]:
        bins = divergence_series[res]
        valid = [b for b in bins if b["divergence_D"] is not None]
        # Filter to 3-6 ka BP window (= approx -4050 to -1050 BCE)
        event_bins = [b for b in valid if 1.0 <= b["bin_mid_ka_bp"] <= 6.0]
        if event_bins:
            mids = [b["bin_mid_ka_bp"] for b in event_bins]
            ds = [b["divergence_D"] for b in event_bins]
            ax.plot(mids, ds, 'o-', label=f"{res}yr bins", markersize=4, alpha=0.7)

    ax.axvline(x=4.2, color='orange', linestyle='--', linewidth=2, label='4.2 ka event')
    ax.set_xlabel("Age (ka BP)")
    ax.set_ylabel("Divergence D (mon_frac - set_frac)")
    ax.set_title("Monument-Settlement Divergence Around the 4.2 ka Event")
    ax.legend(fontsize=8)
    ax.invert_xaxis()

    plt.tight_layout()
    plt.savefig(OUT_DIR / "climate_42ka.png", dpi=150, bbox_inches="tight")
    print(f"  Saved: {OUT_DIR / 'climate_42ka.png'}")
    plt.close()

    # ── Plot 4: Multi-proxy comparison at 500yr resolution ──
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot divergence D
    div_500 = divergence_series[500]
    valid_div = [b for b in div_500 if b["divergence_D"] is not None and b["bin_mid_ka_bp"] <= 9]
    if valid_div:
        mids = [b["bin_mid_ka_bp"] for b in valid_div]
        ds = [b["divergence_D"] for b in valid_div]
        ax.plot(mids, ds, 's-', color='red', linewidth=2, markersize=6, label='Divergence D', zorder=5)

    ax.set_ylabel("Divergence D", color='red')
    ax.tick_params(axis='y', labelcolor='red')

    # Overlay climate on secondary axis
    ax2 = ax.twinx()
    for clim_key, clim_rec in climate_records.items():
        clim_bins = bin_climate_data(clim_rec, 0.5, time_range_ka=(0, 9))
        if clim_bins:
            mids = [b["bin_mid_ka"] for b in clim_bins]
            vals = [b["value"] for b in clim_bins]
            ax2.plot(mids, vals, 'o--', color=colors.get(clim_key, "gray"),
                    label=clim_rec["name"][:25], alpha=0.7, markersize=3)

    ax2.set_ylabel("δ18O (‰ VPDB)")

    ax.set_xlabel("Age (ka BP)")
    ax.set_title("Multi-Proxy Climate vs Monument-Settlement Divergence (500yr bins)")
    ax.axvline(x=4.2, color='orange', linestyle='--', linewidth=1.5, alpha=0.5)
    ax.text(4.2, ax.get_ylim()[1]*0.95, "4.2 ka", ha='center', fontsize=8, color='orange')
    ax.invert_xaxis()

    # Combined legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=7, loc='upper right')

    plt.tight_layout()
    plt.savefig(OUT_DIR / "climate_multiproxy.png", dpi=150, bbox_inches="tight")
    print(f"  Saved: {OUT_DIR / 'climate_multiproxy.png'}")
    plt.close()


if __name__ == "__main__":
    main()
