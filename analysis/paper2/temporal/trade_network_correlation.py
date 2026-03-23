#!/usr/bin/env python3
"""
Directive B: Trade Network Proxy Correlation

Correlates radiocarbon date density in the Egypt-Iran trade corridor
with monument-settlement divergence (D) from the Great Circle analysis.

The hypothesis: if the Great Circle alignment reflects a trade/communication
network, divergence D should correlate with periods of intensified long-distance
exchange across the Egypt-Iran corridor.
"""

import csv
import json
import math
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
EARTH_R_KM = 6371
POLE_LAT = 59.682122
POLE_LON = -138.646087

# Egypt-Iran corridor bounding box
CORRIDOR_LAT_MIN = 20.0
CORRIDOR_LAT_MAX = 38.0
CORRIDOR_LON_MIN = 25.0
CORRIDOR_LON_MAX = 60.0

# Bin parameters
BIN_START = -5000   # 5000 BCE
BIN_END = 500       # 500 CE
BIN_WIDTH = 500

# Paths (relative to project root)
PROJECT_ROOT = os.path.expanduser("~/megalith_site_research")
P3K14C_PATH = os.path.join(PROJECT_ROOT, "data/p3k14c/p3k14c_data.csv")
DIVERGENCE_PATH = os.path.join(PROJECT_ROOT, "results/temporal_divergence_bins.json")
OUT_PLOT = os.path.join(PROJECT_ROOT, "results/trade_network_correlation.png")
OUT_JSON = os.path.join(PROJECT_ROOT, "results/trade_network_correlation.json")

# ---------------------------------------------------------------------------
# Trade-route timeline (hardcoded markers)
# ---------------------------------------------------------------------------
TRADE_MARKERS = [
    (-3200, "First lapis lazuli in Egypt (Naqada)"),
    (-3100, "Earliest Mesopotamian cylinder seals"),
    (-2600, "First Indus seals in Mesopotamia (Ur)"),
    (-2500, "Peak of Dilmun maritime trade"),
    (-2350, "Akkadian Empire unifies Mesopotamia"),
    (-2200, "Akkadian collapse"),
    (-2000, "Ur III collapse"),
    (-1900, "End of Indus-Mesopotamia trade"),
    (-1650, "Hyksos period in Egypt"),
    (-1500, "Egyptian New Kingdom trade with Punt"),
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _pearson(x, y):
    """Pearson correlation coefficient and two-tailed p-value."""
    n = len(x)
    if n < 3:
        return None, None
    mx = sum(x) / n
    my = sum(y) / n
    sxx = sum((xi - mx) ** 2 for xi in x)
    syy = sum((yi - my) ** 2 for yi in y)
    sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    if sxx == 0 or syy == 0:
        return 0.0, 1.0
    r = sxy / math.sqrt(sxx * syy)
    # t-test for significance
    if abs(r) >= 1.0:
        return r, 0.0
    t = r * math.sqrt((n - 2) / (1 - r * r))
    # Approximate two-tailed p from t distribution (df = n-2)
    p = _t_to_p(t, n - 2)
    return r, p


def _rank(arr):
    """Return ranks (1-based, average ties)."""
    indexed = sorted(enumerate(arr), key=lambda x: x[1])
    ranks = [0.0] * len(arr)
    i = 0
    while i < len(indexed):
        j = i
        while j < len(indexed) and indexed[j][1] == indexed[i][1]:
            j += 1
        avg_rank = (i + j + 1) / 2.0  # 1-based average
        for k in range(i, j):
            ranks[indexed[k][0]] = avg_rank
        i = j
    return ranks


def _spearman(x, y):
    """Spearman rank correlation coefficient and approximate p-value."""
    n = len(x)
    if n < 3:
        return None, None
    rx = _rank(x)
    ry = _rank(y)
    return _pearson(rx, ry)


def _t_to_p(t, df):
    """Approximate two-tailed p-value from t statistic using regularised
    incomplete beta function approximation (good enough for df >= 1)."""
    # Use the relationship: p = I_{df/(df+t^2)}(df/2, 1/2)
    # We'll use a simple numerical approximation via the survival function
    x = df / (df + t * t)
    a = df / 2.0
    b = 0.5
    # Regularised incomplete beta via continued fraction (Lentz)
    ibeta = _regularised_incomplete_beta(x, a, b)
    p = ibeta  # two-tailed
    return min(p, 1.0)


def _regularised_incomplete_beta(x, a, b, max_iter=200, tol=1e-12):
    """Regularised incomplete beta function I_x(a, b) via series expansion."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    # Use the series expansion for I_x(a, b)
    # First compute ln(B(a,b)) via lgamma
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(a * math.log(x) + b * math.log(1 - x) - lbeta) / a

    # Series: sum_{n=0}^{inf} (prod_{m=1}^{n} (m-b)/(a+m)) * x^n / 1
    s = 1.0
    term = 1.0
    for n in range(1, max_iter):
        term *= (n - b) / (a + n) * x
        s += term
        if abs(term) < tol:
            break
    return front * s


def try_scipy_correlations(x, y):
    """Try to use scipy for more accurate p-values; fall back to manual."""
    try:
        from scipy.stats import pearsonr, spearmanr
        pr, pp = pearsonr(x, y)
        sr, sp = spearmanr(x, y)
        return (pr, pp, "scipy"), (sr, sp, "scipy")
    except ImportError:
        pr, pp = _pearson(x, y)
        sr, sp = _spearman(x, y)
        return (pr, pp, "manual"), (sr, sp, "manual")


# ---------------------------------------------------------------------------
# 1. Extract corridor dates from p3k14c
# ---------------------------------------------------------------------------

def load_corridor_dates():
    """Load p3k14c dates within the Egypt-Iran corridor bounding box."""
    dates = []
    with open(P3K14C_PATH, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                lat = float(row["Lat"])
                lon = float(row["Long"])
                age = float(row["Age"])
            except (ValueError, KeyError):
                continue
            if (CORRIDOR_LAT_MIN <= lat <= CORRIDOR_LAT_MAX and
                    CORRIDOR_LON_MIN <= lon <= CORRIDOR_LON_MAX):
                cal_year = 1950 - age  # approximate calendar year
                dates.append({
                    "site_id": row.get("SiteID", ""),
                    "site_name": row.get("SiteName", ""),
                    "lat": lat,
                    "lon": lon,
                    "age_bp": age,
                    "cal_year": cal_year,
                })
    return dates


# ---------------------------------------------------------------------------
# 2. Bin corridor dates
# ---------------------------------------------------------------------------

def bin_dates(dates):
    """Bin corridor dates into 500-year windows."""
    bins = []
    for start in range(BIN_START, BIN_END, BIN_WIDTH):
        end = start + BIN_WIDTH
        count = sum(1 for d in dates if start <= d["cal_year"] < end)
        bins.append({
            "bin_start": start,
            "bin_end": end,
            "label": f"{abs(start)} {'BCE' if start < 0 else 'CE'} - {abs(end)} {'BCE' if end < 0 else 'CE'}",
            "corridor_date_count": count,
        })
    return bins


# ---------------------------------------------------------------------------
# 3. Load divergence data
# ---------------------------------------------------------------------------

def load_divergence():
    """Load temporal divergence bins JSON."""
    with open(DIVERGENCE_PATH, "r") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# 4. Correlate
# ---------------------------------------------------------------------------

def correlate(corridor_bins, divergence_data):
    """Compute correlations between corridor density and divergence D."""
    results = {}
    for db_key in ("pleiades_bins", "p3k14c_bins"):
        db_bins = divergence_data[db_key]
        # Build lookup: bin_start -> D
        d_lookup = {}
        for b in db_bins:
            if b["divergence_D"] is not None:
                d_lookup[b["bin_start"]] = b["divergence_D"]

        # Align
        density_vals = []
        d_vals = []
        aligned_labels = []
        for cb in corridor_bins:
            bs = cb["bin_start"]
            if bs in d_lookup:
                density_vals.append(float(cb["corridor_date_count"]))
                d_vals.append(float(d_lookup[bs]))
                aligned_labels.append(cb["label"])

        if len(density_vals) < 3:
            results[db_key] = {
                "n_aligned_bins": len(density_vals),
                "pearson_r": None,
                "pearson_p": None,
                "spearman_r": None,
                "spearman_p": None,
                "method": "insufficient_data",
                "aligned_bins": aligned_labels,
            }
            continue

        (pr, pp, pm), (sr, sp, sm) = try_scipy_correlations(density_vals, d_vals)

        results[db_key] = {
            "n_aligned_bins": len(density_vals),
            "pearson_r": round(pr, 4) if pr is not None else None,
            "pearson_p": round(pp, 6) if pp is not None else None,
            "spearman_r": round(sr, 4) if sr is not None else None,
            "spearman_p": round(sp, 6) if sp is not None else None,
            "method": pm,
            "aligned_bins": aligned_labels,
            "density_values": density_vals,
            "D_values": d_vals,
        }

    return results


# ---------------------------------------------------------------------------
# 5. Plot
# ---------------------------------------------------------------------------

def make_plot(corridor_bins, divergence_data, correlation_results):
    """Generate dual-axis time series with trade markers."""
    fig, ax1 = plt.subplots(figsize=(16, 8))
    ax2 = ax1.twinx()

    # X-axis: bin midpoints
    midpoints = [(b["bin_start"] + b["bin_end"]) / 2 for b in corridor_bins]
    densities = [b["corridor_date_count"] for b in corridor_bins]

    # Plot corridor date density on right axis
    ax2.bar(midpoints, densities, width=400, alpha=0.25, color="steelblue",
            label="Corridor date density", zorder=1)
    ax2.set_ylabel("Corridor 14C date count (per 500-yr bin)", color="steelblue",
                    fontsize=11)
    ax2.tick_params(axis="y", labelcolor="steelblue")

    # Plot divergence D for both databases on left axis
    for db_key, color, marker, lbl in [
        ("pleiades_bins", "darkred", "o", "Pleiades D"),
        ("p3k14c_bins", "darkorange", "s", "p3k14c D"),
    ]:
        db_bins = divergence_data[db_key]
        xs, ys = [], []
        for b in db_bins:
            if b["divergence_D"] is not None:
                xs.append((b["bin_start"] + b["bin_end"]) / 2)
                ys.append(b["divergence_D"])
        ax1.plot(xs, ys, color=color, marker=marker, linewidth=2,
                 markersize=7, label=lbl, zorder=3)

    ax1.set_ylabel("Divergence D (monument z - settlement z)", fontsize=11)
    ax1.set_xlabel("Calendar Year", fontsize=11)
    ax1.axhline(y=0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)

    # Trade markers as vertical lines
    y_min, y_max = ax1.get_ylim()
    for i, (year, label) in enumerate(TRADE_MARKERS):
        ax1.axvline(x=year, color="green", linestyle=":", linewidth=1.0,
                     alpha=0.7, zorder=2)
        # Stagger label heights to reduce overlap
        y_pos = y_max - (y_max - y_min) * (0.05 + 0.07 * (i % 5))
        ax1.text(year + 20, y_pos, label, fontsize=7, color="darkgreen",
                 rotation=90, va="top", ha="left", alpha=0.9)

    # Add correlation text
    text_lines = []
    for db_key, db_label in [("pleiades_bins", "Pleiades"), ("p3k14c_bins", "p3k14c")]:
        cr = correlation_results.get(db_key, {})
        pr = cr.get("pearson_r")
        pp = cr.get("pearson_p")
        sr = cr.get("spearman_r")
        sp = cr.get("spearman_p")
        n = cr.get("n_aligned_bins", 0)
        if pr is not None:
            text_lines.append(
                f"{db_label} (n={n}): Pearson r={pr:.3f} (p={pp:.4f}), "
                f"Spearman r={sr:.3f} (p={sp:.4f})"
            )
        else:
            text_lines.append(f"{db_label}: insufficient aligned bins")

    ax1.text(0.02, 0.02, "\n".join(text_lines), transform=ax1.transAxes,
             fontsize=8, verticalalignment="bottom",
             bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow",
                       edgecolor="gray", alpha=0.9))

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left",
               fontsize=9, framealpha=0.9)

    ax1.set_title(
        "Directive B: Egypt-Iran Corridor Date Density vs. Monument-Settlement Divergence\n"
        "with Trade Network Markers",
        fontsize=13, fontweight="bold"
    )

    # Format x-axis
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(500))
    ax1.set_xlim(BIN_START - 200, BIN_END + 200)
    plt.tight_layout()
    fig.savefig(OUT_PLOT, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved to {OUT_PLOT}")


# ---------------------------------------------------------------------------
# 6. Save JSON
# ---------------------------------------------------------------------------

def save_json(corridor_dates, corridor_bins, divergence_data, correlation_results):
    """Save all results to JSON."""
    # Summarise corridor dates by site
    site_summary = defaultdict(lambda: {"count": 0, "lat": None, "lon": None, "name": ""})
    for d in corridor_dates:
        sid = d["site_id"] or d["site_name"]
        site_summary[sid]["count"] += 1
        site_summary[sid]["lat"] = d["lat"]
        site_summary[sid]["lon"] = d["lon"]
        site_summary[sid]["name"] = d["site_name"]

    top_sites = sorted(site_summary.items(), key=lambda x: -x[1]["count"])[:20]

    output = {
        "meta": {
            "analysis": "Directive B: Trade Network Proxy Correlation",
            "date": "2026-03-19",
            "great_circle_pole": {"lat": POLE_LAT, "lon": POLE_LON},
            "corridor_bounding_box": {
                "lat_min": CORRIDOR_LAT_MIN,
                "lat_max": CORRIDOR_LAT_MAX,
                "lon_min": CORRIDOR_LON_MIN,
                "lon_max": CORRIDOR_LON_MAX,
            },
            "bin_width_years": BIN_WIDTH,
            "date_conversion": "cal_year = 1950 - Age (approximate; no IntCal calibration)",
            "earth_radius_km": EARTH_R_KM,
        },
        "corridor_summary": {
            "total_dates_in_corridor": len(corridor_dates),
            "unique_sites": len(site_summary),
            "top_20_sites": [
                {"site_id": sid, "name": info["name"], "lat": info["lat"],
                 "lon": info["lon"], "date_count": info["count"]}
                for sid, info in top_sites
            ],
        },
        "corridor_bins": corridor_bins,
        "trade_markers": [
            {"year": y, "event": e} for y, e in TRADE_MARKERS
        ],
        "correlations": {
            "pleiades": correlation_results.get("pleiades_bins", {}),
            "p3k14c": correlation_results.get("p3k14c_bins", {}),
        },
        "interpretation": {
            "question": "Does corridor trade activity correlate with monument-settlement divergence?",
            "note": "Positive correlation would suggest that periods of high long-distance exchange "
                    "coincide with periods when monuments preferentially cluster near the Great Circle "
                    "relative to settlements.",
        },
    }

    # Remove non-serialisable items from correlation results
    for key in ("pleiades", "p3k14c"):
        cr = output["correlations"].get(key, {})
        cr.pop("density_values", None)
        cr.pop("D_values", None)

    with open(OUT_JSON, "w") as f:
        json.dump(output, f, indent=2)
    print(f"JSON saved to {OUT_JSON}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Directive B: Trade Network Proxy Correlation")
    print("=" * 70)

    # 1. Load corridor dates
    print("\n[1] Loading p3k14c dates in Egypt-Iran corridor "
          f"({CORRIDOR_LAT_MIN}-{CORRIDOR_LAT_MAX}N, "
          f"{CORRIDOR_LON_MIN}-{CORRIDOR_LON_MAX}E)...")
    corridor_dates = load_corridor_dates()
    print(f"    Found {len(corridor_dates)} dates in corridor")

    # 2. Bin dates
    print("\n[2] Binning into 500-year windows...")
    corridor_bins = bin_dates(corridor_dates)
    for b in corridor_bins:
        print(f"    {b['label']:>30s}: {b['corridor_date_count']:5d} dates")

    # 3. Load divergence
    print("\n[3] Loading temporal divergence results...")
    divergence_data = load_divergence()
    print(f"    Pleiades bins: {len(divergence_data['pleiades_bins'])}")
    print(f"    p3k14c bins:   {len(divergence_data['p3k14c_bins'])}")

    # 4. Correlate
    print("\n[4] Computing correlations...")
    correlation_results = correlate(corridor_bins, divergence_data)
    for db_key, db_label in [("pleiades_bins", "Pleiades"), ("p3k14c_bins", "p3k14c")]:
        cr = correlation_results.get(db_key, {})
        pr = cr.get("pearson_r")
        sr = cr.get("spearman_r")
        n = cr.get("n_aligned_bins", 0)
        method = cr.get("method", "n/a")
        if pr is not None:
            print(f"    {db_label} (n={n}, method={method}):")
            print(f"      Pearson  r = {pr:+.4f}  p = {cr['pearson_p']:.6f}")
            print(f"      Spearman r = {sr:+.4f}  p = {cr['spearman_p']:.6f}")
        else:
            print(f"    {db_label}: insufficient aligned bins for correlation")

    # 5. Plot
    print("\n[5] Generating plot...")
    make_plot(corridor_bins, divergence_data, correlation_results)

    # 6. Save JSON
    print("\n[6] Saving JSON results...")
    save_json(corridor_dates, corridor_bins, divergence_data, correlation_results)

    print("\nDone.")


if __name__ == "__main__":
    main()
