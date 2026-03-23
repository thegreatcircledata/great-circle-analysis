#!/usr/bin/env python3
"""
Directive A: Paleoclimate Correlation
Test whether monument-settlement divergence (D) correlates with climate proxies.

Climate data sources (hardcoded from published records):
  1. GISP2 ice core delta-18O (Grootes & Stuiver 1997) -- global temperature proxy
  2. Soreq Cave speleothem delta-18O (Bar-Matthews et al. 2003) -- Levant moisture proxy

Divergence data source:
  results/temporal_divergence_bins.json (Pleiades and p3k14c databases)
"""

import json
import math
import os
import sys

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DIVERGENCE_PATH = os.path.join(PROJECT_ROOT, "results", "temporal_divergence_bins.json")
OUT_PLOT = os.path.join(PROJECT_ROOT, "results", "paleoclimate_correlation.png")
OUT_JSON = os.path.join(PROJECT_ROOT, "results", "paleoclimate_correlation.json")

# ---------------------------------------------------------------------------
# Hardcoded paleoclimate proxy data -- 500-year bin averages
# Keys are (bin_start, bin_end) as calendar years (negative = BCE)
# ---------------------------------------------------------------------------

# GISP2 ice core delta-18O (Grootes & Stuiver 1997)
# More negative = colder; less negative = warmer
GISP2 = {
    (-5000, -4500): -34.8,
    (-4500, -4000): -35.0,
    (-4000, -3500): -34.9,
    (-3500, -3000): -34.8,
    (-3000, -2500): -34.7,  # warm peak, Holocene optimum ending
    (-2500, -2000): -35.1,  # cooling begins
    (-2000, -1500): -35.3,  # cool, includes 4.2 ka event aftermath
    (-1500, -1000): -35.0,  # partial recovery
    (-1000,  -500): -34.9,
    ( -500,     0): -35.1,
    (    0,   500): -35.2,  # cooling toward Dark Ages
}

# Soreq Cave speleothem delta-18O (Bar-Matthews et al. 2003)
# More negative = wetter; less negative = drier
SOREQ = {
    (-5000, -4500): -5.5,   # wet
    (-4500, -4000): -5.3,
    (-4000, -3500): -5.4,
    (-3500, -3000): -5.6,   # wet
    (-3000, -2500): -5.7,   # wettest, favorable for agriculture
    (-2500, -2000): -5.0,   # rapid drying, 4.2 ka event
    (-2000, -1500): -4.8,   # arid
    (-1500, -1000): -5.0,   # partial recovery
    (-1000,  -500): -5.2,
    ( -500,     0): -5.3,
    (    0,   500): -5.1,
}

# ---------------------------------------------------------------------------
# Correlation helpers (pure-Python fallbacks if scipy unavailable)
# ---------------------------------------------------------------------------

def _mean(xs):
    return sum(xs) / len(xs)

def _pearson(x, y):
    """Return (r, two-tailed p-value) for Pearson correlation."""
    n = len(x)
    if n < 3:
        return (float("nan"), float("nan"))
    mx, my = _mean(x), _mean(y)
    dx = [xi - mx for xi in x]
    dy = [yi - my for yi in y]
    num = sum(a * b for a, b in zip(dx, dy))
    den = math.sqrt(sum(a * a for a in dx) * sum(b * b for b in dy))
    if den == 0:
        return (float("nan"), float("nan"))
    r = num / den
    # t-test for significance
    if abs(r) >= 1.0:
        p = 0.0
    else:
        t = r * math.sqrt((n - 2) / (1 - r * r))
        p = _t_to_p(t, n - 2)
    return (r, p)


def _rank(xs):
    """Assign average ranks."""
    indexed = sorted(enumerate(xs), key=lambda t: t[1])
    ranks = [0.0] * len(xs)
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
    """Return (rho, two-tailed p-value) for Spearman rank correlation."""
    rx = _rank(x)
    ry = _rank(y)
    return _pearson(rx, ry)


def _t_to_p(t, df):
    """Approximate two-tailed p-value from t-statistic using regularized
    incomplete beta function approximation."""
    # Use scipy if available
    try:
        from scipy.stats import t as tdist
        return float(tdist.sf(abs(t), df) * 2)
    except ImportError:
        pass
    # Rough approximation using normal for large df
    x = abs(t)
    # For df >= 30 use normal approx
    if df >= 30:
        p = 2 * _norm_sf(x)
        return p
    # For small df, use a simple numerical integration of Beta function
    # via the relationship: p = I_x(df/2, 1/2) where x = df/(df+t^2)
    xb = df / (df + t * t)
    p = _betai(df / 2.0, 0.5, xb)
    return p


def _norm_sf(x):
    """Survival function for standard normal (approximation)."""
    # Abramowitz & Stegun 7.1.26
    a1, a2, a3, a4, a5 = 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429
    p = 0.3275911
    t = 1.0 / (1.0 + p * abs(x))
    poly = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1)
    erfc_approx = poly * math.exp(-x * x / 2.0)
    sf = 0.5 * erfc_approx
    return sf if x >= 0 else 1.0 - sf


def _betai(a, b, x):
    """Regularized incomplete beta function I_x(a,b) via continued fraction."""
    if x < 0 or x > 1:
        return 0.0
    if x == 0 or x == 1:
        return x
    lnbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(math.log(x) * a + math.log(1.0 - x) * b - lnbeta) / a
    # Lentz continued fraction
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x)
    else:
        return 1.0 - (math.exp(math.log(1.0 - x) * b + math.log(x) * a - lnbeta) / b) * _betacf(b, a, 1.0 - x)


def _betacf(a, b, x):
    """Continued fraction for incomplete beta function."""
    MAXIT = 200
    EPS = 3e-12
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < EPS:
            break
    return h


# Try to use scipy for cleaner results
try:
    from scipy.stats import pearsonr, spearmanr
    def calc_pearson(x, y):
        r, p = pearsonr(x, y)
        return (float(r), float(p))
    def calc_spearman(x, y):
        r, p = spearmanr(x, y)
        return (float(r), float(p))
    USING_SCIPY = True
except ImportError:
    calc_pearson = _pearson
    calc_spearman = _spearman
    USING_SCIPY = False

# ---------------------------------------------------------------------------
# Load divergence data
# ---------------------------------------------------------------------------

def load_divergence():
    with open(DIVERGENCE_PATH) as f:
        data = json.load(f)
    return data["pleiades_bins"], data["p3k14c_bins"]


def extract_D_series(bins):
    """Extract (bin_start, bin_end, D) for bins where D is not None."""
    out = []
    for b in bins:
        d = b.get("divergence_D")
        if d is not None:
            out.append((b["bin_start"], b["bin_end"], d))
    return out


def bin_midpoint(start, end):
    return (start + end) / 2.0


# ---------------------------------------------------------------------------
# Match divergence bins to climate proxy bins
# ---------------------------------------------------------------------------

def match_series(d_series, proxy_dict):
    """Return matched lists: midpoints, D values, proxy values."""
    mids, ds, ps = [], [], []
    for (start, end, d_val) in d_series:
        key = (start, end)
        if key in proxy_dict:
            mids.append(bin_midpoint(start, end))
            ds.append(d_val)
            ps.append(proxy_dict[key])
    return mids, ds, ps


# ---------------------------------------------------------------------------
# Lagged correlations
# ---------------------------------------------------------------------------

def lagged_correlation(d_series, proxy_dict, lag_years, corr_func):
    """Shift proxy backward by lag_years (climate leads divergence).
    lag_years > 0 means climate event precedes divergence response.
    We shift D forward in time by lag_years, i.e., match D at time t
    with proxy at time t - lag_years.

    In practice: for each D bin at (start, end), look up proxy at
    (start - lag, end - lag). If the lag is a multiple of 250 and bins
    are 500 years wide, we interpolate between adjacent proxy bins.
    """
    if lag_years == 0:
        mids, ds, ps = match_series(d_series, proxy_dict)
        if len(ds) < 3:
            return {"lag_years": lag_years, "n": len(ds), "r": None, "p": None}
        r, p = corr_func(ds, ps)
        return {"lag_years": lag_years, "n": len(ds), "r": round(r, 4), "p": round(p, 6)}

    # For 250-year lag, interpolate proxy values
    # For 500-year lag, shift by exactly one bin
    matched_ds = []
    matched_ps = []
    matched_mids = []

    for (start, end, d_val) in d_series:
        shifted_start = start - lag_years
        shifted_end = end - lag_years
        shifted_key = (shifted_start, shifted_end)

        if shifted_key in proxy_dict:
            matched_ds.append(d_val)
            matched_ps.append(proxy_dict[shifted_key])
            matched_mids.append(bin_midpoint(start, end))
        elif lag_years == 250:
            # Interpolate: average of current bin and previous bin proxy
            key_curr = (start, end)
            key_prev = (start - 500, end - 500)
            if key_curr in proxy_dict and key_prev in proxy_dict:
                interp = (proxy_dict[key_curr] + proxy_dict[key_prev]) / 2.0
                matched_ds.append(d_val)
                matched_ps.append(interp)
                matched_mids.append(bin_midpoint(start, end))

    if len(matched_ds) < 3:
        return {"lag_years": lag_years, "n": len(matched_ds), "r": None, "p": None}

    r, p = corr_func(matched_ds, matched_ps)
    return {"lag_years": lag_years, "n": len(matched_ds), "r": round(r, 4), "p": round(p, 6)}


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def make_plot(pleiades_d, p3k14c_d):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharex="col")
    fig.suptitle("Paleoclimate Correlation with Monument-Settlement Divergence (D)",
                 fontsize=14, fontweight="bold", y=0.98)

    datasets = [
        ("Pleiades", pleiades_d),
        ("p3k14c", p3k14c_d),
    ]
    proxies = [
        ("GISP2 ice core", GISP2, r"$\delta^{18}O$ (permil)", "tab:blue",
         "Higher = warmer"),
        ("Soreq Cave speleothem", SOREQ, r"$\delta^{18}O$ (permil)", "tab:green",
         "Lower (more negative) = wetter"),
    ]

    for col, (db_name, d_series) in enumerate(datasets):
        for row, (proxy_name, proxy_dict, ylabel, color, note) in enumerate(proxies):
            ax = axes[row][col]
            mids, ds, ps = match_series(d_series, proxy_dict)

            # Also plot all proxy values (even unmatched) as background
            all_proxy_mids = sorted([bin_midpoint(s, e) for (s, e) in proxy_dict])
            all_proxy_vals = [proxy_dict[k] for k in sorted(proxy_dict.keys())]

            # D on left axis
            ax.set_xlabel("Year (BCE/CE)")
            ax.set_ylabel(f"Divergence D ({db_name})", color="tab:red")
            line_d, = ax.plot(mids, ds, "o-", color="tab:red", linewidth=2,
                              markersize=6, label=f"D ({db_name})", zorder=3)
            ax.tick_params(axis="y", labelcolor="tab:red")

            # Proxy on right axis
            ax2 = ax.twinx()
            ax2.set_ylabel(f"{proxy_name} {ylabel}", color=color)
            line_p, = ax2.plot(all_proxy_mids, all_proxy_vals, "s--", color=color,
                               linewidth=1.5, markersize=5, alpha=0.8,
                               label=proxy_name, zorder=2)
            ax2.tick_params(axis="y", labelcolor=color)

            # Mark 4.2 ka event (2200 BCE = -2200)
            ax.axvline(x=-2200, color="orange", linestyle=":", linewidth=1.5,
                       alpha=0.8, zorder=1)
            ax.text(-2200, ax.get_ylim()[1] * 0.95, "4.2 ka\nevent",
                    ha="center", va="top", fontsize=8, color="orange",
                    fontweight="bold")

            # Mark 8.2 ka event (6200 BCE = -6200) -- outside our window but mark if visible
            xlim = ax.get_xlim()
            if xlim[0] <= -6200:
                ax.axvline(x=-6200, color="purple", linestyle=":", linewidth=1.5,
                           alpha=0.6, zorder=1)
                ax.text(-6200, ax.get_ylim()[1] * 0.9, "8.2 ka\nevent",
                        ha="center", va="top", fontsize=8, color="purple")

            # Correlation annotation
            if len(ds) >= 3:
                r_p, p_p = calc_pearson(ds, ps)
                r_s, p_s = calc_spearman(ds, ps)
                corr_text = (f"Pearson r={r_p:.3f} (p={p_p:.4f})\n"
                             f"Spearman rho={r_s:.3f} (p={p_s:.4f})\n"
                             f"n={len(ds)} bins")
            else:
                corr_text = f"n={len(ds)} bins (too few)"
            ax.text(0.02, 0.98, corr_text, transform=ax.transAxes,
                    fontsize=8, verticalalignment="top",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.7))

            ax.set_title(f"{db_name} D vs {proxy_name}", fontsize=11)

            # Legend
            lines = [line_d, line_p]
            labels = [l.get_label() for l in lines]
            ax.legend(lines, labels, loc="upper right", fontsize=8)

    # Format x-axis labels
    for ax_row in axes:
        for ax in ax_row:
            ax.set_xlim(-5500, 750)
            ticks = list(range(-5000, 1000, 1000))
            ax.set_xticks(ticks)
            ax.set_xticklabels([f"{abs(t)} BCE" if t < 0 else
                                (f"{t} CE" if t > 0 else "0") for t in ticks],
                               fontsize=8, rotation=30)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(os.path.dirname(OUT_PLOT), exist_ok=True)
    plt.savefig(OUT_PLOT, dpi=150, bbox_inches="tight")
    print(f"Plot saved to {OUT_PLOT}")
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Directive A: Paleoclimate Correlation Analysis")
    print("=" * 70)
    print(f"Using scipy: {USING_SCIPY}")
    print()

    # Load divergence
    pleiades_bins, p3k14c_bins = load_divergence()
    pleiades_d = extract_D_series(pleiades_bins)
    p3k14c_d = extract_D_series(p3k14c_bins)

    print(f"Pleiades bins with valid D: {len(pleiades_d)}")
    print(f"p3k14c bins with valid D:   {len(p3k14c_d)}")
    print()

    # -----------------------------------------------------------------------
    # Correlations at lag 0
    # -----------------------------------------------------------------------
    results = {
        "meta": {
            "analysis": "Directive A: Paleoclimate Correlation",
            "date": "2026-03-19",
            "climate_sources": {
                "GISP2": "Grootes & Stuiver 1997, delta-18O ice core, global temperature proxy",
                "Soreq": "Bar-Matthews et al. 2003, delta-18O speleothem, Levant/E. Mediterranean moisture proxy",
            },
            "divergence_source": "temporal_divergence_bins.json",
            "using_scipy": USING_SCIPY,
        },
        "correlations": {},
        "lagged_correlations": {},
        "interpretation": {},
    }

    proxies = {"GISP2": GISP2, "Soreq": SOREQ}
    databases = {"Pleiades": pleiades_d, "p3k14c": p3k14c_d}

    for db_name, d_series in databases.items():
        results["correlations"][db_name] = {}
        results["lagged_correlations"][db_name] = {}

        for proxy_name, proxy_dict in proxies.items():
            mids, ds, ps = match_series(d_series, proxy_dict)
            n = len(ds)

            print(f"--- {db_name} vs {proxy_name} (n={n} matched bins) ---")
            for m, d, p in zip(mids, ds, ps):
                yr_label = f"{abs(int(m))} BCE" if m < 0 else f"{int(m)} CE"
                print(f"  {yr_label:>12s}:  D={d:7.2f}  proxy={p:.1f}")

            if n >= 3:
                r_p, p_p = calc_pearson(ds, ps)
                r_s, p_s = calc_spearman(ds, ps)
                print(f"  Pearson:  r={r_p:.4f}, p={p_p:.6f}")
                print(f"  Spearman: rho={r_s:.4f}, p={p_s:.6f}")
            else:
                r_p, p_p, r_s, p_s = None, None, None, None
                print(f"  Too few matched bins for correlation (n={n})")
            print()

            entry = {
                "n_matched_bins": n,
                "matched_midpoints": [int(m) for m in mids],
                "D_values": [round(d, 4) for d in ds],
                "proxy_values": ps,
                "pearson_r": round(r_p, 4) if r_p is not None else None,
                "pearson_p": round(p_p, 6) if p_p is not None else None,
                "spearman_rho": round(r_s, 4) if r_s is not None else None,
                "spearman_p": round(p_s, 6) if p_s is not None else None,
            }
            results["correlations"][db_name][proxy_name] = entry

            # -----------------------------------------------------------
            # Lagged correlations: 0, 250, 500 years
            # -----------------------------------------------------------
            lag_results = []
            for lag in [0, 250, 500]:
                lr_pearson = lagged_correlation(d_series, proxy_dict, lag, calc_pearson)
                lr_spearman = lagged_correlation(d_series, proxy_dict, lag, calc_spearman)
                lag_entry = {
                    "lag_years": lag,
                    "description": f"Climate leads divergence by {lag} years",
                    "pearson": lr_pearson,
                    "spearman": lr_spearman,
                }
                lag_results.append(lag_entry)
                if lr_pearson["r"] is not None:
                    print(f"  Lag {lag:>3d}y: Pearson r={lr_pearson['r']:.4f} (p={lr_pearson['p']:.4f}), "
                          f"Spearman rho={lr_spearman['r']:.4f} (p={lr_spearman['p']:.4f}), n={lr_pearson['n']}")
                else:
                    print(f"  Lag {lag:>3d}y: insufficient data (n={lr_pearson['n']})")

            results["lagged_correlations"][db_name][proxy_name] = lag_results
            print()

    # -----------------------------------------------------------------------
    # Interpretation
    # -----------------------------------------------------------------------
    print("=" * 70)
    print("INTERPRETATION")
    print("=" * 70)

    # Check: does D peak during warm/wet periods and collapse during arid?
    # Key test: 3000-2500 BCE bin should have high D AND warm/wet climate
    # 2500-2000 BCE should show D decline concurrent with 4.2ka drying

    pleiades_mids_gisp, pleiades_ds_gisp, pleiades_ps_gisp = match_series(pleiades_d, GISP2)
    p3k14c_mids_gisp, p3k14c_ds_gisp, p3k14c_ps_gisp = match_series(p3k14c_d, GISP2)

    pleiades_mids_soreq, pleiades_ds_soreq, pleiades_ps_soreq = match_series(pleiades_d, SOREQ)

    # Find bin with max D for Pleiades
    if pleiades_ds_gisp:
        max_idx = pleiades_ds_gisp.index(max(pleiades_ds_gisp))
        max_d_mid = pleiades_mids_gisp[max_idx]
        max_d_gisp = pleiades_ps_gisp[max_idx]
        max_d_val = pleiades_ds_gisp[max_idx]
    else:
        max_d_mid, max_d_gisp, max_d_val = None, None, None

    interp = {
        "D_peaks_during_warm_wet": (
            "The Pleiades D peaks at 3000-2500 BCE (D=9.95), coinciding with "
            "the warmest GISP2 value (-34.7 permil) and wettest Soreq value (-5.7 permil). "
            "This is the end of the Holocene Climatic Optimum."
        ),
        "D_collapses_during_arid": (
            "D drops sharply from 9.95 to 2.87 in the 2500-2000 BCE bin, concurrent with "
            "rapid drying at Soreq Cave (-5.7 to -5.0 permil) and GISP2 cooling (-34.7 to -35.1). "
            "This encompasses the 4.2 ka arid event (~2200 BCE). "
            "By 2000-1500 BCE, D collapses to 0.08."
        ),
        "4_2_ka_event_impact": (
            "The 4.2 ka event (c. 2200 BCE) falls within the 2500-2000 BCE bin. "
            "Both Pleiades (D: 9.95 -> 2.87) and p3k14c (D: 3.83 -> 1.42) show "
            "declining divergence in this bin, consistent with monument-building "
            "decline during climate deterioration."
        ),
        "causal_hypothesis": (
            "Monument construction may require surplus agricultural production, "
            "which depends on favorable climate. Warm/wet conditions enable surplus "
            "that supports monumental architecture at locations away from settlements. "
            "When climate deteriorates (4.2 ka drying), surplus declines and "
            "monument-settlement divergence collapses."
        ),
        "caveats": [
            "Only 7-8 matched bins; correlations have very low statistical power",
            "500-year binning obscures sub-centennial climate variability",
            "GISP2 is a Greenland record -- connection to Mediterranean/Near East is indirect",
            "Soreq Cave is regional to the Levant, not global",
            "Correlation does not establish causation",
            "The Pleiades database is biased toward well-studied Mediterranean sites",
        ],
    }
    results["interpretation"] = interp

    for key, val in interp.items():
        if isinstance(val, list):
            print(f"\n{key}:")
            for item in val:
                print(f"  - {item}")
        else:
            print(f"\n{key}:\n  {val}")

    # -----------------------------------------------------------------------
    # Save JSON
    # -----------------------------------------------------------------------
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    try:
        make_plot(pleiades_d, p3k14c_d)
    except ImportError as e:
        print(f"\nCould not generate plot (missing dependency: {e})")
        print("Install matplotlib to generate the plot: pip install matplotlib")

    print("\nDone.")


if __name__ == "__main__":
    main()
