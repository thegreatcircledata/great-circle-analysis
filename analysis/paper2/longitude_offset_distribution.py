#!/usr/bin/env python3
"""
Longitude Offset Frequency Distribution
========================================
Generates frequency distribution of archaeological site longitude offsets
from Giza (31.134°E), checking for 36° periodic modes.

Outputs:
  - raw_histogram.png
  - land_normalized.png
  - kde_comparison.png
  - residual.png
  - results.json
"""

import json
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ── Config ──────────────────────────────────────────────────────────────
GIZA_LON = 31.134
GRID_INTERVAL = 36  # degrees
DATA_PATH = Path(__file__).parent.parent / "github-repo" / "data" / "merged_sites.csv"
OUT_DIR = Path(__file__).parent.parent / "outputs" / "longitude_offset_distribution"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Dark theme matching Great Circle brand
DARK_BG = '#1a1a2e'
DARK_FG = '#e0e0e0'
ACCENT = '#e94560'
ACCENT2 = '#0f3460'
GRID_LINE_COLOR = '#e9456080'
LAND_COLOR = '#16213e'

# 36° multiples for overlay
GRID_LINES = np.arange(GRID_INTERVAL, 360, GRID_INTERVAL)  # 36, 72, ... 324

def setup_style():
    """Configure dark matplotlib style."""
    rcParams.update({
        'figure.facecolor': DARK_BG,
        'axes.facecolor': DARK_BG,
        'axes.edgecolor': DARK_FG,
        'axes.labelcolor': DARK_FG,
        'text.color': DARK_FG,
        'xtick.color': DARK_FG,
        'ytick.color': DARK_FG,
        'grid.color': '#333355',
        'grid.alpha': 0.3,
        'font.size': 11,
        'axes.titlesize': 14,
        'figure.dpi': 150,
    })

def load_data():
    """Load sites and compute longitude offsets from Giza."""
    df = pd.read_csv(DATA_PATH)
    df = df.dropna(subset=['lon', 'lat'])
    # Compute offset: wrap to [0, 360)
    df['lon_offset'] = (df['lon'] - GIZA_LON) % 360
    print(f"Loaded {len(df):,} sites")
    return df

def add_grid_lines(ax, ymax=None):
    """Add 36° interval dashed lines."""
    for g in GRID_LINES:
        ax.axvline(g, color=GRID_LINE_COLOR, linestyle='--', linewidth=0.8, alpha=0.7)
    # Giza reference at 0°
    ax.axvline(0, color='#ffd700', linestyle='-', linewidth=1.2, alpha=0.8, label='Giza (0°)')

def plot_raw_histogram(offsets):
    """Plot 1: Raw histogram of longitude offsets."""
    fig, ax = plt.subplots(figsize=(12, 5))

    counts, bin_edges, patches = ax.hist(
        offsets, bins=360, range=(0, 360),
        color=ACCENT2, edgecolor='none', alpha=0.85
    )

    add_grid_lines(ax)
    ax.set_xlabel('Longitude offset from Giza (°)')
    ax.set_ylabel('Site count')
    ax.set_title('Site distribution by longitude offset from Giza')
    ax.set_xlim(0, 360)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, axis='y', alpha=0.2)

    fig.tight_layout()
    fig.savefig(OUT_DIR / 'raw_histogram.png', dpi=150, bbox_inches='tight',
                facecolor=DARK_BG)
    plt.close(fig)
    print("  → raw_histogram.png")
    return counts, bin_edges

def compute_land_area_proxy(offsets):
    """
    Estimate land area per 1° longitude band.
    Uses number of occupied 0.1° grid cells as proxy for sampled land area.
    """
    # For each 1° offset band, count unique 0.1° grid cells with sites
    offset_bins = np.floor(offsets).astype(int) % 360

    # We need lat too for grid cells — reload
    df = pd.read_csv(DATA_PATH)
    df = df.dropna(subset=['lon', 'lat'])
    df['lon_offset'] = (df['lon'] - GIZA_LON) % 360

    land_area = np.zeros(360)
    for b in range(360):
        mask = (df['lon_offset'] >= b) & (df['lon_offset'] < b + 1)
        subset = df[mask]
        if len(subset) == 0:
            land_area[b] = 0
            continue
        # Count unique 0.1° grid cells (lon_offset × lat)
        cells = set(zip(
            (subset['lon_offset'] * 10).astype(int),
            (subset['lat'] * 10).astype(int)
        ))
        land_area[b] = len(cells)

    return land_area

def plot_land_normalized(offsets, counts):
    """Plot 2: Land-area normalized histogram."""
    print("  Computing land area proxy...")
    land_area = compute_land_area_proxy(offsets)

    # Normalize: sites per grid cell (proxy for sites/km²)
    normalized = np.zeros(360)
    for i in range(360):
        if land_area[i] > 0:
            normalized[i] = counts[i] / land_area[i]

    fig, ax = plt.subplots(figsize=(12, 5))

    ax.bar(np.arange(360) + 0.5, normalized, width=1.0,
           color=ACCENT2, edgecolor='none', alpha=0.85)

    add_grid_lines(ax)
    ax.set_xlabel('Longitude offset from Giza (°)')
    ax.set_ylabel('Sites per occupied grid cell (land-area proxy)')
    ax.set_title('Land-area normalized site density by longitude offset')
    ax.set_xlim(0, 360)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, axis='y', alpha=0.2)

    fig.tight_layout()
    fig.savefig(OUT_DIR / 'land_normalized.png', dpi=150, bbox_inches='tight',
                facecolor=DARK_BG)
    plt.close(fig)
    print("  → land_normalized.png")
    return normalized, land_area

def generate_land_random_points(df, n=10000):
    """
    Generate random points on land by sampling from the convex hull of existing sites.
    Uses existing site locations to define 'land' — resample with replacement,
    then add small jitter to avoid exact duplication.
    """
    # Sample from existing site locations with jitter
    rng = np.random.default_rng(42)
    idx = rng.choice(len(df), size=n, replace=True)
    lons = df['lon'].values[idx] + rng.normal(0, 0.5, n)
    lats = df['lat'].values[idx] + rng.normal(0, 0.5, n)
    # Wrap longitude
    lons = ((lons + 180) % 360) - 180
    offsets = (lons - GIZA_LON) % 360
    return offsets

def plot_kde_comparison(df, offsets):
    """Plot 3: KDE comparison — sites vs random land points."""
    x_grid = np.linspace(0, 360, 1000)

    # Site KDE
    site_kde = stats.gaussian_kde(offsets, bw_method=5/offsets.std())
    site_density = site_kde(x_grid)

    # Random land points KDE
    land_offsets = generate_land_random_points(df, n=50000)
    land_kde = stats.gaussian_kde(land_offsets, bw_method=5/land_offsets.std())
    land_density = land_kde(x_grid)

    fig, ax = plt.subplots(figsize=(12, 5))

    ax.fill_between(x_grid, site_density, alpha=0.4, color=ACCENT, label='Archaeological sites')
    ax.plot(x_grid, site_density, color=ACCENT, linewidth=1.5)

    ax.fill_between(x_grid, land_density, alpha=0.3, color='#4ecdc4', label='Random land points')
    ax.plot(x_grid, land_density, color='#4ecdc4', linewidth=1.5, linestyle='--')

    add_grid_lines(ax)
    ax.set_xlabel('Longitude offset from Giza (°)')
    ax.set_ylabel('Density')
    ax.set_title('KDE comparison: archaeological sites vs. random land baseline')
    ax.set_xlim(0, 360)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, axis='y', alpha=0.2)

    fig.tight_layout()
    fig.savefig(OUT_DIR / 'kde_comparison.png', dpi=150, bbox_inches='tight',
                facecolor=DARK_BG)
    plt.close(fig)
    print("  → kde_comparison.png")
    return site_density, land_density, x_grid

def plot_residual(normalized, land_area):
    """Plot 4: Residual after land-area correction."""
    # Expected density = mean of normalized (flat if no pattern)
    valid = normalized[normalized > 0]
    expected = np.mean(valid) if len(valid) > 0 else 0
    residual = normalized - expected

    fig, ax = plt.subplots(figsize=(12, 5))

    colors = [ACCENT if r > 0 else ACCENT2 for r in residual]
    ax.bar(np.arange(360) + 0.5, residual, width=1.0,
           color=colors, edgecolor='none', alpha=0.85)

    ax.axhline(0, color=DARK_FG, linewidth=0.5, alpha=0.5)
    add_grid_lines(ax)

    ax.set_xlabel('Longitude offset from Giza (°)')
    ax.set_ylabel('Residual (observed − expected)')
    ax.set_title('Residual site density after land-area correction')
    ax.set_xlim(0, 360)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, axis='y', alpha=0.2)

    fig.tight_layout()
    fig.savefig(OUT_DIR / 'residual.png', dpi=150, bbox_inches='tight',
                facecolor=DARK_BG)
    plt.close(fig)
    print("  → residual.png")
    return residual, expected

def compute_summary_stats(offsets, counts, normalized, residual, expected):
    """Compute summary statistics and check 36° periodicity."""
    results = {
        "total_sites": int(len(offsets)),
        "giza_longitude": GIZA_LON,
        "grid_interval": GRID_INTERVAL,
        "bin_width_degrees": 1,
    }

    # Check density at 36° multiples vs elsewhere
    grid_bins = [int(g) % 360 for g in GRID_LINES]

    # Raw counts at grid lines (±2° window)
    grid_counts = []
    non_grid_counts = []
    for b in range(360):
        near_grid = any(abs(b - g) <= 2 or abs(b - g - 360) <= 2 for g in grid_bins)
        if near_grid:
            grid_counts.append(float(counts[b]))
        else:
            non_grid_counts.append(float(counts[b]))

    results["raw_counts"] = {
        "mean_at_grid_lines": float(np.mean(grid_counts)),
        "mean_away_from_grid": float(np.mean(non_grid_counts)),
        "ratio": float(np.mean(grid_counts) / np.mean(non_grid_counts)) if np.mean(non_grid_counts) > 0 else None,
    }

    # Normalized density at grid lines
    grid_norm = []
    non_grid_norm = []
    for b in range(360):
        if normalized[b] == 0:
            continue
        near_grid = any(abs(b - g) <= 2 or abs(b - g - 360) <= 2 for g in grid_bins)
        if near_grid:
            grid_norm.append(float(normalized[b]))
        else:
            non_grid_norm.append(float(normalized[b]))

    results["normalized_density"] = {
        "mean_at_grid_lines": float(np.mean(grid_norm)) if grid_norm else None,
        "mean_away_from_grid": float(np.mean(non_grid_norm)) if non_grid_norm else None,
        "ratio": float(np.mean(grid_norm) / np.mean(non_grid_norm)) if grid_norm and non_grid_norm and np.mean(non_grid_norm) > 0 else None,
    }

    # Residual at grid lines
    grid_resid = [float(residual[int(g) % 360]) for g in GRID_LINES]
    results["residual_at_grid_lines"] = {
        "values": {f"{int(g)}°": float(residual[int(g) % 360]) for g in GRID_LINES},
        "mean_residual": float(np.mean(grid_resid)),
        "expected_if_null": 0.0,
        "overall_residual_std": float(np.std(residual[residual != 0])) if np.any(residual != 0) else 0,
    }

    # Is the mean residual at grid lines significantly different from 0?
    valid_residuals = residual[residual != 0]
    if len(valid_residuals) > 0:
        z_score = float(np.mean(grid_resid) / (np.std(valid_residuals) / np.sqrt(len(grid_resid))))
        p_value = float(2 * (1 - stats.norm.cdf(abs(z_score))))
        results["residual_at_grid_lines"]["z_score"] = z_score
        results["residual_at_grid_lines"]["p_value"] = p_value

    # Fourier analysis: is there a peak at period=36°?
    fft_vals = np.fft.rfft(counts)
    freqs = np.fft.rfftfreq(360, d=1.0)  # cycles per degree
    power = np.abs(fft_vals)**2
    # Period of 36° = frequency of 1/36 ≈ 0.02778
    target_freq = 1.0 / 36
    freq_idx = np.argmin(np.abs(freqs - target_freq))

    results["fourier_analysis"] = {
        "target_period_degrees": 36,
        "power_at_36deg": float(power[freq_idx]),
        "mean_power": float(np.mean(power[1:])),  # exclude DC
        "power_ratio": float(power[freq_idx] / np.mean(power[1:])),
        "rank_among_all_frequencies": int(np.sum(power[1:] > power[freq_idx]) + 1),
        "total_frequencies": int(len(power) - 1),
    }

    # Continental peaks
    peak_bins = np.argsort(counts)[-10:][::-1]
    results["top_10_peak_offsets"] = [
        {"offset_deg": int(b), "count": int(counts[b]),
         "approx_longitude": float((b + GIZA_LON) % 360 - 180 if (b + GIZA_LON) % 360 > 180 else (b + GIZA_LON) % 360)}
        for b in peak_bins
    ]

    results["conclusion"] = (
        "The raw histogram shows peaks corresponding to continental land masses. "
        "After land-area normalization, no significant periodic signal at 36° intervals "
        "emerges above the noise floor." if results.get("residual_at_grid_lines", {}).get("p_value", 1) > 0.05
        else "After land-area normalization, there is a statistically significant signal "
        "at 36° intervals — further investigation warranted."
    )

    return results

def main():
    setup_style()

    print("Loading data...")
    df = load_data()
    offsets = df['lon_offset'].values

    print("\nPlot 1: Raw histogram")
    counts, bin_edges = plot_raw_histogram(offsets)
    # counts has 360 values (1° bins)

    print("\nPlot 2: Land-area normalized")
    normalized, land_area = plot_land_normalized(offsets, counts)

    print("\nPlot 3: KDE comparison")
    site_density, land_density, x_grid = plot_kde_comparison(df, offsets)

    print("\nPlot 4: Residual")
    residual, expected = plot_residual(normalized, land_area)

    print("\nComputing summary statistics...")
    results = compute_summary_stats(offsets, counts, normalized, residual, expected)

    with open(OUT_DIR / 'results.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("  → results.json")

    print(f"\nDone. All outputs in {OUT_DIR}/")
    print(f"\nKey findings:")
    print(f"  Raw count ratio (grid/non-grid): {results['raw_counts']['ratio']:.3f}")
    if results['normalized_density']['ratio']:
        print(f"  Normalized ratio (grid/non-grid): {results['normalized_density']['ratio']:.3f}")
    print(f"  Fourier 36° power rank: {results['fourier_analysis']['rank_among_all_frequencies']}/{results['fourier_analysis']['total_frequencies']}")
    if 'p_value' in results['residual_at_grid_lines']:
        print(f"  Residual z-score: {results['residual_at_grid_lines']['z_score']:.3f}")
        print(f"  Residual p-value: {results['residual_at_grid_lines']['p_value']:.4f}")

if __name__ == '__main__':
    main()
