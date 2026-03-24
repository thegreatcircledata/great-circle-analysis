#!/usr/bin/env python3
"""
Generate rose diagram of Nazca line orientations with Great Circle bearing marked.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

# Data from Richter et al. 2021 (axial: 0-180°, mirrored to full circle)
BINS_CENTER = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115,
               125, 135, 145, 155, 165, 175]
BIN_COUNTS = [91, 126, 103, 100, 153, 127, 177, 116, 129, 76,
              104, 86, 165, 140, 172, 177, 124, 142]

GC_BEARING = 63.14
JUNE_SOLSTICE_SR = 65.6

# Mirror to full 360° for the rose diagram
full_centers = BINS_CENTER + [c + 180 for c in BINS_CENTER]
full_counts = BIN_COUNTS + BIN_COUNTS

# Convert to radians (compass: N=0, clockwise → math: E=0, counterclockwise)
# For polar plot: theta = 90° - azimuth (in radians)
theta = np.radians(90 - np.array(full_centers))
width = np.radians(10)

fig, ax = plt.subplots(1, 1, figsize=(10, 10), subplot_kw=dict(projection='polar'))

# Set up as compass (N at top, clockwise)
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)

# Plot bars
bars = ax.bar(np.radians(full_centers), full_counts, width=width,
              bottom=0, alpha=0.7, color='steelblue', edgecolor='navy', linewidth=0.5)

# Color the GC-aligned bins
for i, (center, count) in enumerate(zip(full_centers, full_counts)):
    if abs(center - GC_BEARING) < 5 or abs(center - (GC_BEARING + 180)) < 5:
        bars[i].set_facecolor('#ff4444')
        bars[i].set_alpha(0.9)

# Mark the Great Circle bearing
gc_rad = np.radians(GC_BEARING)
max_count = max(full_counts)
ax.annotate('', xy=(gc_rad, max_count * 1.25), xytext=(gc_rad, 0),
            arrowprops=dict(arrowstyle='->', color='red', lw=2.5))
ax.annotate('', xy=(np.radians(GC_BEARING + 180), max_count * 1.25),
            xytext=(np.radians(GC_BEARING + 180), 0),
            arrowprops=dict(arrowstyle='->', color='red', lw=2.5))

# Label GC bearing
ax.text(gc_rad, max_count * 1.35, f'GC: {GC_BEARING:.1f}°',
        ha='center', va='bottom', fontsize=11, fontweight='bold', color='red')
ax.text(np.radians(GC_BEARING + 180), max_count * 1.35, f'GC: {GC_BEARING + 180:.1f}°',
        ha='center', va='bottom', fontsize=11, fontweight='bold', color='red')

# Mark June solstice sunrise
ss_rad = np.radians(JUNE_SOLSTICE_SR)
ax.plot([ss_rad, ss_rad], [0, max_count * 1.15], '--', color='orange', lw=2, alpha=0.8)
ax.text(ss_rad + np.radians(3), max_count * 1.18, 'Jun solstice\nsunrise',
        ha='left', va='bottom', fontsize=9, color='orange', fontweight='bold')

# Mark Rayleigh mean direction
mean_dir = 154.9
mean_rad = np.radians(mean_dir)
ax.plot([mean_rad, mean_rad], [0, max_count * 1.1], ':', color='green', lw=2, alpha=0.8)
ax.text(mean_rad, max_count * 1.15, f'Mean: {mean_dir:.0f}°',
        ha='center', va='bottom', fontsize=9, color='green', fontweight='bold')

# Uniform expectation circle
expected = sum(BIN_COUNTS) / 18
circle_theta = np.linspace(0, 2*np.pi, 100)
ax.plot(circle_theta, [expected]*100, '--', color='gray', alpha=0.5, lw=1)
ax.text(np.radians(315), expected + 8, f'Uniform\n({expected:.0f}/bin)',
        fontsize=8, color='gray', ha='center')

# Title and labels
ax.set_title('Nazca Line Orientations vs Great Circle Bearing\n'
             '(2,308 lines, Richter et al. 2021)',
             fontsize=14, fontweight='bold', pad=30)
ax.set_rlabel_position(225)

plt.tight_layout()
plt.savefig('nazca_orientation_rose_diagram.png', dpi=150, bbox_inches='tight',
            facecolor='white')
print("Saved nazca_orientation_rose_diagram.png")

# ─── Also make a simpler half-rose (0-180° axial) ──────────────

fig2, ax2 = plt.subplots(figsize=(12, 5))

colors = ['#ff4444' if abs(c - GC_BEARING) < 5 else 'steelblue' for c in BINS_CENTER]
bars2 = ax2.bar(BINS_CENTER, BIN_COUNTS, width=8, color=colors, edgecolor='navy',
                linewidth=0.5, alpha=0.8)

ax2.axhline(y=expected, color='gray', linestyle='--', alpha=0.5, label=f'Uniform ({expected:.0f})')
ax2.axvline(x=GC_BEARING, color='red', linewidth=2.5, label=f'GC bearing ({GC_BEARING:.1f}°)')
ax2.axvline(x=JUNE_SOLSTICE_SR, color='orange', linewidth=2, linestyle='--',
            label=f'June solstice sunrise ({JUNE_SOLSTICE_SR}°)')
ax2.axvline(x=mean_dir, color='green', linewidth=2, linestyle=':',
            label=f'Rayleigh mean ({mean_dir:.0f}°)')

ax2.set_xlabel('Axial Orientation (degrees from North)', fontsize=12)
ax2.set_ylabel('Number of Lines', fontsize=12)
ax2.set_title('Nazca Line Orientations — Axial Distribution (0-180°)\n'
              '2,308 lines from Richter et al. 2021', fontsize=13, fontweight='bold')
ax2.legend(loc='upper right', fontsize=10)
ax2.set_xlim(-5, 185)
ax2.set_xticks(range(0, 181, 15))
ax2.grid(axis='y', alpha=0.3)

# Annotate the GC bin
gc_idx = BINS_CENTER.index(65)
ax2.annotate(f'Peak: {BIN_COUNTS[gc_idx]} lines\n(rank 1/18)',
             xy=(65, BIN_COUNTS[gc_idx]),
             xytext=(65, BIN_COUNTS[gc_idx] + 15),
             fontsize=9, fontweight='bold', color='red',
             ha='center', va='bottom',
             arrowprops=dict(arrowstyle='->', color='red', lw=1.5))

plt.tight_layout()
plt.savefig('nazca_orientation_histogram.png', dpi=150, bbox_inches='tight',
            facecolor='white')
print("Saved nazca_orientation_histogram.png")
