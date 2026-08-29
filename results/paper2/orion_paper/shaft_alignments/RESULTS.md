# O3: Great Pyramid Shaft Alignment Test

> **Note added 2026-08-29.** The figures in this file are the output of the runs of
> 2026-03-22. Some of them were carried into Allan, E., 2026, *A comprehensive statistical
> evaluation of the Orion–Giza correlation*, Journal of Astronomical History and Heritage
> 29(2), 367–378, doi:10.3724/SP.J.1440-2807.2026.02.06, and are superseded by a correction
> in process at the journal. The figures affected are:
>
> - the catalogue size of 162 stars, in the N column of the Phase 4 table (lines 88–91).
>
> For this study the corrected figure is not simply the value used elsewhere in the paper.
> The correction record notes that the shaft study ran on a different star count from the
> rest of the analysis, one more than the catalogue size used elsewhere, both before and
> after the correction. See the DOI above.
>
> The accompanying `results.json` records the same run and is unchanged.
> The values above are left as the runs produced them.

## Overview

Tests whether the four Great Pyramid internal shafts align to specific
bright stars at ~2560 BCE, as claimed by Bauval (1994).

## Phase 1: Target Declinations

| Shaft | Angle | Direction | Target Dec | With Refraction | Uncertainty | Claimed Star |
|-------|-------|-----------|------------|-----------------|-------------|-------------|
| kings_south | 45.0° | south | -15.021° | -15.038° | ±0.500° | Alnitak (ζ Ori) |
| kings_north | 32.5° | north | +87.479° | +87.505° | ±0.500° | Thuban (α Dra) |
| queens_south | 39.5° | south | -20.521° | -20.541° | ±0.500° | Sirius (α CMa) |
| queens_north | 39.0° | north | +80.979° | +81.000° | ±0.500° | Kochab (β UMi) |

## Phase 2: Star Matching at 2560 BCE

### kings_south

Target altitude: 45.0° (south), target dec: -15.021°

| Star | Mag | Dec (2560 BCE) | Transit Alt | Offset | Claimed |
|------|-----|----------------|-------------|--------|---------|
| Alnilam | 1.69 | -15.057° | 44.964° | -0.036° |  |
| Larawag | 2.29 | -15.315° | 44.706° | -0.294° |  |
| Alnitak | 1.77 | -15.349° | 44.672° | -0.328° | YES |
| Sadalsuud | 2.91 | -14.619° | 45.402° | +0.402° |  |
| Mintaka | 2.23 | -14.576° | 45.445° | +0.445° |  |
| Dabih | 3.08 | -15.806° | 44.215° | -0.785° |  |
| Algedi | 3.57 | -13.444° | 46.577° | +1.577° |  |

### kings_north

Target altitude: 32.5° (north), target dec: +87.479°

| Star | Mag | Dec (2560 BCE) | Transit Alt | Offset | Claimed |
|------|-----|----------------|-------------|--------|---------|
| Thuban | 3.65 | +88.713° | 31.266° | -1.234° | YES |

### queens_south

Target altitude: 39.5° (south), target dec: -20.521°

| Star | Mag | Dec (2560 BCE) | Transit Alt | Offset | Claimed |
|------|-----|----------------|-------------|--------|---------|
| Lesath | 2.69 | -20.694° | 39.327° | -0.173° |  |
| Shaula | 1.63 | -20.717° | 39.304° | -0.196° |  |
| Tureis | 2.78 | -21.017° | 39.004° | -0.496° |  |
| Cursa | 2.79 | -21.113° | 38.908° | -0.592° |  |
| Iota Ori | 2.77 | -19.529° | 40.492° | +0.992° |  |
| Ascella | 2.60 | -21.536° | 38.485° | -1.015° |  |
| Menkar | 2.53 | -19.454° | 40.567° | +1.067° |  |
| Saiph | 2.09 | -21.904° | 38.117° | -1.383° |  |
| Kaus Australis | 1.85 | -22.173° | 37.848° | -1.652° |  |

### queens_north

Target altitude: 39.0° (north), target dec: +80.979°

| Star | Mag | Dec (2560 BCE) | Transit Alt | Offset | Claimed |
|------|-----|----------------|-------------|--------|---------|
| Kochab | 2.08 | +80.245° | 39.734° | +0.734° | YES |
| Pherkad | 3.05 | +80.228° | 39.751° | +0.751° |  |

## Phase 3: Exhaustive Comparison

How many bright stars (mag < 4.0) match each shaft within ±1°?

| Shaft | N matching | N total | p_shaft |
|-------|-----------|---------|--------|
| kings_south | 6 | 162 | 0.0370 |
| kings_north | 0 | 162 | 0.0000 |
| queens_south | 5 | 162 | 0.0309 |
| queens_north | 2 | 162 | 0.0123 |

## Phase 4: Combined Significance

### Fisher's Method

- Individual p-values: [0.037037, 0.006173, 0.030864, 0.012346]
- Chi² = 32.512 (df=8)
- Combined p = 0.000075

### Monte Carlo

- 10,000 trials, 4 random angles in [25.0°, 55.0°]
- Tolerance: ±1.0°
- p(all 4 shafts match a bright star) = 0.216300

| Shafts matching | Count | Fraction |
|-----------------|-------|----------|
| 0/4 | 0 | 0.0000 |
| 1/4 | 94 | 0.0094 |
| 2/4 | 2743 | 0.2743 |
| 3/4 | 5000 | 0.5000 |
| 4/4 | 2163 | 0.2163 |

## Phase 5: Epoch Sweep

Best-fit epoch per shaft (where claimed star's transit altitude
most closely matches shaft angle):

| Shaft | Claimed Star | Best Epoch (BCE) | Offset at Best | Offset at 2560 BCE |
|-------|-------------|------------------|----------------|--------------------|
| kings_south | Alnitak | 2500 | -0.0345° | N/A |
| kings_north | Thuban | 3225 | -0.0390° | N/A |
| queens_south | Sirius | 1900 | -0.0337° | N/A |
| queens_north | Kochab | 2350 | +0.0030° | N/A |

**Convergence:** Mean = 2494 BCE, Std = 476 yr, Spread = 1325 yr

## Interpretation

The shaft alignments do NOT show statistically significant correspondence
with claimed stellar targets beyond chance (MC p = 0.2163, Fisher p = 0.0001).
Multiple bright stars can match any given shaft angle at some epoch.

The four shafts do NOT converge tightly (spread = 1325 years).
