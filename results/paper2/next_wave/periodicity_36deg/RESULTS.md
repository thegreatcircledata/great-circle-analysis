# Study 6: 36-Degree Periodicity Validation on Data-Driven Clusters

**Date:** 2026-03-22
**Verdict:** FALSIFIED

## Overview

Previous analysis found 36-degree periodicity in inter-site arc distances
among 6 hand-picked sites along the proposed great circle (GC). This study
tests whether that periodicity persists in data-driven clusters identified
from 676 archaeological sites within 50 km of the GC.

## Phase 1: Cluster Detection

- Sites within 50 km of GC: 676
- KDE bandwidth: ~2 degrees
- Clusters identified (density > 2x mean): **6**
- Cluster arc positions (degrees from Giza, eastward):
  21.1, 102.8, 109.3, 292.4, 342.2, 357.0

## Phase 2: 36-Degree Rayleigh Test

Pairwise arc distances between 6 cluster centers (15 pairs)
tested for 36-degree periodicity.

| Metric | Value |
|--------|-------|
| Rayleigh Z | 0.3866 |
| p-value | 0.686501 |
| Mean resultant length (R) | 0.1605 |
| Mean direction | 135.2 deg |

## Phase 3: Multi-Period Scan

Tested periods: 10, 12, 15, 18, 20, 22.5, 24, 30, 36, 40, 45, 60, 72, 90, 120 degrees.
BH FDR correction applied.

| Period | Rayleigh Z | p (raw) | p (BH) | Survives? |
|--------|-----------|---------|--------|-----------|
| 10 | 0.3361 | 0.721269 | 0.800014 |  |
| 12 | 0.7936 | 0.459624 | 0.800014 |  |
| 15 | 0.3334 | 0.723136 | 0.800014 |  |
| 18 | 0.3477 | 0.713092 | 0.800014 |  |
| 20 | 1.0480 | 0.356704 | 0.800014 |  |
| 22.5 | 1.1707 | 0.315405 | 0.800014 |  |
| 24 | 0.3006 | 0.746679 | 0.800014 |  |
| 30 | 0.4490 | 0.645767 | 0.800014 |  |
| 36 | 0.3866 | 0.686501 | 0.800014 |  |
| 40 | 1.1757 | 0.313825 | 0.800014 |  |
| 45 | 0.4377 | 0.652940 | 0.800014 |  |
| 60 | 0.3026 | 0.745218 | 0.800014 |  |
| 72 | 0.8941 | 0.415956 | 0.800014 |  |
| 90 | 1.1859 | 0.310604 | 0.800014 |  |
| 120 | 0.0481 | 0.954525 | 0.954525 |  |

**Periods surviving BH correction:** None

## Phase 4: Monte Carlo Null Distribution

- 1000 random cluster sets (uniform on 0-360 deg arc, same N=6)
- Observed 36-deg Rayleigh Z: 0.3866
- Null distribution: mean = 0.9984, std = 1.3320
- Percentile rank of observed Z: 35.1%
- Monte Carlo p-value: 0.6490

## Verdict: FALSIFIED

The 36-degree periodicity is NOT statistically significant on data-driven clusters. The original finding appears to depend on hand-picked site selection.

## Outputs

- `results.json` — Full numerical results
- `kde_arc_positions.png` — KDE density plot of site arc positions with cluster centers
- `periodogram.png` — Multi-period scan bar chart
