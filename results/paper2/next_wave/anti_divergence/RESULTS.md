# Study 9: Anti-Divergence Circle Search

## Overview
This study searches for the great circle that maximizes settlement clustering
while minimizing monument clustering -- the inverse of the Alison pattern.
By identifying the 'anti-circle', we can understand what geographic and
archaeological properties distinguish the Alison circle.

## Phase 1: Pole Scan Results

Scanned 684 pole candidates (lat 0-90, step 5; lon -180 to 180, step 10).

### Top 10 Anti-Divergence Poles

| Rank | Lat | Lon | Monuments | Settlements | Inv. Divergence |
|------|-----|-----|-----------|-------------|------------------|
| 1 | 10.0 | -30.0 | 0 | 72 | 72.00 |
| 2 | 40.0 | -60.0 | 2 | 80 | 40.00 |
| 3 | 0.0 | -20.0 | 0 | 39 | 39.00 |
| 4 | 0.0 | 160.0 | 0 | 39 | 39.00 |
| 5 | 35.0 | -20.0 | 0 | 36 | 36.00 |
| 6 | 65.0 | -100.0 | 1 | 33 | 33.00 |
| 7 | 70.0 | -120.0 | 0 | 27 | 27.00 |
| 8 | 70.0 | -90.0 | 1 | 26 | 26.00 |
| 9 | 50.0 | 140.0 | 2 | 50 | 25.00 |
| 10 | 70.0 | -130.0 | 1 | 25 | 25.00 |

## Phase 2: Anti-Circle Characterization

**Anti-circle pole:** (10.0, -30.0)

**Regions crossed:** Middle East / Central Asia, Africa, North America

- Settlements within 50km: 72
- Monuments within 50km: 0
- Inverse divergence: 72.00

### Monte Carlo Validation

Settlement enrichment: Z = -0.22, p = 0.2800 (200 trials)

Monument depletion: Z = -0.50, p(low) = 0.0000 (200 trials)

## Phase 3: Alison vs Anti-Circle Comparison

| Metric | Alison Circle | Anti-Circle |
|--------|--------------|-------------|
| Pole (lat, lon) | (59.682122, -138.646087) | (10.0, -30.0) |
| Monuments within 50km | 83 | 0 |
| Settlements within 50km | 122 | 72 |
| Mon/Set ratio | 0.6803 | 0.0000 |
| Set/Mon ratio | 1.4699 | 72.0000 |

## Phase 4: Orthogonality Test

Angular separation between Alison pole and anti-circle pole: **90.52 degrees**

The poles are approximately orthogonal (80-100 deg). This means the two
circles are nearly perpendicular, which would be a geometric constraint
difficult to explain by chance alone.

## Plots

- `pole_heatmap.png` -- Inverse divergence across all scanned poles
- `gc_comparison.png` -- Alison circle vs anti-divergence circle
- `divergence_comparison.png` -- Bar chart comparison of counts and ratios
