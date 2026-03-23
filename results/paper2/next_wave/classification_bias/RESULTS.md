# Study 2: Pleiades Classification Bias Test

**Verdict: CLEAR**

Monument ratio is NOT significantly correlated with distance-to-GC at the regional level (Spearman p > 0.1), OR within-region divergence persists, indicating the pattern is not a classification artifact.

## Phase 1: Regional Classification Profile

Spearman rank correlation of monument_ratio with mean distance-to-GC:
- rho = -0.2848, p = 0.4250 (n = 10 regions)
- NO significant correlation — classification ratio not driven by GC proximity

| Region | Mon | Set | Both | Unc | Total | MonRat | MeanDist (km) |
|--------|----:|----:|-----:|----:|------:|-------:|--------------:|
| Egypt/Nubia | 253 | 616 | 25 | 393 | 1287 | 0.291 | 286.4 |
| Levant | 195 | 782 | 23 | 445 | 1445 | 0.200 | 335.8 |
| Anatolia | 481 | 2437 | 44 | 1772 | 4734 | 0.165 | 949.7 |
| Greece/Balkans | 744 | 2584 | 114 | 2095 | 5537 | 0.224 | 1020.7 |
| Italy | 1276 | 3831 | 134 | 1983 | 7224 | 0.250 | 1505.1 |
| North Africa (non-E) | 490 | 3093 | 37 | 960 | 4580 | 0.137 | 945.6 |
| Iberia | 225 | 1770 | 26 | 544 | 2565 | 0.113 | 1818.2 |
| Mesopotamia | 171 | 526 | 10 | 630 | 1337 | 0.245 | 479.3 |
| Iran/Central Asia | 90 | 582 | 6 | 360 | 1038 | 0.134 | 654.4 |
| South Asia | 32 | 303 | 2 | 261 | 598 | 0.096 | 812.7 |

## Phase 2: Within-Region Divergence

On-corridor: <50 km from GC | Off-corridor: >200 km from GC

| Region | On_Mon | On_Set | On_Ratio | Off_Mon | Off_Set | Off_Ratio | Fisher p |
|--------|-------:|-------:|---------:|--------:|--------:|----------:|---------:|
| Egypt/Nubia | 58 | 78 | 0.426 | 170 | 327 | 0.342 | 0.0707 |
| Levant | 12 | 39 | 0.235 | 154 | 578 | 0.210 | 0.7229 |
| Anatolia | 0 | 0 | n/a | 525 | 2481 | 0.175 | n/a |
| Greece/Balkans | 0 | 0 | n/a | 858 | 2698 | 0.241 | n/a |
| Italy | 0 | 0 | n/a | 1410 | 3965 | 0.262 | n/a |
| North Africa (non-E) | 0 | 4 | 0.000 | 527 | 3093 | 0.146 | 1.0000 |
| Iberia | 0 | 0 | n/a | 251 | 1796 | 0.123 | n/a |
| Mesopotamia | 0 | 3 | 0.000 | 164 | 466 | 0.260 | 0.5722 |
| Iran/Central Asia | 25 | 34 | 0.424 | 41 | 506 | 0.075 | 0.0000 |
| South Asia | 0 | 4 | 0.000 | 27 | 276 | 0.089 | 1.0000 |

Testable regions (both groups >= 5): 3
Significant divergence (p < 0.05): 1

## Phase 3: Unclassified Site Analysis

- On-corridor:  101 / 302 = 0.3344
- Off-corridor: 11308 / 32954 = 0.3431
- Fisher exact: OR = 0.962, p = 0.8076

## Decision Criteria

- **CLEAR**: monument_ratio NOT correlated with distance-to-GC (Spearman p > 0.1), OR divergence persists within-region
- **CONCERN**: Regional classification varies AND within-region sample sizes too small
- **FATAL**: monument_ratio strongly correlates with GC proximity AND within-region divergence disappears
