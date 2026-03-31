# Coastal Bias Correction: Corrected Statistics for Paper Revision

**Date:** 2026-03-24

## Summary of Corrections

The original Paper 1 distribution-matched Monte Carlo baseline (Approach A) did not constrain synthetic site coordinates to land. When jittered with sigma=2 degrees (~220 km), coastal and island sites frequently landed in the ocean, deflating the expected on-corridor count and inflating z-scores.

**Monument (ancient, pre-2000 BCE) — 50 km corridor:**
- Original unconstrained (N=1000): Z = 11.14, enrichment = 4.91x
- Land-constrained (N=10,000): Z = 6.74, enrichment = 2.52x
- Fallback rate: 0.00%

**Settlement (ancient, pre-2000 BCE) — 50 km corridor:**
- Land-constrained (N=10,000): Z = -2.91, enrichment = 0.5x


## Random Circles Cross-Validation (Approach B)

Approach B (random great circles, fixed sites) completely sidesteps the coastal bias issue because real sites stay at their actual coordinates.

- Monument (ancient) random circles Z = 3.45, enrichment = 12.1x


## Easter Island and Coastal Isolation


- full: 406 sites, Z = 7.08, enrichment = 2.54x

- no_easter_island: 406 sites, Z = 7.08, enrichment = 2.54x

- no_small_islands: 290 sites, Z = 6.91, enrichment = 2.59x

- no_coastal_egypt: 387 sites, Z = 3.98, enrichment = 2.11x

- no_all_coastal: 271 sites, Z = 4.15, enrichment = 2.15x


## MC Convergence


- N = 200: Z = 6.93

- N = 500: Z = 7.11

- N = 1000: Z = 7.08

- N = 5000: Z = 6.83

- N = 10000: Z = 6.74


## Recommendations for Paper Revision


Monument-settlement divergence under land-constrained correction: 9.65
Monument-settlement divergence under random circles: 2.55


These corrected values should replace the original Paper 1 statistics. Paper 2 should cite the corrected values and note the methodological improvement. The core finding — that monumental sites cluster near the Great Circle significantly more than contemporaneous settlements — should be evaluated based on the corrected z-scores above.
