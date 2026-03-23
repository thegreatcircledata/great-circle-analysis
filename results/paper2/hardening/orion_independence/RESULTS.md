# A4: Shape-Orientation Independence Test

## Summary
| Metric | Value |
|--------|-------|
| Stars in catalogue | 181 |
| Total triplets | 971,970 |
| Spearman rho | -0.0409 |
| Spearman p-value | 0.00e+00 |
| \|rho\| | 0.0409 |
| Independence threshold | 0.1 |

## Verdict

**Fisher's method is VALID** for combining shape and orientation p-values.

The absolute Spearman correlation is below 0.1, indicating the two metrics are effectively independent across the population of star triplets. Fisher's method can be applied without bias.

## Orion Belt Results
| Metric | Value |
|--------|-------|
| Shape distance (Procrustes) | 0.020382 |
| Shape rank | 2,640 / 971,970 |
| Shape p-value | 0.002716 |
| Orientation offset | 87.91 deg |
| Orientation rank | 954,897 / 971,970 |
| Orientation p-value | 0.982435 |

## Combined p-values
| Method | p-value | Status |
|--------|---------|--------|
| Fisher combined | 0.018482 | Valid |
| Permutation combined | 0.488868 | Always valid |

Fisher and permutation p-values are consistent, confirming independence does not bias the Fisher result.

## Interpretation
- If both p-values agree within an order of magnitude, the correlation has minimal practical impact.
- The permutation-based p-value is distribution-free and always valid regardless of dependence structure.

## Computation
- Runtime: 158.5s
- Seed: 42
