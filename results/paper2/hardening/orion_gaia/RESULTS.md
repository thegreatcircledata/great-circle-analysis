# Test A1: Orion Gaia DR3 Proper Motion Upgrade

## Purpose
Confirm the Orion epoch sweep using Gaia DR3 proper motions in place of BSC values.

## Pyramid Axis
Bearing (Menkaure to Khufu): **39.60 deg**

## Results

| Metric | BSC | Gaia DR3 | Delta |
|--------|-----|----------|-------|
| Best-fit epoch | -11625 | -11625 | 0 yr (0.0%) |
| Min offset (deg) | 0.017 | 0.088 | 0.071 |
| Window <5 deg (yr) | 1225 | 1250 | 25 (2.0%) |
| Window <10 deg (yr) | 2550 | 2550 | 0 (0.0%) |
| p-value (5 deg) | 0.0475 | 0.0485 | |
| Offset @ -2560 (deg) | 33.493 | 33.555 | 0.062 |

## Procrustes Shape Stability

| Epoch | Distance |
|-------|----------|
| -11750 | 0.000616 |
| -10000 | 0.000553 |
| -5000 | 0.000003 |
| -2560 | 0.000143 |
| 0 | 0.000039 |
| 2000 | 0.000048 |

Max Procrustes distance: **0.000616** -- STABLE

## Verdict
BSC proper motions: **CONFIRMED** by Gaia DR3 upgrade.

- Best-epoch delta: 0 yr (0.0%)
- Window delta (5 deg): 25 yr (2.0%)
- Shape stability: STABLE (max Procrustes = 0.000616)
