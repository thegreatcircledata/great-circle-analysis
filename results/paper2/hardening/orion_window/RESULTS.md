# A2: Alignment Window Sensitivity Sweep

## Purpose
Test whether the combined shape + orientation p-value is robust
to the choice of alignment-window threshold.

## Method
- Pyramid axis bearing (Khufu to Menkaure): 219.61 deg
- Minimum offset between ground-mirrored Belt PA and axis: 0.05 deg at epoch -11800
- Swept Orion Belt PA across epochs -13000 to +2026 (step 25 yr)
- Precession computed via astropy FK5 equinox transformation
- E-W mirror applied for sky-to-ground mapping
- Offset computed mod 180 deg (undirected line comparison)
- Shape p = 0.0025 (from Study 12)
- Combined via Fisher's method (df=4)

## Results

| Threshold (deg) | Window (yr) |   Orient p |   Combined p |
|-----------------|-------------|------------|--------------|
|               1 |         250 |   0.009700 |  0.000281969 |
|               2 |         475 |   0.018431 |  0.000506166 |
|               3 |         725 |   0.028131 |  0.000742831 |
|               4 |        1000 |   0.038802 |  0.000993399 |
|               5 |        1225 |   0.047532 |    0.0011928 |
|               6 |        1500 |   0.058203 |    0.0014311 |
|               7 |        1750 |   0.067903 |   0.00164345 |
|               8 |        2025 |   0.078574 |   0.00187304 |
|               9 |        2250 |   0.087304 |   0.00205815 |
|              10 |        2375 |   0.092154 |   0.00216004 |
|              15 |        3250 |   0.126106 |   0.00285696 |
|              20 |        4775 |   0.185279 |   0.00401932 |
|              30 |        8425 |   0.326905 |   0.00662763 |
|              45 |       14800 |   0.574267 |    0.0108337 |

## Critical Thresholds
- Combined p first exceeds **0.01** at **45 deg**
- Combined p never exceeds **0.05** in tested range

## Pre-registration Argument
Archaeoastronomical precision benchmarks: Schaefer (1993) ~0.5 deg, Ruggles (1999) 1-2 deg, Aveni (2001) 2-3 deg. At 2 deg threshold combined p = 0.0005062; at 3 deg threshold combined p = 0.0007428; at 5 deg threshold combined p = 0.001193. Result is ROBUST: combined p < 0.01 at literature-standard 2-3 deg windows.

## Files
- `results.json` -- machine-readable results
- `window_sensitivity.png` -- combined p vs threshold plot
