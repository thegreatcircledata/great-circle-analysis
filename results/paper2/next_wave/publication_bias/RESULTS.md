# Study 1: Publication & Excavation Density Bias Test

**Verdict: CLEAR**

XRONOS Mann-Whitney p=0.5834 (>0.05) and median ratio=1.000 (<1.5). No significant publication density bias detected. Note: p3k14c has zero on-corridor sites (<50 km) due to geographic sampling bias toward North America and Europe; the great circle passes through North Africa, South America, and the South Pacific where p3k14c has minimal coverage. XRONOS serves as the primary radiocarbon proxy.

## Phase 1: Radiocarbon Intensity

### p3k14c

- skipped (insufficient on-corridor coverage in p3k14c, on=0, off=12815)

- p3k14c geographic coverage does not reach the great circle corridor.

### XRONOS (primary)

- median ratio=1.000, p=0.5834, CI=[0.750, 1.000] (on=316, off=27402)


| Country | On (n) | Off (n) | Median Ratio | p-value |
|---------|--------|---------|-------------|---------|
| BR | 53 | 816 | 2.000 | 0.01554 |
| CL | 45 | 323 | 1.500 | 0.5899 |
| EG | 98 | 149 | 1.000 | 0.1692 |
| IL | 36 | 85 | 0.333 | 0.002588 |
| IR | 20 | 91 | 0.667 | 0.2979 |
| JO | 31 | 42 | 0.750 | 0.4075 |
| LY | 5 | 28 | 0.500 | 0.2339 |
| ML | 6 | 111 | 1.167 | 0.4185 |
| PE | 18 | 191 | 1.000 | 0.99 |


## Phase 2: Pleiades Attention Proxy

- median ratio=inf, p=0.1973, CI=[nan, nan] (on=302, off=32954)


## Phase 3: Monument vs Settlement Decomposition

### Pleiades

- **monument**: median ratio=inf, p=0.4342, CI=[nan, nan] (on=79, off=3914)
- **settlement**: median ratio=inf, p=0.8763, CI=[nan, nan] (on=118, off=17294)

### XRONOS

- **monument**: skipped (too few sites, on=0, off=337)
- **settlement**: median ratio=1.417, p=0.4098, CI=[0.830, 2.083] (on=56, off=4194)


## Phase 4: Distance-Decile Analysis

### XRONOS

- Spearman rho = -0.1636, p = 0.6515

| Decile | Mean dist (km) | Mean dates | Median dates | N sites |
|--------|---------------|-----------|-------------|---------|
| 1 | 459 | 12.54 | 4.0 | 2813 |
| 2 | 1283 | 9.56 | 3.0 | 2813 |
| 3 | 1830 | 7.43 | 3.0 | 2812 |
| 4 | 2123 | 8.33 | 4.0 | 2813 |
| 5 | 2445 | 10.38 | 4.0 | 2813 |
| 6 | 2824 | 8.19 | 3.0 | 2812 |
| 7 | 3065 | 7.28 | 3.0 | 2813 |
| 8 | 3278 | 8.49 | 4.0 | 2812 |
| 9 | 3531 | 7.29 | 3.0 | 2813 |
| 10 | 5707 | 19.83 | 6.0 | 2813 |

### Pleiades

- Spearman rho = -0.1879, p = 0.6032

| Decile | Mean dist (km) | Mean conn. | Median conn. | N sites |
|--------|---------------|-----------|-------------|---------|
| 1 | 221 | 0.23 | 0.0 | 3447 |
| 2 | 599 | 0.56 | 0.0 | 3447 |
| 3 | 873 | 0.41 | 0.0 | 3447 |
| 4 | 1016 | 0.38 | 0.0 | 3447 |
| 5 | 1138 | 0.37 | 0.0 | 3447 |
| 6 | 1353 | 0.48 | 0.0 | 3447 |
| 7 | 1587 | 0.80 | 0.0 | 3447 |
| 8 | 1797 | 0.33 | 0.0 | 3447 |
| 9 | 2117 | 0.23 | 0.0 | 3460 |
| 10 | 2849 | 0.28 | 0.0 | 3434 |

## Parameters

- Pole: (59.682122, -138.646087)
- On-corridor threshold: 50 km
- Off-corridor minimum: 200 km
- Monte Carlo trials: 1000
- Random seed: 42
