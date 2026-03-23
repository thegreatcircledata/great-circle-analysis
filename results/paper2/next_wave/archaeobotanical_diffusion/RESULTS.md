# Study 8: Archaeobotanical Diffusion Rate Test

**Overall Verdict: COUNTER-EVIDENCE**

At least one dataset shows perpendicular spread to GC.

---

## Old World: Wheat-Barley-Lentil Package

**Verdict: NULL**

Old World mean offset=41.9°, V-test p=0.0000, percentile=24.6%. No significant alignment with GC corridor.

### Phase 1: Spread Bearings

- 16 consecutive pairs sorted chronologically

| From | To | Offset (deg) | Distance (km) | Time (yrs) |
|------|----|-------------|---------------|------------|
| Ohalo II | Tell Abu Hureyra | 51.3 | 436.1 | 10000 |
| Tell Abu Hureyra | Jericho | 56.1 | 522.0 | 2000 |
| Jericho | Çayönü | 59.3 | 806.2 | 500 |
| Çayönü | Ali Kosh | 40.8 | 960.3 | 1000 |
| Ali Kosh | Çatalhöyük | 24.6 | 1460.5 | 100 |
| Çatalhöyük | Jarmo | 12.9 | 1106.3 | 400 |
| Jarmo | Knossos | 3.1 | 1790.4 | 0 |
| Knossos | Franchthi Cave | 61.1 | 297.7 | 0 |
| Franchthi Cave | Mehrgarh | 8.1 | 4186.4 | 0 |
| Mehrgarh | Nabta Playa | 15.9 | 3751.7 | 0 |
| Nabta Playa | Karanovo | 84.4 | 2276.0 | 800 |
| Karanovo | Lepenski Vir | 47.9 | 363.8 | 200 |
| Lepenski Vir | Linearbandkeramik (avg) | 48.0 | 1088.0 | 500 |
| Linearbandkeramik (avg) | Fayum A | 61.6 | 2863.7 | 300 |
| Fayum A | Merimde | 72.6 | 101.2 | 400 |
| Merimde | Sarazm | 21.8 | 3460.7 | 1300 |

### Phase 2: Circular Statistics

- Mean angular offset: 41.85 deg
- Median angular offset: 47.95 deg
- Std: 24.56 deg
- V-test vs 0 deg (along GC): V=10.9068, p=0.000003
- V-test vs 90 deg (perp GC): V=9.7943, p=0.000010
- Rayleigh test: Z=13.4304, p=0.000001

### Phase 3: Speed Analysis

- Mean speed along GC: 2.3108 km/yr
- Mean speed perp to GC: 2.0331 km/yr
- Median speed along GC: 1.0213 km/yr
- Median speed perp to GC: 1.1698 km/yr
- Wilcoxon signed-rank (along > perp): stat=32.0, p=0.715332

### Phase 4: Null Model

- Alison GC mean offset: 41.85 deg
- Random GCs mean: 45.04 deg (std=4.09)
- 5th-95th percentile: 38.24 - 51.02 deg
- Alison percentile rank: 24.60%

---

## New World: Maize Diffusion

**Verdict: COUNTER-EVIDENCE**

New World mean offset=60.5° (>60°). Maize diffusion is predominantly perpendicular to GC, counter to corridor hypothesis.

### Phase 1: Spread Bearings

- 7 consecutive pairs sorted chronologically

| From | To | Offset (deg) | Distance (km) | Time (yrs) |
|------|----|-------------|---------------|------------|
| Tehuacan Valley | Guilá Naquitz | 83.7 | 195.9 | 750 |
| Guilá Naquitz | Oaxaca | 46.4 | 34.2 | 250 |
| Oaxaca | Panama (Aguadulce) | 52.5 | 2007.9 | 1000 |
| Panama (Aguadulce) | Coastal Ecuador | 61.5 | 1138.8 | 1000 |
| Coastal Ecuador | Coastal Peru | 84.8 | 1192.1 | 1000 |
| Coastal Peru | SW United States | 84.0 | 6071.7 | 900 |
| SW United States | Eastern Woodlands | 10.7 | 2074.8 | 300 |

### Phase 2: Circular Statistics

- Mean angular offset: 60.49 deg
- Median angular offset: 61.45 deg
- Std: 27.16 deg
- V-test vs 0 deg (along GC): V=3.0653, p=0.013287
- V-test vs 90 deg (perp GC): V=5.5652, p=0.000849
- Rayleigh test: Z=5.7667, p=0.000620

### Phase 3: Speed Analysis

- Mean speed along GC: 1.3577 km/yr
- Mean speed perp to GC: 1.7331 km/yr
- Median speed along GC: 0.5442 km/yr
- Median speed perp to GC: 1.1873 km/yr
- Wilcoxon signed-rank (along > perp): stat=6.0, p=0.921875

### Phase 4: Null Model

- Alison GC mean offset: 60.49 deg
- Random GCs mean: 45.14 deg (std=14.04)
- 5th-95th percentile: 22.78 - 64.63 deg
- Alison percentile rank: 79.60%

---

## Parameters

- Pole: (59.682122, -138.646087)
- Monte Carlo trials: 1000
- Random seed: 42

## Criteria

- SUPPORTED: Mean angular offset < 30 deg AND V-test p < 0.05 AND Alison ranks in bottom 10% of random GCs
- NULL: Mean angular offset ~45 deg AND V-test not significant
- COUNTER-EVIDENCE: Mean angular offset > 60 deg
