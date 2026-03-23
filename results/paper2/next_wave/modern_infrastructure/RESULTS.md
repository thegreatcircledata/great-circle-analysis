# Study 10: Modern Infrastructure Echo Test

**Date:** 2026-03-22
**Seed:** 42

## Question
If the Alison great circle traces a deep habitability/connectivity corridor,
do modern cities — proxies for infrastructure — cluster near it more than
expected by chance?

## Method
- **City database:** 99 of the world's largest cities (population ≥ ~3M)
- **Metric:** Count cities within 50 / 100 / 200 / 500 km of each great circle
- **Null model:** 1000 random great circles (uniformly random poles, upper hemisphere)
- **Correlation:** Pearson r between archaeological-site density and city density
  along 360 equally-spaced points on the Alison GC
- **Land fraction:** Fraction of GC crossing land, Alison vs random

## Results

### City Proximity (Alison vs Random)

| Threshold | Alison | Random Mean ± SD | Z-score | Percentile |
|-----------|--------|------------------|---------|------------|
| 50 km     |      1 |   0.8 ± 1.0 |  +0.20 |  80.5% |
| 100 km    |      2 |   1.6 ± 1.5 |  +0.27 |  76.7% |
| 200 km    |      4 |   3.2 ± 2.4 |  +0.34 |  74.1% |
| 500 km    |     12 |   7.7 ± 4.8 |  +0.89 |  82.3% |

### Cities Within 200 km of Alison GC

| City                      | Lat      | Lon       | Pop(M) | Dist(km) |
|---------------------------|----------|-----------|--------|----------|
| Cairo                     |    30.04 |     31.24 |  21.3 |    12.9 |
| Kolkata                   |    22.57 |     88.36 |  15.1 |    86.1 |
| Alexandria                |    31.20 |     29.92 |   5.3 |   152.7 |
| Yangon                    |    16.87 |     96.20 |   5.3 |   176.7 |

### Ancient vs Modern Density Correlation

| Metric | Value |
|--------|-------|
| Pearson r (all 360 pts) | -0.0276 |
| p-value | 6.01e-01 |
| Pearson r (nonzero, n=114) | -0.1767 |
| p-value (nonzero) | 5.99e-02 |

### Land Fraction

| Circle | Land Fraction |
|--------|---------------|
| Alison | 0.1861 |
| Random mean ± SD | 0.2536 ± 0.0929 |
| Z-score | -0.73 |
| Percentile | 24.9% |

## Interpretation

The Alison great circle's city proximity is compared against a null distribution
of 1000 random great circles. Z-scores above +2 would indicate the Alison
circle passes near significantly more major cities than expected by chance,
supporting the "habitability corridor" hypothesis.

The ancient-vs-modern density correlation tests whether the same segments of the
arc that are rich in archaeological sites also tend to be near modern cities —
suggesting persistent geographic attractors across millennia.

## Outputs
- `results.json` — full numeric results
- `city_count_distributions.png` — histograms at each threshold
- `ancient_modern_density.png` — dual-axis density along the arc
- `land_fraction_distribution.png` — land fraction comparison
- `gc_cities_map.png` — map of GC with cities
