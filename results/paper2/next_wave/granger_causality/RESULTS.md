# Study 5: Directional Granger Causality Along the Arc

## Verdict: **NULL**

Only 6% of pairs show unidirectional causality (p_null=0.2500), insufficient for diffusion claim

## Data Note

The p3k14c radiocarbon database has only **0 records within 100 km** of this
great circle (the arc passes through the Sahara, West Africa, South America, and
the Pacific where p3k14c has minimal coverage). Pleiades (648 sites
within 100 km) provides the primary dataset for this analysis, with
p3k14c run supplementarily using a wider 500 km corridor
(19 sites).

## Method

Tests whether archaeological activity propagates directionally along the
great-circle corridor using Granger causality on adjacent arc segments.

- **Corridor**: 100 km (Pleiades) / 500 km (p3k14c)
- **Arc**: 36 segments of 10 degrees each (reference: Giza at 0 degrees)
- **Time**: 500-year bins from 12000-0 BP (primary analysis)
- **Granger test**: F-test comparing restricted (AR(1)) vs unrestricted (AR(1) + lagged neighbor) model

## Phase 1: Arc-Time Matrix (Pleiades, 500yr)

- Matrix shape: 36 arc segments x 24 time bins
- Non-zero cells: 49
- Occupied segments: 9/36
- Total site counts: 648

## Phase 2: Autocorrelation

| Metric | Mean | Median |
|--------|------|--------|
| Temporal ACF (lag-1) | 0.291 | 0.316 |
| Spatial Moran's I | 0.194 | 0.110 |

## Phase 3: Granger Causality (Pleiades, 500yr bins)

| Category | Count | Fraction |
|----------|-------|----------|
| Forward only (i->i+1) | 1 | 2.9% |
| Reverse only (i+1->i) | 1 | 2.9% |
| Both directions | 2 | 5.7% |
| Neither | 31 | 88.6% |

- **Fraction unidirectional**: 0.057
- **Null model p-value**: 0.2500 (1,000 permutations)
- **Dominant direction**: increasing_arc

### p3k14c Supplementary (500km corridor)

| Metric | Value |
|--------|-------|
| Sites used | 19 |
| Occupied segments | 3/36 |
| Frac unidirectional | 0.000 |
| Null p-value | 1.0000 |

## Phase 4: Diffusion Rate

| Lag | Duration | Fraction Significant | Avg F |
|-----|----------|---------------------|-------|
| 1 | 500 yr | 8.6% (3/35) | 6.544 |
| 2 | 1000 yr | 5.7% (2/35) | 0.523 |
| 3 | 1500 yr | 5.7% (2/35) | 0.701 |

- **Best lag**: 1 time bins = 500 years
- **Arc segment width**: 1112 km
- **Implied diffusion rate**: 2.22 km/yr
- **Comparison**: Neolithic farming spread ~1 km/yr (Ammerman & Cavalli-Sforza)

## Phase 5: Sensitivity

| Bin Width | Frac Unidirectional | Null p-value | Dominant Dir |
|-----------|--------------------:|-------------:|--------------|
| 250 yr | 0.029 | 0.4800 | decreasing_arc |
| 500 yr | 0.057 | 0.2500 | increasing_arc |
| 1000 yr | 0.057 | 0.2310 | increasing_arc |

## Criteria Applied

- **STRONG**: >60% adjacent pairs show significant unidirectional Granger causality, null p < 0.01
- **SUGGESTIVE**: 30-60% significant with some directional bias
- **NULL**: <30% significant or no directional bias

## Output Files

- `arc_time_heatmap.png` -- Arc position vs time heatmap
- `granger_summary.png` -- Per-pair p-values and null distribution
- `diffusion_lag_plot.png` -- Significance by lag
- `sensitivity_binwidth.png` -- Sensitivity to bin width
- `results.json` -- Full numerical results
