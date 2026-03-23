# ML Feature Importance Test — Results

**Date:** 2026-03-21
**Test:** Does "Distance to Great Circle" predict monument placement after geographic controls?

---

## Verdict

**WEAK: The Great Circle is largely a geographic proxy**

After controlling for measurable geographic features, the Great Circle adds minimal predictive power.

---

## Data Summary

| Metric | Value |
|--------|-------|
| Total sites | 21779 |
| Monuments | 4012 |
| Settlements | 17767 |
| Features | 10 |
| Corridor sites | 4158 |

Features used: dist_to_gc_km, elevation_m, dist_to_coast_km, dist_to_river_km, cloud_cover, seismic_pga, abs_latitude, longitude, site_density, cross_type_dist_km

---

## Model Performance (10-fold CV AUC)

| Model | Global AUC | Corridor AUC |
|-------|-----------|-------------|
| Random Forest | 0.8893 ± 0.0060 | 0.8834683125563616 |
| XGBoost | 0.8890 ± 0.0050 | 0.8744555751327063 |
| Logistic Regression | 0.8398 ± 0.0066 | 0.8349587192188268 |

---

## SHAP Feature Importance Ranking

### Global Model

| Rank | Feature | Mean |SHAP| |
|------|---------|-------------|
| 1 | cross_type_dist_km | 2.3442 |
| 2 | site_density | 0.5043 |
| 3 | elevation_m | 0.2712 |
| 4 | dist_to_coast_km | 0.2283 |
| 5 | abs_latitude | 0.1891 |
| 6 | dist_to_river_km | 0.1630 |
| 7 | longitude | 0.1532 |
| 8 | seismic_pga | 0.1201 |
| 9 | cloud_cover | 0.1093 |
| 10 | dist_to_gc_km | 0.0976 **<<<** |

**GC Distance Rank: #10 of 10**

### Corridor Model

**GC Distance Rank: #6 of 10**

---

## Ablation Test

| Metric | Global | Corridor |
|--------|--------|----------|
| AUC with GC | 0.8890 | 0.8744555751327063 |
| AUC without GC | 0.8888 | 0.8740545702932856 |
| Delta AUC | 0.0002 | 0.0004010048394207377 |

**Interpretation:** MINIMAL: GC is largely a geographic proxy

---

## Interpretation Table

| Scenario | Threshold | This Result |
|----------|-----------|-------------|
| GC ranks #1-3, ablation > 0.02 | Circle captures something geography doesn't | NO |
| GC ranks #4-6, ablation 0.005-0.02 | Circle adds modest information | NO |
| GC ranks #7+, ablation < 0.005 | Circle is a geographic proxy | YES |
| GC high SHAP in corridor only | Circle matters regionally | NO |

---

## Imputation Notes

- **soil_productivity**: DROPPED (100% missing)
- **water_table_m**: DROPPED (100% missing)
- **seismic_pga**: median imputed (53 values, 0.2%)

---

## Output Files

- `feature_matrix.csv` — all sites with all features
- `model_results.json` — full model results
- `shap_summary.png` — SHAP beeswarm plot (global)
- `shap_summary_corridor.png` — SHAP beeswarm plot (corridor)
- `ablation_test.json` — AUC with/without GC distance
- `partial_dependence.png` — monument probability vs GC distance
- `feature_interactions.png` — GC distance interactions with other features
- `logistic_coefficients.json` — interpretable model coefficients
- `feature_importance_table.json` — all importance metrics per feature
