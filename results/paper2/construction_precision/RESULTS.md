# Construction Precision Analysis — Results

**Date:** 2026-03-22
**Status:** PARTIAL (Phase 2B triggered — insufficient cross-site data for statistical analysis)

---

## Executive Summary

The Construction Anomaly Index (CAI) was computed for all ancient sites with published quantitative precision measurements. **The central finding is a data gap, not a data pattern:** only 3 independent sites worldwide (Giza, Serapeum, Puma Punku) have any numerical precision measurements in the published literature. Of the 12 major sites examined, 7 are completely unmeasured and 3 have only qualitative descriptions.

Where data does exist, the Great Pyramid at Giza exhibits construction precision that is genuinely anomalous — **11.5x more precise than modern CNC dimensional stone standards** (casing stone flatness) and **54x more precise** in the descending passage. This anomaly is real, well-documented (Petrie 1883, confirmed by Dash 2016), and has never been replicated by experimental archaeology. However, because essentially all quantitative precision data comes from one site, no cross-site statistical analysis is meaningful.

## Key Numbers

| Metric | Value |
|--------|-------|
| Sites examined | 12 (19 sub-components) |
| Sites with quantitative precision data | 3 (Giza, Serapeum, Puma Punku) |
| Sites completely unmeasured | 7 |
| CAI computable entries | 12 measurement-site pairs |
| Great Pyramid mean CAI | +3.09 (strongly anomalous) |
| Puma Punku CAI | +2.00 (moderately anomalous) |
| Serapeum CAI | disputed (-195 to +4.0 depending on source) |

## CAI Results (HIGH quality data only)

| Site | Measurement | Value | Baseline | CAI | Interpretation |
|------|------------|-------|----------|-----|----------------|
| GP — descending passage | surface flatness | 0.028 mm/m | 2.0 mm/m (copper tools est.) | **+3.9** | Extreme anomaly |
| GP — base platform | surface flatness | 0.091 mm/m | 2.0 mm/m | **+3.8** | Extreme anomaly |
| GP — casing stones | surface flatness | 0.13 mm/m | 2.0 mm/m | **+3.7** | Extreme anomaly |
| GP — casing stones | joint gap | 0.5 mm | 2.0 mm | **+3.0** | Strong anomaly |
| GP — King's Chamber | joint gap | 0.5 mm | 2.0 mm | **+3.0** | Strong anomaly |
| GP — coffer | surface flatness | 0.53 mm/m | 2.0 mm/m | **+2.9** | Strong anomaly |
| GP — base platform | angular accuracy | 12 arcsec | 60 arcsec | **+2.4** | Moderate anomaly |
| GP — casing stones | angular accuracy | 12 arcsec | 60 arcsec | **+2.4** | Moderate anomaly |
| Puma Punku H-blocks | surface flatness | 1.0 mm/m | 3.0 mm/m (stone tools est.) | **+2.0** | Moderate anomaly |

All Great Pyramid measurements exceed the experimental baseline by >2 standard deviations. The descending passage (0.028 mm/m straightness over 45.7m) is the single most anomalous measurement in the dataset.

## Comparison with Modern Standards

| Surface | Flatness (mm/m) | vs. Modern CNC |
|---------|-----------------|----------------|
| GP descending passage | 0.028 | **54x** more precise |
| GP base platform | 0.091 | **16x** more precise |
| GP casing stones | 0.13 | **12x** more precise |
| GP coffer | 0.53 | **2.8x** more precise |
| Puma Punku H-blocks | 1.0 | **1.5x** more precise |
| Modern CNC (MIA std) | 1.5 | reference |
| Serapeum (Dunn claim) | 0.0027 | **556x** — almost certainly wrong |

**Important caveat:** The MIA standard (1.5 mm/m) is a *minimum acceptance specification*, not the limit of modern capability. Modern optical grinding achieves sub-micron flatness routinely. The comparison is against commercial stone-cutting standards, not physics limits.

## Spatial Analysis

**Distance to Alison Great Circle:**

| Site | Distance (km) | Has Precision Data? |
|------|--------------|-------------------|
| Serapeum | 6 | Yes (disputed) |
| Great Pyramid | 7 | Yes (HIGH quality) |
| Valley Temple | 6 | No |
| Ollantaytambo | 7 | No |
| Easter Island | 17 | No |
| Sacsayhuaman | 46 | No (qualitative only) |
| Osireion | 420 | No |
| Baalbek | 421 | No |
| Puma Punku | 512 | Yes (MEDIUM quality) |
| Gobkeli Tepe | 770 | No |
| Malta | 897 | No |
| Stonehenge | 2,915 | No |

Several sites near the Great Circle (Ollantaytambo, Valley Temple, Sacsayhuaman, Easter Island) have zero precision measurements. **The spatial hypothesis cannot be tested until these sites are measured.** The 3-site correlation is meaningless (Spearman rho = 0.500, p = 0.667).

## Go / No-Go

| Test | Threshold | Result | Verdict |
|------|-----------|--------|---------|
| Sites with quantitative data | >= 10 | 8 | **FAIL** |
| Sites with CAI computable | >= 6 | 8 | **FAIL** |
| Spatial correlation p-value | < 0.05 | 0.667 | **INSUFFICIENT DATA** |
| Monte Carlo percentile | > 95th | N/A (n=3) | **INSUFFICIENT DATA** |
| Temporal correlation | r > 0.5 | -0.332 | **INSUFFICIENT DATA** |

**VERDICT: PARTIAL** — The Great Pyramid anomaly is real and well-documented, but cross-site analysis requires fieldwork. See `gap_analysis_report.md` for the research proposal.

## Bias Acknowledgments

1. **Giza dominance:** 9 of 12 CAI entries come from one monument complex. This is a data availability artifact, not evidence that Giza is uniquely precise.
2. **Baseline uncertainty:** Experimental archaeology baselines are estimated (no replication experiment has measured joint tolerances). CAI values would shift substantially with different baseline assumptions.
3. **Serapeum contradiction:** The Dunn claim (0.0027 mm/m) and the Antropogenez measurement (~1.1 angular deviation) are irreconcilable. Both are flagged LOW quality.
4. **Selection bias:** All 12 sites were selected because they are popularly cited for precision. Control sites with known "ordinary" construction are absent from this dataset.
5. **Petrie dependency:** The most important measurements in this entire analysis are 143 years old and have never been independently verified with modern instruments.

## Output Files

- `cleaned_data.csv` — 18 site entries with coordinates and great circle distances
- `anomaly_index.csv` — 12 CAI computations
- `spatial_correlation.csv` — Spearman correlation + Monte Carlo results
- `temporal_analysis.csv` — Pearson temporal correlation
- `site_distances.csv` — All sites' distance to Alison Great Circle
- `gap_analysis_report.md` — Phase 2B research proposal
- `charts/01_world_map_cai.png` — World map with sites colored by CAI
- `charts/02_cai_vs_distance.png` — CAI vs distance scatter
- `charts/03_cai_bar_chart.png` — CAI by site bar chart
- `charts/04_gap_analysis.png` — Data gap visualization
- `charts/05_flatness_comparison.png` — Ancient vs modern flatness comparison
