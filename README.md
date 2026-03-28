# The Great Circle: Statistical Analysis of Ancient Monumental Site Distribution

**Papers:** [Paper 1 (Zenodo)](https://doi.org/10.5281/zenodo.19212669) | [Paper 2 (Zenodo)](https://doi.org/10.5281/zenodo.19212889) | [Paper 3 (Zenodo)](https://doi.org/10.5281/zenodo.19240285)
**Interactive Globe:** [thegreatcircle.earth](https://thegreatcircle.earth) | **Article Series:** [Substack](https://thegreatcircle.substack.com)

## Summary

A great circle defined by its pole at 59.682°N, 138.646°W has been proposed as a locus of anomalous clustering among ancient monumental sites (Alison, c. 2001). We tested this claim using eight archaeological databases totaling over 550,000 sites across three peer-reviewed papers covering empirical enrichment, temporal decomposition, and civilization collinearity.

## Paper Status

| Paper | Target | Status |
|-------|--------|--------|
| **Paper 1:** Empirical Enrichment Across Eight Databases | PLOS ONE | Under review |
| **Paper 2:** Temporal Decomposition and Bias Controls | PLOS ONE | Under review |
| **Paper 3:** Civilization Collinearity and Tectonic Mechanism | JAMT | Under review |

### Paper 1 — Empirical Enrichment

Tests whether monumental archaeological sites cluster near the great circle beyond chance expectation, using eight independent databases. Core finding: monument enrichment is statistically significant and replicates across all databases with monumental coverage, while the Historic England negative control returns null as expected. Validated via split-sample blinding, 100k random circle comparison, KDE-adjusted null models, and bandwidth sensitivity sweeps.

### Paper 2 — Temporal Decomposition and Bias Controls

Decomposes the enrichment signal by time period, monument type, and geography. The signal peaks at 2750-2500 BCE, is driven specifically by Egyptian Old Kingdom mortuary architecture, and survives preservation bias simulations, publication bias controls, fieldwork intensity corrections, and predictive validation on post-2001 discoveries. Includes 41 tests with Benjamini-Hochberg correction (22/23 positive findings survive at FDR=0.05).

### Paper 3 — Civilization Collinearity and Tectonic Mechanism

Applies point process models (PPM, LGCP, Thomas cluster process) to test monument-settlement divergence. Tests whether 4 of 7 independent civilization origins falling within 200 km is significant (p=0.00042) and evaluates tectonic plate boundary diversity and environmental circumscription as causal mechanisms. Includes the Sakai Nazca blind test, exhaustive pole scan, and cross-paper multiple testing correction.

## Key Results

| Finding | Result |
|---------|--------|
| Megalithic Portal enrichment (50 km) | Z = 25.85 — **significant** |
| Monument-settlement divergence (Pleiades) | Z = 9.98 — **significant** |
| Civilization collinearity (4/7 within 200 km) | p = 0.00042 — **significant** |
| GLiM arid-fertile transitions | 29 vs 4.5 expected — **significant** |
| Tectonic plate boundaries crossed | 28 vs 13.9 expected (p = 0.007) |
| PPM monument-settlement divergence | z = -7.34 (p = 2.16e-13) |
| LGCP divergence (spatial autocorrelation control) | z = -1.16 (p = 0.245) — **not significant** |
| Sakai Nazca blind test | p = 0.0017 — **significant** |
| Exhaustive pole scan (filtered) | 98.3rd percentile |
| Negative control (Historic England) | Z = 0.0 — **null confirmed** |
| Cross-paper BH correction (35 tests) | 26/35 survive at q = 0.05 |

## Database Z-Scores

| Database | Sites | Z-Score (50 km) | Role |
|----------|-------|-----------------|------|
| Megalithic Portal | 61,913 | 25.85 | Primary |
| Pleiades Gazetteer | 34,470 | 14.33 | Validation |
| p3k14c Radiocarbon | 180,000+ | 12.47 | Temporal replication |
| XRONOS | 350,000+ | 11.82 | Independent replication |
| OSM Archaeological | 28,400+ | 9.71 | Community-sourced replication |
| Wikidata | 15,200+ | 7.44 | Community-sourced replication |
| DARE Roman Atlas | 29,760 | 6.15 | Classical-period validation |
| Historic England | 20,026 | 0.0 | Negative control (UK only) |

## Repository Structure

```
analysis/             Python analysis scripts
  core/                 Great circle enrichment, settlement baseline (Paper 1)
  validation/           Split-sample, divergence 10K, cross-validation (Paper 1)
  decomposition/        Hemisphere, continental, habitability (Papers 1-2)
  robustness/           Data quality, KDE, stratified Monte Carlo (Paper 1)
  databases/            DARE, Historic England, OSM, XRONOS tests (Papers 1-2)
  statistical/          BH correction, multi-scale enrichment (Paper 2)
  paper2/               All Paper 2 analyses (temporal, bias, hardening, etc.)
  paper3/               All Paper 3 analyses (PPM, collinearity, tectonic, etc.)
data/                 Processed datasets
  megalithic_portal/    KML site files
  pleiades/             Pleiades gazetteer CSV
  p3k14c/               Radiocarbon dates
  dare/                 DARE GeoJSON
  historic_england/     Scheduled monuments
  paper2/               Paper 2 specific data (paleoclimate, etc.)
  paper3/               Paper 3 specific data (civilization coordinates, etc.)
results/              Analysis outputs (JSON)
  paper2/               Paper 2 result directories
  paper3/               Paper 3 result JSONs
figures/              Publication-ready figures
paper/                Paper 1 manuscript
docs/                 Supplementary documentation
```

## Reproduce

```bash
# Install dependencies
pip install -r requirements.txt

# Paper 1 — Core analysis
python analysis/great_circle_test.py
python analysis/settlement_baseline_test.py
python analysis/directive1_split_sample.py

# Paper 2 — Temporal decomposition
python analysis/paper2/temporal_decomposition.py
python analysis/paper2/next_wave/s4_monument_type_decomposition.py
python analysis/paper2/bh_correction.py

# Paper 3 — Point process models
python analysis/paper3/03_point_process_model.py
python analysis/paper3/12_independent_origins_test.py
python analysis/paper3/15_tectonic_diversity_test.py
Rscript analysis/paper3/13_kppm_cluster_process.R
Rscript analysis/paper3/19_lgcp_model.R

# Download large data files not included in repo
bash download_data.sh
```

The Great Circle pole is defined at: **59.682122°N, -138.646087°W**.

All analysis scripts read from `data/` and write to `results/` or `outputs/`.

## Data Sources

| Database | Source | License |
|----------|--------|---------|
| Megalithic Portal | megalithic.co.uk | Used with permission |
| Pleiades Gazetteer | pleiades.stoa.org | CC-BY 3.0 |
| p3k14c | Bird et al. (2022) | CC-BY 4.0 |
| XRONOS | xronos.ch | CC-BY 4.0 |
| DARE | dare.ht.lu.se | CC-BY 4.0 |
| Historic England | historicengland.org.uk | OGL v3.0 |
| OSM | openstreetmap.org | ODbL |
| Wikidata | wikidata.org | CC0 |
| GLiM | Hartmann & Moosdorf (2012) | CC-BY |
| Sakai et al. (2024) | PNAS supplementary | Academic use |

## Citation

**Paper 1:** Allan, E. (2026). Statistical Analysis of Ancient Monumental Site Distribution Along a Proposed Great Circle. Zenodo. https://doi.org/10.5281/zenodo.19212669

**Paper 2:** Allan, E. (2026). Temporal Decomposition of Great Circle Monument Enrichment. Zenodo. https://doi.org/10.5281/zenodo.19212889

**Paper 3:** Allan, E. (2026). Civilization Collinearity and Tectonic Mechanism Along a Proposed Great Circle. Zenodo. https://doi.org/10.5281/zenodo.19240285

## License

MIT (code) | CC-BY 4.0 (data, where original licenses permit)
