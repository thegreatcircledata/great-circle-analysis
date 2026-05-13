# Land-constrained Monte Carlo baselines for corridor enrichment analysis

**Institute:** [Deep Time Research Institute](https://deeptime-research.org)
**Author:** Elliot Allan — elliot@deeptime-research.org — ORCID [0009-0008-8541-0944](https://orcid.org/0009-0008-8541-0944)
**Interactive globe:** [thegreatcircle.earth](https://thegreatcircle.earth) | **Article series:** [Substack](https://thegreatcircle.substack.com)

> ⚠️ This repository is under active revision. The analysis described here corresponds to the manuscript currently under review. Some statistics may differ from earlier preprint versions. Files in `archive/` retain pre-correction numbers and are kept as historical record only.

## Current Status

| Paper | Title | Target | Status | Preprint |
|---|---|---|---|---|
| Merged GC | Land-constrained Monte Carlo baselines for corridor enrichment analysis: correcting coastal jitter bias and testing a proposed transcontinental great circle alignment | Scientific Reports | Under review | [zenodo.19046175](https://doi.org/10.5281/zenodo.19046175) |
| Quantifying TCC | Quantifying tectonic–climatic circumscription: a reproducible spatial-statistical test at primary-state locations | Humanities and Social Sciences Communications | Under review | [zenodo.19240284](https://doi.org/10.5281/zenodo.19240284) |

Last updated: 2026-05-13.

## Summary

A great circle defined by its pole at 59.682122°N, 138.646087°W has been proposed as a locus of anomalous clustering among ancient monumental sites (Alison, c. 2001). This repository contains the code and data supporting a land-constrained Monte Carlo framework that corrects coastal jitter bias in corridor enrichment tests and applies it across eight independently compiled archaeological databases (~259,000 unique sites) to evaluate the alignment claim, together with a companion collinearity analysis of independent civilization origins along the same circle.

After land-constrained correction, aggregate site enrichment is fully explained by the circle's path through habitable low-latitude terrain (78.5th percentile of habitability-matched circles; not significant). The residual finding is a monument-settlement divergence — monuments cluster while contemporaneous settlements are depleted — that localises to a single 250-year construction window at the Memphis necropolis during the Egyptian Old Kingdom. South American, Anatolian, astronomical, and pre-Younger Dryas alternative hypotheses are explicitly falsified.

All headline statistics here are drawn from the verified numbers summary at `outputs/merged_paper_prep/SUMMARY.md` (2026-03-30, FINAL).

## Key Results

| Finding | Result | Source |
|---|---|---|
| Megalithic Portal enrichment (50 km, land-constrained) | Z = 8.77 (N = 10,000); 1.62× enrichment | `outputs/merged_paper_prep/portal_10k_land_constrained.json` |
| Pleiades ancient monument enrichment (50 km, LC) | Z = 6.74; 2.52× enrichment | `outputs/merged_paper_prep/land_constrained_all_databases.json` |
| Monument–settlement divergence (Pleiades, LC) | D = 9.65; monument Z = 6.74, settlement Z = −2.91 | `outputs/merged_paper_prep/SUMMARY.md` |
| Random-pole baseline | 0 / 10,000 random great circles exceed the observed divergence; max random D = 0.90 | `outputs/merged_paper_prep/SUMMARY.md` |
| Habitability-adjusted divergence | D = 10.44, 99.63rd percentile of 4,323 habitability-matched circles | `outputs/merged_paper_prep/SUMMARY.md` |
| XRONOS enrichment (LC) | Z = 5.91 (N = 1,000); 1.40× | `outputs/merged_paper_prep/land_constrained_all_databases.json` |
| p3k14c enrichment (LC) | Z = 2.09 (p = 0.02; marginal); 1.16× | `outputs/merged_paper_prep/land_constrained_all_databases.json` |
| Civilization collinearity (Paper 3) | 4 of 7 independent origins within 200 km; p = 0.00042 | `analysis/paper3/` |
| Historic England negative control | D = 0.0; null as geometry predicts | `results/historic_england_negative_control.json` |
| 108° pair hypothesis | Z = −8.01 — falsified | `outputs/merged_paper_prep/SUMMARY.md` |
| Cross-paper BH correction (q = 0.05) | [needs verification — pre-LC BH table awaiting re-run] | Merged paper §4 |

## Database Inventory

| Database | Sites | Role | LC Z (50 km) |
|---|---|---|---|
| Megalithic Portal | 61,913 | Primary | 8.77 (N = 10,000) |
| Pleiades Gazetteer | 34,470 total; 406 ancient monuments pre-2000 BCE | Divergence analysis | 6.74 (monuments); divergence D = 9.65 |
| p3k14c Radiocarbon | 36,693 sites from 173,946 dates | Temporal replication | 2.09 (N = 1,000; p = 0.02) |
| XRONOS | 28,127 sites from 305,400 records | Independent replication | 5.91 (N = 1,000) |
| DARE Roman Atlas | 29,760 | Classification replication | Divergence D(LC) = 5.94 |
| Peru Ministry of Culture | 17,465 | South American null test | — |
| LuwianSiteAtlas | 483 | Anatolian null (0/500 km) | — |
| Historic England | 20,026 | Negative control (UK only) | 0.0 |
| OpenStreetMap archaeological | 28,400+ | Excluded (post-2001 contamination risk) | — |
| Wikidata | 15,200+ | Excluded (Pleiades cross-linkage) | — |

Land-constrained Monte Carlo baselines reduce observed Z-scores by 34–73% versus standard coastal-jittered baselines across the four databases where the correction was computed. See `submission/` and `outputs/merged_paper_prep/SUMMARY.md` for full methodology.

## Repository Structure

```
analysis/             Python analysis scripts
  core/                 Great circle enrichment, settlement baseline
  validation/           Split-sample, divergence 10K, cross-validation
  decomposition/        Hemisphere, continental, habitability
  robustness/           Data quality, KDE, stratified Monte Carlo
  databases/            DARE, Historic England, OSM, XRONOS tests
  statistical/          BH correction, multi-scale enrichment
  paper3/               Collinearity paper (PPM, LGCP, Thomas cluster)
data/                 Processed datasets (see download_data.sh for large files)
outputs/              Computed JSON outputs
  merged_paper_prep/    Verified headline-number source of truth (SUMMARY.md)
results/              Analysis outputs (JSON)
figures/              Publication-ready figures
submission/           Prior submission packages — historical; current submission at Scientific Reports
archive/              Earlier draft versions — contain pre-correction numbers; kept for record only
docs/                 Supplementary documentation
```

## Reproduce

```bash
# Install dependencies
pip install -r requirements.txt

# Download large data files not included in repo
bash download_data.sh

# Core analysis (merged paper)
python analysis/great_circle_test.py
python analysis/settlement_baseline_test.py
python analysis/directive1_split_sample.py

# Collinearity paper
python analysis/paper3/03_point_process_model.py
python analysis/paper3/12_independent_origins_test.py
python analysis/paper3/15_tectonic_diversity_test.py
Rscript analysis/paper3/13_kppm_cluster_process.R
Rscript analysis/paper3/19_lgcp_model.R
```

The Great Circle pole is defined at **59.682122°N, −138.646087°W**. This pole definition was not optimised or adjusted against any database in this study.

All analysis scripts read from `data/` and write to `results/` or `outputs/`.

## Data Sources

| Database | Source | License |
|---|---|---|
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

**Merged GC paper:**
Allan, E. (2026). *Land-constrained Monte Carlo baselines for corridor enrichment analysis: correcting coastal jitter bias and testing a proposed transcontinental great circle alignment.* Zenodo (v9.3, concept DOI). https://doi.org/10.5281/zenodo.19046175

**Paper 3 (Quantifying TCC):**
Allan, E. (2026). *Quantifying tectonic–climatic circumscription: a reproducible spatial-statistical test at primary-state locations.* Deep Time Research Institute. Zenodo (concept DOI). https://doi.org/10.5281/zenodo.19240284

## AI disclosure

Claude (Anthropic) was used as a computational research assistant for code generation, Monte Carlo simulation design, and manuscript drafting. All analyses were conceived, directed, and verified by the author, who takes full responsibility for the content.

## License

MIT (code) | CC-BY 4.0 (data, where original licenses permit). The corresponding author is affiliated with the Deep Time Research Institute, a registered 501(c)(3) nonprofit research organisation.
