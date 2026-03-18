# The Great Circle: Statistical Analysis of Ancient Monumental Site Distribution

**Paper:** [Zenodo (v3.0)](https://zenodo.org/records/19081718) | **Interactive Globe:** [thegreatcircle.earth](https://thegreatcircle.earth) | **Article Series:** [Substack](https://thegreatcircle.substack.com)

## Summary

A great circle defined by its pole at 59.682°N, 138.646°W has been proposed as a locus of anomalous clustering among ancient monumental sites (Alison, c. 2001). We tested this claim using five archaeological databases totaling over 150,000 sites.

**Key finding:** The overall site count is explained by the circle's path through habitable corridors. But within those corridors, monuments cluster at 5× the expected rate while settlements fall below random — a divergence that replicates across four independent databases, survives habitability adjustment (99.63rd percentile), and is not produced by any of 10,000 random great circles.

## Results at a Glance

| Finding | Result |
|---------|--------|
| Overall enrichment (habitability-adjusted) | 78.5th percentile — **not significant** |
| Monument-settlement divergence (Pleiades) | 9.98 Z-units — **significant** |
| Divergence replication (p3k14c) | 6.18 — **significant** |
| Divergence replication (DARE) | 6.15 — **significant** |
| Split-sample validation (100 splits) | Mean Z = 9.45, min = 7.31 — **passes** |
| 10,000 random circles matching divergence | 0/10,000 |
| Habitability-matched circles matching divergence | 16/4,323 (0.37%) |
| Negative control (Historic England) | Z = 0.0 — **null confirmed** |
| Benjamini-Hochberg FDR correction | 31/42 tests survive at q = 0.05 |
| Multi-scale peak significance | 10–20 km |

## Databases

| Database | Sites | Type | Role |
|----------|-------|------|------|
| Megalithic Portal | 61,913 | Community catalog | Primary |
| Pleiades Gazetteer | 34,470 | Academic gazetteer | Validation |
| p3k14c | 36,693 | Radiocarbon dates | Replication |
| DARE | 29,760 | Roman atlas | Validation |
| Historic England | 20,026 | Government register | Negative control |

## Repository Structure

```
paper/              Paper v3 (markdown + PDF)
data/               Processed datasets
  megalithic_portal/  KML files
  pleiades/           Pleiades CSV
  p3k14c/             Radiocarbon data
  dare/               DARE GeoJSON
  historic_england/   Scheduled monuments
analysis/           Python analysis scripts
  core/               Great circle test, settlement baseline
  validation/         Split-sample, divergence 10K, cross-validation
  decomposition/      Hemisphere, continental, habitability
  robustness/         Data quality checks, KDE, stratified MC
  databases/          DARE, Historic England, new database tests
  statistical/        BH correction, multi-scale enrichment
outputs/            Results organized by directive
website/            Interactive globe (thegreatcircle.earth)
```

## Reproduce

```bash
# Install dependencies
pip install numpy scipy pandas geopandas pyproj statsmodels

# Run core analysis
python analysis/core/great_circle_analysis.py

# Run settlement baseline
python analysis/core/settlement_baseline.py

# Run split-sample validation
python analysis/validation/directive1_split_sample.py
```

All analysis scripts read from `data/` and write to `outputs/`. The Great Circle pole is hardcoded: (59.682122°N, -138.646087°W).

## Citation

Allan, E. (2026). Statistical Analysis of Ancient Monumental Site Distribution Along a Proposed Great Circle: Evidence from Five Archaeological Databases and a Hemisphere Decomposition. Zenodo. https://zenodo.org/records/19081718

## License

MIT
