# Deep Time Atlas: Statistical Analysis of the Great Circle Hypothesis

A rigorous statistical test of whether ancient archaeological sites cluster near a specific great circle, as proposed by Jim Alison.

## Key Findings

| Dataset | Sites | Z @ 25 km | Z @ 50 km | Z @ 100 km | Z @ 200 km |
|---|---|---|---|---|---|
| UNESCO Cultural | 871 | 1.2 | 1.9 | 0.2 | 2.7 |
| Megalithic Portal | 61,870 | 23.2 | 23.7 | 19.9 | 18.5 |
| Merged (Portal + Supplement) | 61,913 | 24.6 | 25.9 | 23.3 | 21.3 |
| Pleiades (all) | 34,470 | 3.8 | 0.2 | 2.9 | 8.6 |
| Pleiades (pre-2000 BCE) | 778 | **10.7** | **6.5** | 4.0 | 5.3 |

- **108° angular separation hypothesis**: FALSIFIED (Z = -1.4)
- **Great circle proximity**: CONFIRMED (Z > 18 at all thresholds, merged dataset)
- **Independent validation**: Pleiades pre-2000 BCE sites show Z = 10.7 at 25 km
- **Temporal pattern**: Prehistoric sites (Z ≈ 21) show 2.5× stronger signal than later sites (Z ≈ 8)
- **Geographic coincidence ruled out:** Ancient monumental sites show Z = 11.83 (5× enrichment) while ordinary settlements in the same regions show Z = -0.95 (below random). The line selects for monuments, not geography.
- **Methodology**: Distribution-matched Monte Carlo with 200 independent trials

## Repository Structure

```
analysis/
  great_circle_test.py      # Standalone CLI tool — test any CSV against any great circle
  run_all_tests.py          # Full pipeline: parse, merge, run all 9 tests, generate results
data/
  supplement_sites.json     # 114 supplementary archaeological sites
  circle_coordinates.json   # 360 points along the great circle
  circle_coordinates.csv    # Same, in CSV format
  merged_sites.csv          # 61,913 sites (Portal + supplement, deduplicated)
  unesco_cultural_sites.json # 871 UNESCO Cultural Heritage sites
  README.md                 # Data source descriptions and download links
results/
  108_falsification.json    # Test 1: 108° angular separation (FALSIFIED)
  great_circle_escalation.json  # Test 2: Escalation across datasets and thresholds
  temporal_analysis.json    # Test 3: Prehistoric vs later sites
  multiple_circles.json     # Test 4: 1,000 random great circles comparison
  density_profile.json      # Test 5: Site density along the circle
  type_enrichment.json      # Test 6: Which site types cluster most
  age_analysis.json         # Test 7: Age comparison on-line vs off-line
  pleiades_validation.json  # Test 8: Independent validation with Pleiades data
  supplement_only.json      # Test 9: Hand-curated sites bias test
  settlement_baseline_test.json  # Test 10: Settlement vs monument baseline
  summary.json              # Combined summary of all results
docs/
  paper.md                  # Full methodology and results writeup
  blog.md                   # Accessible summary for general audience
```

## Quick Start

**Requirements**: Python 3.7+ (stdlib only — no pip install needed)

### Test any CSV against the great circle

```bash
python analysis/great_circle_test.py \
  --input data/merged_sites.csv \
  --output my_results.json \
  --trials 200
```

Your CSV must have `lat` and `lon` columns. Optional: `name`, `type`.

### Reproduce all results from scratch

```bash
# Full run (~2-3 hours, 200 Monte Carlo trials per test)
python analysis/run_all_tests.py --trials 200

# Fast run (~30 minutes, 50 trials, good for verification)
python analysis/run_all_tests.py --fast

# Skip Pleiades download (if you don't want to fetch the ~15MB file)
python analysis/run_all_tests.py --trials 200 --skip-pleiades
```

The script will:
1. Parse all KML files from the parent directory
2. Merge with supplementary sites (1 km dedup radius)
3. Download the Pleiades gazetteer (unless `--skip-pleiades`)
4. Run all 9 statistical tests
5. Generate summary JSON

## Methodology

### The Great Circle

Defined by its pole at **59.682122°N, 138.646087°W**. Every point on the circle is exactly one quarter of Earth's circumference (~10,007.5 km) from this pole. The circle passes through Egypt, Iran, northern India, Southeast Asia, the Pacific, Peru, and the North Atlantic.

### Distance Metric

```
gc_distance = |haversine(site, pole) - 10,007.5 km|
```

### Distribution-Matched Monte Carlo

The null hypothesis must account for the highly non-uniform geographic distribution of known archaeological sites (concentrated in Western Europe). For each trial:

1. For each real site, generate a random point by independently selecting a random site's latitude and longitude, then adding Gaussian jitter (σ = 2°)
2. Count random points within each distance threshold
3. Repeat 200 times to build a null distribution

**Z-score** = (observed - baseline_mean) / baseline_std

### Nine Tests

1. **108° Falsification** — Do sites form pairs at 108° more than chance? (No)
2. **Escalation Series** — Z-scores across 3 datasets × 4 thresholds
3. **Temporal Analysis** — Prehistoric vs later construction periods
4. **Multiple Circles** — Rank Alison's circle among 1,000 random circles
5. **Density Profile** — Site density at each degree along the circle
6. **Type Enrichment** — Which site types cluster most near the circle
7. **Age Analysis** — Mean/median age of on-line vs off-line sites
8. **Pleiades Validation** — Independent dataset from a separate project
9. **Supplement Bias Test** — Do hand-curated famous sites inflate the signal?
10. **Settlement Baseline** — Do ordinary settlements also cluster, or only monuments?

## Data Sources

| Source | Sites | License | Included |
|---|---|---|---|
| [Megalithic Portal](https://www.megalithic.co.uk/) | ~61,870 | Copyright (with permission for derived data) | Derived CSV only |
| Supplementary Sites | 114 | Public domain coordinates | Yes |
| [Pleiades Gazetteer](https://pleiades.stoa.org/) | ~34,000 | CC BY 3.0 | Auto-downloaded |
| [UNESCO WHC](https://whc.unesco.org/) | 871 | Public data | Yes |

Raw KML files from the Megalithic Portal are **not included** (copyright). The `merged_sites.csv` contains only coordinates and type classifications (derived data).

## Replication

To fully replicate from raw data:

1. Download KML files from the Megalithic Portal (one per site type)
2. Place them in the parent directory of this repo (or use `--kml-dir`)
3. Run `python analysis/run_all_tests.py --trials 200`

The script uses `random.seed(42)` for reproducibility. Results should match within Monte Carlo variance (~1-2 Z-score points between runs with different seeds).

## Citation

If you use this analysis or code, please cite:

```
[Anonymous]. (2026). Statistical Analysis of Ancient Monumental Site Distribution 
Along a Proposed Great Circle Alignment. The Great Circle Project.
https://doi.org/10.5281/zenodo.19046176
```

## Acknowledgments

- **Andy Burnham** and contributors to [The Megalithic Portal](https://www.megalithic.co.uk/) for the world's most comprehensive megalithic site database
- **Jim Alison** for the original great circle observation
- The **Pleiades** community (pleiades.stoa.org) for the ancient places gazetteer
- **UNESCO** World Heritage Centre for cultural site data
