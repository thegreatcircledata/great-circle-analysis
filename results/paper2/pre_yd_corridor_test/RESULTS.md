# Pre-Younger Dryas Corridor Test — Results

**Date:** 2026-03-20
**Databases:** p3k14c (170,117) + XRONOS (278,936) + ROAD v32 (18,755)
**After deduplication:** 94,181 dates (16,483 pre-YD)
**Corridor threshold:** 50km | Random circles: 1000
**Runtime:** 19.9s

---

## Verdict: CONVENTIONAL — No pre-YD anomaly detected

---

## Test 2 — THE KEY RESULT: Pre-YD Corridor Enrichment

| Metric | Observed | Random Mean ± SD | Z-score | Percentile |
|--------|----------|-----------------|---------|------------|
| Pre-YD (>12800 BP) | 10 | 127.1 ± 248.9 | **-0.47** | 32.8% |
| Deep pre-YD (>15000 BP) | 3 | 110.8 ± 220.2 | -0.49 | 21.8% |
| Late Glacial (12800-15000) | 7 | 16.4 ± 31.0 | -0.30 | 60.1% |
| Total (all periods) | 807 | 728.4 ± 931.5 | 0.08 | — |
| Habitability-matched pre-YD | 10 | 113.4 ± 104.7 | -0.99 | 11.7% |

## Test 3 — YD Disruption

| Metric | Corridor | Global | Z-score |
|--------|----------|--------|---------|
| Crash ratio (post-YD onset / pre-YD) | 1.500 | 1.316 | 0.11 |
| Recovery ratio | 5.500 | 4.481 | — |

## Test 4 — Regional Decomposition

| Egypt/Levant | 6 pre-YD | 443 post-YD | 1.3% pre-YD |
| Iran/Central Asia | 4 pre-YD | 83 post-YD | 4.6% pre-YD |
| Other | 0 pre-YD | 181 post-YD | 0.0% pre-YD |
| Peru/Amazon | 0 pre-YD | 44 post-YD | 0.0% pre-YD |
| North Africa/Mediterranean | 0 pre-YD | 46 post-YD | 0.0% pre-YD |

## Test 5 — Continuity

- On-corridor continuity rate: 2.11%
- Off-corridor continuity rate: 5.57%
- Continuity sites spanning the YD: 4 on-corridor, 1271 off-corridor

## Test 6 — Atlantis Regions

- **Atlantic west of Gibraltar**: 1112 dates, 328 pre-YD, ratio vs global = 1.97×
- **Richat Structure, Mauritania**: 17 dates, 0 pre-YD, ratio vs global = 0.00×
- **Persian Gulf basin**: 23 dates, 0 pre-YD, ratio vs global = 0.00×
- **Sunda Shelf, SE Asia**: 0 dates, 0 pre-YD, ratio vs global = 0.00×
- **Doggerland, North Sea**: 720 dates, 1 pre-YD, ratio vs global = 0.01×

## Sensitivity Checks

- **25km**: Z = -0.4534
- **100km**: Z = -0.3397
- **exclude_north_america**: Z = -0.4334

## Caveats

1. **Calibration:** Tests 1-6 use uncalibrated 14C ages. The YD-era calibration plateau may compress some dates.
2. **p3k14c North America coordinates:** Fuzzed to county centroids. Sensitivity test excluding NA: Z = -0.43
3. **Submerged regions:** Persian Gulf, Sunda Shelf, Doggerland have minimal data due to submersion, not necessarily absence of activity.
4. **Sampling bias:** The Natufian, PPNA, and Clovis are heavily sampled. Matched-corridor comparison controls for this globally but not regionally.
5. **Multiple comparisons:** 40 temporal bins tested; apply Bonferroni or FDR correction for individual bin significance.

## Files

- `test1_temporal_profile.json` — Enrichment time series
- `test2_pre_yd_enrichment.json` — Pre-YD Z-scores and random distributions
- `test3_yd_disruption.json` — YD crash ratio analysis
- `test4_regional_decomposition.json` — Regional breakdown
- `test5_continuity.json` — Site continuity across YD
- `test6_atlantis_regions.json` — Candidate region analysis
- `test7_spd.json` — Simplified SPD
- `figures/` — Visualizations
