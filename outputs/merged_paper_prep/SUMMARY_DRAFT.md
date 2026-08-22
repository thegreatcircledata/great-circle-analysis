# Merged Paper 1+2: Verified Headline Numbers
## Generated: 2026-03-30 (DRAFT — awaiting Task 1 MC results)

---

## Status of All 11 Tasks

| Task | Status | Output File |
|------|--------|-------------|
| 1. Land-constrained MC (Portal, p3k14c, XRONOS) | **RUNNING** | `land_constrained_all_databases.json` |
| 2. SAAID divergence | DATA NOT FOUND | `saaid_divergence.json` |
| 3. Corridor width profile | COMPLETED | `corridor_width_profile.json` |
| 4. Circle family (66/100k) | PARTIALLY VERIFIED | `circle_family_qualifying.json` |
| 5. LGM connectivity Z-scores | VERIFIED | `lgm_connectivity_zscores.json` |
| 6. Egypt-removed Z=4.51 | VERIFIED | `egypt_removal_clarification.json` |
| 7. 550K aggregate | RESOLVED | `database_inventory.json` |
| 8. XRONOS 1000-trial | MERGED WITH TASK 1 | `land_constrained_all_databases.json` |
| 9. Ahramat crossing angle | CLARIFIED | `ahramat_crossing_angle.json` |
| 10. BH paleoclimate correction | COMPUTED | `paleoclimate_bh_correction.json` |
| 11. Circle family consistency | CLARIFIED | `circle_family_consistency.json` |

---

## Headline Numbers for the Merged Paper

### Pleiades (land-constrained, N=10,000 — DEFINITIVE)
- Monument (ancient) Z = **6.74**, enrichment = **2.52×**
- Settlement (ancient) Z = **-2.91**
- Monument-settlement divergence D = **9.65**
- Source: `outputs/coastal_bias_correction/05_all_results.json`

### Megalithic Portal (land-constrained — AWAITING COMPUTATION)
- **PENDING**: Task 1 MC running
- Uncorrected: Z = 13.3 (1,000 trials) or 25.85 (200 trials)
- Expected corrected Z: approximately 7-8 (if ~40% reduction similar to Pleiades)

### p3k14c (land-constrained — AWAITING COMPUTATION)
- **PENDING**: Task 1 MC running
- Uncorrected: Z = 7.86 (1,000 trials, n=25,557) or Z = 7.21 (n=36,693)
- Site count discrepancy: Z=7.86 used 25,557 sites, not 36,693

### XRONOS (land-constrained + 1,000 trials — AWAITING COMPUTATION)
- **PENDING**: Task 1 MC running
- Uncorrected: Z = 24.45 (200 trials, 28,127 sites)

### DARE
- Z = -4.77 (1,000 trials, 29,760 sites)
- Coastal correction: not needed (negative Z)
- Source: `outputs/new_databases/dare_analysis.json`

### Historic England (negative control)
- Z = 0.0 (20,026 sites, 0 within 50 km)
- Source: `outputs/new_databases/historic_england_analysis.json`

---

## Monument-Settlement Divergence

| Database | D | Source |
|----------|---|--------|
| Pleiades (pre-2000 BCE, land-constrained) | **9.65** | coastal_bias_correction/05_all_results.json |
| Pleiades (pre-2000 BCE, unconstrained, 1k MC) | 9.98 | paper/paper_v2.md |
| Megalithic Portal (non-EU) | 6.74 | paper/paper_v2.md |
| p3k14c | 6.18 | paper/paper_v2.md |
| DARE | 6.15 | dare_analysis.json |
| XRONOS (full) | +3.28 | xronos_great_circle_test.json |
| XRONOS (without Egypt) | -14.05 | xronos_great_circle_test.json |
| Peru Ministry | -1.25 | peru_expanded_classification.json |
| SAAID | 0.19 | **UNVERIFIED — no JSON source** |
| CNSA/IPHAN | 12.46 | cnsa_brazil_divergence.json |

### Robustness
- Habitability-adjusted D = 10.44 (99.63rd percentile of 4,323 circles)
- 0/10,000 random circles exceeded D = 9.91 (unconstrained baseline)
- Egypt removed: enrichment Z = 4.51, divergence D = -0.64

---

## Temporal Decomposition

| Period | D (Pleiades) | D (p3k14c) |
|--------|-------------|------------|
| 4000-3500 BCE | 0.13 | — |
| 3500-3000 BCE | -0.59* | — |
| **3000-2500 BCE** | **9.95** | **3.83** |
| 2500-2000 BCE | 2.87 | — |
| 2000-1500 BCE | 0.08 | — |
| 1500-1000 BCE | -1.31* | — |
| **500 BCE-0 CE** | **5.98** | — |

*Values marked with * are below the JSON min_sites_per_bin threshold; paper used a lower threshold.

### Peak Detail
- 2750-2500 BCE (250-yr): Monument Z = 11.26, 31/38 monuments on circle
- All 38 monuments are Egyptian Old Kingdom pyramids (Memphis necropolis)
- 500 BCE mechanism: settlement depletion (Z = -5.36), not monument enrichment

---

## Preservation Test
- Baseline D = 12.83
- +58 buried settlements: D = 13.10 (1000-MC) / 13.49 (200-MC)
- +350 MC settlements: D > 2 in 100% of 1,000 iterations, mean D = 10.21 ± 0.72

---

## Deep-Time
- Early Holocene (10,500-8,500 BP): 4.52× mean enrichment, Z = 5.26 above FC mean
- 0/9 Fertile Crescent circles match
- NERD replication: 1.10× (null)

## Pre-Younger Dryas
- 94,181 merged dates, 10 on corridor before 12,800 BP
- Z = -0.47, 33rd percentile (null)
- Continuity: 2.1% on-corridor vs 5.6% off-corridor

---

## Multi-Scale Profile
- Peak enrichment: **9.23×** at 7 km
- Peak divergence: **D = 2.18** at 18 km
- Half-life: **32 km**
- Corridor width 1.43 ratio: **UNVERIFIED** (closest match: 1.42 at minDate≤0, 5-15 km band)
- Ahramat crossing angle: ~70° acute / ~110° obtuse at 20 km scale

---

## Circle Family
- 60 qualifying circles on 1° grid (paper claims 66/100,000 random — unverified)
- Alison D = 8.26 (rank 5/60, 93.2nd percentile) — self-consistent triple
- Best alternative triplet: Athens+Persepolis+Varanasi, D = 6.08
- All 8 alternative D values verified

---

## Geophysical Scan
- All 8 Z-scores verified; none survive Bonferroni (threshold Z = 2.94)

---

## Null Results (all verified)
- 108° grid: falsified (Z = -8.01)
- Astronomical alignments: null after Bonferroni
- Monument orientations: Rayleigh p = 0.75
- Trade network: r = 0.04
- Luwian: 0/483 within 500 km
- South America: null on 3 databases
- Agricultural: p = 0.57, Lithology: p = 0.74, Groundwater: p = 0.37
- LGM connectivity: connected Z = 12.56, disconnected Z = 12.58 (no difference)

---

## Predictive Validation
- Pre-2001 p3k14c: Z = 0.06
- Post-2001 p3k14c: Z = **16.74** (NOT 17.24 as paper claims — CORRECTED)
- Five grid cells: 0.85% of data eliminates signal

---

## Paleoclimate
- Soreq Cave δ18O at 500-yr lag (p3k14c): r = -0.77, p = 0.025
- BH-adjusted (k=3 lags): p = 0.076 (marginal, does not survive at q=0.05)

---

## Paper 3 Cross-References
- Tectonic diversity: p = 0.011 (99.2nd percentile) ✓
- Circumscription: p = 0.014 (37.2% vs 10.9%) ✓
- Agricultural: p = 0.57, Lithology: p = 0.74, Groundwater: p = 0.37 ✓

---

## Aggregate
- 8 primary databases, 228,937 site records, over 640,000 raw entries/dates
- "640,000" refers to entries/dates, NOT site records
- 228,937 is a sum of per-database counts after within-database deduplication, not a cross-database union. The ~259,000 figure was the TEN-database inventory total (258,948), which includes SAAID and CNSA/IPHAN — neither is in the paper's corpus.

---

## Critical Items Still Outstanding

1. **Portal land-constrained Z** — MC running, result pending
2. **p3k14c land-constrained Z** — MC running, result pending
3. **XRONOS land-constrained Z** — MC running, result pending
4. **SAAID data** — must be downloaded and analyzed from scratch
5. **10k random circle test with land-constrained baselines** — should be rerun
6. **p3k14c site count / Z mismatch** — must choose: 25,557 sites + Z=7.86 OR 36,693 + Z=7.21
