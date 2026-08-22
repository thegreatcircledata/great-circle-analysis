# Merged Paper 1+2: Verified Headline Numbers
## Generated: 2026-03-30 (FINAL — all computations complete)

---

## Land-Constrained Z-Scores (DEFINITIVE)

All three previously-uncorrected databases now have land-constrained baselines.
Sources: `land_constrained_all_databases.json`, `portal_10k_land_constrained.json`

| Database | N sites | Observed (50km) | Uncorr Z | **LC Z** | Enrichment | Z reduction |
|----------|---------|-----------------|----------|----------|------------|-------------|
| **Megalithic Portal** | 61,913 | 319 | 13.2 (N=200) | **8.77** ★ | **1.62×** | 33.6% |
| **p3k14c** | 36,693 | 187 | 7.77 (N=1k) | **2.09** | **1.16×** | 73.1% |
| **XRONOS** | 28,127 | 316 | 12.4 (N=1k) | **5.91** | **1.40×** | 52.3% |
| Pleiades (ancient mon) | 406 | 40 | 11.14 (N=1k) | **6.74** (N=10k) | **2.52×** | 39.5% |

★ Portal: N=10,000 definitive run. N=200 gave 9.04 (3% higher — expected convergence noise).

**Notes:**
- p3k14c 36,693 sites = SiteID-deduplicated full CSV. Paper 1 used 25,557 (lat/lon dedup); same 187 observed.
- p3k14c Z=2.09 is marginal (p=0.02). Report honestly.
- XRONOS at 1,000 trials with land constraint: Z=5.91 (vs Z=24.45 from original 200-trial unconstrained).
- Fallback rates all <0.25% — land mask is well-resolved.

---

## All Headline Numbers for the Merged Paper

### Primary Database Z-Scores (land-constrained where available)

| Database | N sites | LC Z (50km) | LC Enrichment | Source |
|----------|---------|-------------|---------------|--------|
| Megalithic Portal | 61,913 | **8.77** ★ | 1.62× | portal_10k_land_constrained.json |
| Pleiades (all) | 34,470 | 0.19 (uncorr) | 1.01× | pleiades_validation.json |
| Pleiades (pre-2000 BCE, 25km) | 778 | 10.68 (uncorr) | 3.44× | pleiades_validation.json |
| Pleiades (ancient monument) | 406 | **6.74** | 2.52× | coastal_bias_correction/05_all_results.json |
| p3k14c | 36,693 | **2.09** (p=0.02) | 1.16× | land_constrained_all_databases.json |
| XRONOS | 28,127 | **5.91** | 1.40× | land_constrained_all_databases.json |
| DARE | 29,760 | -4.77 (uncorr) | 0.69× | dare_analysis.json |
| Peru Ministry | 17,465 | — | — | peru_expanded_classification.json |
| CNSA/IPHAN | 27,190 | — | — | cnsa_brazil_divergence.json |
| Historic England | 20,026 | 0.0 | 0.0× | historic_england_analysis.json |
| Luwian | 483 | — (0/500km) | — | luwian_anatolia_test.json |
| SAAID | 2,429 | **UNVERIFIED** | — | **NO DATA FILE** |

### Monument-Settlement Divergence

| Database | D | Source |
|----------|---|--------|
| Pleiades (LC, N=10k) | **9.65** | coastal_bias_correction |
| Pleiades (uncorr, N=1k) | 9.98 | paper_v2.md |
| Megalithic Portal (non-EU) | 6.74 | paper_v2.md |
| p3k14c | 6.18 | paper_v2.md |
| DARE | 6.15 | dare_analysis.json |
| XRONOS (full) | +3.28 | xronos_great_circle_test.json |
| XRONOS (no Egypt) | -14.05 | xronos_great_circle_test.json |
| Peru Ministry | -1.25 | peru_expanded_classification.json |
| CNSA/IPHAN | 12.46 | cnsa_brazil_divergence.json |
| SAAID | 0.19 | **UNVERIFIED** |

### Divergence Robustness
- Habitability-adjusted: D = 10.44, 99.63rd percentile of 4,323 circles
- 0/10,000 random circles exceeded D = 9.91 (unconstrained baseline); max random D = 0.90
  - LC Pleiades D = 9.65 is also exceeded by 0/10,000 (since 9.65 >> 0.90 max)
- Egypt removed: enrichment Z = 4.51, divergence D = -0.64

### Temporal
- Peak: 3000-2500 BCE, D = 9.95 (Pleiades) / 3.83 (p3k14c)
- 2750-2500 BCE: Mon Z = 11.26, 31/38 pyramids on circle
- 500 BCE: D = 5.98 (settlement depletion Z = -5.36)

### Preservation
- Baseline D = 12.83; +58 buried: D = 13.10; +350 MC: 100% > D=2, mean = 10.21 ± 0.72

### Deep-Time
- Early Holocene: 4.52× enrichment, Z = 5.26 above FC mean, 0/9 FC match
- NERD: 1.10× (null)
- Pre-YD: Z = -0.47, 33rd percentile (null)

### Multi-Scale
- Peak enrichment: 9.23× at 7 km; Peak D: 2.18 at 18 km; Half-life: 32 km
- Corridor 1.43 ratio: UNVERIFIED (closest: 1.42 at minDate≤0)
- Ahramat angle: ~70° acute / ~110° obtuse at 20 km scale

### Circle Family
- 60 qualifying circles (1° grid); paper's 66/100k random unverified
- Alison D = 8.26 (rank 5/60, 93.2nd percentile) — self-consistent triple
- Best alternative: Athens+Persepolis+Varanasi D = 6.08

### Geophysical
- All 8 Z-scores verified; none survive Bonferroni (threshold 2.94)

### Null Results (all verified)
- 108° grid: falsified (Z = -8.01)
- Astronomical: null after Bonferroni
- Orientations: Rayleigh p = 0.75
- Trade: r = 0.04
- Luwian: 0/483
- South America: null ×3
- Agricultural p = 0.57, Lithology p = 0.74, Groundwater p = 0.37
- LGM: connected Z = 12.56, disconnected Z = 12.58

### Predictive
- Pre-2001: Z = 0.06; Post-2001: Z = **16.74** (NOT 17.24)
- 5 grid cells: 0.85% eliminates signal

### Paleoclimate
- Soreq r = -0.77, p = 0.025 (p3k14c, lag 500yr)
- BH-adjusted: p = 0.076 (marginal)

### Paper 3
- Tectonic p = 0.011, Circumscription p = 0.014
- Nulls: agricultural 0.57, lithology 0.74, groundwater 0.37

### Aggregate
- 228,937 site records across the 8 primary databases (sum of per-database counts after within-database deduplication — NOT a cross-database union; no cross-database dedup was performed)
- 643,420 raw entries/dates across those 8; state as "over 640,000" (use "entries" not "sites"; "550,000" in the paper was an undercount)
- Note: the 258,948 / ~637,000 totals in `database_inventory.json` cover TEN databases, including SAAID (2,429, excluded by the paper) and CNSA/IPHAN (27,582, never a primary database). That file is correct for what it measures; do not join it to the paper's eight-database framing.

---

## Critical Corrections for Merged Paper

1. **Portal Z = 8.77** (land-constrained, N=10,000 definitive) replaces Z = 13.3/25.85/9.04
2. **Post-2001 p3k14c Z = 16.74** (not 17.24)
3. **p3k14c site count**: Use 36,693 sites with Z = 2.09 (land-constrained) or 7.77 (uncorrected)
4. **All Z-scores must use land-constrained baselines** where computed
5. **XRONOS Z = 5.91** (land-constrained, N=1000) replaces Z = 24.45
6. **p3k14c Z = 2.09** (land-constrained, N=1000) replaces Z = 7.86; p=0.02 (marginal)
7. **Enrichment 2.52×** (Pleiades ancient monuments) replaces 5.05×
8. **D = 9.65** (land-constrained) replaces D = 12.78/9.98
9. **640,000+ is entries/dates** (not site records; the 8-database record count is 228,937, not ~259K — that total was the ten-database inventory; "550K" in paper was undercount)
10. **SAAID D = 0.19 is UNVERIFIED** — data file not found anywhere

---

## Files Produced

```
outputs/merged_paper_prep/
├── land_constrained_all_databases.json    ← Task 1 (includes definitive Portal N=200 + N=10k)
├── portal_10k_land_constrained.json       ← Portal N=10,000 DEFINITIVE ★
├── p3k14c_site_count_resolution.json      ← p3k14c 25,557 vs 36,693 resolution
├── saaid_divergence.json                   ← Task 2 (data not found)
├── corridor_width_profile.json             ← Task 3
├── circle_family_qualifying.json           ← Task 4
├── lgm_connectivity_zscores.json           ← Task 5
├── egypt_removal_clarification.json        ← Task 6
├── database_inventory.json                 ← Task 7
├── ahramat_crossing_angle.json             ← Task 9
├── paleoclimate_bh_correction.json         ← Task 10
├── circle_family_consistency.json          ← Task 11
├── corridor_width_profile.py               ← Task 3 script
├── task1_land_constrained_mc.py            ← Task 1 original (slow) script
├── task1_fast_land_constrained.py          ← Task 1 optimized script (N=200/1000)
├── task1_portal_10k.py                     ← Portal N=10,000 script
└── SUMMARY.md                              ← This file
```

---

## Remaining Unresolved Items

1. **SAAID data** — must be downloaded from Pandora/PACHAMAMA platform (no local file exists)
2. **10k random circle test with land-constrained baselines** — current test used unconstrained MC; LC Pleiades D=9.65 still exceeds max random D=0.90 so the conclusion holds, but a proper rerun would give a corrected D threshold for the comparison
3. **Corridor 1.43 ratio** — not reproducible from Pleiades data with stated parameters (closest: 1.42)
4. **Circle family 66/100k** — only 60 found on grid; random sampling output not located

All four items are documented in the JSON files above. Items 2–4 are minor: the scientific conclusions are unaffected by the missing precision. SAAID (item 1) is the only substantive gap.
