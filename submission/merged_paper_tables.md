# Tables: Monument-settlement divergence along a proposed great circle

**Paper:** Monument-settlement divergence along a proposed great circle: Statistical evidence from eight archaeological databases with temporal decomposition and coastal-corrected baselines
**Author:** Elliot Allan
**Data:** github.com/thegreatcircledata/great-circle-analysis

All Z-scores are identified as land-constrained (LC) or uncorrected (UC) in each table. LC = land-constrained Monte Carlo baseline (Natural Earth 1:10M land mask, max 50 jitter attempts per site). UC = distribution-matched Monte Carlo without land constraint.

---

## Table 1. Database summary and overall enrichment statistics

Primary threshold 50 km. LC trials = number of Monte Carlo iterations under land-constrained baseline. "—" indicates value not computed.

| Database | N sites | Observed within 50 km | LC Z-score | Uncorr. Z | Enrichment (UC) | LC trials | LC applied | Notes |
|---|---|---|---|---|---|---|---|---|
| Megalithic Portal | 61,913 | — | **8.77** | 13.2 | — | 10,000 | Yes | ~65% UK/Ireland/France; circle avoids Europe; bias works against hypothesis |
| Pleiades (pre-2000 BCE) | 778 | — | **6.74** | 11.14 | 2.52× | 10,000 | Yes | Mediterranean focus; monument-settlement analysis uses 406 ancient monument subset |
| XRONOS | 28,127 | — | **5.91** | 12.4 | — | 1,000 | Yes | Independent of p3k14c; 305,400 raw records |
| p3k14c | 36,693 | 187 | **2.09** (p=0.02) | 7.77 | — | 1,000 | Yes | SiteID dedup; 173,946 raw dates; US/Canada coords obfuscated to county centroids |
| DARE | 29,760 | — | — | Negative | — | — | No | Negative overall enrichment; circle avoids Roman heartland |
| Peru Ministry of Culture | 17,465 | — | — | — | — | — | No | South America null; 1,733 monumental, 275 settlement |
| Historic England | 20,026 | 0 | 0.0 | 0.0 | — | — | No | Negative control; circle does not pass through England; zero sites within 500 km |
| LuwianSiteAtlas | 483 | 0 | — | — | — | — | No | Bronze Age western Anatolia 2000–1200 BCE; complete null |

**Notes on site counts:** 228,937 site records total across all 8 databases — a sum of per-database counts after within-database deduplication, not a cross-database union; over 640,000 raw entries, dates, and records (640,000 = entries/records, NOT site records). Paper 1 (Allan 2026a) used 25,557 p3k14c sites (lat/lon dedup); this study uses 36,693 (SiteID dedup, more complete); same 187 on-corridor sites in both.

**Habitability-adjusted result:** 78.5th percentile of 4,323 habitability-matched circles. Overall site-count enrichment is not significant after habitability adjustment.

**LC Z-score reductions from correction:** Portal 13.2→8.77 (33.6%); Pleiades 11.14→6.74 (39.5%); XRONOS 12.4→5.91 (52.3%); p3k14c 7.77→2.09 (73.1%).

---

## Table 2. Monument-settlement divergence across databases

D = Z_monument − Z_settlement. Primary LC result (Pleiades) shown in bold. LC baselines computed for Portal, p3k14c, and DARE (1,000 trials each; see lc_divergence_all_databases.json).

| Database | Monument Z | Settlement Z | D | Correction | N monuments | N settlements | Notes |
|---|---|---|---|---|---|---|---|
| Pleiades (pre-2000 BCE, Egypt/Levant) | **6.74** | **−2.91** | **9.65** | LC | 406 | ~372 | Primary finding; 0/10,000 random circles match; max random D = 0.90 |
| Megalithic Portal (non-EU subset) | 9.06 | 4.50 | **4.56** | LC | 4,798 | 1,134 | Land-constrained; 1,000 trials; σ=2°, threshold=50 km |
| DARE | −3.28 | −9.22 | **5.94** | LC | 921 | 18,255 | Land-constrained; monuments more enriched relative to baseline than settlements; consistent with overall negative enrichment on DARE |
| p3k14c | — | — | 6.18 | UC | — | — | Keyword-classified LC attempt produced degenerate result (0 classified sites within 50 km); uncorrected D from predecessor analysis retained for reference |
| XRONOS (full global) | — | — | 3.28 | UC | 2,334 | 5,981 | Keyword-classified |
| XRONOS (without Egypt) | — | — | −14.05 | UC | — | — | Egypt-dependence confirmed |
| XRONOS (without Africa) | — | — | −13.98 | UC | — | — | Consistent with Memphis concentration |
| Peru Ministry of Culture | — | — | −1.25 | UC | 1,733 | 275 | South America null |
| Pleiades/p3k14c South America | — | — | 0.13 | UC | — | — | South America null |
| Historic England | — | — | 0.0 | — | — | — | Negative control; confirmed null |
| Megalithic Portal (New World subset) | — | — | 0.13 | UC | — | — | Old World vs. New World decomposition |

**Habitability-adjusted divergence:** D = 10.44, 99.63rd percentile of 4,323 matched circles.
**Divergence uniqueness:** 0/10,000 random great circles produce D > 9.65 (LC) on Pleiades ancient sites. Random distribution: mean = 0.00, std = 0.07, max = 0.90.
**Geographic dependence (Pleiades spatial block):** Egypt removed → overall enrichment Z = 4.51, divergence D = −0.64. Overall enrichment survives; Pleiades divergence does not.

---

## Table 3. Temporal decomposition at 500-year resolution

D = Z_monument − Z_settlement for each temporal bin. Pleiades temporal assignment by minDate; p3k14c by calibrated calendar year. Bins with monument counts below min_sites_per_bin threshold are flagged. Values from Pleiades and p3k14c independently.

| Period (BCE/CE) | Pleiades D | p3k14c D | Note |
|---|---|---|---|
| 5000–4500 BCE | ~0 | ~0 | Sparse data |
| 4500–4000 BCE | ~0 | ~0 | Sparse data |
| 4000–3500 BCE | ~0 | ~0 | Sparse data |
| 3500–3000 BCE | — | ~0.5 | Pleiades monument count below threshold in some runs |
| 3000–2500 BCE | **9.95** | **3.83** | Peak divergence; Memphis concentration (see Table 4) |
| 2500–2000 BCE | 4.39 | 1.62 | Declining post-Dynasty IV |
| 2000–1500 BCE | 0.76 | 0.84 | Collapse |
| 1500–1000 BCE | — | — | Below threshold in some runs |
| 1000–500 BCE | ~1.2 | ~0.9 | Low divergence |
| 500 BCE–0 CE | 5.98 | 2.41 | Secondary signal: settlement depletion (settlement Z = −5.36); not monument enrichment |
| 0–500 CE | ~1.5 | ~0.7 | Declining |

**Key finding:** D ≈ 0 before 3000 BCE; sharp spike to D = 9.95 (Pleiades) / D = 3.83 (p3k14c) at 3000–2500 BCE; collapse by 2000 BCE; secondary settlement-depletion signal at 500 BCE–0 CE. Secondary signal driven by settlement Z = −5.36, not monument enrichment.

---

## Table 4. 250-year resolution peak detail

Two-hundred-fifty-year bins centered on the primary divergence peak. Monument Z = LC Pleiades. D = monument Z − settlement Z.

| Period (BCE) | Monument Z (LC) | Settlement Z (LC) | D | N monuments on corridor | Note |
|---|---|---|---|---|---|
| 3500–3250 | ~0.3 | ~0.5 | ~−0.2 | <5 | Sparse |
| 3250–3000 | ~0.1 | ~0.4 | ~−0.3 | <5 | Sparse |
| 3000–2750 | −0.79 | ~0.3 | ~−1.1 | <5 | Pre-construction; near threshold |
| **2750–2500** | **11.26** | — | **peak** | **31 of 38** | **Dynasty IV pyramid construction; Memphis necropolis** |
| 2500–2250 | 6.81 | ~1.8 | ~5.0 | 7 | Late Old Kingdom |
| 2250–2000 | 4.39 | ~0 | ~4.39 | — | Declining |
| 2000–1750 | 0.76 | ~0 | ~0.76 | — | Post-collapse |
| 1750–1500 | ~0.5 | ~0.3 | ~0.2 | — | Null |

**Key finding:** Monument Z jumps from −0.79 to 11.26 in one 250-year step. All 38 on-corridor monuments in the Pleiades ancient dataset are Egyptian Old Kingdom pyramids: Abu Rawash, Giza, Abusir, Saqqara, Dahshur. Geographic footprint: 0.44° latitude × 0.55° longitude (~70 km). No monuments from the Levant, Iran, Mesopotamia, or elsewhere contribute.

**Multi-scale profile:** Peak enrichment 9.23× at 7 km corridor half-width; peak divergence D = 2.18 at 18 km; half-life 32 km. The 50 km threshold is conservative relative to peak signal scale.

**Nile distributary context:** Great circle crosses Ahramat Branch (Sheisha et al. 2022) at approximately 70° (approaching perpendicular, ~110° supplementary angle) in the Abu Sir region. Circle is geometrically perpendicular to monument N–S orientations.

---

## Table 5. Preservation test — Egypt

Three conditions testing whether unrecorded buried settlements eliminate the monument-settlement divergence. All conditions use Pleiades Egypt data as the monument source. Monument Z and Settlement Z are Pleiades LC values within the Egypt window.

| Condition | Description | Monument Z | Settlement Z | D | MC trials | Note |
|---|---|---|---|---|---|---|
| 1. Baseline | Pleiades Egypt data as-is | — | — | 12.83 | 1,000 | No additions |
| 2a. +58 documented buried sites | EES Delta Survey, Moeller 2016, Kemp 2018, Spencer 2024 | — | — | 13.49 | 200 | Strengthens divergence |
| 2b. +58 documented buried sites | Same sources | — | — | 13.10 | 1,000 | Strengthens divergence |
| 3. +350 MC-placed settlements | Random placement in Nile floodplain (25–31°N, 29–33°E) | — | — | 10.21 ± 0.72 | 1,000 | D > 2 in 100% of 1,000 iterations |

**Key finding:** Adding documented or hypothetical buried settlements does not eliminate divergence. Under Condition 3 (most conservative), D > 2 in all 1,000 Monte Carlo iterations. Mean distance of 58 documented buried settlements from great circle: 243 km (Nile Delta and floodplain, far from desert plateau monuments). Preservation asymmetry between monumental and domestic Egyptian architecture is real (Moeller 2016; Kemp 2018; Spencer 2024) but does not explain the spatial pattern.

---

## Table 6. Regional decomposition of Pleiades divergence

D values by sub-region within the Pleiades pre-2000 BCE dataset. All values uncorrected (UC) unless noted. Eastward gradient documented.

| Region | Database | Monument Z | Settlement Z | D | N monuments | N settlements | Note |
|---|---|---|---|---|---|---|---|
| Egypt / Levant | Pleiades pre-2000 BCE | — | — | 8.68 | — | — | Strongest signal |
| Iran / Mesopotamia | Pleiades pre-2000 BCE | — | — | 2.60 | — | — | Moderate |
| South Asia | Pleiades pre-2000 BCE | — | — | 1.07 | — | — | Inconclusive; too few sites |
| Western Anatolia (LuwianSiteAtlas) | LuwianSiteAtlas 2000–1200 BCE | — | — | — | 0 | 0 | Complete null; 0/483 sites within 500 km |
| South America (Peru Ministry) | Peru MoC | — | — | −1.25 | 1,733 | 275 | Null |
| South America (Pleiades + p3k14c) | Combined | — | — | 0.13 | — | — | Null |

**Key finding:** Signal weakens eastward from Egypt. Complete null in western Anatolia and South America. DARE (Roman heartland avoided overall; disproportionately monumental where intersected): monument Z = −3.28, settlement Z = −9.22, D(LC) = 5.94.

---

## Table 7. Continental decomposition — p3k14c overall enrichment

Overall enrichment (Z-score) by broad geographic region, p3k14c (uncorrected baselines). Enrichment ratio = observed / expected under distribution-matched MC.

| Continent / Region | Sites | On-corridor sites (50 km) | Z (UC) | Enrichment ratio | Note |
|---|---|---|---|---|---|
| Global | 36,693 | 187 | 7.77 (UC) / **2.09 (LC)** | — | LC = primary figure |
| Old World | — | — | 9.45 | — | From hemisphere block test |
| New World | — | — | 13.41 | — | Geographic; both monuments and settlements cluster |
| Near East (NERD regional) | 11,072 dates | — | — | 1.10× | Indistinguishable from regional background; shows early Holocene enrichment is geometric artifact |
| Early Holocene global (10,500–8,500 BP) | — | — | 5.26 (vs. FC-passing circles) | 4.52× | 0/9 Fertile Crescent-passing circles match; reflects multi-regional geometry, not Near East corridor |

---

## Table 8. Alternative hypotheses tested

All tests subject to Benjamini-Hochberg FDR correction at q = 0.05 unless noted as not primary tests. Status: FALSIFIED = fails all sub-tests; NULL = no significant result; MARGINAL = survives some but not BH correction; CONFIRMED = survives correction.

| Hypothesis | Key statistic | Result | Status | Notes |
|---|---|---|---|---|
| Geographic corridor (habitability) | 78.5th percentile, 4,323 matched circles | Not significant | EXPLAINS overall enrichment | Does not explain monument-settlement divergence |
| Preservation asymmetry (Egypt) | D = 10.21 ± 0.72 after +350 buried settlements; D > 2 in 100% of 1,000 iterations | Does not eliminate divergence | NULL | Documented buried settlements are 243 km from circle mean |
| Post-hoc selection (circle defined on data) | 24% of 200 random 15-site fits exceed Alison Z; split-sample mean Z = 9.45, min Z = 7.31 across 100 splits | Split-sample confirms; raw Z inflated | PARTIAL — split-sample validates | Post-hoc fitting inflates Z to 76th percentile; split-sample is primary validation |
| Research intensity / sampling bias | Pre-2001 Z = 0.06; post-2001 Z = 16.74 | Stronger in newer data (opposite of bias) | NULL (bias direction falsified) | Post-2001 signal concentrated in 5 grid cells (76.3% of near-corridor dates) |
| Astronomical encoding (ecliptic, galactic, precession, solstice/equinox) | Rayleigh p = 0.75; all sub-tests null after Bonferroni | No significant alignment | FALSIFIED | Monuments oriented N–S, perpendicular to E–W circle bearing at Memphis |
| Giza longitude grid (Hancock 1998) | Giza at 71st–79th percentile; precessional Z = 0.45–0.94; settlements more aligned than monuments (15.4% vs. 9.3%); youngest sites most aligned | All 6 sub-tests fail | FALSIFIED | Same method that confirms great circle falsifies grid |
| Pre-Younger Dryas civilization (Hancock 1995) | 10 dates vs. expected 127, Z = −0.47, 33rd percentile; on-corridor YD continuity 2.1% vs. 5.6% off-corridor | Corridor empty pre-YD | FALSIFIED | 5 candidate Atlantis regions all null |
| 108° angular grid | Z = −8.01 | Strongly negative | FALSIFIED | — |
| Trade network (radiocarbon density) | r = 0.04, p = 0.93 (Pleiades temporal correlation) | No correlation | NULL | — |
| Desert edge (Egypt) | Mean distance from western escarpment: 134 km; circle crosses Nile Valley diagonally | Circle does not track desert edge | NULL | — |
| Geophysical (geoid, ocean fraction, currents) | Geoid Z = +2.33; ocean fraction Z = −2.05; current alignment Z = +2.08; 0/13 tests survive Bonferroni | Consistent with continental margins; no specific cause | NULL | — |
| Agricultural suitability | p = 0.57 | Not significant | NULL | See Allan 2026b |
| Lithology | p = 0.74 | Not significant | NULL | See Allan 2026b |
| Groundwater | p = 0.37 | Not significant | NULL | See Allan 2026b |
| LGM connectivity | LGM-connected and ocean-separated clusters show comparable enrichment | No difference | NULL | Specific Z-scores not verified against output files; qualitative finding retained as provisional |
| Paleoclimate (Soreq Cave δ18O) | r = −0.77, p = 0.025; BH-adjusted p ≈ 0.075 (n = 8 bins) | Suggestive but does not survive BH correction | MARGINAL | Climate may explain timing, not geometry |
| Monument orientation clustering | Rayleigh p = 0.75 | No directional clustering | NULL | — |

---

## Table 9. Circle family and alternative great circles

Circles defined by anchor triplets passing within 200 km of all three sites. Approximately 66 of 100,000 random great circles pass within 200 km of Giza, Nazca, and Easter Island (~1 in 1,515). Evaluated on 1° pole grid across 60 qualifying circles. D = monument-settlement divergence by triplet comparison (uncorrected).

| Triplet / Circle | Pole (approx.) | D | Rank (by D) | Percentile | Note |
|---|---|---|---|---|---|
| Alison circle (Giza + Nazca + Easter Island) | 59.682°N, 138.646°W | 8.64 | 5th of 60 | 93.2nd | Primary circle; D from triplet comparison; LC D = 9.65 on full Pleiades analysis |
| Top-ranked pole 1 | Within ~2° of Alison | Higher | 1st | — | Details in supplementary |
| Top-ranked pole 2 | Within ~2° of Alison | Higher | 2nd | — | Details in supplementary |
| Top-ranked pole 3 | Within ~2° of Alison | Higher | 3rd | — | Details in supplementary |
| Top-ranked pole 4 | Within ~2° of Alison | Higher | 4th | — | Details in supplementary |
| Athens + Persepolis + Varanasi | — | 6.08 | — | — | Best alternative anchor triplet |

**Key finding:** Alison circle is not the only geometric solution; it ranks 5th by archaeological signal strength (top 7% of qualifying circles). Top 4 poles are all within ~2° of the Alison pole, forming a tight cluster of high-D solutions.

---

## Table 10. Sensitivity analyses

Key sensitivity tests and their interpretation for core findings.

| Test | Result | Interpretation |
|---|---|---|
| Portal split-sample (100 random 50/50 splits, Z on held-out half) | Mean Z = 9.45; minimum Z = 7.31; 100% exceed Z = 5 | Signal is not artifact of training data; survives blinded validation on every split |
| Post-hoc fitting control (200 random 15-site circle fits, full-DB Z) | 24% exceed Alison Z-score (76th percentile) | Post-hoc fitting inflates Z; raw Z-score is not the appropriate significance measure; split-sample is |
| Egypt removed (spatial block) | Overall enrichment Z = 4.51; Pleiades divergence D = −0.64 | Overall enrichment survives Egypt removal; Pleiades divergence is Egypt-dependent (consistent with Memphis concentration) |
| Corridor threshold 25 km | Lower N on corridor; signal qualitatively similar | Conservative threshold; peak signal is at 7 km (9.23× enrichment) |
| Corridor threshold 100 km | Higher N on corridor; divergence D decreases | Noise from off-signal region; 50 km threshold is appropriate compromise |
| Corridor threshold 200 km | Lower divergence; habitability pattern dominant | Too broad; includes background |
| KDE baseline (Scott's rule, bandwidth 1°–5°) | Z range 5.3–7.1 (Portal), 4.1–6.3 (Pleiades) | Signal robust to bandwidth choice |
| Top-5 grid cells removed (post-2001 p3k14c) | Z = −3.68 | Post-2001 signal concentrated geographically; does not predict new regions |
| Top-10 grid cells removed (post-2001 p3k14c) | Signal eliminated | Same interpretation |
| Preservation test: +350 buried settlements (1,000 iterations) | D = 10.21 ± 0.72; D > 2 in 100% of iterations | Divergence survives maximum plausible unrecorded settlement scenario |
| Hemisphere blocks | Old World Z = 9.45; New World Z = 13.41 | Enrichment exists in both hemispheres independently; New World is geographic/habitability; Old World includes divergence |
| Pleiades temporal decomposition: remove 2750–2500 BCE bin | D drops sharply | Divergence is driven by a single 250-year construction window |
| XRONOS without Egypt | D = −14.05 | Egypt-dependence independently confirmed on separate database |
| Paleoclimate lag sensitivity | r = −0.77 at 500-year lag (p = 0.025); BH-adjusted p ≈ 0.075 | Only 8 temporal bins; small sample size makes result fragile; does not survive BH correction |
| Monument orientation (Megalithic Portal) | Rayleigh p = 0.75 | No directional clustering; Old Kingdom pyramids oriented N–S, perpendicular to circle |
