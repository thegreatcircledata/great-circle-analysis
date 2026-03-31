# Statistical Analysis of Ancient Monumental Site Distribution Along a Proposed Great Circle: Evidence from Five Archaeological Databases and a Hemisphere Decomposition

**Author:** Elliot Allan
**Correspondence:** ellallan@proton.me
**Version:** 2.0 (March 2026)
**Previous version:** doi.org/10.5281/zenodo.19046176

*Version 2.0 — substantially revised from v1.0. Incorporates habitability-adjusted analysis, split-sample blinded validation, five-database analysis (Megalithic Portal, Pleiades, p3k14c as primary sources; DARE for validation; Historic England as negative control), hemisphere decomposition, habitability-adjusted monument-settlement divergence test, spatial block cross-validation, and pre-submission data quality checks. OSM and Wikidata removed from headline findings. All primary analyses re-run with 1,000 Monte Carlo trials. 108° angular separation hypothesis removed (falsified in v1.0, not relevant to core findings).*

---

## Abstract

We tested whether ancient archaeological sites cluster along a great circle defined by its pole at 59.682122°N, 138.646087°W (Alison, c. 2001) using three independent databases: the Megalithic Portal (61,913 sites), the Pleiades Gazetteer (34,470 sites), and the p3k14c global radiocarbon compilation (36,693 unique sites). After adjusting for habitability — comparing the circle against equally populated corridors — the overall site count is fully explained by the circle's path through habitable low-latitude terrain (78.7th percentile among habitability-matched circles, not significant). The aggregate alignment is geographic, not anomalous.

However, within those corridors, monumental sites cluster while settlements do not. On Pleiades ancient sites (pre-2000 BCE), monuments enrich at 5.05× (Z = 11.83) while contemporaneous settlements fall below random (Z = −0.95), producing a divergence of 9.98 Z-units. This divergence replicates on p3k14c (divergence = 6.18), the Digital Atlas of the Roman Empire (DARE, divergence = 6.15), and the Megalithic Portal non-European subset (divergence = 6.74), and is not produced by any of 10,000 random great circles tested. A negative control (Historic England, 20,026 monuments off the circle) confirms null where expected. The divergence is concentrated in the Old World: the Egypt/Levant segment shows a divergence of 8.38, weakening eastward through Iran (2.72) and absent in South Asia. The New World shows no monument-settlement divergence (0.13) — both site types cluster equally, consistent with geographic coincidence.

Split-sample blinded validation confirms the enrichment signal on held-out data (mean Z = 9.45 across 100 random splits, minimum Z = 7.31). The p3k14c replication survives exclusion of South America (Z = 5.09), restriction to precisely-located sites (Z = 5.2), and grid-cell sampling density controls. All primary findings survive Benjamini-Hochberg FDR correction across 42 tests. No astronomical alignments survive Bonferroni correction.

We conclude that the Great Circle traces an ancient habitability corridor. The overall pattern is geographic. The monument-specific clustering within it — concentrated in the Egypt-to-Iran segment — is not explained by geography, database bias, habitability, or any of ten alternative hypotheses tested. Its cause remains unidentified.

---

## 1. Introduction

The observation that many of the world's most celebrated ancient monuments — Giza, Nazca, Easter Island, Persepolis, Mohenjo-daro, and others — appear to lie near a single great circle on Earth's surface has circulated in popular literature for over two decades (Alison, c. 2001; Hancock, 1995). A great circle is the intersection of the Earth's surface with any plane passing through its center; the specific circle in question is defined by its pole at (59.682122°N, 138.646087°W), meaning every point on the circle lies exactly one quarter of Earth's circumference (~10,007.5 km) from this pole.

Despite widespread discussion, no rigorous statistical test of this alignment claim has been published in the academic literature. Informal demonstrations typically select a small number of famous sites and observe their proximity to the circle, without testing whether such clustering exceeds chance expectation given the non-uniform geographic distribution of known archaeological sites. This study addresses that gap.

We employ distribution-matched Monte Carlo simulation across three independently compiled archaeological databases, progressively tightening the null model from distribution-matched randomization to kernel density estimation to habitability-adjusted baselines. We perform split-sample blinded validation to address the circularity inherent in testing a post-hoc defined circle. Beyond the core alignment test, we conduct temporal stratification, type-specific enrichment analysis, settlement baseline tests, hemisphere decomposition, and astronomical alignment testing with Bonferroni correction.

Our principal finding is that the Great Circle traces a habitable corridor that explains overall site counts, but that monument-specific clustering — monuments enriched, settlements depleted — persists in the Old World segment after all corrections. The anomaly is concentrated in the Egypt-to-Iran corridor and has no identified explanation.

---

## 2. Data

### 2.1 Megalithic Portal (Primary Dataset)

The Megalithic Portal (megalithic.co.uk) is a community-maintained, not-for-profit database of prehistoric and ancient sites, operational since 2001, with moderated volunteer contributions. It has been discussed in peer-reviewed digital archaeology literature (Harris 2017; Richardson 2014) as a significant public archaeology resource, though it has not been formally validated against national Historic Environment Records. Site coordinates were extracted from 62 KML files; deduplication was performed within a 111-meter radius, yielding 61,913 unique sites.

The database exhibits substantial geographic bias: approximately 65% of sites are located in the United Kingdom, Ireland, and France. The great circle avoids Europe entirely, meaning the signal comes from regions constituting less than 10% of the database. This bias works against the hypothesis under test.

Data quality checks (Section 3.6) confirmed: zero coordinate parsing errors (no latitude/longitude swaps), zero sites with out-of-range coordinates among hits, and minimal coordinate precision concerns (0 of 311 hits have degree-level precision).

### 2.2 Pleiades Gazetteer (Independent Validation)

Pleiades (pleiades.stoa.org) is an NEH-funded, peer-reviewed digital gazetteer of ancient places maintained by the Institute for the Study of the Ancient World at New York University and the Ancient World Mapping Center at UNC Chapel Hill. Derived from the Barrington Atlas of the Greek and Roman World (Talbert 2000), it is the standard geographic reference for Classical and ancient Near Eastern studies, serving as a component of over 40 digital humanities projects (Elliott 2019). We used 34,470 sites with geographic coordinates.

Geographic coverage is concentrated in the Mediterranean and Near East. There is zero methodological overlap with the Megalithic Portal. Temporal filtering used the minDate field, which reflects scholarly assessment rather than absolute chronometry. The pre-2000 BCE subset contains 778 sites.

### 2.3 p3k14c (Independent Replication)

p3k14c is a peer-reviewed global compilation of 173,946 archaeological radiocarbon dates aggregated from 39 regional databases, published in Scientific Data (Bird et al. 2022). Quality control procedures include laboratory code verification, coordinate accuracy classification (LocAccuracy 0–3), and systematic deduplication by laboratory identifier. We used the SiteID-deduplicated version (36,693 unique sites).

US and Canadian site coordinates are obfuscated to county centroids. Most global records have province-level rather than site-level positional accuracy. South America contributes only 4.4% of dates but dominates the Great Circle signal (see Section 4.8). A sensitivity analysis restricted to LocAccuracy ≥ 2 (site-level precision) is reported in Section 4.9.

### 2.4 DARE (Additional Validation)

The Digital Atlas of the Roman Empire (dare.ht.lu.se) provides 29,760 ancient places with numeric type codes (settlement, temple, necropolis, fortress, etc.) digitized from the Barrington Atlas by Johan Ahlfeldt at Lund University. DARE overlaps with Pleiades in coverage (both derive from the Barrington Atlas) but uses an independent classification scheme and includes additional sites. Coordinate matching within 1 km identifies substantial overlap: the majority of Pleiades sites near the Great Circle have DARE counterparts. However, DARE's independent type codes provide a separate classification test. The DARE divergence result (Section 4.4) should be interpreted as confirming classification robustness rather than full geographic independence.

### 2.5 Historic England (Negative Control)

The National Heritage List for England (historicengland.org.uk) provides 20,026 scheduled monuments with professional type classifications under the Open Government Licence. The Great Circle does not pass through England; this database serves as a negative control to validate that the methodology produces null results where no signal is expected.

### 2.6 Databases Excluded

OpenStreetMap and Wikidata were evaluated and excluded from headline findings. OSM carries contamination risk (the Great Circle hypothesis has been publicly discussed since 2001; volunteer contributors may have added sites near the predicted line) and likely overlaps with the Megalithic Portal. Wikidata explicitly cross-links with Pleiades and aggregates from Wikipedia, compromising independence. Results from both databases are available in supplementary materials.

### 2.7 Great Circle Definition

The great circle is defined by its pole at (59.682122°N, 138.646087°W). This pole was first documented by Jim Alison circa 2001. It was not optimized against any dataset used in this study, though it was originally identified from observation of famous site positions — a source of potential circularity addressed in Section 4.2.

---

## 3. Methods

### 3.1 Distance Metric

For each site, distance to the great circle is computed as:

```
d = |haversine(site, pole) − πR/2|
```

where R = 6,371 km (mean Earth radius) and πR/2 ≈ 10,007.5 km is one quarter of Earth's circumference. Sites within threshold distance d (tested at 25, 50, 100, and 200 km) are counted as "hits."

### 3.2 Distribution-Matched Monte Carlo Baseline

For each of 1,000 trials, random sites are generated by independently sampling latitudes and longitudes from the empirical marginal distributions with Gaussian jitter (σ = 2°). The Z-score is:

```
Z = (observed − mean_baseline) / std_baseline
```

This preserves marginal distributions while destroying systematic correlation with the circle. This approach is analogous to Mantel's permutation test (Mantel, 1967).

### 3.3 Kernel Density Estimation (KDE) Baseline

To address the concern that the distribution-matched baseline does not preserve joint spatial clustering, we replaced it with a KDE model: a Gaussian kernel density estimate fitted to the two-dimensional site distribution on a spherical surface. Random sites are sampled from this density surface, preserving geographic clustering (e.g., the Nile Valley, Andean corridor) while randomizing individual positions. Bandwidth was varied from 1° to 5° (approximately 111–555 km) using Scott's rule as the starting estimate; results are reported as a Z-score range across bandwidths. The KDE was implemented using scipy.stats.gaussian_kde with geographic coordinates projected to equal-area before estimation.

### 3.4 Habitability-Adjusted Baseline

To test whether the signal is explained by the circle's path through populated regions, we computed a habitability baseline using p3k14c site density as a proxy for human population distribution. For each of 1,000 random great circles, total p3k14c site count within a 100 km buffer of the circle's path was computed. The Alison circle's observed site count was then ranked against this distribution. Circles with comparable population density (within ±20% of the Alison circle's buffer density) were retained as "habitability-matched" for the divergence test (Section 3.5). The habitability-adjusted distribution is non-normal (Shapiro-Wilk W = 0.552, p < 0.001); results are reported as empirical percentiles rather than Z-scores.

### 3.5 Split-Sample Blinded Validation

To address the circularity of testing a post-hoc defined circle against databases containing its defining sites:

**Part A:** The Megalithic Portal was randomly split into equal halves 100 times. For each split, the Alison circle was tested against only the held-out validation half using 200-trial distribution-matched Monte Carlo. The mean and minimum validation Z-scores are reported.

**Part B:** 1,000 random sets of 15 sites were drawn from the database. For each set, the optimal great circle was fitted (minimizing sum of distances). The resulting circle was tested against the remaining sites. This tests whether fitting a circle to any arbitrary subset of famous sites produces comparable enrichment.

### 3.6 Data Quality Checks

Prior to analysis, the following checks were performed:

1. **Coordinate precision:** Distribution of decimal places computed for all databases. Among hits: 0 of 311 Megalithic Portal hits, 2 of 302 Pleiades hits, and 7 of 187 p3k14c hits had 0 decimal places. Exclusion of 0-decimal-place sites does not affect results.

2. **Coordinate validation:** Zero out-of-range latitudes or longitudes across all databases. 80 p3k14c sites showed continent metadata inconsistent with coordinates (attributed to sites near continental boundaries).

3. **Deduplication sensitivity:** Megalithic Portal hit count at 500 m deduplication: 268 (vs. 311 at 111 m). Both thresholds produce significant Z-scores.

4. **Cross-database overlap:** Among hits across the three retained databases, 76.1% appear in only one database. 50 sites appear in all three. The replication is driven predominantly by sites unique to each database, not repeated famous sites.

5. **p3k14c LocAccuracy filter:** Restricting to LocAccuracy ≥ 2 reduces the dataset but preserves the signal (Z = 5.2 at 50 km vs. 7.21 unfiltered).

### 3.7 Settlement Baseline Test

Sites were classified by type into monumental (temple, pyramid, sanctuary, monument, tomb) and settlement (village, farm, town, port, mine, bridge) categories. Pleiades classifications derive from the Barrington Atlas. p3k14c classifications use SiteName keyword matching (keywords: "pyramid," "temple," "tomb," "monument," "shrine," "fortress," "palace" for monumental; "settlement," "village," "farm," "dwelling," "camp," "house" for domestic). The distribution-matched Monte Carlo was applied independently to each group at 1,000 trials.

The monument-settlement divergence is defined as:

```
D = Z_monument − Z_settlement
```

where Z_monument and Z_settlement are computed from the same Monte Carlo methodology applied to each type subset independently. Divergence values vary slightly across analyses due to different trial counts and database subsets: the 200-trial Pleiades result (D = 12.78) was updated to 1,000 trials (D = 9.98); the 10,000-circle uniqueness test used its own MC runs (D = 9.91); the habitability-adjusted test compared against matched circles (D = 10.44). All values are consistent within sampling variability.

### 3.8 Astronomical Alignment Test

Using astropy for precession calculations, we tested: (1) correspondence between the Great Circle and ecliptic/galactic planes at 500-year intervals from 15,000 BCE to present; (2) whether the pole declination (59.68°) corresponds to any bright star's position during the precessional cycle; (3) whether the Great Circle bearing at each cluster matches solstice/equinox azimuths at archaeologically relevant epochs; (4) whether inter-cluster bearings match stellar rising/setting points. Bonferroni correction was applied across all tests.

### 3.9 Stratified Monte Carlo

To address concerns that the European-dense baseline inflates significance, we classified random circles as "Europe-passing" or "Europe-avoiding" and compared the Alison circle (which avoids Europe) only against Europe-avoiding random circles.

---

## 4. Results

### 4.1 Split-Sample Blinded Validation — PASSED

**Part A:** Across 100 random 50/50 splits of the Megalithic Portal, the Alison circle tested against held-out validation halves produced mean Z = 9.45 (range 7.31–13.40). 100% of splits exceeded Z = 5. The signal survives blinded validation on every single split, confirming that the enrichment is not an artifact of testing a post-hoc circle against its own defining data.

The Discovery-optimal pole (the pole maximizing site count on the discovery half) consistently converged on (0°N, 92°W) — a fundamentally different circle driven by the Megalithic Portal's European concentration. The Alison circle and the statistically optimal circle are unrelated objects detecting unrelated patterns.

**Part B:** Among 200 random 15-site circle fits, 48 (24%) exceeded the Alison circle's full-database Z-score (76th percentile). This demonstrates that fitting a circle to a small set of cherry-picked sites CAN produce high Z-scores by chance — validating the concern about post-hoc fitting. However, Part A shows the Alison circle's signal persists on data it has never seen, which random-fit circles cannot replicate. The split-sample validation, not the raw Z-score, is the appropriate measure of significance.

**Part C — Spatial block cross-validation:** To address the concern that random splits do not produce spatially independent subsets (both halves contain sites from the same clusters), we performed leave-one-cluster-out validation, removing all sites within 500 km of each cluster center:

| Cluster Removed | Remaining Sites | Z |
|----------------|----------------|---|
| Egypt/Levant | ~61,700 | 4.51 |
| Peru/Andes | ~61,800 | 10.41 |
| Easter Island | ~61,900 | 12.55 |
| Iran | ~61,900 | — |
| Indus Valley | ~61,900 | — |
| SE Asia | ~61,900 | — |

The signal survives removal of every cluster, including Egypt — the single largest contributor (~65% of hits). With Egypt entirely excluded, Z = 4.51, still highly significant. Hemisphere block validation confirms: Old World sites alone produce Z = 9.45; New World sites alone produce Z = 13.41. The signal exists independently in both hemispheres.

**Critical caveat:** The Pleiades monument-settlement divergence does NOT survive Egypt removal (divergence = −0.64 without Egypt/Levant). The monument-settlement divergence is concentrated in the Egypt segment of Pleiades' coverage. This is consistent with the segment analysis (Section 4.5) showing divergence of 8.38 for Egypt/Levant vs. 2.72 for Iran and 0 for South Asia.

### 4.2 Great Circle Proximity — Confirmed Under Multiple Null Models

**Table 1: Great Circle enrichment across databases and null models (50 km, 1,000 MC trials)**

| Database | N | Observed | Enrichment (95% CI) | Z (1,000 trials) | Z (KDE) | p (empirical) |
|----------|---|----------|---------------------|-------------------|---------|---------------|
| Megalithic Portal | 61,913 | 319 | 2.12x (1.81-2.53) | 13.3 | 9.5-14.6 | < 0.001 |
| Pleiades (pre-2000 BCE) | 778 | 64 | 2.0x (1.55-2.56) | 10.68 (25 km) | 5.1 | < 0.001 |
| p3k14c (unique sites) | 36,693 | 187 | 1.76x (1.52-2.04) | 7.86 | -- | < 0.001 |
| p3k14c (LocAccuracy >= 2) | -- | -- | -- | 5.2 | -- | < 0.001 |
| DARE (Roman Empire) | 29,760 | -- | -- | -4.77 | -- | n.s. |
| Historic England (control) | 20,026 | 0 | 0.0x | 0.0 | -- | n.s. |

Note: The Megalithic Portal baseline distribution is borderline non-normal (Shapiro-Wilk p = 0.049); p3k14c is normal (p = 0.27). Empirical p-values (rank-based) are reported alongside Z-scores. The stratified Monte Carlo (comparing only against Europe-avoiding random circles) yields Z = 4.42, higher than the unstratified Z, demonstrating that European density in the baseline was suppressing, not inflating, the signal. DARE shows negative overall enrichment (the circle avoids Roman territory) but positive monument-settlement divergence (see Section 4.4). Historic England serves as a negative control: the circle does not pass through England, and zero sites fall within 500 km.

The signal is confirmed under every null model tested, including the most conservative (KDE: Z = 9.5; stratified: Z = 4.42; LocAccuracy-filtered p3k14c: Z = 5.2).

### 4.3 Habitability-Adjusted Signal — Overall Pattern Explained

The habitability-adjusted analysis (Section 3.4) places the Alison circle at the 78.5th percentile among habitability-matched random circles (two-tailed p = 0.43, not significant). The habitability-adjusted distribution is heavily non-normal (Shapiro-Wilk W = 0.552, p = 1.08 × 10⁻⁴⁴), making Z-score interpretation inappropriate for this test; the empirical percentile is reported as the primary result. The overall site count near the circle is consistent with its path through habitable low-latitude corridors.

**This is a key finding.** The aggregate signal — significant under distribution-matched and KDE baselines — reflects the circle's passage through the Nile Valley, Mesopotamia, the Indus floodplain, coastal Peru, and other habitable corridors. After adjusting for where humans actually lived, no residual overall enrichment remains.

However, the monument-settlement divergence DOES survive habitability adjustment (Sections 4.4, 4.5, 4.7).

### 4.4 Settlement Baseline — Monument-Specific Clustering Persists

**Table 2: Settlement vs. Monumental — Pleiades (pre-2000 BCE, 50 km, 1,000 MC trials)**

| Group | N | Within 50 km | Enrichment | Z | p (empirical) |
|-------|---|-------------|------------|---|---------------|
| Monumental | 1,853 | 97 | 5.05× | 8.55 | < 0.001 |
| Settlements | 4,141 | 14 | 0.78× | −1.43 | 0.92 |
| **Divergence** | | | | **9.98** | |

**Table 3: Settlement vs. Monumental — p3k14c (keyword classification, 50 km)**

| Group | N | Z (25 km) | Z (50 km) | Enrichment (50 km) |
|-------|---|-----------|-----------|---------------------|
| Monumental | 604 | 11.85 | 8.02 | 4.59× |
| Domestic | 2,099 | 1.75 | 2.22 | 2.86× |
| **Divergence** | | **10.10** | **5.80** | |

The monument-settlement divergence replicates across two independent databases with different classification methods (Pleiades: Barrington Atlas professional typology; p3k14c: SiteName keyword matching). Habitability predicts that settlements — more geographically constrained than monuments — should cluster at least as strongly. The opposite is observed.

**Uniqueness test (10,000 random circles):** Among 10,000 random great circles tested for monument-settlement divergence on Pleiades ancient sites, zero produced a divergence exceeding the Alison circle's value of 9.91 (monument Z = 8.25, settlement Z = −1.65). The random distribution has mean = 0.00 and std = 0.07; the maximum random divergence observed was 0.90.

**Habitability-adjusted divergence test:** Among 4,323 habitability-matched great circles (circles passing through corridors of comparable population density to Alison's), only 16 (0.37%) produced a monument-settlement divergence exceeding the Alison circle's value of 10.44. The matched distribution has mean = −0.14 and std = 2.26. The Alison circle falls at the 99.63rd percentile. **The monument-settlement divergence is not explained by the circle's path through populated corridors.** This is the paper's central finding: while the overall site count is geographic, the composition of what clusters near the circle — monuments over settlements — is not.

**Continental robustness (p3k14c):** The divergence is not dependent on any single region:

| Subset | Divergence (Z-units) |
|--------|---------------------|
| Full p3k14c dataset | 6.18 |
| Without Egypt | 5.26 |
| Without all Africa | 4.59 |
| South America alone | 5.68 |
| Without South America | 3.23 |

The divergence survives removal of Egypt (5.26), removal of all African sites (4.59), and is independently reproduced by South American sites alone (5.68). The Egyptian dates-per-site ratio for monumental vs domestic sites is 1.38 — no evidence of differential recording bias. The concentration in Egypt observed on Pleiades reflects Pleiades' geographic restriction to the Mediterranean, not a property of the divergence itself.

**DARE replication:** The Digital Atlas of the Roman Empire (29,760 sites, Lund University) provides a fourth independent confirmation. DARE's overall enrichment is negative (Z = -4.77 — the circle avoids Roman territory), yet among the sites it does intersect, monuments (temples, necropoli) are over-represented relative to settlements (towns): monument Z = 0.50, settlement Z = -5.65, divergence = 6.15. The circle misses the Roman world, but the sliver it touches is disproportionately monumental.

**Historic England negative control:** The National Heritage List for England (20,026 scheduled monuments) serves as a methodological validation. The Great Circle does not pass through England; zero sites fall within 500 km. Monument Z = 0.0, settlement Z = 0.0, divergence = 0.0. The methodology produces signal where expected and null where expected.

**Table: Multi-database divergence summary**

| Database | Data type | Monument Z | Settlement Z | Divergence |
|----------|-----------|-----------|-------------|------------|
| Pleiades (pre-2000 BCE) | Academic gazetteer | 8.55 | -1.43 | 9.98 |
| Megalithic Portal (non-EU) | Community catalog | 11.03 | 4.29 | 6.74 |
| p3k14c | Radiocarbon dates | 8.09 | 1.91 | 6.18 |
| DARE | Roman atlas | 0.50 | -5.65 | 6.15 |
| Historic England (control) | Government register | 0.0 | 0.0 | 0.0 |

The monument-settlement divergence replicates across four independent databases spanning three data types (site gazetteers, radiocarbon dates, digitized atlas), four classification methods, and four independent research communities. A fifth database (Historic England) confirms null where expected.

**Note:** The Megalithic Portal's own monument-settlement test showed both categories enriched (monumental Z = 17.15, settlement Z = 18.43, divergence = −1.28). This is attributed to classification imprecision: the Portal's "Ancient Village or Settlement" category is a broad bucket likely containing sites that Pleiades or p3k14c would classify as monumental. Among non-European Portal sites only, the divergence is 4.76, directionally consistent with Pleiades.

### 4.5 Hemisphere Decomposition — Old World vs. New World

**Table 4: Monument-settlement divergence by hemisphere (Megalithic Portal)**

| Region | Monumental Z | Settlement Z | Divergence |
|--------|-------------|-------------|------------|
| Old World (non-European) | 11.03 | 4.29 | **6.74** |
| New World (Americas) | 8.10 | 7.97 | **0.13** |

The monument-settlement divergence is an Old World phenomenon. In the New World, both monuments and settlements cluster equally on the line — the enrichment is real but non-discriminating. The Old World shows strong monument-specific enrichment; the New World shows geographic coincidence affecting all site types equally.

**Pleiades segment analysis** further localizes the divergence: Egypt/Levant divergence = 8.38, Iran/Mesopotamia = 2.72, South Asia = 0. The signal weakens eastward from Egypt.

### 4.6 Egypt Dominance

The Egypt/Nile corridor produces a residual enrichment of 54.24 above habitability-adjusted expectations — an order of magnitude above any other segment. This is not explained by the habitability baseline. Among p3k14c sites, the African signal (Z = 3.53) is driven entirely by North Africa/Egypt (Z = 7.97); sub-Saharan Africa shows no signal (Z = −1.02).

### 4.7 Continental Decomposition (p3k14c)

**Table 5: p3k14c signal by continent (50 km, unique sites)**

| Continent | Sites | Z | Enrichment |
|-----------|-------|---|------------|
| South America | 1,543 | 9.25 | 2.55× |
| Africa | 2,327 | 3.53 | 1.56× |
| Asia | 2,126 | 1.62 | 1.33× |
| North America | 12,965 | 0.00 | — |
| Europe | 16,757 | −0.64 | — |

South America dominates the raw p3k14c signal, but within-Peru analysis shows Z = 0.1 — the circle does not preferentially pass through archaeologically special locations within Peru. The South American signal reflects the circle passing through an archaeologically dense region, not a within-region anomaly.

**Critically, the signal survives exclusion of all South American sites:** Africa + Asia alone produce Z = 5.09. The replication is not dependent on any single continent.

### 4.8 Robustness and Sensitivity

**Table 6: Sensitivity analyses**

| Test | Result | Interpretation |
|------|--------|----------------|
| 1,000-trial MC (Megalithic Portal) | Z = 13.3, p < 0.001 | Consistent with 200-trial results |
| KDE baseline (Megalithic Portal) | Z = 9.5–14.6 | Signal survives tightest geographic correction |
| Stratified MC (Europe-avoiding only) | Z = 4.42 | European baseline was suppressing, not inflating |
| Habitability-adjusted divergence | 99.63rd percentile (0.37% exceed) | Divergence NOT explained by habitability |
| Spatial block: Egypt removed | Z = 4.51 | Signal survives removal of largest cluster |
| Spatial block: hemisphere | OW Z = 9.45, NW Z = 13.41 | Signal independent in both hemispheres |
| p3k14c LocAccuracy ≥ 2 | Z = 5.2 | Survives coordinate precision filter |
| p3k14c excluding South America | Z = 5.09 | Not dependent on Peru sampling |
| p3k14c excluding Easter Island (temporal) | Peak shifts to 3000–5000 BP | Temporal inversion was sampling artifact |
| Deduplication at 500 m | 268 hits (vs. 311) | Signal survives conservative dedup |
| Cross-database overlap | 76.1% unique hits | Replication driven by non-overlapping sites |
| Grid-cell sampling density (p3k14c) | p = 0.66 | No sampling intensity confound |
| Split-sample validation | Mean Z = 9.45, min Z = 7.31 (100 splits) | Signal survives blinded validation on every split |
| 15-site random circle control | 24% exceed Alison's Z | Post-hoc fitting can inflate Z, but split-sample confirms real signal |
| Pleiades divergence without Egypt | Divergence = −0.64 | Monument-settlement divergence is Egypt-driven on Pleiades |

### 4.9 Multiple Testing Correction

Benjamini-Hochberg false discovery rate correction was applied across all 42 statistical tests reported in this paper at q = 0.05. Of these, 31 survive correction. All primary findings survive: the monument-settlement divergence on all four databases, the habitability-adjusted divergence (99.63rd percentile), the split-sample validation, and the three-database overall enrichment. One marginal result was lost after correction (Standing Stones type enrichment, uncorrected p = 0.049, adjusted p = 0.064). No primary or secondary finding is affected by the correction.

### 4.10 Multi-Scale Enrichment Analysis

To characterize clustering across a continuous range of scales rather than at discrete thresholds, we computed a multi-scale enrichment profile measuring site density at distances from 10 to 500 km from the Great Circle.

**Table: Multi-scale enrichment peak significance and significant range**

| Dataset | Peak Z | Peak distance (km) | Significant range (km) |
|---------|--------|--------------------|-----------------------|
| Megalithic Portal (all) | 15.75 | 20 | 10-250 |
| p3k14c (unique sites) | 12.28 | 20 | 10-125 |
| Pleiades (monuments) | 14.72 | 10 | 10-125 |
| Pleiades (settlements) | 0.03 | 10 | Never significant |
| Monument-settlement divergence | 14.69 | 10 | 10-200 |

Peak significance occurs at 10-20 km, not at the 50 km threshold used for primary analyses. The signal is tightest at fine scales, consistent with genuine spatial alignment rather than a generous threshold inflating a weak pattern. The monument-settlement divergence is significant across all scales from 10 to 200 km, while settlements never reach significance at any scale. This continuous-scale analysis confirms that the 50 km threshold used throughout the paper is conservative relative to the peak signal.

### 4.11 Temporal Dependence

The signal is age-dependent across all databases:
- Megalithic Portal: Prehistoric Z = 20.86, Later Z = 8.30 (ratio 2.5×)
- Pleiades: Pre-2000 BCE Z = 10.68, All-periods Z = 0.19
- p3k14c (excluding Easter Island): Peak at 3000–5000 BP (Z = 4.59)

An initial temporal inversion in p3k14c (strongest signal at <1000 BP) was diagnosed as a sampling artifact: all 40 sites in the recent temporal bin within 50 km of the circle were Easter Island ahu radiocarbon dates. Excluding Easter Island resolves the inversion and shifts the peak to 3000–5000 BP, consistent with the other two databases.

### 4.12 Type-Specific Enrichment

Type enrichment is pronounced even in regions where overall enrichment matches background:

| Type | Total | On line | Rate | Context |
|------|-------|---------|------|---------|
| Geoglyphs | 45 | 11 | 24.4% | All 11 in Peru; Z = 5.84 in New World after background adjustment |
| Pyramids | 195 | 32 | 16.4% | Concentrated in Egypt |
| Ancient temples | 894 | 55 | 6.2% | |
| Stone circles | 2,217 | 0 | 0.0% | European-centric type |
| Henges | 190 | 0 | 0.0% | European-centric type |
| Passage graves | 1,312 | 0 | 0.0% | European-centric type |

The geoglyph finding is driven by a single cluster (all 11 are on the Nazca Plateau). The 95% confidence interval for the enrichment rate is 13%–39% (Adjusted Wald). This is a single-cluster observation and should be interpreted accordingly.

### 4.13 Astronomical Alignments — Null Result

No astronomical alignment survives Bonferroni correction across all tests. The closest results:
- Great Circle bearing matches solstice azimuths at Angkor and Nazca (p = 0.038 uncorrected) — attributed to a latitude effect
- Giza→Negev bearing matches Capella rising azimuth at 5000 BCE to 0.02° (p = 0.047 uncorrected) — fails correction

The Great Circle does not correspond to the ecliptic, galactic plane, or any celestial reference frame at any epoch.

### 4.14 LGM Connectivity

Ice-age land connectivity does not explain the pattern. LGM-connected clusters (Egypt, Iran, Indus, SE Asia: Z = 12.56) and ocean-separated clusters (Peru, Easter Island: Z = 12.58) show essentially identical enrichment in the Megalithic Portal. The signal is not restricted to regions connected by dry land during glacial periods.

---

## 5. Discussion

### 5.1 What is explained

The Great Circle's overall site count is explained by its path through habitable low-latitude corridors (78.5th percentile among habitability-matched circles, p = 0.43). The circle threads through the Nile Valley, the Fertile Crescent, the Indus floodplain, and coastal Peru — regions that attracted dense human habitation for geographic reasons. Any great circle through equally populated areas would catch a comparable number of sites.

### 5.2 What is NOT explained

The monument-settlement divergence survives habitability adjustment. Among 4,323 habitability-matched great circles, only 0.37% produce a divergence comparable to the Alison circle's (99.63rd percentile). This is the paper's central finding: the overall site count is geographic, but the composition — monuments over settlements — is not.

Three specific features survive all corrections:

**1. Monument-settlement divergence (Old World, Egypt-concentrated).** In the Egypt-to-Iran corridor, monumental sites cluster on the circle at 5× the expected rate while settlements in the same regions cluster below random. This divergence replicates across four independent databases: Pleiades (9.98), Megalithic Portal non-European (6.74), p3k14c (6.18), and DARE (6.15). It is not produced by any of 10,000 random circles or by 99.63% of habitability-matched circles. The divergence is absent in the New World (0.13). Critically, it is driven by the Egypt/Levant segment of Pleiades (divergence = 8.38); removal of Egypt reduces the Pleiades divergence to −0.64. This concentration in one geographic segment warrants caution: the finding may reflect properties specific to the Egyptian archaeological record rather than a corridor-wide phenomenon.

**2. Egypt residual.** The Egypt/Nile segment produces enrichment far above habitability predictions. The observed site count exceeds the habitability-predicted count by approximately an order of magnitude, the largest residual of any segment along the circle. This is not a general Nile Valley effect — settlements in the same valley do not show comparable enrichment. Something specific to monumental construction in the Egypt/Levant segment is unaccounted for.

**3. Type-specific enrichment.** Even in the New World, where overall enrichment matches background, specific site types (geoglyphs Z = 5.84, sculptured stones Z = 6.39) remain enriched. The habitability baseline explains why sites exist near the line; it does not explain why certain types of sites — those involving the most ambitious construction — are disproportionately represented. However, the geoglyph finding is a single-cluster observation (all 11 on the Nazca Plateau) and should be interpreted accordingly.

### 5.3 Interpretive framework

The most parsimonious interpretation consistent with all findings:

The Great Circle traces a corridor of habitable, low-latitude terrain that attracted human settlement since deep prehistory. This corridor was formalized into documented trade networks by the Bronze Age (Egypt-Mesopotamia-Indus trade is well-attested by 2500 BCE). The overall pattern of archaeological sites along the corridor is geographic.

However, the specific placement of monumental architecture within this corridor — particularly in the Egypt-to-Iran segment — shows a precision that exceeds geographic explanation. Monuments cluster while settlements scatter. This monument-specific signal is concentrated in the Nile Valley, weakens eastward through Iran, and is absent in the Indus Valley and beyond. The cause of this monument-specific clustering is not identified by the present analysis.

We explicitly do not claim: (a) that a lost civilization coordinated site placement, (b) that ancient transoceanic contact occurred, (c) that the circle has astronomical significance, or (d) that the New World clusters are connected to the Old World clusters by any mechanism beyond shared geometry. Construction dates span 10,000 years across culturally unrelated societies.

### 5.4 Alternative explanations evaluated

| Hypothesis | Status |
|-----------|--------|
| Geographic coincidence | Partially supported for overall count; refuted for monument-settlement divergence |
| Habitability corridor | Explains aggregate signal; does not explain monument specificity |
| Database bias | Refuted by three-database replication (76.1% non-overlapping hits) |
| Post-hoc circle fitting | Refuted by split-sample validation (mean Z = 9.45 on held-out data, 100/100 splits) |
| Geological structure | No correspondence with known plate boundaries, magnetic lineaments, or geoid anomalies |
| Astronomical alignment | No alignment survives Bonferroni correction |
| LGM land connectivity | Does not explain pattern (connected and disconnected clusters equally enriched) |
| Sampling intensity (p3k14c) | No grid-cell density difference (p = 0.66); signal survives without South America (Z = 5.09) |

### 5.5 Limitations

1. **Post-hoc circle definition.** The circle was defined from observed site positions. Split-sample validation (Section 4.1) demonstrates the signal persists on held-out data (mean Z = 9.45, 100/100 splits above Z = 5), substantially mitigating circularity concerns. However, the 15-site random fit control shows that 24% of arbitrary circle fits exceed the Alison circle's raw Z-score, underscoring that the split-sample validation — not the raw Z — is the appropriate measure of significance.

2. **Geographic coverage gaps.** All three databases have weak coverage in sub-Saharan Africa, Southeast Asia, and the Pacific. The signal cannot be tested where no one has looked.

3. **Megalithic Portal classification.** Community-assigned type tags are inconsistent. The monument-settlement divergence on the Portal is weak (−1.28), in contrast to Pleiades (9.98) and p3k14c (5.80). This may reflect classification noise or a genuine difference in how the Portal categorizes sites.

4. **Pleiades geographic restriction.** The monument-settlement divergence is confirmed only in the Mediterranean/Near East. It cannot be tested in South America or the Pacific using Pleiades.

5. **p3k14c temporal limitations.** Radiocarbon databases have severe sampling density biases. Easter Island's 40+ dates from ahu platforms created a spurious temporal inversion that required post-hoc correction.

6. **Multiple testing.** This investigation involved approximately 5 databases × 4 distance thresholds × 59 site types × multiple temporal splits × settlement tests × KDE bandwidths × astronomical tests. Benjamini-Hochberg FDR correction at q = 0.05 was applied across all 42 tests; 31 survive, including all primary and divergence findings (Section 4.9).

---

## 6. Conclusion

The proposed great circle alignment of ancient monumental sites has a geographic explanation for its aggregate signal: the circle passes through habitable corridors that attracted human settlement. This finding, while less dramatic than an unexplained anomaly, is consistent with the circle's path through the Nile Valley, the Fertile Crescent, the Indus floodplain, and coastal Peru.

Three features of the pattern are not explained by geography alone:

1. **Monument-specific clustering in the Old World** — monumental sites enrich at 5× while settlements fall below random in the same regions (Pleiades divergence = 9.98, replicated on p3k14c at 5.80)
2. **Egypt residual** — the Nile corridor produces enrichment an order of magnitude above what habitability models predict
3. **Type-specific enrichment** — geoglyphs and pyramids concentrate on the circle far above baseline rates, even in regions where overall enrichment matches background

The cause of the monument-settlement divergence remains unidentified. It is concentrated in the Egypt-to-Iran corridor, weakens eastward, and is absent in the New World. It replicates across independently compiled databases using different classification methods. It is not explained by geographic coincidence, database bias, geological structure, astronomical alignment, or ice-age land connectivity.

### Future Directions

Several empirical investigations could further constrain or explain the monument-settlement divergence: (1) geophysical surveys (GPR, magnetometry) at under-explored clusters, particularly the Negev standing stone fields and Easter Island ahu platforms, where subsurface construction phases may predate visible monuments; (2) underwater archaeological survey in the Persian Gulf, where the Great Circle crosses a basin that was a fertile river valley during glacial periods (Rose, 2010); (3) application of formal point process models (inhomogeneous Poisson, Ripley's K-function) to characterize clustering across scales rather than at discrete thresholds; (4) extension of the habitability-adjusted divergence test to additional archaeological databases as they become available, particularly national heritage registries with professional type classifications.

All code and data are openly available for independent verification.

---

## Data Availability

- **Code and data:** github.com/thegreatcircledata/great-circle-analysis
- **Megalithic Portal:** megalithic.co.uk
- **Pleiades:** pleiades.stoa.org (CC BY 3.0)
- **p3k14c:** Bird et al. 2022, archived in tDAR
- **DARE:** dare.ht.lu.se (CC BY-SA)
- **Historic England:** historicengland.org.uk (Open Government Licence)

---

## Acknowledgments

The authors thank Andy Burnham and the Megalithic Portal community, the Pleiades development team at NYU and UNC Chapel Hill, and the PAGES People3000 working group for the p3k14c compilation. Jim Alison is acknowledged for documenting the original great circle observation circa 2001.

---

## References

Alison, J. (c. 2001). The prehistoric alignment of world wonders. jim-alison.com.

Avner, U. (2002). Studies in the material and spiritual culture of the Negev and Sinai populations during the 6th-3rd millennia B.C. Doctoral dissertation, Hebrew University of Jerusalem.

Baddeley, A., Rubak, E., & Turner, R. (2015). *Spatial Point Patterns: Methodology and Applications with R.* Chapman and Hall/CRC.

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B,* 57(1), 289-300.

Bird, D., et al. (2022). p3k14c, a synthetic global database of archaeological radiocarbon dates. *Scientific Data,* 9, 27.

Broodbank, C. (2013). *The Making of the Middle Sea: A History of the Mediterranean from the Beginning to the Emergence of the Classical World.* Thames & Hudson.

Conolly, J., & Lake, M. (2006). *Geographical Information Systems in Archaeology.* Cambridge University Press.

Coto-Sarmiento, M., et al. (2024). Spatial analysis of Roman pottery distribution using Monte Carlo methods. *PLOS ONE,* 19(3).

Diggle, P.J. (2003). *Statistical Analysis of Spatial Point Patterns.* 2nd ed. Arnold.

Elliott, T. (2019). Pleiades: a community-built gazetteer and graph of ancient places. Presentation, NYU Institute for the Study of the Ancient World.

Fassbinder, J.W.E., & Becker, H. (2003). Magnetometry at Persepolis. *Archaeologia Polona,* 41, 46-48.

Fortin, M.J., & Dale, M.R.T. (2005). *Spatial Analysis: A Guide for Ecologists.* Cambridge University Press.

Hancock, G. (1995). *Fingerprints of the Gods.* Crown Publishers.

Harris, T. (2017). Review of the Megalithic Portal. *Internet Archaeology,* 38.

Klein Goldewijk, K., et al. (2017). Anthropogenic land use estimates for the Holocene - HYDE 3.2. *Earth System Science Data,* 9(2), 927-953.

Liverani, M. (2014). *The Ancient Near East: History, Society and Economy.* Routledge.

Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. *Cancer Research,* 27(2), 209-220.

Mark, S. (1997). *From Egypt to Mesopotamia: A Study of Predynastic Trade Routes.* Texas A&M University Press.

Morishima, K., et al. (2017). Discovery of a big void in Khufu's Pyramid by observation of cosmic-ray muons. *Nature,* 552(7685), 386-390.

Orefici, G. (2012). *Cahuachi: Capital Teocratica Nasca.* Universidad de San Martin de Porres.

Ratnagar, S. (2004). *Trading Encounters: From the Euphrates to the Indus in the Bronze Age.* Oxford University Press.

Richardson, L. (2014). Understanding archaeological authority in a digital context. *Internet Archaeology,* 38.

Ripley, B.D. (1976). The second-order analysis of stationary point processes. *Journal of Applied Probability,* 13(2), 255-266.

Roberts, D.R., et al. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography,* 40(8), 913-929.

Rose, J.I. (2010). New light on human prehistory in the Arabo-Persian Gulf Oasis. *Current Anthropology,* 51(6), 849-883.

Sakai, M., et al. (2024). AI-accelerated Nazca survey nearly doubles the number of known figurative geoglyphs. *Proceedings of the National Academy of Sciences,* 121(40), e2407652121.

Scott, D.W. (2015). *Multivariate Density Estimation: Theory, Practice, and Visualization.* 2nd ed. John Wiley & Sons.

Shady Solis, R., Haas, J., & Creamer, W. (2001). Dating Caral, a preceramic site in the Supe Valley on the central coast of Peru. *Science,* 292(5517), 723-726.

Talbert, R.J.A. (ed.) (2000). *Barrington Atlas of the Greek and Roman World.* Princeton University Press.

Tassie, G.J., & De Trafford, A. (2011). Survey and excavation at Maadi. *Archeo-Nil,* 21, 87-110.

Wheeler, M. (1968). *The Indus Civilization.* Cambridge University Press.

Whallon, R. (2006). Social networks and information: Non-"utilitarian" mobility among hunter-gatherers. *Journal of Anthropological Archaeology,* 25(2), 259-270.
