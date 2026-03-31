# The Alison Great Circle: Temporal Decomposition, Geographic Resolution, and 41 Tests of the World's Most Famous Archaeological Alignment

**Author:** Elliot Allan
**Correspondence:** ellallan@proton.me
**Preprint (Paper I):** doi.org/10.5281/zenodo.19046176
**Data:** github.com/thegreatcircledata/great-circle-analysis

---

## Abstract

A companion study (Allan 2026; under review, PLOS ONE, PONE-D-26-13610) established that ancient monumental sites cluster along a specific great circle (pole: 59.68°N, 138.65°W) at 5.05× the expected rate while contemporaneous settlements do not, replicated across seven databases. This study conducts 41 follow-up tests to decompose and explain the phenomenon.

The monument-settlement divergence onsets sharply at 2750–2500 BCE, driven entirely by 38 Egyptian Old Kingdom pyramids in the Memphis necropolis. It collapses without Egypt (XRONOS D: +3.28 → −14.05). A preservation test adding 350 estimated buried settlements confirms the divergence is not a recording artifact. Geographic features are sufficient to predict monument vs. settlement classification without Great Circle distance as a predictor (XGBoost, 10 features, 21,779 sites; GC distance SHAP rank #10/10).

The approximate collinearity of Giza, Nazca, and Easter Island is scaffolded by plate tectonics: the geological features enabling these three civilizations define a pole 0.8° from the Alison pole. Easter Island is the only inhabited Pacific island within 200 km of the resulting circle. The sub-km precision is a 1-in-3,000 rarity among equivalent-site combinations.

The corridor was empty before the Younger Dryas (10 of 94,181 radiocarbon dates, 33rd percentile). The proposed Giza longitude grid (Hancock 1998) is not supported (six tests, 508,000 sites). All code and data are open-source.

**Keywords:** spatial statistics, monument-settlement divergence, great circle, Egyptian Old Kingdom, machine learning, plate tectonics

---

## 1. Introduction

Allan (2026, hereafter Paper I) tested a long-standing but previously unexamined claim: that many of the world's ancient monumental sites — Giza, Nazca, Easter Island, Persepolis, Mohenjo-daro — lie near a single great circle on Earth's surface (Alison, c. 2001). Using distribution-matched Monte Carlo simulation across six independent databases totaling over 550,000 entries, Paper I found statistically significant clustering (Z = 25.85 at 50 km on the primary dataset) and, more critically, a monument-settlement divergence: ancient monumental sites cluster at 5.05× the expected rate (Z = 11.83) while contemporaneous settlements do not (Z = −0.95). This divergence was not matched by any of 10,000 random great circles.

Paper I left several questions unresolved. First, was the circle fit to the data? Alison proposed it based on visual inspection of a small number of famous sites; if the pattern depends entirely on those sites, it may reflect post-hoc selection rather than a genuine spatial regularity. Second, the monument-settlement divergence was reported at millennial resolution; finer temporal analysis could reveal whether it onsets gradually (consistent with recording bias) or sharply (consistent with a historical event). Third, Egypt's well-documented preservation asymmetry — stone monuments survive while mudbrick settlements dissolve — was raised as a potential confound but not quantitatively tested. Fourth, the geographic extent of the divergence remained unclear.

This study addresses these questions through five analyses: predictive validation on post-2001 data, temporal decomposition at 250-year resolution, a quantitative preservation test, paleoclimate proxy correlation, and regional decomposition including three new databases. Two additional analyses extend the temporal scope: a pre-Younger Dryas corridor test using 94,181 merged radiocarbon dates to determine whether the corridor shows anomalous activity before 12,800 BP, and a test of the proposed Giza-centered precessional longitude grid (Hancock 1998).

---

## 2. Data

### 2.1 Databases from Paper I

**Pleiades Gazetteer** (pleiades.stoa.org): 34,470 georeferenced ancient sites with minDate fields. Sites classified as monumental (temple, sanctuary, pyramid, monument, amphitheatre, aqueduct; N = 1,853 ancient) or settlement (village, town, settlement, farm, city, port; N = 4,141 ancient). Pleiades records include creation dates enabling a pre/post-2001 split: 20,530 records created 2009–2011 (Barrington Atlas digitization, representing pre-2001 knowledge) and 13,940 records added subsequently.

**p3k14c radiocarbon database** (Bird et al. 2022; building on Martindale et al. 2016): 170,150 radiocarbon dates with geographic coordinates. Reference publication years enable temporal splitting: 45,171 dates from references published before 2001; 102,567 from references published 2001 or later.

### 2.2 New Databases

**XRONOS** (Palmisano et al. 2025): 305,400 radiocarbon records from 28,127 unique sites across 125 countries. Monument-settlement classification by keyword matching yields 2,334 monumental and 5,981 settlement sites.

**Peru Ministry of Culture** (geoportal.cultura.gob.pe; cf. Wernke et al. 2024): 17,465 georeferenced archaeological sites. Expanded keyword classification (Spanish and Quechua terms) yields 1,733 monumental and 275 settlement sites.

**SAAID** (South American Archaeological Isotopic Database, v2.0): 2,429 unique sites across South America with researcher-assigned site type and function fields. Classification yields 209 monumental, 530 settlement, and 180 cemetery sites.

**CNSA/IPHAN** (Cadastro Nacional de Sítios Arqueológicos, Brazil): 27,582 sites with description-field classification yielding 5,260 monumental and 12,420 settlement sites.

**LuwianSiteAtlas** (Scientific Data, 2025): 483 Bronze Age settlement sites in western Anatolia (2000–1200 BCE) with professional classification.

### 2.3 Buried Settlement Sites

58 ancient Egyptian settlement sites (pre-1000 BCE) known to exist but inaccessible to full excavation, compiled from the EES Delta Survey (Spencer 2024), Moeller (2016), Kemp (2018), Bietak (1996), and site-specific publications. 36 Delta sites and 22 Nile Valley sites, geocoded from published coordinates or nearest modern settlement. Full list in Supplementary Table S1.

### 2.4 Pre-Younger Dryas Merged Database

For the pre-YD corridor test, three radiocarbon databases were merged and deduplicated: p3k14c (170,117 dates), XRONOS (278,936 records), and ROAD v32 (18,755 European Palaeolithic dates; Vermeersch 2024). Deduplication by coordinates (within 0.01°) and age (within 200 years) yielded 94,181 unique dates, of which 16,483 predate 12,800 BP. ROAD v32 provides critical European Palaeolithic coverage where p3k14c is thin.

### 2.5 Hancock Longitude Grid Sites

For the longitude grid test, 508,000 sites were compiled from the merged databases used throughout this study. The grid was defined as longitudes at multiples of 36° from Giza (31.13°E), following Hancock's claim (Heaven's Mirror, 1998) that sites at these longitudes encode precessional knowledge.

### 2.6 Paleoclimate Proxies

**Soreq Cave speleothem** (Bar-Matthews et al. 2003): δ18O record from central Israel providing continuous Holocene climate proxy for the Eastern Mediterranean. More negative values indicate wetter conditions.

**GISP2 ice core** (Grootes et al. 1993): δ18O record from central Greenland providing Northern Hemisphere temperature proxy.

---

## 3. Methods

### 3.1 Distribution-Matched Monte Carlo Baseline

Z-scores computed following Paper I methodology (cf. Ripley 1976; Bevan 2012): for each of 1,000 trials, random sites generated by independently selecting latitude and longitude from the empirical distribution with Gaussian jitter (σ = 2°). Z = (observed − μ_baseline) / σ_baseline.

### 3.2 Predictive Validation

Databases split by temporal proxy: Pleiades by record creation date (pre-2012 vs. post-2012); p3k14c by reference publication year (pre-2001 vs. post-2001). Z-scores computed independently for each split at 25 km and 50 km thresholds.

To control for geographic concentration of fieldwork, near-circle dates were analyzed by 5° × 5° grid cell. The signal was tested after removing (a) the five cells with the highest near-circle date counts, and (b) the ten cells with the highest total date counts regardless of circle proximity.

### 3.3 Temporal Divergence Analysis

The divergence metric D = Z_monument − Z_settlement, computed independently for each temporal bin. Two resolutions: 500-year bins (5000 BCE to 500 CE) and 250-year bins (3500 BCE to 500 CE). Pleiades sites binned by minDate; p3k14c dates by calibrated calendar year.

### 3.4 Preservation Test

Three variants: (1) Baseline Pleiades Egypt data; (2) +58 documented buried settlements; (3) +350 settlements randomly placed within the Nile floodplain (25–31°N, 29–33°E), 1,000 Monte Carlo iterations.

### 3.5 Paleoclimate and Trade Network Correlation

Divergence D values at 500-year resolution correlated with contemporaneous and lagged (0, 250, 500 year) Soreq Cave and GISP2 δ18O values. Radiocarbon date density within the Egypt-Iran corridor (20–38°N, 25–60°E; 4,400 p3k14c dates) used as trade/activity proxy.

### 3.6 Deep-Time Radiocarbon Analysis

To test whether the circle captures a pattern of human activity older than the monument-building era, all 170,150 p3k14c radiocarbon dates (Bird et al. 2022) were binned by 500-year windows from 20,000 BP to 5,000 BP, following established radiocarbon-as-demographic-proxy methods (Rick 1987; Surovell et al. 2009; Crema 2022; cf. Palmisano et al. 2017). For each bin, the ratio of dates within 200 km of the circle to the global total was computed and compared to the expected ratio under the null hypothesis. This analysis was repeated for 100 random great circles to determine whether any deep-time enrichment is specific to the Alison circle. A separate test examined the Younger Dryas period (12,800–11,600 BP) specifically, as this interval has been proposed as a catastrophic disruption to a putative Ice Age civilization (Hancock 1995).

### 3.7 Geophysical Property Scan

To test whether the circle corresponds to any measurable physical property of the Earth, four geophysical datasets were sampled at 360 points along the circle (1 per degree of arc) and compared to 100 random great circles: (1) EGM2008 geoid heights (mean height and gradient); (2) EMAG2v3 crustal magnetic anomalies (mean and variance); (3) ETOPO1 global relief (ocean fraction, coastal fraction, elevation variance); (4) ocean surface current alignment (fraction of ocean-crossing segments within 30° of a major current). Z-scores computed against the random circle distribution.

### 3.8 Regional and Site-Level Decomposition

Monument-settlement divergence computed for five regions. All Pleiades monuments within 50 km of the circle in the 2750–2500 BCE peak bin individually identified by name, type, coordinates, and geographic location to determine whether the spike represents distributed corridor activity or geographic clustering.

### 3.9 Pre-Younger Dryas Corridor Test

To test whether the corridor shows anomalous human activity before the Younger Dryas (>12,800 BP), the merged database (Section 2.4; 94,181 dates) was classified as on-corridor (within 50 km of the Great Circle) or off-corridor. Pre-YD date counts were compared against 1,000 random great circles. Seven sub-tests were conducted: (1) temporal enrichment profile in 500-year bins from 25,000 to 5,000 BP; (2) pre-YD enrichment Z-score vs. random and habitability-matched corridors; (3) YD disruption profile comparing the crash ratio (dates in 500 years post-YD onset / 500 years pre-YD) on-corridor vs. off-corridor; (4) regional decomposition of pre-YD dates; (5) continuity test identifying site clusters with dates spanning both sides of the YD boundary; (6) candidate "Atlantis" region analysis testing five proposed locations; (7) simplified summed probability distribution (SPD) comparing on-corridor vs. off-corridor temporal density. Sensitivity checks at 25 km and 100 km thresholds and excluding North America (where p3k14c coordinates are fuzzed to county centroids) were also conducted.

### 3.10 Giza Longitude Grid Test

Hancock (1998) proposed that ancient sites are placed at precessionally significant longitude intervals from Giza — specifically multiples of 36° (with 72°, 108°, and 144° highlighted as precessional numbers). Six tests were conducted against 508,000 sites: (1) site clustering at 36° multiples vs. continental geography control; (2) Giza vs. 1,000 random reference longitudes; (3) precessional angles vs. rotated grids and random angle sets; (4) precision test of claimed site pairs (Angkor, Nazca, Persepolis); (5) monument vs. settlement grid alignment; (6) age-alignment correlation.


### 3.11 Machine Learning Geographic Decomposition

Three classifiers (XGBoost, Random Forest, Logistic Regression) were trained on 21,779 Pleiades sites (4,012 monuments, 17,767 settlements) to predict monument vs. settlement classification using 10 geographic features: distance to the Great Circle (km), elevation (ETOPO1 1-arcmin), coastline distance (Natural Earth 10m), river distance (Natural Earth 10m), absolute latitude, longitude, site density (count within 100 km), seismic hazard (GEM Global Seismic Hazard Map 2023.1), cloud cover (EarthEnv MODCF mean annual), and cross-type distance (distance to nearest opposite-type site). Two candidate features (soil productivity, water table depth) were excluded due to incomplete global coverage (>30% missing values).

Train/test split: 80/20 stratified by site type. XGBoost hyperparameters: max_depth = 6, n_estimators = 200, learning_rate = 0.1, min_child_weight = 5, subsample = 0.8. Five-fold cross-validation was used to confirm stability (AUC SD < 0.01 across folds). SHAP values (Lundberg & Lee 2017) computed on the test set using TreeExplainer. Ablation testing: model retrained with and without the Great Circle distance feature; AUC difference reported.

A secondary corridor-only model was trained on sites within the Egypt-Levant-Iran corridor (20-38°N, 25-60°E; 8,412 sites) using the same features and hyperparameters.

### 3.12 Systematic Great Circle Search

All possible great circle pole positions were scanned on a 2° latitude-longitude grid (16,200 poles). For each pole, the monument-settlement divergence D was computed using the same Monte Carlo methodology as Section 3.1, with 100 trials per pole. Poles with D > 10 were flagged as hot spots; 137 hot spots were refined to 0.5° resolution with 1,000-trial Monte Carlo each.

Independent clusters were identified by grouping poles with D > 5 and angular separation < 20°. For each cluster, the peak pole was identified and the sites captured within 50 km were characterized. An anchor precision test computed the minimum distance from each scanned pole to Giza, Nazca, and Easter Island simultaneously.

### 3.13 Exhaustive Triplet Scan

All C(149,3) = 540,274 triplets from a curated list of 149 globally famous archaeological sites (drawn from UNESCO World Heritage Sites and widely cited monumental complexes) were evaluated. For each triplet, the best-fit great circle was computed via eigendecomposition of the 3×3 Gram matrix of unit position vectors; the eigenvector corresponding to the smallest eigenvalue defines the pole. Geometric fit was measured as the mean great-circle distance from each site to the best-fit circle. Monument-settlement divergence (Poisson D) was computed for each triplet's circle using the Pleiades database. A combined metric D/(fit + 1) was used to rank triplets by both geometric precision and archaeological signal strength.

### 3.14 Terrain and Visibility Tests

Five mechanistic tests were conducted to identify specific landscape properties underlying the geographic decomposition:

1. **Terrain transition zone.** Elevation gradients were computed at 1° intervals along the Alison circle using ETOPO1 data and compared to 200 random great circles.
2. **Visibility/prominence.** Topographic Position Index (TPI; elevation minus mean elevation within 10 km radius) was computed for all Pleiades monuments and settlements within 100 km of the circle and within a 500 km offset control corridor.
3. **Navigational waypoints.** Coastal crossing points were identified where the circle intersects continental coastlines (Natural Earth 10m); topographic prominence at each crossing was compared to 100 random coastal points.
4. **Monument elevation differential.** TPI values for near-circle monuments vs. near-circle settlements were compared using a Mann-Whitney U test.
5. **Nile Valley constriction.** Valley width was measured as the distance between the 100 m contours on either side of the Nile at 0.1° latitude intervals from 29.0°N to 30.5°N. The Great Circle's crossing latitude was compared to the latitude of minimum valley width.

### 3.15 Use of Generative AI Tools

All computational analysis for this study was designed, executed, and iterated using Claude Code (Anthropic, Claude Opus 4 model), an AI-assisted coding tool that generates and runs Python scripts from natural language directives. The author wrote research directives specifying each analysis (hypothesis, datasets, statistical methods, expected outputs); Claude Code implemented the corresponding Python code, executed it, and returned results. All code was reviewed by the author for correctness and is openly available (github.com/thegreatcircledata/great-circle-analysis). Perplexity AI was used for literature search and fact-checking of archaeological claims. Intermediate results were reviewed through iterative AI-assisted critique, in which the author prompted Claude (Anthropic) to identify methodological weaknesses, propose additional controls, and flag potential confounds — functioning as a structured adversarial review process. All hypotheses, interpretations, conclusions, and scientific judgments in this paper are the author's own. The AI tools were used as computational instruments, not as intellectual contributors.

---

## 4. Results

### 4.1 Predictive Validation

**Table 1: Z-scores by database split (pre- vs. post-2001)**

| Database | Split | N | Z (50 km) | Z (25 km) | Enrichment (25 km) |
|----------|-------|---|-----------|-----------|-------------------|
| Pleiades | Pre-2001 | 20,530 | −1.69 | +1.12 | 1.08× |
| Pleiades | Post-2001 | 13,940 | +2.55 | +5.30 | 1.72× |
| p3k14c | Pre-2001 | 45,171 | +0.06 | +3.97 | 1.37× |
| p3k14c | Post-2001 | 102,567 | +17.24 | +27.26 | 3.11× |

Both databases show stronger clustering in post-2001 data — the opposite of what post-hoc selection predicts. Pre-2001 Pleiades data (Barrington Atlas sites) shows Z = −1.69 at 50 km, slightly below random. Pre-2001 p3k14c data shows Z = 0.06, indistinguishable from random.

However, fieldwork intensity analysis reveals the post-2001 p3k14c signal is geographically concentrated. Of 583 near-circle dates, 76.3% fall in five grid cells: Easter Island (206 dates, 100% near-circle), Egypt (82), Levant (63), Peru (60), and Iran (34). Removing these five cells (871 of 102,567 dates) eliminates the signal (Z = −3.68). Removing the ten most research-intensive cells globally (37,308 dates, predominantly UK and North America) does not affect it (Z = +17.74).

**Table 2: Continental decomposition of post-2001 p3k14c signal**

| Continent | N | Near-circle | Z | Enrichment |
|-----------|---|-------------|---|------------|
| Easter Island | 701 | 206 | +28.60 | 5.59× |
| Middle East | 2,151 | 235 | +10.53 | 1.88× |
| South America | 3,417 | 105 | +7.38 | 1.98× |
| Europe | 46,502 | 0 | −1.64 | 0.00× |
| North America | 41,565 | 0 | +0.00 | 0.00× |
| Africa | 3,709 | 37 | −2.15 | 0.70× |

The signal is driven by three regions — Easter Island, the Middle East, and coastal Peru — which together constitute 6.1% of post-2001 dates but 93.3% of near-circle dates. Europe and North America (88% of all data) contribute zero near-circle dates.

The predictive validation confirms the alignment is not solely dependent on sites known in 2001, but it does not demonstrate that the circle predicts archaeological activity in previously unexplored regions. The post-2001 signal reflects continued fieldwork in the specific regions the circle traverses.

### 4.2 XRONOS Replication

XRONOS independently replicates the overall alignment: Z = 24.45 at 50 km (28,127 sites). Monument-settlement divergence: D = +3.28 (full dataset), D = −14.05 (without Egypt), D = −13.98 (without Africa). The divergence is Egypt-dependent on this database.

### 4.3 Temporal Divergence at 500-Year Resolution

**Table 3: Monument-settlement divergence by 500-year bin**

| Period | Pleiades Mon Z | Pleiades Set Z | Pleiades D | p3k14c D |
|--------|---------------|---------------|------------|----------|
| 4000–3500 BCE | 0.47 | 0.34 | 0.13 | — |
| 3500–3000 BCE | −0.38 | 0.21 | −0.59 | −0.42 |
| **3000–2500 BCE** | **9.87** | **−0.08** | **9.95** | **3.83** |
| 2500–2000 BCE | 2.44 | −0.43 | 2.87 | 1.42 |
| 2000–1500 BCE | 0.31 | 0.23 | 0.08 | 0.87 |
| 1500–1000 BCE | −0.89 | 0.42 | −1.31 | −2.51 |
| 500 BCE–0 CE | 5.98 | — | 5.98 | 0.08 |

The divergence is absent before 3000 BCE (D ≈ 0), spikes to D = 9.95 at 3000–2500 BCE, and collapses to D = 0.08 by 2000–1500 BCE. Both databases independently confirm this timing (Figure 1).

### 4.4 Temporal Divergence at 250-Year Resolution

Higher temporal resolution (Figure S2) narrows the signal further.

**Table 4: 250-year temporal resolution (Pleiades)**

| Period | Mon Z | D | Note |
|--------|-------|---|------|
| 3000–2750 BCE | −0.79 | −0.56 | No signal |
| **2750–2500 BCE** | **11.26** | **—** | **31/38 monuments on circle** |
| 2500–2250 BCE | −1.20 | −0.73 | Collapse |
| 2250–2000 BCE | 3.63 | 4.39 | Residual |
| 2000–1750 BCE | −0.55 | 0.76 | Gone |

The onset is sharp: monument Z jumps from −0.79 to 11.26 in a single 250-year step. The collapse is equally sharp, dropping from D = 4.39 to 0.76 between 2250–2000 and 2000–1750 BCE.

### 4.5 Site-Level Identification of the 2750–2500 BCE Spike

All 38 Pleiades monuments within 50 km of the circle in the 2750–2500 BCE bin are Egyptian Old Kingdom pyramids and associated structures in the Memphis necropolis: Abu Rawash, Giza, Abusir, Saqqara, and Dahshur. They span 0.44° of latitude and 0.55° of longitude (~70 km). No monuments from the Levant, Iran, Mesopotamia, or any other region contribute to the spike.

The monument-settlement divergence during its peak period is not a corridor phenomenon. It is the Memphis pyramid fields — the densest concentration of monumental construction in the ancient world — sitting within 18–44 km of the great circle.

### 4.6 Preservation Test

**Table 5: Egypt preservation test (1,000 Monte Carlo iterations)**

| Test | Monument Z | Settlement Z | D |
|------|-----------|-------------|---|
| Baseline (Pleiades) | 9.93 | −2.90 | 12.83 |
| +58 buried settlements | 10.41 | −3.08 | 13.49 |
| +350 estimated (mean ± SD) | 9.93 ± 0.41 | −0.33 ± 0.52 | 10.21 ± 0.72 |

Adding 58 documented buried settlements increases D because the buried sites cluster in the Nile Delta and floodplain (mean distance from circle: 243 km), far from the desert plateau edge where the circle passes. Under the Monte Carlo model, D exceeds 2 in 100% of 1,000 iterations (Figure 4).

**Asymmetry note.** This test addresses one direction of the preservation bias: missing settlements. It does not test the complementary possibility of missing monuments — destroyed or unexcavated monumental sites away from the circle. However, destroyed Egyptian monuments (provincial temples, Upper Egyptian structures) are predominantly in the same regions as surviving ones (Memphis, Thebes, Abydos) rather than in new locations, so missing monuments would not substantially alter the geometry.

The preservation asymmetry between Egyptian monumental and domestic architecture is well-documented (Moeller 2016; Kemp 2018; Spencer 2024). It does not explain the spatial divergence. The buried settlements are geographically distant from the great circle; adding them to the analysis reinforces rather than weakens the monument-settlement divergence.

**Multi-region preservation simulation.** The preservation test was extended beyond Egypt to three additional riverine regions with known unexcavated settlement concentrations: Mesopotamia (+200 estimated settlements, Tigris-Euphrates alluvial plain), the Indus Valley (+500 estimated, Punjab/Cholistan), and coastal Peru (+200 estimated, Lima region). Under realistic geographic distribution (Gaussian placement around known settlement centers), D drops from 8.64 to 5.20 ± 0.69, with D > 2 in 100% of 1,000 iterations. The Indus Valley contributes the most near-circle settlements (23.4% within 50 km) because the Indus floodplain runs closer to the circle than the Egyptian floodplain. Peru is negligible (0.2% within 50 km, mean distance 351 km).

Under adversarial placement (all additional settlements placed within the 50 km band), break-even occurs at approximately 214 settlements (0.47× the known on-corridor settlement count). The true preservation bias lies between these bounds: realistic geographic distribution preserves the divergence comfortably, while adversarial placement makes it fragile. The paper reports both.

### 4.6b Publication and Excavation Bias Control

A critical confound: are on-corridor sites simply more excavated? Three databases were tested for research intensity bias:

| Database | Metric | On-corridor | Off-corridor | p-value |
|----------|--------|-------------|--------------|---------|
| p3k14c | Median dates/site | 3.0 | 2.0 | 0.142 (n.s.) |
| XRONOS | Median dates/site | 4.0 | 4.0 | 0.595 (n.s.) |
| Pleiades | Median edits/site | 3 | 4 | <0.001 (reversed) |

On-corridor sites are NOT more researched. Pleiades actually shows fewer edits for on-corridor sites — the bias runs against the signal. Latitude-matched controls show only 2 of 8 latitude bands with significant differences, confined to the 30-32° band (Egypt/Nile Delta). After downweighting heavily-researched areas, enrichment increases from 3.10× to 3.33×. The signal strengthens, not weakens, when controlling for research intensity.

### 4.6c Monument Type Decomposition

The monument category was decomposed into sub-types using Pleiades featureType classifications:

| Sub-type | Enrichment | Examples |
|----------|-----------|---------|
| Mortuary | 4.30× | tomb, pyramid, necropolis, cemetery |
| Ceremonial | 2.06× | temple, sanctuary |
| Administrative | 1.21× | palace, fort |
| Entertainment | 0.22× | theatre, amphitheatre |
| Infrastructure | 0.75× | aqueduct, bridge |
| Settlements | 0.78× | town, village |

The 2750 BCE spike comprises 49 sites, of which 35 (71%) are mortuary: Giza pyramids, Saqqara mastabas, Abusir necropolis — all within 15 km of the Great Circle. Seven are ceremonial (temples). The divergence is primarily a mortuary phenomenon: the circle captures where ancient Egyptians buried their dead, not where they worshipped or lived.

### 4.7 Desert Edge Test

The hypothesis that the circle tracks the Nile Valley desert-floodplain boundary (the escarpment where pyramids were preferentially built) was tested using SRTM elevation data. Of five circle sample points in the Egypt window, only one (longitude ~31.0°E) lies near the western escarpment (4.1 km). The remaining four are 107–226 km from the escarpment. Mean distance: 134 km. The circle crosses the Nile Valley at a steep diagonal; it does not track the desert edge.

### 4.8 Paleoclimate Correlation

**Table 6: Initial correlation between divergence D and paleoclimate proxies (1000-year bins)**

| Proxy | Database | Lag | Pearson r | p-value |
|-------|----------|-----|-----------|---------|
| Soreq Cave δ18O | Pleiades | 0 yr | −0.70 | 0.080 |
| Soreq Cave δ18O | Pleiades | 250 yr | −0.76 | 0.049 |
| Soreq Cave δ18O | p3k14c | 500 yr | −0.77 | 0.025 |
| GISP2 δ18O | Pleiades | 250 yr | +0.70 | 0.083 |

The initial analysis on 7–8 temporal bins suggested a correlation between D and Soreq Cave δ18O at 250–500 year lag (p = 0.025–0.049). However, subsequent high-resolution analysis using actual downloaded NOAA paleoclimate data at 100-year, 250-year, 500-year, and 1000-year resolution across four independent proxies (Soreq Cave 2003, Grant 2012 updated Soreq, LC21 marine core, Jeita Cave Lebanon) found no significant correlation at any resolution or lag.

**Table 6b: High-resolution climate-divergence correlation (4 proxies × 4 resolutions)**

| Resolution | N bins | Best |r| | Best p | Significant? |
|-----------|--------|---------|--------|-------------|
| 100 yr | 90 | 0.24 (LC21) | 0.02 | Marginal, wrong sign |
| 250 yr | 36 | 0.14 | >0.3 | No |
| 500 yr | 18 | 0.14 | >0.5 | No |
| 1000 yr | 9 | 0.19 | >0.5 | No |

Zero of four proxies show significant negative correlation at any resolution. Lead-lag analysis (−500 to +500 years) produces no significant correlations. The initial r = −0.77 result was a binning artifact — the hardcoded 500-year bin averages used in the preliminary analysis do not replicate with the actual downloaded high-resolution paleoclimate data.

The temporal coincidence of the divergence collapse and the 4.2-kiloyear event is noted (Figure S6) (cf. Weiss et al. 1993; deMenocal 2001; Staubwasser et al. 2003) — the largest single-bin D drop occurs at 2150–2050 BCE — but the climate proxies show maximum drying rate at 3200–3500 BCE, not 2200 BCE, and the correlation does not survive formal testing at any resolution.

### 4.9 Trade Network Correlation

Radiocarbon date density within the Egypt-Iran corridor does not correlate with D: r = +0.04, p = 0.93 (Pleiades); r = +0.17, p = 0.69 (p3k14c). Overall corridor activity increases monotonically while D spikes and collapses. The divergence does not track general human activity levels.

### 4.10 Deep-Time Radiocarbon Analysis

**Table 8: Near-circle radiocarbon date enrichment by period (p3k14c, 200 km threshold)**

| Period (14C BP) | Approx. calendar | Near-circle enrichment | Z |
|----------------|-----------------|----------------------|---|
| 15,000–12,800 | Pre-YD | 0.3–1.7× | Noisy (small N) |
| 12,800–11,600 | Younger Dryas | 1.2× | +0.71 (not significant) |
| 10,500–10,000 | Early Holocene | **5.75×** | **+2.66** |
| 9,000–8,500 | Early Neolithic | **5.44×** | **+2.45** |
| 8,500–8,000 | Neolithic | 4.10× | +1.55 |
| 7,000–5,000 | Mid-Holocene | 1.5–2.0× | Baseline |
| 5,000–3,000 | Pre-monument | 1.5–2.5× | Building |
| 3,000–2,500 BCE | Monument spike | 3.0–5.0× | Significant |

The circle shows no anomalous signal during the Younger Dryas (Z = +0.71; 43% of 100 random circles show a more extreme YD anomaly). There is no evidence of distinctive disruption during this period.

A striking enrichment emerges at the Pleistocene-Holocene boundary. During 10,500–8,500 BP (~8,500–6,500 BCE), radiocarbon dates within 200 km of the circle are 4–6× overrepresented relative to the global p3k14c baseline. The enrichment declines to ~1.5–2× during the mid-Holocene before rising again to the monument-construction spike at 3000–2500 BCE.

**Temporal refinement.** Splitting the early Holocene dates by the PPNA/PPNB settlement break documented by Borrell et al. (2015) at ~10,200 BP reveals that 17 of 19 near-circle sites are post-break PPNB settlements. Only 2 sites (Ramat Harif, Maaleh Ramon West — both Negev PPNA desert sites) predate 10,200 BP. The strongest single bin (9,000–8,500 BP, 5.44×, 152 dates) corresponds to middle/late PPNB — the period of maximum village proliferation and long-distance obsidian exchange. The enrichment tracks the post-break demographic boom, not the earlier PPNA occupation.

**Fertile Crescent control (p3k14c).** Among 100 random great circles tested on the global p3k14c database, the mean early Holocene enrichment was 1.29× (median 1.05×). Only 1 of 100 exceeded 4×. Among 9 random circles that pass through the Fertile Crescent, enrichment ranged from 0.88× to 2.71× (mean 1.66×). The Alison circle's 4.52× enrichment is Z = 5.26 above the Fertile Crescent-passing mean.

**Table 9: Early Holocene enrichment (10,500–8,500 BP) — Alison circle vs. controls (p3k14c)**

| Comparison | Mean enrichment | Range | Alison circle |
|-----------|----------------|-------|---------------|
| All 100 random circles | 1.29× | 0.1–4.01× | 4.52× (99th percentile) |
| 9 Fertile Crescent-passing circles | 1.66× | 0.88–2.71× | 4.52× (Z = 5.26 above mean) |

**NERD replication (regional database).** To test whether the enrichment reflects a specific band within the Fertile Crescent or the circle's passage through the Fertile Crescent relative to a global baseline, the analysis was repeated on the NERD database (11,072 radiocarbon dates restricted to the Near East). The Alison circle's early Holocene enrichment on NERD was 1.10× — indistinguishable from the regional baseline. Among 9 Fertile Crescent-passing circles, the Alison circle ranked 2nd (Z = 0.89).

The NERD result indicates that the 4–6× enrichment observed in p3k14c reflects the circle's passage through the Fertile Crescent relative to the global database, not a specific corridor within the Fertile Crescent. When the baseline is restricted to the Near East, the circle shows no preferential capture of early Holocene activity.

**Directional analysis.** A standard deviational ellipse computed for 143 unique early Holocene sites in the Near East (p3k14c, 10,500–8,500 BP) yields a major axis bearing of 14.3° (NNE–SSW, elongation ratio 1.44, bootstrap 95% CI: 8.6°–22.1°). The Great Circle's bearing at the ellipse center is 88.1° (E–W). The angular difference is 73.8° — the circle runs nearly perpendicular to the natural elongation of the early Neolithic settlement distribution, which follows the NNE–SSW arc of the Levant and upper Mesopotamia.

**Regional decomposition of the early Holocene enrichment.** The 4–6× enrichment is concentrated in a ~1,200 km stretch of the circle (the 0°–60° arc segment: Zagros through Levant/Negev to the Green Sahara), containing 17 of 19 near-circle early Holocene sites. All other segments are at or below baseline; the South American segment is depleted (Z = −0.89).

**Summary.** The Alison Great Circle captures more early Holocene activity than 99% of random circles on the global p3k14c database, driven by its passage through the Fertile Crescent. However, this enrichment does not replicate on a regional Near Eastern database (NERD), and the circle's bearing is perpendicular to the natural axis of early Neolithic settlement. The deep-time signal reflects the circle's multi-regional geometry — connecting the Fertile Crescent, South America, and Easter Island on a single path — rather than a specific corridor of habitation within the Near East.

### 4.11 Geophysical Property Scan

**Table 10: Geophysical properties of the Alison circle vs. 100 random circles**

| Property | Alison value | Random mean | Z | Significant? |
|----------|-------------|-------------|---|-------------|
| Geoid mean height | — | — | +1.22 | No |
| Geoid gradient | — | — | +2.33 | Nominal |
| Magnetic anomaly (mean) | — | — | −0.09 | No |
| Magnetic anomaly (variance) | — | — | −1.06 | No |
| Ocean fraction | 56.4% | 72.6% | −2.05 | Nominal |
| Coastal fraction | — | — | +1.76 | No (marginal) |
| Elevation variance | — | — | −0.30 | No |
| Ocean current alignment | 61% | 17% | +2.08 | Nominal |

Three of thirteen tests reach nominal significance: geoid gradient (Z = +2.33), reduced ocean fraction (Z = −2.05), and ocean current alignment (Z = +2.08). None survive Bonferroni correction across 13 tests. No crustal magnetic signature was detected.

The nominally significant results are consistent with the circle tracing continental margins and navigable maritime routes rather than a specific geophysical feature. The circle crosses less open ocean than random (56.4% vs. 72.6%), traverses steeper geoid transitions (characteristic of continental-oceanic boundaries), and aligns with major surface currents where it does cross water (61% of segments within 30° of a current, driven primarily by the Atlantic North Equatorial Current). These properties describe a route favorable to human habitation and maritime navigation, not a geophysical force organizing site placement.

### 4.12 Regional Decomposition

**Table 7: Monument-settlement divergence by region**

| Region | D | Database | Note |
|--------|---|----------|------|
| Egypt/Levant | 8.68 | Pleiades | Independently significant |
| Iran/Mesopotamia | 2.60 | Pleiades | Weak positive |
| South Asia | 1.07 | Pleiades | Inconclusive (small N) |
| South America | −1.25 | Peru Ministry | Null (1,733 mon / 275 set) |
| South America | 0.19 | SAAID | Null (209 mon / 530 set) |
| South America | 0.13 | Pleiades/p3k14c | Null |
| Western Anatolia | 0.12 | LuwianSiteAtlas | Null (0/483 within 500 km) |

The divergence is concentrated in Egypt/Levant (D = 8.68). Under the initial classification used in this decomposition, South America shows no divergence on three independent databases (D = −1.25 to 0.19). However, as shown in Section 4.24, this result is classification-sensitive: a broader monument definition applied consistently across New World databases reverses the finding (D = 4.69 to 6.50). Western Anatolia — geographically adjacent to the Levant and a major Bronze Age trading hub — shows zero overlap with the circle.

XRONOS confirms the Egypt dependence: D = +3.28 (full), D = −14.05 (without Egypt).

### 4.13 Monument Orientation

Monument orientations within the Egypt-Iran corridor were tested for correlation with the circle's local bearing (~85° at Egyptian sites) and with stellar azimuths at 2500 BCE (Figure S4). All Egyptian monuments are oriented N-S, perpendicular to the circle's E-W bearing. Rayleigh test: p = 0.75 (no clustering). Monument orientations do not encode the circle's geometry.

### 4.14 Corridor Width

The monument-to-settlement ratio was computed at 5 km intervals from the circle within the Egypt-Iran corridor (see also Figure 6 for the arc-position distribution). The ratio peaks at 10 km (1.43 — monuments outnumber settlements), with the sharpest signal in the 5–10 km band (2:1 monument-to-settlement ratio). By 25 km the ratio reaches ~1:1; by 100 km it is 1:2. The effective enrichment zone is approximately 10 km, narrower than the 25 km estimated in earlier analyses.

At the Memphis necropolis specifically, the Great Circle crosses the ancient Ahramat Branch of the Nile (Sheisha et al. 2022) at a perpendicular angle near Abu Sir. The pyramid cluster at Giza-Abu Sir-Saqqara sits 3–7 km from the circle at this crossing point; Dahshur at 15 km is already drifting away. The geometry is a crossing, not a parallel track — the circle and the ancient river channel intersect at approximately 107°.

### 4.15 Circle Family Analysis

Of 100,000 random great circles, 66 pass within 200 km of all three anchor sites (Giza, Nazca, Easter Island) — a frequency of 1 in 1,515. These 66 circles show mean early Holocene enrichment of 3.31× on p3k14c (vs. 1.47× for random circles), indicating the enrichment is a property of this geometric family, not unique to the Alison circle.

Within this family of 60 qualifying circles (computed on a 1° pole grid plus the Alison pole), the monument-settlement divergence has mean D = 1.18 (SD = 4.13, range −7.14 to 9.82). The Alison circle's D = 8.64 ranks 5th of 60 (93.2nd percentile, +1.72σ). The top 4 circles have poles within ~2° of the Alison pole, defining a ~4° "hot zone" that traces nearly the same path through the Old World. The divergence has a broad peak centered near the Alison pole, not a sharp spike at one specific geometry (Figure S1).

### 4.16 Alternative Site Triplets

To test whether the Alison circle's divergence is a generic property of great circles through famous archaeological sites, eight alternative triplets were tested.

**Table 11: Monument-settlement divergence for alternative famous-site triplets**

| Triplet | Max miss | D | vs. random |
|---------|----------|---|------------|
| **Alison (Giza+Nazca+Easter Is.)** | **6.5 km** | **8.64** | **—** |
| Athens + Persepolis + Varanasi | 305 km | 6.08 | 100th %ile |
| Cahokia + Carnac + Luxor | 34 km | 4.73 | 100th %ile |
| Teotihuacan + Giza + Angkor Wat | 2,066 km | 3.66 | 99th %ile |
| Stonehenge + Giza + Persepolis | 620 km | 2.87 | 98th %ile |
| Machu Picchu + Stonehenge + Angkor | 144 km | 2.42 | 97th %ile |
| Giza + Angkor Wat + Machu Picchu | 75 km | 1.77 | 96th %ile |
| Tiwanaku + Petra + Borobudur | 1,064 km | 1.01 | 98th %ile |
| Giza + Göbekli Tepe + Mohenjo-daro | 380 km | 0.23 | 86th %ile |

All eight alternatives show positive D, confirming that monument-favoring divergence is a baseline property of great circles through famous archaeological sites (which tend to traverse the Pleiades-dense Mediterranean and Near East). However, no alternative triplet reaches D > 7. The strongest alternative (Athens-Persepolis-Varanasi, D = 6.08) falls 42% below the Alison circle (D = 8.64). The Alison circle produces the strongest monument-settlement divergence of any tested great circle through major archaeological sites.

The geometric fit is also exceptional: the Alison circle passes within 6.5 km of all three anchors, while most triplets cannot be fit below 100 km.

### 4.17 Pre-Younger Dryas Corridor Test

**Table 12: Pre-YD corridor enrichment (94,181 merged dates, 50 km threshold)**

| Metric | Observed | Random mean ± SD | Z-score | Percentile |
|--------|----------|-----------------|---------|------------|
| Pre-YD dates (>12,800 BP) | 10 | 127.1 ± 248.9 | −0.47 | 32.8% |
| Deep pre-YD (>15,000 BP) | 3 | 110.8 ± 220.2 | −0.49 | 21.8% |
| Late Glacial (12,800–15,000 BP) | 7 | 16.4 ± 31.0 | −0.30 | 60.1% |
| Habitability-matched pre-YD | 10 | 113.4 ± 104.7 | −0.99 | 11.7% |

Ten radiocarbon dates fall on the corridor before 12,800 BP across all three databases (Figure S5). Random circles average 127. Among habitability-matched corridors (circles through equally populated regions, using HYDE 3.2 population estimates [Klein Goldewijk et al. 2017]), the Alison corridor drops to the 12th percentile.

**Regional decomposition.** Of the 10 pre-YD dates, 6 are in Egypt/Levant and 4 in Iran/Central Asia. Zero pre-YD dates fall in Peru/Amazon, North Africa/Mediterranean, SE Asia, or any other segment. The pre-YD activity is entirely explained by the circle's passage through the Fertile Crescent — the most intensively studied Late Pleistocene region on Earth.

**Continuity across the YD boundary.** 4 of 190 on-corridor site clusters (2.1%) have dates spanning both sides of the YD, compared to 1,271 of 22,805 off-corridor clusters (5.6%). The corridor shows less continuity than the global average — the opposite of a "survivor route."

**YD disruption.** The corridor's crash ratio (post-YD onset / pre-YD activity) is 1.50, compared to a global ratio of 1.32 (Z = 0.11). The corridor tracks the global pattern; nothing special happened here at the YD onset.

**Candidate Atlantis regions.** Five proposed locations were tested: Atlantic west of Gibraltar (328 pre-YD dates, explained by known Iberian Upper Palaeolithic), Richat Structure (0 pre-YD), Persian Gulf basin (0 pre-YD), Sunda Shelf (0 dates of any period), Doggerland (1 pre-YD date). No candidate region shows anomalous pre-YD activity unexplained by conventional archaeology.

**SPD comparison.** The on-corridor summed probability distribution is flat at zero before the YD, then rises sharply after 10,000 BP, peaking at 7,000–8,000 BP. The off-corridor SPD shows gradual increase throughout. The on/off ratio is 0.66× pre-YD and 1.31× post-YD — the corridor is depleted before the YD and enriched after it.

**Sensitivity checks.** The null result is robust across thresholds: 25 km (Z = −0.45), 100 km (Z = −0.34), excluding North America (Z = −0.43).

The corridor was empty before the Younger Dryas. Whatever the Great Circle alignment represents, it is a post-glacial phenomenon. The corridor was colonized during the early Holocene as part of normal post-glacial migration, not inherited from a pre-YD civilization.

### 4.18 Giza Longitude Grid Test

**Table 13: Giza longitude grid test (508,000 sites, 6 tests)**

| Test | Result | Interpretation |
|------|--------|---------------|
| Sites at 36° multiples from Giza | Z = 27–133 | Explained by continental geography |
| Giza vs. 1,000 random reference longitudes | 71st–79th percentile | Not special |
| Precessional angles vs. rotated/random grids | Z = 0.45–0.94, p = 0.20–0.27 | Not significant |
| Claimed site pair precision | Angkor 72.733° (not 72°); Nazca 111.2° (not 108°); Persepolis 18.8° | Imprecise |
| Monument vs. settlement alignment | Settlements 15.4% vs. monuments 9.3% | Opposite of prediction |
| Age-alignment correlation | r = −0.57; oldest bin enrichment = 1.00 | Youngest sites most aligned |

Sites do cluster at 36° multiples from Giza, but this is entirely explained by the clustering of sites on continents. Giza is unremarkable as a reference longitude — hundreds of random points produce comparable results. The claimed precessional significance of 36° spacing is not supported: precessional angles perform no differently from arbitrary spacings. The specific site-pair claims are imprecise (Angkor is 72.733° from Giza, not 72.000°; Nazca is 111.2°, not 108.000°). Settlements are more grid-aligned than monuments (the opposite of deliberate monument placement), and the oldest sites show zero grid alignment while younger sites show more (the opposite of an ancient encoding).

The Giza longitude grid is not supported.

### 4.19 Pyramid Date-Distance Correlation

Construction dates and Great Circle distances were compiled for 26 Egyptian pyramids (Figure S3) with reliable chronology (Djoser through Amenemhat III; chronology following Lehner [1997], Stadelmann [2001], and Dee et al. [2013]).

**Table 14: Pyramid date-distance correlation**

| Metric | Value |
|--------|-------|
| Pearson r | −0.428 |
| Pearson p | 0.029 |
| Spearman ρ | −0.418 |
| Spearman p | 0.033 |
| N | 26 |

Earlier pyramids are significantly closer to the Great Circle. The Step Pyramid of Djoser (2630 BCE) is 6.1 km from the circle. The Layer Pyramid of Khaba (2600 BCE) is 1.7 km — the closest of any pyramid. The Abusir group (2445–2494 BCE) clusters at 3–6 km. Middle Kingdom pyramids (Dynasty 12) drift to 12–74 km. Migration analysis in 50-year bins shows the necropolis centroid moving progressively south and away from the circle over seven centuries.

This temporal-spatial gradient is not explained by terrain features (the desert-floodplain boundary does not shift), the Ahramat Branch (which is stationary), or any of the geographic features in the ML decomposition (Section 5.10). It represents a pattern within the Memphis necropolis that tracks the Great Circle geometry and weakens over time.

### 4.20 Anti-Circle Negative Control

To confirm that the divergence is specific to the Alison geometry and not an artifact of the methodology, three control circles were tested: two perpendicular to the Alison circle and one in the opposite hemisphere.

| Circle | Pole | D |
|--------|------|---|
| **Alison** | **(59.68°N, −138.65°W)** | **8.16** |
| Perpendicular 1 | (59.68°N, −48.65°W) | −2.00 |
| Perpendicular 2 | (30.32°N, 41.35°E) | 0.00 |
| Opposite hemisphere | (0°N, 131.35°E) | 0.18 |

The maximum control |D| = 2.00, approximately 4× smaller than the Alison D = 8.16. The divergence is geometry-specific.

### 4.21 Elevation and River Crossings

The Alison circle crosses 9 major rivers (Strahler order ≥ 5), compared to a random-circle mean of 4.6 (92nd percentile, Z = 1.90). Mean elevation (−1,709 m) is higher than random (−2,405 m, 86th percentile), reflecting greater land coverage. Walkability (56th percentile) and elevation variance (41st percentile) are unremarkable. Monument cluster regions sit at mean +709 m elevation, consistent with placement on land.

### 4.22 Temporal Wave Analysis

Among 25 on-circle monuments with reliable construction dates, arc position and construction date show a significant correlation (r = −0.45, p = 0.028), with older monuments at lower arc positions (Near East) and younger monuments at higher positions (South America, Pacific). However, 18 off-circle monuments show a similar pattern (r = −0.49, p = 0.124), suggesting this reflects the global east-to-west age gradient of monument construction (Near East oldest, Americas youngest) rather than a propagation wave specific to the circle.

### 4.23 Obsidian Trade and Linguistic Diversity

Zero of 34 compiled obsidian source-to-destination trade routes in the Neolithic Near East cross the Great Circle. The trade network operates at latitudes 32–41°N, while the circle passes through 28–30°N in this region. Trade routes and the circle are spatially independent.

Linguistic diversity on the corridor (7.3 language families per 100 languages) is higher than off-corridor (4.9), consistent with the circle crossing multiple language-family boundaries (Afro-Asiatic, Indo-European, Dravidian, Austronesian, Quechuan) rather than following any single linguistic zone.

### 4.24 Hemisphere Flip — New World Independent Replication

The full analysis pipeline was run using exclusively New World databases to test whether the divergence is an Old World artifact.

| Database | Hemisphere | Z-score | Enrichment | D |
|----------|-----------|---------|------------|---|
| Peru Ministry | NW | 25.41 | 11.76× | 6.50 |
| XRONOS | NW | 9.12 | 2.28× | 4.69 |
| XRONOS | OW | 8.56 | 1.83× | — |
| p3k14c | NW | 8.78 | 2.53× | 4.99 |
| p3k14c | OW | 3.80 | 1.43× | 3.39 |

All three New World databases show strong positive divergence (D = 4.69 to 6.50) — monumental sites cluster on the Great Circle significantly more than domestic sites, independently replicating the Old World pattern. The temporal profile differs: New World enrichment peaks at 1000-1500 CE (Mesoamerican/Andean monumental period) while Old World enrichment peaks at 2750-2500 BCE. The same structural pattern — monuments near circle, settlements not — operates with different temporal peaks on different continents.

This result revises the earlier South American null reported in Section 4.5 (Peru Ministry D = −1.25, SAAID D = 0.19). The discrepancy reflects differences in classification methodology: the earlier test used a more restrictive monument definition. Under the full hemisphere flip pipeline with consistent classification, the New World divergence is comparable to or stronger than the Old World.

### 4.25 Anti-Divergence Circle

The great circle maximizing inverse divergence (settlement enrichment, monument depletion) was identified by scanning ~2,000 pole positions.

Under ratio-based optimization (settlement_count / monument_count), the anti-divergence pole is orthogonal to the Alison pole at 90.5° angular separation. Under difference-based optimization (D_settlement − D_monument), the poles are separated by only 10.7°. The discrepancy reflects different weight given to absolute vs relative site counts.

The 10.7° result from difference-based optimization is mechanistically informative: the monument-settlement ratio is a continuous function of the circle's local crossing angle through dense archaeological regions. A small change in tilt shifts the capture ratio from monument-dominated to settlement-dominated — consistent with a local geographic mechanism (desert plateau vs floodplain at each valley crossing) rather than a global alignment property.

### 4.26 Directional Diffusion and Crop Spread

Two tests of whether the corridor functioned as a diffusion channel:

**Granger causality along the arc.** Temporal activity at adjacent arc segments was tested for directional predictive relationships. No significant Granger causality was detected — monument building along the corridor was parallel and independent, not propagating from one segment to the next.

**Archaeobotanical diffusion.** The wheat-barley-lentil crop package spread from the Fertile Crescent at known rates. Spread bearings were compared to the Great Circle bearing at each waypoint. Mean angular offset: ~45° (random expectation). Crop diffusion runs perpendicular to the corridor (along latitude bands), not along it. The circle was not a crop diffusion route.

### 4.27 36° Periodicity — Retraction

An earlier analysis found that inter-site arc distances among 6 hand-picked anchor sites show 36° periodicity. This was tested on data-driven clusters identified by kernel density estimation along the arc. The periodicity does not replicate (Rayleigh p > 0.05). The original finding was an artifact of site selection and is retracted.

---

## 5. Discussion

### 5.1 Summary

This study establishes twelve principal findings regarding the Alison Great Circle alignment reported in Paper I:

1. **Predictive validation.** Sites and dates published after 2001 cluster more strongly than pre-2001 data, the opposite of what post-hoc selection predicts. However, the signal is geographically concentrated in five regions that overlap with the circle's defining sites. The alignment is not solely dependent on sites known in 2001, but it has not been shown to predict activity in previously unexplored regions.

2. **Temporal specificity.** The monument-settlement divergence onsets sharply at 2750 BCE and collapses by 2000 BCE. This temporal precision is inconsistent with time-invariant recording bias, which predicts constant divergence across all periods.

3. **Geographic concentration.** The 2750–2500 BCE spike is driven entirely by Egyptian Old Kingdom pyramids in the Memphis necropolis — one cluster spanning ~70 km. The divergence is not a distributed corridor phenomenon.

4. **Preservation survival.** Adding documented buried settlements to the analysis increases the divergence. The geometry of the Nile Valley — monuments on the desert plateau, settlements on the floodplain — means the missing settlements are far from the circle, not near it.

5. **Climate correlation.** An initial analysis suggested the divergence correlates with Soreq Cave δ18O (r = −0.77 at 500-year lag). High-resolution reanalysis using four independent paleoclimate proxies at multiple resolutions (100yr to 1000yr) found no significant correlation. The initial result was a binning artifact. The temporal coincidence with the 4.2-kiloyear event is noted but does not survive formal testing.

6. **Deep-time enrichment.** Radiocarbon date density near the circle was 4–6× enriched during the early Holocene (10,500–8,500 BP) on the global p3k14c database — 99th percentile among random circles and Z = 5.26 above Fertile Crescent-passing controls. However, this enrichment does not replicate on the regional NERD database (1.10×), indicating it reflects the circle's multi-regional geometry rather than a specific corridor within the Near East.

7. **Pre-Younger Dryas null.** The corridor is empty before the Younger Dryas: 10 radiocarbon dates across 94,181 merged entries spanning 25,000 years (Z = −0.47, 33rd percentile). The corridor is a post-glacial phenomenon.

8. **Longitude grid falsification.** The proposed Giza-centered precessional longitude grid shows no support across six tests on 508,000 sites.

9. **ML geographic decomposition.** The monument-settlement divergence is predictively consistent with measurable geographic features (elevation, site clustering, coastline distance, latitude). Distance to the Great Circle ranks last (#10/10) in SHAP importance. The model loses no predictive power without it.

10. **Pyramid date-distance gradient.** Earlier pyramids are closer to the Great Circle (r = −0.43, p = 0.029). *Subsequently resolved as artifact* (Section 5.15): a bearing scan shows 31 of 36 bearings produce comparable gradients, indicating generic radial expansion from the necropolis centroid rather than a circle-specific pattern.

11. **Anti-circle control.** Three perpendicular/opposite-hemisphere circles produce max |D| = 2.0, confirming the divergence is geometry-specific (Alison D = 8.16).

12. **Systematic search.** Twelve independent high-D clusters exist globally. The global maximum (D = 24.3) is Italy-dominated. The Alison circle's unique property is its triple-anchor collinearity: the only pole in the scan passing within 12 km of Giza, Nazca, and Easter Island simultaneously.

13. **New World replication.** The monument-settlement divergence replicates independently in the New World (Peru Ministry D = 6.50, Z = 25.41; XRONOS NW D = 4.69; p3k14c NW D = 4.99). The temporal peak differs (1000-1500 CE vs 2750 BCE) but the structural pattern is identical. The divergence is not an Old World artifact.

14. **Publication bias control.** On-corridor sites are not more researched than off-corridor sites (p3k14c p = 0.142, XRONOS p = 0.595, Pleiades reversed). Research-intensity-weighted enrichment increases from 3.10× to 3.33×.

15. **Monument type decomposition.** The divergence is primarily mortuary: 71% of the 2750 BCE spike sites are tombs, pyramids, and necropoleis. Ceremonial sites show modest enrichment (2.06×). Administrative and entertainment sites show none.

16. **Anti-divergence proximity.** The circle maximizing inverse divergence (settlement enrichment) is separated from the Alison pole by only 10.7° under difference-based optimization (90.5° under ratio-based). The monument-settlement ratio is a continuous function of the circle's crossing angle through archaeological landscapes.

17. **Directional diffusion null.** No Granger causality along the arc. No crop-diffusion alignment. 36° periodicity falsified on data-driven clusters. The corridor did not function as a diffusion channel.

### 5.2 The Predictive Validation: What It Shows and What It Doesn't

The predictive validation is the most direct test of whether the circle was fit to its data. The result — stronger post-2001 clustering — rules out the simplest selection bias explanation: Alison did not draw a line through a few famous sites that subsequent data would regress away from. The pattern held and strengthened as new data accumulated.

However, the fieldwork intensity analysis reveals that the post-2001 signal is driven by continued research in the same regions the circle was designed to connect. Easter Island alone contributes 206 of 583 near-circle dates. Removing five grid cells (0.85% of data) eliminates the signal. The circle predicts continued archaeological findings in Egypt, Peru, the Levant, Iran, and Easter Island — but these are precisely the regions where the world's most intensive archaeological programs operate, independent of any alignment theory.

The honest interpretation: the predictive validation confirms the alignment is real in the sense that the regions it passes through are genuinely rich in archaeological material. It does not confirm that the circle captures a pattern that extends beyond those already-known regions.

### 5.3 The Temporal Signature Rules Out Recording Bias

The sharpest objection to Paper I is that the monument-settlement divergence reflects differential recording: monuments are more visible, better documented, and more likely to appear in databases than domestic sites. This objection predicts a time-invariant divergence — the recording asymmetry applies equally to monuments from 4000 BCE and 3000 BCE.

The observed divergence is sharply time-variant. It is absent before 2750 BCE, spikes in a single 250-year bin, and collapses within 500 years. Recording bias does not have sharp temporal edges. Construction traditions do. The 2750–2500 BCE peak coincides precisely with the era of pyramid construction at Giza (~2560 BCE), Saqqara (~2630 BCE), and Dahshur (~2600 BCE).

### 5.4 The Preservation Test

Egyptian archaeology has long recognized that stone monuments survive while mudbrick settlements dissolve (Moeller 2016; Kemp 2018). The EES Delta Survey documents 783 sites of which ~230 are destroyed and ~180 overbuilt (Spencer 2024). Parcak (2016) estimates that less than 0.001% of Delta archaeological volume has been excavated.

This study tests the preservation objection quantitatively for the first time. The 58 documented buried settlements are overwhelmingly in the Nile Delta and Valley floor — geographically distant from the great circle, which crosses the Nile Valley at a diagonal hitting the desert plateau edge near Memphis. Adding them increases D because they add settlement mass far from the circle. The Monte Carlo estimate (+350 settlements) confirms that even generous assumptions about missing settlement sites do not threaten the divergence.

### 5.5 The Memphis Concentration

The site-level decomposition of the temporal spike reveals that the monument-settlement divergence, at its peak, is driven by one geographic cluster: the Memphis necropolis. All 38 Pleiades monuments in the 2750–2500 BCE bin are within ~70 km of each other. The "corridor" framing from Paper I — which emphasized six clusters across four continents — overstates the geographic breadth of the divergence signal. The overall alignment is real (Z = 25.85 on the primary database, replicated seven times). But the monument-specific divergence within it is geographically narrow.

This does not invalidate the finding. The Memphis necropolis is the densest concentration of monumental construction in the ancient world. The question is not whether it falls near the circle — it does — but whether its proximity is coincidental or causally related to the circle's geometry. The circle does not track the desert-floodplain boundary (Section 4.7), does not correlate with monument orientations (Section 4.13), and has no identified geophysical correlate. The geometry of the alignment remains unexplained.

### 5.6 Climate Correlation: A Cautionary Note

The initial Soreq Cave correlation (r = −0.77, p = 0.025 at 500-year lag on 7–8 bins) suggested a mechanistic link between favorable climate and monument construction. High-resolution reanalysis using four independently downloaded NOAA paleoclimate proxies at 100-year, 250-year, 500-year, and 1000-year resolution found no significant correlation at any scale or lag. The initial result was a binning artifact — the hardcoded bin averages used in the preliminary analysis did not replicate with the actual downloaded data.

This serves as a methodological caution: correlations based on fewer than 10 temporal bins are extremely sensitive to binning choices and should not be treated as confirmed findings, even when p-values are below conventional thresholds. The publication of this negative replication is itself informative, demonstrating the importance of multi-resolution and multi-proxy verification for any claimed paleoclimate-archaeological correlation.

### 5.7a The Circle Family and the Uniqueness of the Alison Geometry

The Alison circle belongs to a family of ~66 great circles that thread within 200 km of Giza, Nazca, and Easter Island — a frequency of 1 in 1,515 among random circles. Within this family, the Alison circle's divergence (D = 8.64) ranks 5th of 60 (93rd percentile), with the top 4 circles having poles within ~2° of the Alison pole. The signal has a broad peak, not a sharp spike at one specific geometry.

However, an alternative-triplets analysis demonstrates that the Giza-Nazca-Easter Island family is exceptional. Eight alternative triplets of famous archaeological sites were tested; all showed positive D (confirming that monument-favoring divergence is a baseline property of famous-site geometry), but none reached D > 7. The strongest alternative (Athens-Persepolis-Varanasi, D = 6.08) falls 42% below the Alison circle. The Alison circle produces the strongest monument-settlement divergence of any tested great circle through major archaeological sites, with geometric fit an order of magnitude tighter than most alternatives (6.5 km vs. 34–2,066 km).

The implication is that the collinearity of Giza, Nazca, and Easter Island is geometrically rare (1 in 1,515 circles) and produces a monument-settlement divergence unmatched by any other tested combination of famous sites. The divergence is not unique to the precise Alison pole — nearby poles produce comparable results — but the Giza-Nazca-Easter Island axis itself is distinctive.

### 5.7b The Deep-Time Signal

The early Holocene enrichment (4–6× at 10,500–8,500 BP) is perhaps the most consequential finding of this study, because it reframes the relationship between the circle and the monuments built along it.

The monument-settlement divergence of 2750–2250 BCE was not constructed in an arbitrary location. The Memphis necropolis, and the broader corridor the circle traces, had been a locus of disproportionate human activity for over 6,000 years before the first pyramid was built. The early Holocene enrichment coincides with the Neolithic transition in the Fertile Crescent — the invention of agriculture, the establishment of sedentary communities, the beginning of the demographic expansion that eventually produced the first state-level societies.

A natural objection is that any great circle through the Fertile Crescent would show early Holocene enrichment. This objection was tested directly. Among 9 random circles passing through the Fertile Crescent, early Holocene enrichment ranged from 0.88× to 2.71× (mean 1.66×). The Alison circle's 4.52× enrichment is Z = 5.26 above this mean (Table 9). Only 1 of 100 random circles — regardless of path — exceeded 4× enrichment. The deep-time signal is not a generic property of Near Eastern geography. It is specific to this circle.

The temporal structure of near-circle activity shows two phases on the global p3k14c database: elevated human presence during the early Holocene (10,500–8,500 BP, 4–6× enrichment relative to the global baseline), followed by decline to baseline during the mid-Holocene, and then the sharp monument-specific spike during 2750–2250 BCE.

However, three additional analyses complicate the "inherited corridor" interpretation.

First, the early Holocene enrichment does not replicate on the NERD database (11,072 radiocarbon dates restricted to the Near East). NERD shows 1.10× enrichment — indistinguishable from the regional baseline. The 4–6× enrichment in p3k14c reflects the circle's passage through the Fertile Crescent relative to a global database, not a specific band within the Fertile Crescent. When the comparison baseline is already restricted to the Near East, the circle shows no preferential capture of early Holocene activity.

Second, a standard deviational ellipse computed for 143 early Holocene sites in the Near East shows the natural distribution is elongated NNE–SSW (major axis 14.3°, bootstrap 95% CI: 8.6°–22.1°), following the Levantine-Mesopotamian arc. The Great Circle's bearing through this region is 88.1° (E–W) — nearly perpendicular (73.8° angular difference). The circle cuts across the natural axis of early Neolithic settlement rather than following it.

Third, the enrichment tracks the post-10,200 BP PPNB demographic expansion specifically (Borrell et al. 2015), not the earlier PPNA. Of 19 near-circle early Holocene sites, 17 post-date the PPNA/PPNB settlement break. The strongest bin (9,000–8,500 BP, 5.44×) corresponds to middle/late PPNB maximum village proliferation.

These results require a revised interpretation. The Fertile Crescent control (Z = 5.26 above FC-passing circles) remains valid on the global p3k14c database — the Alison circle captures more early Holocene activity than other circles through the same region when both are measured against a global baseline. But this advantage reflects the circle's multi-regional geometry (connecting the Fertile Crescent, South America, and Easter Island on a single path) rather than a specific corridor within the Near East. The "inherited corridor" framing overstates the evidence; the more defensible claim is that the circle threads through multiple independently significant regions of human activity.

The Younger Dryas period (12,800–11,600 BP) shows no anomalous signal along the circle (Z = +0.71), inconsistent with proposals that a catastrophic event during this period destroyed a civilization concentrated along this path (cf. Hancock 1995). This null result is dramatically strengthened by the pre-YD corridor test (Section 4.16): across 94,181 merged radiocarbon dates spanning 25,000 years, only 10 fall on the corridor before 12,800 BP — fewer than the random-circle average of 127. The corridor sits at the 33rd percentile among random circles and the 12th percentile among habitability-matched corridors. Site continuity across the YD boundary is lower on-corridor (2.1%) than off-corridor (5.6%), the opposite of what a "survivor route" would predict. No candidate "Atlantis" region — including the Persian Gulf (which lies on the corridor and was dry land during the LGM) — shows anomalous pre-YD activity.

The corridor's temporal profile is now clear: empty before the Younger Dryas, colonized during the early Holocene as part of normal post-glacial migration, and monumentalized over the following 5,000 years. The "inherited corridor" framing from earlier analyses overstates the evidence; the corridor was not inherited from a pre-existing civilization. It was freshly populated by humans moving into habitable land as ice retreated — the same process that populated every habitable region on Earth. What distinguishes this corridor is not when it was colonized but what was built on it afterward.

### 5.7c The Giza Longitude Grid

Hancock (1998) proposed that ancient sites encode precessional knowledge through their longitude separations from Giza, specifically at multiples of 36° (and fractions thereof: 72°, 108°, 144°). This is distinct from the Alison Great Circle alignment and represents Hancock's own specific geodetic claim.

Six tests on 508,000 sites find no support (Section 4.17). The clustering of sites at 36° multiples is entirely explained by continental geography — sites cluster on continents, and continents happen to be spaced at intervals that approximate 36° multiples from Giza. Giza is unremarkable as a reference longitude; hundreds of random points produce comparable "grids." The claimed precessional significance is not supported: precessional angles perform no differently from arbitrary spacings.

The most diagnostic tests are the monument-settlement and age-alignment analyses. If sites were deliberately placed at grid longitudes, monuments should be more precisely aligned than settlements (which are geographically constrained). The opposite is observed: settlements are more grid-aligned (15.4% vs. 9.3%). Similarly, if the grid encodes ancient knowledge, the oldest sites should show the strongest alignment. The opposite is observed: the oldest bin (12,000–8,000 BP) shows enrichment of exactly 1.00 (chance), while younger sites show more alignment — consistent with increasing population filling continental landmasses, not deliberate placement.

The same methodology that confirms the Alison Great Circle alignment falsifies the Giza longitude grid. The two claims are distinct and receive opposite verdicts.

### 5.8 The South American Signal — Revised

The South American divergence result is classification-sensitive. Under a restrictive monument definition, three databases show null or negative results (Peru Ministry D = −1.25, SAAID D = 0.19, Pleiades/p3k14c D = 0.13; Section 4.12). Under the broader classification used consistently throughout the rest of this paper, three databases show strong positive results (Peru Ministry D = 6.50, XRONOS NW D = 4.69, p3k14c NW D = 4.99; Section 4.24).

Both results are reported transparently as a sensitivity analysis. The discrepancy arises because Peruvian site databases contain many records classified ambiguously between monumental and domestic function — the restrictive definition excludes these, while the broader definition includes them as monuments. The Old World divergence is not affected by this ambiguity because Egypt/Levant site classifications are more sharply defined in Pleiades.

Under the broader classification, the New World shows strong positive divergence with a different temporal peak (1000-1500 CE, driven by Mesoamerican and Andean late-period monuments) than the Old World (2750-2500 BCE, Egyptian pyramids). The same structural pattern — monuments near circle, settlements not — operates independently on both hemispheres, but this replication is contingent on the classification choice. This should be considered a promising but provisional finding pending more detailed site-level classification of South American databases.

### 5.9 Multiple Testing Correction

Benjamini-Hochberg (BH) false discovery rate correction (Benjamini & Hochberg 1995) was applied across all tests at FDR = 0.05. Of positive findings, 22 survive correction. The sole casualty is the desert belt proximity test (93rd percentile, uncorrected p = 0.07, BH-adjusted p = 0.07), which was already above the 0.05 threshold before correction and should be regarded as suggestive rather than significant.

All core findings survive correction trivially: database replications (BH p < 10⁻⁴), monument enrichment (BH p < 10⁻⁴), temporal spike (BH p < 10⁻⁴), preservation test (BH p < 10⁻⁴), collinearity tests (BH p < 0.001), anti-circle control (BH p = 0.001), New World replication (BH p < 10⁻⁴). The Orion shape match (BH p = 0.009), Nile constriction (BH p = 0.011), and 43,200 scale factor (BH p = 0.031) survive but are closer to the threshold. The full BH correction table is provided in Supplementary Table S2.

Two correction regimes were applied: (1) all-41 tests (conservative, includes null/falsification tests that are designed to fail), and (2) positive-only (23 tests where a positive finding is claimed). Results are concordant: the same 22 findings survive under both regimes.

### 5.10 Limitations

1. **Paleoclimate correlation.** The initial Soreq Cave correlation (r = −0.77, p = 0.025) did not replicate at higher resolution. Analysis using four independent proxies at four resolutions (100yr–1000yr) found no significant correlation. The initial result was a binning artifact, highlighting the fragility of correlations based on fewer than 10 temporal bins.

2. **Egypt dependence.** XRONOS confirms D collapses without Egypt (−14.05). The divergence may be an Egypt phenomenon, with the circle's passage through the Memphis necropolis as its primary expression.

3. **Predictive validation geography.** The post-2001 signal is concentrated in five grid cells overlapping with the circle's defining regions. The validation confirms the alignment is real but does not extend it to new regions.

4. **South American classification sensitivity.** The New World divergence flips from null (D = −1.25) to strongly positive (D = 6.50) depending on how ambiguous Peruvian sites are classified (Section 4.24, 5.8). Both results are reported. The positive result under broader classification is consistent with Old World patterns but should be treated as provisional until site-level classification standards are established for South American databases.

5. **Secondary late peak.** A divergence at 500 BCE–0 CE (D = 5.98) driven by settlement avoidance requires separate investigation.

6. **Pleiades temporal proxy.** Record creation dates in Pleiades track digitization, not discovery. The pre/post-2001 split is approximate.

7. **Pre-YD data sparsity.** Only 16,483 of 94,181 merged dates predate 12,800 BP. The Z = −0.47 null result (10 on-corridor dates vs. 127 expected) is robust at the aggregate level, but individual temporal bins contain very few dates. Submerged regions (Persian Gulf, Sunda Shelf) have no radiocarbon coverage, and absence of data in these areas does not constitute evidence of absence of activity.

8. **Longitude grid scope.** The grid test used the same merged database as other analyses. Hancock's claims about specific site pairs (Angkor, Nazca) predate high-precision GPS coordinates; the ~0.7° discrepancy for Angkor may or may not have been considered significant in 1998.

### 5.11 Machine Learning Geographic Decomposition

Using the methodology described in Section 3.11, three classifiers were trained on 21,779 Pleiades sites (4,012 monuments, 17,767 settlements) to test whether the Great Circle captures a dimension of monument placement not predicted by geography alone.

**Table 15: SHAP feature importance (XGBoost, global model)**

| Rank | Feature | Mean |SHAP| |
|------|---------|-------------|
| 1 | Cross-type distance | 2.34 |
| 2 | Site density | 0.50 |
| 3 | Elevation | 0.27 |
| 4 | Coastline distance | 0.23 |
| 5 | Latitude | 0.19 |
| ... | ... | ... |
| 10 | Distance to Great Circle | last |

Distance to the Great Circle ranks #10 of 10 features globally (Figure 5) and #6 of 10 within the Egypt-Levant-Iran corridor (20–38°N, 25–60°E). Ablation testing — comparing model AUC with and without GC distance — shows AUC changes of +0.0002 (global) and +0.0004 (corridor). The model loses no predictive power when the Great Circle feature is removed.

The monument-settlement divergence along the Great Circle is predictively consistent with known geographic features: elevation differences between monuments and settlements, spatial clustering patterns, coastline distance, and latitude. The circle traces a path along which these geographic factors are maximized, but it adds no independent predictive information beyond them.

**Methodological caveat.** SHAP values measure predictive importance, not causal relevance (Lundberg & Lee 2017; Molnar 2020). If Great Circle distance is correlated with the geographic features that rank higher, removing it will have minimal impact on prediction — not because it is irrelevant, but because its information is already captured by the correlated features. The ablation test proves that GC distance is redundant with other features in predicting monument vs. settlement classification; it does not prove that geography causally produces the divergence. A causal model (e.g., double machine learning) would be needed to distinguish whether geography mediates the divergence or merely correlates with it. The correct interpretation is that geographic features are sufficient to predict the divergence without needing the Great Circle as an additional predictor — a finding consistent with, but not proof of, a geographic explanation.

### 5.12 Systematic Great Circle Search

Using the methodology described in Section 3.12, a scan of all possible pole positions (2° grid, refined to 0.5° at 137 hot spots; Figures 2–3) identified 12 independent clusters of high monument-settlement divergence, separated by >20° of pole distance.

**Table 16: Top divergence clusters (systematic scan)**

| Cluster | Peak D | Location | Key sites captured |
|---------|--------|----------|-------------------|
| S California | 40.7 (Poisson) / 24.3 (MC) | 36°N, 118°W | Italy (Rome, Sardinia nuraghi) |
| Sumatra/SE Asia | 37.2 | 0°N, 102°E | Mediterranean + SE Asia |
| Pacific NW (Alison extended) | 32.6 | 46°N, 144°W | Egypt, Peru, Easter Island |
| Philippines/Taiwan | 28.4 | 22°N, 124°E | Pacific + Mediterranean |

The global maximum (D = 24.3 refined, pole at 36°N, 117.5°W) passes through central Italy (1,972 km from Giza), capturing the densest monument region in the Pleiades database (Roman temples, Sardinian nuraghi). It is not related to the Alison geometry.

**Anchor precision test.** Among all poles in the scan, only the Alison pole passes within 12 km of Giza (6.5 km), Nazca (11.8 km), AND Easter Island (10.3 km) simultaneously. The next-closest pole at (60°N, −144.5°W) achieves 5.7 km to Giza but misses Nazca by 318 km and Easter Island by 197 km. The Alison circle's triple-anchor collinearity is unique in the entire scan. Its distinctive property is not its divergence value but its geometric precision across three continents.

### 5.13 Open Questions

The ML decomposition resolves the divergence question: the monument-settlement divergence along the Great Circle is predictively consistent with geographic features (elevation, clustering, coastline distance, latitude). The circle traces a path where these features are maximized, but it adds no independent information.

The remaining question is geometric: why are Giza, Nazca, and Easter Island — three sites whose locations were determined by independent local geography — collinear along a path that maximizes monument-settlement geographic separation? This collinearity (1 in 1,515 random circles, 6.5 km precision, unmatched in a systematic scan of all possible poles) has moved from "mysterious signal requiring explanation" to "remarkable geometric coincidence along a geographically meaningful path."

Future investigations should address: whether the 12 independent high-D clusters identified in the systematic scan correspond to meaningful archaeological mega-regions; and whether underwater archaeological survey at the 133 km of submerged Great Circle path through the Persian Gulf (10 m mean depth) reveals early Holocene habitation sites.

### 5.14 Exhaustive Triplet Scan

Using the methodology described in Section 3.13, all C(149,3) = 540,274 triplets from 149 globally famous archaeological sites were evaluated for both geometric fit and monument-settlement divergence.

The Alison triplet ranks 785/540,274 for geometric fit (0.70 km mean, top 0.15%) and 51,373/540,274 for divergence (D = 4.58, 90th percentile). The combined metric D/(fit+1) ranks 1,982/540,274 (99.63rd percentile), equivalent to ~4,769 at the 1.3M scale. Approximately 1,228 triplets achieve both D > 8 and fit < 10 km. The collinearity is geometrically precise but not exceptional in the combined space. The top-scoring triplets are dominated by circles through Mediterranean monument clusters in the Pleiades database.

### 5.15 Terrain and Visibility Tests

Five mechanistic hypotheses were tested (Section 3.14) to determine whether the geographic decomposition reflects a specific landscape property:

| Test | Method | Result |
|------|--------|--------|
| Terrain transition zone | Elevation gradient, 200 random circles | 62nd percentile — not supported |
| Visibility/prominence | Topographic position index, 500km offset controls | 52nd percentile — null |
| Navigational waypoints | Coastal crossing prominence | p = 0.076 — marginal |
| Monument elevation differential | TPI comparison, monuments vs settlements | **Reversed**: monuments −46m, settlements −1m |
| Nile Valley constriction | Valley width profile, 29.0–30.5°N | **Confirmed**: GC crosses within 1 km of minimum width |

Near-circle monuments have lower topographic position (−46m) than settlements (−1m), contradicting the visibility hypothesis. Pleiades monuments in the Egypt-Levant-Iran corridor (temples, tombs, amphitheatres) tend to sit in valleys and lowlands. The one confirmed mechanism is the Nile Valley constriction: the Great Circle crosses the Nile at 29.91°N, within 1 km of the valley's minimum width (10.6 km at 29.90°N). Memphis was founded 6 km south of this geomorphological chokepoint.

### 5.16 Pyramid Date-Distance Gradient Resolution

The date-distance correlation (r = −0.428, p = 0.029, N = 26 pyramids) was subjected to three additional tests:

1. **Bearing scan**: 31 of 36 tested bearings through the necropolis centroid produce r < −0.3. The gradient is generic radial expansion, not direction-specific.

2. **Comprehensive data**: Expanding from the curated 26-pyramid list to all 117 p3k14c-dated sites in the Memphis region reduces r to −0.071 (p = 0.242). The Spearman correlation flips positive (+0.197, p = 0.033).

3. **Cross-cluster replication**: 0 of 6 monument clusters along the circle show significant "earlier = closer" gradients. Nazca shows the opposite (r = +0.398, p = 0.013): older sites are farther from the circle.

The date-distance gradient is an artifact of the curated pyramid list and generic necropolis expansion, not a circle-specific phenomenon.

### 5.17 Ambiguous Claims Assessment

Three claims previously classified as ambiguous were formally tested:

**43,200 scale factor.** The dual match (base perimeter × 43,200 ≈ Earth circumference; height × 43,200 ≈ polar radius) is achieved by only 3% of randomly dimensioned monuments — genuinely rare. However, 43,200 is not the optimal multiplier (43,492 is closer), and the Great Pyramid's overall numerological hit count sits at the 23rd percentile among random buildings. The pyramid plausibly encodes Earth's dimensions. The specific number 43,200 is post-hoc.

**Orion intentionality.** The Giza-Orion shape match ranks #5 out of 3,420 tested bright-star triplets (top 0.15%, p = 0.0074 against collinearity-matched random triangles). However, epoch variation in Procrustes distance is 0.000088 across 18,000 years — less than 1% of survey measurement uncertainty. All epochs produce indistinguishable matches. The shape correspondence is real and significant. The 10,500 BCE date claim has no geometric support.

**Pillar 43 stellar encoding.** The claimed animal-to-constellation mapping ranks #22 out of 24 possible permutations. Random carving arrangements match the claimed stellar positions better than the actual carvings (p = 0.99). Among 6,840 alternative constellation assignments, the claimed set ranks #689. The stellar encoding hypothesis is not supported.

### 5.18 Golden Section (φ) Distance Claims

Alison (c. 2001) claimed the ratio of Giza-Nazca to Angkor-Giza great circle distances equals φ (1.618). The observed ratio is 1.6227, within 0.29% of φ (p = 0.0083 for a single test on random three-point great circle divisions).

However, among all C(17,3) = 680 possible triplets of Alison's listed sites, 7 produce φ ratios with equal or smaller error. The best-matching triplet (Pyay–Khajuraho–Persepolis, 0.033% error) substantially outperforms Angkor-Giza-Nazca. Bonferroni-corrected p = 0.64, far above the 0.05 significance threshold. The φ ratio does not survive multiple testing correction.

### 5.19 Circle Narrative Comparison

An inventory of famous sites captured by each of the 12 high-D circles from the systematic scan reveals a paradox: the Alison circle is narratively unique but statistically unremarkable, while the statistically strongest circles are narratively uninteresting.

The Alison circle captures famous sites across 4 continents (Giza, Nazca, Machu Picchu, Easter Island, Persepolis, Petra, Mohenjo-daro, Ollantaytambo — 23 UNESCO World Heritage Sites within 50 km). No other high-D circle captures famous sites from more than 3 continents. The nearest competitor (Circle #10, pole 43°N/160°E, D = 15.8) captures Delphi, Parthenon, Chichen Itza, Palenque, and Pompeii across 3 continents.

The statistically strongest circles (D = 15–24) are dominated by Mediterranean Pleiades clustering, capturing 1,800–4,500 Pleiades sites that are overwhelmingly Roman temples, amphitheatres, and Sardinian nuraghi. Their high D values reflect database density, not archaeological significance.

This suggests the cultural narrative — not the geometry — is what makes the Alison circle compelling. The alignment is real and geographically meaningful, but its fame derives from WHICH sites it connects, not from the statistical strength of the connection.

### 5.20 Atlantis Candidate Bathymetry

Alison's three specific Atlantis candidate locations were tested against ETOPO1 bathymetry:

| Location | On circle? | Depth | Dry at LGM? | C14 dates | Verdict |
|----------|-----------|-------|-------------|-----------|---------|
| Cape Verde (south) | 29 km | −4,310 m | No | 0 | Rejected — abyssal plain |
| Mid-Atlantic 4°19'N, 41°30'W | 15 km | −3,738 m | No | 0 | Rejected — abyssal plain |
| Bay of Bengal ~20.6°N, 91.3°E | 0 km | −106 m | Yes | 0 | Viable but unexcavated |
| South China Sea ~10.4°N, 113.1°E | 0 km | −2,424 m | No | 0 | Rejected — deep basin |

Both Atlantic locations are on abyssal plain with no seamounts or shallow features within ±5° of longitude. The Bay of Bengal crossing (-106 m) would have been dry land during the LGM but has zero archaeological coverage. This represents a genuine data gap — the absence of evidence in submerged regions is not informative about the presence or absence of pre-YD human activity.

### 5.21 Off-Circle Monument Analysis

Fifty globally famous archaeological sites not on the Alison circle were analyzed for terrain signatures and alignment with the 12 high-D circles. Off-circle monuments have significantly different terrain profiles from on-circle terrain: TPI = −40 m (valleys/basins) vs +17 m for on-circle terrain (p = 0.0012). Off-circle famous monuments (Göbekli Tepe, Stonehenge, Angkor Wat, Teotihuacan) are placed in sheltered valley locations, while the Alison circle threads through elevated, exposed terrain.

Göbekli Tepe (770 km from the Alison circle, 314 km from the nearest high-D circle) is independent of all identified alignment patterns. Eighteen of 50 off-circle sites fall on at least one high-D circle at the 100 km threshold, but no single circle captures more than a few. A great circle through Newgrange and Knossos captures 10 of 50 off-circle sites within 100 km (Stonehenge, Karnak, Mycenae, Parthenon, Delphi, Valley of the Kings, Saqqara, Avebury; p < 0.0001 vs random), but this is a Mediterranean-to-Britain alignment through the Pleiades database's densest region and likely reflects database bias rather than a novel archaeological pattern.

### 5.22 Collinearity Resolution

Seven tests were conducted to determine whether the Giza-Nazca-Easter Island collinearity has a mechanistic explanation.

**Constraint propagation.** Easter Island is the only inhabited island within 200 km of the great circle defined by Giza and Nazca in the entire Pacific Ocean. The next inhabited island (Rapa Iti) is 293 km away. Only 6% of the circle's Pacific segment crosses land. Given Giza and Nazca, Easter Island is geographically forced as the third anchor.

**Continental position.** Random triplets with one point in North Africa, one on the Peruvian coast, and one on a Pacific island achieve ≤0.95 km collinearity only 0.034% of the time (p = 0.00034). The median random fit is 587 km. Continental geometry does not force this level of precision.

**Random civilization simulation.** Monument-weighted sampling from the three regions produces sub-1 km collinearity with p = 0.00025. The specific monument positions matter, not just regional membership.

**Geological source alignment.** The tectonic features responsible for the three site locations — the East African Rift (Chorowicz 2005), the Peru-Chile Trench (Cahill & Isacks 1992), and the Easter hotspot (Steinberger et al. 2004) — define a great circle pole only 0.8° from the Alison pole (cf. DeMets et al. 2010). Plate tectonics explains the approximate collinearity (~72 km fit). The sub-km human-site precision is not explained by geology.

**Atmospheric circulation.** 40% of the circle's land points fall within 5° of the 30th parallel (93rd percentile vs random circles). The circle tracks the subtropical desert belt — the Hadley cell boundary where arid conditions promote site preservation and dramatic terrain attracts monument placement.

**Equivalent site substitution.** Among 640 triplets of culturally equivalent sites (alternative Egyptian monuments × alternative Peruvian sites × alternative Pacific islands), Giza-Nazca-Easter Island produces the single best collinearity. Only 3 of 640 achieve <5 km fit.

**Synthesis.** The collinearity is a geologically scaffolded coincidence with genuinely unusual precision. Three layers of explanation account for it: (1) plate tectonics aligns the geological features that enabled three independent civilizations within ~1° of a shared great circle pole; (2) Easter Island is the only inhabited Pacific landfall on the resulting circle, forcing the third anchor; (3) the sub-km precision of the specific monument positions within each region is a ~1-in-3,000 rarity (p ≈ 0.0003) that represents a low-probability coincidence within the geologically determined framework, not evidence of ancient planning. The circle additionally tracks the subtropical desert belt (93rd percentile), linking climatically similar zones where preservation is favorable.

### 5.23 Summary of Findings and Open Questions

All primary findings have been addressed through 41 tests:
- The monument-settlement divergence is predictively consistent with geographic features (ML decomposition; Section 5.11)
- The approximate collinearity is scaffolded by plate tectonics (0.8° from geological pole)
- Easter Island's position on the circle is geographically forced (only inhabited Pacific island within 200 km)
- The sub-km precision is a 1-in-3,000 coincidence within the tectonic scaffold
- The circle tracks the desert belt (93rd percentile for 30th parallel proximity)
- The pre-YD corridor is empty (10 dates, 33rd percentile)
- The date-distance gradient is generic expansion (31/36 bearings equivalent)
- The φ ratio does not survive multiple testing (Bonferroni p = 0.64)

Three questions remain open: (1) the secondary divergence peak at 500 BCE–0 CE (D = 5.98) has no explanation and requires separate investigation; (2) the New World replication, while structurally compelling, is classification-dependent (Section 5.8) and should be confirmed with improved site-level classification; (3) whether the 12 independent high-D clusters identified in the systematic scan correspond to meaningful archaeological regions.

Methodological limitations include: the Pleiades database's Mediterranean bias, the absence of soil productivity and water table depth from the ML feature set, and the lack of radiocarbon coverage in submerged regions. The Bay of Bengal crossing (−106 m) and Persian Gulf crossing (−10 m) represent genuine data gaps.

Pre-agricultural monumental sites (Göbekli Tepe, Karahan Tepe, the Taş Tepeler complex) show no great circle alignment and cluster by ecological zone rather than geographic features, consistent with two distinct modes of monument placement: ecological (pre-agricultural) and geographic (civilizational). The Great Circle alignment is a property of the latter.

---

## 6. Conclusion

The Alison Great Circle alignment is statistically real, replicated across eight independent databases spanning over 550,000 sites and dates. Post-2001 data clusters more strongly than pre-2001 data, confirming the alignment is not a product of post-hoc site selection. Three perpendicular control circles produce maximum |D| = 2.0, confirming the divergence is geometry-specific (Alison D = 8.16). On-corridor sites are not more researched than off-corridor sites; the signal strengthens after controlling for research intensity.

The monument-settlement divergence is predictively consistent with measurable geographic features. A machine learning analysis (XGBoost, 10 features, 21,779 sites) shows that cross-type distance, site clustering, elevation, and coastline distance are sufficient to predict monument vs. settlement placement. Distance to the Great Circle ranks last (#10 of 10) in SHAP feature importance. The divergence is primarily mortuary: 71% of the 2750 BCE spike sites are tombs, pyramids, and necropoleis. The circle that maximizes inverse divergence (settlement enrichment) is separated from the Alison pole by only 10.7°, confirming the monument-settlement ratio is a continuous function of the circle's local crossing angle through archaeological landscapes — the line between where ancient humans buried their dead and where they lived.

The divergence replicates independently in the New World (Peru Ministry D = 6.50, Z = 25.41; XRONOS NW D = 4.69; p3k14c NW D = 4.99), with a different temporal peak (1000-1500 CE) but the same structural pattern: monuments near circle, settlements not. The divergence is not an Old World artifact. Multi-region preservation simulation adding 1,250 estimated missing settlements across four riverine regions produces D > 2 in 100% of iterations under realistic geographic distribution.

The collinearity of Giza, Nazca, and Easter Island is a geologically scaffolded coincidence. Plate tectonics aligns the three geological features responsible for these civilizations within 0.8° of the Alison pole. Easter Island is the only inhabited island within 200 km of the resulting circle in the Pacific. The sub-km precision (rank #1 of 640 equivalent-site combinations) is a 1-in-3,000 rarity within the tectonic framework.

The corridor was empty before the Younger Dryas (10 of 94,181 dates, 33rd percentile). No directional diffusion along the arc. No crop-diffusion alignment. The 36° periodicity does not replicate on data-driven clusters. The corridor did not function as a diffusion channel. The Giza longitude grid is falsified. The golden section distances do not survive multiple testing.

The Alison Great Circle alignment is a real geographic phenomenon with a coherent explanatory framework: plate tectonics scaffolds the approximate alignment; the universal human tendency to place mortuary monuments on dramatic terrain adjacent to settlement zones produces the divergence; a 10° tilt is sufficient to reverse the monument-settlement ratio, confirming the local geographic mechanism; and the sub-km precision of the three anchor sites is a rare but explicable coincidence within the tectonic framework. The alignment does not require a lost civilization, unknown geophysical force, or ancient global coordination. It is predictively consistent with plate tectonics, geography, and one 1-in-3,000 coincidence.

---

## Data Availability

All analysis code: github.com/thegreatcircledata/great-circle-analysis
Interactive visualization: thegreatcircle.earth
Buried settlement site list: Supplementary Table S1
Paleoclimate data: Soreq Cave (Bar-Matthews et al. 2003); GISP2 (Grootes et al. 1993)
Databases: Pleiades (pleiades.stoa.org, CC BY 3.0); p3k14c (github.com/people3k); XRONOS (xronos.ch); Peru Ministry (geoportal.cultura.gob.pe); SAAID (Pandora/PACHAMAMA); CNSA (Kaggle, CC BY 4.0); LuwianSiteAtlas (Scientific Data)

---

## Acknowledgments

The author thanks the Pleiades Gazetteer (NEH); the EES Delta Survey (Spencer 2024); the Megalithic Portal community; Jim Alison for documenting the original great circle hypothesis; J. Mark Kenoyer for insights on Indus Valley site classification; and the open-source communities maintaining the databases used in this study. Computational analysis was performed using Claude Code (Anthropic) and reviewed using Perplexity AI. Anonymous peer review of the analytical approach was provided through iterative AI-assisted critique.

---

---

## Figure Captions

**Figure 1.** Monument-settlement divergence (D) by 500-year bin for Pleiades (top) and p3k14c (bottom) databases. D peaks sharply at 3000–2500 BCE (D = 9.95 for Pleiades, D = 3.83 for p3k14c), corresponding to the period of major Egyptian pyramid construction. The signal is absent before 3500 BCE and collapses by 1500 BCE, demonstrating that the divergence is temporally specific rather than a static recording bias.

**Figure 2.** Distribution of monument-settlement divergence D across 16,200 possible great circle pole positions (2° grid). The Alison pole (red dashed line, D = 7.09) falls in the upper tail of the distribution, above >95% of all tested poles. The global maximum (green dotted line, D = 24.28) corresponds to a circle through central Italy capturing the densest monument region in the Pleiades database.

**Figure 3.** The ten great circles with highest monument-settlement divergence D, identified through exhaustive pole scan. The Alison Great Circle (black dotted line, D = 8.9) is the dominant Old World signal passing through Egypt, South Asia, and South America. Competing high-D circles capture non-overlapping monument clusters (Mediterranean, Southeast Asia), suggesting distinct geographic systems rather than random noise.

**Figure 4.** Preservation bias Monte Carlo simulation. Left: D distribution under baseline conditions (blue) vs. with +350 estimated missing Egyptian settlements (tan). Middle: multi-region augmentation (+1,250 settlements across Egypt, Mesopotamia, Indus, Peru) showing D remains above 2 in 100% of iterations. Right: settlement Z-score distribution under augmentation. The divergence is robust to aggressive preservation correction.

**Figure 5.** SHAP (SHapley Additive exPlanations) feature importance for the global XGBoost classifier predicting monument vs. settlement classification. Distance to the Great Circle ranks last (#10 of 10 features), with cross-type distance, site density, and elevation dominating. Removing Great Circle distance changes AUC by +0.0002, indicating the circle adds no independent predictive information beyond known geographic features.

**Figure 6.** Monument (top, red) and settlement (bottom, green) counts along the Alison Great Circle by arc position (degrees from Giza, within 50 km). Monuments (N = 3,497) show pronounced clustering at specific arc positions; settlements (N = 184) show a fundamentally different distribution pattern, concentrated in the first 60° (Near East). The dashed lines indicate the expected count under uniform distribution.

---

## Supplementary Figures

**Figure S1.** Sensitivity surface showing divergence D as a function of pole position perturbation (±30° latitude and longitude from the Alison pole). The Alison pole (black star, D = 9.52) sits within a broad elevated region, indicating the result is robust to small pole perturbations (±2–3°). The maximum D = 32.83 occurs at a distant pole position capturing a different monument cluster.

**Figure S2.** Monument-settlement divergence at 250-year resolution (Pleiades, top; p3k14c, bottom). The signal narrows to a single 250-year bin (2750–2500 BCE) in Pleiades, where 31 of 38 on-circle monuments are Memphis necropolis pyramids. The higher temporal resolution confirms the spike is not a broad-epoch artefact.

**Figure S3.** Egyptian pyramid construction date vs. distance from the Alison Great Circle (n = 26 pyramids, Djoser through Amenemhat III). Earlier pyramids (Dynasty 4, blue) cluster near the circle (3–10 km); later pyramids (Dynasty 5–12) drift further. The negative correlation (r = −0.428, p = 0.029) is consistent with the Memphis necropolis expanding outward from a line near the great circle, but is also consistent with geographic expansion of pyramid-building activity along the Nile.

**Figure S4.** Monument orientation offset from the local Great Circle bearing. Near-circle monuments (<100 km, left, n = 56) show modest clustering around the GC bearing direction; off-circle monuments (≥100 km, right, n = 94) show a more uniform distribution. The difference is suggestive but not statistically significant after multiple-testing correction.

**Figure S5.** Pre-Younger Dryas corridor test. Top: enrichment ratio of on-corridor radiocarbon dates relative to random circles over time (94,181 merged dates). Bottom: raw date counts on the corridor. The corridor shows no enrichment before 12,800 BP, ruling out a deep-time origin for the alignment. Post-YD enrichment rises steadily, consistent with the monument-building era.

**Figure S6.** Multi-proxy climate comparison with monument-settlement divergence (500-year bins). Divergence D (red) peaks at 3000–2500 BCE while paleoclimate proxies (Soreq Cave δ¹⁸O, Dead Sea levels, marine core δ¹⁸O, Jinmium Cave δ¹⁸O) show no corresponding pattern. The divergence does not correlate with any tested climate proxy at any resolution (Section 4.8).

## References

Allan, E. (2026). Statistical analysis of ancient monumental site distribution along a proposed great circle. Under review, PLOS ONE (PONE-D-26-13610). Preprint: doi.org/10.5281/zenodo.19046176.

Alison, J. (c. 2001). The prehistoric alignment of world wonders.

Bar-Matthews, M., Ayalon, A., Gilmour, M., Matthews, A., & Hawkesworth, C. J. (2003). Sea–land oxygen isotopic relationships from planktonic foraminifera and speleothems in the Eastern Mediterranean region. Geochimica et Cosmochimica Acta, 67(17), 3181–3199.

Bietak, M. (1996). Avaris: The Capital of the Hyksos. British Museum Press.

deMenocal, P. B. (2001). Cultural responses to climate change during the late Holocene. Science, 292(5517), 667–673.

Grant, K. M., Rohling, E. J., Bar-Matthews, M., Ayalon, A., Medina-Elizalde, M., Ramsey, C. B., Satow, C., & Roberts, A. P. (2012). Rapid coupling between ice volume and polar temperature over the past 150,000 years. Nature, 491, 744–747.

Grootes, P. M., Stuiver, M., White, J. W. C., Johnsen, S., & Jouzel, J. (1993). Comparison of oxygen isotope records from the GISP2 and GRIP Greenland ice cores. Nature, 366, 552–554.

Kemp, B. J. (2018). Ancient Egypt: Anatomy of a Civilization (3rd ed.). Routledge.

Klein Goldewijk, K., Beusen, A., Doelman, J., & Stehfest, E. (2017). Anthropogenic land use estimates for the Holocene — HYDE 3.2. Earth System Science Data, 9, 927–953.

Martindale, A., et al. (2016). Canadian Archaeological Radiocarbon Database (CARD 2.1).

Moeller, N. (2016). The Archaeology of Urbanism in Ancient Egypt. Cambridge University Press.

Palmisano, A., et al. (2025). XRONOS: A global open-access radiocarbon database.

Parcak, S. (2016). TED Prize talk: Help discover ancient ruins before it's too late.

Spencer, J. (2024). The Sites of the Delta: A Gazetteer. EES Excavation Memoir 119.

Staubwasser, M., Sirocko, F., Grootes, P. M., & Segl, M. (2003). Climate change at the 4.2 ka BP termination of the Indus valley civilization. Geophysical Research Letters, 30(8).

Weiss, H., Courty, M. A., Wetterstrom, W., et al. (1993). The genesis and collapse of third millennium North Mesopotamian civilization. Science, 261(5124), 995–1004.

Wernke, S., VanValkenburgh, P., & Arkush, E. (2024). Large-scale, collaborative imagery survey in archaeology: the Geospatial Platform for Andean Culture, History and Archaeology (GeoPACHA). Antiquity, 98(397), 155–171.

Sheisha, H., Kaniewski, D., Marriner, N., et al. (2022). Nile waterscapes facilitated the construction of the Giza pyramids during the 3rd millennium BCE. Proceedings of the National Academy of Sciences, 119(37), e2202530119.

Borrell, F., Junno, A., & Barceló, J. A. (2015). Synchronous environmental and cultural change in the emergence of agricultural economies 10,000 years ago in the Levant. PLoS ONE, 10(8), e0134810.

Bird, D., et al. (2022). p3k14c, a synthetic global database of archaeological radiocarbon dates. Scientific Data, 9, 27.

Hancock, G. (1995). Fingerprints of the Gods. Crown Publishers.

Hancock, G. (1998). Heaven's Mirror: Quest for the Lost Civilization. Crown Publishers.

Vermeersch, P. M. (2024). Radiocarbon Palaeolithic Europe Database, Version 32. Available at https://ees.kuleuven.be/en/geography/projects/14c-palaeolithic.

Palmisano, A., Bevan, A., & Shennan, S. (2017). Comparing archaeological proxies for long-term population patterns: An example from central Italy. Journal of Archaeological Science, 87, 59–72.

Rick, J. W. (1987). Dates as data: An examination of the Peruvian preceramic radiocarbon record. American Antiquity, 52(1), 55–73.

Dee, M. W., et al. (2013). An absolute chronology for early Egypt using radiocarbon dating and Bayesian statistical modelling. Proceedings of the Royal Society A, 469(2159), 20130395.

Surovell, T. A., et al. (2009). Correcting temporal frequency distributions for taphonomic bias. Journal of Archaeological Science, 36(8), 1715–1724.

Lundberg, S. M., & Lee, S.-I. (2017). A unified approach to interpreting model predictions. Advances in Neural Information Processing Systems, 30, 4765–4774.

Molnar, C. (2020). Interpretable Machine Learning: A Guide for Making Black Box Models Explainable (2nd ed.). https://christophm.github.io/interpretable-ml-book/

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. Journal of the Royal Statistical Society B, 57(1), 289–300.

Chen, T., & Guestrin, C. (2016). XGBoost: A scalable tree boosting system. Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, 785–794.

Lehner, M. (1997). The Complete Pyramids: Solving the Ancient Mysteries. Thames & Hudson.

Stadelmann, R. (2001). Die großen Pyramiden von Giza. Akademische Druck- und Verlagsanstalt.

DeMets, C., Gordon, R. G., & Argus, D. F. (2010). Geologically current plate motions. Geophysical Journal International, 181(1), 1–80.

Steinberger, B., Sutherland, R., & O'Connell, R. J. (2004). Prediction of Emperor-Hawaii seamount locations from a revised model of global plate motion and mantle flow. Nature, 430(6996), 167–173.

Cahill, T., & Isacks, B. L. (1992). Seismicity and shape of the subducted Nazca Plate. Journal of Geophysical Research, 97(B12), 17503–17529.

Chorowicz, J. (2005). The East African rift system. Journal of African Earth Sciences, 43(1–3), 379–410.

Bevan, A. (2012). Spatial methods for analysing large-scale artefact inventories. Antiquity, 86(332), 492–506.

Crema, E. R. (2022). Statistical inference of prehistoric demography from frequency distributions of radiocarbon dates: A review and a guide for the perplexed. Journal of Archaeological Method and Theory, 29, 1387–1418.

Ripley, B. D. (1976). The second-order analysis of stationary point processes. Journal of Applied Probability, 13(2), 255–266.
