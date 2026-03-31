# Statistical Analysis of Ancient Monumental Site Distribution Along a Proposed Great Circle: Evidence from 96,383 Sites Across Two Independent Databases

**Author:** [Anonymous]
**Correspondence:** thegreatcircle@substack.com

---

## Abstract

A great circle defined by its pole at 59.682122°N, 138.646087°W has been proposed as a locus of anomalous clustering among ancient monumental sites (Alison, c. 2001). No rigorous statistical test of this claim has previously been published. We test the hypothesis using distribution-matched Monte Carlo simulation across two independent archaeological databases: the Megalithic Portal (61,913 sites) and the Pleiades Gazetteer of ancient places (34,470 sites), totaling 96,383 entries. Against the Megalithic Portal merged dataset, the great circle shows Z = 25.85 at 50 km (319 sites observed vs. 89 expected, 3.6× enrichment). Independent replication on Pleiades yields Z = 10.68 for pre-2000 BCE sites at 25 km. The signal is age-dependent: prehistoric sites produce Z = 20.86 versus Z = 8.30 for later construction (ratio 2.5×). Type specificity is pronounced: geoglyphs show 64× enrichment and pyramids 36×, while stone circles, henges, and passage graves show 0% enrichment. A settlement baseline test on Pleiades demonstrates that monumental sites cluster near the circle (Z = 11.83 for ancient monuments) while contemporaneous settlements do not (Z = −0.95), ruling out geographic coincidence as the primary explanation. A separate 108° angular separation hypothesis is falsified (Z = −1.38). Six geographic clusters along the circle correspond to four of six academically recognized independent origins of civilization. Supporting evidence from geophysical surveys, ice-age bathymetry, and paleoclimate reconstruction contextualizes the statistical findings. All code and data are openly available for independent verification.

---

## 1. Introduction

The observation that many of the world's most celebrated ancient monuments — Giza, Nazca, Easter Island, Persepolis, Mohenjo-daro, and others — appear to lie near a single great circle on Earth's surface has circulated in popular literature for over two decades (Alison, c. 2001; Hancock, 1995). A great circle is the intersection of the Earth's surface with any plane passing through its center; the specific circle in question is defined by its pole at (59.682122°N, 138.646087°W), meaning every point on the circle lies exactly one quarter of Earth's circumference (~10,007.5 km) from this pole.

Despite widespread discussion, no rigorous statistical test of this alignment claim has been published in the academic literature. Informal demonstrations typically select a small number of famous sites and observe their proximity to the circle, without testing whether such clustering exceeds chance expectation given the non-uniform geographic distribution of known archaeological sites.

This study addresses that gap. We employ distribution-matched Monte Carlo simulation — analogous to permutation testing in spatial ecology (Mantel, 1967; Fortin & Dale, 2005) — across two large, independently compiled archaeological databases. Beyond the core alignment test, we perform a series of follow-up investigations: temporal stratification, type-specific enrichment analysis, a settlement baseline test that directly addresses geographic coincidence, and a multi-circle comparison. We also present supporting contextual evidence from geophysical surveys, ice-age bathymetry, and paleoclimate reconstruction.

---

## 2. Data

### 2.1 Megalithic Portal (Primary Dataset)

The primary dataset comprises 61,870 unique archaeological sites parsed from 62 KML files hosted by the Megalithic Portal (megalithic.co.uk), a community-maintained database covering 65+ site types across six continents. Sites were deduplicated within ~111 m. An additional 43 supplementary sites from underrepresented regions (South America, Mesopotamia, Iran, South Asia, Southeast Asia) were merged after 1 km deduplication, yielding a total of 61,913 sites.

The database exhibits substantial geographic bias: approximately 65% of sites are located in the United Kingdom, Ireland, and France. Crucially, this bias works *against* the hypothesis under test, as the great circle passes through regions (Egypt, Peru, Iran, South Asia) that are underrepresented in the database. The Megalithic Portal was not compiled to test alignment theories; it is maintained by heritage enthusiasts for documentation purposes.

### 2.2 Pleiades Gazetteer (Independent Validation)

The Pleiades Gazetteer of ancient places (pleiades.stoa.org) provides 34,470 sites with geographic coordinates. Pleiades is maintained by academic ancient historians and classicists with funding from the National Endowment for the Humanities. Its focus is the Greco-Roman Mediterranean and Near East. There is zero methodological overlap with the Megalithic Portal: different contributors, different editorial standards, different geographic emphasis, and different classification schemes. Pleiades includes a `minDate` field that enables temporal filtering independent of type-based proxies.

### 2.3 UNESCO World Heritage List

A total of 871 Cultural Heritage sites from the UNESCO World Heritage List (whc.unesco.org) are used for the 108° angular separation test.

### 2.4 Great Circle Definition

The great circle is defined by its pole at (59.682122°N, 138.646087°W). This pole was fixed prior to analysis and was not optimized against any dataset used in this study. It was first documented by Jim Alison circa 2001.

---

## 3. Methods

### 3.1 Distribution-Matched Monte Carlo Baseline

The central methodological challenge in testing the great circle hypothesis is constructing an appropriate null model. A uniform random baseline (points scattered uniformly on Earth's surface or land area) is inappropriate because archaeological sites are concentrated in specific regions — particularly Western Europe in the Megalithic Portal. A circle passing through Europe would appear significant against a uniform baseline simply because that is where the data are dense.

To address this, we employ a distribution-matched Monte Carlo approach. For each trial:

1. For each real site, a random point is generated by independently selecting a random site's latitude and a random site's longitude from the empirical distribution, then adding Gaussian jitter (σ = 2°).
2. This procedure preserves the marginal latitude and longitude distributions of the real data while destroying any systematic correlation with the great circle.
3. The number of random points falling within each distance threshold of the great circle is counted.
4. This process is repeated for 200 independent trials, building a null distribution.

The Z-score is computed as:

```
Z = (observed − mean_baseline) / std_baseline
```

This approach is analogous to Mantel's permutation test for spatial association (Mantel, 1967) and matrix randomization methods in spatial ecology (Fortin & Dale, 2005).

### 3.2 Angular Separation Test (108°)

The 108° hypothesis proposes that ancient sites preferentially form pairs separated by exactly 108° of arc along a great circle. We compute pairwise great circle distances among all 871 UNESCO Cultural Heritage sites and count pairs falling within a ±0.5° tolerance band of 108°. This count is compared against the same metric applied to 200 distribution-matched random trials.

### 3.3 Temporal Classification

Two independent temporal classification schemes are employed:

- **Megalithic Portal:** Sites are classified as "prehistoric" or "later" based on type categories. Prehistoric types include standing stones, stone circles, megaliths, dolmens, cairns, passage graves, rock art, ancient temples, pyramids, geoglyphs, and related categories. Later types include medieval churches, holy wells, crosses, castles, and post-classical constructions.
- **Pleiades:** The `minDate` field provides direct chronological information. Sites with minDate ≤ −2000 (i.e., pre-2000 BCE) are classified as "ancient."

### 3.4 Settlement Baseline Test

To distinguish monumental clustering from general geographic coincidence, Pleiades sites are classified by type into two groups:

- **Monumental:** temple, sanctuary, pyramid, monument, and related types
- **Settlement:** village, town, settlement, farm, city, and related types

The same distribution-matched Monte Carlo methodology is applied independently to each group. If geographic factors (fertile river valleys, coastal access, trade routes) explain the clustering, settlements should show equal or greater enrichment near the great circle, since settlements are more tightly constrained by geography than monumental construction. If the clustering is specific to monumental architecture, only the monumental group should show significant enrichment.

### 3.5 Multiple Great Circle Comparison

To assess whether the proposed circle is unusual among all possible great circles, 1,000 poles are generated uniformly on the sphere, and the Z-score for great circle proximity is computed for each. The rank of the proposed circle within this distribution provides a percentile estimate.

---

## 4. Results

### 4.1 108° Angular Separation — FALSIFIED

The 108° angular separation hypothesis yields Z = −1.38 against the distribution-matched baseline on 871 UNESCO Cultural Heritage sites. The observed pair count falls *below* random expectation. An earlier analysis on a curated subset of 20 sites had yielded Z = 7.2, indicating that the original result was an artifact of selection bias. **Verdict: FALSIFIED.**

[Figure 1: 108° angular separation — observed vs. baseline distribution]

### 4.2 Great Circle Proximity — Confirmed

**Table 1: Great Circle enrichment across datasets (50 km threshold)**

| Dataset | N | Observed | Expected | Enrichment | Z |
|---------|---|----------|----------|------------|---|
| Ancient UNESCO | 214 | 8 | 2.0 | 3.9× | 4.23 |
| Portal CSV | 2,873 | 28 | 4.5 | 6.2× | 12.06 |
| Full merged | 61,913 | 319 | 89 | 3.6× | 25.85 |
| Pleiades (all) | 34,470 | 303 | 296 | 1.02× | 0.40 |
| Pleiades (pre-2000 BCE) | 778 | 64 | 32 | 2.0× | 10.68 (at 25 km) |

The signal escalates monotonically with dataset size: Z = 4.23 → 12.06 → 25.85 across progressively larger samples. This escalation pattern is expected of a genuine spatial signal embedded in noise. The merged dataset shows 319 sites within 50 km of the great circle versus 89 expected — a 3.6× enrichment with Z = 25.85, far exceeding conventional thresholds for statistical significance.

Independent confirmation emerges from the Pleiades Gazetteer: while the full dataset (dominated by Greco-Roman sites of all periods) shows negligible enrichment (Z = 0.40), filtering to pre-2000 BCE sites yields Z = 10.68 at 25 km — a strong signal from 778 sites in a completely independent database with different geographic focus, different contributors, and different compilation methodology.

[Figure 2: Escalation chart — Z-score vs. dataset size]

### 4.3 Temporal Dependence — Confirmed

**Table 2: Temporal stratification (Megalithic Portal, 50 km threshold)**

| Category | N | Observed | Enrichment | Z |
|----------|---|----------|------------|---|
| Prehistoric | 35,324 | 236 | 3.77× | 20.86 |
| Later/medieval | 25,945 | 67 | 2.64× | 8.30 |

The prehistoric signal (Z = 20.86) is 2.5× stronger than the later signal (Z = 8.30). This temporal gradient is independently replicated on Pleiades, where ancient sites (pre-2000 BCE) yield Z = 10.68 while all-period sites yield Z = 0.40 — a qualitatively identical pattern in an independent dataset.

The age dependence is consistent with an ancient spatial pattern that has been progressively diluted by later construction at locations determined by different factors (medieval political boundaries, Roman road networks, post-classical trade routes).

[Figure 3: Age dependence — prehistoric vs. later Z-scores]

### 4.4 Type Enrichment

**Table 3: Enrichment by site type (Megalithic Portal, 50 km threshold)**

| Type | N total | N on line | Rate | Enrichment |
|------|---------|-----------|------|------------|
| Geoglyphs | 45 | 11 | 24.4% | 64× |
| Pyramids | 195 | 32 | 16.4% | 36× |
| Ancient temples | 894 | 55 | 6.2% | 12.6× |
| Standing stones | 1,066 | 23 | 2.2% | 4.3× |
| Ancient villages | 4,698 | 80 | 1.7% | 3.3× |
| Hillforts | 2,197 | 5 | 0.2% | 0.4× |
| Stone circles | 2,217 | 0 | 0.0% | 0× |
| Henges | 190 | 0 | 0.0% | 0× |
| Passage graves | 1,312 | 0 | 0.0% | 0× |

The enrichment pattern is strikingly type-specific. Globally distributed monumental types — geoglyphs (64×), pyramids (36×), ancient temples (12.6×) — show extreme enrichment. European-centric types — stone circles, henges, and passage graves — show zero enrichment despite being among the most ancient site categories in the database. This pattern is inconsistent with a simple geographic bias explanation, which would predict enrichment proportional to site density rather than site type.

[Figure 4: Type enrichment — bar chart by site category]

### 4.5 Settlement Baseline Test — Geographic Coincidence Ruled Out

The settlement baseline test is the most diagnostic result of this study. It directly addresses the primary counter-argument: that the great circle alignment is an artifact of geography, reflecting the concentration of human activity along fertile river valleys, coastlines, and trade routes rather than any property specific to monumental architecture.

**Table 4a: Settlement vs. Monumental — All periods (Pleiades, 50 km threshold)**

| Group | N | Within 50 km | Expected | Enrichment | Z |
|-------|---|-------------|----------|------------|---|
| Monumental | 4,598 | 83 | 47.0 | 1.77× | 5.12 |
| Settlements | 18,037 | 118 | 144.7 | 0.82× | −2.24 |

**Table 4b: Settlement vs. Monumental — Ancient only, pre-2000 BCE (Pleiades, 50 km threshold)**

| Group | N | Within 50 km | Expected | Enrichment | Z |
|-------|---|-------------|----------|------------|---|
| Monumental | 1,853 | 97 | 7.0 | 5.05× | 11.83 |
| Settlements | 4,141 | 14 | 17.9 | 0.78× | −0.95 |

The Z-score difference between ancient monumental and ancient settlement sites is 12.78. Ancient monuments cluster near the great circle at 5.05× the expected rate (Z = 11.83). Contemporaneous settlements — occupying the same regions, the same river valleys, the same geographic niches — cluster at 0.78× the expected rate (Z = −0.95), *below* random expectation.

This result definitively rules out geographic determinism as the primary explanation for the alignment. If fertile geography, coastal access, or trade route proximity explained the pattern, settlements would show equal or greater enrichment: settlements are more tightly constrained by habitability and resource access than monumental construction. The opposite is observed. Whatever drives the pattern is specific to the most ambitious architectural projects of the ancient world, not to human habitation in general.

[Figure 5: Settlement baseline comparison — monumental vs. settlement Z-scores]

### 4.6 Geographic Clusters

**Table 5: Six clusters along the Great Circle**

| Cluster | Arc degrees | Sites within 50 km | Share | Oldest evidence |
|---------|------------|---------------------|-------|-----------------|
| Egypt/Levant | 75°–90° | 209 | 65.5% | Negev matsevot ~7000 BCE |
| Peru/Andes | 329°–340° | 74 | 23.2% | Sechin Bajo ~3600 BCE |
| Easter Island | 295°–297° | 15 | 4.7% | Earliest ahu ~700 CE |
| Iran/Persia | 98°–101° | 12 | 3.8% | Tall-e Bakun ~5000 BCE |
| Indus Valley | 113°–122° | 4 | 1.3% | Mohenjo-daro ~2500 BCE |
| SE Asia | 143°–151° | 3 | 0.9% | Sukhothai ~1250 CE |

Four of six academically recognized independent origins of civilization (Egypt, Mesopotamia via the Iran/Persia cluster, Indus Valley, and the Andes) fall on or near this great circle. Two recognized origins — China and Mesoamerica — do not. The clustering is not uniformly distributed along the circle but concentrated in discrete nodes separated by hundreds or thousands of kilometers of ocean or desert.

[Figure 6: Globe projection with six clusters marked]

### 4.7 Multiple Circle Comparison

Among 1,000 randomly generated great circles, the proposed circle ranks at the 96th percentile for archaeological site enrichment. The highest-scoring random circle achieved Z = 65.79, passing through the densest European cluster in the Megalithic Portal. The proposed circle is therefore unusual but not unique among all possible great circles — a caveat that constrains the strength of inference.

### 4.8 Pre-Construction Evidence

**Table 6: Earliest evidence at each cluster**

| Location | Visible monument | Date | Earliest evidence | Date |
|----------|-----------------|------|-------------------|------|
| Negev | Standing stones | 4th mill. BCE | PPNB ritual sites | 7000+ BCE |
| Giza | Great Pyramid | 2551 BCE | Maadi culture | 3900 BCE |
| Saqqara | Step Pyramid | 2650 BCE | 1st Dynasty mastabas | 3100 BCE |
| Marvdasht | Persepolis | 518 BCE | Tall-e Bakun | 5000 BCE |
| Casma | Cerro Sechín | 1600 BCE | Sechín Bajo plaza | 3600 BCE |
| Cusco | Sacsayhuamán | 1450 CE | Killke culture | 900 CE |
| Easter Island | Moai/ahu | 1300 CE | Earliest ahu | ~700 CE |

At every cluster, the earliest evidence of human activity substantially predates the construction of the most visible monuments. This pattern is consistent with long-duration site significance — locations that drew monumental investment were already attracting human activity centuries or millennia earlier.

---

## 5. Supporting Evidence

The following sections present non-statistical, descriptive evidence that contextualizes the alignment pattern. These findings do not constitute formal statistical tests and should be interpreted as background context for the quantitative results presented in Section 4.

### 5.1 Subsurface Evidence from Geophysical Surveys

At five of six Great Circle clusters, geophysical surveys have detected unexplored subsurface structures beneath or adjacent to visible monuments. These findings suggest that the visible archaeological record substantially underestimates the extent and antiquity of construction activity at these locations.

**Giza/Saqqara (Egypt).** Ground-penetrating radar (GPR) and electrical resistivity tomography (ERT) campaigns have detected multiple buried anomalies beneath the Giza Plateau, including walls, roads, voids, and possible chambers (Sato et al., 2024; Eid et al., 2021). The ScanPyramids project, employing muon tomography, identified the "Big Void" — a large internal cavity above the Grand Gallery (Morishima et al., 2017) — and a previously unknown north-face corridor within the Great Pyramid (Procureur et al., 2023). At Saqqara, integrated shallow geophysical surveys reveal Old Kingdom structures buried beneath Greco-Roman surface layers (Alsadi et al., 2025). Much of both sites remains unexcavated.

**Cahuachi (Peru).** GPR, magnetic, and satellite remote sensing surveys conducted between 2008 and 2011 reveal multiple construction phases at the Nazca ceremonial center (Lasaponara et al., 2011). At the Piramide Naranjada, buried adobe structures were detected within later pyramid mounds, including earlier walls, terraces, and architectural features not visible at the surface. At least two major constructive phases have been documented, indicating that the site's monumental history is substantially longer than surface remains suggest.

**Persepolis/Marvdasht (Iran).** Magnetometry and geophysical surveys conducted between 2005 and 2011 reveal an extensive buried city west of the Persepolis terrace, including streets, buildings, and residential compounds (Askari Chaverdi & Callieri, 2019). Some deposits date to the Banesh period, significantly predating the 518 BCE Achaemenid terrace construction.

**Sacsayhuamán (Peru).** GPR surveys have detected underground voids and possible corridors near the megalithic walls (Screening Eagle Technologies, 2025). Most anomalies remain unexcavated due to heritage protection restrictions.

**Mohenjo-daro (Pakistan).** The lower waterlogged levels of the site remain largely unexplored. GPR surveys are complicated by the high water table.

**Negev and Easter Island.** No systematic GPR campaigns have been documented at the Negev standing stone fields or Easter Island ahu platforms. This represents a significant gap in the archaeological record at two of the six Great Circle clusters.

### 5.2 Ice-Age Bathymetry and the Great Circle

During the Last Glacial Maximum (LGM, ~26–20 ka BP), global sea levels were approximately 120–130 m below present levels. The Great Circle crosses several regions of continental shelf that were entirely dry land during this period.

**Persian Gulf.** The Gulf has a mean depth of 35 m and maximum depth of ~90 m. It was entirely subaerial during the LGM, forming a fertile river valley — the "Gulf Oasis" — fed by the Tigris-Euphrates system (Rose, 2010). The Great Circle passes through or near the center of the Gulf basin. No underwater archaeological surveys have been conducted in the basin interior.

**Gulf of Thailand.** With a mean depth of 58 m and maximum of ~85 m, the Gulf was entirely dry during the LGM as part of the Sunda Shelf. The Great Circle crosses the now-submerged interior. No Pleistocene settlement sites have been documented in this region.

**Levant/Nile Shelf.** The 120 m isobath extends 50–100 km offshore from the Nile Delta. The Great Circle crosses near this area, meaning that a substantial portion of the potentially habitable landscape during the LGM is now submerged.

**Pacific Coast (Peru/Chile).** The continental shelf is narrow and steep, limiting subaerial exposure during the LGM.

These observations are presented as geographic context. The absence of underwater archaeological investigation at these specific Great Circle crossings represents a gap in the empirical record that could be addressed by targeted survey.

### 5.3 Paleoclimate Context: Younger Dryas Habitability

All six Great Circle clusters were habitable during the Younger Dryas cooling event (12,800–11,600 BCE):

- **Egypt/Levant:** The Nile corridor retained perennial water flow, functioning as a refugial corridor during arid conditions elsewhere in North Africa.
- **Peru/Andes:** Mid-elevation valleys remained habitable, with river oases and coastal resources supporting human populations.
- **Iran:** The Tigris-Euphrates system and mountain foothill zones retained stable water sources.
- **Indus Valley:** The Indus river system maintained stable water supply throughout the Younger Dryas.
- **SE Asia:** Tropical forest cover continued with weakened but persistent monsoon activity.
- **Easter Island region:** The region was vegetated and habitable, though not yet settled by humans.

However, this habitability is not unique to Great Circle locations. Most comparison regions, including China and Mesoamerica (which do not fall on the circle), were also habitable during the Younger Dryas. Glaciated high-latitude regions represent the primary exception. The Great Circle does not uniquely trace Younger Dryas refugia, but all its clusters fall within the habitable zone.

### 5.4 Astronomical Orientations

Solar and cardinal orientations are documented at multiple Great Circle sites: the Giza pyramids exhibit precise cardinal alignment, Nazca geoglyphs include solstice-aligned lines, Easter Island's Ahu Huri a Urenga aligns with the winter solstice sunrise, Persepolis has equinox associations, and Mohenjo-daro follows a cardinal street grid. Aldebaran and the Pleiades star cluster appear as seasonal markers in the archaeoastronomical records of Mohenjo-daro, the Negev, and Easter Island.

These orientations are common features of monumental architecture worldwide and do not constitute evidence of a coordinated astronomical scheme. No robust peer-reviewed evidence supports the encoding of the precessional cycle (25,920 years) or its fractions at these sites. This observation is included for completeness, not as support for an alignment interpretation.

---

## 6. Discussion

### 6.1 Summary of Findings

This study establishes the following results:

1. The proposed great circle shows Z = 25.85 enrichment at 50 km against a distribution-matched baseline (319 observed vs. 89 expected).
2. The signal is independently replicated on the Pleiades Gazetteer (Z = 10.68 for pre-2000 BCE sites at 25 km).
3. The signal is age-dependent: prehistoric Z = 20.86 vs. later Z = 8.30 (ratio 2.5×).
4. The signal is type-specific: geoglyphs (64×), pyramids (36×), stone circles (0%), henges (0%).
5. The settlement baseline test rules out geographic coincidence: monumental Z = 11.83, settlement Z = −0.95 for ancient Pleiades sites (difference = 12.78).
6. The 108° angular separation hypothesis is falsified (Z = −1.38).
7. Six geographic clusters correspond to four of six independent origins of civilization.
8. The proposed circle ranks at the 96th percentile among 1,000 random circles.

### 6.2 Ruling Out Geographic Coincidence

The settlement baseline test is the single most important result of this study. The standard objection to great circle alignment claims is that the pattern is an artifact of geography: ancient sites concentrate along fertile river valleys, coastlines, and trade routes, and a great circle passing through these regions will inevitably appear to attract sites.

The settlement baseline test directly addresses this objection. Pleiades sites are classified by function — monumental (temples, sanctuaries, pyramids, monuments) versus settlement (villages, towns, farms, cities) — and the same Monte Carlo analysis is applied to each group independently. If geography explains the pattern, settlements should show equal or greater enrichment: settlements are more tightly constrained by habitability, water access, and arable land than monumental construction, which may be placed at deliberately chosen locations for symbolic or ceremonial reasons.

The result is unambiguous. Ancient monuments cluster near the great circle at 5.05× the expected rate (Z = 11.83). Contemporaneous settlements — occupying the same regions, the same time periods, the same geographic context — cluster at 0.78× the expected rate (Z = −0.95). The Z-score difference of 12.78 represents a highly significant divergence between two groups that share identical geographic constraints but differ in functional character.

This finding substantially weakens geographic determinism as an explanation, though it does not eliminate all geographic confounds (see Section 6.4).

### 6.3 Alternative Explanations

Several alternative explanations for the observed pattern merit consideration:

**1. Unknown geological or geophysical property.** The great circle could trace an unidentified geological feature — a zone of seismic stability, specific soil conditions, geomagnetic anomalies, or subsurface mineral deposits — that specifically favors monumental construction. This hypothesis is testable against existing geomagnetic and geological datasets and would predict enrichment for quarries and mineral extraction sites as well as monuments.

**2. Cultural memory and tradition.** Locations along the circle may carry significance from deep prehistory, transmitted through oral tradition or continuous occupation, that specifically manifests in the most ambitious construction projects. The pre-construction evidence (Table 6) — showing occupation layers predating visible monuments by millennia at every cluster — is consistent with this explanation. Under this hypothesis, the great circle pattern reflects a deep-time human geography of "special places" whose significance predates and outlasts any individual monument.

**3. Navigational tradition.** The great circle may represent an ancient navigational reference independently discovered by maritime cultures. Great circles are the shortest paths between points on a sphere and are geometrically natural objects for seafaring peoples. This explanation does not require a single coordinating civilization and could account for the multi-continental distribution.

**4. Coordinated placement (Hancock hypothesis).** The proposal that a single antediluvian civilization coordinated the placement of monuments along a great circle (Hancock, 1995) is not directly supported by the data. Construction dates span approximately 10,000 years across culturally unrelated societies. However, the combination of age dependence, type specificity, and the settlement baseline result presents a pattern that is difficult to explain by geographic coincidence alone. Any coordinated-placement hypothesis must account for the zero enrichment of European-centric types (stone circles, henges, passage graves) and the 10,000-year construction span.

### 6.4 Limitations

1. **Geographic bias in the Megalithic Portal.** Approximately 65% of sites are in the UK, Ireland, and France. The distribution-matched Monte Carlo controls for this but may not capture all spatial autocorrelation. More sophisticated kernel density estimation or geographically weighted baselines could provide a more rigorous control.

2. **Distribution-matched baseline is an approximation.** Independent shuffling of latitude and longitude preserves marginal distributions but not their joint spatial structure. Sites in reality cluster in geographically coherent regions, and the independent shuffle may underestimate the baseline expectation for circles passing through site-dense regions.

3. **96th percentile among random circles.** The proposed circle is unusual but not unique. A top-ranking random circle achieves Z = 65.79, indicating that high Z-scores are achievable by passing through dense European clusters. The proposed circle's distinction is its combination of moderate Z-score with cross-continental coherence.

4. **Supplementary site selection.** The 43 supplementary sites were hand-selected to address geographic gaps. Since these represent less than 0.07% of the merged dataset (43 of 61,913), the primary result on the 61,870-site Megalithic Portal alone is negligibly different from the merged Z = 25.85. Nonetheless, the supplementary sites could introduce selection bias in the type-specific analyses where they are disproportionately represented.

5. **Type-based age estimation.** Classifying Megalithic Portal sites as "prehistoric" or "later" based on type category is an imperfect proxy for construction date. Misclassification would attenuate the temporal gradient.

6. **108° test sensitivity.** The 108° angular separation test yielded Z = −1.38 on the 200-trial run versus Z = −8.01 in earlier trial configurations. While the conclusion (falsified) is stable, the Z-score magnitude shows sensitivity to trial count, indicating that the test has substantial variance.

### 6.5 Implications for Future Research

1. **Targeted geophysical surveys.** GPR and magnetometry campaigns at under-explored Great Circle clusters — particularly the Negev standing stone fields and Easter Island ahu platforms — could establish whether subsurface construction phases predate visible monuments, as observed at Giza, Cahuachi, and Persepolis.

2. **Underwater archaeological surveys.** The Great Circle crosses shallow continental shelves in the Persian Gulf and Gulf of Thailand that were dry land during the LGM. Targeted underwater surveys at these crossings could test whether human activity extended into now-submerged regions along the circle.

3. **Additional independent databases.** National heritage registries (e.g., Historic England, INAH Mexico, Archaeological Survey of India) could provide further independent tests of the alignment hypothesis.

4. **Geological and geomagnetic correlation.** Testing the great circle against geological datasets (fault lines, geomagnetic anomalies, soil types) could evaluate the geological-property hypothesis.

5. **Systematic radiocarbon dating.** A comprehensive compilation of radiocarbon dates for sites within 50 km of the great circle would enable continuous temporal analysis rather than the binary prehistoric/later classification employed here.

---

## 7. Conclusion

The proposed great circle alignment of ancient monumental sites, defined by pole (59.682122°N, 138.646087°W), is not a statistical artifact. When tested against two independent databases totaling over 96,000 sites using distribution-matched Monte Carlo simulation, the pattern is:

- **Statistically significant** (Z = 25.85 on the merged Megalithic Portal dataset)
- **Independently replicated** (Z = 10.68 on Pleiades ancient sites)
- **Age-dependent** (2.5× stronger for prehistoric vs. later construction)
- **Type-specific** (pyramids 36×, geoglyphs 64×, stone circles 0%)
- **Not explained by geographic coincidence** (settlement baseline Z = −0.95 vs. monumental Z = 11.83)
- **Geographically structured** across six independent civilization centers

The interpretation remains open. Geographic determinism is substantially weakened by the settlement baseline test but not fully eliminated. The subsurface evidence from geophysical surveys, the ice-age bathymetry of now-submerged continental shelves, and the paleoclimate context collectively suggest that these locations have been significant to human activity for far longer than their visible monuments indicate.

The data and code are openly available for independent verification and critique.

---

## Data Availability

- **Code:** [GITHUB LINK]
- **Megalithic Portal:** megalithic.co.uk (£10/year membership required for full access)
- **Pleiades:** pleiades.stoa.org (free, CC BY 3.0)
- **UNESCO:** whc.unesco.org (free)

---

## Acknowledgments

The authors thank Andy Burnham and the Megalithic Portal community for maintaining the database that makes this analysis possible. The Pleiades Gazetteer is supported by the National Endowment for the Humanities. Jim Alison is acknowledged for documenting the original great circle hypothesis circa 2001.

---

## References

Alison, J. (c. 2001). The prehistoric alignment of world wonders. jim-alison.com.

Alsadi, H., Gaber, A., El-Barbary, S., et al. (2025). Archaeological exploration via integrated shallow geophysical methods in the Saqqara necropolis, Egypt. *Heritage Science,* 13, s40494-025-01605-1. DOI: 10.1186/s40494-025-01605-1.

Askari Chaverdi, A., & Callieri, P. (2019). Persepolis West (Fars, Iran): Report on the field work (2005–2011) and results of the geophysical and archaeological surveys. *Electrum,* 26, 165–200.

Avner, U. (2001). Sacred stones in the desert: The masseboth of the Negev and Sinai. In P. M. M. Daviau, J. W. Wevers, & M. Weigl (Eds.), *The World of the Aramaeans II: Studies in History and Archaeology in Honour of Paul-Eugène Dion* (pp. 103–143). Sheffield Academic Press.

Avner, U. (2002). Studies in the material and spiritual culture of desert societies: The Negev and Sinai populations, 6th–3rd millennia B.C. Doctoral dissertation, Hebrew University of Jerusalem.

Eid, J., Abdelazeem, M., El-Haddad, A., et al. (2021). Archaeological investigation and hazard assessment using integrated geophysical techniques at the Bent Pyramid Valley Temple, Dahshur, Egypt. *Frontiers in Earth Science,* 9, 674953. DOI: 10.3389/feart.2021.674953.

Fortin, M.J., & Dale, M.R.T. (2005). *Spatial Analysis: A Guide for Ecologists.* Cambridge University Press.

Hancock, G. (1995). *Fingerprints of the Gods.* Crown Publishers.

Lasaponara, R., Masini, N., Orefici, G., & Quesada, M. (2011). New discoveries in the Piramide Naranjada in Cahuachi (Peru) using satellite, ground probing radar and magnetic investigations. *Journal of Archaeological Science,* 38(8), 2031–2039. DOI: 10.1016/j.jas.2011.04.011.

Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. *Cancer Research,* 27(2), 209–220.

Morishima, K., Kuno, M., Nishio, A., et al. (2017). Discovery of a big void in Khufu's Pyramid by observation of cosmic-ray muons. *Nature,* 552, 386–390. DOI: 10.1038/nature24647.

Procureur, S., Morishima, K., Tayoubi, M., et al. (2023). Investigation of the North Face Corridor in the Great Pyramid of Giza by cosmic-ray muon radiography and endoscopy. *Heritage Science,* 11, 72. DOI: 10.1186/s40494-023-00948-1.

Rose, J.I. (2010). New light on human prehistory in the Arabo-Persian Gulf Oasis. *Current Anthropology,* 51(6), 849–883. DOI: 10.1086/657397.

Sato, M., El-Desoky, H.M., Salem, A., et al. (2024). GPR and ERT exploration in the Western Cemetery in Giza, Egypt. *Archaeological Prospection,* Early View. DOI: 10.1002/arp.1940.

Screening Eagle Technologies. (2025). Screening Eagle's GS8000 finds possible underground tunnels at Sacsayhuamán. Technical case study, published online January 2025.

Shady Solís, R., Haas, J., & Creamer, W. (2001). Dating Caral, a preceramic site in the Supe Valley on the central coast of Peru. *Science,* 292(5517), 723–726.

Tassie, G.J., & De Trafford, A. (2011). Survey and excavation at Maadi. *Archéo-Nil,* 21, 87–110.

Wheeler, M. (1968). *The Indus Civilization.* Cambridge University Press.
