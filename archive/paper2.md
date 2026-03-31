# From Migration to Monumentality: The Alison Great Circle as a 60,000-Year Human Corridor

**Author:** Elliot Allan
**Affiliation:** Independent Researcher, New York, United States
**Correspondence:** ellallan@proton.me
**Preprint (companion study):** doi.org/10.5281/zenodo.19046176
**Data and code:** github.com/thegreatcircledata/great-circle-analysis

---

## Abstract

A great circle defined by its pole at 59.682122°N, 138.646087°W connects ancient monumental sites — Giza, Nazca, Easter Island, Persepolis, Mohenjo-daro — at statistically significant concentrations, while contemporaneous settlements do not cluster on the same path (monument enrichment 5.05×, Z = 11.83; settlement enrichment 0.78×, Z = −0.95; divergence replicated across seven independent databases; Allan 2026, preprint). This study investigates the mechanism. We present twelve independent analyses organized in three categories.

**Temporal depth.** The corridor is enriched for archaeological sites dating to the initial Out-of-Africa southern dispersal (~60,000 BP; Z = 4.42, p = 0.001), including Mololo Cave, Waigeo Island, West Papua, Indonesia (55,000–50,000 BP, 12 km from circle; Gaffney et al. 2024). The circle passes closer to independent centers of plant and animal domestication than 99.8% of random great circles (p = 0.002), with nine origins within 500 km across three continents. It approximates the optimal early Holocene overland route between the Levant and the Indus via a southern arc bypassing Mesopotamia (p = 0.008). Bronze Age monumental sites are significantly co-located with early Holocene occupation sites along the corridor (p = 0.004). Five of six tested epochs from 130,000 BP to the present show corridor enrichment, with the single gap coinciding with the Last Glacial Maximum when the Arabian corridor was impassable.

**Geometric precision in Egypt.** The circle crosses the Nile within 10 km of the Memphis necropolis — the densest ancient monument complex on Earth — at a probability of 0.003 among random circles through Egypt. It intersects the recently discovered Ahramat Branch of the Nile near-perpendicularly at the pyramid field (intersection angle 70–86° depending on scale; p = 0.00016). Directional spatial analysis using Ripley's K-function reveals that pyramids are not merely near the circle but arranged along its axis (parallel/perpendicular clustering ratio = 0.08, p < 0.0001) — geometric evidence beyond simple proximity. Egyptian administrative capitals (nome centers) show no enrichment (p = 0.54), confirming the signal is monument-specific, not geographic.

**Hypothesis elimination.** Eleven alternative explanations return null results: astronomical alignment (four sub-tests), geomagnetic correlation (modern and historical field), gravitational anomaly, oral tradition clustering, Nazca geoglyph orientation (confounded with solstice sunrise), Polynesian navigation alignment, subterranean architecture, pre-Younger Dryas corridor activity, Giza precessional longitude grid, and Bronze Age trade routes. Underground architecture follows geology, not the circle. The Great Circle is less aligned with the galactic plane than 92% of random circles.

We propose that the Alison Great Circle is consistent with the southern coastal route of the Out-of-Africa dispersal — a geographic corridor first used ~60,000 years ago that attracted recurrent human activity through every subsequent epoch. The route-fit test cannot statistically distinguish the southern from the northern dispersal route in isolation (p = 0.44 for the differential), but the convergence of deep-time site enrichment, least-cost path analysis, and Holocene continuity evidence all favor the southern arc. The monument-settlement divergence reflects the final and most visible layer of this history: communities building their most permanent structures at locations that had been culturally significant for tens of millennia. The corridor was not designed. It was inherited.

**Keywords:** spatial statistics, Out-of-Africa migration, southern dispersal route, monument-settlement divergence, great circle, Monte Carlo simulation, archaeological corridor, agricultural origins, Egyptian Old Kingdom, Ahramat Branch, human migration

---

## 1. Introduction

### 1.1 A Pattern Without an Explanation

A great circle defined by its pole at (59.682122°N, 138.646087°W) passes within 50 km of Giza, the Nazca plateau, Easter Island, Persepolis, and Mohenjo-daro — five of the most studied ancient monumental complexes on Earth. This geometric coincidence was first documented by Jim Alison circa 2001 and has circulated widely in popular literature, attracting both enthusiastic speculation and scholarly dismissal.

A companion study (Allan, 2026, preprint; hereafter the "companion study") subjected this claim to rigorous statistical testing for the first time. Using distribution-matched Monte Carlo simulation across seven independent archaeological databases totaling over 550,000 entries, it found that the overall clustering of sites near the circle is explained by the circle's path through habitable, populated terrain — not anomalous in itself. But within those populated corridors, a specific pattern emerges: ancient monumental sites cluster at 5.05× the expected rate (Z = 11.83) while contemporaneous settlements fall below random (0.78×, Z = −0.95). This divergence — monuments enriched, settlements depleted, in the same geographic regions — is replicated on seven independent databases and is not produced by any of 10,000 random great circles tested. The companion study concluded: "The anomaly is concentrated in the Egypt-to-Iran corridor and has no identified explanation."

This paper addresses that gap. The monument-settlement divergence is real and robust. What causes it?

### 1.2 Two Categories of Explanation

Proposed mechanisms for any great circle alignment of ancient sites fall into two broad categories.

The first category invokes a **physical property of Earth**: astronomical alignment (the circle corresponds to an ecliptic, a precession cycle, or a significant stellar configuration), geomagnetic correlation (the circle follows a magnetic field contour or anomaly), gravitational structure, or geological feature. If any physical property explains the pattern, we would expect it to predict not just the monument clustering but the specific geographic distribution of the signal.

The second category invokes **human activity**: the circle captures a corridor that was geographically significant to ancient peoples — a migration route, a trade network, a pilgrimage path, or a travel corridor whose importance accumulated over time. This category predicts that the corridor should show evidence of human activity that predates the monuments, increasing in density and permanence as time progresses.

We test both categories systematically.

### 1.3 The Out-of-Africa Southern Dispersal Route

The leading model for human dispersal from Africa proposes two waves: a northern route through the Levant into Eurasia, and a southern coastal route from East Africa across the southern Arabian Peninsula, along the Iranian coast, through South Asia, into Southeast Asia, and ultimately into Sahul — the combined Australia-Papua New Guinea landmass connected by lower Pleistocene sea levels (Reyes-Centeno et al. 2014; Malaspinas et al. 2016; Field & Lahr 2006; Groucutt et al. 2015). The southern route is supported by genetic evidence (deep divergence of Aboriginal Australian and Papua New Guinean lineages from other non-African populations; Rasmussen et al. 2011; Malaspinas et al. 2016), archaeological evidence (Jebel Faya, UAE, ~125,000–100,000 BP; Jwalapuram, India, ~74,000 BP; Mololo Cave, Waigeo Island, West Papua, Indonesia, ~55,000–50,000 BP; Gaffney et al. 2024; O'Connell et al. 2018), and cranial morphological data (Reyes-Centeno et al. 2014).

The geographic arc of the southern route — East Africa → southern Arabia → coastal Iran → South Asia → Southeast Asia → Sahul — closely resembles the Great Circle's path through the Old World. This resemblance has never been quantitatively tested. We do so here.

### 1.4 The Early Holocene and Agricultural Origins

Between approximately 12,000 and 7,000 BP, the global climate warmed following the Last Glacial Maximum (LGM), sea levels rose ~120 m as ice sheets melted, and the southern Arabian Peninsula — previously a hyperarid barrier — experienced a period of increased humidity now termed "Green Arabia" (Groucutt & Petraglia 2012; Rosenberg et al. 2011). During this period, the southern corridor through Arabia and Persia was traversable overland for the first time since the early Upper Paleolithic. Concurrently, multiple independent centers of plant and animal domestication emerged across Eurasia and the Americas (Larson et al. 2014; Purugganan & Fuller 2009). We test whether these agricultural origins cluster along the Great Circle corridor.

### 1.5 Study Overview

We present twelve analyses: four addressing temporal depth (how old is the corridor?), four addressing geometric precision in Egypt (how specific is the Egyptian signal?), and four addressing physical-property and cultural-transmission hypotheses (what doesn't explain it?). All analyses use distribution-matched Monte Carlo simulation with latitude-profile-matched random circles unless otherwise specified. All code and data are publicly available (github.com/thegreatcircledata/great-circle-analysis).

---

## 2. Data

### 2.1 Databases from the Companion Study

The following databases are described in full in the companion study; brief descriptions follow.

**Megalithic Portal** (megalithic.co.uk): 61,913 georeferenced prehistoric and ancient sites, predominantly from the United Kingdom and France (~65%). The Great Circle avoids Europe entirely, so the signal originates from the ~10% of the database covering non-European sites — a bias working against the hypothesis under test.

**Pleiades Gazetteer** (pleiades.stoa.org): 34,470 georeferenced ancient places from the NEH-funded peer-reviewed gazetteer maintained by New York University and the University of North Carolina. Sites classified as monumental (temple, sanctuary, pyramid, monument, amphitheatre, aqueduct; N = 1,853 pre-2000 BCE) or settlement (village, town, city, farm, port; N = 4,141 pre-2000 BCE).

**p3k14c radiocarbon database** (Bird et al. 2022): 180,070 radiocarbon dates with geographic coordinates aggregated from 39 regional databases (170,150 retained after filtering for valid coordinates and calibrated date ranges).

**XRONOS** (Palmisano et al. 2025): 305,400 radiocarbon records from 28,127 unique sites across 125 countries.

**Peru Ministry of Culture** (geoportal.cultura.gob.pe): 17,465 georeferenced Peruvian archaeological sites.

### 2.2 New Databases

**ArchaeoGLOBE** (Stephens et al. 2019): Global land-use reconstruction at multiple time slices (10,000, 8,000, 6,000, 4,000, 2,000 BP) across 10,000 grid cells, used to test whether corridor enrichment correlates with early agricultural land use.

**Independent domestication centers** (Larson et al. 2014): Twenty-four independently identified centers of plant and animal domestication with coordinates and approximate dates, spanning Eurasia, Africa, and the Americas.

**Published Out-of-Africa site compilations**: Sites from Groucutt et al. (2015), O'Connell et al. (2018), Petraglia et al. (2007), and Armitage et al. (2011), comprising 40 archaeological sites dated >40,000 BP along or adjacent to the proposed southern dispersal route. Key sites: Jebel Faya, UAE (~125,000–100,000 BP; Armitage et al. 2011); Jwalapuram, India (~74,000 BP); Tam Pa Ling, Laos (~86,000–68,000 BP; Freidline et al. 2023); Mololo Cave, Waigeo Island, West Papua, Indonesia (~55,000–50,000 BP; Gaffney et al. 2024); Madjedbebe, Australia (~65,000 BP; Clarkson et al. 2017).

**Ahramat Branch reconstruction** (Ghoneim et al. 2024): Reconstructed centerline of a now-extinct Nile tributary running approximately N-S through the Memphis necropolis, derived from synthetic aperture radar and sediment core analysis. Digitized from published figures at approximately 1 km waypoint resolution.

**IGRF-13**: International Geomagnetic Reference Field, 13th generation, sampled at 360 points along the Great Circle at 1° intervals.

**EMAG2v3**: Earth Magnetic Anomaly Grid, 2-arc-minute resolution crustal magnetic anomaly map (National Centers for Environmental Information).

**EGM2008**: Earth Gravitational Model 2008, free-air gravity anomaly sampled at 360 points along the Great Circle.

**CALS10k.2** (Constable et al. 2016): Spherical harmonic model of Earth's geomagnetic field for the past 10,000 years.

**Berezkin Analytical Catalogue**: Georeferenced mythological motif distributions across ~1,000 ethnic groups, compiled from published sources and Berezkin (2015). Seven target motif categories analyzed: flood/deluge, earth-diver, world axis, travelling transformer, cosmic hunt, celestial navigation, and giants/ancient builders.

**Sakai et al. 2024**: 2,308 Nazca geoglyph orientations from AI-assisted satellite analysis (doi:10.1073/pnas.2407652121).

**Egyptian Predynastic sites**: Compiled from Hendrickx (1999/2006) and Pleiades, classified by period: Badarian (4400–4000 BCE), Naqada I–III (4000–3100 BCE), Early Dynastic (3100–2686 BCE), Old Kingdom (2686–2181 BCE).

**Buried Egyptian settlements**: 58 pre-1000 BCE Egyptian settlement sites documented but inaccessible to full excavation (EES Delta Survey; Spencer 2024; Moeller 2016; Kemp 2018; 36 Delta, 22 Nile Valley).

### 2.3 Great Circle Definition

The great circle is defined by its pole at (59.682122°N, 138.646087°W), first documented by Alison (c. 2001). This pole was not optimized against any dataset in this study. Distance to the circle for each site is computed as d = |haversine(site, pole) − πR/2| where R = 6,371 km.

---

## 3. Methods

### 3.1 Distribution-Matched Monte Carlo Baseline

Following the companion study, baseline enrichment is computed by generating 1,000 synthetic site distributions by independently sampling latitude and longitude from the empirical distribution with Gaussian jitter (σ = 2°). For each synthetic distribution, the proportion of sites within threshold distance d of the Great Circle is computed. The Z-score is (observed proportion − mean baseline proportion) / standard deviation of baseline proportions. Primary threshold: 50 km. Reported at 25, 100, and 200 km for sensitivity.

For all new analyses involving a direct circle-vs-random-circle comparison (rather than circle vs. random sites), 10,000 random great circles are generated by sampling poles uniformly on the sphere, then filtered to match the Great Circle's latitude profile (fraction of the circle in each 10° latitude band ± 5%).

### 3.2 Out-of-Africa Route Fit

The consensus southern dispersal route was digitized from Reyes-Centeno et al. (2014), Tassi et al. (2015), and Field & Lahr (2006) as a series of geographic waypoints spanning East Africa to Papua New Guinea. Mean distance between the Great Circle and route waypoints was computed. Significance tested against 10,000 latitude-matched random great circles. The northern route (through Anatolia and Central Asia) was digitized from the same sources and tested identically.

### 3.3 Deep-Time Site Enrichment

For each of 40 sites >40,000 BP compiled from published sources, distance to the Great Circle was computed. Sites were tested at 300 km and 500 km thresholds. Z-score and p-value computed against 10,000 random great circles through the same geographic regions.

### 3.4 Early Holocene Campsite Characterization

p3k14c radiocarbon dates in the epoch 12,500–8,500 BP (10,500–6,500 BCE) were classified by material dated: charcoal/wood (campfire indicator), grain/seed (agriculture indicator), faunal remains (hunting/domestic indicator). On-corridor dates (≤50 km) were compared to off-corridor dates in the same temporal window by chi-square test. Z-scores computed for each material category.

### 3.5 Agricultural Origin Center Test

Mean distance from the Great Circle to 24 independent domestication centers (Larson et al. 2014) was computed and compared to 10,000 latitude-matched random great circles. Additionally, the number of centers within 500 km was computed and ranked against the random circle distribution.

### 3.6 Least-Cost Path Analysis

A cost surface was constructed for the early Holocene (~10,000 BP) using SRTM elevation data (slope cost), a sea-level adjustment of −40 m (ICE-6G_C model; Peltier et al. 2015), and aridity masks from published Green Arabia reconstructions (Groucutt & Petraglia 2012). The least-cost overland path was computed between three nodes: the Levant (32°N, 35°E), Persepolis (30°N, 52.9°E), and Mohenjo-daro (27.3°N, 68.1°E). Mean distance between the Great Circle and the least-cost path was compared to 10,000 random circles.

### 3.7 Holocene-Bronze Age Continuity

For each Bronze Age monument within 50 km of the Great Circle (3000–2000 BCE, from Pleiades), distance to the nearest early Holocene radiocarbon date (p3k14c, 12,500–6,500 BP) was computed. Monte Carlo null: 10,000 random shuffles of early Holocene site positions along the corridor, preserving the distance-to-circle distribution.

### 3.8 Temporal Decomposition

Monument-settlement divergence D = Z_monument − Z_settlement computed at 250-year resolution from 3500 BCE to 500 CE using Pleiades minDate and p3k14c calibrated dates. Sites in each temporal bin were tested independently using the standard Monte Carlo baseline.

### 3.9 Preservation Bias Test

Three variants of the Egyptian settlement database: (1) Pleiades baseline; (2) +58 documented buried settlements; (3) +350 settlements randomly placed within the Nile floodplain (25–31°N, 29–33°E), 1,000 Monte Carlo iterations.

### 3.10 Ahramat Branch Intersection Geometry

The Ahramat Branch centerline (Ghoneim et al. 2024) was digitized at ~1 km resolution. The Great Circle's intersection point with the branch was computed using line intersection on the sphere. Intersection angle was computed at four scales: local (2 km), necropolis (20 km), regional (50 km), and full branch (~200 km). Monte Carlo: 100,000 random great circles, each tested for (a) crossing the Ahramat Branch at all, (b) crossing within 10 km of a pyramid, and (c) crossing at angle >80° within 10 km of a pyramid.

### 3.11 Directional Ripley's K-Function

For all known Egyptian pyramids (N ≈ 138; coordinates from Lehner 1997 and Verner 2001), positions were decomposed into coordinates parallel and perpendicular to the Great Circle's local bearing. K(r) was computed separately in each direction at scales r = 0.5, 1, 2, 5, 10, 20, 50 km. The parallel/perpendicular ratio was compared to 10,000 random distributions of the same number of points within the same geographic bounding box using the Ripley edge correction.

### 3.12 "Why Memphis" Compound Probability

(1) Among 100,000 random great circles with poles within 20° of the Alison pole (i.e., circles through the Egypt latitude band), the fraction crossing the Nile within 10 km of Memphis (29.85°N, 31.25°E) was computed. (2) For each crossing, monument density within 20 km was computed from the Pleiades and Megalithic Portal. The compound probability of (crossing within 10 km) × (monument density at or above Memphis level) was computed.

### 3.13 Predynastic Site Distribution

All Egyptian sites classified as Predynastic or Early Dynastic were binned by period (Badarian, Naqada I, Naqada II, Naqada III, Early Dynastic, Old Kingdom) and their distances to the Great Circle computed. Enrichment (proportion within 50 km) was computed at each period and tested against the distribution-matched Monte Carlo baseline.

### 3.14 Hypothesis Elimination Tests

Each hypothesis was tested with an appropriate statistical framework and compared to a Monte Carlo null of 10,000 random great circles:

- *Astronomical*: (a) Frozen equator — true polar wander rate from Steinberger & Torsvik 2008 applied to required angular displacement; (b) Stellar horizon events — bearing at each of 5 key sites compared to rise/set azimuths of all stars with magnitude < 3 at 10 epochs using astropy precession; (c) Celestial pole projection — angular distance from Great Circle pole to ecliptic pole, galactic north pole, and all bright stars computed across the precession cycle; (d) Galactic plane alignment — mutual inclination between Great Circle and galactic plane minimized across all epochs.
- *Geomagnetic*: IGRF-13 field components (F, D, I, H) sampled along the circle; variance compared to random circles. CALS10k.2 magnetic pole distance from Great Circle pole computed at 8 epochs. EMAG2 crustal anomaly mean and variance along circle compared to random.
- *Gravitational*: EGM2008 geoid height variance along circle compared to random circles.
- *Oral tradition*: Berezkin motif enrichment ratio (on-corridor/off-corridor) at 200 km band tested against 10,000 latitude-matched random circles for each of 7 categories.
- *Nazca orientations*: Great Circle bearing at Nazca (63.14°) compared to geoglyph orientation distribution using Rayleigh test and V-test. Solstice sunrise azimuth at Nazca latitude (−14.7°) computed for context.
- *Polynesian navigation*: Marae/ahu enrichment within 50–500 km of circle tested against 10,000 random circles through the Pacific. Colonization date-distance correlation (Spearman) compared to random.
- *Subterranean*: Global underground site enrichment from Wikidata, Pleiades, Megalithic Portal, and UNESCO cave/underground heritage sites. Monumental vs. utilitarian underground divergence computed.
- *Pre-Younger Dryas*: As in companion study; 94,181 merged radiocarbon dates from p3k14c, XRONOS, and ROAD v32.
- *Giza longitude grid*: As in companion study; 508,000 sites tested against 36° multiples from Giza.
- *Trade routes*: As in companion study; p3k14c corridor density as trade-activity proxy, correlation with divergence D.

---

## 4. Results

### 4.1 The Corridor Is 60,000 Years Old

#### 4.1.1 Out-of-Africa Southern Route Fit

The Great Circle fits the consensus Out-of-Africa southern dispersal route at the 96th percentile of latitude-matched random great circles (mean distance to route waypoints: 1,186 km vs. random circle mean of 3,627 km; p = 0.042). The northern dispersal route — through the Levant, Anatolia, and Central Asia — has a mean distance from the Great Circle of 1,377 km; the southern route outperforms the northern at the 56th percentile of the difference distribution, though this differential is not individually significant (p = 0.44). Both routes fit far better than the random baseline.

We note an important interpretive limitation: the route-fit test alone cannot statistically distinguish the Great Circle's association with the southern route from a similar association with the northern route. Both dispersal routes pass through overlapping geographic regions across much of the Old World, and both fit the circle substantially better than chance. The southern-arc interpretation receives independent support from the least-cost path analysis (Section 4.1.4), which finds that the optimal early Holocene route between the Levant and the Indus follows a southern arc rather than the northern Fertile Crescent path, and from the deep-time site distribution (Section 4.1.2), which shows enrichment consistent with the southern coastal arc. Taken together, these analyses favor a southern-route interpretation, though this cannot be established by the route-fit test in isolation.

#### 4.1.2 Deep-Time Site Enrichment

Of 40 archaeological sites dated >40,000 BP compiled from published Out-of-Africa dispersal literature, 16 (40%) fall within 500 km of the Great Circle. The Z-score against 10,000 random great circles is 4.42 (p = 0.001). The 300 km threshold produces Z = 4.35 (p = 0.0025).

Two well-dated sites fall within 50 km of the circle: Mololo Cave, Waigeo Island, West Papua, Indonesia (~55,000–50,000 BP; 12 km; Gaffney et al. 2024) and Kebar Valley, Bird's Head Peninsula, West Papua, Indonesia (~26,000 BP; 16 km; Bartstra 1998). Five additional well-dated sites fall within 200 km.

Recurrent occupation across the full time range is demonstrated by epoch binning (Table 1): five of six tested epochs show corridor enrichment. The single gap — the LGM period (26,000–14,000 BP) — is expected from paleoclimate modeling: the Arabian segment of the corridor was hyperarid and effectively impassable during the glacial maximum. Corridor occupation resumes in the early Holocene as the climate ameliorated. We note that the ~12,000-year LGM interruption is consistent with the corridor hypothesis rather than problematic for it: it represents a period of geographic impassability rather than cultural discontinuity. The corridor's resumption in the early Holocene, at the same locations and with the same material signature, is itself evidence that the route reflects persistent geographic optimization — the same topographic and resource logic that made the corridor favorable before the LGM made it favorable again once conditions allowed. A corridor that showed continuous occupation through the LGM, when Arabia was uninhabitable, would be more difficult to explain than one that paused and resumed.

**Table 1. Temporal layering of corridor activity**

| Epoch | Indicator | Enrichment | p-value |
|-------|-----------|------------|---------|
| >40,000 BP | Out-of-Africa sites (N=40) | Z = 4.42 | 0.001 |
| 40,000–26,000 BP | Dispersal-era sites | Z = 3.18 | 0.008 |
| 26,000–14,000 BP (LGM) | p3k14c dates | Below random | expected |
| 12,500–8,500 BP | Early Holocene occupation (p3k14c) | +3.91σ charcoal dates | <0.001 |
| 8,000–5,000 BP | Agricultural origins | 9 centers ≤500 km | 0.002 |
| 5,000–2,000 BP | Bronze Age monuments | 5.05×, Z = 11.83 | <0.001 |

#### 4.1.3 Early Holocene Campsite Characterization

Among p3k14c radiocarbon dates in the early Holocene window (12,500–8,500 BP), on-corridor dates (≤50 km; N = 97) show a distinct material profile: significant enrichment for charcoal and wood (+3.91σ relative to off-corridor dates in the same period) and significant depletion of faunal remains (−4.67σ). Grain and seed materials are absent from on-corridor dates in this period.

Geographically, on-corridor early Holocene dates are concentrated at two sites: Beidha (21 km from circle) and Ghwair I (48 km), both substantial Pre-Pottery Neolithic A/B village settlements in Jordan, at the northern end of the southern arc and the junction of the Levant and the Arabian corridor. The charcoal/wood enrichment in the radiocarbon record likely reflects sample selection strategies at these two sites (N_sites = 2) rather than site function — both Beidha and Ghwair I are long-occupied communities with permanent architecture, not transient encampments. The material profile therefore characterizes the dating record, not the nature of occupation. What remains significant is the geographic pattern: the two sites with the densest early Holocene radiocarbon coverage in the region both fall within 50 km of the circle.

#### 4.1.4 Agricultural Origins Along the Corridor

The Great Circle passes closer to independent centers of plant and animal domestication than 99.8% of random great circles (p = 0.002). Nine of twenty-four known independent domestication centers (Larson et al. 2014) fall within 500 km of the circle, spanning three continents: zebu cattle (64 km, South Asia), taro (180 km, SE Asia), llama/alpaca (192 km, Peru), potato (194 km, Peru), and wheat/barley (201 km, Levant), among others.

The least-cost path analysis finds that the Great Circle approximates the optimal early Holocene overland route between the Levant and the Indus Valley better than 99.2% of random great circles (p = 0.008). The optimal route follows a southern arc rather than the northern Fertile Crescent route through Anatolia and Mesopotamia. This southern arc corresponds to the "Green Arabia" corridor: during the early Holocene, when Arabia received substantially more rainfall than today, a passable route ran south of Mesopotamia through what is now the Empty Quarter.

#### 4.1.5 Holocene-Bronze Age Continuity

Bronze Age monuments within 50 km of the Great Circle (N = 180 sites) sit a mean of 87 km and median of 81 km from the nearest early Holocene radiocarbon date (p3k14c, 12,000–8,000 BP; N = 124 on-corridor dates). The Monte Carlo null — generated by shuffling early Holocene site positions while preserving their distance-to-circle distribution — produces a null mean distance of 4,406 km (σ = 3,090 km; p = 0.004).

A methodological note: the null shuffle redistributes early Holocene positions along the full global extent of the corridor, which means null comparisons can place a Levantine date adjacent to a Peruvian monument, inflating the apparent ratio. The result is best interpreted as a test of whether Bronze Age monuments are significantly co-located with the specific early Holocene dates on their own continental segment of the corridor, and the p-value (0.004) is the more reliable summary statistic. The corridor shows strong spatial continuity between Holocene occupation loci and subsequent monument locations within each regional cluster (Egypt/Levant, Iran/Persia, Peru).

This result rules out the hypothesis that the Bronze Age monument clustering is a purely Bronze Age phenomenon with no antecedent. The corridor was not "discovered" by ancient monument-builders — it was already well-used. The clustering is heavily concentrated in the Egypt/Levant cluster (143 Bronze Age sites and 108 early Holocene dates at mean separation 80 km) with secondary clusters in Iran/Persia and Peru. They built at the same places where people had stopped, camped, and returned to for thousands of years.

### 4.2 Agricultural and Demographic Context

The corridor-agriculture relationship cannot be reduced to a simple "civilizations arise where farming originated" narrative. The monument-settlement divergence shows that while monumental sites cluster strongly (5.05×), settlements cluster weakly or not at all (0.78×) in the same regions. If the pattern were merely tracking civilization-level development, both monument types should cluster similarly.

Instead, the pattern suggests that the corridor attracted monumental construction specifically — ritual, commemorative, or territorial marking at locations whose significance derived from their long history of human passage, not from their agricultural productivity. The nome capital null result (p = 0.54; Egyptian administrative centers show no enrichment) supports this interpretation: where Egyptians chose to govern was driven by local resources; where they built their most permanent structures was driven by something else.

### 4.3 Geometric Precision in Egypt: The Memphis Complex

#### 4.3.1 Temporal Decomposition

At 250-year temporal resolution, the monument-settlement divergence onsets sharply at 2750–2500 BCE (monument Z = 11.26, settlement Z = −1.47; divergence = 12.73) and collapses within a single bin. Thirty-one of 38 Pleiades monuments driving this peak are Egyptian Old Kingdom pyramids in the Memphis necropolis (Giza to Dahshur, ~70 km spread). The remaining 7 are mortuary temples and royal tombs in the same complex.

This temporal sharpness was noted in the companion study as potentially anomalous — a "spike" rather than a sustained signal. The Ahramat Branch intersection (Section 4.3.2) provides a geometric explanation: the spike reflects the convergence of two alignment systems at one point, not a diffuse corridor-wide phenomenon.

#### 4.3.2 Ahramat Branch Intersection

Ghoneim et al. (2024) identified the Ahramat Branch, a now-extinct Nile tributary running approximately N-S through the Memphis necropolis that provided waterfront access to Old Kingdom pyramid construction sites. The Great Circle intersects this branch at 29.921°N, 31.149°E — near modern Zawyet el-Aryan, between Giza and Abu Sir — at the following angles depending on measurement scale:

- Local (2 km segment): 54° (oblique)
- Necropolis axis (20 km): 70°
- Regional (50 km): 81°
- Full branch bearing: 86°

Fifteen pyramids fall within 10 km of the intersection point; 27 within 20 km. The intersection occurs at the densest segment of the entire pyramid field.

Monte Carlo significance: among 100,000 random great circles, only 0.47% cross the Ahramat Branch at all. Of those that cross the branch, the probability of a crossing within 10 km of a pyramid at an angle >70° (necropolis-scale perpendicularity) is 0.00082. Conditioning on the circle passing through the Egypt latitude band, the compound probability of the observed configuration is 0.00016 (approximately 1 in 6,250).

We note a statistical caveat: this compound probability treats the three conditions (crossing the Ahramat Branch; crossing within 10 km of a pyramid; crossing at angle >70°) as approximately independent. In practice, these are not fully independent — the Ahramat Branch runs through the pyramid field, so a circle that crosses the branch near a pyramid is more likely to do so simply by passing through that geographic region. The figure should therefore be treated as an approximation. A more conservative framing: the probability of any random circle crossing the branch within 10 km of a pyramid at the observed angle is p ≈ 0.00016 among all great circles through the Egypt latitude band.

We note that the "perpendicular crossing" claim is scale-dependent and should be qualified: at the local 2 km scale, the crossing is oblique (54°); at the regional scale relevant to the entire pyramid field, it is near-perpendicular (81–86°). Both values are reported for transparency.

#### 4.3.3 Directional Ripley's K

Applying Ripley's K-function to the positions of all known Egyptian pyramids (N ≈ 138) with decomposition into circle-parallel and circle-perpendicular coordinates reveals extreme directional clustering. At all scales from 0.5 to 50 km, pyramids cluster along the direction of the Great Circle's local bearing (parallel K-function Z > 18, p < 0.0001). The ratio of parallel to perpendicular clustering is 0.08 — meaning the pyramid distribution is effectively a razor-thin linear strip aligned with the Great Circle axis. This provides geometric evidence beyond simple proximity: the pyramids are not merely near the circle, they are arranged in its direction.

A Monte Carlo test (10,000 random distributions of 138 points within the Memphis corridor bounding box) produces a parallel/perpendicular ratio exceeding 0.08 in fewer than 0.01% of trials.

#### 4.3.4 "Why Memphis"

Memphis — capital of Egypt for over three millennia and site of the highest concentration of monumental architecture in the ancient world — sits at the specific point where the Great Circle crosses the Nile. Among 100,000 random great circles passing through the Egypt latitude band (24°N–32°N), only 0.75% cross the Nile within 10 km of Memphis. Among those that cross within 10 km of Memphis, monument density within 20 km at the crossing point matches or exceeds the observed Memphis density in fewer than 0.5% of cases. The compound probability — crossing near Memphis AND at the densest monument concentration — is approximately 0.003.

We do not claim that the Great Circle was the reason Memphis was sited where it was; the conventional explanations (apex of the Nile Delta, strategic position between Upper and Lower Egypt, proximity to the Ahramat Branch waterway) are well-supported. Rather, we note that the Great Circle, a 40,026 km circumference great circle defined by a point in Alaska, manages to cross the Nile at the specific location that became Earth's most monumental city — and that this is unlikely by chance.

#### 4.3.5 Predynastic Emergence

The monument-settlement divergence in Egypt does not appear suddenly at the Old Kingdom. Among Predynastic sites compiled from Hendrickx (1999/2006), zero sites from the Badarian (4400–4000 BCE) or Naqada I (4000–3500 BCE) periods fall within 50 km of the Great Circle. The first on-corridor sites appear at Naqada II (~3500–3200 BCE): Ma'adi-culture settlements and proto-urban centers in the Memphis region. Enrichment intensifies through Naqada III (3200–3100 BCE) and the Early Dynastic period (3100–2686 BCE) before reaching full expression in the Old Kingdom (Z = 2.44, p = 0.034 for Predynastic-Early Dynastic combined, N = 23 sites).

The small sample sizes in pre-Old Kingdom periods preclude strong statistical claims. However, the directional trajectory — from zero to gradual accumulation to peak — is consistent with the long-term corridor model: locations on the circle became progressively more significant over centuries before the pyramid-building era.

#### 4.3.6 Egyptian Administrative Geography as Negative Control

Ancient Egypt's ~42 nome (administrative province) capitals (coordinates from Baines & Malek 2000 and Pleiades), which mark where the Egyptian state chose to locate its governance infrastructure, show no enrichment near the Great Circle (p = 0.54). This is the expected result if the signal reflects monument-specific placement rather than general geographic attractiveness. Administrative capitals were located according to local resource logic — river crossings, fertile floodplain segments, defensible positions — independent of any corridor. The monument-settlement divergence is not a property of Egypt's geography in general; it is specific to where Egyptians built their most permanent and ritually significant structures.

### 4.4 Hypothesis Elimination

**Table 2. Hypotheses tested and outcomes**

| # | Hypothesis | Test | Key Statistic | Result |
|---|-----------|------|---------------|--------|
| 1 | Frozen equator (true polar wander) | Angular displacement / TPW rate | Requires ~30 Myr; human timescales impossible | NULL |
| 2 | Stellar horizon alignment | Best star at each of 5 sites, 10 epochs | Best match: Procyon at Giza, 0.1° offset, 11th percentile | NULL |
| 3 | Celestial pole projection | GC pole projected through precession cycle | Min distance to ecliptic pole: 18.7°; no bright star within 3° at any epoch | NULL |
| 4 | Galactic plane alignment | Mutual inclination, minimized over time | Min inclination 37.8°; less aligned than 92% of random circles | NULL |
| 5 | Geomagnetic field (modern) | IGRF-13 variance along circle | Field variance: 11th–12th percentile; 55° from magnetic equator | NULL |
| 6 | Geomagnetic field (historical) | CALS10k.2, 8 epochs | Closest to magnetic equator: 57.45° at 1000 BCE (89th percentile) | NULL |
| 7 | Crustal magnetic anomaly | EMAG2v3 profile | Variance: unremarkable; no contour-following | NULL |
| 8 | Gravity anomaly | EGM2008 geoid heights | Variance: 60th percentile | NULL |
| 9 | Oral tradition clustering | Berezkin motifs, 7 categories | Best: travelling transformer, 1.07×, p = 0.29 | NULL |
| 10 | Nazca geoglyph orientation | Rayleigh/V-test, GC bearing 63.14° | 60–70° bin most populated but GC bearing = solstice sunrise ± 2.5° | CONFOUNDED |
| 11 | Polynesian navigation alignment | Marae enrichment, colonization wavefront | Significant only at Easter Island (9.6 km); not Pacific-wide | NULL (local) |
| 12 | Subterranean architecture | Global underground sites | Z ≈ 0 across all thresholds; geology predicts, not circle | NULL |
| 13 | Pre-Younger Dryas corridor | 94,181 merged radiocarbon dates | Z = −0.47, 33rd percentile before 12,800 BP | NULL |
| 14 | Giza longitude grid | 508,000 sites vs. 36° multiples | Giza: 71st–79th percentile; settlements MORE aligned than monuments | NULL |
| 15 | Bronze Age trade routes | Corridor density correlation | Pearson r = +0.04, p = 0.93; trade-relevant sites show D = −4.00 | NULL |

**Physical Earth properties** (Hypotheses 1–8): No tested physical property of Earth accounts for the pattern. The Great Circle is not a frozen equator (this would require true polar wander at geological rather than human timescales — approximately 30 million years to displace the pole by the required 30°). No bright star aligns with the circle's bearing at multiple independent sites at the same epoch. The circle's projected pole traces no path near the ecliptic pole, galactic north pole, or any notable bright star across the full 26,000-year precessional cycle. Geomagnetic field variance along the circle is in the 11th percentile — not unusually low or high. Crustal magnetic anomalies show no correlation. Gravity anomaly variance is entirely unremarkable (60th percentile). The circle is less aligned with the galactic plane than 92% of random circles.

**Cultural transmission** (Hypotheses 9–12): Mythological motif distributions do not cluster along the corridor. The Nazca geoglyph result is a confounded null: the Great Circle's bearing at Nazca (63.14°) falls 2.5° from the June solstice sunrise azimuth at Nazca's latitude (65.6°), making these two hypotheses statistically indistinguishable. The Polynesian signal is entirely driven by Easter Island (9.6 km from circle, ~900 moai), with no significant enrichment across the broader Pacific or correlation with colonization sequence — Easter Island was settled late relative to other Polynesian islands despite being on the circle. Subterranean architecture follows local geology (volcanic tuff in Cappadocia, soluble limestone in the Levant) rather than the circle; Cappadocia's famous underground cities are 905–991 km from the circle.

**Catastrophist and alternative hypotheses** (Hypotheses 13–15): The corridor shows no anomalous activity before 12,800 BP — fewer radiocarbon dates than the random-circle mean before the Younger Dryas (Z = −0.47, 33rd percentile). There is no pre-YD civilization signal. The Giza precessional longitude grid is not supported. Bronze Age trade routes, when tested as an alternative explanation for the divergence, produce the opposite pattern: sites relevant to trade activity are anti-correlated with the divergence (D = −4.00).

---

## 5. Discussion

### 5.1 The Corridor Hypothesis

The convergence of independent lines of evidence points toward a single explanation: the Alison Great Circle is consistent with the Out-of-Africa southern dispersal route — the geographic corridor first used by anatomically modern humans as they colonized the eastern hemisphere approximately 60,000–70,000 years ago. As discussed in Section 4.1.1, the route-fit test alone cannot statistically distinguish southern from northern route association. The southern interpretation is supported by the independent convergence of the least-cost path (Section 4.1.4), deep-time site enrichment along the coastal arc (Section 4.1.2), and Holocene material patterns in the Arabia-to-Levant segment.

The evidence is layered across five independent time scales. The earliest layer (>40,000 BP) consists of the Out-of-Africa dispersal sites themselves — Mololo Cave in West Papua at 12 km, Kebar Valley (West Papua) at 16 km, Jwalapuram through a southern arc connecting the same regions as the Great Circle. The Holocene layer (12,500–6,500 BP) shows charcoal-dominant radiocarbon dates along the corridor, consistent with occupation nodes rather than agricultural settlements. The agricultural layer (8,000–5,000 BP) shows that nine independent centers of food production arose within 500 km of the corridor — a pattern more significant than 99.8% of random circles. The continuity layer (Bronze Age monuments significantly co-located with early Holocene occupation sites; p = 0.004) shows that specific locations mattered across millennia. And the monumental layer (the pyramid field itself, arranged along the circle's axis) is the final, most visible expression of what had been a significant geographic corridor for tens of thousands of years.

No single analysis in this chain is decisive alone. The route fit (p = 0.042) is marginal and would not survive stringent correction for multiple comparisons. But the deep-time enrichment (p = 0.001), domestication center proximity (p = 0.002), least-cost path (p = 0.008), and Holocene-Bronze Age continuity (p = 0.004) are independent tests using different data and different methodologies — and they all point the same direction. The probability that five independent tests all support the corridor hypothesis by chance decreases with each additional positive result.

### 5.2 Why Monuments, Not Settlements?

The monument-settlement divergence — the central puzzle from the companion study — is interpretable in this framework. Monuments are not placed where people live; they are placed at locations of cultural significance. Along a 60,000-year migration corridor, the nodes where travelers historically stopped, rested, and gathered would accumulate cultural significance long before any pyramid was laid.

An analogy: medieval European cathedrals cluster along the Camino de Santiago pilgrimage route; towns and villages fill the surrounding countryside according to agricultural and commercial logic. The cathedrals mark the corridor; the towns do not. Similarly, the monumental sites along the Great Circle may mark the corridor — ritual intensification at locations whose significance derived from thousands of years of human passage — while settlements clustered according to local resource availability, independent of the corridor.

The nome capital result supports this interpretation precisely: where Egyptians governed followed local geographic logic (p = 0.54, null). Where they built their most permanent structures followed something else — something that consistently pointed toward the corridor.

We anticipate the objection that early Holocene occupation sites cluster along the corridor (Section 4.1.3) while later settlements do not — and that this temporal shift requires explanation. The corridor hypothesis provides one: as populations grew and agricultural land use intensified through the Neolithic and Bronze Age, settlement siting became increasingly governed by local resource availability (arable land, water, defensible terrain) rather than by corridor proximity. Monument siting, by contrast, remained anchored to locations of accumulated cultural significance — the nodes along the corridor where human activity had been concentrated for millennia. The absence of settlement clustering is therefore not a weakness of the corridor hypothesis but a core prediction: settlements follow resources, monuments follow memory (cf. Bradley 1998 on the role of accumulated cultural significance in monument placement).

### 5.3 The Memphis Intersection as Crossroads

The concentration of the signal at the Memphis necropolis — where the Great Circle intersects the Ahramat Branch near-perpendicularly at the densest part of the pyramid field — is the most geometrically striking feature of the entire pattern. Two independent alignment systems converge at one point: the N-S Ahramat Branch (a practical waterway that explained pyramid siting for millennia) and the E-W Great Circle (whose mechanism is the subject of this paper).

The Ahramat Branch explains HOW the pyramids were built (water access for material transport). The Great Circle corridor explains WHERE the entire complex was located in Egypt's landscape. The intersection of the two at Memphis — a 1-in-6,250 geometric coincidence — is consistent with the corridor hypothesis: the Great Circle traces a route that was geographically significant long before the Ahramat Branch was a relevant factor, and Memphis's position at the junction of Upper and Lower Egypt placed it at the intersection of both.

### 5.4 The Levant as the Corridor's Northern Hinge

The early Holocene occupation concentration at Beidha and Ghwair I in Jordan — Pre-Pottery Neolithic A/B village settlements at the junction of the Levant and the Arabian corridor — provides a geographic anchor for the human movement interpretation. Both are substantial, permanently occupied communities with architectural sequences spanning centuries; the charcoal-dominant radiocarbon record at these sites reflects dating strategies rather than site function. Their position at the northern end of the Green Arabia corridor, where travelers would have entered or left the Arabian peninsula, is nonetheless consistent with a role as key nodes on a travel route, whether or not the communities themselves were transient.

This also provides context for why the signal is strongest in the Egypt-to-Iran segment rather than the full global circle. The Old World overland portion of the southern dispersal route — from the Levant through Arabia to Persia — was the most geographically constrained section, where the corridor was narrowest. In Africa (west of the Levant), the continent's width offered many alternative routes. In South Asia and beyond, coastal and riverine options dispersed human movement. But through Arabia, the route was funneled — a geographic bottleneck that produced the highest concentration of significance per unit distance.

### 5.5 What This Does Not Claim

This paper does not claim that ancient peoples consciously placed monuments along a great circle alignment. It does not claim that the Egyptians, the Andeans, or the builders of Easter Island knew of each other or shared a cosmological tradition. It does not claim that a lost civilization mapped the Earth's surface and transmitted coordinates to successor cultures. All such interpretations are inconsistent with the evidence: the pre-Younger Dryas corridor is empty (ruling out any Ice Age antecedent civilization), the oral traditions null result rules out cultural transmission of the corridor's geometry through mythology, and the monument-settlement divergence is geographically concentrated in the Old World, not a global phenomenon.

The claim is more parsimonious: a geographic corridor that was optimal for human movement in the early Holocene, and that attracted the earliest human movement through the same region 60,000 years ago, became the location of the most intensive monumental construction the ancient world produced. The corridor shaped human activity across deep time. The monuments are its most recent signature.

We also address a natural counter-hypothesis: that the Great Circle merely passes through the most archaeologically productive latitude band on Earth (~25–35°N), and that the apparent signal is a property of that latitude rather than of a specific great circle. The latitude-profile matching of our 10,000 random circles controls for this directly: every random circle in the null distribution captures a similar fraction of each 10° latitude band (±5%), ensuring that circles with comparable continental and latitude coverage are included. The enrichment signal survives this control, meaning it cannot be reduced to latitude-band effects alone. Nevertheless, we acknowledge that a complementary test — comparing enrichment within ±50 km of a latitude line (e.g., 30°N) to enrichment within ±50 km of the Great Circle through the same region — would further strengthen the case. This test is included in the companion study's supplementary materials, where the Great Circle outperforms the nearest latitude line at all tested thresholds.

### 5.6 Limitations

The Alison Great Circle was identified post-hoc from the observation that it passes through well-known monumental sites. This raises the concern of post-hoc pattern-matching — the "Texas Sharpshooter" fallacy, in which a target is drawn around an existing cluster and the resulting alignment treated as meaningful. We address this in three ways. First, the null-hypothesis framework tests the circle against 10,000 latitude-matched random great circles, not against a baseline of zero; the enrichment signal must exceed what any comparably positioned circle produces, not merely exist. Second, the monument-settlement divergence is not predicted by post-hoc selection: if the circle were merely passing through populated areas, both monuments and settlements should cluster, yet settlements fall below random (0.78×, Z = −0.95). A sharpshooter cannot paint a target around monuments without also enclosing the settlements beside them. Third, the deep-time enrichment (Z = 4.42, p = 0.001) predates the monuments by 50,000 years and was not part of the original observation that motivated the circle's identification. Nevertheless, this study is exploratory, not confirmatory. The twelve tests were not pre-registered, and the results should be treated as hypothesis-generating rather than hypothesis-confirming. Independent pre-registered replication — ideally using a corridor identified from migration data alone, without reference to monument locations — would substantially strengthen the case.

The route fit p-value (p = 0.042) is marginal and would not survive Bonferroni correction applied across all analyses. Applying a Benjamini-Hochberg false discovery rate correction across the twelve primary analyses, the route fit would not survive at FDR q < 0.10; however, the seven remaining tests (deep-time enrichment, domestication proximity, least-cost path, Holocene continuity, Ahramat intersection, Ripley's K, and "Why Memphis") all survive at q < 0.05. The argument therefore does not depend on any single test but on the convergence of independent lines of evidence, each with its own dataset and methodology.

Deep-time archaeological sites are sparse — 40 sites covering 60,000 years across the entire Old World. The enrichment test is accordingly low-powered in specific regions. The Arabian segment of the southern route is particularly poorly sampled archaeologically, partly because site preservation is poor in hyperarid environments and partly because much of the region was submerged under the Persian Gulf (now confirmed dry during the LGM) or under alluvial deposits.

Genetic evidence was considered but is uninformative for this analysis: post-Neolithic population movements (Indo-European expansion, Austronesian expansion, Bantu expansion) have thoroughly overwritten any deep-time migration signal in modern genomes, and ancient DNA coverage along the southern route remains insufficient for corridor-specific analysis. We include this as a testable prediction (Section 5.7) rather than a current test.

The Predynastic emergence finding (p = 0.034) is based on 23 sites and should be treated as suggestive pending larger-sample analysis.

The monument-settlement divergence pattern is not replicated in South America in the primary analysis (three independent databases; companion study). The Peruvian monuments on the circle — Nazca, Machu Picchu, Cusco — are present, but the divergence D value does not reach significance with standard classification methodology. We note that the companion study's New World analysis was sensitive to site classification approach: under alternative classification schemes, weak positive signals emerged in the Peru Ministry dataset (D = 6.50) and XRONOS New World subset (D = 4.69). The southern dispersal route does not extend to the Americas in any standard model, consistent with the primary absence of the divergence pattern there, and the New World signal is therefore treated as null unless replicated with more complete data.

### 5.7 Testable Predictions

If the corridor hypothesis is correct, several predictions follow:

1. Systematic survey of the Arabian segment of the circle (UAE, Oman, coastal Arabia) should yield early Holocene occupation sites consistent with the campsite material signature identified in the Levant.
2. GPR and magnetometry surveys at Saqqara and other Memphis-area sites should find buried structures closer to the Great Circle's path than to the periphery — a testable prediction for future geophysical survey programs.
3. If the Ahramat Branch reconstruction is extended north and south of its currently documented extent, the Great Circle's intersection point should remain near-perpendicular to the branch regardless of meander, because the branch's overall N-S orientation is geologically constrained.
4. Future ancient DNA from early Holocene sites along the southern corridor should show greater affinity with Pleistocene southern-route populations (Andamanese, Papua New Guinean highland communities) than with northern Fertile Crescent farming populations.

---

## 6. Conclusions

The Alison Great Circle is consistent with the Out-of-Africa southern coastal dispersal route — the geographic corridor along which anatomically modern humans colonized the eastern hemisphere beginning approximately 60,000–70,000 years ago. The route-fit test cannot distinguish southern from northern route association in isolation, but the southern interpretation is independently supported by least-cost path analysis, deep-time site distribution, and Holocene material patterns. Statistical evidence for corridor significance comes from five independent time scales: enrichment of Out-of-Africa dispersal sites (p = 0.001), proximity to independent agricultural origins (p = 0.002), fit to the optimal early Holocene overland route (p = 0.008), and spatial continuity between Holocene occupation sites and Bronze Age monuments (p = 0.004). The corridor shows human occupation in five of six tested epochs spanning 130,000 BP to the present, with the single gap at the LGM consistent with known Arabian aridity.

In Egypt — where the signal is most concentrated — the circle crosses the Nile at the world's densest monument complex at a 1-in-300 probability, intersects the Ahramat Branch near-perpendicularly at the pyramid field (p = 0.00016), and pyramids are arranged along its directional axis rather than merely near it (directional Ripley's K ratio = 0.08, p < 0.0001). Egyptian administrative geography shows no enrichment (p = 0.54), confirming the signal is monument-specific.

Fifteen alternative hypotheses — spanning astronomical, geomagnetic, gravitational, geological, mythological, navigational, catastrophist, and cultural-transmission explanations — all return null results.

The monuments were not placed by design along a global alignment. They were built, over centuries and millennia, at locations on a corridor that had been geographically significant for 60,000 years — locations where people had always stopped, camped, gathered, and returned. The corridor came first. The monuments are its most permanent record.

---

## Acknowledgments

The author thanks Jim Alison for the original geometric observation that motivated this research. Statistical analyses were conducted in Python; code is available at github.com/thegreatcircledata/great-circle-analysis. The author has no conflicts of interest to declare.

**AI disclosure.** The author used AI language model tools (Anthropic Claude) for editorial feedback, reference verification, and manuscript review during preparation of this paper. All statistical analyses, code, data collection, and scientific interpretations are the author's own work.

---

## References

Allan, E. (2026). Statistical analysis of ancient monumental site distribution along a proposed great circle: Evidence from five archaeological databases and a hemisphere decomposition. Zenodo preprint. doi:10.5281/zenodo.19046176

Alison, J. (c. 2001). The prehistoric alignment of world wonders. Retrieved from home.hiwaay.net/~jalison/

Armitage, S.J. et al. (2011). The southern route "out of Africa": Evidence for an early expansion of modern humans into Arabia. Science, 331(6016), 453–456.

Baines, J., & Malek, J. (2000). Cultural Atlas of Ancient Egypt (rev. ed.). Checkmark Books.

Bradley, R. (1998). The Significance of Monuments: On the Shaping of Human Experience in Neolithic and Bronze Age Europe. Routledge.

Bartstra, G.J. (Ed.) (1998). Bird's Head Approaches: Irian Jaya Studies — A Programme for Interdisciplinary Research. Modern Quaternary Research in Southeast Asia, 15. A.A. Balkema.

Berezkin, Y.E. (2015). Folklore and mythology catalogue: Its networking and cognitive studies. In A. Kyriakidis (Ed.), Folklore Studies in the Twenty-First Century. Tartu University Press.

Bird, D. et al. (2022). p3k14c, a synthetic global database of archaeological radiocarbon dates. Scientific Data, 9, 27.

Clarkson, C. et al. (2017). Human occupation of northern Australia by 65,000 years ago. Nature, 547, 306–310.

Constable, C. et al. (2016). Persistent high paleosecular variation activity in Southern hemisphere for at least 10,000 years. Earth and Planetary Science Letters, 453, 78–86.

Field, J.S., & Lahr, M.M. (2006). Assessment of the southern dispersal: GIS-based analyses of potential routes at oxygen isotope stage 4. Journal of World Prehistory, 19(1), 1–45.

Freidline, S.E. et al. (2023). Early presence of Homo sapiens in Southeast Asia by 86–68 kya at Tam Pa Ling, Laos. Nature Communications, 14, 3193.

Gaffney, D. et al. (2024). Waigeo Island occupation at Mololo Cave reveals earliest known human presence in the Raja Ampat archipelago. Antiquity, 98(400), e22.

Ghoneim, E. et al. (2024). Radar topography reveals a massive ancient river system below Egypt's pyramids. Communications Earth & Environment, 5, 233.

Groucutt, H.S., & Petraglia, M.D. (2012). The prehistory of the Arabian Peninsula: Deserts, dispersals, and demography. Evolutionary Anthropology, 21(3), 113–125.

Groucutt, H.S. et al. (2015). Rethinking the dispersal of Homo sapiens out of Africa. Evolutionary Anthropology, 24(4), 149–164.

Hendrickx, S. (1999/2006). Predynastic-Early Dynastic chronology. In Hornung, Krauss, & Warburton (Eds.), Ancient Egyptian Chronology. Brill.

Kemp, B. (2018). Ancient Egypt: Anatomy of a Civilisation (3rd ed.). Routledge.

Larson, G. et al. (2014). Current perspectives and the future of domestication studies. PNAS, 111(17), 6139–6146.

Lehner, M. (1997). The Complete Pyramids. Thames & Hudson.

Malaspinas, A.S. et al. (2016). A genomic history of Aboriginal Australia. Nature, 538, 207–214.

Moeller, N. (2016). The Archaeology of Urbanism in Ancient Egypt. Cambridge University Press.

O'Connell, J.F. et al. (2018). When did Homo sapiens first reach Southeast Asia and Sahul? PNAS, 115(34), 8482–8490.

Palmisano, A. et al. (2025). XRONOS: An open data infrastructure for archaeological chronology. Journal of Computer Applications in Archaeology, 8(1), 242–263. doi:10.5334/jcaa.191

Peltier, W.R. et al. (2015). Space geodesy constrains ice age terminal deglaciation. PNAS, 112(30), 9287–9292.

Petraglia, M. et al. (2007). Middle Paleolithic assemblages from the Indian subcontinent before and after the Toba super-eruption. Science, 317(5834), 114–116.

Petraglia, M. et al. (2010). Out of Africa: New hypotheses and evidence for the dispersal of Homo sapiens along the Indian Ocean rim. Annals of Human Biology, 37(3), 288–311.

Purugganan, M.D., & Fuller, D.Q. (2009). The nature of selection during plant domestication. Nature, 457, 843–848.

Rasmussen, M. et al. (2011). An Aboriginal Australian genome reveals separate human dispersals into Asia. Science, 334(6052), 94–98.

Reyes-Centeno, H. et al. (2014). Genomic and cranial phenotype data support multiple modern human dispersals from Africa and a southern route into Asia. PNAS, 111(20), 7248–7253.

Rosenberg, T.M. et al. (2011). Humid periods in southern Arabia: Windows of opportunity for modern human dispersal. Geology, 39(12), 1115–1118.

Sakai, M. et al. (2024). New discoveries of Nazca geoglyphs using AI-assisted satellite image analysis. PNAS, 121(40), e2407652121.

Spencer, J. (2024). The Archaeological Sites of Lower Egypt: A Gazetteer. Excavation Memoir 119. London: Egypt Exploration Society.

Steinberger, B., & Torsvik, T.H. (2008). Absolute plate motions and true polar wander in the absence of hotspot tracks. Nature, 452, 620–623.

Stephens, L. et al. (2019). Archaeological assessment reveals Earth's early transformation through land use. Science, 365(6456), 897–902.

Tassi, F. et al. (2015). Early modern human dispersal from Africa: Genomic evidence for multiple waves of migration. Investigative Genetics, 6, 13.

Verner, M. (2001). The Pyramids: The Mystery, Culture, and Science of Egypt's Great Monuments. Grove Press.

---

*Word count (main text excluding references): ~9,800 words*
*Figures required: 9 (see outline)*
*Supplementary materials: directive output files, full site catalogs, Monte Carlo code*
