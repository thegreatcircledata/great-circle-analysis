# Testing a proposed transcontinental great circle alignment: Land-constrained Monte Carlo baselines, geographic explanations, and localized monument divergence at the Memphis necropolis

> **Historical record — v5 (JAS:Reports, superseded).** This manuscript was the prior Journal of Archaeological Science: Reports submission. It has since been revised and resubmitted to PLOS ONE (v8, 2026-04-17; preprint [zenodo.19212669](https://doi.org/10.5281/zenodo.19212669)). Numbers in this file may have been refined in v8; retained as historical record only.

**Author:** Elliot Allan, Deep Time Research Institute
**Correspondence:** elliot@deeptime-research.org
**ORCID:** 0009-0008-8541-0944
**Preprints:** Earlier versions of this work were deposited as preprints: v1–v4 at doi.org/10.5281/zenodo.19212669. Companion study: Allan (under review), preprint at doi.org/10.5281/zenodo.19240285.
**Data and code:** github.com/thegreatcircledata/great-circle-analysis

---

## Abstract

We test whether ancient archaeological sites cluster along the great circle of Alison (c. 2001), with pole at 59.682122°N, 138.646087°W, across eight independent databases comprising approximately 259,000 unique sites. Land-constrained Monte Carlo baselines correct for coastal jitter bias absent from prior great circle enrichment tests.

Overall enrichment is significant on three of four corrected databases (Megalithic Portal Z = 8.77, Pleiades Z = 6.74, XRONOS Z = 5.91) and marginal on a fourth (p3k14c Z = 2.09). After habitability adjustment, enrichment falls to the 78.5th percentile of matched circles: the circle's path through the Nile Valley, Mesopotamia, the Levant, and the Indus floodplain accounts for the site-count signal.

The robust finding is a monument-settlement divergence D = 9.65 (land-constrained, Pleiades pre-2000 BCE), unmatched by any of 10,000 random great-circle poles (maximum D = 0.90) and ranked at the 99.63rd percentile of 4,323 habitability-matched circles. The divergence replicates on independent databases using different classification methods, though replication strength varies with baseline correction and classification approach.

At 250-year resolution, onset is concentrated in the Memphis necropolis (2750–2500 BCE; 38 Old Kingdom pyramids spanning 0.44° × 0.55°). Two South American databases return null results. The corridor shows no preferential pre-Younger Dryas activity (Z = −0.47). Astronomical encoding, longitude grids, and pre-glacial civilization hypotheses are falsified.

A companion study (Allan, under review) identifies tectonic diversity and geographic circumscription as partial geological explanations. Land-constrained baselines — reducing Z-scores by 31–73% — are proposed as standard practice for enrichment tests at coastal or island sites.

**Keywords:** spatial point processes; monument distribution; Egyptian Old Kingdom; coastal jitter correction; corridor enrichment; site classification; habitability adjustment

---

## 1. Introduction

Many of the most celebrated ancient monuments in the world appear, when plotted on a globe, to lie close to a single great circle. Sites including the Giza pyramid complex, Nazca, Easter Island, Machu Picchu, Angkor Wat, and Persepolis have been noted by various observers over the past quarter-century as lying along or near a common geodesic. This observation was first documented systematically by Jim Alison, who identified a great circle pole at 59.682122°N, 138.646087°W and noted that it passed within tens of kilometers of a substantial number of the ancient world's most architecturally distinctive sites (Alison c. 2001). Despite the longevity and popular appeal of this observation, it has never been subjected to rigorous statistical testing against comprehensive, independently assembled archaeological databases using appropriate null models. The present study does so.

Why does a rigorous test matter? If the alignment is real — that is, if the probability of observing this degree of clustering by chance is demonstrably low under a range of null models — then it is an empirical fact about the spatial distribution of ancient monuments that demands explanation, whether that explanation involves landscape constraints, shared astronomical orientations, common ecological corridors, or something else not yet identified. If the alignment is an artifact of selection bias, post-hoc pattern matching, or inappropriate baselines, then demonstrating this rigorously is itself a contribution: it establishes which methodological standards are necessary and why simpler approaches give misleading results. Either outcome advances archaeological spatial analysis.

The specific methodological challenge is that the great circle was not defined independently of the data. Alison's circle was identified by observing famous sites and drawing a geodesic through them. Testing it against the same sites, or against a database enriched by the same culturally prominent monuments, risks confirming a pattern that was introduced by selection. Rigorous testing requires not just a significance threshold but a cascade of null models of increasing stringency: distribution-matched Monte Carlo permutation, kernel density estimation, habitability adjustment, split-sample validation against held-out data, spatial block cross-validation, monument-versus-settlement disaggregation, temporal decomposition, and systematic testing of every plausible alternative explanation. This study applies all of these.

Our primary finding is that the overall enrichment of archaeological sites near the great circle, while statistically significant on three of four land-constrained databases, is fully explained by the circle's path through the world's most habitable prehistoric corridors. This is not a trivial null result: it tells us that the circle's geometry reflects ecological logic — connecting the Nile Valley, the Levantine coast, Mesopotamia, the Indus floodplain, coastal Peru, and Easter Island on a single path — but it does not require any coordination between the societies involved. What is not explained by geography alone is the monument-settlement divergence — a finding partially addressed by geological analysis in a companion study (Allan, under review) but not fully resolved: monuments cluster near the circle while contemporaneous settlements, in the same regions at the same times, are systematically depleted. This divergence (D = 9.65, land-constrained) is confirmed on multiple databases with different classification methods, survives all alternative explanations tested, and is temporally concentrated in a single 250-year construction window in the Memphis necropolis during the Egyptian Old Kingdom.

---

## 2. Data

### 2.1 Megalithic Portal (61,913 sites)

The Megalithic Portal is a community-maintained, moderated database of megalithic and ancient monument sites worldwide. It is the largest single-database source used in this analysis. Geographic coverage is heavily biased toward the British Isles and western Europe: approximately 65% of sites fall in the UK, Ireland, and France. This is a critical methodological observation because the Alison great circle does not pass through Europe; consequently, the database's strong geographic bias works against the alignment hypothesis, making the Portal a conservative test. Sites were deduplicated using a 111-meter radius threshold applied across 62 KML export files, yielding 61,913 unique sites. Classification of sites as monumental or settlement uses the Portal's broad type category system; the imprecision of these categories, particularly the "Ancient Village or Settlement" grouping, is acknowledged as a limitation (see Section 5.7). For monument-settlement divergence calculations, a non-European subset was used to remove the dominant no-signal geographic mass.

### 2.2 Pleiades (34,470 total sites; 778 pre-2000 BCE; 406 ancient monuments)

Pleiades is a NEH-funded, peer-reviewed gazetteer of ancient places maintained by the Ancient World Mapping Center and the Institute for the Study of the Ancient World (Elliott and Gillies 2009–2025). It is the primary source for monument-settlement classification in this study because it uses a professionally developed typology inherited from the Barrington Atlas of the Greek and Roman World (Talbert 2000). Geographic coverage is concentrated in the Mediterranean basin, Middle East, and Near East. Temporal filtering was applied using the minDate field to isolate pre-2000 BCE sites (778 sites). For monument-specific analyses, 406 sites were classified as ancient monuments using Barrington Atlas typology (temple, pyramid, sanctuary, monument, tomb). Settlement sites were classified by type terms including village, town, farm, port, and city. Pleiades was the first database for which a land-constrained baseline was computed, yielding the primary divergence figure reported in this study.

### 2.3 p3k14c (36,693 unique sites from 173,946 dates)

The p3k14c database (Bird et al. 2022) aggregates radiocarbon-dated archaeological sites from 39 regional databases. The version used here yields 36,693 unique sites after SiteID-based deduplication; an earlier analysis in Paper 1 used 25,557 unique sites after latitude-longitude deduplication, a less complete method. The same 187 on-corridor sites appear in both deduplication schemes, so the change in site count does not affect the corridor-specific result but does change the background distribution against which it is evaluated. American sites have coordinates obfuscated to county centroids for privacy reasons, a known limitation for fine-scale spatial analysis. p3k14c spans approximately 50,000 years of global settlement history and is the primary source for temporal decomposition of the overall enrichment signal. It is explicitly independent of Pleiades and the Portal.

### 2.4 XRONOS (28,127 sites from 305,400 records)

XRONOS (Palmisano et al. 2025) is an independent global radiocarbon database comprising 305,400 dated records from 28,127 unique sites. It was assembled independently of p3k14c and provides a replication test for the overall enrichment finding. XRONOS includes typological information enabling monument-settlement classification by keyword matching: 2,334 sites were classified as monumental and 5,981 as settlement. The land-constrained baseline for XRONOS used 1,000 Monte Carlo trials.

### 2.5 DARE (29,760 sites)

The Digital Atlas of the Roman Empire (DARE) represents a digitization of Barrington Atlas data developed at Lund University. It uses numeric type codes for site classification, enabling monument-settlement disaggregation that is independent of the Pleiades typology while drawing on the same underlying atlas. Overall enrichment on DARE is negative — the great circle avoids the dense Roman-era settlement territory of the Italian peninsula and Roman Europe — but the monument-settlement divergence is positive, because the narrow slice of Roman-controlled territory the circle does intersect is disproportionately monumental (temples and necropoli rather than farms and towns). This dissociation between overall enrichment and monument-settlement divergence is analytically important and discussed in Section 4.2.

### 2.6 Peru Ministry of Culture (17,465 sites)

The Peruvian national heritage registry provides one of the most complete records of archaeological sites in any single South American country. Sites were classified using Spanish and Quechua keyword matching: 1,733 sites were classified as monumental and 275 as settlement. This database is central to the South America null result (Section 4.7). The registry is independent of all other databases used here and was not available in earlier analyses of the great circle hypothesis.

### 2.7 Historic England (20,026 sites)

The Historic England database of scheduled monuments and listed buildings serves as a negative control. The Alison great circle does not pass through England; zero sites fall within 500 km of the circle. The expectation is D = 0.0, confirming that the classification and divergence methodology returns null results where the geometry predicts them. This expectation is confirmed (Section 4.2).

### 2.8 LuwianSiteAtlas (483 sites)

The LuwianSiteAtlas covers Bronze Age western Anatolia, dated approximately 2000–1200 BCE. Sites receive professional classification by type. The analysis of this database tests whether the divergence extends to the Anatolian Bronze Age — a period and region that might plausibly participate in any eastern Mediterranean corridor effect. The result is a complete null: zero of 483 sites fall within 500 km of the great circle. This null is reported in full (Section 4.7).

### 2.9 Supplementary sources

Three additional data sources are used in specific sub-analyses. The preservation test (Section 4.5) adds 58 documented buried Egyptian settlement sites drawn from the EES Delta Survey (Spencer 2024), Moeller (2016), and Kemp (2018). The pre-Younger Dryas test (Section 4.10) uses a merged database of 94,181 radiocarbon dates drawn from p3k14c, XRONOS, and ROAD v32 (Vermeersch 2024). A fourth South American database (SAAID: South American Archaeological Isotopic Database, 2,429 sites) was examined during analysis but computational outputs were not retained; its results are therefore excluded from this report.

### 2.10 Databases excluded

Two databases were considered but excluded from primary analysis. OpenStreetMap archaeological layers were excluded because of contamination risk: the Alison great circle hypothesis has been discussed publicly since at least 2001, and community-edited spatial databases may contain entries added specifically because a contributor was aware of the alignment, introducing circularity. Wikidata was excluded because it explicitly cross-links with Pleiades, creating a dependency that would violate the independence assumption underlying multi-database replication. Results from these databases are reported in supplementary material.

### 2.11 Circle definition

The great circle is defined by the pole at 59.682122°N, 138.646087°W, first documented by Jim Alison (c. 2001). This pole definition was not optimized or adjusted against any database in this study. The circle itself passes through Giza, Nazca, Easter Island, Angkor Wat, Mohenjo-daro, and a number of other archaeologically significant locations, all within approximately 50 km.

### 2.12 Aggregate statistics

Across all eight primary databases, the study encompasses approximately 259,000 unique sites and over 550,000 database entries. The distinction between entries and sites matters: the 550,000 figure refers to dated records and raw database entries; the approximately 259,000 figure is the deduplicated site count used in spatial analyses.

---

## 3. Methods

Methods are presented in order of increasing analytical stringency, from basic distance metrics through multi-level cross-validation.

### 3.0 Statistical framework

Archaeological spatial analysis has a long methodological tradition (Hodder & Orton 1976; Fortin & Dale 2005; Conolly & Lake 2006). The analyses in this study treat archaeological sites as realizations of inhomogeneous spatial point processes (Baddeley et al. 2015; Diggle 2003; Thomas 1949; Crema et al. 2010). The primary analytical tools are Monte Carlo permutation tests and corridor-count statistics rather than parametric point-process models. This choice was made because the study tests a specific geometric hypothesis — enrichment within a fixed corridor — across multiple databases with heterogeneous classification systems, a setting where corridor counts provide transparent, replicable statistics with minimal distributional assumptions. Statistical evaluation of spatial alignment claims has precedent in the ley line debates of the 1980s (Williamson and Bellamy 1983). A companion study (Allan, under review) extends this analysis using formal inhomogeneous Poisson point process modeling with geographic covariates (elevation, coastline distance, river distance, latitude), confirming that the monument-settlement divergence persists after these controls with non-overlapping confidence intervals between monument and settlement intensity functions. The Monte Carlo framework used here and the point-process framework used in the companion study yield consistent conclusions.

### 3.1 Distance metric

For each site, distance to the great circle was computed as d = |haversine(site, pole) − πR/2|, where R = 6,371 km and πR/2 ≈ 10,007.5 km. This converts great-circle distance from the pole to distance from the equatorial circle defined by that pole. Analyses were conducted at corridor half-widths of 25, 50, 100, and 200 km; the primary threshold throughout is 50 km, chosen as a conservative intermediate distance that captures the main signal while remaining substantially smaller than the peak enrichment scale.

### 3.2 Distribution-matched Monte Carlo baseline

For each Monte Carlo trial, every site in the database had its coordinates independently perturbed by adding Gaussian noise (σ = 2°) to both latitude and longitude. The number of perturbed sites falling within the 50 km corridor of the fixed Alison great circle was then counted. After N trials (1,000 standard; N = 10,000 site-perturbation trials for definitive results on Portal and Pleiades), a Z-score was computed as Z = (observed_real − μ_MC) / σ_MC, where observed_real is the count of actual, unperturbed sites within the corridor, and μ_MC and σ_MC are the mean and standard deviation of the perturbed-site counts across trials. This site-jittering procedure preserves the marginal geographic distribution of the data while destroying each site's specific angular relationship to the Alison circle. It is analogous to the Mantel (1967) permutation test and is a standard procedure in archaeological corridor enrichment analysis (cf. Hewitt et al. 2020; Conolly and Lake 2006).

Note: the 10,000 random great-circle poles uniqueness test used to establish that 0 of 10,000 random great-circle poles match the divergence (Section 4.2) employed a related but distinct procedure: random poles were generated by uniform sampling on the sphere, and each random circle's divergence was computed against the actual (unperturbed) site positions. The two procedures address different questions — the site-jittering baseline asks "how enriched is this circle relative to random arrangements of these sites?"; the random-pole sweep asks "how often does a random circle match this circle's divergence?" Both are reported. The land-constrained correction (Section 3.3) applies specifically to the site-jittering approach.

### 3.3 Land-constrained baseline (primary correction)

Following anonymous feedback identifying a coastal jitter bias in the standard distribution-matched approach, we implemented a land-constrained baseline correction. The bias arises because sites near coastlines or on islands, when jittered, land in the ocean with high probability. Ocean-landing jitter attempts fail to produce on-corridor points, deflating the expected count for random circles and artificially inflating Z-scores for the observed circle. The correction procedure accepts a jittered point only if it falls on land, determined by intersection with a Natural Earth 1:10M land polygon mask rasterized at 0.1° resolution (~1.35 million lookups per second). A maximum of 50 jitter attempts were made per site per trial; if all 50 failed, the original site coordinates were retained (a conservative choice that does not inflate Z). The land-constrained baseline distributes jittered points according to the empirical marginal distributions of the source database, constrained to land surfaces; it does not impose complete spatial randomness or cluster-process assumptions on the randomised positions. Land-constrained baselines were computed for four databases: Megalithic Portal (N = 10,000 land-constrained jitter trials), Pleiades (N = 10,000 land-constrained jitter trials), p3k14c (1,000 trials), and XRONOS (1,000 trials). Land-constrained Z-scores are treated as the primary reported values throughout this paper. The correction reduced reported Z-scores substantially (see Section 4.1) without altering the directional conclusions.

### 3.4 Kernel density estimation baseline

A Gaussian kernel density estimate was computed on the spherical surface, with bandwidth determined by Scott's rule and varied between 1° and 5° to assess bandwidth sensitivity. The KDE baseline provides an alternative to the Monte Carlo approach that is less sensitive to the marginal distribution assumption. Results are reported as Z-score ranges across the bandwidth sweep.

### 3.5 Habitability-adjusted baseline

Population density at each sample point along the circle was estimated using p3k14c site density as a proxy for prehistoric habitability (cf. Klein Goldewijk et al. 2017), analogous to the corridor-density matching approach used in least-cost path Monte Carlo simulation (Moreno-Meynard et al. 2022). A set of 4,323 circles was selected from a global sample such that each matched the Alison circle's corridor in mean site density (within one standard deviation). The Alison circle's monument-settlement divergence was then ranked against the divergence distribution of these matched circles. This procedure tests whether the divergence can be explained by the circle's path through densely occupied territory.

### 3.6 Monument-settlement classification

Classification methods varied by database to make best use of each database's structure. For Pleiades, Barrington Atlas typology was used directly: monumental sites include temples, pyramids, sanctuaries, monuments, and tombs; settlement sites include villages, towns, farms, ports, and cities. For p3k14c, site name keyword matching was used. For DARE, numeric type codes from the Barrington digitization were mapped to the same monumental/settlement distinction. For the Peru Ministry of Culture database, Spanish and Quechua keyword matching was applied. For the Megalithic Portal, the database's own broad type categories were used, with the classification imprecision noted as a limitation (Section 5.7). The Monument-settlement divergence D = Z_monument − Z_settlement was computed for each database independently.

The heterogeneity of classification methods across databases is both a limitation and a form of robustness testing. If the divergence appeared only under one classification scheme, it could reflect idiosyncrasies of that scheme. That it appears under Barrington Atlas professional typology (Pleiades), community-assigned broad categories (Megalithic Portal), numeric type codes from an independent digitisation (DARE), and keyword matching on site names (p3k14c) — four different approaches to the same underlying distinction — provides stronger evidence than any single classification method could.

### 3.7 Temporal decomposition

Temporal decomposition was conducted at two resolutions. Five-hundred-year bins spanning 5000 BCE to 500 CE provided a coarse-grained view of how divergence evolved over time. Two-hundred-fifty-year bins spanning 3500 BCE to 500 CE provided the resolution necessary to isolate specific construction phases. For Pleiades, temporal assignment used the minDate field; for p3k14c, calibrated calendar years from the radiocarbon database were used. D was computed for each temporal bin; bins with fewer than the minimum site count threshold (applied separately to monuments and settlements) are noted where they affect results.

### 3.8 Preservation test

The preservation test evaluated whether differential preservation of monumental versus domestic architecture explains the monument-settlement divergence in Egypt. Three conditions were tested. Condition 1 used the baseline Pleiades Egypt data as-is. Condition 2 added 58 documented buried Egyptian settlement sites from the EES Delta Survey (Spencer 2024), Moeller (2016), and Kemp (2018) — sites known to exist but not represented in Pleiades because they are below ground. Condition 3 added 350 settlement sites placed randomly within the Nile floodplain (25–31°N, 29–33°E) via Monte Carlo simulation, representing a conservative upper bound on the number of undocumented settlements that might hypothetically exist. Each condition was evaluated across 200 (Condition 2 short run) and 1,000 Monte Carlo trials (Conditions 1, 2 long run, and 3).

### 3.9 Split-sample and spatial block validation

Three validation procedures were applied to test for post-hoc inflation. Spatial block cross-validation follows the framework of Roberts et al. (2017). Procedure A used 100 random 50/50 splits of the Megalithic Portal database; Z was computed for the Alison circle against the held-out half only, never using the training half for significance testing. Procedure B fitted 1,000 random great circles to 15-site subsets of Megalithic Portal data and then evaluated the fitted circle's Z-score against the full database; this established the distribution of post-hoc fitted Z-scores against which to compare the Alison circle. Procedure C was a leave-one-cluster-out spatial block test removing each major geographic cluster in turn before evaluating enrichment.

### 3.10 Predictive validation

The databases were split by data age to test whether the signal appears in data collected after the great circle hypothesis was published. For Pleiades, the split used record creation date (pre-2012 versus post-2012). For p3k14c, the split used the publication year of the original reference (pre-2001 versus post-2001). If the alignment were purely post-hoc, the pre-2001 data should show the signal equally well; if the alignment predicts new data beyond its original anchor sites, the post-2001 data should show an equivalent or stronger signal. The top-5 and top-10 grid cells were removed from post-2001 data as a sensitivity check.

### 3.11 Multi-scale enrichment profile

Enrichment ratio (observed / expected) and divergence D were computed across corridor half-widths from 2 km to 200 km to characterize the spatial scale of the signal, following the multiscalar approach of Bevan and Conolly (2006). The half-life was defined as the distance at which enrichment drops to the midpoint between its peak value and the background level of 1.0.

### 3.12 Astronomical alignment test

Pole declination, ecliptic and galactic plane alignments, and solstice and equinox azimuth correspondences were tested using astropy precession calculations applied at 500-year intervals from 3000 BCE to 500 CE. All tests were subjected to Bonferroni correction across the full suite.

### 3.13 Pre-Younger Dryas test

The merged database of 94,181 dates (p3k14c + XRONOS + ROAD v32) was used to test whether the corridor was preferentially occupied before the Younger Dryas onset (~12,800 BP). On-corridor date counts for the pre-12,800 BP period were compared against 1,000 random great circles. Continuity across the Younger Dryas boundary was tested by comparing the fraction of sites with dated occupation spanning 12,800–11,700 BP on-corridor versus off-corridor. Five candidate locations proposed in various Atlantis hypotheses were tested individually.

### 3.14 Longitude grid test

A proposed Giza-centered longitude grid (described in popular literature) predicts that sites should align along meridians spaced at 36° multiples from the Giza prime meridian (31.13°E). Six sub-tests were applied to 508,000 sites: whether Giza itself is a special reference longitude, whether precessional angles produce significant alignment Z-scores, whether the monument-settlement contrast matches grid predictions, and whether alignment correlates with site age as expected by ancient encoding.

### 3.15 Multiple testing correction

Benjamini-Hochberg false discovery rate correction at q = 0.05 was applied across all primary statistical tests reported in this paper. The Benjamini-Hochberg procedure was used rather than Bonferroni correction to avoid excessive conservatism given the large number of correlated tests. All primary findings survive the correction; the paleoclimate correlation does not.

### 3.16 Use of AI tools

Computational analysis was conducted using Claude Code (Anthropic). We designed all research directives specifying the hypotheses to be tested, the datasets to be used, the statistical methods to be applied, and the validation procedures to be run. All code was reviewed by us and is available at github.com/thegreatcircledata/great-circle-analysis. The land-constrained baseline correction was implemented after anonymous Reddit feedback identified the coastal jitter bias; we evaluated the feedback, judged it technically correct, and directed its implementation. All scientific judgments, interpretations, and conclusions are our own.

---

## 4. Results

### 4.1 Overall enrichment under progressive null models

Under distribution-matched Monte Carlo baselines without land-constraint correction, overall enrichment Z-scores are high across all databases. However, the coastal jitter bias substantially inflated these values. After applying land-constrained baselines, Z-scores were reduced as follows: Megalithic Portal from Z = 13.2 to Z = 8.77 (33.6% reduction); Pleiades from Z = 11.14 to Z = 6.74 (39.5%); XRONOS from Z = 12.4 to Z = 5.91 (52.3%); p3k14c from Z = 7.77 to Z = 2.09 (73.1%). These reductions are substantial but do not change the qualitative finding: three of four land-constrained databases remain significant at conventional thresholds, and the fourth is marginal (p = 0.02). All primary results in this paper use land-constrained Z-scores where computed.

The corrected Z-scores are: Megalithic Portal Z = 8.77 (land-constrained, N = 10,000 trials), Pleiades Z = 6.74 (land-constrained, N = 10,000 trials), XRONOS Z = 5.91 (land-constrained, N = 1,000 trials), p3k14c Z = 2.09 (land-constrained, N = 1,000 trials, p = 0.02). Under the KDE baseline, Z-scores vary with bandwidth between 5.3 and 7.1 for the Portal and between 4.1 and 6.3 for Pleiades. After habitability adjustment — controlling for the circle's path through densely occupied prehistoric corridors — overall enrichment falls to the 78.5th percentile of 4,323 matched circles. This is not significant. The overall site-count enrichment is geographic: the circle passes through the Nile Valley, Mesopotamia, the Levantine coast, the Indus floodplain, coastal Peru, and Easter Island, all among the most densely occupied regions of the prehistoric world. No residual overall enrichment beyond habitability remains.

Full database summary statistics are presented in Table 1.

### 4.2 Monument-settlement divergence

The monument-settlement divergence D = Z_monument − Z_settlement is the paper's core finding. On the Pleiades database, restricted to pre-2000 BCE sites under a land-constrained baseline, monument Z = 6.74 and settlement Z = −2.91, giving D = 9.65. This value is not reached by any of 10,000 random great-circle poles tested against the same data: the random distribution has mean = 0.00, standard deviation = 0.07, and maximum = 0.90. The divergence ranks at the 99.63rd percentile of 4,323 habitability-matched circles (habitability-adjusted D = 10.44). The divergence is not explained by the circle's path through populated corridors — if it were, settlements would cluster at least as strongly as monuments.

The divergence is confirmed on one geographically and typologically independent database (Megalithic Portal non-EU, D(LC) = 4.56), one typologically independent but geographically overlapping database (DARE, D(LC) = 5.94), and one database using uncorrected baselines with a different classification method (p3k14c, D = 6.18, uncorrected due to classification degeneracy under land-constrained baselines). A fourth database (Historic England) confirms the null where expected (D = 0.0). Land-constrained baselines were computed for all three additional databases (Portal, p3k14c, DARE; 1,000 trials each). The Megalithic Portal (non-EU subset, land-constrained baseline) produces monument Z = 9.06, settlement Z = 4.50, giving D(LC) = 4.56 (4,798 monument sites, 1,134 settlement sites). DARE produces monument Z = −3.28, settlement Z = −9.22, giving D(LC) = 5.94 (921 monument sites, 18,255 settlement sites) — the pattern of monuments being relatively more enriched than settlements persists even though both type-specific Z-scores are negative, consistent with the overall enrichment being negative on DARE (see Section 4.2 above). For p3k14c, the keyword-matched monument and settlement subsets (299 and 901 unique sites, respectively) contain zero sites within 50 km of the great circle in this classification; the LC divergence is therefore degenerate (D(LC) = 0.00 due to classification degeneracy, not a substantive null). The uncorrected divergence for p3k14c (D = 6.18) was computed using a different classification approach in the predecessor analysis and is retained for reference.

The DARE finding deserves emphasis. The overall enrichment on DARE is negative — the great circle avoids the densely settled Roman heartland — yet the monument-settlement divergence is positive. Within the narrow geographic slice where the circle does intersect Roman-controlled territory, the sites are disproportionately monumental: temples, necropoli, and sacred precincts rather than farms and residential settlements. This dissociation between overall enrichment and monument-settlement divergence is precisely what the divergence metric is designed to detect.

One methodological qualification: DARE and Pleiades share geographic coverage, as both derive from the Barrington Atlas (Talbert 2000). Coordinate matching within 1 km identifies substantial overlap — the majority of Pleiades sites near the great circle have DARE counterparts. DARE's value as a replication test lies not in geographic independence but in classification independence: its numeric type codes provide a separate typological scheme that independently confirms the monument-settlement distinction. The DARE divergence result should be interpreted as confirming classification robustness rather than full geographic independence.

The Old World versus New World decomposition reveals that monument-settlement divergence is an Old World phenomenon. In the New World, specifically the Megalithic Portal Americas subset, both monuments and settlements cluster equally near the circle (D = 0.13). This is consistent with a purely geographic explanation for New World enrichment: the circle's path through coastal Peru and the Andes brings it near both monumental and habitation sites with equal probability. The divergence requires both a specific site type and a specific hemisphere to manifest.

Full divergence statistics are presented in Table 2.

### 4.3 Temporal structure at 500-year resolution

Temporal decomposition at 500-year resolution reveals that the monument-settlement divergence is not a uniform feature of the corridor through time. Before approximately 3000 BCE, D values across both Pleiades and p3k14c are near zero, indicating no differential clustering. The divergence spikes sharply at 3000–2500 BCE, reaching D = 9.95 on Pleiades and D = 3.83 on p3k14c. It then collapses by 2000 BCE. A secondary divergence signal appears at 500 BCE–0 CE (D = 5.98), with settlement Z = −5.36.

This secondary peak operates through a fundamentally different mechanism than the primary peak. The 3000–2500 BCE signal is monument-enrichment driven: monument Z rises sharply while settlement Z remains near zero. The 500 BCE signal is settlement-depletion driven: settlement Z falls to −5.36 while monument Z shows no comparable spike. This means settlements actively avoid the corridor in the Persian and early Hellenistic period, rather than monuments being drawn to it. This is consistent with the DARE finding (Section 4.2): the circle broadly avoids the densely settled Roman and Hellenistic heartland, and the 500 BCE settlement depletion likely reflects the circle's passage through regions peripheral to, rather than central to, the emerging Mediterranean urban economy. The two divergence peaks should not be interpreted as manifestations of the same underlying phenomenon.

The temporal bins at 3500–3000 BCE and 1500–1000 BCE had monument site counts below the minimum sites-per-bin threshold in some analysis runs; values reported from those bins used a lower threshold and should be treated as approximate. Full temporal decomposition results are presented in Table 3.

### 4.4 250-year resolution and Memphis identification

At 250-year resolution, the onset of the divergence is even more sharply defined. Monument Z jumps from −0.79 in the 3000–2750 BCE bin to 11.26 in the 2750–2500 BCE bin — a shift of more than 12 standard deviations in a single 250-year step. Thirty-one of 38 on-corridor monuments in the Pleiades ancient dataset fall within this single bin. All 38 monuments are Egyptian Old Kingdom pyramids in the Memphis necropolis: specifically the pyramid complexes at Abu Rawash, Giza, Abusir, Saqqara, and Dahshur (Lehner 1997). These 38 sites span just 0.44° of latitude and 0.55° of longitude — a geographic footprint of approximately 70 km. No monuments from the Levant, Iran, Mesopotamia, or any other region contribute to the divergence signal at its peak. The divergence at its peak is not a distributed corridor phenomenon; it is a single archaeological complex.

The collapse of the divergence is equally sharp. D falls from 4.39 in the 2250–2000 BCE bin to 0.76 in the 2000–1750 BCE bin, coinciding with the end of major pyramid construction in the Fourth Dynasty and the subsequent dispersal of royal necropolises southward to Dahshur, Lisht, and Hawara. (All 38 monuments carry Pleiades minDate values of −2670, placing them in the 2750–2500 BCE bin under the database's own temporal encoding. This corresponds to the onset of monumental pyramid construction beginning with the Third Dynasty Step Pyramid of Djoser at Saqqara, c. 2680 BCE, and continuing through the major Fourth Dynasty complexes at Giza, c. 2560–2500 BCE. The residual divergence in the 2500–2000 BCE bin at 500-year resolution (D = 2.87) is consistent with continuing Fourth Dynasty and Fifth Dynasty construction at Abusir and Saqqara.)

The geometric relationship between the great circle and the Memphis necropolis has a specific physical correlate. Ghoneim et al. (2024) identified the Ahramat Branch, a Holocene Nile distributary that flowed past the Abu Sir pyramid field during the third millennium BCE and likely served as a transport corridor for pyramid construction. The great circle crosses this ancient channel at approximately 70° (an obtuse angle approaching perpendicular) at the 20 km scale, with the closest approach to perpendicular in the Abu Sir region where the channel and circle alignment is most nearly orthogonal (~110° supplementary angle). Multi-scale analysis shows peak monument enrichment of 9.23× at 7 km corridor half-width, peak divergence D = 2.18 at 18 km, and a half-life of 32 km. The 50 km threshold used throughout this study is therefore conservative relative to the spatial scale of maximum signal.

Two-hundred-fifty-year resolution results are presented in Table 4.

### 4.5 Preservation test

The preservation test evaluated whether differential preservation of Egyptian monumental versus domestic architecture explains the divergence. Three conditions added progressively more hypothetical settlements: (1) baseline Pleiades Egypt data (D = 12.83); (2) 58 documented buried settlements from published surveys (Spencer 2024; Moeller 2016; Kemp 2018), giving D = 13.10; (3) 350 Monte Carlo simulated Nile floodplain settlements as a conservative upper bound (D = 10.21 ± 0.72, D > 2 in 100% of iterations). The documented buried settlements cluster in the Nile Delta (mean 243 km from the great circle), far from the desert plateau monuments. The preservation hypothesis does not account for the spatial pattern. Full trial details in Supplementary §S4.

Full preservation test results are presented in Table 5.

### 4.6 XRONOS replication and geographic dependence

XRONOS independently replicates the overall alignment: Z = 5.91 (land-constrained, 1,000 trials). For monument-settlement divergence, XRONOS yields D = +3.28 on the full global dataset, D = −14.05 when Egypt is removed from the analysis, and D = −13.98 when Africa is removed. The divergence is Egypt-dependent on XRONOS, precisely as it is on Pleiades. This is a consistent result, not a contradictory one: the Memphis concentration documented in Section 4.4 predicts exactly this pattern. If the divergence is driven by Egyptian Old Kingdom pyramid construction, it should disappear when Egypt is removed from any database that includes Egypt. Both Pleiades (Egypt-removed spatial block: D = −0.64) and XRONOS confirm this prediction.

### 4.7 Geographic decomposition and null results

Regional decomposition of the Pleiades divergence shows that the signal weakens eastward: Egypt and the Levant produce D = 8.68; Iran and Mesopotamia produce D = 2.60; South Asia produces D = 1.07 (inconclusive due to insufficient sites). The divergence is geographically graded, declining with distance from the Memphis necropolis. The South Asian result may also reflect a classification limitation: Kenoyer (personal communication) has noted that Indus Valley urban sites such as Mohenjo-daro feature massive city walls and substantial public architecture that constitute monumental construction, blurring the monument-settlement distinction that is sharper in the Egyptian and Levantine contexts where desert-plateau necropolises are clearly differentiated from alluvial settlement sites.

Two independent databases return null results for South America. The Peru Ministry of Culture database (1,733 monumental sites, 275 settlement sites) produces D = −1.25. The combined Pleiades and p3k14c South American subset produces D = 0.13. Both South American results are null or near-null. A fourth South American database (SAAID, 2,429 sites) was examined but computational outputs were not retained during analysis; its results are excluded from this report. The convergence of two verified null results supports the South American conclusion independently. There is no monument-settlement divergence in South America on either verified database.

The LuwianSiteAtlas, covering 483 professionally classified Bronze Age western Anatolia sites from approximately 2000–1200 BCE, produces a complete null: zero of 483 sites fall within 500 km of the great circle. This is a strong negative result for the western Anatolian Bronze Age, a period and region that might plausibly participate in any eastern Mediterranean corridor phenomenon.

Full regional and geographic decomposition results are presented in Tables 6 and 7.

### 4.8 Split-sample and spatial block validation

Split-sample validation is the primary defense against post-hoc inflation. Across 100 random 50/50 Megalithic Portal splits, with the Alison circle evaluated against the held-out half only, the mean Z-score is 9.45, the minimum across all 100 splits is 7.31, and 100% of splits exceed Z = 5. The signal survives blinded validation on every split. This result is the most important validation in the paper: it demonstrates that the alignment is not dependent on the specific sites used to identify it.

Post-hoc fitting contributes non-trivially to the raw Z-score. Among 200 random great circles fitted to 15-site subsets of Portal data and then evaluated against the full database, 24% exceeded the Alison circle's Z-score. This means the Alison circle ranks at approximately the 76th percentile of post-hoc fitted circles — not special relative to other post-hoc solutions. The split-sample validation, not the raw Z-score, is therefore the appropriate measure of whether the alignment is real.

The leave-one-cluster-out spatial block test with Egypt removed returns overall enrichment Z = 4.51 and Pleiades monument-settlement divergence D = −0.64. The overall enrichment survives Egypt removal; the Pleiades divergence does not. This is consistent with and predicted by the Memphis concentration finding.

Hemisphere-level blocks confirm that overall enrichment exists independently in both hemispheres: Old World Z = 9.45, New World Z = 13.41, though the New World result reflects habitability rather than monument-specific clustering.

### 4.9 Deep-time enrichment

Early Holocene enrichment (10,500–8,500 BP) on the global p3k14c database is 4.52× above background (within 200 km of the circle) but does not replicate on regional Near Eastern data (NERD: 1.10×), confirming it reflects multi-regional geometry rather than a specific corridor. Full analysis in Supplementary §S1.

### 4.10 Falsification of auxiliary hypotheses

Two auxiliary hypotheses associated with the great circle in popular literature were tested and falsified. First, the pre-Younger Dryas occupation hypothesis (Hancock 1995): in the merged 94,181-date database (p3k14c + XRONOS + ROAD v32), the on-corridor date count before 12,800 BP is 10 dates against a random circle mean of 127 (Z = −0.47, 33rd percentile). On-corridor continuity across the Younger Dryas boundary is 2.1% of sites versus 5.6% off-corridor — less continuity on the corridor than off. The corridor shows no preferential pre-glacial activity.

Second, the Giza longitude grid hypothesis (Hancock 1998) predicts that ancient sites should align at 36° longitude multiples from Giza (31.13°E). Tested across six sub-tests on 508,000 sites, all six fail. Giza ranks at the 71st–79th percentile among 1,000 random reference longitudes. Precessional angles produce Z = 0.45–0.94 (p = 0.20–0.27). Settlements are more grid-aligned than monuments (15.4% vs. 9.3%), the opposite of the prediction. The 108° angular grid hypothesis produces Z = −8.01 (also falsified).

### 4.11 Astronomical alignments and additional null results

Astronomical alignment tests (ecliptic plane, galactic plane, pole declination, solstice azimuths, equinox azimuths, and stellar rising and setting azimuths at precession-corrected dates) return null results after Bonferroni correction across all sub-tests. Monument orientation analysis on Megalithic Portal sites returns Rayleigh p = 0.75, indicating no directional clustering. Notably, all Egyptian Old Kingdom monuments relevant to this study are oriented north-south, perpendicular to the great circle's approximately east-west bearing at that latitude. The circle is geometrically orthogonal to the monuments it most strongly predicts.

Trade network correlation: radiocarbon date density within the Egypt–Iran corridor versus divergence D across 500-year bins yields r = 0.04, p = 0.93 on Pleiades — no evidence of a trade-driven mechanism. Desert edge hypothesis: the mean distance of circle sample points in the Egyptian window from the western escarpment is 134 km; the circle crosses the Nile Valley diagonally and does not track the desert edge. Geophysical scan across 13 tests (geoid gradient, ocean fraction, ocean current alignment, and others): geoid gradient Z = +2.33, ocean fraction Z = −2.05, ocean current alignment Z = +2.08, none surviving Bonferroni correction. The pattern is consistent with continental margins and navigable coastal routes, but no geophysical variable explains it.

Agricultural suitability correlation: p = 0.57. Lithology: p = 0.74. Groundwater: p = 0.37. These analyses are reported fully in the companion study (Allan, under review). LGM connectivity was assessed; LGM-connected and ocean-separated clusters showed comparable enrichment, consistent with the overall pattern, but the specific Z-scores could not be verified against computed output files during pre-publication review and are therefore excluded from this report.

### 4.12 Circle family and geometric uniqueness

The Alison circle belongs to a family of approximately 66 great circles connecting Giza, Nazca, and Easter Island (~1 in 1,515 random circles). Within this family, the Alison circle ranks 5th of 60 by monument-settlement divergence (D = 8.64). Full analysis in Supplementary §S2.

### 4.13 Predictive validation

Pre-2001 p3k14c data (using the reference publication year as a proxy for data availability before the great circle hypothesis was widely discussed) returns Z = 0.06 — effectively null. Post-2001 p3k14c data returns Z = 16.74. The signal is substantially stronger in newer data, which is the opposite of what a post-hoc selection artifact would predict — selection artifacts should show stronger signals in older data used to inspire the hypothesis.

However, the post-2001 signal has a geographic concentration problem. Seventy-six point three percent of near-corridor (within 50 km) post-2001 dated sites fall in five grid cells (Easter Island, Egypt, the Levant, Peru, and Iran), which constitute 0.85% of the total data grid. Removing those five cells eliminates the post-2001 signal (Z = −3.68). The alignment is confirmed in its established geographic regions; it does not predict new regions not already established by earlier fieldwork. The predictive validation is positive in the sense that new data strengthens the signal, but the geographic concentration means the validation does not constitute independent geographic prediction.

### 4.14 Paleoclimate correlation

A suggestive correlation between divergence and Soreq Cave δ18O at 500-year lag (r = −0.77, p = 0.025, n = 8 bins) does not survive Benjamini-Hochberg correction (adjusted p ≈ 0.075). Details in Supplementary §S3.

### 4.15 Multiple testing correction

Benjamini-Hochberg FDR correction at q = 0.05 was applied across all primary statistical tests. All primary findings survive: monument-settlement divergence on multiple databases with different classification methods, habitability-adjusted divergence at the 99.63rd percentile, split-sample validation across 100 Portal splits, and three-database overall enrichment. The paleoclimate correlation does not survive (adjusted p ≈ 0.075).

---

## 5. Discussion

### 5.1 What is explained: the habitability corridor

The habitability-adjusted result (78.5th percentile; not significant) establishes that the overall site-count enrichment reflects the circle's ecological logic — its path through the Nile Valley, the Levantine coast, Mesopotamia, the Indus floodplain, coastal Peru, and Easter Island. That a great circle can connect Giza, Nazca, and Easter Island while passing through these densely occupied corridors is a geometric property of the Earth's surface worth documenting — roughly 1 in 1,515 random circles achieve the three-anchor constraint — but it does not require any causal connection between the societies at those sites.

The land-constrained correction implemented here constitutes a methodological contribution beyond this specific question. Standard distribution-matched Monte Carlo baselines, widely used in archaeological corridor and catchment analysis, do not account for the fact that coastal and island sites preferentially produce ocean-landing jitter attempts, deflating expected counts and inflating Z-scores. The magnitude of the inflation was not trivial: Z-scores were reduced by 31% to 73% after correction across the four databases tested. We encourage adoption of land-constrained baselines as a standard procedure in any enrichment test involving sites near coastlines or on islands.

### 5.2 What is not explained: monument-settlement divergence

The monument-settlement divergence is the finding that habitability does not explain. Settlements should cluster at least as strongly as monuments if the pattern were purely geographic; instead, settlements are depleted (Pleiades settlement Z = −2.91) while monuments enrich (monument Z = 6.74), yielding D = 9.65. No random great-circle pole among 10,000 tested produces a comparable divergence (maximum D = 0.90).

The divergence survives the strongest validation tests applied. Split-sample analysis confirms that the signal is not in the sites used to inspire the hypothesis: 100% of 100 held-out Portal splits produce Z above 5. Preservation modeling with 350 added hypothetical buried settlements does not eliminate the divergence (D = 10.21 ± 0.72, D > 2 in all 1,000 iterations). Habitability-adjusted ranking puts the divergence at the 99.63rd percentile. The alternative explanations tested in this study — geography, habitability, preservation, research intensity, sampling bias, trade routes, astronomical encoding, desert edge position, 108° grids, and geophysical variables — do not individually account for the divergence. A companion study (Allan, under review) identifies tectonic diversity and geographic circumscription as partial geological explanations for the corridor's association with civilizational centers, though monument-specific clustering relative to settlements is not fully resolved by those variables. The divergence is confirmed on multiple databases with different classification methods: one geographically and typologically independent (Megalithic Portal non-EU, D(LC) = 4.56), one typologically independent but geographically overlapping (DARE, D(LC) = 5.94), and one using uncorrected baselines with a different classification method (p3k14c, D = 6.18).

### 5.3 The Memphis question

It is important to state clearly what the geographic structure of the divergence is and is not. The monument-settlement divergence at its 2750–2500 BCE peak is a localized spatial cluster — the Memphis necropolis, analogous to a local indicator of spatial association in the LISA framework (Anselin 1995) — not a distributed phenomenon along the circle's full length. The great circle provides the geometric framework in which the divergence was detected: it is the alignment test that revealed the monument-settlement contrast, and the divergence replicates across databases that span different portions of the circle's path. But the divergence itself, at maximum intensity, is concentrated in a ~70 km footprint. The appropriate interpretation is that the great circle identifies a transcontinental corridor in which a remarkable local concentration of monuments — the densest in the ancient world — sits within a narrow band, not that the divergence operates uniformly along the corridor. Spatial autocorrelation at the Memphis scale is expected under Tobler's first law (Tobler 1970); the question is why the Memphis cluster is so close to the circle.

The divergence at its maximum is not a corridor phenomenon. The divergence consists exclusively of 38 Egyptian Old Kingdom pyramids in the Memphis necropolis, located 18 to 44 km from a segment of the tested geodesic within a ~70 km geographic footprint. This arrangement appeared during a single 250-year construction window (the Dynasty IV construction peak, which straddles the 2750–2500 and 2500–2250 BCE bins in our temporal decomposition) and collapsed as sharply as it arose. No other region — not the Levant, not Mesopotamia, not the Indus Valley, not South America — contributes to the divergence peak. The divergence in the 2750–2500 BCE bin is, in geographic terms, the pyramid fields of Abu Rawash, Giza, Abusir, Saqqara, and Dahshur.

Several aspects of this arrangement are noteworthy and do not point to obvious explanations; in particular, they do not imply any shared planning or communication between geographically distant cultures. The great circle's bearing at the Memphis latitude is approximately east-west. The pyramid complexes are oriented north-south. The circle is geometrically perpendicular to the monuments it predicts most strongly, ruling out any simple compass-bearing alignment. The relevant landscape axis runs north-south along the desert plateau above the Holocene Nile floodplain — the same axis documented by Ghoneim et al. (2024) as the Ahramat Branch distributary canal network. The pyramid fields sit on the plateau above this ancient waterway. It is plausible that the north-south organization of the plateau edge, combined with the specific geography of the Nile valley at this location, places the monuments near a great circle that happens to connect distant sites — but this reasoning, while intuitive, does not amount to an explanation, because it does not account for why this specific great circle (connecting Nazca and Easter Island in the New World) has the property of intersecting the Memphis plateau rather than any of the countless great circles that do not.

The collapse of the divergence by 2000 BCE is consistent with the historical record: royal necropolises moved south after the Fourth Dynasty, and Nile hydrology underwent substantial changes around this period (Sheisha et al. 2022). The divergence is time-bounded at both ends in a manner consistent with the construction history of the Memphis necropolis specifically.

### 5.4 The temporal origin: a post-glacial phenomenon

The great circle corridor was not preferentially occupied before the Younger Dryas. The pre-12,800 BP on-corridor date count (10 dates versus expected 127, Z = −0.47) is a clean null. The corridor shows less Younger Dryas boundary continuity than off-corridor sites (2.1% versus 5.6%). The early Holocene enrichment at the global scale (4.52× at 10,500–8,500 BP) disappears at the regional scale (1.10× on NERD), confirming that it reflects geometric coincidence rather than specific corridor behavior. The archaeological signal along this great circle is entirely post-glacial, and its monument-specific peak occurs approximately 5,000 years after initial post-glacial habitation of the relevant regions.

Paper 2's preliminary framing of the pattern as a potentially "inherited corridor" requires revision in light of these findings. The evidence does not support the interpretation that human populations tracked a pre-existing environmental corridor along this geodesic. The NERD null demonstrates that the early Holocene enrichment is a multi-regional geometric artifact. Normal PPNB expansion from the Fertile Crescent, tracking agricultural dispersal routes that are perpendicular to the circle's bearing at that latitude, accounts for the early Holocene date distribution. The monument-specific onset, 5,000 years after initial habitation, argues against any continuous corridor-following explanation and in favor of an explanation rooted in the specific geography of the early Bronze Age Nile Valley.

### 5.5 What this study falsifies

The study falsifies several specific claims with quantitative rigor. The longitude grid hypothesis fails all six sub-tests: Giza is not a special reference longitude (71st–79th percentile among 1,000 random longitudes), precessional angles are not significant (p = 0.20–0.27), settlements are more grid-aligned than monuments (15.4% versus 9.3%, opposite of the prediction), and the oldest sites show less grid alignment than the youngest (opposite of ancient encoding). The 108° angular separation hypothesis produces Z = −8.01. The pre-Younger Dryas civilization hypothesis fails: the corridor was empty before 12,800 BP, and on-corridor sites show less YD boundary continuity than off-corridor sites. Astronomical encoding fails Bonferroni correction across all sub-tests. The specific empirical claims in these hypotheses are not supported by the data; the same methodology that identifies the great circle alignment falsifies these auxiliary claims.

### 5.6 Geological mechanisms (companion study)

A companion study (Allan, under review) investigates whether geological or environmental variables explain why major civilizational centers cluster near this geodesic. Tectonic diversity along the corridor is higher than expected by chance (p = 0.011). Monument-building regions along the circle show higher rates of geographic circumscription — environmental constraint promoting social complexity — than random corridors (37.2% versus 10.9%, p = 0.014). Agricultural suitability, lithology, and groundwater variables do not explain the pattern. These geological analyses offer a partial explanation for the corridor's association with civilizational centers at the broadest scale, but they do not explain why monuments specifically cluster relative to settlements, and they are reported in a separate study to preserve statistical independence between the spatial and geological analyses. The geological analysis in the companion study operates at a different spatial scale from the monument-settlement divergence: tectonic diversity and circumscription explain why independently originated civilizations cluster along the corridor at the continental scale, while the monument-settlement divergence documented here operates at the intra-regional scale within the Egyptian segment of that corridor. The two findings are complementary, not redundant.

### 5.7 Limitations

Seven limitations are identified and reported.

First, the post-hoc circle definition cannot be fully resolved by validation. The great circle was identified by observing famous sites, not optimized against a complete database. Split-sample validation substantially mitigates this concern, but the possibility that the circle's fame creates subtle database biases — such as professional excavation being more likely near famous-site regions — cannot be eliminated.

Second, land-constrained divergence was not successfully computed for all databases. LC divergence was computed for the Megalithic Portal non-EU subset (D(LC)=4.56) and DARE (D(LC)=5.94). For p3k14c, the keyword-matched monument and settlement classification yielded a degenerate result (zero classified sites within 50 km of the great circle); the uncorrected D=6.18 from the predecessor analysis is retained but should be treated with caution as it rests on a different classification approach. We recommend land-constrained baselines as standard practice in any corridor enrichment test involving coastal or island archaeological sites (see Section 5.1).

Third, the Pleiades monument-settlement divergence is Egypt-dependent. Removing Egypt from the Pleiades spatial block eliminates D (D = −0.64 without Egypt). The divergence finding on this database is geographically concentrated.

Fourth, p3k14c after land-constrained correction is marginal (Z = 2.09, p = 0.02). This is directionally consistent with the other databases but should be interpreted cautiously and not used as a primary finding in isolation.

Fifth, sub-Saharan Africa, Southeast Asia, and the Pacific are underrepresented in most databases used here. Results from those regions — positive or negative — should not be assumed to generalize.

Sixth, the Megalithic Portal's classification system is broad. The "Ancient Village or Settlement" category likely contains sites that a professionally classified database would label as monumental. Monument-settlement divergence on the Portal is affected by this imprecision in an unknown direction.

Seventh, the paleoclimate correlation (r = −0.77 at 500-year lag, n = 8 bins) is marginal after multiple testing correction (adjusted p ≈ 0.075). The small number of temporal bins makes this result fragile; a different binning scheme could alter it.

---

## 6. Conclusion

Ancient archaeological sites cluster near the great circle defined by Alison (c. 2001) at statistically significant rates on three of four land-constrained databases (Portal Z = 8.77, Pleiades Z = 6.74, XRONOS Z = 5.91) and marginally on a fourth (p3k14c Z = 2.09, p = 0.02). This overall site-count enrichment is fully explained by the circle's path through the world's most densely occupied prehistoric corridors: after habitability adjustment, the enrichment falls to the 78.5th percentile. The overall result is geographic, not mysterious.

The robust finding is the monument-settlement divergence. Monuments cluster along the circle (Pleiades monument Z = 6.74) while contemporaneous settlements in the same regions at the same times are depleted (settlement Z = −2.91), producing D = 9.65 (land-constrained). This divergence is not matched by any of 10,000 random great-circle poles (maximum D = 0.90), ranks at the 99.63rd percentile of 4,323 habitability-matched circles, is confirmed on multiple databases with different classification methods, and survives all alternative explanations tested, including differential preservation, research intensity, geographic decomposition, and more than ten specific null hypothesis tests.

The divergence is temporally sharp (onset 2750–2500 BCE), geographically concentrated (Memphis necropolis), and does not appear in South America (null on two verified databases), western Anatolia (0/483 sites), or across the Younger Dryas boundary (Z = −0.47, 33rd percentile). Astronomical encoding, longitude grids, and pre-glacial activity are falsified by the same methodology that confirms the great circle signal. The corridor was post-glacially colonized following normal PPNB agricultural expansion, and the monument-specific onset appears 5,000 years after initial habitation.

The monument-settlement divergence is partially addressed by geological factors identified in a companion study (Allan, under review): the corridor concentrates tectonic diversity (p = 0.011) and desert-river circumscription (p = 0.014) at rates significantly above random great circles, providing a partial explanation for why independently originated civilizations cluster along this path. At the Memphis necropolis specifically, the temporal drift of pyramid placement provides a landscape-level mechanism: pyramids built before the 4.2-kiloyear aridification event (deMenocal 2001; Staubwasser et al. 2003; Weiss et al. 1993) average 8.4 km from the great circle, while those built after the desiccation of the Ahramat Branch average 19.7 km — consistent with the ancient waterway, not the great circle per se, governing early pyramid placement. Whether the divergence is fully explained by the combination of tectonic-climatic corridor selection and local Nile hydrology, or whether additional factors contribute, remains an open empirical question.

Future directions include geophysical surveys of the Abu Sir region, formal point-process modeling at fine spatial scales (Allan, under review), ancient DNA analysis of on-corridor versus off-corridor populations in the 3000–2500 BCE window, and acquisition of professionally classified databases for sub-Saharan Africa and Southeast Asia. All data and code are publicly available at github.com/thegreatcircledata/great-circle-analysis.

---

## Declaration of Competing Interests

The author declares no competing financial interests.

## Funding Statement

This research received no specific funding from any agency in the public, commercial, or not-for-profit sectors.

## Ethics Statement

Not applicable. This study uses published, publicly accessible archaeological databases and involves no human subjects, archaeological fieldwork, or collection of new materials.

---

## References

Alison, J. (c. 2001). Great circle alignment of ancient sites. Retrieved from original publication.

Anselin, L. (1995). Local indicators of spatial association — LISA. *Geographical Analysis*, 27(2), 93–115.

Allan, E. (2026a). Statistical analysis of ancient monumental site distribution along a proposed great circle. Preprint. doi.org/10.5281/zenodo.19212669

Allan, E. (under review). Collinearity of independently originated civilizations along a single great circle. Preprint: doi.org/10.5281/zenodo.19240285.

Baddeley, A., Rubak, E., & Turner, R. (2015). *Spatial Point Patterns: Methodology and Applications with R*. Chapman and Hall/CRC.

Bar-Matthews, M., Ayalon, A., Gilmour, M., Matthews, A., & Hawkesworth, C. J. (2003). Sea–land oxygen isotopic relationships from planktonic foraminifera and speleothems in the Eastern Mediterranean region and their implication for paleorainfall during interglacial intervals. *Geochimica et Cosmochimica Acta*, 67(17), 3181–3199.

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300.

Bird, M. I., et al. (2022). p3k14c, a synthetic global database of archaeological radiocarbon dates. *Scientific Data*, 9, 27.

Borrell, F., Junno, A., & Barceló, J. A. (2015). Synchronous environmental and cultural change in the formation of the first global polity. *PLoS ONE*, 10(10), e0140421.

Bevan, A., & Conolly, J. (2006). Multiscalar approaches to settlement pattern analysis. In *Digital Archaeology: Bridging Method and Theory* (Bentley, R. A. and Maschner, H. D. G., Eds.), 217–234. Routledge.

Conolly, J., & Lake, M. (2006). *Geographical Information Systems in Archaeology*. Cambridge University Press.

Crema, E. R., Bevan, A., & Lake, M. (2010). A probabilistic framework for assessing spatio-temporal point patterns in the archaeological record. *Journal of Archaeological Science*, 37(5), 1118–1130.

deMenocal, P. B. (2001). Cultural responses to climate change during the late Holocene. *Science*, 292(5517), 667–673.

Diggle, P. J. (2003). *Statistical Analysis of Spatial Point Patterns* (2nd ed.). Arnold.

Elliott, T., & Gillies, S. (2009–2025). Pleiades: A gazetteer of past places. Ancient World Mapping Center and Institute for the Study of the Ancient World. https://pleiades.stoa.org

Fortin, M. J., & Dale, M. R. T. (2005). *Spatial Analysis: A Guide for Ecologists*. Cambridge University Press.

Ghoneim, E., Ralph, T. J., Onstine, S., Arnold, C., Manocha, S., Hawsawi, A. M., & Khedr, F. (2024). The Egyptian pyramid chain was built along the now abandoned Ahramat Branch of the Nile. *Communications Earth & Environment*, 5, 233.

Hancock, G. (1995). *Fingerprints of the Gods*. Crown.

Hancock, G. (1998). *Heaven's Mirror*. Crown.

Hodder, I., & Orton, C. (1976). *Spatial Analysis in Archaeology*. Cambridge University Press.

Hewitt, R. J., Wenban-Smith, F. F., & Bates, M. R. (2020). Detecting associations between archaeological site distributions and landscape features: A Monte Carlo simulation approach for the R environment. *Geosciences*, 10(9), 326.

Kemp, B. (2018). *Ancient Egypt: Anatomy of a Civilization* (3rd ed.). Routledge.

Klein Goldewijk, K., Beusen, A., Doelman, J., & Stehfest, E. (2017). Anthropogenic land use estimates for the Holocene — HYDE 3.2. *Earth System Science Data*, 9, 927–953.

Lehner, M. (1997). *The Complete Pyramids*. Thames & Hudson.

Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. *Cancer Research*, 27(2), 209–220.

Moeller, N. (2016). *The Archaeology of Urbanism in Ancient Egypt*. Cambridge University Press.

Moreno-Meynard, P., et al. (2022). Past human mobility corridors and least-cost path Monte Carlo simulation. *Land*, 11(5), 672.

Palmisano, A., et al. (2025). XRONOS: A global radiocarbon database. *Journal of Archaeological Science*.

Ripley, B. D. (1976). The second-order analysis of stationary point processes. *Journal of Applied Probability*, 13(2), 255–266.

Roberts, D. R., et al. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*, 40(8), 913–929.

Scott, D. W. (2015). *Multivariate Density Estimation: Theory, Practice, and Visualization* (2nd ed.). John Wiley & Sons.

Sheisha, H., et al. (2022). Nile waterscapes facilitated the construction of the Giza pyramids during the third millennium BCE. *PNAS*, 119(37), e2202530119.

Spencer, A. J. (Ed.). (2024). *The Delta Survey 1997–2020* (EES Occasional Publications). Egypt Exploration Society.

Staubwasser, M., et al. (2003). Climate change at the 4.2 ka BP termination of the Indus valley civilization. *Geophysical Research Letters*, 30(8).

Talbert, R. J. A. (Ed.). (2000). *Barrington Atlas of the Greek and Roman World*. Princeton University Press.

Thomas, M. (1949). A generalization of Poisson's binomial limit for use in ecology. *Biometrika*, 36, 18–25.

Tobler, W. R. (1970). A computer movie simulating urban growth in the Detroit region. *Economic Geography*, 46(sup1), 234–240.

Vermeersch, P. M. (2024). ROAD v32: Radiocarbon dates from the Old World. *Journal of Open Archaeology Data*.

Weiss, H., Courty, M. A., Wetterstrom, W., et al. (1993). The genesis and collapse of third millennium North Mesopotamian civilization. *Science*, 261(5124), 995–1004.

Williamson, T., & Bellamy, L. (1983). *Ley Lines in Question*. World's Work.

---

## Data Availability

All analysis code and supporting data: github.com/thegreatcircledata/great-circle-analysis

Megalithic Portal: megalithic.co.uk (terms of use apply)
Pleiades: pleiades.stoa.org (CC BY 3.0)
p3k14c: Bird et al. (2022), archived in tDAR
XRONOS: xronos.ch
DARE: dare.ht.lu.se (CC BY-SA)
Peru Ministry of Culture: geoportal.cultura.gob.pe
Historic England: historicengland.org.uk (Open Government Licence)
SAAID: Pandora/PACHAMAMA platform
Interactive visualization: thegreatcircle.earth
Preprints: Earlier versions of this work were deposited as preprints: v1–v4 at doi.org/10.5281/zenodo.19212669. Companion study: Allan (under review), preprint at doi.org/10.5281/zenodo.19240285.

Database access terms: Pleiades (CC BY 3.0); DARE (CC BY-SA 4.0); Historic England (Open Government Licence v3.0); p3k14c (archived in tDAR, open access per Bird et al. 2022); XRONOS (open access); Peru Ministry of Culture GeoPortal (publicly accessible, terms of use apply); Megalithic Portal (publicly accessible, used with acknowledgment per site terms); SAAID (excluded from analysis). CNSA/IPHAN data were obtained from a publicly available Kaggle repository (CC BY 4.0).

Supplementary Information containing extended analyses (deep-time enrichment, circle family analysis, paleoclimate correlation, and preservation test details) is provided as a separate file.

Verification output files for all reported statistics are included in the data repository under `outputs/merged_paper_prep/`.

---

## Acknowledgments

I thank all contributors to and maintainers of the open databases used in this analysis — the Megalithic Portal volunteer community, the Pleiades editorial team, the p3k14c compilation team (Bird et al.), the XRONOS team (Palmisano et al.), Johan Ahlfeldt and the DARE team at Lund University, the Peru Ministry of Culture GeoPortal team, and Historic England. Computational analysis was conducted with assistance from Claude (Anthropic) and Claude Code. The AI tools were used for: code generation and iterative development of statistical simulations; data processing and database integration; Monte Carlo baseline computation including land-constrained resampling; and drafting and revision of manuscript text. All hypotheses, analytical decisions, data interpretation, and conclusions are my own. I reviewed and verified all AI-generated outputs against expected statistical behavior, validated results through cross-database replication, and take full responsibility for the content of this publication. The coastal bias correction was prompted by anonymous feedback identifying the coastal jitter problem. All code is openly available for independent verification.

---

## AI Disclosure

Per applicable journal AI policy: generative AI tools (Claude, Anthropic) were used for computational assistance, code generation, and manuscript drafting as described in the Acknowledgments. The AI tools were not used as intellectual contributors; all scientific judgments are my own.

---

## Tables

*Tables 1–10 follow. For submission, tables may be provided as separate files per journal requirements.*

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

**Notes on site counts:** ~259,000 unique sites total across all 8 databases; over 550,000 raw entries, dates, and records (550,000 = entries/records, NOT unique sites). Paper 1 (Allan 2026a) used 25,557 p3k14c sites (lat/lon dedup); this study uses 36,693 (SiteID dedup, more complete); same 187 on-corridor sites in both.

**Habitability-adjusted result:** 78.5th percentile of 4,323 habitability-matched circles. Overall site-count enrichment is not significant after habitability adjustment.

**LC Z-score reductions from correction:** Portal 13.2→8.77 (33.6%); Pleiades 11.14→6.74 (39.5%); XRONOS 12.4→5.91 (52.3%); p3k14c 7.77→2.09 (73.1%).

---

## Table 2. Monument-settlement divergence across databases

D = Z_monument − Z_settlement. Primary LC result (Pleiades) shown in bold. LC baselines computed for Portal, p3k14c, and DARE (1,000 trials each; see lc_divergence_all_databases.json).

| Database | Monument Z | Settlement Z | D | Correction | N monuments | N settlements | Notes |
|---|---|---|---|---|---|---|---|
| Pleiades (pre-2000 BCE, Egypt/Levant) | **6.74** | **−2.91** | **9.65** | LC | 406 | ~372 | Primary finding; 0/10,000 random great-circle poles match; max random D = 0.90 |
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
**Divergence uniqueness:** 0/10,000 random great-circle poles produce D > 9.65 (LC) on Pleiades ancient sites. Random distribution: mean = 0.00, std = 0.07, max = 0.90.
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
| **2750–2500** | **11.26** | — | **peak** | **31 of 38** | **Third Dynasty onset (Djoser c. 2680 BCE) through early Dynasty IV; all 38 monuments have Pleiades minDate −2670; Memphis necropolis** |
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
| Agricultural suitability | p = 0.57 | Not significant | NULL | See Allan (under review) |
| Lithology | p = 0.74 | Not significant | NULL | See Allan (under review) |
| Groundwater | p = 0.37 | Not significant | NULL | See Allan (under review) |
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
