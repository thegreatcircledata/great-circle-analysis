# TEMPORAL PROPAGATION & WAVEFRONT ANALYSIS — Research Directive v1.0

**Date:** 2026-03-23
**Author:** Claude (for Ell)
**Objective:** Determine whether monument construction propagated sequentially along the Great Circle (implying knowledge transfer or migration) or emerged simultaneously at independent points (implying shared response to a common stimulus or coincidence).
**Repository:** great-circle-analysis/
**Output Directory:** outputs/temporal_propagation/
**Estimated Runtime:** 1–2 hours
**Parallel Safe:** Yes

---

## Background & Motivation

The temporal decomposition (paper2-v3) established that the monument-settlement divergence peaks at 2750–2500 BCE, concentrated in the Memphis necropolis. But the analysis binned temporally, not spatially-temporally. We haven't asked: **did monumentality originate at one point on the circle and spread along it?** This is a fundamentally different question from "when was the peak?" — it asks about *directionality* and *causality*.

If monuments appear simultaneously at multiple widely separated points → independent emergence (the corridor is a latent attractor, activated independently).
If they propagate sequentially → knowledge transfer, migration, or cultural diffusion along the corridor.

The Indus synchrony is particularly tantalizing: Mohenjo-daro's urban/monumental phase (2600–2500 BCE) overlaps almost exactly with the Memphis pyramid-building peak.

---

## Analysis 1: Space-Time Wavefront Plot

### Data
- **All dated monumental sites within 50km of the Great Circle** — compile from:
  - Pleiades (with start_date field)
  - Megalithic Portal (with type-based age estimates from age_analysis)
  - p3k14c (radiocarbon dates with lat/lon)
  - Peru Ministry of Culture (with period assignments)
- **Key fields needed:** lat, lon, estimated_date_BCE, site_name, type/category, database_source

### Method
1. For each dated site within 50km of the circle:
   a. Compute its position along the circle in degrees (0–360°, starting from an arbitrary reference point — suggest starting at Giza = 0°)
   b. Record its estimated date (use midpoint of date range if a range)
2. **Create the wavefront plot:** X-axis = position along circle (degrees), Y-axis = date (BCE, oldest at top)
   - Each point = one dated site
   - Color by cluster region (Egypt, Iran, Indus, Peru, Easter Island, SE Asia)
3. **Visual inspection:** Does the plot show:
   a. A horizontal line (simultaneous emergence at multiple points) → independent origin
   b. A diagonal line (propagation in one direction) → diffusion
   c. A V-shape or expanding wavefront (originating from one point) → radiation from a center
   d. No pattern (scattered) → no temporal-spatial structure

4. **Formal test — Spearman rank correlation:**
   - Compute Spearman's ρ between position-along-circle and construction date
   - If ρ is significantly positive or negative → evidence of directional propagation
   - Monte Carlo: shuffle the dates among the sites 10,000 times, compute ρ each time → p-value

5. **Regional synchrony test:**
   - For each pair of cluster regions (Egypt-Iran, Egypt-Indus, Egypt-Peru, Iran-Indus, etc.):
     a. Compute the mean construction date for monumental sites in each cluster
     b. Compute the temporal overlap between the two clusters' active periods
     c. Test: given the geographic separation, is the temporal overlap expected under a diffusion model? (Use estimated diffusion rates from known cultural expansions: Neolithic ~1 km/year, Bronze Age trade ~5 km/year)
     d. If overlap is too large for diffusion to explain → independent emergence

### Output
- `wavefront_plot.png` — the space-time visualization (this will be visually striking either way)
- `spearman_correlation.json` — ρ, p-value, interpretation
- `cluster_synchrony_matrix.json` — pairwise temporal overlap between all cluster pairs, diffusion-expected vs. observed
- `dated_sites_along_circle.csv` — master dataset used

---

## Analysis 2: Egypt-Indus Synchrony Deep Dive

### Background
Memphis pyramid-building peak: 2750–2500 BCE (31 of 38 sites in the temporal bin).
Mohenjo-daro/Harappa urban phase: 2600–2300 BCE (Mature Harappan period).
These overlap almost exactly despite ~3,500 km separation along the circle.

### Data
- **Egyptian sites:** Use the dated pyramids and monuments from the Megalithic Portal + Pleiades
- **Indus sites:** Compile from Pleiades (filter to Indus region) + any Kenoyer datasets available
  - If no detailed Indus dataset: use the known chronology from Possehl 2002 and Kenoyer 1998 (compiled list of major Indus sites with founding dates)

### Method
1. Compile all dated monumental/urban sites in both clusters (within 100km of circle)
2. For each cluster, compute:
   a. Earliest monument date
   b. Peak monument construction rate (construction starts per century)
   c. Duration of monument-building phase
3. **Synchrony test:**
   a. Compute the temporal offset between the two peaks
   b. Monte Carlo: if both clusters drew their construction dates from a uniform distribution spanning 4000–1000 BCE, how often would the peaks fall within X years of each other?
   c. Report: is the observed synchrony statistically significant?
4. **Diffusion feasibility:**
   a. Distance along the circle between Memphis and Mohenjo-daro: compute in km
   b. If diffusion rate = 1 km/year (Neolithic expansion rate): expected delay = distance/rate
   c. If diffusion rate = 5 km/year (Bronze Age trade): expected delay = distance/rate
   d. Compare to observed delay (which appears to be ~100–150 years at most)

### Output
- `egypt_indus_synchrony.json` — peak dates, temporal offset, Monte Carlo p-value
- `diffusion_feasibility.json` — expected delay vs. observed delay at different diffusion rates
- `dual_construction_timeline.png` — overlaid construction rate curves for both clusters

---

## Analysis 3: Early Holocene → Bronze Age Transition

### Background
The early Holocene enrichment (8500–6500 BCE) predates monuments by ~4,000–5,000 years. Is there continuity? Did the same locations that were active in the early Holocene become monumental in the Bronze Age?

### Method
1. Identify all early Holocene sites on the corridor (from p3k14c, 10500–6500 BCE, within 50km)
2. Identify all Bronze Age monumental sites on the corridor (3000–2000 BCE, within 50km)
3. For each Bronze Age monumental site, compute the distance to the nearest early Holocene site
4. **Monte Carlo:** For 10,000 random shuffles of the early Holocene site positions (keeping them on the corridor but randomizing their positions along it), compute the mean distance to Bronze Age sites
5. If Bronze Age monuments preferentially sit near early Holocene occupation sites → continuity (the corridor was used consistently for 6,000+ years)
6. If no spatial correlation between the two epochs → the corridor was "rediscovered" rather than continuously occupied

### Output
- `holocene_bronze_age_continuity.json` — mean distance to nearest precursor, Monte Carlo percentile
- `continuity_map.png` — early Holocene sites (blue dots) and Bronze Age monuments (red dots) along the circle

---

## Lookahead / Bias Warnings
- Date estimates from the Megalithic Portal type-based system have enormous uncertainties (±1000 years in some cases). Use the caveats from age_analysis_results.json.
- Pleiades start_date is often "traditional" rather than radiometric. Prefer radiocarbon-dated sites where possible.
- The Egypt cluster has far more precisely dated sites than other regions — this creates an apparent Egypt-centrism that may be a data quality artifact, not a real pattern.
- Diffusion models assume continuous land routes — the Egypt-Indus corridor has major geographic barriers (Arabian desert). Account for this in feasibility estimates.

---

## Deliverables
1. `outputs/temporal_propagation/RESULTS.md` — narrative summary with clear conclusion
2. All CSVs, JSONs, and figures listed above
3. The wavefront plot alone could be a standalone Substack graphic regardless of result
