# EARLY HOLOCENE CORRIDOR CHARACTERIZATION — Research Directive v1.0

**Date:** 2026-03-23
**Author:** Claude (for Ell)
**Objective:** Determine what kind of human activity was concentrated along the Great Circle during the early Holocene (10,500–6,500 BCE), whether it correlates with independent centers of agricultural origin, and whether genetically related populations cluster along the corridor.
**Repository:** great-circle-analysis/
**Output Directory:** outputs/early_holocene_corridor/
**Estimated Runtime:** 2–4 hours
**Parallel Safe:** Yes

---

## Background & Motivation

Paper 2 found 4–6× enrichment of radiocarbon dates near the Great Circle during 8,500–6,500 BCE — thousands of years before any monuments were built. This is one of the most provocative findings in the entire project: the corridor was already "special" in the early Holocene, long before the 2750–2500 BCE monument spike. But we don't know *what* people were doing there. Were these campsites? Early agricultural settlements? Ritual sites? Understanding the *character* of early Holocene activity on the corridor could distinguish between competing explanations (migration route vs. resource corridor vs. something else entirely).

---

## Analysis 1: p3k14c Site Context Classification

### Data
- **Primary:** p3k14c radiocarbon database (170,150+ dates with site metadata)
- Download: https://doi.org/10.1038/s41597-022-01118-7 (or use existing local copy)
- **Fields needed:** lat, lon, age_BP, site_name, site_type (if available), material_dated, context/phase fields

### Method
1. Filter to dates between 12,500–8,500 BP (10,500–6,500 BCE)
2. Compute distance from each date's lat/lon to the Great Circle (pole: 59.682122°N, 138.646087°W)
3. Classify as "on-corridor" (≤50km) vs. "off-corridor" (>50km)
4. For on-corridor dates:
   a. Extract all available context fields (material_dated, site descriptions, archaeological phase)
   b. Classify into categories: **campsite/ephemeral**, **settlement/domestic**, **ritual/funerary**, **agricultural**, **unknown**
   c. Use material_dated as proxy: charcoal/wood = campfire/occupation; grain/seed = agriculture; bone = hunting/domestic; shell = coastal foraging
5. Compare the category distribution of on-corridor vs. off-corridor dates using chi-square test
6. Specifically test: is the on-corridor population enriched for agricultural-indicator materials (grain, seed, domesticated animal bone) relative to the global baseline?

### Output
- `early_holocene_site_types.csv` — every on-corridor date with classification
- `category_comparison.json` — chi-square test results, category distributions on/off corridor
- `temporal_profile.png` — 250-year bin histogram of on-corridor dates by category

### Success Criteria
- If on-corridor dates are significantly enriched for agricultural indicators → supports "corridor as agricultural highway" hypothesis
- If enriched for ephemeral/campsite materials → supports "migration route" hypothesis
- If no category difference → the corridor enrichment is type-agnostic (activity-neutral)

---

## Analysis 2: Agricultural Origin Center Overlay

### Data
- **ArchaeoGLOBE dataset** (Stephens et al. 2019, doi:10.1126/science.aax1192) — global land-use reconstruction at 10,000, 8,000, 6,000, 4,000, 2,000 BP
  - Download: https://doi.org/10.7910/DVN/CQWUBI (Harvard Dataverse)
  - Grid-cell format with land-use intensity scores
- **Alternative:** Larson et al. 2014 domestication origin points (compile from paper: doi:10.1073/pnas.1312787110)
  - ~24 independent domestication centers for plants/animals with lat/lon and approximate dates

### Method
1. **ArchaeoGLOBE overlay:**
   a. For the 10,000 BP and 8,000 BP time slices, extract all grid cells within 50km of the Great Circle
   b. Compute mean land-use intensity score for on-corridor vs. global average
   c. Z-score the on-corridor mean against distribution-matched Monte Carlo (1,000 random great circles, same latitude profile)
   d. Repeat for 6,000 BP and 4,000 BP to track temporal evolution

2. **Domestication center test:**
   a. Compile the ~24 independent plant/animal domestication centers with coordinates
   b. Compute each center's distance to the Great Circle
   c. Monte Carlo: generate 10,000 random great circles, compute mean distance to domestication centers
   d. Report: what percentile is the Great Circle's mean distance? Is it closer to domestication centers than expected?

### Output
- `archaeoglobe_corridor_overlay.json` — land-use intensity on/off corridor by epoch
- `domestication_center_distances.csv` — each center, distance to circle, and significance
- `domestication_monte_carlo.json` — percentile ranking, p-value

### Success Criteria
- If the Great Circle passes significantly closer to independent domestication centers than random → the corridor traces the Neolithic agricultural network
- If ArchaeoGLOBE shows early high land-use intensity on corridor → early agricultural adoption on the corridor

---

## Analysis 3: Ancient DNA Spatial Correlation (Exploratory)

### Data
- **Allen Ancient DNA Resource (AADR)** v54+ — ~15,000 published ancient genomes
  - Download: https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes
  - Key fields: lat, lon, date_BP, genetic_cluster/ancestry_component, publication
- If AADR is too complex to parse, use the **supplementary table** from Lazaridis et al. 2022 (doi:10.1038/s41586-022-05085-2) which covers the Southern Arc (exactly the Egypt-Iran corridor) with ~700 genomes

### Method
1. Filter AADR to dates 12,000–4,000 BP (the full Holocene corridor window)
2. Compute distance to Great Circle for each genome
3. For genomes within 200km of the circle, extract ancestry components
4. Test: do genetically similar populations (same primary ancestry component) cluster along the circle more than expected?
   a. For each pair of genetically similar genomes, compute: (i) geographic distance, (ii) angular separation along the circle
   b. If the circle is a migration corridor, genetically similar populations should show lower angular separation along the circle than perpendicular to it
5. Specifically test the Egypt-Iran corridor: do Egyptian and Iranian early Holocene genomes share more ancestry than expected given their geographic distance?

### Output
- `adna_corridor_analysis.json` — ancestry component distributions on/off corridor
- `genetic_similarity_along_vs_across.json` — directional clustering test results

### Lookahead / Bias Warnings
- aDNA sampling is heavily biased toward Europe and the Near East — this overlaps the corridor's strongest signal region. Must normalize for sampling density.
- Ancestry components are model-dependent (ADMIXTURE K choice matters). Report results at multiple K values if possible.
- This analysis is EXPLORATORY — it may not have enough statistical power with current aDNA sample sizes outside Europe.

---

## Analysis 4: Least-Cost-Path Comparison

### Data
- **SRTM elevation data** (already used in geological correlation test)
- **Paleoclimate reconstruction:** Use CHELSA-TraCE21k (doi:10.1038/s41597-023-02549-6) for 10,000 BP temperature and precipitation, OR use the simpler approach of HYDE 3.3 population density as a passability proxy
- **Coastline reconstruction:** Use ICE-6G_C sea level model (Peltier et al. 2015) for LGM and early Holocene coastlines

### Method
1. Build a cost surface for early Holocene (~10,000 BP):
   a. Elevation → slope cost (steeper = more expensive)
   b. Water bodies → impassable (except known shallow crossings)
   c. Desert/extreme aridity → high cost (use precipitation data)
   d. Sea level adjusted to 10,000 BP (~-40m)
2. Compute least-cost-path connecting five nodes: Egypt (Giza), Persepolis, Mohenjo-daro, and optionally Easter Island and Nazca (require ocean crossings — handle separately)
3. For the overland portion (Egypt → Persepolis → Mohenjo-daro):
   a. Compute the least-cost path
   b. Measure its mean distance from the Great Circle
   c. Monte Carlo: compute 10,000 random great circles, measure each's mean distance from the least-cost path
   d. Report: what percentile is the Great Circle? How close is it to the "optimal" overland route?
4. Visualize: overlay the Great Circle, the least-cost path, and the actual locations of early Holocene radiocarbon dates on a single map

### Output
- `least_cost_path_vs_circle.json` — mean separation, percentile ranking
- `cost_surface_map.png` — the cost surface with circle and LCP overlaid
- `holocene_corridor_map.png` — combined visualization

### Success Criteria
- If the Great Circle approximates the least-cost path → parsimonious explanation: the corridor is geographically optimal for overland travel
- If the Great Circle diverges significantly from the LCP → the corridor is NOT just the easiest route; something else is at play

---

## Deliverables Summary
1. `outputs/early_holocene_corridor/RESULTS.md` — narrative summary with go/no-go for each sub-analysis
2. All CSVs and JSONs listed above
3. Figures: temporal profile, domestication center map, cost surface overlay
4. Recommendation on which findings merit inclusion in Paper v2 or a standalone Substack piece
