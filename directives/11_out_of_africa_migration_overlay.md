# OUT-OF-AFRICA MIGRATION CORRIDOR OVERLAY — Research Directive v1.0

**Date:** 2026-03-24
**Author:** Claude (for Ell)
**Objective:** Test whether the Great Circle traces the Out-of-Africa southern coastal migration route (~70,000–45,000 BP), connecting the deepest-time archaeological sites on the corridor (PNG at 55,000 BP) to the Holocene campsites and Bronze Age monuments through a single persistent geographic corridor spanning 60,000 years.
**Repository:** great-circle-analysis/
**Output Directory:** outputs/out_of_africa_overlay/
**Estimated Runtime:** 2–3 hours
**Parallel Safe:** Yes

---

## Background & Motivation

Previous directives have established:
- **Directive 01:** The corridor is enriched for early Holocene campsites (10,500–6,500 BCE), passes closer to independent agricultural origins than 99.8% of random circles (p = 0.002), and approximates the optimal overland route through the early Holocene landscape (p = 0.008) via a **southern arc** (Negev → Arabia → Persia) rather than the northern Fertile Crescent.
- **Directive 04:** Bronze Age monuments sit 50× closer to early Holocene sites than chance (p = 0.004) — the corridor has been continuously occupied for 8,000+ years.

But the corridor contains sites far older than the Holocene:
- **Mololo Cave, PNG:** 54km from circle, **55,000 BP** — among the earliest human presence in Sahul
- **Kebar Valley, PNG:** 16km from circle, **26,000 BP** — deep ice age occupation
- **Caution Bay, PNG:** 40km from circle — significant early occupation
- **That Nang Ing, SE Asia:** 0.4km from circle

These sites predate the Holocene by tens of thousands of years. The agricultural corridor explanation cannot account for them. But the **southern dispersal route** — the leading model for how anatomically modern humans left Africa and reached Australia/PNG between ~70,000–50,000 BP — traces almost exactly the same geographic arc.

**The hypothesis:** The Great Circle doesn't just trace a Holocene agricultural corridor. It traces the original Out-of-Africa southern coastal migration route — a geographic corridor that has been used by humans for 60,000+ years. The Holocene campsites, the agricultural origins, and the Bronze Age monuments are successive layers of activity on a corridor that was first walked by the earliest humans to leave Africa.

If confirmed, this would unify ALL the significant findings into a single temporal narrative:
1. ~60,000 BP: First humans traverse the corridor (southern coastal route) → PNG, SE Asia sites
2. ~10,000 BP: Post-glacial populations reoccupy the corridor as a campsite/transit route → p3k14c enrichment
3. ~8,000 BP: Agriculture independently develops at multiple points along the corridor → domestication center proximity (p = 0.002)
4. ~5,000 BP: Monuments are built at the same locations → monument-settlement divergence, continuity (p = 0.004)

---

## Data Sources

### Published Southern Dispersal Route Reconstructions

1. **Reyes-Centeno et al. (2014)** — doi:10.1073/pnas.1323666111
   - "Genomic and cranial phenotype data support multiple modern human dispersals from Africa and a southern route into Asia"
   - Provides a reconstructed southern route with geographic waypoints based on craniometric and genomic data
   - Key waypoints extractable from figures and supplementary materials

2. **Tassi et al. (2017)** — doi:10.1038/s41598-017-01802-y
   - "Early modern human dispersal from Africa: genomic evidence for multiple waves of migration"
   - Models both southern and northern routes with dates and waypoints

3. **Petraglia et al. (2010)** — doi:10.1073/pnas.0911799107
   - "Middle Paleolithic assemblages from the Indian subcontinent before and after the Toba super-eruption"
   - Maps archaeological sites along the southern route through India

4. **Groucutt et al. (2015)** — doi:10.1016/j.quascirev.2015.01.006
   - "Stone tool assemblages and models for the dispersal of Homo sapiens out of Africa"
   - Comprehensive review with mapped sites along both dispersal routes

5. **Field & Lahr (2006)** — doi:10.1016/j.jas.2005.07.002
   - "Assessment of the southern dispersal: GIS-based analyses of potential routes at Oxygen Isotope Stage 4"
   - GIS-based least-cost path analysis of the southern route during OIS-4 (~74,000–59,000 BP)
   - This is directly comparable to our Directive 01 least-cost path analysis but for a much earlier epoch

6. **Bulbeck (2007)** — "Where River Meets Sea: A parsimoniously calibrated SE Asian model for early human migration"
   - Maps the coastal route through SE Asia to Sahul

### Archaeological Sites Along the Southern Route

7. **Armitage et al. (2011)** — doi:10.1126/science.1199113
   - Jebel Faya, UAE — anatomically modern human tools at 125,000 BP, on the Arabian Peninsula
   - One of the key waypoints for the southern route

8. **Rose et al. (2011)** — doi:10.1086/658436
   - "The Nubian Complex of Dhofar, Oman" — African-style stone tools in southern Arabia

9. **Petraglia et al. (2007)** — doi:10.1126/science.1141564
   - Jwalapuram, India — Middle Paleolithic tools below and above the Toba ash layer (~74,000 BP)
   - Direct evidence of humans on the southern route through India

10. **O'Connell et al. (2018)** — doi:10.1073/pnas.1807724115
    - "When did Homo sapiens first reach Southeast Asia and Sahul?"
    - Comprehensive review of earliest dates for human arrival across the southern route

### Genetic/aDNA Evidence for the Southern Route

11. **Rasmussen et al. (2011)** — doi:10.1126/science.1211177
    - Aboriginal Australian genome shows deep divergence from other non-African populations, consistent with early southern route dispersal

12. **Malaspinas et al. (2016)** — doi:10.1038/nature18299
    - "A genomic history of Aboriginal Australia"
    - Confirms southern route and early divergence (~72,000 BP)

13. **Pagani et al. (2016)** — doi:10.1038/nature19792
    - "Genomic analyses inform on migration events during the peopling of Eurasia"
    - Models multiple dispersal waves with route reconstructions

---

## Analysis 1: Route Fit Comparison

### Method
1. **Digitize the southern dispersal route** from the published sources above:
   - Extract waypoints from Reyes-Centeno (2014), Tassi (2017), and Field & Lahr (2006)
   - If figures are the only source, digitize from maps using georeferenced coordinates
   - Create a composite "consensus southern route" by averaging waypoints across publications
   - The route should run approximately: East Africa → Bab el-Mandeb/southern Arabia → coastal Iran/Makran coast → coastal India → SE Asia → Sahul/PNG

2. **Compute the fit between the Great Circle and the southern route:**
   a. For each waypoint on the consensus route, compute distance to the Great Circle
   b. Compute mean distance, median distance, max distance
   c. Compute the RMS deviation of the route from the Great Circle

3. **Monte Carlo significance:**
   a. Generate 10,000 random great circles
   b. For each, compute mean distance to the consensus southern route waypoints
   c. Percentile rank the Great Circle
   d. Report: how many random circles fit the southern route as well as or better than the Great Circle?

4. **Also compute fit to the NORTHERN route** (through the Levant → Anatolia → Central Asia → East Asia):
   - If the circle fits the southern route significantly better than the northern route, that confirms the directional specificity already found in Directive 01

### Output
- `southern_route_waypoints.csv` — digitized consensus route
- `route_fit_comparison.json` — mean distances, RMS, Monte Carlo percentile for both routes
- `route_overlay_map.png` — world map showing Great Circle + southern route + northern route

---

## Analysis 2: Deep-Time Archaeological Site Enrichment

### Method
1. **Compile all known archaeological sites >40,000 BP along the Out-of-Africa southern route:**
   - From the papers above and supplementary databases
   - Key sites to include:
     - **Africa (departure):** Blombos Cave, Border Cave, Klasies River Mouth
     - **Arabia:** Jebel Faya (~125,000 BP), Dhofar Nubian Complex sites, Al Wusta (~85,000 BP)
     - **India:** Jwalapuram (~74,000 BP), Attirampakkam, Bhimbetka
     - **SE Asia:** Tam Pa Ling (Laos, ~63,000 BP), Niah Cave (Borneo, ~40,000 BP), That Nang Ing
     - **Sahul/PNG:** Mololo Cave (55,000 BP), Kebar Valley (26,000 BP), Madjedbebe (Australia, ~65,000 BP), Lake Mungo (~42,000 BP)

2. **For each site:** compute distance to the Great Circle

3. **Enrichment test:**
   a. What proportion of sites >40,000 BP are within 200km of the circle? (Wider band appropriate for this epoch — coastlines were different, data is sparse)
   b. Monte Carlo: 10,000 random great circles, same enrichment metric
   c. Report percentile

4. **Temporal layering test (THE KEY ANALYSIS):**
   a. Bin all on-corridor sites by epoch:
      - >50,000 BP (initial dispersal)
      - 50,000–25,000 BP (early occupation)
      - 25,000–12,000 BP (Last Glacial Maximum)
      - 12,000–8,000 BP (early Holocene — connects to Directive 01)
      - 8,000–5,000 BP (Neolithic — connects to domestication centers)
      - 5,000–2,000 BP (Bronze Age — connects to Directive 04)
   b. Compute enrichment at each epoch
   c. **Test for continuity:** Is there evidence of continuous occupation at all epochs, or are there gaps?
   d. This directly tests whether the corridor shows a continuous 60,000-year occupation sequence

### Output
- `deep_time_sites.csv` — all compiled sites with distances and dates
- `temporal_layering.json` — enrichment by epoch
- `continuity_test.json` — gap analysis across epochs
- `deep_time_timeline.png` — timeline visualization showing corridor occupation from 60,000 BP to present

---

## Analysis 3: Paleoclimate Corridor Persistence

### Concept
The southern route was viable during OIS-4 (~74,000–59,000 BP) when sea levels were lower and Arabia was periodically greener. But it also needed to remain traversable during the Last Glacial Maximum (~26,000–20,000 BP) and the Younger Dryas (~12,800–11,600 BP). **Was the corridor continuously habitable for 60,000 years?**

### Data
- **Sea level curve:** Lambeck et al. (2014) — global sea level over the past 135,000 years
  - doi:10.1073/pnas.1411762111
- **Paleoclimate:**
  - CHELSA-TraCE21k for the last 21,000 years (already used in Directive 01)
  - For deeper time: use the Singarayer & Valdes (2010) HadCM3 paleoclimate simulations
  - Or simpler: use the HYDE 3.3 population estimates (already in dataset) as a habitability proxy for the Holocene portion
- **Green Arabia periods:** Reviewed in Groucutt & Petraglia (2012, doi:10.1016/j.quascirev.2012.09.011)
  - Multiple humid phases: MIS 5 (~130,000–71,000 BP), parts of MIS 3 (~60,000–27,000 BP)
  - These correspond to when the southern route through Arabia was passable

### Method
1. Sample the Great Circle at 1° intervals through the overland portion (Africa → SE Asia)
2. For each sample point, extract:
   a. Modern elevation and distance to coast
   b. LGM (~20,000 BP) coastline distance (sea level -120m)
   c. OIS-4 (~65,000 BP) coastline distance (sea level -80m)
   d. Estimated precipitation at each paleoclimate epoch (from published reconstructions)
3. **Corridor viability assessment:**
   a. At each epoch, is the corridor traversable? (Not blocked by ice, desert, or deep water)
   b. Identify bottlenecks: sections where the corridor narrows to a single passable route
   c. Map the corridor's width at each epoch
4. **The Persian Gulf factor:**
   a. During LGM and OIS-4, the Gulf was dry land — the circle passes through the dry basin
   b. During "Green Arabia" phases, the Arabian segment was passable
   c. When was the corridor at its widest (most habitable) vs. narrowest (barely passable)?
5. **Visualization:** A time-lapse style figure showing the corridor's habitability from 70,000 BP to present, with archaeological sites appearing at their respective dates

### Output
- `corridor_viability_by_epoch.json` — habitability scores along the corridor at each epoch
- `bottleneck_analysis.json` — narrowest passable points at each epoch
- `paleoclimate_corridor_timelapse.png` — multi-panel figure showing corridor evolution
- `persian_gulf_dry_periods.json` — when the Gulf crossing was passable

---

## Analysis 4: Genetic Distance Along vs. Across the Corridor

### Concept
If the corridor has been a migration highway for 60,000 years, populations along the corridor should be more genetically similar to each other (relative to geographic distance) than populations perpendicular to the corridor. This is testable with published genome-wide data.

### Data
- **HUGO Pan-Asian SNP Consortium (2009)** — doi:10.1126/science.1177074
  - 73 Asian/Oceanian populations with genome-wide SNP data
  - Genetic distances (Fst) between all population pairs published in supplementary materials
- **1000 Genomes Project** populations along and perpendicular to the corridor
- **Malaspinas et al. (2016)** — Aboriginal Australian genomes with published Fst to other populations

### Method
1. For each population pair in the dataset:
   a. Compute geographic distance
   b. Compute the pair's "corridor alignment" — are both populations close to the Great Circle? Or is one on-corridor and one off-corridor?
   c. Classify pairs as: both on-corridor (within 500km), both off-corridor, or mixed
2. **Isolation-by-distance residual test:**
   a. Fit a standard isolation-by-distance regression (Fst ~ geographic distance)
   b. Compute residuals for each pair
   c. Test: do on-corridor pairs have more negative residuals (genetically MORE similar than distance predicts) than off-corridor pairs?
   d. This is the key test: the corridor should reduce effective genetic distance between populations along it
3. **Monte Carlo:** Replace the Great Circle with 10,000 random great circles, recompute the residual test for each. Percentile rank.

### Output
- `genetic_distance_pairs.csv` — population pairs with Fst, geographic distance, corridor classification
- `isolation_by_distance_residuals.json` — residual test results
- `genetic_corridor_effect.json` — Monte Carlo percentile
- `ibd_plot.png` — isolation-by-distance plot with on-corridor pairs highlighted

### Caveats
- Genetic similarity is shaped by many forces (drift, selection, admixture), not just migration corridors
- Post-Neolithic migrations (Indo-European expansion, Bantu expansion, Austronesian expansion) may overwrite deep-time signals
- The test may be underpowered for the Pacific segment where populations are small and isolated
- Frame results carefully: "consistent with" a migration corridor, not "proves" one

---

## Analysis 5: The Full 60,000-Year Timeline (Synthesis)

### Method
This is a synthesis visualization combining all results from this directive AND Directives 01 and 04:

1. Create a master timeline of corridor activity from 60,000 BP to 2,000 BP
2. Plot on a dual-axis figure:
   - X-axis: time (60,000 BP → present)
   - Y-axis (top): number of dated sites on the corridor per 1,000-year bin
   - Y-axis (bottom): corridor viability/habitability score
3. Annotate key events:
   - First human arrival in Arabia (~125,000 BP — Jebel Faya)
   - Southern coastal migration (~70,000–50,000 BP)
   - First arrival in Sahul/PNG (~65,000–50,000 BP)
   - Toba eruption (~74,000 BP)
   - Last Glacial Maximum (~26,000–20,000 BP)
   - Green Arabia phases
   - Persian Gulf flooding (~8,000 BP)
   - Early Holocene campsites (Directive 01)
   - Agricultural origins (Directive 01)
   - Bronze Age monument peak (Directive 04)
4. Overlay the Great Circle's enrichment Z-score at each epoch (where data permits)

### Output
- `master_timeline_60k.png` — THE capstone visualization for the entire project
- `timeline_data.json` — all data points underlying the figure
- `narrative_summary.md` — written synthesis tying all directives together

---

## Lookahead / Bias Warnings
- The southern dispersal route is the **leading hypothesis** but not the only model — some geneticists argue for a single northern+southern wave, others for multiple waves. Present the southern route as one model being tested, not established fact.
- Archaeological sites >40,000 BP are extremely sparse — there may be fewer than 50 sites total along the entire route. Statistical power will be limited. Report effect sizes and confidence intervals, not just p-values.
- Coastlines were radically different during the dispersal period — sites that were coastal are now underwater (especially in the Persian Gulf and Sunda Shelf). The corridor may have more sites than we can ever find.
- The temporal layering test may show gaps (especially during the LGM when human populations contracted). Gaps are informative, not failures — they show when the corridor was less used.
- The genetic distance analysis may be confounded by post-Neolithic population movements. This is the weakest analysis and should be framed as exploratory.

---

## Go / No-Go Interpretation

| Finding | Implication |
|---------|------------|
| Circle fits southern route (p < 0.05) | The corridor IS the Out-of-Africa southern route |
| Deep-time sites enriched on corridor | The corridor was used from the very beginning of human dispersal |
| Continuous temporal layering (no major gaps) | 60,000 years of persistent corridor use |
| Genetic similarity enhanced along corridor | The corridor functioned as a gene-flow highway |
| Paleoclimate shows persistent habitability | The corridor remained viable across glacial cycles |

If 3+ of these return positive → the Great Circle traces the oldest continuously used human migration corridor on Earth. That's a unifying explanation for everything: PNG, the Holocene campsites, the agricultural origins, and the monuments.

---

## Deliverables
1. `outputs/out_of_africa_overlay/RESULTS.md` — full narrative
2. All CSVs, JSONs, and figures listed above
3. The master 60,000-year timeline figure is the single most important visualization in the entire project
4. If the route fit is significant → this reframes the entire paper from "mysterious monument alignment" to "the world's oldest human corridor, still marked by the monuments of every civilization that walked it"
