# POLYNESIAN NAVIGATION & SACRED SITE ALIGNMENT — Research Directive v1.0

**Date:** 2026-03-23
**Author:** Claude (for Ell)
**Objective:** Test whether traditional Polynesian navigation routes and sacred sites (marae) show spatial statistical alignment with the Great Circle or with great circle arcs in general, leveraging the fact that the Great Circle crosses the Pacific through Easter Island and near other Polynesian islands.
**Repository:** great-circle-analysis/
**Output Directory:** outputs/polynesian_navigation/
**Estimated Runtime:** 2–3 hours
**Parallel Safe:** Yes

---

## Background & Motivation

The Great Circle passes directly through Easter Island (Rapa Nui) — one of the most isolated inhabited places on Earth and home to ~900 moai statues (60+ ahu platforms are within the analysis band). It also passes through the Pacific with segments near the Marquesas, Tuamotus, and other island groups.

Polynesian navigators are documented to have used star paths (kaveinga in Tongan, Te Ala o Te Manu in Tuvaluan) — great circle routes on the celestial sphere that were projected onto ocean paths. The star compass systems of Polynesian wayfinding effectively divide the horizon into directional sectors associated with specific stars. Some scholars (Gladwin 1970, Lewis 1972, Finney 1994) have mapped these traditional routes.

**Key questions:**
1. Do documented Polynesian inter-island routes follow great circle arcs more closely than expected?
2. Does the Great Circle specifically correlate with any documented star path or migration route?
3. Do marae (Polynesian temple platforms) cluster along the Great Circle where it passes through the Pacific?

This connects the Great Circle research to one of the most remarkable feats of navigation in human history.

---

## Data Sources

### Polynesian Sacred Sites (Marae)
1. **Easter Island ahu:** Already in dataset (60+ ahu within analysis band)
2. **Society Islands marae:** Compiled in Emory & Sinoto (1965), Green et al. (1967)
   - ~300+ known marae on Tahiti, Raiatea, Huahine, Moorea
   - If no digital database: use the major ones listed in UNESCO/ICOMOS documentation for Taputapuatea (Raiatea) — a UNESCO World Heritage Site
3. **Hawaiian heiau:** State Historic Preservation Division database
   - ~600+ known heiau locations on Hawaiian islands
   - Some publicly available through Hawaii SHPD GIS
4. **Marquesas me'ae:** Archaeological surveys compiled in Suggs (1961), Linton (1925)
5. **Tongan langi (burial mounds):** Documented in McKern (1929) and more recent surveys
6. **OSM "archaeological_site" + "place_of_worship" tags in Polynesia** — already partly captured in the OSM dataset

### Navigation Routes
1. **Traditional star paths:** Compiled in Gladwin (1970) "East is a Big Bird" and Lewis (1972) "We, the Navigators"
   - Routes between major island groups with approximate waypoints
2. **Simulated voyaging routes:** Fitzpatrick & Callaghan (2009, doi:10.1016/j.jas.2008.09.026) — computer-simulated optimal sailing routes between Pacific islands
   - These account for winds, currents, and island visibility
3. **Polynesian Voyaging Society routes:** Documented routes of Hōkūle'a and other recreated voyaging canoes
   - Available from PVS website and published accounts

### Migration/Colonization Routes
- Kirch (2000) "On the Road of the Winds" — canonical Polynesian expansion model
- Estimated colonization dates: Fiji/Tonga/Samoa ~1000 BCE → Society Islands ~800 CE → Easter Island ~1200 CE → Hawaii ~1000 CE → New Zealand ~1250 CE
- These routes with dates allow a "wavefront" analysis similar to Directive 04

---

## Analysis 1: Marae Distribution Along the Great Circle

### Method
1. Compile all available marae/ahu/heiau locations across Polynesia
2. Compute distance from each to the Great Circle
3. Given the Pacific's vastness, use a wider band: test at 50km, 100km, 200km, 500km
4. **Monte Carlo:** Generate 10,000 random great circles passing through the Pacific (poles constrained to produce Pacific-crossing circles). Compute marae enrichment for each.
5. **Key sub-question:** The Great Circle passes near specific island groups. Which islands are within 200km of the circle?
   - Easter Island: ON the circle
   - Pitcairn Islands: distance?
   - Marquesas: distance?
   - Tuamotus: distance?
   - Society Islands (Tahiti/Raiatea): distance?
   - Cook Islands: distance?
6. Compute: on islands that the circle passes near, are marae more numerous or more monumental than on islands the circle misses?

### Output
- `marae_enrichment.json` — enrichment at each distance threshold with Monte Carlo significance
- `island_distances.csv` — every major Polynesian island group with distance to circle
- `pacific_circle_map.png` — Great Circle through the Pacific with island groups and marae marked

---

## Analysis 2: Navigation Route Geometry

### Method
1. Digitize 10+ traditional/reconstructed Polynesian voyaging routes
2. For each route:
   a. Fit a great circle arc to the route waypoints
   b. Compute RMS deviation from the best-fit great circle
   c. Compare: is the route closer to a great circle than a rhumb line (constant bearing)?
   d. Monte Carlo: generate random paths between the same endpoints, compute RMS from great circle
3. **Specifically test:** Do any documented routes parallel the Great Circle?
   a. Compute the angular separation between each route's best-fit great circle and the Great Circle
   b. If any route's great circle is within 5° of the Great Circle's orientation → significant
4. **Star path overlay:** The traditional Polynesian star compass assigns specific stars to specific bearings. At Easter Island, what star corresponds to the Great Circle's bearing? Is it a star of known navigational importance?

### Output
- `route_geometries.csv` — digitized routes with great circle fits
- `route_circle_comparison.json` — angular separation from Great Circle for each route
- `star_compass_bearing.json` — Great Circle bearing at key islands mapped to star compass directions

---

## Analysis 3: Easter Island Ahu Orientation vs. Circle Bearing

### Method
This is a focused sub-test similar to the Nazca directive (06), but for Easter Island:
1. Easter Island ahu have well-documented orientations (most face inland, some face the sea)
2. Compile ahu orientations from published surveys (Martinsson-Wallin 1994, Van Tilburg 1994)
3. Compute the Great Circle's bearing at Easter Island
4. Test: do ahu orientations correlate with the circle bearing?
   - Rayleigh test on angular offsets
   - V-test directed at the circle bearing
5. **Alternative test:** Ahu are distributed around the island's perimeter. Does ahu density peak where the Great Circle intersects the coastline?

### Output
- `ahu_orientation_test.json` — Rayleigh/V-test results
- `ahu_density_by_bearing.json` — coastal density relative to circle crossing points
- `easter_island_circle_overlay.png` — island map with ahu, their orientations, and circle path

---

## Analysis 4: Colonization Wavefront and the Circle

### Method
1. Using the established Polynesian colonization chronology (Kirch 2000, Wilmshurst et al. 2011):
   - Plot each island group's colonization date against its distance from the Great Circle
   - Test: were islands closer to the circle colonized earlier than islands farther away?
   - Spearman correlation between distance-to-circle and colonization date
2. This tests whether the Great Circle traces a preferential colonization route through the Pacific

### Output
- `colonization_vs_distance.json` — correlation test results
- `colonization_timeline.png` — timeline plot colored by circle distance

---

## Caveats
- Polynesian sacred site databases are fragmentary — many marae have been destroyed or are unrecorded
- The Pacific is so vast that even "near" the circle means hundreds of km for some island groups
- Polynesian navigation routes were constrained by winds and currents, not just geometry — a great circle route is not always the fastest sailing route
- Easter Island's extreme isolation makes it a statistical outlier regardless of any alignment
- The colonization chronology is debated (especially early dates for some islands)

---

## Deliverables
1. `outputs/polynesian_navigation/RESULTS.md`
2. All data files, figures, and maps
3. The Pacific map with the Great Circle and island sacred sites would be visually stunning
4. Any connection between Polynesian star navigation and the Great Circle geometry would be extraordinary
