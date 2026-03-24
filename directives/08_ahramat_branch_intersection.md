# AHRAMAT BRANCH PERPENDICULAR INTERSECTION — Research Directive v1.0

**Date:** 2026-03-23
**Author:** Claude (for Ell)
**Objective:** Quantitatively test the geometric relationship between the Great Circle and the recently discovered Ahramat Branch of the Nile, focusing on: (1) the intersection angle, (2) whether the intersection predicts undiscovered sites, and (3) the probability of this geometric coincidence.
**Repository:** great-circle-analysis/
**Output Directory:** outputs/ahramat_intersection/
**Estimated Runtime:** 2–3 hours
**Parallel Safe:** Yes

---

## Background & Motivation

Ghoneim et al. (2024, Communications Earth & Environment, doi:10.1038/s43247-024-01379-7) discovered that 31 Old Kingdom pyramids and major mortuary temples align along a now-extinct branch of the Nile they named the "Ahramat Branch." This branch ran roughly N-S from the Faiyum region to the Memphis area, providing waterfront access to pyramid construction sites.

The Great Circle passes through the Memphis pyramid field nearly **perpendicularly** to this reconstructed branch. This means the densest monument cluster on the entire Great Circle sits at the intersection of two alignment systems:
- **N-S axis:** The Ahramat Branch (fluvial, practical — provided water transport for construction)
- **E-W axis:** The Great Circle (mechanism unknown)

This "crossroads" geometry has not been quantitatively tested. The key questions:
1. How perpendicular is the intersection actually? (Is it 85°? 90°? 70°?)
2. What are the odds of a random great circle intersecting this specific river branch at this angle and at this location?
3. Do the newest discoveries (El-khteeb 2025 Saqqara GPR) fall closer to the intersection point?

---

## Data

### Ahramat Branch Reconstruction
- **Primary:** Ghoneim et al. 2024 supplementary data
  - The paper includes a reconstructed centerline of the Ahramat Branch derived from radar remote sensing (SIR-C, SRTM, Sentinel-1) and soil sample analysis
  - Extract waypoints from the supplementary figures/data
  - If digital data unavailable: digitize from the paper's Figure 1 (the branch route is clearly mapped)
- **Branch extent:** Approximately from 29.3°N to 30.0°N, running roughly N-S with meanders
- **Width estimates:** 200–500m at various points (paper provides cross-sections)

### Pyramid Locations
- **Use existing Megalithic Portal data** for pyramids within 50km of the Great Circle
- **Supplement with:** Lehner (1997) "The Complete Pyramids" — authoritative positions for all ~138 Egyptian pyramids
- **Key subset:** The 31 pyramids identified by Ghoneim et al. as aligned with the Ahramat Branch

### El-khteeb 2025 Saqqara Discoveries
- **Paper:** El-khteeb et al. (2025, already referenced in project docs)
- Key findings: buried walls, road networks, structures 8m deep near the Saqqara pyramid complex
- Extract coordinates of newly discovered features from the paper

### Great Circle Parameters
- Pole: 59.682122°N, 138.646087°W
- Circle bearing at Memphis (29.85°N, 31.25°E): compute precisely

---

## Analysis 1: Intersection Geometry

### Method
1. Digitize the Ahramat Branch centerline as a series of waypoints
2. Compute the Great Circle's bearing at each point along the branch
3. Find the intersection point(s): where does the Great Circle cross the branch?
   - Use line-intersection on the sphere (the branch is approximately linear at the scale of the intersection)
4. At the intersection point:
   a. Compute the branch's local bearing (from the two nearest waypoints)
   b. Compute the Great Circle's bearing
   c. Compute the intersection angle: |circle_bearing - branch_bearing|
   d. Report: how close to 90° is it?
5. Compute the intersection point's proximity to the densest part of the pyramid field
   - Distance to Giza, Saqqara, Dahshur, Abu Sir, and Abusir

### Output
- `intersection_geometry.json` — intersection point, angle, distances to major pyramid sites
- `intersection_map.png` — map showing the Ahramat Branch, the Great Circle, and pyramid positions with the intersection point highlighted

---

## Analysis 2: Monte Carlo Probability

### Method
1. **Test 1: Random Circle × Fixed Branch**
   a. Generate 100,000 random great circles (uniformly distributed poles on the sphere)
   b. For each, compute: does it cross the Ahramat Branch within 10km of a pyramid?
   c. If yes, compute the intersection angle
   d. Count: how many random circles produce a near-perpendicular (>80°) crossing within 10km of a pyramid?
   e. Report: probability of the observed configuration by chance

2. **Test 2: Random Circle × Fixed Branch × Monument Density**
   a. For circles that DO cross the branch, compute: how many pyramids are within 10km of the intersection point?
   b. Compare to the Great Circle's count
   c. This tests whether the Great Circle hits the branch at the *densest* part of the pyramid field

3. **Test 3: Conditional on passing through Egypt**
   a. Of random circles that pass through the Memphis-to-Faiyum latitude band (29°N–30.5°N):
      - What fraction cross the Ahramat Branch?
      - Of those, what fraction are near-perpendicular?
      - Of those, what fraction intersect at the pyramid-dense segment?
   b. This isolates the geometric coincidence from the "the circle passes through Egypt" coincidence

### Output
- `monte_carlo_intersection.json` — probabilities for each test level
- `random_circle_angle_histogram.png` — distribution of intersection angles for circles that cross the branch

---

## Analysis 3: Predictive Test — Do New Discoveries Fall on the Intersection?

### Method
1. The El-khteeb 2025 discoveries at Saqqara include GPR-detected features not visible from the surface
2. Extract coordinates of newly discovered buried structures
3. Compute their distances to:
   a. The Great Circle
   b. The Ahramat Branch
   c. The intersection point
4. Compare to the distances of previously known structures at Saqqara
5. **Prediction test:** If the intersection geometry is meaningful, newly discovered (i.e., deeper/older) structures should cluster CLOSER to the intersection point than visible surface structures
6. Also check: Hawass's 2020 Hermopolis Magna discoveries and any other very recent Egyptian archaeological finds — do any fall near the Great Circle?

### Output
- `new_discoveries_proximity.json` — distances of El-khteeb features to circle, branch, and intersection
- `prediction_test.json` — comparison of new vs. known structure distances

---

## Analysis 4: The "Crossroads" Hypothesis — Are Other Crossroads Sites Also Monument-Dense?

### Concept
If the Memphis intersection of the Great Circle and the Ahramat Branch is significant, we should ask: do OTHER intersections of the Great Circle with major river systems also show monument clustering?

### Method
1. Identify all major river crossings of the Great Circle:
   - Nile (already analyzed — Ahramat Branch)
   - Indus (near Mohenjo-daro)
   - Amazon tributaries (near Tapajós, 0.2km)
   - Mekong/Chao Phraya (SE Asia segment)
   - Tigris/Euphrates (if circle crosses)
2. At each crossing point:
   a. Compute the intersection angle between the circle and the river
   b. Count monuments within 10km of the intersection
   c. Compare to the average monument density along the circle
3. **Test:** Are river-crossing points enriched for monuments relative to the rest of the circle?
4. Monte Carlo: for random great circles, are river crossings enriched for monuments?

### Output
- `river_crossings.csv` — all crossings with angles, monument counts
- `crossing_enrichment.json` — enrichment test results
- `crossroads_map.png` — all river-circle intersections with monument counts

---

## Lookahead / Bias Warnings
- The Ahramat Branch reconstruction has uncertainties — it's derived from remote sensing, not excavation. The centerline may shift with future data.
- The intersection angle depends on which segment of the branch is used (the branch meanders). Report the angle at multiple scales (1km, 5km, 10km averaging).
- Digitization error from paper figures introduces ~100m uncertainty. This is small relative to the scales being tested.
- The "prediction test" (Analysis 3) is weak if El-khteeb's features are only at one location — it becomes a single data point.

---

## Deliverables
1. `outputs/ahramat_intersection/RESULTS.md`
2. All data files, figures, and maps
3. The intersection map alone is publication-quality material
4. If the perpendicular intersection is confirmed and Monte Carlo probability is low → standalone paper or killer Substack piece
5. If the crossroads hypothesis generalizes to other river crossings → major finding
