# NAZCA GEOGLYPH ORIENTATION ANALYSIS — Research Directive v1.0

**Date:** 2026-03-23
**Author:** Claude (for Ell)
**Objective:** Test whether Nazca geoglyph orientations (the directions figures "face" or the bearing of linear geoglyphs) correlate with the Great Circle's bearing at Nazca, using the Sakai et al. 2024 PNAS dataset of 733 figurative geoglyphs.
**Repository:** great-circle-analysis/
**Output Directory:** outputs/nazca_orientation/
**Estimated Runtime:** 1–2 hours
**Parallel Safe:** Yes

---

## Background & Motivation

Sakai et al. (2024, PNAS, doi:10.1073/pnas.2407652121) identified 733 figurative geoglyphs on the Nazca Pampa and surrounding foothills using AI-assisted satellite imagery analysis. Many of these have clear orientational axes (the direction a figure faces, the long axis of a linear geoglyph, or the "entry" direction of a path-type geoglyph). The Great Circle passes directly through the Nazca region — one of its strongest cluster points.

The existing astronomical alignment test (paper v3) checked whether monument orientations correlate with the circle bearing across ALL sites globally. It came back null. But that global test diluted any site-specific signal. This directive focuses exclusively on Nazca, where the density of oriented features is unmatched and the Great Circle proximity is tight.

Additionally, the linear Nazca lines (not the figurative geoglyphs) are famously oriented in specific directions. Some researchers (Aveni, Ruggles) have proposed astronomical alignments; others (Johnson) have proposed hydrological alignments. The Great Circle's bearing at Nazca provides a third hypothesis to test.

---

## Data

### Primary: Sakai et al. 2024 Dataset
- **Paper:** doi:10.1073/pnas.2407652121
- **Supplementary data:** Check if geoglyph coordinates and orientations are in the supplementary materials
- If not publicly available: extract from the paper's figures and tables
- **Key fields needed:** lat, lon, type (figurative/linear), orientation_angle (degrees from north), size, period (if available)

### Secondary: Aveni's Line Database
- Aveni (1990, 2000) catalogued ~800 linear features (lines, trapezoids, ray centers) with bearings
- Published in: "The Lines of Nazca" (1990) and subsequent papers
- If digital data unavailable: the key result (bimodal distribution peaking at ~10° and ~80°) is well-documented and can be used as a comparison

### Tertiary: Clarkson & Dorn Line Surveys
- Additional orientation data from more recent surveys
- Check: doi:10.1016/j.jasrep.2016.03.025

---

## Analysis 1: Great Circle Bearing at Nazca

### Method
1. Compute the Great Circle's bearing at Nazca (approximately 14.7°S, 75.1°W):
   - The bearing is the azimuth of the Great Circle as it passes through this point
   - This gives two values: the bearing in each direction (B and B+180°)
   - Use the standard formula: bearing from point P along great circle with pole at (φ_p, λ_p)
2. Also compute the bearing at 5 specific points across the Nazca Pampa (the plateau is ~50km across; the bearing may vary slightly)

### Expected Result
- The bearing should be roughly E-W through Nazca (the circle runs roughly east-west through Peru at this latitude)
- Compute the exact value to use as the test angle

---

## Analysis 2: Geoglyph Orientation vs. Circle Bearing

### Method (if Sakai orientation data available)
1. For each geoglyph with a measurable orientation:
   a. Compute the angular offset between the geoglyph's orientation and the circle's bearing at that point
   b. Use the minimum of |offset| and |180° - offset| (since most geoglyphs are symmetric — a figure facing NE is also "aligned" at NE+180°=SW)
2. **Rayleigh test on offsets:**
   - If geoglyphs are randomly oriented relative to the circle → offsets are uniformly distributed (0–90°)
   - If aligned → offsets cluster near 0°
   - If perpendicular → offsets cluster near 90°
   - Compute Rayleigh statistic on the doubled angles (2θ to handle axial data)
3. **V-test (directed):**
   - Test specifically against the hypothesis that offsets cluster at 0° (alignment with circle)
   - More powerful than Rayleigh if the expected direction is known
4. **Monte Carlo control:**
   - Replace the Great Circle bearing with 10,000 random bearings
   - For each, compute the Rayleigh statistic on offsets
   - Percentile rank the actual Great Circle's Rayleigh statistic

### Method (if only Aveni's aggregate data available)
1. Use the known bimodal distribution of Nazca line orientations (peaks at ~10° and ~80° from north)
2. Compute the Great Circle bearing at Nazca
3. Test: is the circle bearing within 5° of either modal peak?
4. Report qualitatively

---

## Analysis 3: Spatial Gradient of Alignment

### Method
1. If geoglyph orientations are available with precise lat/lon:
   a. Divide the Nazca Pampa into 5km × 5km grid cells
   b. In each cell, compute the mean angular offset between geoglyphs and the circle bearing
   c. Test: does alignment improve closer to the circle? (i.e., is there a distance-dependent effect?)
   d. Spearman correlation between distance-to-circle and mean angular offset
2. This tests whether geoglyphs *near* the circle are more aligned than those farther away — a much stronger signal than just testing all geoglyphs together

---

## Analysis 4: Line Centers ("Ray Centers") and the Circle

### Background
The Nazca lines radiate from ~60 "ray centers" — hilltops or natural features where multiple lines converge. Some researchers believe these are the organizing nodes of the Nazca line system.

### Method
1. Compile ray center locations from Aveni (1990) or subsequent surveys
2. Compute each ray center's distance to the Great Circle
3. Test: do ray centers cluster closer to the circle than the general Nazca Pampa surface?
4. Specifically: does any ray center sit ON the circle (within 1km)?

---

## Output
- `circle_bearing_at_nazca.json` — precise bearing values
- `orientation_analysis.json` — Rayleigh test, V-test, Monte Carlo results
- `spatial_gradient.json` — distance-dependent alignment test
- `ray_center_distances.csv` — ray center locations and distances to circle
- `nazca_orientation_rose_diagram.png` — rose diagram of geoglyph orientations with circle bearing marked
- `nazca_alignment_map.png` — map showing geoglyphs colored by alignment degree with circle overlaid

---

## Caveats
- "Orientation" of figurative geoglyphs is subjective — a hummingbird's "facing direction" is interpretive
- Linear geoglyphs have unambiguous bearings; figurative ones are fuzzier
- The circle bearing at Nazca may happen to align with the dominant wind direction or the valley axis — test these alternative explanations
- Nazca geoglyph creation spans ~800 BCE to ~600 CE (1,400 years). If alignment is found, test whether it's epoch-dependent.

---

## Deliverables
1. `outputs/nazca_orientation/RESULTS.md`
2. All data files and figures
3. Rose diagram is highly visual — good Substack/social media material regardless of result
