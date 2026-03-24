# ASTRONOMICAL EPOCH ANALYSIS — Research Directive v1.0

**Date:** 2026-03-23
**Author:** Claude (for Ell)
**Objective:** Test whether the Great Circle has astronomical significance at any historical epoch — as a frozen equator, via stellar horizon events at key sites, via its pole's celestial projection, or via galactic plane intersection.
**Repository:** great-circle-analysis/
**Output Directory:** outputs/astronomical_epoch/
**Estimated Runtime:** 2–3 hours
**Parallel Safe:** Yes

---

## Background & Motivation

The existing astronomical alignment test (paper v3) checked whether individual monument *orientations* correlate with the circle's *bearing* at each site. Result: null (Rayleigh p = 0.75). But that test was narrow — it asked "do monuments point along the circle?" rather than "is the circle itself astronomically meaningful?" There are four completely untested astronomical questions, each addressing a different geometric relationship between the Great Circle and the sky.

---

## Analysis 1: The Circle as a Frozen Equator

### Concept
Earth's rotational axis precesses with a ~25,772-year period (axial precession). Additionally, true polar wander (TPW) shifts the geographic poles over geological timescales. At some past epoch, the orientation of Earth's axis relative to the mantle was different — meaning the equator was in a different position.

**Question:** At what epoch (if any) would the Great Circle have coincided with Earth's equator?

### Method
1. **Axial precession only (short timescale, <100,000 years):**
   - The Great Circle's pole is at 59.682°N, 138.646°W
   - Earth's current rotational pole is at 90°N
   - The angular distance between the circle's pole and Earth's pole is: arccos(sin(59.682°)) ≈ 30.3°
   - For the Great Circle to be the equator, the rotational pole would need to be at the circle's pole position
   - Axial precession changes the direction of the pole relative to the *stars* but NOT relative to the mantle/geography — the geographic pole stays at ~90°N
   - **Therefore, axial precession alone cannot make the Great Circle into the equator.** Document this clearly.

2. **True Polar Wander (geological timescale):**
   - TPW shifts the geographic pole relative to the mantle at ~1°/Myr (Steinberger & Torsvik 2008)
   - To move the pole 30.3° would require ~30 million years
   - This is well outside human timescales
   - **Check:** Are there any proposed rapid TPW events (excursions) in the Quaternary? Review Kirschvink 1997 and more recent literature.
   - **Conclusion expected:** The Great Circle was NOT Earth's equator at any point in human history. Document this rigorously — it's a common fringe claim worth preemptively debunking.

3. **Obliquity variation:**
   - Earth's axial tilt oscillates between 22.1° and 24.5° with a ~41,000-year period
   - Current tilt: ~23.44°; at 10,000 BP: ~24.2°
   - This changes where the tropics and arctic circles are, not the equator
   - But compute: does the circle's inclination (30.3° from equator) match any obliquity-related angle at any epoch? (e.g., 2× obliquity, obliquity + some factor)
   - This is exploratory — likely null but worth documenting.

### Output
- `frozen_equator_analysis.json` — angular distances, TPW timescales, obliquity comparison
- Clear narrative conclusion in RESULTS.md

---

## Analysis 2: Stellar Horizon Events at Key Sites

### Concept
At each site on the Great Circle, the circle crosses the horizon at a specific azimuth (bearing). At historical epochs, specific bright stars rose or set at specific azimuths. **Question:** At the key sites (Giza, Nazca, Easter Island, Persepolis, Mohenjo-daro), did a prominent star rise or set along the circle's bearing at the epoch of monument construction?

### Data
- **Bright star catalog:** Use the Hipparcos catalog or Yale Bright Star Catalog (magnitude < 3.0 gives ~170 stars)
- **Proper motion:** For stars within 500 light-years, proper motion matters over 5,000 years. Use proper motion corrections from Hipparcos.
- **Precession:** Stellar coordinates shift due to precession. Use standard precession matrices (Lieske 1979 or IAU 2006).

### Method
For each of the 5 key sites:
1. Compute the Great Circle's azimuth (bearing) at the site — two values: the direction along the circle in each direction (A and A+180°)
2. For each target epoch (3000 BCE, 2500 BCE, 2000 BCE, 10000 BCE, 6500 BCE):
   a. Precess all bright star coordinates to that epoch
   b. Compute the rising and setting azimuths of each star at the site's latitude using: `Az = arccos(-sin(δ)/cos(φ))` where δ = star declination, φ = site latitude
   c. Find the closest star to the circle's bearing at rise/set
   d. Record the angular offset and the star's identity/magnitude
3. **Monte Carlo significance:**
   a. Generate 10,000 random great circles
   b. For each, compute bearing at Giza and find closest bright-star rise/set azimuth
   c. Report: what percentile is the Great Circle's best match?
4. **Special attention to:**
   - Sirius (α CMa) — most important star in Egyptian astronomy
   - Orion's Belt stars — associated with Giza in popular literature
   - Vega (α Lyr) — former pole star
   - Canopus (α Car) — second brightest star, important in navigation

### Python Libraries
- `astropy.coordinates` for precession, coordinate transforms
- `astropy.time` for epoch handling
- `skyfield` as alternative for high-precision positional astronomy

### Output
- `stellar_horizon_events.csv` — for each site × epoch × direction: closest star, angular offset, magnitude
- `stellar_significance.json` — Monte Carlo percentile ranking
- `horizon_diagram_giza.png` — visual showing circle bearing vs. star rise/set positions at 2500 BCE

### Success Criteria
- If a bright star (mag < 2) aligns within 1° of the circle's bearing at multiple independent sites at the same epoch → significant astronomical correlation
- If no star aligns better than expected by chance → astronomical hypothesis remains null (consistent with existing findings)

---

## Analysis 3: The Circle's Pole on the Celestial Sphere

### Concept
The Great Circle's pole is at 59.682°N, 138.646°W. If we project this point onto the celestial sphere (treating it as a fixed direction relative to Earth's surface), it traces a circle on the sky due to precession over 26,000 years. **Question:** At any epoch, does this projected point coincide with a significant celestial object or location (a bright star, the north ecliptic pole, the north galactic pole)?

### Method
1. Convert the circle's pole (geographic) to celestial coordinates at each epoch:
   a. At epoch T, compute the Local Sidereal Time at the pole's longitude
   b. The pole's RA = LST, Dec = latitude = 59.682°N (fixed)
   c. Actually: the geographic pole projects to Dec = geographic latitude, but RA depends on epoch due to Earth's rotation — this gives a *fixed* declination but varying RA
   d. More precisely: the geographic location 59.682°N, 138.646°W has a fixed relationship to the celestial sphere only at a given sidereal time. Over precession cycles, the celestial coordinates of "straight up from this geographic point" change.
   e. Use the standard geographic-to-celestial transform accounting for precession.

2. For epochs from 25,000 BCE to present (in 100-year steps):
   a. Compute the celestial coordinates (RA, Dec) of the zenith at the circle's pole location
   b. Compute angular distance to:
      - North Ecliptic Pole (λ=270°, β=+66.56°, essentially fixed)
      - North Galactic Pole (RA 12h51m, Dec +27.13°, essentially fixed in galactic coordinates)
      - Every bright star (mag < 2) at that epoch
      - The precessional pole itself (the ecliptic pole)
   c. Find the epoch of closest approach to each reference point

3. **Key question:** Does the circle's pole ever coincide with the ecliptic pole? If so, that would mean the circle is the ecliptic at that epoch — which would be remarkable.
   - The ecliptic pole is at Dec ≈ +66.56° in J2000 coordinates
   - The circle's pole is at Dec ≈ +59.68°
   - Angular separation ≈ 7°
   - Does precession ever close this gap? (Precession moves the celestial pole, not the ecliptic pole — so the answer depends on how the circle's pole projects through the precession cycle)

### Output
- `pole_celestial_track.csv` — epoch, RA, Dec of circle's pole projection
- `pole_proximity.json` — closest approach to ecliptic pole, galactic pole, and each bright star, with epoch
- `pole_track_sky_chart.png` — the circle's pole traced on a star chart over 26,000 years

---

## Analysis 4: Galactic Plane Intersection

### Concept
The Milky Way (galactic plane) traces a great circle on the celestial sphere. The Great Circle, projected onto the sky at a given epoch, also traces a great circle. Two great circles on a sphere always intersect at two points. **Question:** At what epoch(s) do these two great circles coincide most closely, and is there any epoch where they are particularly aligned?

### Method
1. The galactic plane has a fixed orientation: the galactic north pole is at RA 12h51m26.3s, Dec +27°07'42" (J2000)
2. The Great Circle projected onto the sky at epoch T has a pole whose celestial coordinates we computed in Analysis 3
3. Compute the angle between the two great circles' poles at each epoch — this gives the mutual inclination
4. Find the epoch of minimum mutual inclination (maximum alignment)
5. Monte Carlo: for 10,000 random great circles on Earth's surface, compute the same minimum mutual inclination with the galactic plane across all epochs. What percentile is the Great Circle?

### Output
- `galactic_alignment.json` — mutual inclination vs. epoch, minimum epoch, Monte Carlo percentile

---

## Lookahead / Bias Warnings
- Astronomical calculations must use proper precession models (not simplified linear extrapolation)
- "Significance" of a stellar alignment at one site is meaningless without Monte Carlo control — there are ~170 bright stars and the sky is full of them
- The key test is whether the SAME star aligns at MULTIPLE independent sites — that's what separates signal from look-elsewhere effect
- Do NOT cherry-pick epochs. Report the full sweep and let the Monte Carlo determine significance.

---

## Deliverables
1. `outputs/astronomical_epoch/RESULTS.md` — narrative summary
2. All CSVs, JSONs, and figures listed above
3. Clear go/no-go on each sub-analysis for inclusion in future papers
