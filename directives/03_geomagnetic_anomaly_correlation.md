# GEOMAGNETIC ANOMALY CORRELATION — Research Directive v1.0

**Date:** 2026-03-23
**Author:** Claude (for Ell)
**Objective:** Test whether the Great Circle correlates with Earth's geomagnetic field — either the modern field, the historical field (past 10,000 years), or the paleomagnetic field — and additionally test gravity anomaly correlation as a quick secondary analysis.
**Repository:** great-circle-analysis/
**Output Directory:** outputs/geomagnetic_correlation/
**Estimated Runtime:** 1–2 hours
**Parallel Safe:** Yes

---

## Background & Motivation

The geological correlation test (already completed) showed the circle doesn't follow plate boundaries, magnetic lineaments, or other geophysical features. But that test used SRTM elevation data and qualitative assessment. It did NOT quantitatively test the geomagnetic field itself. Some fringe hypotheses (and some serious ones — see Kirschvink's magnetoreception work) propose that ancient peoples were sensitive to magnetic fields or that sacred sites correlate with geomagnetic anomalies. This has never been tested for the Great Circle with actual geomagnetic data.

Additionally, Earth's gravitational field varies measurably. The GRACE satellite data provides detailed gravity anomaly maps. A correlation with gravity gradients would be unexpected but is trivial to test and ruling it out adds rigor.

---

## Analysis 1: Modern Geomagnetic Field Correlation

### Data
- **IGRF-13** (International Geomagnetic Reference Field, 13th generation)
  - Python package: `pyIGRF` or compute from coefficients
  - Provides: total intensity (F), declination (D), inclination (I), horizontal intensity (H) at any lat/lon/date
- **EMAG2v3** (Earth Magnetic Anomaly Grid) — 2-arc-minute resolution crustal magnetic anomaly map
  - Download: https://www.ncei.noaa.gov/products/earth-magnetic-model-anomaly-grid-2-arc-minute
  - This captures crustal anomalies (local geological magnetism), distinct from the core field

### Method

#### A. Core Field Test
1. Sample the Great Circle at 1° intervals (360 points)
2. At each point, compute IGRF field components (F, D, I, H) for epoch 2025
3. Compute the same components at 360 random points at the same latitude distribution (to control for latitude dependence of the geomagnetic field)
4. Test: does the circle follow a contour of constant field intensity, declination, or inclination?
   - Compute the variance of each field component along the circle
   - Monte Carlo: compute variance along 10,000 random great circles
   - If the Great Circle has unusually LOW variance in any component → it follows a magnetic contour
5. Repeat at epoch 3000 BCE using the CALS10k.2 model (see Analysis 2)

#### B. Crustal Anomaly Test
1. Extract EMAG2v3 anomaly values at 1° intervals along the circle
2. Compute mean absolute anomaly on-corridor (within 50km) vs. global mean
3. Test: does the circle follow a crustal magnetic anomaly ridge or valley?
4. Specifically: at the 6 cluster locations, are there magnetic anomalies? (Some studies link sacred sites to local magnetic anomalies)
5. Monte Carlo: 10,000 random circles, compare anomaly profiles

### Output
- `core_field_profile.csv` — field components along circle
- `crustal_anomaly_profile.csv` — EMAG2 values along circle
- `geomag_variance_test.json` — Monte Carlo percentile for field variance
- `geomag_cluster_anomalies.json` — anomaly values at the 6 cluster sites

---

## Analysis 2: Historical Geomagnetic Field (CALS10k.2)

### Data
- **CALS10k.2** (Constable et al. 2016, doi:10.1016/j.epsl.2016.06.005)
  - Spherical harmonic model of Earth's magnetic field for the past 10,000 years
  - Python implementation: `pymagsv` package, or download coefficients from:
    https://earthref.org/ERDA/2207/
  - Alternative: **SHAWQ-Iron Age** model for 1000 BCE–0 CE (Osete et al. 2020) for higher resolution in the critical period

### Method
1. For epochs: 8000 BCE, 6000 BCE, 4000 BCE, 3000 BCE, 2500 BCE, 2000 BCE, 1000 BCE, 0 CE
   a. Compute the geomagnetic field (F, D, I) at each 1° sample along the Great Circle
   b. Compute the position of the geomagnetic north pole at that epoch (the dipole axis)
   c. Compute the angular distance between the Great Circle's pole and the geomagnetic pole
   d. Compute the angular distance between the Great Circle and the geomagnetic equator
2. **Key question:** At any epoch, does the Great Circle's pole coincide with the geomagnetic pole? Or does the circle approximate the magnetic equator?
   - The geomagnetic pole wanders significantly over millennia (it was near 80°N, 90°W at 3000 BCE)
   - If the circle ever coincided with the magnetic equator, that would be remarkable
3. Monte Carlo: for 1,000 random great circles, compute minimum angular separation from the magnetic equator across all epochs. Percentile rank the Great Circle.

### Output
- `historical_field_profile.json` — field components along circle at each epoch
- `magnetic_pole_distance.csv` — angular distance from circle's pole to magnetic pole, by epoch
- `magnetic_equator_separation.json` — separation between circle and magnetic equator, by epoch, with Monte Carlo percentile

---

## Analysis 3: Gravity Anomaly Overlay (Quick Test)

### Data
- **EGM2008** (Earth Gravitational Model 2008) — free-air or Bouguer gravity anomaly
  - 2.5-arc-minute resolution
  - Download: http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/
  - Or use `boule` Python package for reference field calculations
- **Alternative shortcut:** Use the GRACE-derived global gravity anomaly grid (lower resolution but easier to work with)
  - GGMplus: http://ddfe.curtin.edu.au/models/GGMplus/

### Method
1. Extract gravity anomaly values at 1° intervals along the Great Circle
2. Compute mean anomaly on-corridor vs. global mean
3. Test: does the circle follow a gravity gradient or anomaly ridge?
4. Monte Carlo: 10,000 random circles, compare gravity anomaly profiles
5. Specifically test: do the 6 cluster locations sit on gravity anomalies?

### Output
- `gravity_anomaly_profile.csv` — gravity values along circle
- `gravity_monte_carlo.json` — percentile ranking vs. random circles

### Expected Result
Almost certainly null — but documenting it adds one more ruled-out explanation.

---

## Deliverables
1. `outputs/geomagnetic_correlation/RESULTS.md` — narrative summary
2. All data files and figures listed above
3. If ANY test shows p < 0.05 → flag for deeper investigation and independent replication
4. If all null → add to the "ruled out" list in the next paper revision
