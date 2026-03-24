# Nazca Geoglyph Orientation Analysis — Results

**Date:** 2026-03-24
**Directive:** 06_nazca_geoglyph_orientation.md v1.0
**Data source:** Richter, Teichert & Pavelka 2021, *Applied Sciences* 11(4):1637 — 2,308 straight line azimuths from NascaGIS database

---

## Data Availability

The directive specified three potential data sources:

| Source | Status |
|--------|--------|
| **Sakai et al. 2024** (PNAS, 733 figurative geoglyphs) | **Unavailable** — coordinates restricted by Peru Ministry of Culture; paper contains no orientation data whatsoever |
| **Richter et al. 2021** (MDPI, 2,308 straight lines) | **Partial** — aggregate histogram available from Figure 13; individual azimuths in NascaGIS database (not public) |
| **Aveni 1990** (ray centers, ~800 lines) | **Print only** — bimodal peaks at ~10° and ~80° cited in secondary literature; exact values not digitally accessible |
| **OpenStreetMap** | Only 3 geoglyph polygons found — insufficient |

Analysis proceeds with the Richter histogram (2,308 lines, 10° bins).

---

## Analysis 1: Great Circle Bearing at Nazca

The Great Circle (pole: 59.682°N, 138.646°W) bearing at the Nazca Pampa:

| Location | Bearing (E) | Bearing (W) | Distance to Circle |
|----------|-------------|-------------|-------------------|
| Pampa center (-14.735, -75.130) | **63.14°** | **243.14°** | 11.8 km |
| Pampa NE (-14.60, -75.00) | 63.11° | 243.11° | 4.7 km |
| Hummingbird geoglyph (-14.692, -75.149) | 63.14° | 243.14° | 6.6 km |
| Cahuachi temple (-14.817, -75.120) | 63.14° | 243.14° | 20.4 km |

- **Bearing:** 63.14° ± 0.02° (ENE-WSW; 26.9° north of due east)
- **Direction:** Roughly NE–SW
- The bearing is essentially constant across the 50km Pampa
- The circle passes within **4.7 km** of the NE corner of the Pampa

---

## Analysis 2: Geoglyph Orientations vs Circle Bearing

### 2a. Bin-level comparison

The GC bearing (63.14°) falls in the **60–70° bin**, which contains **177 lines** — the **most populated bin** of all 18 bins (tied with 155°).

| Metric | Value |
|--------|-------|
| GC bin count | 177 lines |
| Uniform expectation | 128.2 lines/bin |
| Ratio to expected | **1.38×** |
| Bin rank | **1 of 18** (tied) |
| Chi² (this bin alone) | 18.56 |

The overall distribution is **significantly non-uniform** (χ² = 133.04, df = 17, p < 10⁻⁶).

### 2b. Full distribution

```
  Bin     Count   Ratio
   5°:     91    0.71×
  15°:    126    0.98×
  25°:    103    0.80×
  35°:    100    0.78×
  45°:    153    1.19×
  55°:    127    0.99×
  65°:    177    1.38×  ← GC bearing + June solstice sunrise
  75°:    116    0.90×
  85°:    129    1.01×
  95°:     76    0.59×  (trough)
 105°:    104    0.81×
 115°:     86    0.67×  (trough)
 125°:    165    1.29×
 135°:    140    1.09×
 145°:    172    1.34×
 155°:    177    1.38×  ← 2nd peak (tied)
 165°:    124    0.97×
 175°:    142    1.11×
```

### 2c. Circular statistics

**Rayleigh test (axial, doubled angles):**
- Mean resultant length R̄ = 0.0399
- Mean preferred direction = **154.9°** (NNW–SSE)
- Rayleigh Z = 3.67, p = 0.025
- → The distribution has a **significant** preferred direction, but it points at **155°, not 63°**

**V-test (directed, against GC bearing):**
- V = -0.0398, u = -2.71, p = 0.997
- → Lines do **NOT** cluster toward the GC bearing
- The negative V-statistic means lines actually point *away* from the GC direction (toward the perpendicular)

### 2d. Monte Carlo comparison

10,000 random test bearings:
- GC V-statistic percentile: **1.9%** (lower is worse)
- The GC bearing is one of the *least* aligned directions relative to the overall distribution's mean
- Bin density: 10.8% of random bearings fall in bins ≥ as dense as the GC bin (so the GC bin is in the top ~11% by density — notable but far from unique)

---

## Critical Confound: Solar Alignment

The Great Circle bearing at Nazca (**63.14°**) is only **2.5°** from the June solstice sunrise azimuth at this latitude (**65.6°**).

| Reference Azimuth | Value | Δ from GC | Bin Count | Rank |
|---|---|---|---|---|
| **Great Circle** | 63.1° | — | 177 | 1/18 |
| **June solstice sunrise** | 65.6° | 2.5° | 177 | 1/18 |
| Equinox sunrise | 90.0° | 26.9° | 76 | 18/18 |
| December solstice sunrise | 114.4° | 51.3° | 86 | 17/18 |

**The GC bearing and the June solstice sunrise fall in the same 10° bin.** Any observed correlation between Nazca line orientations and the GC bearing cannot be distinguished from a solar alignment hypothesis. In fact, the Richter et al. paper specifically investigates solar alignment at Nazca and finds the ~65° bin enriched.

---

## Analysis 3: Spatial Gradient

**Blocked.** Requires individual geoglyph coordinates with orientations. The Sakai dataset is restricted and the NascaGIS database is not publicly accessible. Cannot test whether alignment improves closer to the circle.

---

## Analysis 4: Ray Centers

**Blocked.** Aveni's ray center coordinates are in print-only publications (1990 book). Only 3 Nazca features found in OpenStreetMap. Cannot test ray center proximity to the circle.

---

## Peak Structure

The distribution has 8 peaks (bins higher than both neighbors):

| Peak | Count | Δ from GC |
|------|-------|-----------|
| 65° | 177 | 1.9° |
| 155° | 177 | 88.1° |
| 125° | 165 | 61.9° |
| 45° | 153 | 18.1° |
| 175° | 142 | 68.1° |
| 85° | 129 | 21.9° |
| 15° | 126 | 48.1° |
| 105° | 104 | 41.9° |

The two highest peaks (65° and 155°) are **90° apart**, suggesting a dominant cross-axis pattern in the Nazca line system. The 65° direction is the June solstice sunrise; the 155° direction has no obvious astronomical correlate but may relate to topography or hydrology (water flows roughly NW→SE across the Pampa).

---

## Interpretation

### What was found:
1. The GC bearing at Nazca is **63.14°** (NE-SW), and the circle passes within **4.7 km** of the Pampa
2. The 60-70° bin is the **most populated** orientation bin for Nazca lines (177 of 2,308, 1.38× uniform expectation)
3. The distribution is significantly non-uniform (p < 10⁻⁶)

### Why this is NOT evidence for a Great Circle correlation:
1. **Solar confound:** The GC bearing is indistinguishable from the June solstice sunrise (Δ = 2.5°). The enrichment at ~65° has a straightforward astronomical explanation that was proposed decades ago
2. **Wrong mean direction:** The overall preferred direction is 155° (Rayleigh test), roughly perpendicular to the GC bearing. The V-test against the GC bearing gives p = 0.997 — actively anti-correlated
3. **Multi-peak structure:** The distribution has many peaks; 65° is one of two co-dominant orientations. The pattern is better explained by a 65°/155° cross-axis (solstice sunrise + topographic drainage) than by a single preferred direction

### What would change this assessment:
- Access to the Sakai figurative geoglyph orientations (restricted) — if the *figurative* geoglyphs show a different pattern from the *linear* features
- Individual geoglyph coordinates to test the spatial gradient (Analysis 3) — do features *closer* to the circle show stronger alignment?
- Access to ray center coordinates to test whether they cluster on the circle (Analysis 4)

---

## Files

| File | Description |
|------|-------------|
| `01_circle_bearing.py` | Analysis 1: GC bearing computation |
| `circle_bearing_at_nazca.json` | Bearing results at 7 reference points |
| `03_orientation_analysis.py` | Analysis 2: full statistical tests |
| `orientation_analysis.json` | All test results |
| `04_rose_diagram.py` | Visualization code |
| `nazca_orientation_rose_diagram.png` | Rose diagram (full 360°) |
| `nazca_orientation_histogram.png` | Axial histogram (0-180°) |

---

## Data Sources

- Richter, C., Teichert, B., & Pavelka, K. (2021). Astronomical investigation to verify the calendar theory of the Nasca lines. *Applied Sciences*, 11(4), 1637. doi:10.3390/app11041637
- Sakai, M., et al. (2024). AI-accelerated Nazca survey. *PNAS*, 121(40), e2407652121. doi:10.1073/pnas.2407652121 [data restricted]
- Aveni, A.F. (1990). Order in the Nazca lines. In *The Lines of Nazca*. American Philosophical Society. [print only]
- Great Circle pole: (59.682122°N, 138.646087°W) per Jim Alison
