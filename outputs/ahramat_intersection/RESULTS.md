# AHRAMAT BRANCH PERPENDICULAR INTERSECTION — Results

**Date:** 2026-03-24
**Directive:** 08_ahramat_branch_intersection.md

---

## Analysis 1: Intersection Geometry

The Great Circle crosses the Ahramat Branch at approximately **29.9213°N, 31.1487°E**.

### Intersection Angle (multi-scale)

The directive warned that the branch meanders, so we report the angle at multiple measurement scales:

| Scale | Branch Bearing | GC Bearing | Intersection Angle | Deviation from 90° | Branch Points |
|-------|---------------|-----------|-------------------|-------------------|--------------|
| local_2km | 139.1° | 84.9° | **54.2°** | 35.8° | 1 |
| 5km | 139.1° | 84.9° | **54.2°** | 35.8° | 3 |
| 10km | 144.1° | 84.9° | **59.3°** | 30.7° | 7 |
| 20km | 154.9° | 84.9° | **70.1°** | 19.9° | 12 |
| 50km | 166.2° | 84.9° | **81.3°** | 8.7° | 15 |
| full_branch | 178.8° | 84.9° | **86.1°** | 3.9° | 17 |

**Key finding:** The intersection angle varies significantly with measurement scale due to local meanders:
- **Local (2km):** 54.2° — dominated by a local SSE meander at Zawyet el-Aryan
- **Representative (20km):** 70.1° — captures the Memphis necropolis axis
- **Full branch:** 86.1° — includes the distant Meidum-Faiyum bend

At the 20km scale (the Memphis necropolis pyramid field), the deviation from perpendicular is **19.9°**.

### Pyramid Proximity

- Pyramids within 10km of intersection: **15**
- Pyramids within 20km of intersection: **27**

Nearest pyramids to the intersection point:
- Layer Pyramid: 1.9km (Dynasty 3)
- Unfinished N. Pyramid: 2.2km (Dynasty 4)
- Sahure: 5.9km (Dynasty 5)
- Raneferef: 5.9km (Dynasty 5)
- Neferirkare: 6.0km (Dynasty 5)

---

## Analysis 2: Monte Carlo Probability

**100,000** random great circles tested.

### Test 1: Random Circle × Fixed Branch

| Metric | Count | Probability |
|--------|-------|------------|
| Cross the branch | 471 | 0.0047 |
| Cross near pyramid (<10km) | 460 | 0.00460 |
| Perpendicular (>80°) + near pyramid | 82 | **0.000820** |

### Test 2: Monument Density at Crossing

- Mean pyramids within 10km of random crossing: 5.46
- Observed (Great Circle): **15**
- P(≥observed): **0.16348**

### Test 3: Conditional on Passing Through Egypt

- Circles passing through Egypt's latitude band: 87,363
- Of those, crossing the branch: 471
- Of those, near-perpendicular: 83
- Of those, at dense pyramid segment: 14
- **P(perpendicular + dense | passes Egypt): 0.000160**

### Crossing Angle Distribution

Mean crossing angle for random circles: 58.47° (median: 60.74°).
See `random_circle_angle_histogram.png`.

---

## Analysis 3: Predictive Test — El-khteeb 2025 Saqqara Discoveries

| Metric | New Discoveries | Known Pyramids | Prediction Supported? |
|--------|----------------|---------------|---------------------|
| Mean dist to GC | 6.21 km | 7.16 km | Yes |
| Mean dist to intersection | 8.46 km | 9.36 km | Yes |

**Caveat:** El-khteeb features are all at one location (Saqqara), making this a weak single-site test.

---

## Analysis 4: Crossroads Hypothesis — Global River Crossings

| River | Intersection Angle | Notable Monuments Nearby |
|-------|-------------------|-------------------------|
| Nile (Ahramat Branch) | 70.1° | 15 pyramids within 10km |
| Indus | 83.2° | Mohenjo-daro (214km) |
| Mekong | 58.1° | Angkor Wat (248km) |
| Amazon (Tapajós confluence) | 20.1° | Tapajós earthworks (80km) |
| Ganges | 9.6° | Varanasi (Kashi Vishwanath) (2km) |

**Note:** River bearings are approximate. Monument lists are illustrative, not exhaustive.

---

## Files Generated

- `intersection_geometry.json` — intersection point, multi-scale angles, pyramid distances
- `monte_carlo_intersection.json` — probability analysis (100k trials)
- `new_discoveries_proximity.json` — El-khteeb 2025 proximity test
- `river_crossings.json` — global river crossing analysis
- `intersection_map.png` — annotated map of the Memphis region
- `random_circle_angle_histogram.png` — angle distribution from Monte Carlo
