# Directive 09: Polynesian Navigation & Sacred Site Alignment — RESULTS

**Date:** 2026-03-24
**Status:** Complete
**Scripts:** 00–05 (6 scripts, all ran successfully)

---

## Executive Summary

The Great Circle passes directly through Easter Island (9.6 km from island center) but is otherwise remarkably isolated from the rest of Polynesia. Only Rapa Iti (286 km) falls within 500 km; the next closest group (Tonga) is 503 km away. The major Polynesian cultural centers — Society Islands, Marquesas, Hawaii — are 1,300–5,500 km from the circle.

**Key findings:**

| Analysis | Result | Significance |
|----------|--------|-------------|
| Marae enrichment (50 km) | 18 sites, 36.2x enrichment | p = 0.011 * |
| Marae enrichment (500 km) | 21 sites, 4.4x enrichment | p = 0.033 * |
| Navigation route alignment | Fiji-Tonga route within 6.3° | 1 of 14 routes |
| Ahu clustering (Rayleigh) | Bimodal clustering at 49° axial | p = 0.038 * |
| Ahu clustering (V-test at GC) | V = 0.25 | p = 0.065 (marginal) |
| Colonization wavefront | ρ = 0.62 (farther = later) | p = 0.032 * |
| Star compass at Easter Island | GC bearing = 75.7° (ENE/WSW) | Hokule'a/Arcturus sector |

---

## Analysis 1: Island Distances & Marae Enrichment

### Island Distances

The Great Circle's Pacific crossing is dominated by open ocean. Easter Island is the only major island group within 200 km:

| Island Group | Distance to GC |
|---|---|
| Easter Island | 9.6 km |
| Rapa Iti | 286 km |
| Tonga | 503 km |
| Pitcairn | 551 km |
| Fiji | 558 km |
| Society Islands | 1,386 km |
| Marquesas | 2,322 km |
| Hawaii | 5,461 km |

### Marae Enrichment (Monte Carlo)

All 18 Easter Island ahu fall within 50 km of the circle (the island is ~15 km from center to coast). Against 10,000 random Pacific-crossing great circles:

- **50 km band:** 18 actual vs 0.5 ± 2.1 expected → **36.2x enrichment** (p = 0.011)
- **100 km band:** 18 vs 1.0 ± 3.0 → **18.8x** (p = 0.023)
- **200 km band:** 18 vs 2.0 ± 4.3 → **9.2x** (p = 0.049)
- **500 km band:** 21 vs 4.8 ± 6.3 → **4.4x** (p = 0.033) — adds 3 Tongan sites

**Caveat:** The enrichment is almost entirely driven by Easter Island. The 18 ahu on an island that sits directly on the circle dominate all narrow-band counts. The 500 km band adds only 3 Tongan sites. This is a "one island drives the signal" situation — the circle happens to pass through the most monument-dense island in the Pacific, but this is essentially one geographic event (the circle intersects Easter Island), not a distributed pattern of alignment across Polynesia.

---

## Analysis 2: Navigation Route Geometry

### Great Circle vs Rhumb Line

All 14 voyaging routes are better fit by great circle arcs than by rhumb lines (constant bearing), confirming that Polynesian navigators followed great circle paths — as expected for navigators using star courses.

### Route Alignment with THE Great Circle

Most routes are not aligned with the Great Circle. Angular separations from the circle:

| Route | Angular Sep | Notes |
|---|---|---|
| Fiji-Tonga | **6.3°** | Within 10° — closest alignment |
| Samoa-Society Islands | 11.3° | Moderate alignment |
| W→E Polynesia (simulated) | 13.2° | |
| Society Is.-Mangareva-Easter Is. | 20.0° | Route to Easter Island |
| Society Is.-Cook Is. | 31.1° | |
| Marquesas-Easter Island | 36.6° | |
| Hawaii-Tahiti routes | 72–74° | Nearly perpendicular |

The **Fiji-Tonga** route (6.3°) is the closest match, but this is the shortest route (1,024 km) and connects two islands that happen to both be ~500 km from the circle. The **Society Islands-Mangareva-Easter Island** route (20.0°) is more interesting — this is the reconstructed colonization route to Easter Island and has moderate angular alignment.

No routes are within 5° of the Great Circle.

### Star Compass Analysis

At Easter Island, the Great Circle bearing is **75.7° / 255.7°** (ENE-WSW). In the traditional Polynesian star compass, this corresponds to the **Hokule'a (Arcturus) rising sector** — the star that is the zenith star of Hawaii and namesake of the famous voyaging canoe. This is a notable cultural connection, though Arcturus is one of the most prominent navigational stars generally.

At other Pacific locations, the circle bearing maps to:
- **Tahiti/Raiatea:** ~96° — Hikianalia (Spica) rising sector (roughly E-W)
- **Marquesas:** ~91° — nearly due E-W
- **Tongatapu:** ~107° — ESE-WNW

---

## Analysis 3: Easter Island Ahu Distribution

### Ahu Positions Relative to GC Crossing

The Great Circle crosses Easter Island on an ENE-WSW axis (75.7° / 255.7°). Ahu are distributed unevenly around the island's perimeter:

- **GC crossing sectors** (60-90° and 240-270°): 6 ahu → **2.0x enrichment** over uniform expectation
- The 240-270° sector (WSW, near Ahu Akivi / Ahu Vai Teka / Tahai complex) has the highest density of any sector (5 ahu)

### Statistical Tests

- **Rayleigh test (bimodal):** R̄ = 0.423, Z = 3.23, **p = 0.038** — significant clustering with mean axial direction 49°
- **V-test (directed at GC bearing):** V = 0.252, u = 1.51, **p = 0.065** — marginally significant
- **Monte Carlo (random crossings):** mean offset 36.3° vs expected 44.8° ± 10.9° → p = 0.307, **not significant**

The Rayleigh test detects significant bimodal clustering, but the preferred axis (49°) is offset ~27° from the GC crossing axis (75.7°). The V-test is marginal. The Monte Carlo test shows the clustering pattern is not unusual compared to random crossing directions.

**Interpretation:** Ahu show non-uniform coastal distribution (clustering in the WSW-W and E-SE sectors), but this pattern is more likely driven by geological factors (suitable coastal platforms, protected bays) than by the Great Circle. Easter Island's coastline is not uniform — the south and southwest coasts have more accessible flat terrain for ahu construction.

---

## Analysis 4: Colonization Wavefront

### Spearman Correlation

Distance to the Great Circle vs colonization date: **ρ = 0.62, p = 0.032** — significant positive correlation, meaning islands farther from the circle were colonized later.

Against 10,000 random great circles: actual ρ is at the **95.3rd percentile** (MC p = 0.048).

### East Polynesia Only

Restricting to post-AD 1000 (removing the West Polynesian early colonization): ρ = 0.60, **p = 0.086** — not significant. The signal weakens substantially when the early Fiji/Tonga/Samoa dates are removed.

### Interpretation

The positive correlation is driven primarily by:
1. **West Polynesia** (Fiji, Tonga, Samoa) — colonized ~1000-800 BCE, at moderate distances (500-1300 km)
2. **Hawaii** — colonized ~AD 1245, at extreme distance (5,461 km)
3. **Easter Island** — colonized ~AD 1230, directly on the circle (but LATE, not early)

The pattern is: the earliest settlements are at moderate distances, not on the circle. Easter Island, despite being on the circle, was among the last islands colonized (~AD 1230). This is opposite to what a "circle as ancient migration corridor" hypothesis would predict. The significant correlation appears to reflect the geographic reality that West Polynesia (near the circle) was settled first and the most remote islands (Hawaii, New Zealand) were settled last — a pattern driven by accessibility and distance from the western Pacific origin, not by alignment with the Great Circle.

---

## Deliverables

| File | Description |
|---|---|
| `00_polynesian_site_data.py` | Master data module (41 sacred sites, 20 island groups, 14 routes, 12 colonization entries) |
| `01_island_distances.py` | Island distance computation and Pacific map |
| `02_marae_enrichment.py` | Monte Carlo enrichment analysis |
| `03_navigation_routes.py` | Route geometry, GC comparison, star compass |
| `04_ahu_orientation.py` | Ahu distribution analysis with circular statistics |
| `05_colonization_wavefront.py` | Colonization timeline correlation |
| `island_distances.csv` | 20 island groups with distances |
| `island_distances.json` | Full distance analysis results |
| `pacific_circle_map.png` | Great Circle across Pacific with island groups |
| `marae_enrichment.json` | Monte Carlo enrichment at 4 thresholds |
| `route_geometries.csv` | 14 routes with GC fits |
| `route_circle_comparison.json` | Route-to-circle angular separations |
| `star_compass_bearing.json` | GC bearing mapped to star compass at 7 locations |
| `ahu_orientation_test.json` | Rayleigh/V-test results |
| `ahu_density_by_bearing.json` | Ahu count per 30° sector |
| `easter_island_circle_overlay.png` | Easter Island map with ahu and GC path |
| `colonization_vs_distance.json` | Spearman correlation and MC results |
| `colonization_timeline.png` | Date vs distance scatter with MC distribution |

---

## Caveats & Limitations

1. **Database completeness:** Our 41-site dataset is a fraction of the hundreds of known marae across Polynesia. Easter Island's 18 ahu dominate. A complete database would include ~300+ Society Islands marae, ~600+ Hawaiian heiau, etc. The enrichment analysis would change substantially with a fuller dataset.

2. **Easter Island dominance:** The circle passes through Easter Island — the most monument-dense island in the Pacific. This single geographic fact drives nearly all the statistical significance. It is one data point, not a distributed pattern.

3. **Pacific vastness:** The Pacific Ocean is so large that "near the circle" means very different things than in continental contexts. At 500 km threshold, only Easter Island and Rapa Iti qualify. The circle is essentially alone in the eastern Pacific.

4. **Navigation constraints:** Polynesian voyaging was constrained by winds, currents, and island visibility — not pure geometry. Great circle routes are not always the fastest or safest sailing routes. The trade winds, ITCZ, and western Pacific monsoon system dictated real-world routing.

5. **Colonization chronology:** The "short chronology" (Wilmshurst et al. 2011) is used here. Some researchers argue for earlier East Polynesian dates. The West-East "long pause" (~1800 years between West Polynesian settlement and East Polynesian expansion) is the dominant feature in the data.

6. **Star compass significance:** The GC bearing at Easter Island (75.7°, Hokule'a/Arcturus sector) is culturally interesting but not statistically surprising — any bearing at a tropical location will correspond to some named star in the Polynesian compass.

---

## Bottom Line

The Great Circle's relationship with Polynesia is essentially a story about **Easter Island** — the circle passes through the most isolated, most monument-dense island in the Pacific. This produces impressive enrichment statistics but represents a single geographic coincidence, not a distributed Pacific-wide pattern. The broader Polynesian cultural world (Society Islands, Hawaii, Marquesas, Samoa, Tonga) is far from the circle, and documented navigation routes do not preferentially align with it. The colonization chronology shows a weak positive correlation that disappears when restricted to East Polynesia alone. The star compass connection (Hokule'a/Arcturus at Easter Island) is evocative but not statistically distinctive.
