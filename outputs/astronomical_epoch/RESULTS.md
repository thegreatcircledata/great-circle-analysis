# Astronomical Epoch Analysis — Results

**Date:** 2026-03-24
**Directive:** 02_astronomical_epoch_analysis.md
**Analyst:** Claude (for Ell)

---

## Executive Summary

Four independent astronomical hypotheses were tested for the Great Circle (pole: 59.682°N, 138.646°W). **All four returned null or near-null results.** The Great Circle has no detectable astronomical significance at any historical epoch.

| Analysis | Core Question | Result | Verdict |
|----------|--------------|--------|---------|
| 1. Frozen Equator | Was the GC ever Earth's equator? | Requires 30 Myr of TPW; impossible in human timescales | **Definitively null** |
| 2. Stellar Horizons | Do bright stars rise/set along the GC at key sites? | Best match (Procyon at Giza, 0.1°) is at 11th percentile vs. random | **Null** |
| 3. Pole on Celestial Sphere | Does the GC pole coincide with a significant sky point? | Closest approach to ecliptic pole: 18.7°; no coincidence | **Null** |
| 4. Galactic Plane | Does the GC align with the Milky Way at any epoch? | Min inclination 37.8°; 91.8th percentile (less aligned than average) | **Null** |

---

## Analysis 1: The Circle as a Frozen Equator

### Question
At what epoch (if any) would the Great Circle have been Earth's equator?

### Results

**Axial precession** changes the orientation of Earth's spin axis relative to the stars but does NOT move the geographic pole on Earth's surface. The North Pole stays at 90°N geographic regardless of precession epoch. Therefore, precession alone **cannot** make the Great Circle into the equator — ever.

**True polar wander (TPW)** does move the geographic pole relative to the mantle, but at rates of ~1°/Myr (Steinberger & Torsvik 2008). The Great Circle's pole is 30.3° from the geographic pole, requiring ~30 million years of TPW at typical rates. Even at the fastest proposed Quaternary rates (5°/Myr), this would take ~6 Myr. Anatomically modern humans have existed for ~300,000 years, during which the pole has moved <0.3°.

**Obliquity variation** (22.1° to 24.5° over 41,000 years) changes the tropics and arctic circles but not the equator. The GC's inclination (30.3°) does not match any simple obliquity-related angle; closest match is 2× max obliquity (off by 5.8°), with no physical mechanism to link them.

### Conclusion
**The Great Circle was NOT Earth's equator at any point in human history.** The "ancient equator" claim, which circulates widely online, is definitively refuted by basic geophysics. Axial precession doesn't move the equator; TPW is far too slow; obliquity doesn't apply.

### Files
- `frozen_equator_analysis.json`

---

## Analysis 2: Stellar Horizon Events at Key Sites

### Question
At the key sites (Giza, Nazca, Easter Island, Persepolis, Mohenjo-daro, Angkor Wat, Machu Picchu, Ollantaytambo), did a prominent star rise or set along the Great Circle's bearing at the epoch of monument construction?

### Method
For each site, the Great Circle's azimuth (bearing) was computed. For each target epoch, all bright stars (mag < 2.5, ~46 stars) were precessed to that epoch using astropy, and their rising/setting azimuths were computed. The closest match was recorded. A 10,000-trial Monte Carlo test was run at Giza (2500 BCE) using random great circles.

### Results — Individual Sites

| Site | Epoch | Best Star | Direction | Offset |
|------|-------|-----------|-----------|--------|
| Giza | 2500 BCE | **Procyon** | setting | **0.1°** |
| Giza | 2600 BCE | Procyon | setting | 0.2° |
| Giza | 3000 BCE | Alhena | setting | 0.5° |
| Giza | 10000 BCE | Menkalinan | rising | 1.1° |
| Nazca | 200 BCE | Nunki | rising | 0.4° |
| Nazca | 500 BCE | Wezen | rising | 0.6° |
| Easter Island | 1000 CE | Aldebaran | setting | 1.4° |
| Persepolis | 500 BCE | Mintaka | setting | 0.7° |
| Mohenjo-daro | 2600 BCE | Alnath | rising | 0.6° |
| Angkor Wat | 1150 CE | Wezen | setting | 0.5° |
| Machu Picchu | 1450 CE | Nunki | rising | 0.1° |
| Ollantaytambo | 1450 CE | Nunki | rising | 0.0° |

The GC bearings at Giza are 84.9° and 264.9° — essentially due east and west with a slight northward tilt. This means the circle's bearing at Giza is close to the equinox sunrise/sunset direction, which is where many stars rise and set. The Procyon alignment at 0.1° looks impressive in isolation.

### Monte Carlo Significance

| Metric | Value |
|--------|-------|
| GC best offset at Giza (2500 BCE) | 0.13° |
| Random GC mean best offset | 2.41° ± 4.15° |
| Random GC median best offset | 0.76° |
| **GC percentile** | **11.1%** |

The Great Circle's best stellar alignment at Giza is at the 11th percentile — better than average, but not significant by any conventional threshold (5% or 1%). Roughly 1 in 9 random great circles achieves an equally good or better stellar alignment at Giza.

### Multi-Site Consistency

Stars aligning within 5° at 2+ sites:
- **Nunki** (σ Sgr, mag 2.0): 3 sites — Nazca, Machu Picchu, Ollantaytambo. However, these three sites are geographically clustered in Peru (within ~200 km), so they are not independent — the GC bearing is nearly identical at all three (~62.5°).
- **Wezen** (δ CMa, mag 1.8): Nazca and Angkor Wat — but at different epochs (500 BCE vs. 1150 CE), and the same star at different epochs is expected given how slowly bright star declinations change.
- **Procyon** (α CMi, mag 0.3): Giza and Persepolis — at different epochs (2500 BCE vs. 2000 BCE), offset 0.1° and 1.0° respectively.

**No star aligns within 1° at multiple truly independent sites at the same epoch.** The multi-site test — the one that would separate signal from look-elsewhere noise — is null.

### Conclusion
Individual tight alignments are not unusual given ~46 bright stars × 2 rise/set azimuths × 2 GC bearings = ~184 opportunities per site-epoch. The Monte Carlo confirms the GC is not an outlier. The multi-site consistency check finds no common stellar alignment across independent sites. **Astronomical hypothesis: null.**

### Files
- `stellar_horizon_events.csv` — all site × epoch results
- `stellar_significance.json` — Monte Carlo percentiles
- `stellar_horizon_full.json` — complete results with multi-site analysis

---

## Analysis 3: The Circle's Pole on the Celestial Sphere

### Question
Does the GC's pole, projected onto the celestial sphere, coincide with any significant celestial object at any epoch within the precession cycle?

### Method
The GC pole (geographic: 59.682°N, 138.646°W) was projected to celestial coordinates at each epoch from 23,000 BCE to 3,000 CE in 100-year steps using astropy's precession model. Angular separations to the North Ecliptic Pole, North Galactic Pole, and bright reference stars were computed.

### Results

**Ecliptic Pole:** Closest approach = **18.7°** at epoch 19,900 BCE. The GC pole never comes close to the ecliptic pole. The Great Circle is NOT a projection of the ecliptic at any epoch.

**Galactic Pole:** Closest approach = **37.8°** at epoch 1,600 CE. No alignment.

**Bright stars (closest approaches):**

| Star | Closest Approach | Epoch |
|------|-----------------|-------|
| Thuban (α Dra) | **1.7°** | 3000 CE |
| Vega (α Lyr) | **4.7°** | 6200 BCE |
| Polaris (α UMi) | 4.7° | 18,400 BCE |
| Deneb (α Cyg) | 6.4° | 10,100 BCE |

Thuban, the pole star of ancient Egypt (~2700 BCE), comes within 1.7° of the GC pole projection at 3000 CE — not at the epoch when Thuban was the pole star. Vega, which was a rough pole star ~12,000 BCE, comes within 4.7° at 6200 BCE. These are mild curiosities but not tight enough alignments to be significant, especially given that bright stars are distributed across the sky.

### Conclusion
The GC pole does not coincide with any significant celestial reference point (ecliptic pole, galactic pole, or bright star) at any epoch within the precession cycle. **Null result.**

### Files
- `pole_celestial_track.csv` — 261 epoch positions
- `pole_proximity.json` — closest approaches to all reference points

---

## Analysis 4: Galactic Plane Intersection

### Question
At what epoch does the Great Circle, projected onto the sky, most closely align with the Milky Way's plane? Is this alignment unusual compared to random great circles?

### Method
The mutual inclination between the GC (projected via its pole's celestial coordinates) and the galactic plane was computed at 100-year steps from 23,000 BCE to 3,000 CE. A 10,000-trial Monte Carlo test compared the GC's minimum inclination against random great circles.

### Results

| Metric | Value |
|--------|-------|
| Minimum mutual inclination | **37.80°** at epoch 1500 CE |
| Maximum mutual inclination | 83.44° at epoch 11,800 BCE |
| Current (2000 CE) | 37.96° |
| MC mean minimum inclination | 18.23° ± 12.21° |
| MC median | 16.48° |
| **GC percentile** | **91.8%** |

The Great Circle is *less* aligned with the galactic plane than the typical random great circle. Its minimum inclination of 37.8° is well above the Monte Carlo median of 16.5°. This means the GC is actually anti-correlated with galactic alignment — it avoids the Milky Way more than a random circle would.

The inclination oscillates between 37.8° and 83.4° over the precession cycle, never approaching 0° (perfect alignment) or even the random-circle median.

### Conclusion
**Null — emphatically so.** The Great Circle is less aligned with the galactic plane than ~92% of random great circles. There is no galactic alignment at any epoch.

### Files
- `galactic_alignment.json` — full results with Monte Carlo
- `galactic_alignment_timeseries.csv` — inclination vs. epoch

---

## Overall Conclusions

1. **The "ancient equator" hypothesis is physically impossible.** Precession doesn't move the equator; TPW is 5–6 orders of magnitude too slow. This is the most important finding for public communication.

2. **No stellar alignment survives Monte Carlo testing.** Individual site-epoch matches look impressive (Procyon at Giza: 0.1°!) but are expected by chance given the density of bright stars and the degrees of freedom available.

3. **The GC pole traces an unremarkable path on the celestial sphere.** It doesn't approach the ecliptic pole (18.7° minimum separation), the galactic pole (37.8°), or any very bright star at a historically meaningful epoch.

4. **The existing paper's null result (Rayleigh p = 0.75 for monument orientations) is reinforced by four additional null tests.** The Great Circle's relationship to the sky is thoroughly non-special.

### Implications for the Research Program
The astronomical hypotheses were always long shots, and these null results strengthen the paper's core argument: the Great Circle's significance lies in its *geographical* relationship to ancient monuments, not in any astronomical encoding. The frozen equator debunking is valuable public-facing content. The stellar horizon analysis methodology is sound and worth documenting even with null results — it shows the work was done rigorously.

### Suggested Next Steps
- **Include the frozen equator debunking** in the paper as a supplementary note — it preempts a common objection
- **Archive the stellar horizon data** — if new sites are added to the circle, the framework can quickly test them
- Move to Directive 03 (Geomagnetic) or 04 (Temporal Propagation) — both are likely to produce more interesting results
