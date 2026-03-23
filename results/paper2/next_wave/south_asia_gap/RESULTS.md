# Study 3: South Asia Gap Diagnosis

## Question

Why does the great circle show weak signal through northern India/Pakistan,
where the Indus Valley Civilization should produce strong alignment?

## Phase 1: Major IVC Site Distances

| Site | GC Distance (km) | Within 50km | Within 100km |
|------|-------------------|-------------|--------------|
| Mohenjo-daro | 26.2 | Yes | Yes |
| Harappa | 443.2 | No | No |
| Lothal | 443.0 | No | No |
| Dholavira | 348.5 | No | No |
| Kalibangan | 352.1 | No | No |
| Rakhigarhi | 385.6 | No | No |

Reference: Giza ~6.5km, Nazca ~1.0km, Easter Island ~10.3km

## Phase 2: Database Coverage

| Region | Pleiades | p3k14c |
|--------|----------|--------|
| South Asia | 663 | 118 |
| Egypt | 1,287 | 1,412 |
| Anatolia | 4,734 | 2,940 |
| Italy | 7,224 | 3,643 |

## Phase 3: Type Distribution (Pleiades)

| Region | Total | Monumental % | Settlement % |
|--------|-------|-------------|-------------|
| South Asia | 663 | 5.4% | 50.2% |
| Egypt | 1,287 | 21.6% | 47.9% |
| Anatolia | 4,734 | 11.1% | 51.5% |
| Italy | 7,224 | 19.5% | 53.0% |

## Phase 4: Divergence Analysis

### Pleiades + p3k14c

- Monuments: 36 total, 0 within 50km
- Settlements: 451 total, 5 within 50km
- Monument Z-score: -0.613
- Settlement Z-score: +0.613
- **Divergence: -1.227**

### Merged South Asia Dataset

- Monuments: 1613 total, 42 within 50km
- Settlements: 3382 total, 124 within 50km
- **Divergence: -3.856**

## Harappan KMZ Analysis

- Total sites: 1,744
- Within 50km of GC: 49 (2.8%)
- Median distance: 392 km
- Enrichment factor (50km): **2.01x**

Nearest Harappan sites to GC:
| Site | Distance (km) | Period |
|------|--------------|--------|
| Jhukar | 0.0 | Harappan |
| Jhukar | 0.0 | Late_Harappan |
| Naru Waro Dharo | 3.5 | Harappan |
| Lalanji Mari | 4.0 | Early_Harappan |
| Lalanji Mari | 4.0 | Harappan |
| Tupi | 4.3 | Early_Harappan |
| Tupi | 4.3 | Harappan |
| Kander Bhit | 4.3 | Harappan |
| Sabharo | 5.2 | Harappan |
| Zayak North | 9.7 | Early_Harappan |

## Wikipedia IVC Sites

- Total: 52
- Within 50km: 4
- Median distance: 397 km

## Phase 5: GC Segment Through South Asia

- Latitude range: 17.05 to 29.59
- Longitude range: 55.19 to 99.73

## Diagnosis

**Verdict: EXPLAINED**

Major IVC sites are mostly far from the great circle. The GC passes through the region but misses the main centers. The Harappan KMZ dataset shows enrichment near the GC, suggesting the corridor does intersect IVC territory.

### Contributing Factors

- **MAJOR_SITES_FAR**: Only 1/6 within 50km, 1/6 within 100km
- **COVERAGE_ADEQUATE**: SA/Egypt ratio = 0.52
- **LOW_MONUMENT_FRACTION**: SA monumental = 5.4%, Egypt = 21.6%
- **NO_DIVERGENCE**: Monument-settlement divergence = -1.23
- **HARAPPAN_ENRICHED**: Harappan sites 2.0x enriched near GC

## Plots

- `ivc_distances.png` — Major IVC site distances + Harappan histogram
- `coverage_comparison.png` — Database coverage and type distribution
- `gc_south_asia_map.png` — GC segment through South Asia with sites
