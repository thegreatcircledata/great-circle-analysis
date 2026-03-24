# ORAL TRADITION & MYTHOLOGY SPATIAL MAPPING — Research Directive v1.0

**Date:** 2026-03-23
**Author:** Claude (for Ell)
**Objective:** Test whether specific mythological motifs — particularly flood myths, creation narratives, sky-deity traditions, and megalith-origin stories — cluster along the Great Circle more than expected by chance.
**Repository:** great-circle-analysis/
**Output Directory:** outputs/oral_tradition_mapping/
**Estimated Runtime:** 3–5 hours (data acquisition is the bottleneck)
**Parallel Safe:** Yes

---

## Background & Motivation

Yuri Berezkin (Museum of Anthropology and Ethnography, St. Petersburg) has compiled the world's largest database of folklore motifs: ~50,000 entries across ~1,000+ ethnic groups, each georeferenced to a homeland region. The database catalogs the presence/absence of ~2,500 distinct narrative motifs across cultures. This is, in principle, the oral tradition equivalent of the Megalithic Portal — a global georeferenced database of cultural phenomena.

If the Great Circle traces an ancient corridor of human activity and cultural exchange, we might expect specific motifs (especially those associated with deep antiquity — flood myths, world-pillar myths, creator-traveler myths) to cluster along it. This has never been tested with spatial statistics.

Additionally, Julien d'Huy (Paris-Sorbonne) has used phylogenetic methods to reconstruct the deep history of specific myths, dating some (like the Cosmic Hunt) to the Paleolithic. His published reconstructions include geographic dispersal patterns that could be overlaid with the Great Circle.

---

## Data Sources

### Primary: Berezkin Analytical Catalogue of Mythological Motifs
- **Access:** http://www.ruthenia.ru/folklore/berezkin/
- The interface is in Russian but navigable; data can be extracted by motif code
- Each motif has: code, name, description, and a list of ethnic groups where it's attested
- Each ethnic group has an approximate geographic center (lat/lon extractable from the associated maps)
- **Alternative access:** Some of Berezkin's data is compiled in his published papers (doi:10.1515/spiram-2015-0009, doi:10.15463/ie1568.2483) with maps that can be digitized
- **If full database is inaccessible:** Use the motif distribution maps from his 2015 book "Cosmic Hunt" and manually digitize ~20 key motifs

### Secondary: d'Huy's Phylogenetic Myth Reconstructions
- Published in: d'Huy 2012 (Cosmic Hunt, doi:10.1515/spiram-2012-0002), d'Huy 2013 (Polyphemus), d'Huy 2016 (Pygmalion)
- These include reconstructed dispersal routes with geographic waypoints

### Tertiary: Witzel's "Origins of the World's Mythologies" (2012)
- Contains proposed "Laurasian" and "Gondwanan" mythology distributions
- The Laurasian mythos has a proposed dispersal route out of Africa through Eurasia — this route can be compared to the Great Circle

---

## Target Motifs

Focus on motifs with proposed deep antiquity and wide geographic distribution:

1. **Flood/Deluge myths** (Berezkin motifs A28–A33 approximately)
   - Universal near-universal; but do they cluster on the circle more than off it?
2. **World-pillar / Axis Mundi** (sky support, world tree, mountain holding up sky)
   - Relevant because the Great Circle itself could be conceived as a "world axis" projected flat
3. **Creator-traveler / Culture hero who travels and creates landmarks**
   - These myths describe a being who walks across the landscape creating features (cf. Aboriginal Dreamtime, Greek Heracles, Polynesian Maui)
   - If the corridor was a travel route, traveler-creator myths should follow it
4. **Star myths / Celestial navigation**
   - Myths specifically about using stars for travel or orientation
5. **Giants / Ancient builders**
   - Myths about an earlier race of beings who built stone structures
6. **Cosmic Hunt** (d'Huy's best-dated myth — possibly Paleolithic)
   - Well-mapped geographic distribution — can directly overlay with circle

---

## Method

### Step 1: Data Compilation
1. From Berezkin's database, extract the geographic distribution of each target motif
2. For each ethnic group: assign a centroid lat/lon (use Ethnologue coordinates or HRAF region centers as a fallback)
3. Build a master table: motif_code | ethnic_group | lat | lon | present (1/0)

### Step 2: Circle Proximity Test (Per Motif)
For each target motif:
1. Compute the distance from each ethnic group centroid to the Great Circle
2. Classify groups as "on-corridor" (centroid ≤200km from circle — wider band because ethnic territories are large) vs. "off-corridor"
3. Compute the proportion of on-corridor groups that have the motif vs. off-corridor
4. **Enrichment ratio:** (on-corridor prevalence) / (off-corridor prevalence)
5. **Monte Carlo:** Generate 10,000 random great circles. For each, compute the enrichment ratio. Percentile rank the Great Circle.
6. **Latitude control:** Since motif distributions are latitude-dependent (tropical vs. temperate), use distribution-matched random circles (matching the Great Circle's latitude profile)

### Step 3: Multi-Motif Composite Test
1. For each ethnic group, count how many of the target motifs are present (score 0–6)
2. Test whether on-corridor groups have a higher mean motif count than off-corridor
3. Monte Carlo significance as above

### Step 4: Dispersal Route Overlay
1. Digitize d'Huy's reconstructed dispersal routes for the Cosmic Hunt and Polyphemus myths
2. Compute the mean distance between each dispersal route and the Great Circle
3. For each route, compute: does it follow the Great Circle segment by segment, or diverge?
4. Visual overlay: Great Circle + myth dispersal routes on a world map

---

## Output
- `motif_distributions.csv` — master table of all motif × ethnic group × lat/lon data
- `enrichment_by_motif.json` — enrichment ratio, Monte Carlo p-value for each target motif
- `composite_motif_score.json` — multi-motif test results
- `dispersal_route_overlay.png` — visual comparison with d'Huy's reconstructions
- `mythology_corridor_map.png` — world map with motif density heat overlay on/near the circle

---

## Caveats & Bias Warnings
- Ethnic group centroids are approximations — some groups span thousands of km
- The 200km band is necessarily wider than the 50km used for archaeological sites; results must be interpreted with this in mind
- Berezkin's database has Euro-centric collection bias (European myths are better documented than, say, Amazonian ones)
- Motif similarity can reflect common human cognition rather than historical contact (flood myths may be universal because floods are universal)
- This analysis is EXPLORATORY — even a positive result would need careful interpretation
- The biggest risk is the "unfalsifiability" of mythology studies. Frame all results as "consistent with" or "inconsistent with," never as proof of contact.

---

## Deliverables
1. `outputs/oral_tradition_mapping/RESULTS.md` — narrative with honest assessment
2. All data files and figures
3. If positive: flag as Substack piece ("The Myths Along the Circle")
4. If null: equally valuable — document that oral traditions do NOT preferentially follow the corridor
