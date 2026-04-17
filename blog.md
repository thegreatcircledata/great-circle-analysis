# Do Ancient Sites Really Line Up? Testing the Great Circle Hypothesis

> **This post has been revised (2026-04-17).** The original version reported pre-land-constraint Monte Carlo Z-scores that substantially overstated the effect. The corrected, land-constrained values below reproduce the peer-reviewed manuscript currently under review at PLOS ONE.

## The Claim

There's a claim that has circulated in alternative-history writing for two decades: many of the world's most famous ancient sites — the Great Pyramids at Giza, Nazca, Angkor Wat, Easter Island, Persepolis, Mohenjo-daro — lie suspiciously close to a single great circle around the Earth. A great circle is the intersection of Earth's surface with any plane through its centre; this one is tilted roughly 30° from the equator and passes through North Africa, South Asia, the Pacific, and South America.

The claim was proposed by Jim Alison in 2001. Informal demonstrations typically cite a handful of famous sites. Until recently, no formal statistical test had been published. We tested it with eight independent archaeological databases comprising approximately 259,000 unique sites.

## The Test

For each site in each database we calculated its distance from Alison's great circle (pole: 59.682°N, 138.646°W). We then asked: is the number of sites near the circle more than you'd expect by chance?

Defining "chance" is the hard part. Archaeological sites are not uniformly distributed — they cluster in Europe, along coastlines, near rivers. Standard Monte Carlo baselines jitter site coordinates to build an expected distribution, but when corridors intersect coastlines or islands, jittered points frequently land in the ocean, depressing expected on-corridor counts and inflating Z-scores. We introduce a **land-constrained Monte Carlo baseline** that rejects invalid ocean landings. This correction reduces reported Z-scores by 34–73% across the four databases tested, demonstrating that uncorrected coastal jitter bias substantially overestimates enrichment.

## The Results (land-constrained)

| Database | Sites | Z-score (50 km, LC) | Enrichment |
|---|---|---|---|
| Megalithic Portal | 61,913 | 8.77 (N = 10,000) | 1.62× |
| Pleiades ancient monuments | 406 (pre-2000 BCE) | 6.74 | 2.52× |
| XRONOS radiocarbon | 28,127 | 5.91 (N = 1,000) | 1.40× |
| p3k14c radiocarbon | 36,693 | 2.09 (p = 0.02; marginal) | 1.16× |
| Historic England (negative control) | 20,026 | 0.0 | null as expected |

After habitability adjustment (comparing the circle against equally populated corridors), overall enrichment sits at the 78.5th percentile — not significant. **The aggregate alignment is geographic, not anomalous.**

## What Does Survive

Within those corridors, monuments cluster while contemporaneous settlements do not. On Pleiades pre-2000 BCE data, the monument–settlement divergence D = Z_monument − Z_settlement = 9.65 (land-constrained). No random great circle out of 10,000 produces a comparable divergence (max random D = 0.90), and the observed value sits at the 99.63rd percentile among 4,323 habitability-matched circles.

The divergence is concentrated in a single 250-year construction window in the Memphis necropolis during the Egyptian Old Kingdom (2750–2500 BCE, 31/38 pyramids within 50 km). Outside that window, the signal weakens sharply.

## What Is Falsified

- **108° pair hypothesis.** Z = −8.01. Fewer pairs at 108° separation than random chance. Falsified.
- **South American anomaly.** Null on three independent databases (Peru Ministry, SAAID, Pleiades/p3k14c).
- **Anatolian anomaly.** LuwianSiteAtlas: 0 of 483 sites within 500 km of the circle.
- **Astronomical alignments.** Null after Bonferroni correction.
- **Monument orientations.** Rayleigh p = 0.75.
- **Trade network correlation.** r = 0.04.
- **Pre-Younger Dryas corridor.** Z = −0.47 (33rd percentile).

## What It Means

A statistically significant monument–settlement divergence in one 250-year window of Egyptian pharaonic construction does not imply that ancient people intentionally built along a great circle. Some explanations to consider:

1. **Tectonic corridors.** The circle crosses an unusually high number of plate boundaries (28 vs 13.9 expected, p = 0.007) and arid–fertile transitions (29 vs 4.5 expected).
2. **Geographical circumscription.** Old World monument-building civilizations cluster in desert-bounded alluvial corridors.
3. **Residual bias.** Our null models may not fully capture every spatial-correlation structure in ancient site placement.

What we can say is that the Egyptian Old Kingdom monument–settlement divergence is real, is large, survives every control we could construct, and is not explained by the obvious confounds of data density, habitability, or coastal jitter. The *cause* remains open.

## Reproduce

```bash
git clone https://github.com/thegreatcircledata/great-circle-analysis.git
cd great-circle-analysis
pip install -r requirements.txt
bash download_data.sh
python analysis/great_circle_test.py --input data/merged_sites.csv --output my_results.json
```

All verified headline numbers are in [`outputs/merged_paper_prep/SUMMARY.md`](outputs/merged_paper_prep/SUMMARY.md). The current manuscript is under review at PLOS ONE; preprint: [10.5281/zenodo.19212669](https://doi.org/10.5281/zenodo.19212669).

## Acknowledgments

- **Andy Burnham** and the Megalithic Portal community for building the world's largest megalithic site database.
- **Jim Alison** for the original great-circle observation.
- **Pleiades** for the independent ancient-places gazetteer.
- **p3k14c, XRONOS, DARE** for independent global radiocarbon and typology data.

---

Elliot Allan, Deep Time Research Institute — elliot@deeptime-research.org
