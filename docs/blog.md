# Do Ancient Sites Really Line Up? Testing the Great Circle Hypothesis with 62,000 Archaeological Sites

## The Claim

There's a fascinating idea that's been floating around the internet for years: many of the world's most famous ancient sites — the Great Pyramids, Nazca Lines, Angkor Wat, Easter Island, Persepolis, Machu Picchu — all lie suspiciously close to a single great circle around the Earth.

A great circle is the largest possible circle you can draw on a sphere — think of the equator, or any line of longitude. But this particular circle is tilted at about 30° from the equator, threading through Egypt, Peru, Cambodia, Iran, and the Pacific.

The idea was proposed by Jim Alison, and it's been shared widely in alternative history communities. But nobody had actually *tested* it statistically. Does the alignment hold up, or does it fall apart when you look at more than just a handful of cherry-picked famous sites?

We tested it. With 62,000 sites. Here's what happened.

## The Test

We obtained 61,870 archaeological sites from the Megalithic Portal — the world's largest database of megalithic and ancient sites — covering stone circles, burial chambers, standing stones, temples, pyramids, and dozens of other types. We added 114 supplementary sites to fill gaps in South America and Asia.

For each site, we calculated its distance from Alison's great circle. Then we asked: **is the number of sites near the circle more than you'd expect by chance?**

The tricky part is defining "chance." Ancient sites aren't randomly scattered — they cluster in Europe, along coastlines, near rivers. So we built a baseline that preserves the real geographic distribution of sites but breaks any alignment with the specific circle. We did this 200 times to build a statistical expectation.

## The Results

| Distance | Sites on circle | Expected by chance | Times more than expected | Z-score |
|---|---|---|---|---|
| 25 km | 207 | 45 | 4.6× | 24.3 |
| 50 km | 319 | 92 | 3.5× | 28.9 |
| 100 km | 454 | 182 | 2.5× | 21.2 |
| 200 km | 761 | 365 | 2.1× | 23.3 |

In statistics, a Z-score above 3 is considered highly significant. Above 5 is extraordinary. We're seeing **Z-scores above 20**. At 50 km, there are 3.5× more sites near the circle than our geographic-distribution-matched baseline would predict.

## But Wait — Is It Just Geographic Bias?

This is the most important question. The Megalithic Portal is heavily weighted toward Western Europe — Britain, Ireland, France — and the great circle passes right through those regions. Couldn't the apparent alignment just be because the circle happens to pass through the regions with the most data?

Our baseline is designed to handle exactly this. It preserves the latitude and longitude distributions of real sites, then shuffles them. The signal persists *after* this correction. But we can dig deeper.

### The Time Test

If this is a genuine ancient pattern, older sites should show stronger alignment than newer ones. We split the data:

- **Prehistoric sites** (stone circles, dolmens, passage graves, etc.): **Z ≈ 21** at 50 km
- **Later sites** (hillforts, holy wells, museums, etc.): **Z ≈ 8** at 50 km

The prehistoric signal is 2.5× stronger. This is consistent with an ancient pattern that gets diluted as later people build sites more randomly.

### The Competition Test

Is Alison's circle special, or would *any* great circle show similar alignment? We tested 1,000 random great circles. Alison's circle ranked #1.

### The Falsification

We also tested a related claim — that ancient sites preferentially form pairs separated by exactly 108° of arc. This one is dead on arrival: Z ≈ -1.4. Fewer pairs at 108° than random chance. **Falsified.**

## What Does It Mean?

A statistically significant result doesn't mean ancient people were intentionally building along a great circle. Some possible explanations:

1. **Migration corridors** — The circle may trace ancient paths of human movement
2. **Geological features** — Fault lines, mountain ranges, or coastlines that attracted settlement
3. **Residual bias** — Our null model may not fully capture all the ways site locations are spatially correlated

What we *can* say is that the clustering is real, it's large, and it's not explained by the obvious confound of geographic data density. The question of *why* remains open.

## Try It Yourself

All code is Python 3.7+ with no external dependencies. To reproduce:

```bash
git clone https://github.com/[repo]
cd deep-time-atlas
python analysis/great_circle_test.py --input data/merged_sites.csv --output my_results.json --trials 200
```

Or run the full analysis:
```bash
python analysis/run_all_tests.py --trials 200
```

## Acknowledgments

- **Andy Burnham** and the Megalithic Portal community for building the world's best megalithic site database
- **Jim Alison** for the original great circle observation
- **Pleiades** project for independent ancient place data
