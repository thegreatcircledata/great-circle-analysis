# Pillar 43 Stellar Encoding — Verdict

## Claim
Sweatman & Tsikritsis (2017) claim that animal carvings on Pillar 43 at Göbekli Tepe
encode stellar positions corresponding to the summer solstice sky at ~10,950 BCE.

Key identifications:
- Vulture → Sagittarius
- Scorpion → Scorpius
- Bird above scorpion → Libra
- Circle/disk → Sun at solstice

## Test Results

### A. Temporal Sweep (is 10,950 BCE the best epoch?)
| Metric | Value |
|---|---|
| Best matching epoch | -5000 CE |
| Best distance | 0.981332 |
| Distance at 10,950 BCE | 0.981500 |
| Total variation | 0.000280 (0.0% of mean) |

**NOT SUPPORTED — best epoch is -5000, not 10,950 BCE; variation is minimal**

### B. Random Arrangement Test
| Test | p-value | Interpretation |
|---|---|---|
| Fixed epoch (10,950 BCE) | 0.989800 | Not significant |
| Best epoch (temporal optimization) | 0.988000 | Not significant |

**NOT SUPPORTED — random arrangements match just as well (p = 0.9898)**

Note: The "best epoch" test is the fairer comparison — it asks whether the actual
carving arrangement, optimized across all epochs, matches better than random arrangements
similarly optimized.

### C. Permutation Test (is the claimed animal→constellation mapping optimal?)
- Tested all 24 permutations of 4 items
- Claimed mapping rank: **#22**
- Best mapping: {'vulture': 'Sagittarius', 'scorpion': 'Sun_solstice', 'bird_above_scorp': 'Libra', 'circle_disk': 'Scorpius'}

**NOT SUPPORTED — claimed mapping ranks #22/24, not optimal**

### D. Expanded Constellation Test
- Tested 6840 possible constellation assignments (20 constellations, all combos)
- Claimed set (Sgr/Sco/Lib) rank: **#689**

**NOT SUPPORTED — claimed set ranks #689/6840**

## Verdict: **NOT SUPPORTED**

No conditions met. The carving arrangement is consistent with chance.

### Methodological Caveats
1. **Carving positions are approximate** — extracted from published figures, not measured in situ.
   Precise positions could change results, though the random arrangement test partially controls for this.
2. **Constellation identification is subjective** — the animal→constellation mapping is an assumption
   of the hypothesis, not independently derived.
3. **The "Sun" position** — using the solstice Sun as one of the matched points builds in a
   constraint that may artificially improve the match for epochs near the solstice alignment.
4. **Multiple testing** — testing many epochs and finding the best match inflates significance
   if not corrected for. The "best epoch" random test addresses this.
