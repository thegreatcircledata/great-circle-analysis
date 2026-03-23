# O2: Queens' Pyramids Extension Study

## Phase 1: Nine-Pyramid Layout

- G1 queens (G1a/G1b/G1c) collinearity: **180.00 deg** (180 = perfect line)
- G3 queens (G3a/G3b/G3c) collinearity: **180.00 deg**
- G1 queen separations: 100 m, 100 m
- G3 queen separations: 164 m, 164 m

## Phase 2: 6-Point Belt + Sword Match

Tests Bauval's claim that G1 queens correspond to Orion's Sword (theta1 Ori, iota Ori, 42 Ori).

- Observed Procrustes d (fixed correspondence): **0.125292**
- MC null (100,000 random 6-star sextuplets): rank **39** / 100001
- p-value: **0.000390**
- MC median: 0.766007, 5th percentile: 0.407007

**Result:** Significant at p < 0.05.

## Phase 3: 7-Point Saqqara-Dahshur Extension

Tests Orofino's broader claim mapping 7 pyramids to 7 Orion stars.

| Pyramid | Star |
|---------|------|
| Khufu | Alnitak |
| Khafre | Alnilam |
| Menkaure | Mintaka |
| Djedefre | Bellatrix |
| Djoser | Saiph |
| Bent | Rigel |
| Red | Betelgeuse |

- Observed Procrustes d (fixed correspondence): **0.872056**
- MC null (100,000 random 7-star septuplets): rank **49793** / 100001
- p-value: **0.497925**
- MC median: 0.872853, 5th percentile: 0.452199

**Result:** NOT significant at p < 0.05. The 7-point extension is unremarkable.

## Phase 4: Queens-Only Exhaustive Triplet Test

Tested G1a/G1b/G1c against all 695,520 star triplets (all 6 permutations).
Queens collinearity: 180.00 deg.

- **Best match:** Deneb, Sheratan, Alula Bor (d = 0.000000)
- p-value: **0.000001** (rank 1 / 695,520)

### Top 10 Matches

| Rank | Procrustes d | Collinearity | Stars |
|------|-------------|-------------|-------|
| 1 | 0.000000 | 180.0 deg | Deneb, Sheratan, Alula Bor |
| 2 | 0.000001 | 179.9 deg | Muphrid, Algieba, Tianguan |
| 3 | 0.000001 | 179.8 deg | Alhena, Schedar, Kraz |
| 4 | 0.000001 | 179.8 deg | Sirius, Zubenelgenubi, Skat |
| 5 | 0.000002 | 179.8 deg | Schedar, Zosma, Deneb Algedi |
| 6 | 0.000002 | 179.7 deg | Menkalinan, Schedar, Alula Bor |
| 7 | 0.000003 | 179.7 deg | Alpheratz, Phecda, Alrai |
| 8 | 0.000004 | 179.8 deg | Naos, Menkar, Arneb |
| 9 | 0.000004 | 179.7 deg | Spica, Nashira, Cursa |
| 10 | 0.000004 | 180.0 deg | Enif, Kornephoros, Alula Bor |

### Best match containing Iota Ori (sword star)
- Stars: Gienah, Iota Ori, Ascella
- d = 0.000025, rank = 57 / 695,520
- p-value: 0.000082

## Summary

| Test | Procrustes d | p-value | Significant? |
|------|-------------|---------|--------------|
| 6-point Belt+Sword | 0.125292 | 0.000390 | Yes |
| 7-point Orofino extension | 0.872056 | 0.497925 | No |
| Queens-only (best of all) | 0.000000 | 0.000001 | N/A (best possible) |

Note: All MC nulls use fixed correspondence (testing the specific claimed mapping,
not the best of all permutations). This is the appropriate test for an a priori claim.
