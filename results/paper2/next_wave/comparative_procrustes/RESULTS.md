# Study 12: Comparative Procrustes Across Claimed Stellar-Architectural Sites

## Overview
Three sites claimed to mirror Orion's Belt were tested using Procrustes analysis
against all C(167,3) = 762,355 bright-star triplets (mag < 4.0).

## Method
1. Site coordinates converted to local Cartesian (equirectangular at site latitude)
2. Star coordinates (J2000 RA/Dec) converted to Cartesian similarly
3. Procrustes distance D computed via `scipy.spatial.procrustes`
4. Orion Belt triplet ranked against all possible bright-star triplets
5. Collinearity control: restricted to triplets with angle > 160 deg
6. Benjamini-Hochberg correction across all tests

## Results

| Site | Claimed Stars | Procrustes D | Rank (all) | p (all) | Rank (collinear) | p (collinear) |
|------|--------------|-------------|-----------|---------|-----------------|---------------|
| giza | Alnitak, Alnilam, Mintaka | 0.000347 | 1,878/762,355 | 0.002463 | 1,878/123,260 | 0.015236 |
| teotihuacan | Alnitak, Alnilam, Mintaka | 0.080989 | 666,703/762,355 | 0.874531 | 116,575/123,260 | 0.945765 |
| thornborough | Alnitak, Alnilam, Mintaka | 0.005988 | 28,263/762,355 | 0.037073 | 10,840/123,260 | 0.087944 |

## Benjamini-Hochberg Correction

| Test | p (raw) | p (BH-adjusted) | Survives alpha=0.05? |
|------|---------|-----------------|---------------------|
| giza_all | 0.002463 | 0.014778 | Yes |
| giza_collinear | 0.015236 | 0.045708 | Yes |
| teotihuacan_all | 0.874531 | 0.945765 | No |
| teotihuacan_collinear | 0.945765 | 0.945765 | No |
| thornborough_all | 0.037073 | 0.074146 | No |
| thornborough_collinear | 0.087944 | 0.131916 | No |

## Interpretation

After BH correction, 2 test(s) survive at alpha=0.05: giza_all, giza_collinear. These site-star matches show alignment that is better than expected by chance among bright-star triplets.

### Collinearity Effect
Orion's Belt is near-collinear, and all three sites have near-collinear layouts. When restricting comparison to similarly collinear star triplets, the p-values generally increase (or remain high), indicating that the claimed match is driven primarily by shared collinearity rather than specific geometric correspondence.

## Metadata
- Stars in catalogue: 167
- Total triplets: 762,355
- Collinearity threshold: 160 deg
- Random seed: 42
