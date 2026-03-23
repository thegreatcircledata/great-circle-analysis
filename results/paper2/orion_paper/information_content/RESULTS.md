# O5 Information-Theoretic Framing

## Phase 1: Degrees of Freedom Analysis

Three points in 2D = 6 coordinates. Procrustes alignment removes 4 DOF
(2 translation + 1 rotation + 1 scale), leaving **2 independent shape parameters**:
1. **Distance ratio** d(Khufu,Khafre) / d(Khafre,Menkaure)
2. **Bend angle** at the middle point (Khafre / Alnilam)

| Parameter | Giza | Orion Belt | Difference |
|-----------|------|-----------|------------|
| Distance ratio | 1.0130 | 0.9785 | 0.0344 |
| Bend angle | 168.80 deg | 172.48 deg | 3.68 deg |

## Phase 2: Per-Parameter Significance

Evaluated all 690,015 valid star triplets (from 695,520 total, 162 stars mag < 4.0).

| Parameter | Tolerance | Matches | p-value | Bits |
|-----------|-----------|---------|---------|------|
| Distance ratio | +/-0.0344 | 21,336 | 0.0309 | 5.02 |
| Bend angle | +/-3.68 deg | 19,636 | 0.0285 | 5.14 |
| **Joint (both)** | -- | 532 | 0.000771 | 10.34 |

Independence check: observed joint p = 0.000771, expected if independent = 0.000880, ratio = 0.88x.

## Phase 3: Bits of Information

The Procrustes shape match (p = 0.0025) encodes:
- **8.67 bits** of information
- Theoretical maximum with 2 DOF at ~1m precision over ~486m span: **17.9 bits**
- Observed = 49% of theoretical max

## Phase 4: Comparison to Known Encodings

| Claimed Encoding | p-value | Bits of Information |
|-----------------|---------|---------------------|
| Pi ratio (2b/h = 3.1454 vs pi) | 0.002171 | 8.85 |
| 43,200 scale factor | 0.0300 | 5.06 |
| **Orion shape match** | 0.0025 | **8.67** |

Ranking by information content: Pi ratio (8.8) > Orion shape (8.7) > 43,200 scale (5.1).

The Orion shape match and pi ratio carry comparable information content
(~8.7-8.9 bits), both substantially exceeding the 43,200 scale factor (~5.1 bits).
