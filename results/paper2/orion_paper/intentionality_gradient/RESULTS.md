# O6 — Intentionality Gradient

Multi-axis intentionality scoring for Orion Belt correlation claims.

## Axes

| # | Axis | Metric | Range |
|---|------|--------|-------|
| 1 | Shape match | -log10(p_procrustes) | 0+ (higher = better fit) |
| 2 | Brightness-size | Spearman rho of size rank vs star brightness rank | 0-1 |
| 3 | Orientation | 1 - (offset/90 deg) | 0-1 |
| 4 | Textual evidence | Categorical | 0-3 |
| 5 | Pattern extension | Categorical | 0-3 |

## Scoring Table

| Axis | Giza | Teotihuacan | Thornborough | Gobekli Tepe P43 |
|------|------|-------------|--------------|------------------|
| Shape match (-log10 p) | 2.61 | 0.06 | 1.43 | 0.04 |
| Brightness-size (rho) | 0.50 | 0.50 | 0.00 | N/A |
| Orientation (0-1) | 0.44 | 0.33 | 0.33 | N/A |
| Textual evidence (0-3) | 3 | 1 | 0 | 1 |
| Extension (0-3) | TBD | 0 | 0 | 0 |
| **Composite (0-1)** | **0.704** | **0.237** | **0.162** | **0.115** |

## Ranking

1. **Giza** — composite = 0.704 (4 axes)
2. **Teotihuacan** — composite = 0.237 (5 axes)
3. **Thornborough** — composite = 0.162 (5 axes)
4. **Gobekli Tepe P43** — composite = 0.115 (3 axes)

## Notes

- Giza Extension axis = TBD (pending O2 queens' extension results)
- Gobekli Tepe axes 2-3 = N/A (single pillar, not three structures)
- Composite = mean of all applicable normalized axes per site
- Brightness-size for Giza: brightest star (Alnilam, mag 1.69) pairs with second-largest pyramid (Khafre), not the largest — partial match only (rho = 0.5)
- Thornborough henges are near-identical in diameter, producing no meaningful size gradient to correlate with Belt star brightness
