# O1 Builder's Precision Envelope — Results

## Parameters
- Naked-eye star accuracy: sigma = 0.5 deg (Schaefer 1993)
- Ground survey error: 58 mm / 230 m (Petrie) — negligible vs star error
- Monte Carlo trials: 10,000
- Seed: 42

## Intentional-Encoding Envelope
- Median D: 0.066310
- Mean D: 0.091314 (std: 0.081756)
- 90% CI: [0.003142, 0.248746]
- Range: [0.000001, 0.437339]

## Site Comparison

| Site | Observed D | Intentional Median | Within 90% CI? | Percentile | Interpretation |
|------|-----------|-------------------|----------------|------------|----------------|
| Giza | 0.000347 | 0.066310 | No | 0.80% | BETTER than expected — exceeds builder precision limit |
| Teotihuacan | 0.080989 | 0.066310 | Yes | 56.26% | Within intentional-encoding envelope |
| Thornborough | 0.005988 | 0.066310 | Yes | 8.05% | Within intentional-encoding envelope |

## Interpretation

The intentional-encoding envelope (90% CI) spans D = [0.003142, 0.248746],
with median D = 0.066310. This represents the range of Procrustes
distances an Old Kingdom builder would achieve if deliberately encoding the
Orion Belt pattern, limited only by naked-eye observational accuracy.

**Giza** (D = 0.000347, percentile = 0.80%) falls below the 5th percentile of the envelope. This means the Giza layout is *more precise* than what naked-eye encoding would typically produce — the match is suspiciously good, suggesting either extraordinary skill, a different surveying method, or coincidence in the Procrustes metric.

**Teotihuacan** (D = 0.080989, percentile = 56.26%) falls within the envelope.

**Thornborough** (D = 0.005988, percentile = 8.05%) falls within the intentional-encoding envelope, consistent with deliberate encoding.
