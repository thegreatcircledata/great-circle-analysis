# C1 — DARE Independent Type-Code Replication

## Purpose
Replicate the mortuary-enrichment finding (Pleiades Z = +7.584) using
DARE's independent archaeological classification, confirming the signal
is not an artefact of the Pleiades type taxonomy.

## Dataset
- **Source:** DARE places2.geojson (29,760 features)
- **Corridor:** < 50.0 km from great circle (pole 59.682122, -138.646087)
- **On-corridor sites:** 162 (0.54% base rate)
- **MC trials:** 1,000, seed 42

## Category Counts

| Category | DARE Type Codes | N sites |
|----------|----------------|---------|
| mortuary | 32, 63 | 726 |
| ceremonial | 21, 22, 24, 61 | 1,925 |
| settlement | 11, 12, 13, 14, 16, 31, 34, 35, 56, 57, 58, 68 | 20,427 |
| admin | 15, 17, 18, 51, 52, 53, 54, 59, 60, 65, 66, 67, 69, 70, 71, 72, 76 | 5,061 |
| other | (remainder) | 1,621 |

## Results — Comparison Table

| Type | Pleiades Z | DARE shuffle Z | DARE sample Z | Sign agrees? |
|------|-----------|---------------|---------------|-------------|
| mortuary | +7.584 | +5.978 | +6.175 | YES |
| ceremonial | -2.287 | +3.688 | +3.687 | NO |
| settlement | — | -2.386 | -2.420 | — |
| admin | -4.345 | -0.508 | -0.503 | YES |

## Detailed Per-Category Results

### Mortuary
- N sites: 726
- Observed on-corridor: 16
- **Shuffle baseline:** mean 4.07, std 1.996, Z = +5.978, p = 0.0000
- **Sample baseline:** mean 3.91, std 1.958, Z = +6.175, p = 0.0000

### Ceremonial
- N sites: 1,925
- Observed on-corridor: 22
- **Shuffle baseline:** mean 10.50, std 3.118, Z = +3.688, p = 0.0010
- **Sample baseline:** mean 10.48, std 3.124, Z = +3.687, p = 0.0010

### Settlement
- N sites: 20,427
- Observed on-corridor: 97
- **Shuffle baseline:** mean 111.04, std 5.885, Z = -2.386, p = 0.0110
- **Sample baseline:** mean 111.18, std 5.859, Z = -2.420, p = 0.0080

### Admin
- N sites: 5,061
- Observed on-corridor: 25
- **Shuffle baseline:** mean 27.44, std 4.796, Z = -0.508, p = 0.3410
- **Sample baseline:** mean 27.43, std 4.837, Z = -0.503, p = 0.3520

## Verdict

Mortuary shuffle Z = **+5.978** (Pleiades reference: +7.584).

**PASS** — Mortuary enrichment replicates under DARE's independent type taxonomy at high significance. The signal is not an artefact of the Pleiades classification system.

*Elapsed: 2.4s*
