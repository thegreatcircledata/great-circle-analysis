# Hardening Round Verdict

**Date:** 2026-03-22
**Tests completed:** 9 / 9 (A1-A4, B1, C1-C4)

---

## Block A: Orion Combined Significance

- **A1 (Gaia upgrade): CONFIRMED** — Gaia DR3 proper motions reproduce BSC results exactly. Best-fit epoch identical (-11,625), alignment windows within 2%, shape Procrustes stable (max delta 0.000616). Proper motion precision is not a vulnerability.

- **A2 (Window sensitivity): ROBUST** — Combined p-value is remarkably stable: p=0.000282 at 1°, p=0.000743 at 3°, p=0.00119 at 5°. Combined p does not exceed 0.01 until 45° and never exceeds 0.05 in the tested range. At archaeoastronomically justified thresholds (1-3°, per Schaefer/Ruggles/Aveni), p remains well below 0.001. Not threshold-dependent.

- **A3 (Construction epoch): NOT SIGNIFICANT** — At 2560 BCE, the Belt-to-pyramid offset is 66.85° (p=0.74 vs random axis). Among top shape-matched triplets, Orion ranks 5,564th in orientation (p=0.85). The orientation at the construction epoch adds no meaningful evidence beyond the shape match. The combined significance from Study 11 (p=0.003) is driven entirely by the best-fit epoch (11,750 BCE), which has no archaeological support.

- **A4 (Independence): FISHER VALID, BUT MOOT** — Shape and orientation are effectively independent (Spearman rho=-0.041), so Fisher's method is technically valid. However, Orion's orientation p at construction is 0.98 (nearly worst possible), making the combined Fisher p=0.018 — weaker than shape alone. Permutation combined p=0.489 confirms orientation contributes nothing.

- **Revised combined p-value:** The meaningful number is the **shape p-value alone: p=0.0025** (from Study 12, survives BH correction). The orientation/epoch analysis is exploratory context only. The combined p=0.003 from Study 11 should not be used as the headline statistic — it depends on a best-fit epoch with no archaeological basis.

- **Paper-ready?** **Yes, with reframing.** Lead with shape match (p=0.0025, survives BH, confirmed with Gaia DR3). Present orientation/epoch analysis as exploratory: "the shape match is genuine; the orientation sweeps through alignment at ~11,750 BCE due to precession, but this has no known archaeological significance, and the orientation at the construction epoch (2560 BCE) is unremarkable." Drop the combined p=0.003 claim.

---

## Block B: Anti-Divergence Orthogonality

- **B1.1 (Fine grid): POLE SHIFTED** — The original (10°N, -30°W) pole is not the optimum on a finer grid. Metric A yields (30°N, -48°E), metric B yields (52°N, -152°E), metric C yields (0°N, -104°E). The three divergence metrics do not converge.

- **B1.2 (Multi-metric): METRIC-DEPENDENT** — Different definitions of "anti-divergence" yield different poles with different angular separations from Alison. The result is not robust to metric choice.

- **B1.3 (Threshold sensitivity): UNSTABLE** — Orthogonality appears only at exactly 50km. At 25km: 44.9°; at 75km: 72.7°; at 200km: 43.5°. The 90° finding is threshold-specific.

- **B1.4 (Bootstrap CI): WIDE** — 95% CI is [50.88°, 90.52°], indicating the anti-divergence pole location is unstable under resampling.

- **B1.5 (Null model): 90° IS EXPECTED** — 46.5% of random reference poles produce separations ≥90.52° to the anti-divergence pole. Mean null separation is 83.52° ± 41.65°. Orthogonality is a geometric property of the sphere (independently optimized circles tend toward orthogonality because orthogonal great circles sample maximally different regions), not a specific property of the Alison circle.

- **Include in paper?** **No.** The orthogonality finding does not survive hardening on any dimension: the pole is unstable, metric-dependent, threshold-specific, and the ~90° separation is expected for any reference pole. Remove from Paper v2.1.

---

## Block C: Mortuary Narrowing

- **C1 (DARE replication): CONFIRMED** — DARE independently confirms mortuary enrichment (Z=+5.978 vs Pleiades +7.584, p=0.0000). Admin depletion sign agrees (both negative). Ceremonial diverges between databases (DARE enriched, Pleiades depleted — classification boundary differences). The core mortuary finding replicates across independent taxonomies.

- **C2 (Leave-one-out): STABLE, min Z = 13.66** — Removing the most influential site (Helwan necropolis) only drops Z by 1.79. Leave-five-out: 100% of 1000 runs maintain Z>3, minimum observed Z=11.40. The signal is distributed across many sites, not driven by outliers.

- **C3 (Outside Egypt): EGYPT-ONLY** — All 41 on-corridor pre-2000 BCE mortuary sites are in Egypt/Nubia. Zero mortuary enrichment in any other segment (non-Egypt aggregate Z=0.000). 90.2% are pyramid or tomb. The temporal driver is 3000-2000 BCE (Old Kingdom). This is specifically an Egyptian pyramid/tomb corridor.

- **C4 (Ceremonial depletion): EXPLAINED** — The overall ceremonial "depletion" (Study 4's Z=-2.287) is entirely driven by post-Roman **churches** (Z=-1.781), which proliferate in Mediterranean regions the circle avoids. Temples/sanctuaries remain enriched overall (Z=+3.758) and in every time period. Pre-2000 BCE combined sacred sites (mortuary+ceremonial) yield Z=+32.327 vs settlement Z=+6.447. The depletion is a church geography artifact.

- **Recommended framing:** The hardening results support two viable framings depending on scope:

  1. **Narrow (conservative):** "Egyptian Old Kingdom pyramid/tomb corridor" — defensible given C3 (Egypt-only mortuary enrichment). This is the most honest framing: the signal is 41 Egyptian mortuary sites, concentrated 3000-2000 BCE.

  2. **Broader (supported by C4):** "Pre-2000 BCE sacred architecture corridor" — defensible because temples/sanctuaries are ALSO enriched (Z=+3.758), and combined sacred sites yield Z=+32.327. However, C3 shows the mortuary component (the strongest sub-signal) is Egypt-only. The ceremonial enrichment's geographic distribution needs further analysis to confirm it extends beyond Egypt.

  **Recommendation:** Use framing #1 for the core claim (it's bulletproof), and note framing #2 as a broader pattern that warrants further investigation. Do NOT claim a global "sacred corridor" without testing whether ceremonial enrichment also extends beyond Egypt.

---

## Summary Table

| Test | Verdict | Implication |
|------|---------|-------------|
| A1 Gaia upgrade | Confirmed | Astrometry is precise |
| A2 Window sweep | Robust | Not threshold-dependent |
| A3 Construction epoch | Not significant | Orientation adds nothing at 2560 BCE |
| A4 Independence | Fisher valid, moot | Shape alone is the finding |
| B1 Orthogonality | Geometric truism | Remove from paper |
| C1 DARE replication | Confirmed | Mortuary enrichment is real |
| C2 Leave-one-out | Stable (min Z=13.66) | Not driven by outliers |
| C3 Outside Egypt | Egypt-only | Signal is geographically specific |
| C4 Ceremonial depletion | Church artifact | Temples actually enriched |
| C5 Ceremonial geography | Egypt-only | "Sacred corridor" not supported |

## C5 Follow-Up: Ceremonial Geography (Added 2026-03-22)

- **C5 (Ceremonial regional distribution): EGYPT-ONLY** — Only 2 pre-2000 BCE ceremonial sites on-corridor (Temple of Amenhotep II, Giza sanctuary), both Egyptian. Zero non-Egypt regions show ceremonial enrichment (all Z ≈ 0). Egypt-only sacred Z=+52.79; non-Egypt sacred Z=-0.97 with 0 on-corridor sites. All 42 on-corridor pre-2000 BCE sacred sites are Egyptian.

- **Framing recommendation updated:** The C4 finding that "temples are enriched (Z=+3.76)" was driven by the overall category across all periods. When restricted to pre-2000 BCE and excluding churches, the ceremonial component is tiny (2 sites) and entirely Egyptian. The "pre-2000 BCE sacred corridor" framing is **not supported**. The correct framing is **"Old Kingdom sacred architecture alignment."**

---

## Revised Paper Recommendations (updating SYNTHESIS.md)

### Orion Paper
- **Lead with shape p=0.0025** (BH-corrected, Gaia-confirmed, unique among claimed sites)
- Present orientation/epoch sweep as exploratory context, NOT as combined significance
- Drop the combined p=0.003 claim
- Retain: Teotihuacan rejected (p=0.875), Thornborough marginal (p=0.074, doesn't survive BH)

### Great Circle Paper (v2.1)
- **Reframe as Egyptological finding:** "Old Kingdom monumental architecture (pyramids, tombs, temples) is significantly concentrated along this great circle (Z=+7.58 Pleiades, Z=+5.98 DARE), a signal that is robust to site removal (leave-one-out min Z=13.66) and not attributable to publication bias (Study 1) or classification artifacts (Study 2)."
- **C5 narrows the scope:** All 42 on-corridor pre-2000 BCE sacred sites are Egyptian. No non-Egypt region shows enrichment. This is not a global "sacred corridor" — it is a property of Old Kingdom Egyptian monumental architecture and the Nile Valley's geographic relationship to this great circle.
- **Remove:** 36° periodicity claim (falsified), anti-divergence orthogonality (geometric truism), crop diffusion corridor (counter-evidence), "sacred corridor" framing (C5 refutes)
- **Address:** Preservation bias vulnerability (Study 7, break-even at 0.47x) — argue that the DARE replication (independent database, same result) provides cross-validation that the signal is real, even if the magnitude is uncertain.
- **Target journal:** Journal of Egyptian Archaeology or JHAA (not JAS/Antiquity — the finding is too geographically specific for a broad archaeology journal)
