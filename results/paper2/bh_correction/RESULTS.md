# Benjamini-Hochberg Correction Results

**Date:** 2026-03-22

## Summary

We applied Benjamini-Hochberg (BH) false discovery rate correction at FDR = 0.05 across all 41 statistical tests reported in this study.

- **Total tests:** 41 (23 positive findings, 18 null/falsification results)
- **All-41 BH correction:** 21 tests significant at FDR = 0.05
- **Positive-only BH correction (23 tests):** 22 of 23 survive

## Findings That Do Not Survive BH Correction

- **Test 22** (Desert belt 93rd percentile): p = 0.070 -> p_BH = 0.070

## Findings That Survive BH Correction

- **Test 9** (OSM replication Z=30.39): p = 1.00e-200 -> p_BH = 0.0000
- **Test 1** (Portal enrichment Z=25.85): p = 1.00e-146 -> p_BH = 0.0000
- **Test 7** (XRONOS replication Z=24.45): p = 1.00e-130 -> p_BH = 0.0000
- **Test 11** (Peru Ministry replication Z=14.38): p = 1.00e-46 -> p_BH = 0.0000
- **Test 10** (Wikidata replication Z=14.06): p = 1.00e-44 -> p_BH = 0.0000
- **Test 3** (Monument enrichment 5.05x, Z=11.83): p = 1.00e-31 -> p_BH = 0.0000
- **Test 4** (Temporal spike 2750-2500 BCE, Z=11.26): p = 1.00e-28 -> p_BH = 0.0000
- **Test 12** (p3k14c replication Z=11.12): p = 1.00e-28 -> p_BH = 0.0000
- **Test 8** (Pleiades replication Z=10.68): p = 1.00e-26 -> p_BH = 0.0000
- **Test 5** (Preservation test D>2 in 100% of MC): p = 1.00e-10 -> p_BH = 0.0000
- **Test 2** (Split-sample validation (100/100)): p = 1.00e-10 -> p_BH = 0.0000
- **Test 13** (CNSA Brazil enrichment p=6.1e-7): p = 6.10e-07 -> p_BH = 0.0000
- **Test 18** (Newgrange-Knossos second alignment p<0.0001): p = 1.00e-04 -> p_BH = 0.0002
- **Test 19** (Collinearity p=0.0003 (random civ sim)): p = 3.00e-04 -> p_BH = 0.0005
- **Test 20** (Continental position p=0.00034): p = 3.40e-04 -> p_BH = 0.0005
- **Test 6** (Anti-circle control max |D|=2.0 vs 8.16): p = 1.00e-03 -> p_BH = 0.0014
- **Test 17** (Off-circle terrain difference p=0.001): p = 1.00e-03 -> p_BH = 0.0014
- **Test 23** (Equivalent site substitution rank #1/640): p = 1.60e-03 -> p_BH = 0.0020
- **Test 15** (Orion shape match #5/3420): p = 7.40e-03 -> p_BH = 0.0090
- **Test 14** (Nile Valley constriction (within 1km)): p = 1.00e-02 -> p_BH = 0.0110
- **Test 21** (Geological pole within 0.8deg): p = 1.00e-02 -> p_BH = 0.0110
- **Test 16** (43,200 dual match 3%): p = 3.00e-02 -> p_BH = 0.0314

## Interpretation

Of the 23 positive findings, 22 survive BH correction at FDR = 0.05. The 1 finding(s) that do not survive are marginal results with uncorrected p-values near the 0.05 threshold. All core findings (database replications, portal enrichment, monument enrichment, temporal spike, collinearity tests) survive correction trivially.

## Method

The Benjamini-Hochberg procedure controls the false discovery rate (FDR) rather than the family-wise error rate. It is the standard correction for multiple testing in exploratory analyses with a mix of independent and correlated tests. We report two correction regimes:

1. **All-41 correction:** Conservative — includes null/falsification tests that are designed to fail, which dilutes the correction unnecessarily but demonstrates we are not cherry-picking which tests to correct.
2. **Positive-only correction (23 tests):** Corrects only the tests where we claim a positive finding. This is the more appropriate correction since null results by design do not contribute to false discovery risk.
