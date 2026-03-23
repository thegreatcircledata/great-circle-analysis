# C5: Ceremonial Geography Follow-Up

## Context
C3 showed all 41 on-corridor pre-2000 BCE mortuary sites are in Egypt. C4 showed temples/sanctuaries are enriched (Z=+3.76) and combined pre-2000 BCE sacred sites yield Z=+32.3. This test determines whether the ceremonial enrichment extends beyond Egypt.

## Phase 1: On-Corridor Pre-2000 BCE Ceremonial Sites
- Total pre-2000 BCE ceremonial (excl church/monastery/synagogue): 30
- On-corridor (<50km): 2
- Regions represented: Egypt/Nubia

On-corridor sites:

| Name | Lat | Lon | Dist (km) | Region | minDate | featureTypes |
|------|-----|-----|-----------|--------|---------|--------------|
| Temple of Amenhotep II | 29.976 | 31.138 | 6.1 | Egypt/Nubia | -2670 | temple, temple-2 |
| Giza | 29.978 | 31.132 | 6.4 | Egypt/Nubia | -2670 | sanctuary, cemetery, pyramid,  |

## Phase 2: Per-Region Ceremonial Enrichment
MC method: 2° Gaussian jitter, 1000 trials.

| Region | N total | N on-corridor | MC mean | MC std | Z |
|--------|---------|---------------|---------|--------|---|
| Egypt/Nubia | 8 | 2 | 0.6 | 0.7 | +2.002 |
| Anatolia | 2 | 0 | 0.0 | 0.0 | +0.000 |
| North_Africa_nonEgypt | 6 | 0 | 0.0 | 0.0 | +0.000 |
| Other | 11 | 0 | 0.0 | 0.03 | -0.032 |
| Mesopotamia | 2 | 0 | 0.0 | 0.04 | -0.045 |
| Iran/Central_Asia | 1 | 0 | 0.02 | 0.13 | -0.132 |

## Phase 3: Framing Decision
- Regions with on-corridor sites: 1 (Egypt/Nubia)
- Non-Egypt regions with Z > 1: 0 (none)
- Egypt fraction of on-corridor: 100.0%
- **Recommended framing: "Old Kingdom sacred architecture alignment"**

## Phase 4: Egypt vs Non-Egypt Decomposition
- Full sacred pre-2000 BCE: n=177, on=42, Z=+32.987
- Egypt-only: n=71, on=42, Z=+52.786
- Non-Egypt: n=106, on=0, Z=-0.968

Per-region divergence (sacred Z - settlement Z):

| Region | Sacred N | Sacred On | Sacred Z | Settle N | Settle On | Settle Z | Divergence |
|--------|----------|-----------|----------|----------|-----------|----------|------------|
| Egypt/Nubia | 71 | 42 | -0.008 | 43 | 9 | -0.035 | +0.026 |
| Anatolia | 7 | 0 | +0.000 | 61 | 0 | +0.000 | +0.000 |
| East_Africa | 1 | 0 | +0.000 | 3 | 0 | +0.000 | +0.000 |
| Levant | 6 | 0 | +0.000 | 72 | 0 | +0.000 | +0.000 |
| North_Africa_nonEgypt | 6 | 0 | +0.000 | 4 | 0 | +0.000 | +0.000 |
| Other | 75 | 0 | +0.000 | 97 | 0 | +0.000 | +0.000 |
| South_Asia | 1 | 0 | +0.000 | 4 | 0 | +0.000 | +0.000 |
| Iran/Central_Asia | 5 | 0 | +0.000 | 40 | 6 | +0.025 | -0.025 |
| Mesopotamia | 5 | 0 | +0.000 | 112 | 1 | +0.032 | -0.032 |

## Phase 5: Sacred Site Catalog
- Total on-corridor pre-2000 BCE sacred sites: 42
- By category: {'mortuary': 40, 'ceremonial': 1, 'mortuary+ceremonial': 1}
- By region: {'Egypt/Nubia': 42}
- Saved to: sacred_site_catalog.csv

## Conclusion
The pre-2000 BCE ceremonial enrichment is overwhelmingly Egyptian (100% of on-corridor sites). Combined with C3's finding that all 41 mortuary sites are Egyptian, the paper should frame this as an **Old Kingdom sacred architecture alignment** — the corridor signal in sacred architecture is specific to Old Kingdom Egypt, not a pan-regional phenomenon.

