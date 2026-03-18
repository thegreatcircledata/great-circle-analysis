# Data Sources

Raw site data from restricted sources is NOT included in this repository.
To reproduce the analysis, obtain data directly from the original sources listed below.
Analysis scripts in /analysis/ expect these files in the /data/ directory.

## Primary Databases

**Megalithic Portal** — 61,913 sites deduplicated from 62 KML files.
Source: megalithic.co.uk (membership required). Data is NOT redistributable.
To reproduce: join the Megalithic Portal, download KML files, and run the parsing scripts.

**Pleiades Gazetteer** (pleiades-places-latest.csv) — 34,470 ancient places.
Source: pleiades.stoa.org (CC BY 3.0). Download: https://pleiades.stoa.org/downloads

**p3k14c** (p3k14c_data.csv) — 173,946 radiocarbon dates, 36,693 unique sites.
Source: Bird et al. 2022, Scientific Data. Archived in tDAR: https://core.tdar.org/collection/70213/p3k14c-data

## Validation

**DARE** — 29,760 Roman Empire places.
Source: dare.ht.lu.se (CC BY-SA). API: https://imperium.ahlfeldt.se/api/

## Negative Control

**Historic England** — 20,026 scheduled monuments.
Source: historicengland.org.uk (Open Government Licence). Download: https://historicengland.org.uk/listing/the-list/data-downloads/

## Reference

**circle_coordinates.csv/json** — Great Circle coordinates (included in repo).
Pole: 59.682122°N, 138.646087°W
