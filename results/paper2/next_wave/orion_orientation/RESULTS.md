# Study 11 — Best-Fit Orientation Epoch (Orion Paper)

**Date:** 2026-03-22
**Precession method:** astropy FK5

---

## Key Finding

The Giza pyramid diagonal axis (Menkaure to Khufu) has a bearing of
**39.6 degrees** from true north.

Orion's Belt axis angle changes with precession over a ~25,772-year cycle.
The best-fit epoch — where the Belt's position angle most closely matches
the pyramid axis — is **11750 BCE**, with a residual offset of just
**0.0 degrees**.

| Epoch | Angular Offset |
|-------|----------------|
| Best fit (11750 BCE) | 0.0° |
| Construction (~2560 BCE) | 33.8° |
| Bauval's 10,500 BCE | 10.9° |

## Alignment Windows

- Epochs within **5 degrees** of match: **1,225 years** (p = 0.0475)
- Epochs within **10 degrees** of match: **2,400 years** (p = 0.0931)

## Combined Statistical Significance

Using Fisher's method to combine the shape match (p = 0.007) with the
orientation match (p = 0.0475):

- **Fisher chi-squared:** 16.016
- **Combined p-value:** 0.002997

## Interpretation

The best-fit orientation epoch (11750 BCE) falls within
approximately 1,250 years of Bauval's proposed
10,500 BCE date, which itself shows an offset of 10.9 degrees.
This is within a reasonable tolerance.

Critically, the construction-era epoch (~2560 BCE) shows a large angular
offset of 33.8 degrees, indicating that the Belt's
orientation at the time of construction did NOT match the pyramid axis.
This means the pyramid layout cannot reflect the sky as it appeared to
the builders — the shape match occurs at one epoch but the orientation
match occurs at a very different epoch.

The combined shape + orientation p-value (0.0030) is significant
at the 1% level, but the best-fit epoch predates construction by
approximately 9,190 years, making
intentional design an extraordinary claim requiring extraordinary evidence.

## Method Notes

- Pyramid axis: forward azimuth from Menkaure to Khufu (geodesic bearing).
- Belt position angle: PA of Mintaka-to-Alnitak line, computed at each epoch
  after applying proper motion and precession (astropy FK5).
- Epoch sweep: -13000 to 2025, step 25 years (602 points).
- Shape p-value (0.007) from existing Procrustes analysis.
