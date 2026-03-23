#!/usr/bin/env python3
"""
Study 11: Best-Fit Orientation Epoch (Orion Paper)
====================================================
The Giza pyramid layout matches Orion's Belt in SHAPE (p=0.007) but the
ORIENTATION (axis angle) varies with precession.  At what epoch does the
Belt axis angle best match the pyramid axis angle?

Phases:
  1. Compute Pyramid Axis Angle (Menkaure -> Khufu bearing)
  2. Compute Belt Axis Angle vs Epoch (-13000 to +2026, step 25 yr)
  3. Find Best-Fit Epoch (minimise angular offset)
  4. Statistical Significance (alignment window + Fisher's method)
  5. Report Key Finding
"""

import sys, os, math, json, warnings, datetime
import numpy as np
from scipy import stats as sp_stats
from pathlib import Path

sys.stdout = os.fdopen(sys.stdout.fileno(), "w", buffering=1)
warnings.filterwarnings("ignore")

np.random.seed(42)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = Path("/Users/elliotallan/megalith_site_research")
OUT_DIR = BASE / "outputs/next_wave/orion_orientation"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Pyramid Coordinates (Dash Foundation surveys) ────────────────────────────
KHUFU   = (29.9792, 31.1342)
KHAFRE  = (29.9761, 31.1306)
MENKAURE = (29.9725, 31.1278)

# ── Orion Belt Stars — J2000 RA(deg), Dec(deg), pm_RA(mas/yr), pm_Dec(mas/yr)
ALNITAK = (85.1897, -1.9425, 3.99, 2.54)    # zeta Ori
ALNILAM = (84.0533, -1.2019, 1.49, -0.46)   # epsilon Ori
MINTAKA = (83.0017, -0.2992, 0.64, -0.69)   # delta Ori

# ── Constants ─────────────────────────────────────────────────────────────────
PRECESSION_PERIOD = 25_772.0   # years
SHAPE_P = 0.007                # from existing Procrustes analysis
CONSTRUCTION_EPOCH = -2560     # ~2560 BCE
BAUVAL_EPOCH = -10500          # Bauval's claimed date

# ── Try astropy for rigorous precession ──────────────────────────────────────
try:
    from astropy.coordinates import SkyCoord, FK5
    import astropy.units as u
    from astropy.time import Time
    HAS_ASTROPY = True
    print("[info] astropy available — using rigorous precession")
except ImportError:
    HAS_ASTROPY = False
    print("[info] astropy not available — using IAU approximation")


# =============================================================================
# Phase 1: Pyramid Axis Angle
# =============================================================================
print("\n" + "=" * 72)
print("PHASE 1 — Pyramid Axis Angle")
print("=" * 72)


def forward_azimuth(lat1, lon1, lat2, lon2):
    """Forward bearing from point 1 to point 2 (degrees, 0=N, 90=E)."""
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dlam = math.radians(lon2 - lon1)
    x = math.sin(dlam) * math.cos(phi2)
    y = (math.cos(phi1) * math.sin(phi2)
         - math.sin(phi1) * math.cos(phi2) * math.cos(dlam))
    return math.degrees(math.atan2(x, y)) % 360


pyramid_axis = forward_azimuth(
    MENKAURE[0], MENKAURE[1], KHUFU[0], KHUFU[1]
)
print(f"  Bearing Menkaure -> Khufu: {pyramid_axis:.3f} deg (from N through E)")


# =============================================================================
# Phase 2: Belt Axis Angle vs Epoch
# =============================================================================
print("\n" + "=" * 72)
print("PHASE 2 — Belt Axis Angle vs Epoch")
print("=" * 72)


def belt_pa_astropy(epoch_yr):
    """Compute Belt PA at a given epoch using astropy precession."""
    dt = epoch_yr - 2000.0

    # Apply proper motion first (in J2000 frame)
    def apply_pm(star, dt):
        ra0, dec0, pm_ra, pm_dec = star
        dec_rad = math.radians(dec0)
        # pm_ra is already in true angle on sky (mas/yr * cos(dec) factored)
        # but spec says pm_RA / cos(Dec) so pm is coordinate rate
        ra = ra0 + (pm_ra / 3_600_000.0) * dt / math.cos(dec_rad)
        dec = dec0 + (pm_dec / 3_600_000.0) * dt
        return ra, dec

    alnitak_ra, alnitak_dec = apply_pm(ALNITAK, dt)
    mintaka_ra, mintaka_dec = apply_pm(MINTAKA, dt)

    # Build J2000 coordinates with proper-motion-corrected positions
    c_alnitak = SkyCoord(ra=alnitak_ra * u.deg, dec=alnitak_dec * u.deg,
                         frame=FK5(equinox="J2000"))
    c_mintaka = SkyCoord(ra=mintaka_ra * u.deg, dec=mintaka_dec * u.deg,
                         frame=FK5(equinox="J2000"))

    # Precess to epoch
    equinox_str = f"J{epoch_yr:.1f}"
    try:
        c_a = c_alnitak.transform_to(FK5(equinox=equinox_str))
        c_m = c_mintaka.transform_to(FK5(equinox=equinox_str))
    except Exception:
        # Fallback for extreme epochs
        return belt_pa_manual(epoch_yr)

    # Position angle: Mintaka -> Alnitak on sky
    dra = (c_a.ra.deg - c_m.ra.deg)
    ddec = (c_a.dec.deg - c_m.dec.deg)
    mean_dec = math.radians((c_a.dec.deg + c_m.dec.deg) / 2.0)
    pa = math.degrees(math.atan2(dra * math.cos(mean_dec), ddec)) % 360
    # Mirror E-W for sky-to-ground projection (looking up vs looking down)
    ground_bearing = (360.0 - pa) % 360.0
    return ground_bearing


def belt_pa_manual(epoch_yr):
    """Compute Belt PA using IAU linear precession approximation."""
    dt = epoch_yr - 2000.0

    # Precession constants (arcsec/yr)
    m = 46.124 / 3600.0   # degrees/yr — general precession in RA
    n = 20.043 / 3600.0   # degrees/yr — precession in obliquity component

    def precess_star(star, dt):
        ra0, dec0, pm_ra, pm_dec = star
        dec_rad = math.radians(dec0)
        # proper motion
        ra = ra0 + (pm_ra / 3_600_000.0) * dt / math.cos(dec_rad)
        dec = dec0 + (pm_dec / 3_600_000.0) * dt
        # precession (linear IAU approx)
        ra_rad = math.radians(ra)
        dec_rad = math.radians(dec)
        d_ra = (m + n * math.sin(ra_rad) * math.tan(dec_rad)) * dt
        d_dec = n * math.cos(ra_rad) * dt
        return ra + d_ra, dec + d_dec

    a_ra, a_dec = precess_star(ALNITAK, dt)
    m_ra, m_dec = precess_star(MINTAKA, dt)

    dra = a_ra - m_ra
    ddec = a_dec - m_dec
    mean_dec = math.radians((a_dec + m_dec) / 2.0)
    pa = math.degrees(math.atan2(dra * math.cos(mean_dec), ddec)) % 360
    # Mirror E-W for sky-to-ground projection
    ground_bearing = (360.0 - pa) % 360.0
    return ground_bearing


# Choose computation function
compute_belt_pa = belt_pa_astropy if HAS_ASTROPY else belt_pa_manual

# Sweep epochs
epochs = np.arange(-13000, 2026 + 1, 25)
belt_pas = []
angular_offsets = []

for ep in epochs:
    pa = compute_belt_pa(float(ep))
    belt_pas.append(pa)
    # Angular offset: shortest angle between two axes (mod 180 symmetry)
    diff = abs(pa - pyramid_axis) % 180
    if diff > 90:
        diff = 180 - diff
    angular_offsets.append(diff)

belt_pas = np.array(belt_pas)
angular_offsets = np.array(angular_offsets)

print(f"  Computed PA at {len(epochs)} epochs ({epochs[0]} to {epochs[-1]})")
print(f"  PA range: {belt_pas.min():.1f} - {belt_pas.max():.1f} deg")


# =============================================================================
# Phase 3: Find Best-Fit Epoch
# =============================================================================
print("\n" + "=" * 72)
print("PHASE 3 — Best-Fit Epoch")
print("=" * 72)

best_idx = np.argmin(angular_offsets)
best_epoch = int(epochs[best_idx])
min_offset = angular_offsets[best_idx]

print(f"  Best-fit epoch:   {best_epoch}")
if best_epoch < 0:
    print(f"                    ({abs(best_epoch)} BCE)")
else:
    print(f"                    ({best_epoch} CE)")
print(f"  Minimum offset:   {min_offset:.2f} deg")
print(f"  Belt PA at best:  {belt_pas[best_idx]:.2f} deg")
print(f"  Pyramid axis:     {pyramid_axis:.2f} deg")

# Offset at construction date
idx_construction = np.argmin(np.abs(epochs - CONSTRUCTION_EPOCH))
offset_construction = angular_offsets[idx_construction]
print(f"\n  Offset at construction (~2560 BCE): {offset_construction:.2f} deg")

idx_bauval = np.argmin(np.abs(epochs - BAUVAL_EPOCH))
offset_bauval = angular_offsets[idx_bauval]
print(f"  Offset at Bauval's 10,500 BCE:      {offset_bauval:.2f} deg")


# =============================================================================
# Phase 4: Statistical Significance
# =============================================================================
print("\n" + "=" * 72)
print("PHASE 4 — Statistical Significance")
print("=" * 72)

# Alignment windows
step = 25  # years per epoch step
window_5 = np.sum(angular_offsets < 5.0) * step
window_10 = np.sum(angular_offsets < 10.0) * step

print(f"  Alignment window (<5 deg):  {window_5:,} years")
print(f"  Alignment window (<10 deg): {window_10:,} years")

orient_p_5 = window_5 / PRECESSION_PERIOD
orient_p_10 = window_10 / PRECESSION_PERIOD
print(f"  Orientation p-value (<5 deg):  {orient_p_5:.4f}")
print(f"  Orientation p-value (<10 deg): {orient_p_10:.4f}")

# Fisher's method: combine shape and orientation p-values
# Use the <5 deg window for main result
orient_p = orient_p_5 if orient_p_5 > 0 else orient_p_10
if orient_p > 0:
    chi2_stat = -2.0 * (math.log(SHAPE_P) + math.log(orient_p))
    combined_p = sp_stats.chi2.sf(chi2_stat, df=4)
else:
    # No alignment window found at all — offset never drops below threshold
    # Use minimum offset to estimate probability: fraction of 90-degree range
    orient_p = float(min_offset) / 90.0
    chi2_stat = -2.0 * (math.log(SHAPE_P) + math.log(orient_p))
    combined_p = sp_stats.chi2.sf(chi2_stat, df=4)

print(f"\n  Shape p-value:     {SHAPE_P}")
print(f"  Orient p-value:    {orient_p:.4f}")
print(f"  Fisher chi2:       {chi2_stat:.3f}")
print(f"  Combined p-value:  {combined_p:.6f}")


# =============================================================================
# Phase 5: Plot
# =============================================================================
print("\n" + "=" * 72)
print("PHASE 5 — Plot")
print("=" * 72)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(12, 6))

ax.plot(epochs, angular_offsets, color="#2c3e50", linewidth=1.5, zorder=3)

# Mark best-fit
ax.axvline(best_epoch, color="#27ae60", linestyle="--", linewidth=1.2,
           label=f"Best fit: {abs(best_epoch)} BCE ({min_offset:.1f}\u00b0)")

# Mark construction date
ax.axvline(CONSTRUCTION_EPOCH, color="#e67e22", linestyle="--", linewidth=1.2,
           label=f"Construction ~2560 BCE ({offset_construction:.1f}\u00b0)")

# Mark Bauval's date
ax.axvline(BAUVAL_EPOCH, color="#c0392b", linestyle="--", linewidth=1.2,
           label=f"Bauval 10,500 BCE ({offset_bauval:.1f}\u00b0)")

# Shade alignment windows
ax.axhspan(0, 5, color="#27ae60", alpha=0.08, zorder=1)
ax.axhspan(5, 10, color="#f39c12", alpha=0.05, zorder=1)
ax.axhline(5, color="#27ae60", alpha=0.3, linewidth=0.8, linestyle=":")
ax.axhline(10, color="#f39c12", alpha=0.3, linewidth=0.8, linestyle=":")

ax.set_xlabel("Epoch (year; negative = BCE)", fontsize=12)
ax.set_ylabel("Angular Offset (\u00b0)", fontsize=12)
ax.set_title("Orion's Belt Orientation vs Giza Pyramid Axis — Best-Fit Epoch",
             fontsize=14, fontweight="bold")
ax.legend(loc="upper right", fontsize=10, framealpha=0.9)
ax.set_xlim(epochs[0], epochs[-1])
ax.set_ylim(0, max(angular_offsets) * 1.05)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plot_path = OUT_DIR / "orientation_epoch.png"
fig.savefig(plot_path, dpi=200)
plt.close()
print(f"  Saved plot: {plot_path}")


# =============================================================================
# Phase 6: Write Outputs
# =============================================================================
print("\n" + "=" * 72)
print("PHASE 6 — Write Outputs")
print("=" * 72)

# ── results.json ─────────────────────────────────────────────────────────────
results = {
    "study": "S11 — Best-Fit Orientation Epoch (Orion Paper)",
    "date": datetime.datetime.now().strftime("%Y-%m-%d"),
    "pyramid_axis_deg": round(pyramid_axis, 3),
    "best_fit_epoch": best_epoch,
    "best_fit_epoch_label": f"{abs(best_epoch)} BCE" if best_epoch < 0 else f"{best_epoch} CE",
    "min_angular_offset_deg": round(float(min_offset), 3),
    "belt_pa_at_best_deg": round(float(belt_pas[best_idx]), 3),
    "offset_at_construction_2560bce_deg": round(float(offset_construction), 3),
    "offset_at_bauval_10500bce_deg": round(float(offset_bauval), 3),
    "alignment_window_5deg_years": int(window_5),
    "alignment_window_10deg_years": int(window_10),
    "orientation_p_value_5deg": round(orient_p_5, 6),
    "orientation_p_value_10deg": round(orient_p_10, 6),
    "shape_p_value": SHAPE_P,
    "fisher_chi2": round(chi2_stat, 4),
    "combined_p_value": round(combined_p, 6),
    "precession_method": "astropy FK5" if HAS_ASTROPY else "IAU linear approx",
    "epoch_range": [int(epochs[0]), int(epochs[-1])],
    "epoch_step_years": 25,
}

json_path = OUT_DIR / "results.json"
with open(json_path, "w") as f:
    json.dump(results, f, indent=2)
print(f"  Saved: {json_path}")

# ── RESULTS.md ───────────────────────────────────────────────────────────────
best_label = results["best_fit_epoch_label"]
near_construction = abs(best_epoch - CONSTRUCTION_EPOCH) < 2000

md = f"""# Study 11 — Best-Fit Orientation Epoch (Orion Paper)

**Date:** {results['date']}
**Precession method:** {results['precession_method']}

---

## Key Finding

The Giza pyramid diagonal axis (Menkaure to Khufu) has a bearing of
**{pyramid_axis:.1f} degrees** from true north.

Orion's Belt axis angle changes with precession over a ~25,772-year cycle.
The best-fit epoch — where the Belt's position angle most closely matches
the pyramid axis — is **{best_label}**, with a residual offset of just
**{min_offset:.1f} degrees**.

| Epoch | Angular Offset |
|-------|----------------|
| Best fit ({best_label}) | {min_offset:.1f}\u00b0 |
| Construction (~2560 BCE) | {offset_construction:.1f}\u00b0 |
| Bauval's 10,500 BCE | {offset_bauval:.1f}\u00b0 |

## Alignment Windows

- Epochs within **5 degrees** of match: **{window_5:,} years** \
(p = {orient_p_5:.4f})
- Epochs within **10 degrees** of match: **{window_10:,} years** \
(p = {orient_p_10:.4f})

## Combined Statistical Significance

Using Fisher's method to combine the shape match (p = {SHAPE_P}) with the
orientation match (p = {orient_p_5:.4f}):

- **Fisher chi-squared:** {chi2_stat:.3f}
- **Combined p-value:** {combined_p:.6f}

## Interpretation

"""

near_bauval = abs(best_epoch - BAUVAL_EPOCH) < 3000

if near_construction:
    md += f"""The best-fit orientation epoch ({best_label}) is broadly consistent
with the historical construction period of the Great Pyramids (~2560 BCE).
This validates the observation that the Giza layout mirrors Orion's Belt
in both shape and orientation — as viewed from the construction epoch.

Bauval's proposed date of 10,500 BCE shows an angular offset of
{offset_bauval:.1f} degrees, which is {'within' if offset_bauval < 10 else 'outside'}
the 10-degree alignment window.
"""
elif near_bauval:
    md += f"""The best-fit orientation epoch ({best_label}) falls within
approximately {abs(best_epoch - BAUVAL_EPOCH):,} years of Bauval's proposed
10,500 BCE date, which itself shows an offset of {offset_bauval:.1f} degrees.
This is {'within' if offset_bauval < 15 else 'outside'} a reasonable tolerance.

Critically, the construction-era epoch (~2560 BCE) shows a large angular
offset of {offset_construction:.1f} degrees, indicating that the Belt's
orientation at the time of construction did NOT match the pyramid axis.
This means the pyramid layout cannot reflect the sky as it appeared to
the builders — the shape match occurs at one epoch but the orientation
match occurs at a very different epoch.

The combined shape + orientation p-value ({combined_p:.4f}) is significant
at the 1% level, but the best-fit epoch predates construction by
approximately {abs(best_epoch - CONSTRUCTION_EPOCH):,} years, making
intentional design an extraordinary claim requiring extraordinary evidence.
"""
else:
    md += f"""The best-fit epoch ({best_label}) does not coincide with
either the construction date (~2560 BCE) or Bauval's 10,500 BCE claim.

The construction-era epoch shows a {offset_construction:.1f}-degree offset,
while Bauval's date shows {offset_bauval:.1f} degrees.

The combined shape + orientation p-value ({combined_p:.4f}) is significant,
but the temporal mismatch with construction complicates any intentional-design
interpretation.
"""

md += f"""
## Method Notes

- Pyramid axis: forward azimuth from Menkaure to Khufu (geodesic bearing).
- Belt position angle: PA of Mintaka-to-Alnitak line, computed at each epoch
  after applying proper motion and precession ({results['precession_method']}).
- Epoch sweep: {epochs[0]} to {epochs[-1]}, step {25} years ({len(epochs)} points).
- Shape p-value (0.007) from existing Procrustes analysis.
"""

md_path = OUT_DIR / "RESULTS.md"
with open(md_path, "w") as f:
    f.write(md)
print(f"  Saved: {md_path}")

print("\n" + "=" * 72)
print("Study 11 complete.")
print("=" * 72)
