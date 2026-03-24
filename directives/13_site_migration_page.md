# Directive 13 — Site: Migration Page & Figures Integration

**Project:** thegreatcircle.earth
**Repository path:** `/Users/elliotallan/megalith_site_research/site`
**Author:** Elliot Allan
**Priority:** High — supports Paper 2 submission
**Depends on:** All 12 research directives complete

---

## Objective

Add a new **"Migration"** page to the site that presents the 60,000-year corridor narrative from Paper 2, and integrate the five paper figures into the existing Results page. This directive is self-contained — do not modify any analysis code or data files.

---

## Context You Need to Know

The site at `thegreatcircle.earth` currently has five pages:

| Page | File | Nav label |
|------|------|-----------|
| Globe explorer | `index.html` | Globe |
| Methodology | `methodology.html` | Methodology |
| Results | `results.html` | Results |
| Paper (companion) | `paper.html` | Paper |
| Timeline | `timeline.html` | Timeline |

The site uses **Mapbox GL JS v3.3.0** for all interactive maps. The Mapbox token (already in globe.js) is:
```
MAPBOX_TOKEN_REDACTED
```

**Design tokens** (copy exactly — do not deviate):
```css
--bg: #08080D;
--bg-card: #12121E;
--bg-card-hover: #1A1A2E;
--text: #E8E6E1;
--text-muted: #8A8A9A;
--accent: #C9A84C;
--accent-light: #E2C76A;
--accent-dim: rgba(201, 168, 76, .15);
--accent-teal: #14B8A6;
--border: #252538;
--border-accent: rgba(201, 168, 76, .25);
--font-heading: 'Cinzel', 'Georgia', serif;
--font: 'Space Grotesk', -apple-system, sans-serif;
--mono: 'Space Mono', 'Consolas', monospace;
```

**Google Fonts already used site-wide:**
```html
<link href="https://fonts.googleapis.com/css2?family=Cinzel:wght@400;600;700;900&family=Space+Grotesk:wght@300;400;500;600;700&family=Space+Mono:wght@400;700&display=swap" rel="stylesheet">
```

**Nav pattern:** Every page has a `<nav>` with the same structure. The active page gets `class="active"` on its `<a>` tag. Copy the nav exactly from `index.html`, changing only which link has `class="active"` and adding the new Migration link.

---

## Tasks

### Task 1 — Copy Figures to Site Images Folder

Copy the five paper figures from:
```
/Users/elliotallan/The Great Circle/figure1_world_map.png
/Users/elliotallan/The Great Circle/figure2_divergence.png
/Users/elliotallan/The Great Circle/figure3_temporal.png
/Users/elliotallan/The Great Circle/figure4_egypt.png
/Users/elliotallan/The Great Circle/figure5_ripleys_k.png
```

Into:
```
/Users/elliotallan/megalith_site_research/site/images/
```

Name them exactly as above (figure1_world_map.png through figure5_ripleys_k.png).

---

### Task 2 — Create `migration.html`

Create `/Users/elliotallan/megalith_site_research/site/migration.html`.

This is a **full-screen interactive Mapbox map** (same approach as `timeline.html`) that shows the 60,000-year corridor. Below is the full specification.

#### 2a. Page structure

```
<head>
  - GA4 tag (copy from any existing page)
  - meta charset, viewport, OG tags for migration page
  - og:title = "Migration — The Great Circle"
  - og:description = "The Great Circle traces the Out-of-Africa southern dispersal route — a corridor first walked 60,000 years ago that attracted continuous human activity through every subsequent epoch."
  - og:url = "https://thegreatcircle.earth/migration.html"
  - Google Fonts (same link as other pages)
  - Mapbox GL JS CSS + JS (v3.3.0, same as timeline.html)
  - Inline <style> block (see 2b)
</head>
<body>
  - Minimal nav bar (same structure as timeline.html's .tl-nav)
  - #map div (full screen)
  - Epoch sidebar panel (left, see 2c)
  - Stats panel (right, see 2d)
  - Timeline controls (bottom, see 2e)
  - Site info popup (see 2f)
  - <script> block (see 2g)
</body>
```

#### 2b. CSS

Full-screen layout — map fills the viewport. Nav floats on top. Use the same approach as `timeline.html`. Key rules:

```css
html, body { height: 100%; overflow: hidden; background: var(--bg); }
#map { width: 100%; height: 100%; }
```

The nav (`.tl-nav`) should match `timeline.html`'s nav exactly (logo left, title center, back-link right).

**Epoch panel** (left sidebar, 280px wide, overlapping the map):
```css
.epoch-panel {
  position: absolute; top: 52px; left: 0; bottom: 140px;
  width: 280px; z-index: 15;
  background: rgba(8, 8, 13, .88); backdrop-filter: blur(12px);
  border-right: 1px solid var(--border);
  overflow-y: auto; padding: 1rem;
}
```

**Epoch button** (one per epoch — 5 total plus "All"):
```css
.epoch-btn {
  display: block; width: 100%; text-align: left;
  padding: .6rem .8rem; margin-bottom: .35rem;
  background: transparent; border: 1px solid var(--border);
  color: var(--text-muted); cursor: pointer;
  font-family: var(--font); font-size: .8rem;
  transition: all .2s; border-radius: 0;
}
.epoch-btn:hover { border-color: var(--accent); color: var(--text); }
.epoch-btn.active { background: var(--accent-dim); border-color: var(--accent); color: var(--accent); }
.epoch-btn .epoch-dot { display: inline-block; width: 8px; height: 8px; border-radius: 50%; margin-right: .5rem; }
.epoch-btn .epoch-stat { font-family: var(--mono); font-size: .7rem; color: var(--accent); display: block; margin-top: .2rem; padding-left: 1rem; }
```

**Stats card** (right side, fixed width 220px):
```css
.stats-card {
  position: absolute; top: 60px; right: 12px; width: 220px; z-index: 15;
  background: rgba(8, 8, 13, .88); backdrop-filter: blur(12px);
  border: 1px solid var(--border-accent); padding: 1rem;
}
.stats-card h4 { font-family: var(--font-heading); font-size: .75rem; color: var(--text-muted); letter-spacing: .12em; text-transform: uppercase; margin-bottom: .75rem; }
.stat-row { display: flex; justify-content: space-between; padding: .3rem 0; border-bottom: 1px solid var(--border); font-size: .75rem; }
.stat-row:last-child { border-bottom: none; }
.stat-val { font-family: var(--mono); color: var(--accent); }
```

**Bottom timeline controls** — same structure as `timeline.html`'s `.tl-controls`.

**Site popup** (appears on click):
```css
.site-popup {
  position: absolute; z-index: 25; min-width: 240px; max-width: 300px;
  background: var(--bg-card); border: 1px solid var(--border-accent);
  padding: 1rem 1.1rem; pointer-events: auto;
}
.site-popup h3 { font-family: var(--font-heading); font-size: .95rem; color: var(--accent); margin-bottom: .4rem; }
.site-popup .tag { display: inline-block; font-size: .7rem; font-family: var(--mono); padding: .15rem .5rem; margin-bottom: .5rem; }
.site-popup p { font-size: .78rem; color: var(--text-muted); line-height: 1.5; }
.site-popup .close-btn { position: absolute; top: .5rem; right: .6rem; background: none; border: none; color: var(--text-muted); cursor: pointer; font-size: 1rem; }
```

#### 2c. Epoch Sidebar Content

Six buttons (in order):

| Button label | Epoch key | Dot colour | Stat line |
|---|---|---|---|
| All Epochs | `all` | `#C9A84C` | 60,000 BP → Present |
| Out-of-Africa (>50k BP) | `ooa` | `#8B0000` | Z = 4.42 · p = 0.001 |
| Upper Paleolithic (50–25k BP) | `upper` | `#CC4400` | Z = 2.21 · marginal |
| Early Holocene (12–8k BP) | `holo` | `#FF8C00` | +3.91σ campsite signal |
| Agricultural Origins (8–5k BP) | `agri` | `#DAA520` | 9 of 24 centers · p = 0.002 |
| Bronze Age Monuments (5–2k BP) | `bronze` | `#1a6e3c` | 5.05× enrichment · Z = 11.83 |

Below the buttons, add a text block:

```html
<div style="margin-top:1.2rem; font-size:.75rem; color:var(--text-muted); line-height:1.6; border-top:1px solid var(--border); padding-top:.8rem;">
  <p>The corridor shows human occupation in <span style="color:var(--accent)">5 of 6 tested epochs</span> spanning 130,000 BP to the present. The single gap (LGM, 26–14k BP) is expected — Arabia was impassable during peak glaciation.</p>
  <p style="margin-top:.6rem;">Bronze Age monuments sit <span style="color:var(--accent)">~55× closer</span> to early Holocene campsites than chance (p = 0.004), confirming the corridor was inherited, not invented.</p>
</div>
```

Add a "Play Epochs" button below:
```html
<button id="btn-play" class="epoch-btn" style="margin-top:.8rem;border-color:var(--accent);color:var(--accent);">▶ Animate Epochs</button>
```

#### 2d. Stats Panel Content

```html
<div class="stats-card">
  <h4>Key Statistics</h4>
  <div class="stat-row"><span>Deep-time enrichment</span><span class="stat-val">Z = 4.42</span></div>
  <div class="stat-row"><span>Monument enrichment</span><span class="stat-val">5.05×</span></div>
  <div class="stat-row"><span>Settlement (null)</span><span class="stat-val">0.78×</span></div>
  <div class="stat-row"><span>Agricultural origins</span><span class="stat-val">p = 0.002</span></div>
  <div class="stat-row"><span>Ahramat intersection</span><span class="stat-val">p = 0.00016</span></div>
  <div class="stat-row"><span>Bronze Age continuity</span><span class="stat-val">~55× closer</span></div>
  <div class="stat-row"><span>Route fit (southern)</span><span class="stat-val">p = 0.042</span></div>
  <div class="stat-row"><span>Epochs occupied</span><span class="stat-val">5 of 6</span></div>
</div>
```

#### 2e. Timeline Controls

Match `timeline.html` bottom controls exactly. Strip the year slider (not needed here — epoch filter replaces it). Replace with:
- Large epoch label (current epoch name in Cinzel gold)
- One-line description of what this epoch shows
- Paper 2 citation: `Allan (2026) · From Migration to Monumentality · Preprint`

#### 2f. Site Data

Define the sites array in the `<script>` block. Use this exact data:

```javascript
const MIGRATION_SITES = [
  // OUT-OF-AFRICA >50,000 BP
  { name: "Mololo Cave", lat: -6.0, lon: 145.5, age_bp: 55000, epoch: "ooa", dist_km: 12,
    desc: "Papua New Guinea · 55,000 BP · 12 km from circle. Earliest confirmed human occupation in PNG. Sits on the corridor's eastern terminus — the same arc as the Out-of-Africa southern dispersal." },
  { name: "Madjedbebe", lat: -12.45, lon: 132.89, age_bp: 65000, epoch: "ooa", dist_km: 1107,
    desc: "Australia · 65,000 BP. Oldest confirmed site in Australia. Entry point from Sahul during glacial low sea levels." },
  { name: "Jebel Faya", lat: 25.21, lon: 56.02, age_bp: 125000, epoch: "ooa", dist_km: 473,
    desc: "UAE · 125,000 BP. Earliest confirmed modern human site outside Africa. Arabian Peninsula crossing during MIS 5 humid phase." },
  { name: "Jwalapuram", lat: 15.37, lon: 78.13, age_bp: 74000, epoch: "ooa", dist_km: 1031,
    desc: "India · 74,000 BP. Continuity through the Toba super-eruption. Key node on the southern India coastal corridor." },
  { name: "16R Dune (Thar Desert)", lat: 24.8, lon: 73.72, age_bp: 96000, epoch: "ooa", dist_km: 159,
    desc: "Rajasthan, India · 96,000 BP · 159 km from circle. Within 200 km of the Great Circle." },
  { name: "Tam Pa Ling", lat: 20.22, lon: 103.4, age_bp: 63000, epoch: "ooa", dist_km: 487,
    desc: "Laos · 63,000 BP. Key node on the southern SE Asia arc of the dispersal route." },
  { name: "Al Wusta", lat: 26.0, lon: 56.1, age_bp: 85000, epoch: "ooa", dist_km: 385,
    desc: "Arabia · 85,000 BP. Confirms Green Arabia dispersal during MIS 5 humid phase." },

  // UPPER PALEOLITHIC 50–25k BP
  { name: "Kebar Valley", lat: -6.5, lon: 146.2, age_bp: 26000, epoch: "upper", dist_km: 16,
    desc: "Papua New Guinea · 26,000 BP · 16 km from circle. Persistent corridor occupation into Upper Paleolithic." },
  { name: "That Nang Ing", lat: 18.31, lon: 103.75, age_bp: 40000, epoch: "upper", dist_km: 313,
    desc: "Laos · 40,000 BP. On the mainland SE Asia corridor." },
  { name: "Khorat Plateau", lat: 14.5, lon: 100.0, age_bp: 40000, epoch: "upper", dist_km: 243,
    desc: "Thailand · 40,000 BP · 243 km from circle." },
  { name: "Niah Cave", lat: 3.82, lon: 113.77, age_bp: 40000, epoch: "upper", dist_km: 604,
    desc: "Borneo · 40,000 BP. Important early occupation site in island SE Asia." },

  // EARLY HOLOCENE 12–8k BP
  { name: "Beidha", lat: 30.5, lon: 35.5, age_bp: 10000, epoch: "holo", dist_km: 21,
    desc: "Jordan · 10,000 BP · 21 km from circle. Pre-Pottery Neolithic A waystation. Material: charcoal-dominant (transit camp, no grain)." },
  { name: "Ghwair I", lat: 30.3, lon: 35.7, age_bp: 9500, epoch: "holo", dist_km: 48,
    desc: "Jordan · 9,500 BP · 48 km from circle. PPNA/B campsite cluster at the Levant–Arabia junction." },
  { name: "Caution Bay", lat: -9.5, lon: 147.2, age_bp: 7000, epoch: "holo", dist_km: 45,
    desc: "Papua New Guinea · 7,000 BP · 45 km. Corridor occupation continues at the eastern terminus into the Holocene." },

  // AGRICULTURAL ORIGINS 8–5k BP
  { name: "Wheat/Barley Origin (Levant)", lat: 36.0, lon: 37.0, age_bp: 9000, epoch: "agri", dist_km: 201,
    desc: "Levant · 9,000 BP · 201 km from circle. Independent origin of wheat and barley — one of 9 domestication centers within 500 km." },
  { name: "Zebu Cattle (Indus)", lat: 24.0, lon: 68.0, age_bp: 8000, epoch: "agri", dist_km: 64,
    desc: "South Asia · 8,000 BP · 64 km — closest major domestication center to the circle." },
  { name: "Taro Origin (SE Asia)", lat: 10.5, lon: 104.0, age_bp: 7000, epoch: "agri", dist_km: 280,
    desc: "Southeast Asia · 7,000 BP. One of 9 independent domestication centers within 500 km of the circle." },
  { name: "Potato / Llama (Andes)", lat: -13.5, lon: -71.9, age_bp: 7000, epoch: "agri", dist_km: 192,
    desc: "Peru · 7,000 BP · 192 km from circle. Independent Andean agricultural origin." },

  // BRONZE AGE MONUMENTS 5–2k BP
  { name: "Giza Pyramids", lat: 29.979, lon: 31.134, age_bp: 2560, epoch: "bronze", dist_km: 6.5,
    desc: "Egypt · 2,560 BCE · 6.5 km from circle. Ripley's K analysis shows pyramids are arranged along the circle's axis (ratio 0.08, p < 0.0001). Ahramat Branch intersection p = 0.00016." },
  { name: "Saqqara / Memphis", lat: 29.872, lon: 31.216, age_bp: 2700, epoch: "bronze", dist_km: 3.2,
    desc: "Egypt · 2,700 BCE · Step Pyramid of Djoser. Memphis: where the Great Circle crosses the Nile. Compound probability p = 0.003." },
  { name: "Dahshur", lat: 29.795, lon: 31.214, age_bp: 2600, epoch: "bronze", dist_km: 8.1,
    desc: "Egypt · 2,600 BCE. Red Pyramid + Bent Pyramid. Part of the 70-km pyramid strip aligned along the Great Circle axis." },
  { name: "Nazca Lines", lat: -14.73, lon: -75.13, age_bp: 1800, epoch: "bronze", dist_km: 18,
    desc: "Peru · 200 BCE–600 CE · 18 km from circle. Note: geoglyph orientation (63.1°) is confounded with June solstice sunrise (65.6°) — hypothesis is not independently testable." },
  { name: "Easter Island (Rapa Nui)", lat: -27.12, lon: -109.36, age_bp: 1100, epoch: "bronze", dist_km: 9.6,
    desc: "Pacific · settled ~900 CE · 9.6 km from circle. ~900 moai. Polynesian navigation analysis shows this island's alignment is local, not network-wide." },
  { name: "Persepolis", lat: 29.94, lon: 52.89, age_bp: 2500, epoch: "bronze", dist_km: 12,
    desc: "Iran · 515 BCE · 12 km from circle. Achaemenid ceremonial capital at the densest segment of the Old World corridor." },
  { name: "Mohenjo-daro", lat: 27.33, lon: 68.14, age_bp: 4500, epoch: "bronze", dist_km: 31,
    desc: "Pakistan · 2,500 BCE · 31 km. One of the Indus Valley Civilization's largest cities. On the corridor connecting the Levant to SE Asia." },
  { name: "Machu Picchu", lat: -13.16, lon: -72.54, age_bp: 600, epoch: "bronze", dist_km: 48,
    desc: "Peru · 1450 CE · 48 km from circle. Inca royal estate on the Andean arc." },
];

const EPOCH_COLORS = {
  all:    "#C9A84C",
  ooa:    "#8B0000",
  upper:  "#CC4400",
  holo:   "#FF8C00",
  agri:   "#DAA520",
  bronze: "#1a6e3c"
};

const EPOCH_INFO = {
  all:    { label: "All Epochs — 60,000 BP to Present",
            desc:  "The complete corridor, layered across deep time. Five of six tested epochs show significant enrichment." },
  ooa:    { label: ">50,000 BP — Out-of-Africa Dispersal",
            desc:  "Modern humans leave Africa on the southern coastal route. Deep-time enrichment Z = 4.42, p = 0.001. Mololo Cave (PNG) sits 12 km from the circle at 55,000 BP." },
  upper:  { label: "50,000–25,000 BP — Upper Paleolithic",
            desc:  "Persistent corridor use. Z = 2.21. The LGM gap (26–14k BP) follows — Arabia becomes impassable. Corridor occupation resumes when the climate ameliorates." },
  holo:   { label: "12,000–8,000 BP — Early Holocene Campsites",
            desc:  "Green Arabia opens the southern corridor. Radiocarbon dates show charcoal-dominant material (+3.91σ) — transit camps, not settlements. Beidha (Jordan) sits 21 km from circle." },
  agri:   { label: "8,000–5,000 BP — Agricultural Origins",
            desc:  "9 of 24 independent domestication centers fall within 500 km. p = 0.002 vs. 10,000 random great circles. The corridor passes through the origins of wheat, cattle, and taro." },
  bronze: { label: "5,000–2,000 BP — Bronze Age Monuments",
            desc:  "Monument enrichment 5.05× (Z = 11.83). Settlements: 0.78× (null). Bronze Age sites sit ~55× closer to Holocene campsites than chance (p = 0.004). The corridor was inherited, not designed." }
};
```

#### 2g. JavaScript Implementation

Initialize Mapbox with:
```javascript
const map = new mapboxgl.Map({
  container: 'map',
  style: 'mapbox://styles/mapbox/dark-v11',  // same as other pages presumably
  projection: 'globe',
  center: [60, 20],
  zoom: 1.5,
  pitch: 0,
  bearing: 0
});
```

**On map load:**
1. Load the Great Circle GeoJSON from `data/circle.json` and add as a red line layer (same color `#c0392b` as the existing globe uses — check `globe.js` for the exact layer style and replicate it).
2. Add a 50-km corridor band — a buffered version of the circle in red with low opacity (0.08).
3. Add site markers as circles using Mapbox GL's circle layer, one per epoch, color-coded by `EPOCH_COLORS`. Use `circle-radius` expressions scaling by epoch importance (monuments bigger).

**Epoch filter function:**
```javascript
function setEpoch(epochKey) {
  // Update button active states
  // Update bottom label and description text
  // For each site layer:
  //   if epochKey === 'all' → show all
  //   else → filter visibility to matching epoch only
  // Fly camera to relevant region:
  //   ooa:    { center: [115, 5],  zoom: 2.5 }
  //   upper:  { center: [120, 0],  zoom: 2.5 }
  //   holo:   { center: [40, 30],  zoom: 3.0 }
  //   agri:   { center: [60, 25],  zoom: 2.5 }
  //   bronze: { center: [35, 25],  zoom: 2.8 }
  //   all:    { center: [60, 20],  zoom: 1.5 }
}
```

**Click handler:** On clicking a site marker, show the `.site-popup` with the site's `name`, epoch tag (coloured by epoch), `age_bp`, `dist_km`, and `desc`. Position the popup at the click coordinates using Mapbox's `LngLat`.

**Animate button:**
Cycle through epochs `['ooa', 'upper', 'holo', 'agri', 'bronze', 'all']` at 3.5-second intervals, calling `setEpoch()` each step.

---

### Task 3 — Update Nav on All Existing Pages

Add `<li><a href="migration.html">Migration</a></li>` to the nav `<ul>` on **all five existing pages**, positioned between `Timeline` and the Substack icon.

On `migration.html` itself, give the Migration link `class="active"`.

The nav on `migration.html` should use the `.tl-nav` style (minimal floating bar, same as `timeline.html`) not the full nav bar, since this is a full-screen page.

---

### Task 4 — Add Figures Section to `results.html`

Append a new section to `results.html` **after the existing last section** and before `</main>` (or wherever the content ends). Title: **"Paper 2 — From Migration to Monumentality"**.

Use this structure:
```html
<section id="paper2-figures" style="margin-top:3rem;padding-top:2rem;border-top:1px solid var(--border);">
  <h2 style="font-family:var(--font-heading);...">Paper 2 Findings</h2>
  <p class="muted" style="margin-bottom:2rem;max-width:700px;">
    The following figures are from the companion study
    <em>From Migration to Monumentality: The Alison Great Circle as a 60,000-Year Human Corridor</em>
    (Allan 2026, in revision). All analyses use distribution-matched Monte Carlo simulation
    (10,000 iterations) against latitude-profile-matched random great circles.
    <a href="migration.html" style="color:var(--accent);">→ Explore the interactive migration map</a>
  </p>

  <!-- Figure grid: 2 columns where possible -->
  <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(520px,1fr));gap:2rem;">

    <figure>
      <img src="images/figure1_world_map.png" alt="Figure 1: World map showing the Great Circle corridor with sites colour-coded by epoch"
           style="width:100%;border:1px solid var(--border);display:block;">
      <figcaption style="font-size:.78rem;color:var(--text-muted);margin-top:.5rem;font-family:var(--mono);">
        Figure 1. The Alison Great Circle corridor with sites colour-coded by epoch from Out-of-Africa dispersal (>50,000 BP) to Bronze Age monuments.
      </figcaption>
    </figure>

    <figure>
      <img src="images/figure2_divergence.png" alt="Figure 2: Monument-settlement divergence across seven databases"
           style="width:100%;border:1px solid var(--border);display:block;">
      <figcaption style="font-size:.78rem;color:var(--text-muted);margin-top:.5rem;font-family:var(--mono);">
        Figure 2. Monument-settlement divergence replicated across seven independent databases. Monuments: 4–5× enrichment. Settlements: ~0.8× (below random).
      </figcaption>
    </figure>

    <figure style="grid-column:1/-1;">
      <img src="images/figure3_temporal.png" alt="Figure 3: Temporal layering of corridor activity 60,000 BP to Bronze Age"
           style="width:100%;border:1px solid var(--border);display:block;">
      <figcaption style="font-size:.78rem;color:var(--text-muted);margin-top:.5rem;font-family:var(--mono);">
        Figure 3. Temporal layering of corridor activity, 60,000 BP to Bronze Age. ★ = statistically significant enrichment. Grey band = Last Glacial Maximum gap (expected).
      </figcaption>
    </figure>

    <figure>
      <img src="images/figure4_egypt.png" alt="Figure 4: Egypt close-up with Ahramat Branch intersection"
           style="width:100%;border:1px solid var(--border);display:block;">
      <figcaption style="font-size:.78rem;color:var(--text-muted);margin-top:.5rem;font-family:var(--mono);">
        Figure 4. Great Circle–Ahramat Branch intersection in the Memphis necropolis. Compound probability p = 0.00016.
      </figcaption>
    </figure>

    <figure>
      <img src="images/figure5_ripleys_k.png" alt="Figure 5: Directional Ripley's K clustering analysis"
           style="width:100%;border:1px solid var(--border);display:block;">
      <figcaption style="font-size:.78rem;color:var(--text-muted);margin-top:.5rem;font-family:var(--mono);">
        Figure 5. Directional Ripley's K-function. Pyramids are arranged along the Great Circle axis (parallel/perpendicular ratio = 0.08, p < 0.0001).
      </figcaption>
    </figure>

  </div>
</section>
```

---

### Task 5 — Verify

After completing all tasks:

1. Open `migration.html` in a browser. Confirm:
   - Map loads without console errors
   - Great Circle line renders in red
   - All 26 sites appear as coloured markers
   - Clicking "Out-of-Africa" epoch filters to 7 sites and flies to SE Asia
   - Clicking a marker shows the popup with correct text
   - "Animate Epochs" button cycles through all 6 epochs automatically
   - Nav links to all other pages work correctly

2. Open `results.html`. Confirm:
   - All 5 figures render correctly
   - No layout breaks on mobile (check at 375px width)

3. Check all 5 existing HTML files have the Migration nav link.

4. Confirm the figure images are present in `site/images/`.

---

## Files Modified

| File | Action |
|------|--------|
| `site/migration.html` | **CREATE** |
| `site/results.html` | **MODIFY** (append figures section) |
| `site/index.html` | **MODIFY** (add Migration nav link) |
| `site/methodology.html` | **MODIFY** (add Migration nav link) |
| `site/paper.html` | **MODIFY** (add Migration nav link) |
| `site/timeline.html` | **MODIFY** (add Migration nav link) |
| `site/images/figure1_world_map.png` | **COPY** |
| `site/images/figure2_divergence.png` | **COPY** |
| `site/images/figure3_temporal.png` | **COPY** |
| `site/images/figure4_egypt.png` | **COPY** |
| `site/images/figure5_ripleys_k.png` | **COPY** |

---

## Do Not

- Do not modify `js/globe.js`, any existing JS, or any data files in `site/data/`
- Do not change the design of existing pages (only add the nav link and the figures section)
- Do not create a separate CSS file — inline all styles for `migration.html` as done in `timeline.html`
- Do not use Three.js — use Mapbox GL JS exclusively, matching the existing site stack

---

*Directive written by Ell · March 2026*
