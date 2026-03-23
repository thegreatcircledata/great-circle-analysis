# Test A3: Construction-Epoch Orientation Significance

## Purpose
Is the Orion Belt orientation offset at the Giza construction epoch (2560 BCE)
unusual compared to what you'd expect for a random building axis, or compared
to other well-shape-matched star triplets?

## Data
- **Stars:** 156 bright stars (mag < 4.0)
- **Triplets:** 620,620 = C(156, 3)
- **Construction epoch:** -2560 (2560 BCE)
- **Pyramid axis:** 39.604° (Menkaure → Khufu bearing, mod 180 = 39.60°)
- **Orion Belt PA at construction:** 106.46°
- **Orientation offset:** 66.85°
- **Seed:** 42

## Results

### Part 1: Random Axis Monte Carlo (N = 10,000)
A uniform random axis has probability **0.7427** (74.3%) of falling
within 66.85° of the Orion Belt PA at 2560 BCE.

MC estimate: p = 0.7395.

**Not significant.** About 74% of random axes would be at least this close.

### Part 2: Orientation Among Top Shape-Matched Triplets
Among the top 6,530 triplets by Procrustes distance
(shape p ≤ 0.0105):

- Orion orientation rank: **5564 / 6530**
- p(orient | top shape) = **0.8521**

**Not significant.** Orion's orientation at 2560 BCE is unremarkable even among the best shape matches.

### Part 3: Orientation Rank Among ALL Triplets
- Rank: 515,707 / 620,620
- p = **0.8310**

### Combined Shape + Orientation
- p(shape) × p(orient | shape) = 0.010522 × 0.8521 = **0.008965**

## Conclusion
The construction-epoch orientation offset (66.85°) is **not independently significant**.
A random building axis has a ~74% chance of being within this angular distance
of any fixed reference direction. Among the top shape-matched triplets, Orion's
orientation at 2560 BCE does not stand out. The construction-epoch orientation
adds no meaningful evidence beyond the shape match alone.
