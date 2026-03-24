#!/usr/bin/env python3
"""
Analysis A: Nome Capital Geometry
=================================
Tests whether ancient Egyptian nome capitals cluster near the Great Circle.

Steps:
1. Compile 42 nome capital positions (scholarly consensus locations)
2. Compute distances to the Great Circle
3. Monte Carlo: 10,000 random great circles through Egypt
4. Sub-test: Upper vs. Lower Egypt
"""
import sys, os, math, random, json, csv
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'archive', 'great-circle-analysis', 'analysis'))
from utils import haversine_km, gc_distance, POLE_LAT, POLE_LON, EARTH_R_KM, QUARTER_CIRC, save_json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ============================================================
# NOME CAPITAL DATABASE
# ============================================================
# 42 nomes with capital identifications from Baines & Málek (2000),
# Helck (1974), and Pleiades/Trismegistos cross-references.
# Coordinates from modern town identifications.
#
# Upper Egypt: 22 nomes (south to north)
# Lower Egypt: 20 nomes (Memphis region and Delta)

NOME_CAPITALS = [
    # === UPPER EGYPT (South to North) ===
    # Nome 1 (Ta-Seti) - Elephantine/Abu
    {"nome": "UE-01", "name": "Elephantine (Abu)", "capital": "Elephantine", "lat": 24.085, "lon": 32.886, "region": "upper"},
    # Nome 2 (Throne of Horus) - Apollonopolis Magna/Edfu
    {"nome": "UE-02", "name": "Edfu (Behdet)", "capital": "Apollonopolis Magna", "lat": 24.978, "lon": 32.873, "region": "upper"},
    # Nome 3 (Shrine) - Nekhen/Hierakonpolis then El-Kab
    {"nome": "UE-03", "name": "Nekhen/Hierakonpolis", "capital": "Hierakonpolis", "lat": 25.098, "lon": 32.779, "region": "upper"},
    # Nome 4 (Sceptre) - Thebes/Waset
    {"nome": "UE-04", "name": "Thebes (Waset)", "capital": "Thebes", "lat": 25.719, "lon": 32.601, "region": "upper"},
    # Nome 5 (Two Falcons) - Koptos/Gebtu
    {"nome": "UE-05", "name": "Koptos (Gebtu)", "capital": "Coptos", "lat": 25.997, "lon": 32.813, "region": "upper"},
    # Nome 6 (Crocodile) - Iuny/Dendera
    {"nome": "UE-06", "name": "Dendera (Iunet)", "capital": "Tentyris", "lat": 26.142, "lon": 32.670, "region": "upper"},
    # Nome 7 (Sistrum) - Diospolis Parva/Hut-sekhem
    {"nome": "UE-07", "name": "Diospolis Parva (Hu)", "capital": "Diospolis Parva", "lat": 26.241, "lon": 32.291, "region": "upper"},
    # Nome 8 (Great Land) - Thinis/Abydos
    {"nome": "UE-08", "name": "Thinis/Abydos", "capital": "Thinis", "lat": 26.419, "lon": 31.950, "region": "upper"},
    # Nome 9 (Min) - Akhmim/Ipu
    {"nome": "UE-09", "name": "Akhmim (Ipu/Panopolis)", "capital": "Panopolis", "lat": 26.562, "lon": 31.749, "region": "upper"},
    # Nome 10 (Cobra) - Tjebu/Qaw el-Kebir
    {"nome": "UE-10", "name": "Tjebu (Antaeopolis)", "capital": "Antaeopolis", "lat": 26.886, "lon": 31.509, "region": "upper"},
    # Nome 11 (Seth) - Shashotep/Shutb
    {"nome": "UE-11", "name": "Shashotep (Hypselis)", "capital": "Hypselis", "lat": 27.039, "lon": 31.263, "region": "upper"},
    # Nome 12 (Viper Mountain) - Per-Nemty/Antinoë then Deir el-Gebrawi
    {"nome": "UE-12", "name": "Per-Nemty (Hieracon)", "capital": "Hieracon", "lat": 27.180, "lon": 31.044, "region": "upper"},
    # Nome 13 (Viper before) - Asyut/Zawty
    {"nome": "UE-13", "name": "Asyut (Zawty/Lycopolis)", "capital": "Lycopolis", "lat": 27.191, "lon": 31.172, "region": "upper"},
    # Nome 14 (Hathor below) - Qis/Cusae
    {"nome": "UE-14", "name": "Cusae (Qis)", "capital": "Cusae", "lat": 27.440, "lon": 30.905, "region": "upper"},
    # Nome 15 (Hare) - Hermopolis Magna/Khemenu
    {"nome": "UE-15", "name": "Hermopolis Magna (Khemenu)", "capital": "Hermopolis Magna", "lat": 27.783, "lon": 30.800, "region": "upper"},
    # Nome 16 (Oryx) - Hebenu/Beni Hasan area
    {"nome": "UE-16", "name": "Hebenu (Kom el-Ahmar)", "capital": "Hebenu", "lat": 28.039, "lon": 30.883, "region": "upper"},
    # Nome 17 (Jackal) - Saka/Kynopolis
    {"nome": "UE-17", "name": "Kynopolis (Saka)", "capital": "Cynopolis", "lat": 28.210, "lon": 30.762, "region": "upper"},
    # Nome 18 (Anti) - Hut-Nesut/Hipponon
    {"nome": "UE-18", "name": "Hut-Nesut (Hipponon)", "capital": "Hipponon", "lat": 28.390, "lon": 30.697, "region": "upper"},
    # Nome 19 (Two Sceptres) - Per-Medjed/Oxyrhynchus
    {"nome": "UE-19", "name": "Oxyrhynchus (Per-Medjed)", "capital": "Oxyrhynchus", "lat": 28.535, "lon": 30.663, "region": "upper"},
    # Nome 20 (Southern Sycamore) - Henen-Nesut/Herakleopolis Magna
    {"nome": "UE-20", "name": "Herakleopolis Magna (Henen-nesut)", "capital": "Herakleopolis Magna", "lat": 29.081, "lon": 30.932, "region": "upper"},
    # Nome 21 (Northern Sycamore) - Semenhor/Crocodilopolis
    {"nome": "UE-21", "name": "Crocodilopolis (Semenhor)", "capital": "Crocodilopolis", "lat": 29.314, "lon": 30.840, "region": "upper"},
    # Nome 22 (Knife) - Atfih/Aphroditopolis
    {"nome": "UE-22", "name": "Aphroditopolis (Atfih)", "capital": "Aphroditopolis", "lat": 29.571, "lon": 31.252, "region": "upper"},

    # === LOWER EGYPT (Memphis and Delta, south to north) ===
    # Nome 1 (White Wall) - Memphis/Ineb-Hedj
    {"nome": "LE-01", "name": "Memphis (Ineb-Hedj)", "capital": "Memphis", "lat": 29.845, "lon": 31.255, "region": "lower"},
    # Nome 2 (Thigh) - Letopolis/Khem (modern Ausim)
    {"nome": "LE-02", "name": "Letopolis (Khem)", "capital": "Letopolis", "lat": 30.108, "lon": 31.140, "region": "lower"},
    # Nome 3 (West) - Apis/Kom el-Hisn
    {"nome": "LE-03", "name": "Apis (Imu)", "capital": "Apis", "lat": 30.780, "lon": 30.595, "region": "lower"},
    # Nome 4 (Southern Shield) - Prosopitis/modern Tanta area
    {"nome": "LE-04", "name": "Prosopitis", "capital": "Prosopitis", "lat": 30.628, "lon": 30.851, "region": "lower"},
    # Nome 5 (Northern Shield) - Sais/Sa el-Hagar
    {"nome": "LE-05", "name": "Sais (Zau)", "capital": "Sais", "lat": 30.965, "lon": 30.769, "region": "lower"},
    # Nome 6 (Mountain Bull) - Xois/Sakha
    {"nome": "LE-06", "name": "Xois (Khaset)", "capital": "Xois", "lat": 31.094, "lon": 30.949, "region": "lower"},
    # Nome 7 (West Harpoon) - Per-aha/Damanhur area (Hermopolis Parva/Metelis)
    {"nome": "LE-07", "name": "Metelis (Per-aha)", "capital": "Metelis", "lat": 31.030, "lon": 30.470, "region": "lower"},
    # Nome 8 (East Harpoon) - Per-Atum/Pithom (Tell el-Maskhuta)
    {"nome": "LE-08", "name": "Pithom (Per-Atum)", "capital": "Pithom", "lat": 30.553, "lon": 31.879, "region": "lower"},
    # Nome 9 (Andjeti) - Busiris/Per-Usir (Abu Sir Bana)
    {"nome": "LE-09", "name": "Busiris (Per-Usir)", "capital": "Busiris", "lat": 30.911, "lon": 31.135, "region": "lower"},
    # Nome 10 (Black Bull) - Athribis/Hut-Hery-Ib (Tell Atrib)
    {"nome": "LE-10", "name": "Athribis (Hut-Hery-Ib)", "capital": "Athribis", "lat": 30.393, "lon": 31.194, "region": "lower"},
    # Nome 11 (Heseb-cow) - Leontopolis/Taremu (Tell el-Muqdam)
    {"nome": "LE-11", "name": "Leontopolis (Taremu)", "capital": "Leontopolis", "lat": 30.737, "lon": 31.367, "region": "lower"},
    # Nome 12 (Calf and Cow) - Tjeb-Netjer/Sebennytos
    {"nome": "LE-12", "name": "Sebennytos (Tjeb-Netjer)", "capital": "Sebennytos", "lat": 30.946, "lon": 31.118, "region": "lower"},
    # Nome 13 (Prospering Sceptre) - Heliopolis/Iunu
    {"nome": "LE-13", "name": "Heliopolis (Iunu)", "capital": "Heliopolis", "lat": 30.131, "lon": 31.313, "region": "lower"},
    # Nome 14 (Foremost of the East) - Tjaru/Sile (near Tell Hebua)
    {"nome": "LE-14", "name": "Tjaru (Sile)", "capital": "Sile", "lat": 30.931, "lon": 32.237, "region": "lower"},
    # Nome 15 (Ibis) - Hermopolis Parva/Ba'h (el-Baqliya)
    {"nome": "LE-15", "name": "Hermopolis Parva (Ba'h)", "capital": "Hermopolis Parva", "lat": 30.919, "lon": 31.364, "region": "lower"},
    # Nome 16 (Fish) - Djedet/Mendes (Tell el-Rub'a)
    {"nome": "LE-16", "name": "Mendes (Djedet)", "capital": "Mendes", "lat": 30.961, "lon": 31.518, "region": "lower"},
    # Nome 17 (Throne) - Semabehdet/Diospolis Inferior (Tell el-Balamun)
    {"nome": "LE-17", "name": "Diospolis Inferior (Semabehdet)", "capital": "Diospolis Inferior", "lat": 31.257, "lon": 31.553, "region": "lower"},
    # Nome 18 (Prince of the South) - Per-Bastet/Bubastis (Tell Basta)
    {"nome": "LE-18", "name": "Bubastis (Per-Bastet)", "capital": "Bubastis", "lat": 30.571, "lon": 31.515, "region": "lower"},
    # Nome 19 (Prince of the North) - Dja'net/Tanis (San el-Hagar)
    {"nome": "LE-19", "name": "Tanis (Dja'net)", "capital": "Tanis", "lat": 30.976, "lon": 31.881, "region": "lower"},
    # Nome 20 (Plumed Falcon) - Per-Sopdu (Saft el-Hinna)
    {"nome": "LE-20", "name": "Per-Sopdu (Saft el-Hinna)", "capital": "Per-Sopdu", "lat": 30.555, "lon": 31.637, "region": "lower"},
]

assert len(NOME_CAPITALS) == 42, f"Expected 42 nomes, got {len(NOME_CAPITALS)}"


def gc_distance_to_pole(lat, lon, pole_lat, pole_lon):
    """Distance from a point to the great circle defined by its pole."""
    return abs(haversine_km(lat, lon, pole_lat, pole_lon) - QUARTER_CIRC)


def random_pole_through_egypt():
    """Generate a random great circle pole such that the circle passes through Egypt.

    Constraint: the circle must cross the Nile Valley between 24°N and 31°N.
    We generate random poles and accept those whose circle passes through
    the Egypt latitude band (24-31°N) and longitude band (29-33°E).
    """
    while True:
        # Random pole on the sphere
        z = random.uniform(-1, 1)
        pole_lat = math.degrees(math.asin(z))
        pole_lon = random.uniform(-180, 180)

        # Check: does this circle pass through Egypt?
        # Sample points on the circle and check if any are in Egypt
        hit = False
        for bearing in range(0, 360, 5):
            lat_r = math.radians(pole_lat)
            lon_r = math.radians(pole_lon)
            d = math.pi / 2
            brng = math.radians(bearing)
            lat2 = math.asin(math.sin(lat_r) * math.cos(d) +
                             math.cos(lat_r) * math.sin(d) * math.cos(brng))
            lon2 = lon_r + math.atan2(math.sin(brng) * math.sin(d) * math.cos(lat_r),
                                       math.cos(d) - math.sin(lat_r) * math.sin(lat2))
            lat2_d = math.degrees(lat2)
            lon2_d = math.degrees(lon2)
            if 24 <= lat2_d <= 31 and 29 <= lon2_d <= 33:
                hit = True
                break
        if hit:
            return pole_lat, pole_lon


def main():
    print("=== Analysis A: Nome Capital Geometry ===\n")

    # Step 1: Compute distances
    print("Computing distances to Great Circle...")
    for nc in NOME_CAPITALS:
        nc['gc_dist_km'] = gc_distance_to_pole(nc['lat'], nc['lon'], POLE_LAT, POLE_LON)

    # Save CSV
    csv_path = os.path.join(OUTPUT_DIR, 'nome_capital_distances.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['nome', 'name', 'capital', 'lat', 'lon', 'region', 'gc_dist_km'])
        writer.writeheader()
        for nc in NOME_CAPITALS:
            writer.writerow({k: nc[k] for k in ['nome', 'name', 'capital', 'lat', 'lon', 'region', 'gc_dist_km']})
    print(f"  Saved: {csv_path}")

    # Step 2: Summary statistics
    dists = [nc['gc_dist_km'] for nc in NOME_CAPITALS]
    mean_dist = np.mean(dists)
    median_dist = np.median(dists)
    within_50 = sum(1 for d in dists if d <= 50)
    within_100 = sum(1 for d in dists if d <= 100)

    print(f"\n--- Observed Statistics ---")
    print(f"  Mean distance: {mean_dist:.1f} km")
    print(f"  Median distance: {median_dist:.1f} km")
    print(f"  Within 50 km: {within_50}/42 ({within_50/42*100:.1f}%)")
    print(f"  Within 100 km: {within_100}/42 ({within_100/42*100:.1f}%)")

    # Upper vs Lower
    upper = [nc for nc in NOME_CAPITALS if nc['region'] == 'upper']
    lower = [nc for nc in NOME_CAPITALS if nc['region'] == 'lower']
    upper_mean = np.mean([nc['gc_dist_km'] for nc in upper])
    lower_mean = np.mean([nc['gc_dist_km'] for nc in lower])
    print(f"\n  Upper Egypt ({len(upper)} nomes): mean dist = {upper_mean:.1f} km")
    print(f"  Lower Egypt ({len(lower)} nomes): mean dist = {lower_mean:.1f} km")

    # Step 3: Monte Carlo — random great circles through Egypt
    print(f"\nRunning Monte Carlo (10,000 random great circles through Egypt)...")
    random.seed(42)
    n_trials = 10000
    mc_means = []
    mc_within_50 = []
    mc_within_100 = []

    for i in range(n_trials):
        if (i + 1) % 1000 == 0:
            print(f"  Trial {i+1}/{n_trials}...")
        p_lat, p_lon = random_pole_through_egypt()
        qc = EARTH_R_KM * math.pi / 2
        trial_dists = [abs(haversine_km(nc['lat'], nc['lon'], p_lat, p_lon) - qc) for nc in NOME_CAPITALS]
        mc_means.append(np.mean(trial_dists))
        mc_within_50.append(sum(1 for d in trial_dists if d <= 50))
        mc_within_100.append(sum(1 for d in trial_dists if d <= 100))

    mc_means = np.array(mc_means)
    mc_within_50 = np.array(mc_within_50)
    mc_within_100 = np.array(mc_within_100)

    # Percentile ranking
    mean_percentile = np.sum(mc_means >= mean_dist) / n_trials * 100
    p_value_mean = np.sum(mc_means <= mean_dist) / n_trials

    w50_percentile = np.sum(mc_within_50 >= within_50) / n_trials * 100
    p_value_50 = np.sum(mc_within_50 >= within_50) / n_trials

    w100_percentile = np.sum(mc_within_100 >= within_100) / n_trials * 100
    p_value_100 = np.sum(mc_within_100 >= within_100) / n_trials

    print(f"\n--- Monte Carlo Results ---")
    print(f"  Baseline mean distance: {np.mean(mc_means):.1f} ± {np.std(mc_means):.1f} km")
    print(f"  Observed mean distance: {mean_dist:.1f} km")
    print(f"  Percentile (lower is closer): {mean_percentile:.1f}%")
    print(f"  p-value (mean dist): {p_value_mean:.4f}")
    print(f"  Baseline within-50km: {np.mean(mc_within_50):.1f} ± {np.std(mc_within_50):.1f}")
    print(f"  Observed within-50km: {within_50}")
    print(f"  p-value (within 50km): {p_value_50:.4f}")
    print(f"  Baseline within-100km: {np.mean(mc_within_100):.1f} ± {np.std(mc_within_100):.1f}")
    print(f"  Observed within-100km: {within_100}")
    print(f"  p-value (within 100km): {p_value_100:.4f}")

    # Save results
    results = {
        'observed': {
            'mean_dist_km': round(mean_dist, 2),
            'median_dist_km': round(median_dist, 2),
            'within_50km': within_50,
            'within_100km': within_100,
            'n_nomes': 42,
        },
        'upper_egypt': {
            'n_nomes': len(upper),
            'mean_dist_km': round(upper_mean, 2),
        },
        'lower_egypt': {
            'n_nomes': len(lower),
            'mean_dist_km': round(lower_mean, 2),
        },
        'monte_carlo': {
            'n_trials': n_trials,
            'baseline_mean_dist': round(float(np.mean(mc_means)), 2),
            'baseline_std_dist': round(float(np.std(mc_means)), 2),
            'observed_mean_dist': round(mean_dist, 2),
            'percentile_rank': round(mean_percentile, 2),
            'p_value_mean_dist': round(p_value_mean, 4),
            'baseline_mean_within_50': round(float(np.mean(mc_within_50)), 2),
            'observed_within_50': within_50,
            'p_value_within_50': round(p_value_50, 4),
            'baseline_mean_within_100': round(float(np.mean(mc_within_100)), 2),
            'observed_within_100': within_100,
            'p_value_within_100': round(p_value_100, 4),
        }
    }
    save_json(results, os.path.join(OUTPUT_DIR, 'nome_monte_carlo.json'))
    print(f"\n  Saved: nome_monte_carlo.json")

    # Step 4: Map
    print("\nGenerating map...")
    fig, ax = plt.subplots(1, 1, figsize=(10, 14))

    # Plot circle path through Egypt
    circle_lats = []
    circle_lons = []
    for bearing in np.linspace(0, 360, 3600):
        lat_r = math.radians(POLE_LAT)
        lon_r = math.radians(POLE_LON)
        d = math.pi / 2
        brng = math.radians(bearing)
        lat2 = math.asin(math.sin(lat_r) * math.cos(d) +
                         math.cos(lat_r) * math.sin(d) * math.cos(brng))
        lon2 = lon_r + math.atan2(math.sin(brng) * math.sin(d) * math.cos(lat_r),
                                   math.cos(d) - math.sin(lat_r) * math.sin(lat2))
        lat2_d = math.degrees(lat2)
        lon2_d = math.degrees(lon2)
        # Normalize longitude
        lon2_d = ((lon2_d + 180) % 360) - 180
        if 22 <= lat2_d <= 32 and 28 <= lon2_d <= 35:
            circle_lats.append(lat2_d)
            circle_lons.append(lon2_d)

    ax.plot(circle_lons, circle_lats, 'r-', linewidth=2, alpha=0.7, label='Great Circle', zorder=3)

    # Nile course (approximate)
    nile_lats = [24.08, 24.5, 25.0, 25.3, 25.7, 26.0, 26.3, 26.6, 27.0, 27.5,
                 28.0, 28.5, 29.0, 29.3, 29.5, 29.85, 30.1, 30.4, 30.7, 31.0, 31.5]
    nile_lons = [32.90, 32.88, 32.85, 32.78, 32.60, 32.72, 32.30, 31.90, 31.50,
                 31.15, 30.90, 30.70, 30.90, 30.85, 31.10, 31.25, 31.20, 31.15,
                 31.00, 31.00, 31.50]
    ax.plot(nile_lons, nile_lats, 'b-', linewidth=1.5, alpha=0.5, label='Nile (approx.)', zorder=2)

    # Plot nome capitals
    for nc in NOME_CAPITALS:
        color = '#2166ac' if nc['region'] == 'upper' else '#b2182b'
        marker = 's' if nc['region'] == 'upper' else 'D'
        ax.scatter(nc['lon'], nc['lat'], c=color, marker=marker, s=50, zorder=4, edgecolors='black', linewidth=0.5)
        ax.annotate(nc['nome'], (nc['lon'], nc['lat']), fontsize=5, ha='left',
                    xytext=(4, 2), textcoords='offset points')

    # Dummy points for legend
    ax.scatter([], [], c='#2166ac', marker='s', s=50, label=f'Upper Egypt ({len(upper)} nomes)', edgecolors='black', linewidth=0.5)
    ax.scatter([], [], c='#b2182b', marker='D', s=50, label=f'Lower Egypt ({len(lower)} nomes)', edgecolors='black', linewidth=0.5)

    ax.set_xlim(28, 35)
    ax.set_ylim(22, 32)
    ax.set_xlabel('Longitude (°E)')
    ax.set_ylabel('Latitude (°N)')
    ax.set_title(f'42 Nome Capitals and the Great Circle\nMean dist: {mean_dist:.1f} km | p = {p_value_mean:.4f}')
    ax.legend(loc='upper left', fontsize=8)
    ax.set_aspect(1 / math.cos(math.radians(27)))
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'nome_map.png'), dpi=150)
    plt.close()
    print(f"  Saved: nome_map.png")

    # Distance histogram
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.hist(mc_means, bins=60, alpha=0.7, color='gray', label='Random GC mean distance')
    ax.axvline(mean_dist, color='red', linewidth=2, linestyle='--', label=f'Observed: {mean_dist:.1f} km')
    ax.set_xlabel('Mean distance of 42 nome capitals to great circle (km)')
    ax.set_ylabel('Count (of 10,000 random circles)')
    ax.set_title(f'Monte Carlo: Nome Capital Proximity\np = {p_value_mean:.4f} (percentile: {mean_percentile:.1f}%)')
    ax.legend()
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'nome_monte_carlo_hist.png'), dpi=150)
    plt.close()
    print(f"  Saved: nome_monte_carlo_hist.png")

    print("\n=== Analysis A Complete ===")
    return results


if __name__ == '__main__':
    main()
