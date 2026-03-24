#!/usr/bin/env python3
"""
Analysis D: Predynastic Site Distribution
==========================================
Tests whether the Great Circle's monument-settlement divergence
extends back before the Old Kingdom into Predynastic Egypt.

Steps:
1. Compile Predynastic Egyptian sites from Pleiades + supplementary data
2. Classify as monument vs. settlement
3. Compute enrichment by period
4. Monte Carlo significance testing
5. Temporal gradient analysis
"""
import sys, os, math, random, json, csv
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'archive', 'great-circle-analysis', 'analysis'))
from utils import (haversine_km, gc_distance, POLE_LAT, POLE_LON, EARTH_R_KM,
                   QUARTER_CIRC, save_json, random_pole, load_pleiades,
                   classify_pleiades)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'data')

# Period definitions (BCE, negative years in Pleiades)
PERIODS = {
    'Badarian': (-4400, -4000),
    'Naqada I': (-4000, -3500),
    'Naqada II': (-3500, -3200),
    'Naqada III': (-3200, -3000),
    'Early Dynastic': (-3100, -2686),
    'Old Kingdom': (-2686, -2181),
}

# Supplementary Predynastic sites not in Pleiades
# Major sites from Hendrickx (1999/2006) and published excavation reports
SUPPLEMENTARY_PREDYNASTIC = [
    # Badarian culture sites
    {"name": "Badari", "lat": 26.98, "lon": 31.42, "type": "cemetery", "period": "Badarian", "min_date": -4400},
    {"name": "Mostagedda", "lat": 26.94, "lon": 31.40, "type": "cemetery", "period": "Badarian", "min_date": -4400},
    {"name": "Matmar", "lat": 27.02, "lon": 31.36, "type": "cemetery", "period": "Badarian", "min_date": -4400},
    {"name": "Hammamiya", "lat": 26.97, "lon": 31.44, "type": "settlement", "period": "Badarian", "min_date": -4400},
    {"name": "Qau el-Kebir", "lat": 26.89, "lon": 31.51, "type": "cemetery", "period": "Badarian", "min_date": -4400},

    # Naqada I (Amratian) sites
    {"name": "Naqada (main site)", "lat": 25.90, "lon": 32.72, "type": "cemetery", "period": "Naqada I", "min_date": -4000},
    {"name": "Hierakonpolis HK29A", "lat": 25.10, "lon": 32.76, "type": "ceremonial", "period": "Naqada I", "min_date": -4000},
    {"name": "Hierakonpolis HK6", "lat": 25.11, "lon": 32.77, "type": "cemetery", "period": "Naqada I", "min_date": -3800},
    {"name": "El-Mahâsna", "lat": 26.38, "lon": 31.91, "type": "cemetery", "period": "Naqada I", "min_date": -3900},
    {"name": "Abadiya", "lat": 26.00, "lon": 32.60, "type": "cemetery", "period": "Naqada I", "min_date": -3900},
    {"name": "Hu (Diospolis Parva)", "lat": 26.24, "lon": 32.29, "type": "cemetery", "period": "Naqada I", "min_date": -3900},
    {"name": "Armant", "lat": 25.62, "lon": 32.54, "type": "settlement", "period": "Naqada I", "min_date": -3900},
    {"name": "Gebelein", "lat": 25.52, "lon": 32.49, "type": "settlement", "period": "Naqada I", "min_date": -3900},

    # Naqada II (Gerzean) sites
    {"name": "Gerzeh", "lat": 29.38, "lon": 31.07, "type": "cemetery", "period": "Naqada II", "min_date": -3500},
    {"name": "Tarkhan", "lat": 29.50, "lon": 31.18, "type": "cemetery", "period": "Naqada II", "min_date": -3400},
    {"name": "Harageh", "lat": 29.27, "lon": 31.10, "type": "cemetery", "period": "Naqada II", "min_date": -3400},
    {"name": "Maadi", "lat": 29.95, "lon": 31.30, "type": "settlement", "period": "Naqada II", "min_date": -3500},
    {"name": "Wadi Digla", "lat": 29.93, "lon": 31.32, "type": "settlement", "period": "Naqada II", "min_date": -3500},
    {"name": "Naqada South Town", "lat": 25.89, "lon": 32.71, "type": "settlement", "period": "Naqada II", "min_date": -3500},
    {"name": "Abydos Cemetery U", "lat": 26.40, "lon": 31.94, "type": "cemetery", "period": "Naqada II", "min_date": -3400},
    {"name": "Adaima", "lat": 25.23, "lon": 32.75, "type": "settlement", "period": "Naqada II", "min_date": -3500},
    {"name": "El-Amra", "lat": 26.32, "lon": 32.05, "type": "cemetery", "period": "Naqada II", "min_date": -3500},
    {"name": "Buto (Tell el-Fara'in)", "lat": 31.20, "lon": 30.74, "type": "settlement", "period": "Naqada II", "min_date": -3500},
    {"name": "Tell el-Farkha", "lat": 30.78, "lon": 31.37, "type": "settlement", "period": "Naqada II", "min_date": -3400},
    {"name": "Minshat Abu Omar", "lat": 30.93, "lon": 31.95, "type": "cemetery", "period": "Naqada II", "min_date": -3400},
    {"name": "Helwan necropolis", "lat": 29.84, "lon": 31.34, "type": "cemetery", "period": "Naqada II", "min_date": -3200},

    # Naqada III / Dynasty 0 sites
    {"name": "Abydos Cemetery B (Royal)", "lat": 26.41, "lon": 31.95, "type": "ceremonial", "period": "Naqada III", "min_date": -3200},
    {"name": "Hierakonpolis Fort", "lat": 25.09, "lon": 32.76, "type": "ceremonial", "period": "Naqada III", "min_date": -3200},
    {"name": "Elephantine temple precinct", "lat": 24.09, "lon": 32.89, "type": "ceremonial", "period": "Naqada III", "min_date": -3200},
    {"name": "Koptos temple", "lat": 25.99, "lon": 32.81, "type": "ceremonial", "period": "Naqada III", "min_date": -3200},
    {"name": "Tell Ibrahim Awad", "lat": 30.72, "lon": 31.47, "type": "settlement", "period": "Naqada III", "min_date": -3200},

    # Early Dynastic
    {"name": "Saqqara (Dynasty 1 tombs)", "lat": 29.87, "lon": 31.22, "type": "ceremonial", "period": "Early Dynastic", "min_date": -3100},
    {"name": "Abydos (Dynasty 1 royal cemetery)", "lat": 26.42, "lon": 31.95, "type": "ceremonial", "period": "Early Dynastic", "min_date": -3100},
    {"name": "Helwan (Early Dynastic)", "lat": 29.85, "lon": 31.33, "type": "cemetery", "period": "Early Dynastic", "min_date": -3000},
    {"name": "Abu Rawash (Dynasty 1)", "lat": 30.03, "lon": 31.08, "type": "cemetery", "period": "Early Dynastic", "min_date": -2900},
    {"name": "Nag el-Deir", "lat": 26.48, "lon": 31.86, "type": "cemetery", "period": "Early Dynastic", "min_date": -3000},
]

# Classification: monumental/ceremonial vs. settlement
MONUMENTAL_TYPES = {'ceremonial', 'temple', 'pyramid', 'tomb', 'monument',
                     'sanctuary', 'architecturalcomplex'}
SETTLEMENT_TYPES = {'settlement', 'village', 'farmstead', 'findspot', 'port', 'mine'}
CEMETERY_TYPE = {'cemetery'}  # cemeteries are a middle category — test both ways


def classify_predynastic(type_str):
    """Classify a predynastic site as 'monumental', 'settlement', or 'cemetery'."""
    t = type_str.lower().strip()
    if t in MONUMENTAL_TYPES or 'ceremonial' in t or 'temple' in t:
        return 'monumental'
    elif t in SETTLEMENT_TYPES or 'settlement' in t or 'village' in t:
        return 'settlement'
    elif t in CEMETERY_TYPE or 'cemetery' in t:
        return 'cemetery'
    return 'other'


def assign_period(min_date):
    """Assign a site to a period based on its min_date."""
    if min_date is None:
        return None
    for name, (start, end) in PERIODS.items():
        if start <= min_date < end:
            return name
    if min_date < -4400:
        return 'Pre-Badarian'
    if min_date >= -2181:
        return 'Post-OK'
    return None


def compute_enrichment(sites, threshold_km=50):
    """Compute the enrichment ratio: (observed within threshold) / (expected).

    Returns enrichment, observed count, expected fraction.
    """
    if not sites:
        return 0, 0, 0
    dists = [gc_distance(s['lat'], s['lon']) for s in sites]
    observed = sum(1 for d in dists if d <= threshold_km)
    # Expected: fraction of Egypt's area within threshold of the circle
    # Egypt approx 1,000 x 600 km; circle corridor width = 2 * threshold
    # Rough: fraction = 2 * threshold / characteristic_width
    # Better: use the actual site distribution
    total = len(sites)
    if total == 0:
        return 0, 0, 0
    return observed, total, observed / total


def main():
    print("=== Analysis D: Predynastic Site Distribution ===\n")

    # Step 1: Load Pleiades Egyptian sites
    print("Loading Pleiades data...")
    pleiades_path = os.path.join(DATA_DIR, 'pleiades', 'pleiades-places-latest.csv')
    all_pleiades = load_pleiades(pleiades_path)

    # Filter to Egypt and Predynastic/Early Dynastic
    egypt_predynastic = []
    for s in all_pleiades:
        if not (22 <= s['lat'] <= 32 and 25 <= s['lon'] <= 35):
            continue
        if s.get('min_date') is not None and s['min_date'] <= -3000:
            # Classify
            pleiades_class = classify_pleiades(s.get('type', ''))
            s['site_class'] = pleiades_class
            s['period'] = assign_period(s['min_date'])
            s['source'] = 'pleiades'
            egypt_predynastic.append(s)

    print(f"  Pleiades pre-3000 BCE Egyptian sites: {len(egypt_predynastic)}")

    # Step 2: Add supplementary sites
    for s in SUPPLEMENTARY_PREDYNASTIC:
        s['site_class'] = classify_predynastic(s['type'])
        if s.get('period') is None:
            s['period'] = assign_period(s['min_date'])
        s['source'] = 'supplement'
        egypt_predynastic.append(s)

    # Deduplicate by proximity (within 2km)
    deduped = []
    for s in egypt_predynastic:
        is_dup = False
        for d in deduped:
            if haversine_km(s['lat'], s['lon'], d['lat'], d['lon']) < 2:
                is_dup = True
                break
        if not is_dup:
            deduped.append(s)

    sites = deduped
    print(f"  After deduplication: {len(sites)} sites")

    # Compute GC distances
    for s in sites:
        s['gc_dist_km'] = gc_distance(s['lat'], s['lon'])

    # Save CSV
    csv_path = os.path.join(OUTPUT_DIR, 'predynastic_sites.csv')
    fields = ['name', 'lat', 'lon', 'site_class', 'period', 'min_date', 'gc_dist_km', 'source']
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for s in sites:
            writer.writerow({k: s.get(k, '') for k in fields})
    print(f"  Saved: {csv_path}")

    # Step 3: Overall statistics
    print(f"\n--- Overall Statistics ---")
    classes = {}
    for s in sites:
        c = s.get('site_class', 'other')
        classes.setdefault(c, []).append(s)

    for c, ss in sorted(classes.items()):
        dists = [s['gc_dist_km'] for s in ss]
        within_50 = sum(1 for d in dists if d <= 50)
        print(f"  {c}: {len(ss)} sites, {within_50} within 50km ({within_50/len(ss)*100:.1f}%), "
              f"mean dist: {np.mean(dists):.1f} km")

    # Step 4: Enrichment by period
    print(f"\n--- Enrichment by Period ---")
    threshold = 50  # km

    period_results = {}
    for period_name, (start, end) in PERIODS.items():
        period_sites = [s for s in sites if s.get('period') == period_name]
        if not period_sites:
            print(f"  {period_name}: no sites")
            period_results[period_name] = {'n_sites': 0}
            continue

        # All sites in period
        n_total = len(period_sites)
        n_within = sum(1 for s in period_sites if s['gc_dist_km'] <= threshold)
        frac = n_within / n_total if n_total > 0 else 0

        # Split by class
        mon_sites = [s for s in period_sites if s['site_class'] in ('monumental', 'ceremonial')]
        set_sites = [s for s in period_sites if s['site_class'] in ('settlement',)]
        cem_sites = [s for s in period_sites if s['site_class'] in ('cemetery',)]

        mon_within = sum(1 for s in mon_sites if s['gc_dist_km'] <= threshold) if mon_sites else 0
        set_within = sum(1 for s in set_sites if s['gc_dist_km'] <= threshold) if set_sites else 0
        cem_within = sum(1 for s in cem_sites if s['gc_dist_km'] <= threshold) if cem_sites else 0

        mon_frac = mon_within / len(mon_sites) if mon_sites else 0
        set_frac = set_within / len(set_sites) if set_sites else 0
        cem_frac = cem_within / len(cem_sites) if cem_sites else 0

        print(f"  {period_name} ({start} to {end} BCE):")
        print(f"    Total: {n_total} sites, {n_within} within {threshold}km ({frac*100:.1f}%)")
        print(f"    Monumental: {len(mon_sites)} sites, {mon_within} within ({mon_frac*100:.0f}%)")
        print(f"    Settlement: {len(set_sites)} sites, {set_within} within ({set_frac*100:.0f}%)")
        print(f"    Cemetery: {len(cem_sites)} sites, {cem_within} within ({cem_frac*100:.0f}%)")

        period_results[period_name] = {
            'n_sites': n_total,
            'n_within_50km': n_within,
            'fraction_within': round(frac, 4),
            'monumental': {'n': len(mon_sites), 'n_within': mon_within, 'frac': round(mon_frac, 4)},
            'settlement': {'n': len(set_sites), 'n_within': set_within, 'frac': round(set_frac, 4)},
            'cemetery': {'n': len(cem_sites), 'n_within': cem_within, 'frac': round(cem_frac, 4)},
        }

    # Step 5: Monte Carlo — random great circles
    print(f"\nRunning Monte Carlo (10,000 random circles)...")
    random.seed(42)
    n_trials = 10000

    # Compute observed enrichment (monumental vs settlement)
    all_mon = [s for s in sites if s['site_class'] in ('monumental', 'ceremonial')]
    all_set = [s for s in sites if s['site_class'] == 'settlement']

    obs_mon_within = sum(1 for s in all_mon if s['gc_dist_km'] <= threshold)
    obs_set_within = sum(1 for s in all_set if s['gc_dist_km'] <= threshold)
    obs_mon_frac = obs_mon_within / len(all_mon) if all_mon else 0
    obs_set_frac = obs_set_within / len(all_set) if all_set else 0
    obs_divergence = obs_mon_frac - obs_set_frac

    mc_divergence = []
    mc_mon_frac = []
    mc_set_frac = []

    qc = EARTH_R_KM * math.pi / 2
    for trial in range(n_trials):
        if (trial + 1) % 2000 == 0:
            print(f"  Trial {trial+1}/{n_trials}...")

        p_lat, p_lon = random_pole()
        mon_within = sum(1 for s in all_mon
                         if abs(haversine_km(s['lat'], s['lon'], p_lat, p_lon) - qc) <= threshold)
        set_within = sum(1 for s in all_set
                         if abs(haversine_km(s['lat'], s['lon'], p_lat, p_lon) - qc) <= threshold)
        mf = mon_within / len(all_mon) if all_mon else 0
        sf = set_within / len(all_set) if all_set else 0
        mc_mon_frac.append(mf)
        mc_set_frac.append(sf)
        mc_divergence.append(mf - sf)

    mc_divergence = np.array(mc_divergence)
    mc_mon_frac = np.array(mc_mon_frac)
    mc_set_frac = np.array(mc_set_frac)

    z_div = (obs_divergence - np.mean(mc_divergence)) / np.std(mc_divergence) if np.std(mc_divergence) > 0 else 0
    p_div = np.sum(mc_divergence >= obs_divergence) / n_trials

    z_mon = (obs_mon_frac - np.mean(mc_mon_frac)) / np.std(mc_mon_frac) if np.std(mc_mon_frac) > 0 else 0
    p_mon = np.sum(mc_mon_frac >= obs_mon_frac) / n_trials

    print(f"\n--- Monte Carlo Results ---")
    print(f"  Monumental sites within {threshold}km: {obs_mon_within}/{len(all_mon)} = {obs_mon_frac*100:.1f}%")
    print(f"    Baseline: {np.mean(mc_mon_frac)*100:.1f}% ± {np.std(mc_mon_frac)*100:.1f}%")
    print(f"    Z-score: {z_mon:.2f}, p = {p_mon:.4f}")
    print(f"  Settlement sites within {threshold}km: {obs_set_within}/{len(all_set)} = {obs_set_frac*100:.1f}%")
    print(f"  Divergence (mon - set): {obs_divergence*100:.1f}%")
    print(f"    Baseline: {np.mean(mc_divergence)*100:.1f}% ± {np.std(mc_divergence)*100:.1f}%")
    print(f"    Z-score: {z_div:.2f}, p = {p_div:.4f}")

    # Save enrichment results
    enrichment_results = {
        'by_period': period_results,
        'overall': {
            'n_monumental': len(all_mon),
            'n_settlement': len(all_set),
            'monumental_within_50km': obs_mon_within,
            'settlement_within_50km': obs_set_within,
            'monumental_fraction': round(obs_mon_frac, 4),
            'settlement_fraction': round(obs_set_frac, 4),
            'divergence': round(obs_divergence, 4),
        },
        'monte_carlo': {
            'n_trials': n_trials,
            'monumental_z': round(z_mon, 3),
            'monumental_p': round(p_mon, 4),
            'divergence_z': round(z_div, 3),
            'divergence_p': round(p_div, 4),
            'baseline_divergence_mean': round(float(np.mean(mc_divergence)), 4),
            'baseline_divergence_std': round(float(np.std(mc_divergence)), 4),
        }
    }
    save_json(enrichment_results, os.path.join(OUTPUT_DIR, 'predynastic_enrichment_by_period.json'))
    save_json({
        'monumental_fraction': round(obs_mon_frac, 4),
        'settlement_fraction': round(obs_set_frac, 4),
        'divergence': round(obs_divergence, 4),
        'divergence_z': round(z_div, 3),
        'divergence_p': round(p_div, 4),
    }, os.path.join(OUTPUT_DIR, 'predynastic_divergence.json'))
    print(f"\nSaved: predynastic_enrichment_by_period.json, predynastic_divergence.json")

    # Step 6: Timeline plot
    print("\nGenerating timeline plot...")
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    period_names = ['Badarian', 'Naqada I', 'Naqada II', 'Naqada III', 'Early Dynastic', 'Old Kingdom']
    period_midpoints = [-4200, -3750, -3350, -3100, -2893, -2434]
    period_widths = [400, 500, 300, 200, 414, 505]

    # Panel 1: Fraction within 50km by class
    ax = axes[0]
    mon_fracs = []
    set_fracs = []
    cem_fracs = []
    for pn in period_names:
        pr = period_results.get(pn, {})
        mon_fracs.append(pr.get('monumental', {}).get('frac', 0))
        set_fracs.append(pr.get('settlement', {}).get('frac', 0))
        cem_fracs.append(pr.get('cemetery', {}).get('frac', 0))

    x = np.array(period_midpoints)
    w = np.array(period_widths) * 0.25

    ax.bar(x - w, mon_fracs, width=w, color='#d73027', alpha=0.8, label='Monumental/Ceremonial')
    ax.bar(x, set_fracs, width=w, color='#4575b4', alpha=0.8, label='Settlement')
    ax.bar(x + w, cem_fracs, width=w, color='#fee090', alpha=0.8, label='Cemetery')
    ax.set_ylabel(f'Fraction within {threshold}km of Great Circle')
    ax.set_title('Predynastic to Old Kingdom: Site Proximity to Great Circle by Type')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    # Panel 2: Site counts
    ax = axes[1]
    total_counts = [period_results.get(pn, {}).get('n_sites', 0) for pn in period_names]
    within_counts = [period_results.get(pn, {}).get('n_within_50km', 0) for pn in period_names]

    ax.bar(x, total_counts, width=w * 2, color='lightgray', alpha=0.7, label='Total sites')
    ax.bar(x, within_counts, width=w * 2, color='#d73027', alpha=0.7, label=f'Within {threshold}km')
    ax.set_xlabel('Date (BCE)')
    ax.set_ylabel('Number of sites')
    ax.set_title('Site Counts by Period')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    # Add period labels
    for pn, mid in zip(period_names, period_midpoints):
        axes[1].text(mid, -0.5, pn, ha='center', va='top', fontsize=7, rotation=30)

    ax.invert_xaxis()
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'predynastic_timeline.png'), dpi=150)
    plt.close()
    print(f"Saved: predynastic_timeline.png")

    print("\n=== Analysis D Complete ===")
    return enrichment_results


if __name__ == '__main__':
    main()
