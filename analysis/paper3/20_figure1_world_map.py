#!/usr/bin/env python3
"""
Figure 1: World map with Great Circle and 7 civilization origins.
Publication-quality Robinson projection.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os

POLE_LAT = 59.682122
POLE_LON = -138.646087

# Civilization origins
ORIGINS = {
    "Egypt": {"lat": 29.87, "lon": 31.22, "dist": 6.3, "on": True},
    "Mesopotamia": {"lat": 31.33, "lon": 45.64, "dist": 120.2, "on": True},
    "Indus": {"lat": 27.32, "lon": 68.14, "dist": 26.6, "on": True},
    "China": {"lat": 34.27, "lon": 108.94, "dist": 2123, "on": False},
    "Mesoamerica": {"lat": 19.69, "lon": -98.84, "dist": 4558, "on": False},
    "Andes": {"lat": -13.16, "lon": -72.55, "dist": 16.8, "on": True},
    "West Africa": {"lat": 11.50, "lon": -4.00, "dist": 1124, "on": False},
}

# Label offsets (dx, dy in data coords for Robinson)
LABEL_OFFSETS = {
    "Egypt": (-2, 7),
    "Mesopotamia": (3, -8),
    "Indus": (4, 6),
    "China": (4, 4),
    "Mesoamerica": (4, -6),
    "Andes": (-20, -7),
    "West Africa": (4, 5),
}

def to_xyz(lat, lon):
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    return np.array([np.cos(lat_r)*np.cos(lon_r), np.cos(lat_r)*np.sin(lon_r), np.sin(lat_r)])

def gc_points(pole_lat, pole_lon, n=3600):
    """Generate points along the great circle defined by its pole."""
    pole = to_xyz(pole_lat, pole_lon)
    if abs(pole[2]) < 0.9:
        ref = np.array([0, 0, 1])
    else:
        ref = np.array([1, 0, 0])
    u = np.cross(pole, ref)
    u /= np.linalg.norm(u)
    v = np.cross(pole, u)
    v /= np.linalg.norm(v)

    angles = np.linspace(0, 2*np.pi, n, endpoint=False)
    points = np.outer(np.cos(angles), u) + np.outer(np.sin(angles), v)

    lats = np.degrees(np.arcsin(np.clip(points[:, 2], -1, 1)))
    lons = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    return lats, lons

def main():
    # Generate great circle points
    gc_lats, gc_lons = gc_points(POLE_LAT, POLE_LON, 3600)

    # Sort by longitude for cleaner plotting — split into segments to avoid wrap-around lines
    # Group consecutive points and split where longitude jumps > 180°
    segments_lats = []
    segments_lons = []
    seg_lat = [gc_lats[0]]
    seg_lon = [gc_lons[0]]

    for i in range(1, len(gc_lats)):
        if abs(gc_lons[i] - gc_lons[i-1]) > 180:
            # Wrap-around: start new segment
            segments_lats.append(np.array(seg_lat))
            segments_lons.append(np.array(seg_lon))
            seg_lat = [gc_lats[i]]
            seg_lon = [gc_lons[i]]
        else:
            seg_lat.append(gc_lats[i])
            seg_lon.append(gc_lons[i])
    segments_lats.append(np.array(seg_lat))
    segments_lons.append(np.array(seg_lon))

    # Create figure
    fig = plt.figure(figsize=(10, 6), dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())

    # Map features
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='#E8E8E8', edgecolor='none', zorder=1)
    ax.add_feature(cfeature.OCEAN, facecolor='white', edgecolor='none', zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.3, color='#999999', zorder=2)

    # Plot great circle segments
    for seg_lat, seg_lon in zip(segments_lats, segments_lons):
        ax.plot(seg_lon, seg_lat,
                color='#1a5276', linewidth=1.5,
                transform=ccrs.Geodetic(), zorder=3)

    # Plot civilization origins
    for name, info in ORIGINS.items():
        color = '#27ae60' if info['on'] else '#c0392b'
        edgecolor = '#1a6b3a' if info['on'] else '#7b241c'

        ax.plot(info['lon'], info['lat'], 'o',
                color=color, markersize=7, markeredgecolor=edgecolor,
                markeredgewidth=0.8, transform=ccrs.PlateCarree(), zorder=5)

        # Label
        dx, dy = LABEL_OFFSETS[name]
        dist_str = f"{info['dist']:.0f} km" if info['dist'] >= 100 else f"{info['dist']:.1f} km"
        label = f"{name}\n({dist_str})"

        ax.text(info['lon'] + dx, info['lat'] + dy, label,
                fontsize=5.5, fontweight='normal', color='#2c3e50',
                transform=ccrs.PlateCarree(), zorder=6,
                ha='left', va='center',
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                         edgecolor='none', alpha=0.8))

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#27ae60',
               markeredgecolor='#1a6b3a', markersize=7, label='Within 200 km'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#c0392b',
               markeredgecolor='#7b241c', markersize=7, label='Beyond 200 km'),
        Line2D([0], [0], color='#1a5276', linewidth=1.5, label='Alison Great Circle'),
    ]
    ax.legend(handles=legend_elements, loc='lower left', fontsize=6,
              frameon=True, facecolor='white', edgecolor='#cccccc',
              framealpha=0.9, borderpad=0.5)

    plt.tight_layout()

    out_path = os.path.expanduser('~/megalith_site_research/paper3/figures/figure1_world_map_origins.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved to {out_path}")
    print(f"Size: {os.path.getsize(out_path)/1024:.0f} KB")

if __name__ == "__main__":
    main()
