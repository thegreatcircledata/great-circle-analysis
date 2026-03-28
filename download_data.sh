#!/bin/bash
# download_data.sh — Download large data files excluded from the repository
#
# These files are too large for GitHub but are needed to reproduce some analyses.
# Run this script from the repository root after cloning.
#
# Total download size: ~6.5 GB
# Required disk space: ~10 GB (compressed + uncompressed)

set -e

echo "=== Great Circle Analysis: Large Data Download Script ==="
echo ""

# ----------------------------------------------------------------
# ETOPO1 Global Relief Model (Paper 2: desert edge analysis, Paper 4: LCP cost surfaces)
# Source: NOAA National Centers for Environmental Information
# ----------------------------------------------------------------
echo "[1/7] Downloading ETOPO1 global relief model (~890 MB)..."
mkdir -p data/geophysical/etopo
if [ ! -f data/geophysical/etopo/ETOPO1_Ice_g_gmt4.grd ]; then
    curl -L -o data/geophysical/etopo/ETOPO1_Ice_g_gmt4.grd \
        "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz" \
        && gunzip data/geophysical/etopo/ETOPO1_Ice_g_gmt4.grd.gz 2>/dev/null || true
    echo "  NOTE: If the NOAA URL has changed, download ETOPO1 manually from:"
    echo "  https://www.ncei.noaa.gov/products/etopo-global-relief-model"
else
    echo "  Already exists, skipping."
fi

# ----------------------------------------------------------------
# EMAG2 v3 Earth Magnetic Anomaly Grid (Paper 2: geophysical supplement)
# Source: NOAA/NCEI
# ----------------------------------------------------------------
echo "[2/7] Downloading EMAG2 v3 magnetic anomaly grid (~1.1 GB compressed)..."
mkdir -p data/geophysical/emag2
if [ ! -f data/geophysical/emag2/EMAG2_V3_20170530.zip ]; then
    curl -L -o data/geophysical/emag2/EMAG2_V3_20170530.zip \
        "https://www.ngdc.noaa.gov/geomag/EMAG2/data/EMAG2_V3_20170530.zip"
    echo "  NOTE: If URL has changed, search for 'EMAG2 v3' at https://www.ncei.noaa.gov/"
else
    echo "  Already exists, skipping."
fi

# ----------------------------------------------------------------
# HYDE 3.2 Historical Population Density (Paper 2: population baseline)
# Source: Netherlands Environmental Assessment Agency (PBL)
# ----------------------------------------------------------------
echo "[3/7] Downloading HYDE 3.2 population density grids (~410 MB)..."
mkdir -p data/hyde
echo "  HYDE data requires manual download from:"
echo "  https://doi.org/10.17026/dans-25g-gez3"
echo "  Download the population density (popd) grids for: 3000BC, 2000BC, 1000BC, 0AD, 1000AD"
echo "  Place .asc files in data/hyde/"

# ----------------------------------------------------------------
# SRTM elevation tiles (Paper 4: cost surface generation)
# Source: NASA Shuttle Radar Topography Mission via OpenTopography
# ----------------------------------------------------------------
echo "[4/7] SRTM elevation data..."
mkdir -p data/srtm_cache
echo "  SRTM tiles are downloaded on demand by the analysis scripts."
echo "  They will be cached in data/srtm_cache/ (~1 GB for full coverage)."
echo "  No manual download needed — scripts use the 'elevation' or 'srtmgl1' API."

# ----------------------------------------------------------------
# Groundwater data (Paper 3: causal mechanisms)
# Source: Fan et al. (2013) Science
# ----------------------------------------------------------------
echo "[5/7] Groundwater water table depth data..."
mkdir -p data/groundwater
echo "  Download from Fan et al. (2013) supplementary data:"
echo "  https://doi.org/10.1126/science.1229881"
echo "  Files needed: EURASIA_WTD_annualmean.nc, AFRICA_WTD_annualmean.nc"
echo "  Place in data/groundwater/"

# ----------------------------------------------------------------
# XRONOS radiocarbon database (Paper 1-2: independent replication)
# Source: xronos.ch
# ----------------------------------------------------------------
echo "[6/7] Downloading XRONOS radiocarbon database..."
mkdir -p data/xronos
if [ ! -f data/xronos/xronos_raw.json ]; then
    echo "  XRONOS data is accessed via API at https://xronos.ch"
    echo "  Run: curl -o data/xronos/xronos_raw.json 'https://xronos.ch/api/data?format=json'"
    echo "  Note: Full download is ~168 MB"
else
    echo "  Already exists, skipping."
fi

# ----------------------------------------------------------------
# JPL DE441 Ephemeris (Orion paper: precession calculations)
# Source: NASA JPL
# ----------------------------------------------------------------
echo "[7/7] JPL DE441 planetary ephemeris (~3.1 GB)..."
echo "  Only needed for the Jelitto planetary test (exploratory, not in any paper)."
echo "  If needed, Skyfield will download it automatically on first use."
echo "  Or manually: curl -o data/de441.bsp https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441.bsp"

echo ""
echo "=== Download complete ==="
echo "Some files require manual download (HYDE, groundwater). See instructions above."
echo ""
echo "Verify data integrity:"
echo "  python3 -c \"import os; print(sum(os.path.getsize(os.path.join(dp,f)) for dp,dn,fn in os.walk('data') for f in fn) / 1e9, 'GB')\""
