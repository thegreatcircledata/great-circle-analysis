#!/usr/bin/env Rscript
# Paper 3 — kppm() Cluster Process Model
# ========================================
# Fit Thomas cluster process via kppm() with method="mincon" on the
# pair correlation function, following Bilotti et al. (2024) methodology.
# Tests whether GC distance remains significant after accounting for
# spatial clustering.

library(spatstat)
library(jsonlite)

cat("=" , rep("=", 69), "\n", sep="")
cat("KPPM CLUSTER PROCESS MODEL\n")
cat("=" , rep("=", 69), "\n", sep="")

BASE <- path.expand("~/megalith_site_research")
OUT <- file.path(BASE, "paper3/results")

# ============================================================
# LOAD DATA
# ============================================================
cat("\nLoading Pleiades data...\n")

pleiades <- read.csv(file.path(BASE, "data/pleiades/pleiades-places-latest.csv"),
                     stringsAsFactors = FALSE)

MONUMENTAL <- c("pyramid", "temple", "temple-2", "sanctuary", "monument",
                "architecturalcomplex", "tomb", "church", "church-2",
                "amphitheatre", "hippodrome", "aqueduct", "dam", "lighthouse",
                "arch", "theatre", "basilica", "acropolis", "nuraghe", "tumulus",
                "shrine", "cemetery")
SETTLEMENT <- c("settlement", "settlement-modern", "village", "villa", "farmstead",
                "station", "findspot", "port", "mine", "mine-2", "quarry",
                "bath", "cistern", "bridge", "fortified-settlement", "townhouse", "production")

# Parse and classify
pleiades$lat <- as.numeric(pleiades$reprLat)
pleiades$lon <- as.numeric(pleiades$reprLong)
pleiades <- pleiades[!is.na(pleiades$lat) & !is.na(pleiades$lon), ]
pleiades <- pleiades[abs(pleiades$lat) > 0.001 | abs(pleiades$lon) > 0.001, ]

classify <- function(ft_string) {
  fts <- trimws(strsplit(ft_string, ",")[[1]])
  is_m <- any(fts %in% MONUMENTAL)
  is_s <- any(fts %in% SETTLEMENT)
  if (is_m & !is_s) return("monument")
  if (is_s & !is_m) return("settlement")
  return("other")
}

pleiades$type <- sapply(pleiades$featureTypes, classify)
monuments <- pleiades[pleiades$type == "monument", ]
settlements <- pleiades[pleiades$type == "settlement", ]

cat(sprintf("  Monuments: %d, Settlements: %d\n", nrow(monuments), nrow(settlements)))

# ============================================================
# COMPUTE GC DISTANCE
# ============================================================
POLE_LAT <- 59.682122
POLE_LON <- -138.646087
R_EARTH <- 6371.0
QC <- R_EARTH * pi / 2

gc_distance <- function(lats, lons) {
  la1 <- POLE_LAT * pi / 180
  lo1 <- POLE_LON * pi / 180
  la2 <- lats * pi / 180
  lo2 <- lons * pi / 180
  dl <- la2 - la1
  dn <- lo2 - lo1
  a <- sin(dl/2)^2 + cos(la1) * cos(la2) * sin(dn/2)^2
  d <- R_EARTH * 2 * asin(sqrt(pmin(1.0, a)))
  abs(d - QC)
}

monuments$gc_dist <- gc_distance(monuments$lat, monuments$lon)
settlements$gc_dist <- gc_distance(settlements$lat, settlements$lon)

# ============================================================
# RESTRICT TO MEDITERRANEAN/MIDDLE EAST WINDOW
# ============================================================
# Pleiades is concentrated in lat 20-55, lon -12 to 80
# Use a rectangular window in projected coordinates (km from origin)
cat("\nSetting up spatial window...\n")

lat_min <- 20; lat_max <- 55
lon_min <- -12; lon_max <- 80

# Filter sites to window
monuments <- monuments[monuments$lat >= lat_min & monuments$lat <= lat_max &
                       monuments$lon >= lon_min & monuments$lon <= lon_max, ]
settlements <- settlements[settlements$lat >= lat_min & settlements$lat <= lat_max &
                           settlements$lon >= lon_min & settlements$lon <= lon_max, ]

cat(sprintf("  After window filter - Monuments: %d, Settlements: %d\n",
            nrow(monuments), nrow(settlements)))

# Project to km (simple equirectangular from center)
center_lat <- (lat_min + lat_max) / 2
km_per_deg_lat <- 111.32
km_per_deg_lon <- 111.32 * cos(center_lat * pi / 180)

monuments$x <- (monuments$lon - lon_min) * km_per_deg_lon
monuments$y <- (monuments$lat - lat_min) * km_per_deg_lat
settlements$x <- (settlements$lon - lon_min) * km_per_deg_lon
settlements$y <- (settlements$lat - lat_min) * km_per_deg_lat

# Create window
win <- owin(xrange = c(0, (lon_max - lon_min) * km_per_deg_lon),
            yrange = c(0, (lat_max - lat_min) * km_per_deg_lat))

# Create point patterns
pp_mon <- ppp(monuments$x, monuments$y, window = win)
pp_set <- ppp(settlements$x, settlements$y, window = win)

cat(sprintf("  Window: %.0f x %.0f km\n", diff(win$xrange), diff(win$yrange)))
cat(sprintf("  Monument point pattern: n=%d\n", pp_mon$n))
cat(sprintf("  Settlement point pattern: n=%d\n", pp_set$n))

# ============================================================
# ATTACH COVARIATES AS PIXEL IMAGES
# ============================================================
cat("\nCreating covariate images...\n")

# Create covariate functions that can be evaluated at any point
# We'll create pixel images from the site data

# For kppm with covariates, we need im objects
# Create them at ~50 km resolution (sufficient for this scale)
nx <- 50
ny <- 25

# Absolute latitude image
make_ablat_im <- function() {
  mat <- matrix(0, ny, nx)
  x_seq <- seq(win$xrange[1], win$xrange[2], length.out = nx)
  y_seq <- seq(win$yrange[1], win$yrange[2], length.out = ny)
  for (i in 1:ny) {
    lat_val <- lat_min + y_seq[i] / km_per_deg_lat
    mat[i, ] <- abs(lat_val)
  }
  im(mat, xcol = x_seq, yrow = y_seq)
}

# GC distance image
make_gc_dist_im <- function() {
  mat <- matrix(0, ny, nx)
  x_seq <- seq(win$xrange[1], win$xrange[2], length.out = nx)
  y_seq <- seq(win$yrange[1], win$yrange[2], length.out = ny)
  for (i in 1:ny) {
    for (j in 1:nx) {
      lat_val <- lat_min + y_seq[i] / km_per_deg_lat
      lon_val <- lon_min + x_seq[j] / km_per_deg_lon
      mat[i, j] <- gc_distance(lat_val, lon_val)
    }
  }
  im(mat, xcol = x_seq, yrow = y_seq)
}

# Coast distance approximation (use longitude as proxy — sites closer to
# Mediterranean coast are at lower longitude in this window)
# This is a simplification; proper coast distance would need coastline data
make_lon_proxy_im <- function() {
  mat <- matrix(0, ny, nx)
  x_seq <- seq(win$xrange[1], win$xrange[2], length.out = nx)
  y_seq <- seq(win$yrange[1], win$yrange[2], length.out = ny)
  for (i in 1:ny) {
    for (j in 1:nx) {
      lon_val <- lon_min + x_seq[j] / km_per_deg_lon
      mat[i, j] <- lon_val  # raw longitude as proxy
    }
  }
  im(mat, xcol = x_seq, yrow = y_seq)
}

im_ablat <- make_ablat_im()
im_gc_dist <- make_gc_dist_im()
im_lon <- make_lon_proxy_im()

cat("  Covariate images created (abs_lat, gc_dist, lon_proxy)\n")

# ============================================================
# FIT kppm() MODELS
# ============================================================
cat("\n", rep("=", 70), "\n", sep="")
cat("FITTING kppm() MODELS\n")
cat(rep("=", 70), "\n", sep="")

# Model 1: Thomas process WITHOUT GC distance
cat("\n--- Monument model WITHOUT GC distance ---\n")
tryCatch({
  kppm_mon_no_gc <- kppm(pp_mon ~ im_ablat + im_lon,
                          clusters = "Thomas",
                          method = "mincon",
                          statistic = "pcf")
  cat("  Converged!\n")
  print(summary(kppm_mon_no_gc))
  mon_no_gc_ok <- TRUE
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message))
  mon_no_gc_ok <<- FALSE
})

# Model 2: Thomas process WITH GC distance
cat("\n--- Monument model WITH GC distance ---\n")
tryCatch({
  kppm_mon_gc <- kppm(pp_mon ~ im_ablat + im_lon + im_gc_dist,
                       clusters = "Thomas",
                       method = "mincon",
                       statistic = "pcf")
  cat("  Converged!\n")
  print(summary(kppm_mon_gc))
  mon_gc_ok <- TRUE
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message))
  mon_gc_ok <<- FALSE
})

# Model 3: Settlement WITHOUT GC distance
cat("\n--- Settlement model WITHOUT GC distance ---\n")
tryCatch({
  kppm_set_no_gc <- kppm(pp_set ~ im_ablat + im_lon,
                          clusters = "Thomas",
                          method = "mincon",
                          statistic = "pcf")
  cat("  Converged!\n")
  print(summary(kppm_set_no_gc))
  set_no_gc_ok <- TRUE
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message))
  set_no_gc_ok <<- FALSE
})

# Model 4: Settlement WITH GC distance
cat("\n--- Settlement model WITH GC distance ---\n")
tryCatch({
  kppm_set_gc <- kppm(pp_set ~ im_ablat + im_lon + im_gc_dist,
                       clusters = "Thomas",
                       method = "mincon",
                       statistic = "pcf")
  cat("  Converged!\n")
  print(summary(kppm_set_gc))
  set_gc_ok <- TRUE
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message))
  set_gc_ok <<- FALSE
})

# ============================================================
# EXTRACT COEFFICIENTS AND BUILD RESULTS
# ============================================================
cat("\n", rep("=", 70), "\n", sep="")
cat("RESULTS SUMMARY\n")
cat(rep("=", 70), "\n", sep="")

results <- list(
  analysis = "kppm Thomas Cluster Process (spatstat)",
  method = "mincon on pair correlation function",
  window = list(lat_range = c(lat_min, lat_max),
                lon_range = c(lon_min, lon_max)),
  covariates = c("abs_lat", "lon_proxy", "gc_dist"),
  note = "lon_proxy used instead of coast_dist (coastline shapefile not loaded in R)"
)

if (exists("mon_gc_ok") && mon_gc_ok) {
  coefs_mon <- coef(kppm_mon_gc)
  se_mon <- sqrt(diag(vcov(kppm_mon_gc)))
  clust_mon <- kppm_mon_gc$clustpar

  cat("\nMonument model (with GC):\n")
  cat(sprintf("  GC distance coefficient: %.6f\n", coefs_mon["im_gc_dist"]))
  cat(sprintf("  GC distance SE: %.6f\n", se_mon["im_gc_dist"]))
  cat(sprintf("  GC distance z: %.4f\n", coefs_mon["im_gc_dist"] / se_mon["im_gc_dist"]))
  cat(sprintf("  Cluster scale (kappa): %.4f\n", clust_mon["kappa"]))
  cat(sprintf("  Cluster radius (scale): %.4f\n", clust_mon["scale"]))

  z_mon <- coefs_mon["im_gc_dist"] / se_mon["im_gc_dist"]
  p_mon <- 2 * pnorm(-abs(z_mon))

  results$monuments <- list(
    coefficients = as.list(coefs_mon),
    std_errors = as.list(se_mon),
    gc_dist_beta = unname(coefs_mon["im_gc_dist"]),
    gc_dist_se = unname(se_mon["im_gc_dist"]),
    gc_dist_z = unname(z_mon),
    gc_dist_p = unname(p_mon),
    gc_dist_significant = unname(p_mon < 0.05),
    cluster_params = as.list(clust_mon)
  )
} else {
  cat("\nMonument model with GC: FAILED\n")
  results$monuments <- list(status = "failed")
}

if (exists("set_gc_ok") && set_gc_ok) {
  coefs_set <- coef(kppm_set_gc)
  se_set <- sqrt(diag(vcov(kppm_set_gc)))
  clust_set <- kppm_set_gc$clustpar

  cat("\nSettlement model (with GC):\n")
  cat(sprintf("  GC distance coefficient: %.6f\n", coefs_set["im_gc_dist"]))
  cat(sprintf("  GC distance SE: %.6f\n", se_set["im_gc_dist"]))
  cat(sprintf("  GC distance z: %.4f\n", coefs_set["im_gc_dist"] / se_set["im_gc_dist"]))
  cat(sprintf("  Cluster scale (kappa): %.4f\n", clust_set["kappa"]))
  cat(sprintf("  Cluster radius (scale): %.4f\n", clust_set["scale"]))

  z_set <- coefs_set["im_gc_dist"] / se_set["im_gc_dist"]
  p_set <- 2 * pnorm(-abs(z_set))

  results$settlements <- list(
    coefficients = as.list(coefs_set),
    std_errors = as.list(se_set),
    gc_dist_beta = unname(coefs_set["im_gc_dist"]),
    gc_dist_se = unname(se_set["im_gc_dist"]),
    gc_dist_z = unname(z_set),
    gc_dist_p = unname(p_set),
    gc_dist_significant = unname(p_set < 0.05),
    cluster_params = as.list(clust_set)
  )
} else {
  cat("\nSettlement model with GC: FAILED\n")
  results$settlements <- list(status = "failed")
}

# Divergence comparison
if (exists("mon_gc_ok") && mon_gc_ok && exists("set_gc_ok") && set_gc_ok) {
  beta_diff <- unname(coefs_mon["im_gc_dist"] - coefs_set["im_gc_dist"])
  se_diff <- sqrt(unname(se_mon["im_gc_dist"])^2 + unname(se_set["im_gc_dist"])^2)
  z_diff <- beta_diff / se_diff
  p_diff <- 2 * pnorm(-abs(z_diff))

  cat(sprintf("\nDIVERGENCE TEST (kppm):\n"))
  cat(sprintf("  Monument β_gc: %.6f (SE=%.6f)\n", coefs_mon["im_gc_dist"], se_mon["im_gc_dist"]))
  cat(sprintf("  Settlement β_gc: %.6f (SE=%.6f)\n", coefs_set["im_gc_dist"], se_set["im_gc_dist"]))
  cat(sprintf("  Difference: %.6f\n", beta_diff))
  cat(sprintf("  z-test: %.4f\n", z_diff))
  cat(sprintf("  p-value: %.2e\n", p_diff))
  cat(sprintf("  CIs overlap: %s\n",
              ifelse(abs(beta_diff) < 1.96 * se_diff, "YES", "NO")))

  results$divergence <- list(
    mon_beta_gc = unname(coefs_mon["im_gc_dist"]),
    set_beta_gc = unname(coefs_set["im_gc_dist"]),
    difference = beta_diff,
    se_diff = se_diff,
    z_test = z_diff,
    p_value = p_diff,
    significant = p_diff < 0.05,
    cis_overlap = abs(beta_diff) < 1.96 * se_diff
  )

  if (p_diff < 0.05) {
    cat("\n→ DIVERGENCE PERSISTS under Thomas cluster process model\n")
    results$conclusion <- "Divergence persists under kppm Thomas cluster process"
  } else {
    cat("\n→ Divergence does NOT persist under Thomas cluster process model\n")
    results$conclusion <- "Divergence does NOT persist under kppm Thomas cluster process"
  }
}

# ============================================================
# SAVE
# ============================================================
write_json(results, file.path(OUT, "kppm_cluster_process.json"),
           pretty = TRUE, auto_unbox = TRUE)
cat(sprintf("\nResults saved: %s\n", file.path(OUT, "kppm_cluster_process.json")))
cat("Done.\n")
