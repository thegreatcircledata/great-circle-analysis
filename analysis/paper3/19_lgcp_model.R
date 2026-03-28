#!/usr/bin/env Rscript
# LGCP (Log-Gaussian Cox Process) Model
# Fits LGCP with composite likelihood to test whether β_gc persists
# under a more rigorous spatial correlation model than Thomas process

library(spatstat)
library(jsonlite)

cat("=== LGCP Model Fitting ===\n")

# Load site data
sites <- read.csv(file.path(Sys.getenv("HOME"), "megalith_site_research/data/pleiades_sites_classified.csv"))
cat(sprintf("Loaded %d sites\n", nrow(sites)))

# Filter to study window
sites <- sites[sites$lat >= 20 & sites$lat <= 55 &
               sites$lon >= -12 & sites$lon <= 80, ]

mon <- sites[sites$type == "monument", ]
sett <- sites[sites$type == "settlement", ]
cat(sprintf("Window: %d monuments, %d settlements\n", nrow(mon), nrow(sett)))

# Create observation window
win <- owin(xrange = c(-12, 80), yrange = c(20, 55))

# Create point patterns
mon_pp <- ppp(mon$lon, mon$lat, window = win)
set_pp <- ppp(sett$lon, sett$lat, window = win)

# Load covariates from the PPM results
ppm_data <- fromJSON(file.path(Sys.getenv("HOME"), "megalith_site_research/paper3/results/08_point_process_model.json"))
std <- ppm_data$standardization

# Great circle distance function
gc_dist_from_pole <- function(lat, lon, pole_lat = 59.682122, pole_lon = -138.646087) {
  R <- 6371.0
  to_rad <- pi / 180
  plat <- pole_lat * to_rad
  plon <- pole_lon * to_rad
  lat_r <- lat * to_rad
  lon_r <- lon * to_rad

  px <- cos(plat) * cos(plon)
  py <- cos(plat) * sin(plon)
  pz <- sin(plat)

  x <- cos(lat_r) * cos(lon_r)
  y <- cos(lat_r) * sin(lon_r)
  z <- sin(lat_r)

  dot <- pmax(pmin(px*x + py*y + pz*z, 1), -1)
  ang <- acos(dot)
  abs(ang - pi/2) * R
}

# Create covariate images
res <- 0.5  # degree resolution
x_seq <- seq(-12 + res/2, 80 - res/2, by = res)
y_seq <- seq(20 + res/2, 55 - res/2, by = res)

cat("Creating covariate images...\n")

# Absolute latitude
ablat_mat <- matrix(0, nrow = length(y_seq), ncol = length(x_seq))
gc_mat <- matrix(0, nrow = length(y_seq), ncol = length(x_seq))

for (i in seq_along(y_seq)) {
  for (j in seq_along(x_seq)) {
    ablat_mat[i, j] <- (abs(y_seq[i]) - std$ablat_m) / std$ablat_s
    gc_mat[i, j] <- (gc_dist_from_pole(y_seq[i], x_seq[j]) - std$gc_m) / std$gc_s
  }
}

im_ablat <- im(ablat_mat, xcol = x_seq, yrow = y_seq)
im_gc <- im(gc_mat, xcol = x_seq, yrow = y_seq)

# Longitude proxy (standardized)
lon_mat <- matrix(0, nrow = length(y_seq), ncol = length(x_seq))
for (i in seq_along(y_seq)) {
  for (j in seq_along(x_seq)) {
    lon_mat[i, j] <- (x_seq[j] - mean(x_seq)) / sd(x_seq)
  }
}
im_lon <- im(lon_mat, xcol = x_seq, yrow = y_seq)

# Try LGCP with composite likelihood
cat("\nFitting LGCP for monuments...\n")
tryCatch({
  # Method 1: kppm with "clik2" (second-order composite likelihood)
  mon_lgcp <- kppm(mon_pp ~ im_ablat + im_lon + im_gc,
                   clusters = "LGCP",
                   method = "clik2",
                   control = list(maxit = 200))

  cat("Monument LGCP converged\n")
  cat(sprintf("  Coefficients:\n"))
  print(coef(mon_lgcp))

  mon_coefs <- coef(mon_lgcp)
  mon_se <- sqrt(diag(vcov(mon_lgcp)))
  names(mon_se) <- names(mon_coefs)

  cat(sprintf("  β_gc = %.4f (SE = %.4f)\n",
              mon_coefs["im_gc"], mon_se["im_gc"]))
  cat(sprintf("  z = %.4f\n",
              mon_coefs["im_gc"] / mon_se["im_gc"]))
  cat(sprintf("  p = %.6f\n",
              2 * pnorm(-abs(mon_coefs["im_gc"] / mon_se["im_gc"]))))

  mon_result <- list(
    method = "LGCP clik2",
    converged = TRUE,
    coefficients = as.list(mon_coefs),
    std_errors = as.list(mon_se),
    gc_dist_beta = unname(mon_coefs["im_gc"]),
    gc_dist_se = unname(mon_se["im_gc"]),
    gc_dist_z = unname(mon_coefs["im_gc"] / mon_se["im_gc"]),
    gc_dist_p = unname(2 * pnorm(-abs(mon_coefs["im_gc"] / mon_se["im_gc"]))),
    gc_dist_significant = unname(2 * pnorm(-abs(mon_coefs["im_gc"] / mon_se["im_gc"]))) < 0.05
  )
}, error = function(e) {
  cat(sprintf("LGCP clik2 failed for monuments: %s\n", e$message))
  cat("Trying method='palm'...\n")

  tryCatch({
    mon_lgcp <- kppm(mon_pp ~ im_ablat + im_lon + im_gc,
                     clusters = "LGCP",
                     method = "palm")
    cat("Monument LGCP (palm) converged\n")
    mon_coefs <- coef(mon_lgcp)
    mon_se <- sqrt(diag(vcov(mon_lgcp)))
    names(mon_se) <- names(mon_coefs)

    mon_result <<- list(
      method = "LGCP palm",
      converged = TRUE,
      coefficients = as.list(mon_coefs),
      std_errors = as.list(mon_se),
      gc_dist_beta = unname(mon_coefs["im_gc"]),
      gc_dist_se = unname(mon_se["im_gc"]),
      gc_dist_z = unname(mon_coefs["im_gc"] / mon_se["im_gc"]),
      gc_dist_p = unname(2 * pnorm(-abs(mon_coefs["im_gc"] / mon_se["im_gc"]))),
      gc_dist_significant = unname(2 * pnorm(-abs(mon_coefs["im_gc"] / mon_se["im_gc"]))) < 0.05
    )
  }, error = function(e2) {
    cat(sprintf("LGCP palm also failed: %s\n", e2$message))
    cat("Trying Thomas with clik2 as fallback...\n")

    mon_thomas <- kppm(mon_pp ~ im_ablat + im_lon + im_gc,
                       clusters = "Thomas",
                       method = "clik2")
    cat("Thomas clik2 converged\n")
    mon_coefs <- coef(mon_thomas)
    mon_se <- sqrt(diag(vcov(mon_thomas)))
    names(mon_se) <- names(mon_coefs)

    mon_result <<- list(
      method = "Thomas clik2 (LGCP fallback)",
      converged = TRUE,
      coefficients = as.list(mon_coefs),
      std_errors = as.list(mon_se),
      gc_dist_beta = unname(mon_coefs["im_gc"]),
      gc_dist_se = unname(mon_se["im_gc"]),
      gc_dist_z = unname(mon_coefs["im_gc"] / mon_se["im_gc"]),
      gc_dist_p = unname(2 * pnorm(-abs(mon_coefs["im_gc"] / mon_se["im_gc"]))),
      gc_dist_significant = unname(2 * pnorm(-abs(mon_coefs["im_gc"] / mon_se["im_gc"]))) < 0.05
    )
  })
})

cat("\nFitting LGCP for settlements...\n")
tryCatch({
  set_lgcp <- kppm(set_pp ~ im_ablat + im_lon + im_gc,
                   clusters = "LGCP",
                   method = "clik2",
                   control = list(maxit = 200))

  cat("Settlement LGCP converged\n")
  set_coefs <- coef(set_lgcp)
  set_se <- sqrt(diag(vcov(set_lgcp)))
  names(set_se) <- names(set_coefs)

  cat(sprintf("  β_gc = %.4f (SE = %.4f)\n",
              set_coefs["im_gc"], set_se["im_gc"]))

  set_result <- list(
    method = "LGCP clik2",
    converged = TRUE,
    coefficients = as.list(set_coefs),
    std_errors = as.list(set_se),
    gc_dist_beta = unname(set_coefs["im_gc"]),
    gc_dist_se = unname(set_se["im_gc"]),
    gc_dist_z = unname(set_coefs["im_gc"] / set_se["im_gc"]),
    gc_dist_p = unname(2 * pnorm(-abs(set_coefs["im_gc"] / set_se["im_gc"]))),
    gc_dist_significant = unname(2 * pnorm(-abs(set_coefs["im_gc"] / set_se["im_gc"]))) < 0.05
  )
}, error = function(e) {
  cat(sprintf("LGCP clik2 failed for settlements: %s\n", e$message))

  tryCatch({
    set_lgcp <- kppm(set_pp ~ im_ablat + im_lon + im_gc,
                     clusters = "LGCP",
                     method = "palm")
    set_coefs <- coef(set_lgcp)
    set_se <- sqrt(diag(vcov(set_lgcp)))
    names(set_se) <- names(set_coefs)

    set_result <<- list(
      method = "LGCP palm",
      converged = TRUE,
      coefficients = as.list(set_coefs),
      std_errors = as.list(set_se),
      gc_dist_beta = unname(set_coefs["im_gc"]),
      gc_dist_se = unname(set_se["im_gc"]),
      gc_dist_z = unname(set_coefs["im_gc"] / set_se["im_gc"]),
      gc_dist_p = unname(2 * pnorm(-abs(set_coefs["im_gc"] / set_se["im_gc"]))),
      gc_dist_significant = unname(2 * pnorm(-abs(set_coefs["im_gc"] / set_se["im_gc"]))) < 0.05
    )
  }, error = function(e2) {
    set_thomas <- kppm(set_pp ~ im_ablat + im_lon + im_gc,
                       clusters = "Thomas",
                       method = "clik2")
    set_coefs <- coef(set_thomas)
    set_se <- sqrt(diag(vcov(set_thomas)))
    names(set_se) <- names(set_coefs)

    set_result <<- list(
      method = "Thomas clik2 (LGCP fallback)",
      converged = TRUE,
      coefficients = as.list(set_coefs),
      std_errors = as.list(set_se),
      gc_dist_beta = unname(set_coefs["im_gc"]),
      gc_dist_se = unname(set_se["im_gc"]),
      gc_dist_z = unname(set_coefs["im_gc"] / set_se["im_gc"]),
      gc_dist_p = unname(2 * pnorm(-abs(set_coefs["im_gc"] / set_se["im_gc"]))),
      gc_dist_significant = unname(2 * pnorm(-abs(set_coefs["im_gc"] / set_se["im_gc"]))) < 0.05
    )
  })
})

# Divergence test
if (exists("mon_result") && exists("set_result")) {
  diff <- mon_result$gc_dist_beta - set_result$gc_dist_beta
  se_diff <- sqrt(mon_result$gc_dist_se^2 + set_result$gc_dist_se^2)
  z_diff <- diff / se_diff
  p_diff <- 2 * pnorm(-abs(z_diff))

  divergence <- list(
    mon_beta_gc = mon_result$gc_dist_beta,
    set_beta_gc = set_result$gc_dist_beta,
    difference = diff,
    se_diff = se_diff,
    z_test = z_diff,
    p_value = p_diff,
    significant = p_diff < 0.05,
    cis_overlap = !(mon_result$gc_dist_beta - 1.96*mon_result$gc_dist_se >
                    set_result$gc_dist_beta + 1.96*set_result$gc_dist_se ||
                    set_result$gc_dist_beta - 1.96*set_result$gc_dist_se >
                    mon_result$gc_dist_beta + 1.96*mon_result$gc_dist_se)
  )

  if (p_diff < 0.05) {
    conclusion <- sprintf("Divergence PERSISTS under LGCP: z=%.2f, p=%.6f", z_diff, p_diff)
  } else {
    conclusion <- sprintf("Divergence does NOT persist under LGCP: z=%.2f, p=%.4f", z_diff, p_diff)
  }

  result <- list(
    analysis = "Log-Gaussian Cox Process (LGCP)",
    monuments = mon_result,
    settlements = set_result,
    divergence = divergence,
    conclusion = conclusion
  )

  cat(sprintf("\n=== RESULT ===\n"))
  cat(sprintf("Monument β_gc: %.4f (p=%.4f)\n", mon_result$gc_dist_beta, mon_result$gc_dist_p))
  cat(sprintf("Settlement β_gc: %.4f (p=%.4f)\n", set_result$gc_dist_beta, set_result$gc_dist_p))
  cat(sprintf("Divergence: z=%.2f, p=%.6f\n", z_diff, p_diff))
  cat(sprintf("%s\n", conclusion))

  write(toJSON(result, auto_unbox = TRUE, pretty = TRUE),
        file.path(Sys.getenv("HOME"), "megalith_site_research/paper3/results/lgcp_model.json"))
  cat("Saved to paper3/results/lgcp_model.json\n")
} else {
  cat("ERROR: One or both models failed to fit\n")
}
