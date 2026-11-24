# ---------------------------------------------------------------------------
# ATSF2026 SIRS + seasonal-beta + EnKF forecasting driver
# ---------------------------------------------------------------------------
# Model setup:
#   - 11 locations, each forecast INDEPENDENTLY (no spatial coupling)
#   - Each location has its own effective population N_eff
#   - Plain SIRS compartments (S, I, R) with WANING immunity (R -> S at rate omega)
#   - Transmission rate beta evolves via OU on log-scale with a
#     SEASONAL mean (winter-peaked), not a constant
#   - EnKF joint state-parameter update on [S, I, R, beta]
#   - Observation: ILI% mapped to I compartment via I_count = ili% * N_eff / 100
#
# Output:
#   ONE CSV per origin_date (134 files) in output_dir.
#   Each file has 1012 rows = 11 locations x 4 horizons x 23 quantiles.
#   Filename pattern: {origin_date}-syed-SIRSEnKF.csv
#
# Usage (from the repository root):  Rscript R/run_forecast.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(lubridate)
})

source("R/sirs_tau_leap.R")
source("R/forecast_step.R")
source("R/enkf_update.R")

# ---------------- Configuration ----------------------------------------------

CONFIG <- list(
  N_ensemble    = 300,
  steps_per_obs = 7,
  horizons      = 1:4,
  
  N_eff_default = 1e5,
  N_eff_by_loc  = NULL,
  
  # SIRS dynamics
  gamma = 1/7,
  omega = 1/180,            # ~6-month average immunity
  
  # Seasonal OU on log(beta)
  ou_theta  = 1/200,
  beta0     = 0.25,
  a_season  = 0.4,
  doy_peak  = 25,           # Jan 25 (historical US flu peak)
  ou_sigma  = 0.06,
  
  obs_sigma_frac = 5e-4,
  
  # Initial ensemble.
  init_I_frac_range  = c(1e-4, 5e-3),
  init_R_frac_range  = c(0.20, 0.40),
  init_beta_range    = c(0.15, 0.40),
  
  quantile_levels = c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99),
  
  # Output settings
  output_dir   = "model-output/syed-SIRSEnKF",
  model_name   = "syed-SIRSEnKF",
  overwrite    = TRUE
)

# ---------------- Helpers ----------------------------------------------------

ili_to_I_count <- function(ili_perc, pop) ili_perc / 100 * pop
I_count_to_ili <- function(I_count, pop)  I_count / pop * 100

N_eff_for <- function(loc, config) {
  if (!is.null(config$N_eff_by_loc) && loc %in% names(config$N_eff_by_loc)) {
    return(config$N_eff_by_loc[[loc]])
  }
  config$N_eff_default
}

season_start_for <- function(ref_date, all_dates) {
  yr <- year(ref_date)
  season_yr <- if (month(ref_date) >= 8) yr else yr - 1
  season_cutoff <- as.Date(sprintf("%d-08-01", season_yr))
  candidates <- all_dates[all_dates >= season_cutoff & all_dates <= ref_date]
  if (length(candidates) == 0) return(NA)
  min(candidates)
}

init_ensemble <- function(N, pop, config) {
  I_frac <- runif(N, config$init_I_frac_range[1], config$init_I_frac_range[2])
  R_frac <- runif(N, config$init_R_frac_range[1], config$init_R_frac_range[2])
  
  I_vec <- I_frac * pop
  R_vec <- R_frac * pop
  S_vec <- pop - I_vec - R_vec
  
  sir <- array(0, dim = c(1, 3, N))
  sir[1, 1, ] <- S_vec
  sir[1, 2, ] <- I_vec
  sir[1, 3, ] <- R_vec
  
  beta <- matrix(runif(N, config$init_beta_range[1], config$init_beta_range[2]),
                 nrow = 1, ncol = N)
  list(sir = sir, beta = beta)
}

# ---------------- Core filter + forecast routines ---------------------------

run_filter_one_location <- function(obs_series, obs_dates, pop, config) {
  N <- config$N_ensemble
  ens <- init_ensemble(N = N, pop = pop, config = config)
  obs_sigma <- config$obs_sigma_frac * pop
  
  cur_date <- obs_dates[1] - config$steps_per_obs
  
  for (k in seq_along(obs_series)) {
    fc <- forecast_step(ens$sir, ens$beta,
                        population    = pop,
                        current_date  = cur_date,
                        steps_per_obs = config$steps_per_obs,
                        theta    = config$ou_theta,
                        beta0    = config$beta0,
                        a_season = config$a_season,
                        doy_peak = config$doy_peak,
                        sigma    = config$ou_sigma,
                        gamma = config$gamma,
                        omega = config$omega)
    ens$sir  <- fc$sir
    ens$beta <- fc$beta
    
    if (!is.na(obs_series[k])) {
      upd <- enkf_update_step(ens$sir, ens$beta,
                              observation = obs_series[k],
                              fixed_sigma = obs_sigma)
      ens$sir  <- upd$sir
      ens$beta <- upd$beta
    }
    
    cur_date <- obs_dates[k]
  }
  list(ens = ens, last_date = cur_date)
}

project_forward <- function(ens, last_date, pop, horizons, config) {
  N <- dim(ens$beta)[2]
  out <- matrix(0, nrow = length(horizons), ncol = N)
  max_h <- max(horizons)
  
  cur      <- ens
  cur_date <- last_date
  for (h in seq_len(max_h)) {
    fc <- forecast_step(cur$sir, cur$beta,
                        population    = pop,
                        current_date  = cur_date,
                        steps_per_obs = config$steps_per_obs,
                        theta    = config$ou_theta,
                        beta0    = config$beta0,
                        a_season = config$a_season,
                        doy_peak = config$doy_peak,
                        sigma    = config$ou_sigma,
                        gamma = config$gamma,
                        omega = config$omega)
    cur      <- fc
    cur_date <- cur_date + config$steps_per_obs
    if (h %in% horizons) {
      out[which(horizons == h), ] <- fc$sir[1, 2, ]     # I is row 2
    }
  }
  out
}

# ---------------- Per-date forecast (writes one CSV) ------------------------

forecast_one_date <- function(ori_date, raw, all_dates, locs,
                              config, verbose = TRUE) {
  ori_date <- as.Date(ori_date)
  season_start <- season_start_for(ori_date, all_dates)
  if (is.na(season_start)) {
    message("  skip: no observations before ", ori_date)
    return(invisible(NULL))
  }
  
  rows <- list()
  idx  <- 1L
  
  for (loc in locs) {
    pop <- N_eff_for(loc, config)
    
    hist <- raw |>
      filter(location == loc,
             target_end_date >= season_start,
             target_end_date <= ori_date) |>
      arrange(target_end_date)
    if (nrow(hist) < 2) {
      if (verbose) message("  skip location ", loc, " (insufficient history)")
      next
    }
    
    obs_counts <- ili_to_I_count(hist$observation, pop)
    obs_dates  <- hist$target_end_date
    
    filter_out <- run_filter_one_location(obs_counts, obs_dates, pop, config)
    ens        <- filter_out$ens
    last_date  <- filter_out$last_date
    
    fc      <- project_forward(ens, last_date, pop, config$horizons, config)
    fc_perc <- I_count_to_ili(fc, pop)
    
    for (hi in seq_along(config$horizons)) {
      h     <- config$horizons[hi]
      tgt_d <- ori_date + 7 * h
      qs    <- quantile(fc_perc[hi, ], probs = config$quantile_levels,
                        type = 7, na.rm = TRUE)
      rows[[idx]] <- tibble(
        origin_date     = ori_date,
        location        = loc,
        target          = "ili perc",
        horizon         = h,
        target_end_date = tgt_d,
        output_type     = "quantile",
        output_type_id  = config$quantile_levels,
        value           = as.numeric(qs)
      )
      idx <- idx + 1L
    }
  }
  
  if (length(rows) == 0) return(invisible(NULL))
  
  out <- bind_rows(rows) |>
    select(origin_date, location, target, horizon, target_end_date,
           output_type, output_type_id, value) |>
    arrange(location, horizon, output_type_id)
  
  if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)
  fname <- file.path(config$output_dir,
                     sprintf("%s-%s.csv", format(ori_date, "%Y-%m-%d"),
                             config$model_name))
  if (config$overwrite || !file.exists(fname)) {
    write_csv(out, fname)
    if (verbose) message(sprintf("  wrote %s  (%d rows)", fname, nrow(out)))
  } else {
    if (verbose) message("  exists, skipping: ", fname)
  }
  
  invisible(out)
}

# ---------------- Top-level driver ------------------------------------------

run_all_forecasts <- function(data_url, origin_dates, locations = NULL,
                              config = CONFIG, seed = 42, verbose = TRUE) {
  set.seed(seed)
  
  raw <- read_csv(data_url, show_col_types = FALSE) |>
    filter(target == "ili perc") |>
    mutate(target_end_date = as.Date(target_end_date))
  
  all_locs <- sort(unique(raw$location))
  if (verbose) {
    message("Locations in data (", length(all_locs), "): ",
            paste(all_locs, collapse = ", "))
  }
  
  if (!is.null(locations)) {
    missing <- setdiff(locations, all_locs)
    if (length(missing) > 0) {
      warning("Requested locations not in data: ",
              paste(missing, collapse = ", "))
    }
    locs <- intersect(locations, all_locs)
  } else {
    locs <- all_locs
  }
  
  ori_dates <- as.Date(origin_dates)
  all_dates <- sort(unique(raw$target_end_date))
  
  for (i_ref in seq_along(ori_dates)) {
    ori_date <- ori_dates[i_ref]
    if (verbose) {
      message(sprintf("[%d/%d] origin_date=%s",
                      i_ref, length(ori_dates), ori_date))
    }
    forecast_one_date(ori_date, raw, all_dates, locs, config, verbose)
  }
  
  message(sprintf("Done. Files written to: %s/", config$output_dir))
}

# ---------------- Entry point -----------------------------------------------

if (sys.nframe() == 0) {
  origin_dates <- c(
    "2015-10-24", "2015-10-31", "2015-11-07", "2015-11-14", "2015-11-21",
    "2015-11-28", "2015-12-05", "2015-12-12", "2015-12-19", "2015-12-26",
    "2016-01-02", "2016-01-09", "2016-01-16", "2016-01-23", "2016-01-30",
    "2016-02-06", "2016-02-13", "2016-02-20", "2016-02-27", "2016-03-05",
    "2016-03-12", "2016-03-19", "2016-03-26", "2016-04-02", "2016-04-09",
    "2016-04-16", "2016-04-23", "2016-04-30", "2016-05-07",
    "2016-10-29", "2016-11-05", "2016-11-12", "2016-11-19", "2016-11-26",
    "2016-12-03", "2016-12-10", "2016-12-17", "2016-12-24", "2016-12-31",
    "2017-01-07", "2017-01-14", "2017-01-21", "2017-01-28", "2017-02-04",
    "2017-02-11", "2017-02-18", "2017-02-25", "2017-03-04", "2017-03-11",
    "2017-03-18", "2017-03-25", "2017-04-01", "2017-04-08", "2017-04-15",
    "2017-04-22", "2017-04-29", "2017-05-06",
    "2017-10-28", "2017-11-04", "2017-11-11", "2017-11-18", "2017-11-25",
    "2017-12-02", "2017-12-09", "2017-12-16", "2017-12-23", "2017-12-30",
    "2018-01-06", "2018-01-13", "2018-01-20", "2018-01-27", "2018-02-03",
    "2018-02-10", "2018-02-17", "2018-02-24", "2018-03-03", "2018-03-10",
    "2018-03-17", "2018-03-24", "2018-03-31", "2018-04-07", "2018-04-14",
    "2018-04-21", "2018-04-28", "2018-05-05",
    "2018-10-20", "2018-10-27", "2018-11-03", "2018-11-10", "2018-11-17",
    "2018-11-24", "2018-12-01", "2018-12-08", "2018-12-15", "2018-12-22",
    "2018-12-29", "2019-01-05", "2019-01-12", "2019-01-19", "2019-01-26",
    "2019-02-02", "2019-02-09", "2019-02-16", "2019-02-23", "2019-03-02",
    "2019-03-09", "2019-03-16", "2019-03-23", "2019-03-30", "2019-04-06",
    "2019-04-13", "2019-04-20", "2019-04-27", "2019-05-04",
    "2019-10-19", "2019-10-26", "2019-11-02", "2019-11-09", "2019-11-16",
    "2019-11-23", "2019-11-30", "2019-12-07", "2019-12-14", "2019-12-21",
    "2019-12-28", "2020-01-04", "2020-01-11", "2020-01-18", "2020-01-25",
    "2020-02-01", "2020-02-08", "2020-02-15", "2020-02-22", "2020-02-29"
  )
  
  stopifnot(length(origin_dates) == 134)
  
  message(sprintf(
    "Running %d origin_dates; each file = 11 locs x 4 horizons x 23 quantiles = 1012 rows",
    length(origin_dates)))
  
  run_all_forecasts(
    data_url     = "https://raw.githubusercontent.com/sjfox/ATSF2026/refs/heads/main/target-data/time-series.csv",
    origin_dates = origin_dates,
    locations    = NULL
  )
}
