# ---------------------------------------------------------------------------
# ATSF2026 SARIMA(0,1,1)(0,1,1)[52] forecasting driver
# ---------------------------------------------------------------------------
# Statistical baseline for comparison with SIRS+EnKF. For each origin_date
# and each location independently:
#   1. Take all ILI% history for that location strictly BEFORE origin_date
#      (expanding-window: SARIMA needs multi-season data to learn seasonality)
#   2. Convert to weekly ts() with frequency 52
#   3. Fit fixed SARIMA(0,1,1)(0,1,1)[52] on log(ILI%)
#   4. Forecast 4 weeks ahead -> mean and sd of Gaussian predictive
#   5. Map to 23 ATSF2026 quantiles via qnorm() on log scale, exp-back to %
#
# Output: ONE CSV per origin_date in output_dir, 1012 rows each, same
# column layout as the SIRS pipeline. Filename: {origin_date}-syed-SARIMA.csv
#
# Usage (from the repository root):  Rscript R/run_sarima.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(lubridate)
  library(forecast)        # for Arima() and forecast()
})

# ---------------- Configuration ---------------------------------------------

CONFIG <- list(
  horizons        = 1:4,

  # Fixed SARIMA orders (FluSight-standard "airline" model)
  arima_order     = c(0, 1, 1),
  seasonal_order  = c(0, 1, 1),
  seasonal_period = 52,

  # Floor before taking log (avoids log(0); ILI% in data is always > 0 in
  # practice, but this guards against edge cases)
  log_floor       = 0.05,

  quantile_levels = c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99),

  output_dir      = "model-output/syed-SARIMA",
  model_name      = "syed-SARIMA",
  overwrite       = TRUE,

  # Minimum number of observations required to fit. With seasonal period 52
  # and one round of seasonal differencing, you need at least ~52 weeks
  # before the model is well-defined.
  min_obs         = 60
)

# ---------------- Helpers ---------------------------------------------------

#' Find the appropriate weekly frequency for a series of weekly ILI dates.
#' We use frequency 52 (52.18 weeks/year, but 52 is what FluSight conventionally uses).
make_ts <- function(values, dates) {
  # Compute the start year and week-of-year of the first observation
  first <- min(dates)
  yr   <- as.integer(format(first, "%Y"))
  wk   <- as.integer(format(first, "%U"))
  if (wk == 0) wk <- 1
  ts(values, start = c(yr, wk), frequency = 52)
}

#' Forecast 4 weeks ahead from a fitted Arima object on log scale.
#' Returns a [n_horizons x n_quantiles] matrix of predicted ILI%.
sarima_quantile_forecast <- function(fit, horizons, quantile_levels) {
  fc <- forecast(fit, h = max(horizons), level = 95)
  # `forecast` returns mean and lower/upper at one level, but we can
  # recover the predictive sd from the level-95 interval
  z95 <- qnorm(0.975)
  pred_mean_log <- as.numeric(fc$mean)
  # upper95 - mean = z95 * sd  =>  sd = (upper95 - mean) / z95
  pred_sd_log <- (as.numeric(fc$upper[, 1]) - pred_mean_log) / z95

  out <- matrix(0, nrow = length(horizons), ncol = length(quantile_levels))
  for (i in seq_along(horizons)) {
    h <- horizons[i]
    z <- qnorm(quantile_levels)
    # quantiles on log scale, then exp-back to % scale
    log_q   <- pred_mean_log[h] + z * pred_sd_log[h]
    out[i, ] <- exp(log_q)
  }
  out
}

# ---------------- Per-date forecast (writes one CSV) -----------------------

forecast_one_date <- function(ori_date, raw, locs, config, verbose = TRUE) {
  ori_date <- as.Date(ori_date)
  rows <- list()
  idx  <- 1L

  for (loc in locs) {
    hist <- raw |>
      filter(location == loc, target_end_date < ori_date) |>
      arrange(target_end_date)

    if (nrow(hist) < config$min_obs) {
      if (verbose) message(sprintf("  skip %s (only %d obs, need %d)",
                                   loc, nrow(hist), config$min_obs))
      next
    }

    # Take logs after flooring
    y <- pmax(hist$observation, config$log_floor)
    log_y <- log(y)
    ts_log <- make_ts(log_y, hist$target_end_date)

    fit <- tryCatch(
      Arima(ts_log,
            order          = config$arima_order,
            seasonal       = list(order = config$seasonal_order,
                                  period = config$seasonal_period),
            method         = "CSS-ML"),
      error = function(e) {
        message(sprintf("  Arima failed for %s on %s: %s",
                        loc, ori_date, e$message))
        NULL
      }
    )
    if (is.null(fit)) next

    quant_mat <- tryCatch(
      sarima_quantile_forecast(fit, config$horizons, config$quantile_levels),
      error = function(e) {
        message(sprintf("  forecast failed for %s on %s: %s",
                        loc, ori_date, e$message))
        NULL
      }
    )
    if (is.null(quant_mat)) next

    for (hi in seq_along(config$horizons)) {
      h     <- config$horizons[hi]
      tgt_d <- ori_date + 7 * h
      rows[[idx]] <- tibble(
        origin_date     = ori_date,
        location        = loc,
        target          = "ili perc",
        horizon         = h,
        target_end_date = tgt_d,
        output_type     = "quantile",
        output_type_id  = config$quantile_levels,
        value           = as.numeric(quant_mat[hi, ])
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
  }

  invisible(out)
}

# ---------------- Top-level driver ------------------------------------------

run_all_sarima <- function(data_url, origin_dates, locations = NULL,
                           config = CONFIG, verbose = TRUE) {
  raw <- read_csv(data_url, show_col_types = FALSE) |>
    filter(target == "ili perc") |>
    mutate(target_end_date = as.Date(target_end_date))

  all_locs <- sort(unique(raw$location))
  if (verbose) {
    message("Locations in data (", length(all_locs), "): ",
            paste(all_locs, collapse = ", "))
  }

  if (!is.null(locations)) {
    locs <- intersect(locations, all_locs)
  } else {
    locs <- all_locs
  }

  ori_dates <- as.Date(origin_dates)

  for (i_ref in seq_along(ori_dates)) {
    ori_date <- ori_dates[i_ref]
    if (verbose) {
      message(sprintf("[%d/%d] origin_date=%s",
                      i_ref, length(ori_dates), ori_date))
    }
    forecast_one_date(ori_date, raw, locs, config, verbose)
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
    "Running SARIMA on %d origin_dates; each file = 11 locs x 4 horizons x 23 quantiles = 1012 rows",
    length(origin_dates)))

  run_all_sarima(
    data_url     = "https://raw.githubusercontent.com/sjfox/ATSF2026/refs/heads/main/target-data/time-series.csv",
    origin_dates = origin_dates,
    locations    = NULL
  )
}
