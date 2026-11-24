# ---------------------------------------------------------------------------
# Evaluate 134 SIHR+EnKF forecast files against ground truth.
#
# Usage: setwd() to the folder containing this script (your `files/` folder),
#        make sure the forecasts are at the output_dir path below,
#        then source("evaluate_forecasts.R").
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
})

# ---- Config ----------------------------------------------------------------

FORECAST_DIR <- "C:/Users/syed0/INF599/Final_project/forecasts_SARIMA"
TRUTH_URL    <- "https://raw.githubusercontent.com/sjfox/ATSF2026/refs/heads/main/target-data/time-series.csv"
TABLES_DIR   <- "C:/Users/syed0/INF599/Final_project/eval_tables_SARIMA"
PLOTS_DIR    <- "C:/Users/syed0/INF599/Final_project/eval_plots_SARIMA"

dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR,  showWarnings = FALSE, recursive = TRUE)

# ---- Load all forecast files ----------------------------------------------

forecast_files <- list.files(FORECAST_DIR, pattern = "\\.csv$", full.names = TRUE)
cat(sprintf("Found %d forecast files\n", length(forecast_files)))

forecasts <- forecast_files |>
  lapply(read_csv, show_col_types = FALSE) |>
  bind_rows() |>
  mutate(origin_date     = as.Date(origin_date),
         target_end_date = as.Date(target_end_date))

cat(sprintf("Loaded %d total forecast rows\n", nrow(forecasts)))

# ---- Sanity checks ---------------------------------------------------------

cat("\n=== SANITY CHECKS ===\n")
rows_per_file <- forecasts |> count(origin_date) |> pull(n)
cat(sprintf("Rows per file: min=%d max=%d (expected 1012)\n",
            min(rows_per_file), max(rows_per_file)))
cat(sprintf("NaN/NA values: %d\n",
            sum(is.na(forecasts$value) | is.nan(forecasts$value))))
cat(sprintf("Negative values: %d\n", sum(forecasts$value < 0, na.rm = TRUE)))

# ---- Load truth -----------------------------------------------------------

truth <- read_csv(TRUTH_URL, show_col_types = FALSE) |>
  filter(target == "ili perc") |>
  mutate(target_end_date = as.Date(target_end_date)) |>
  select(location, target_end_date, truth = observation)

# ---- Pivot forecasts wide and merge truth ---------------------------------

fc_wide <- forecasts |>
  select(origin_date, location, horizon, target_end_date,
         output_type_id, value) |>
  pivot_wider(names_from = output_type_id, values_from = value,
              names_prefix = "q") |>
  left_join(truth, by = c("location", "target_end_date"))

# ---- WIS with decomposition -----------------------------------------------

# WIS per Bracher et al. 2021. We also return the three components:
#   sharpness   = sum of weighted interval widths (cost of being uncertain)
#   over_pred   = penalty when truth < lower bound (forecast too high)
#   under_pred  = penalty when truth > upper bound (forecast too low)
# For probabilistic forecasts these decomposed pieces are diagnostic gold.
wis_with_components <- function(row) {
  y <- row$truth
  if (is.na(y)) return(c(WIS = NA_real_, sharp = NA_real_,
                         over = NA_real_, under = NA_real_))
  m <- row$q0.5
  
  alphas <- c(0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50,
              0.60, 0.70, 0.80, 0.90)
  qlow_names  <- c("q0.01","q0.025","q0.05","q0.1","q0.15","q0.2","q0.25",
                   "q0.3","q0.35","q0.4","q0.45")
  qhigh_names <- c("q0.99","q0.975","q0.95","q0.9","q0.85","q0.8","q0.75",
                   "q0.7","q0.65","q0.6","q0.55")
  
  K <- length(alphas)
  width <- numeric(K); over <- numeric(K); under <- numeric(K)
  for (i in seq_along(alphas)) {
    a <- alphas[i]
    l <- row[[qlow_names[i]]];  u <- row[[qhigh_names[i]]]
    width[i] <- u - l
    over[i]  <- if (y < l) (2 / a) * (l - y) else 0
    under[i] <- if (y > u) (2 / a) * (y - u) else 0
  }
  
  norm <- 1 / (K + 0.5)
  med_term  <- 0.5 * abs(y - m)
  sharp     <- norm * sum((alphas / 2) * width)
  over_p    <- norm * sum((alphas / 2) * over)
  under_p   <- norm * sum((alphas / 2) * under)
  wis_total <- norm * med_term + sharp + over_p + under_p
  
  c(WIS = wis_total, sharp = sharp, over = over_p, under = under_p)
}

cat("\nComputing WIS with decomposition ...\n")
score_mat <- t(vapply(seq_len(nrow(fc_wide)),
                      function(i) wis_with_components(fc_wide[i, ]),
                      numeric(4)))
fc_wide$WIS         <- score_mat[, "WIS"]
fc_wide$WIS_sharp   <- score_mat[, "sharp"]
fc_wide$WIS_over    <- score_mat[, "over"]
fc_wide$WIS_under   <- score_mat[, "under"]

# Point-forecast errors and coverage
fc_wide <- fc_wide |>
  mutate(
    sq_err     = (truth - q0.5)^2,
    abs_err    = abs(truth - q0.5),
    in_50      = truth >= q0.25  & truth <= q0.75,
    in_95      = truth >= q0.025 & truth <= q0.975
  )

# ---- Summary helpers ------------------------------------------------------

summarize_block <- function(df, ...) {
  df |> filter(!is.na(WIS)) |>
    group_by(...) |>
    summarise(
      n            = n(),
      mean_WIS     = mean(WIS),
      WIS_sharp    = mean(WIS_sharp),
      WIS_over     = mean(WIS_over),
      WIS_under    = mean(WIS_under),
      RMSE         = sqrt(mean(sq_err)),
      MAE          = mean(abs_err),
      coverage_50  = mean(in_50),
      coverage_95  = mean(in_95),
      .groups = "drop"
    )
}

# ---- Compute and write tables ---------------------------------------------

overall      <- summarize_block(fc_wide)
by_horizon   <- summarize_block(fc_wide, horizon)
by_location  <- summarize_block(fc_wide, location) |> arrange(mean_WIS)
by_loc_hor   <- summarize_block(fc_wide, location, horizon)
by_origin    <- summarize_block(fc_wide, origin_date)

write_csv(overall,     file.path(TABLES_DIR, "scores_overall.csv"))
write_csv(by_horizon,  file.path(TABLES_DIR, "scores_by_horizon.csv"))
write_csv(by_location, file.path(TABLES_DIR, "scores_by_location.csv"))
write_csv(by_loc_hor,  file.path(TABLES_DIR, "scores_by_location_horizon.csv"))
write_csv(by_origin,   file.path(TABLES_DIR, "scores_by_origin_date.csv"))

# Also save the per-forecast scores for fine-grained analysis
write_csv(fc_wide |> select(origin_date, location, horizon, target_end_date,
                            truth, q0.025, q0.25, q0.5, q0.75, q0.975,
                            WIS, WIS_sharp, WIS_over, WIS_under,
                            abs_err, sq_err, in_50, in_95),
          file.path(TABLES_DIR, "scores_per_forecast.csv"))

# ---- Print to console ------------------------------------------------------

cat("\n=== OVERALL SCORES ===\n")
print(overall)

cat("\n=== BY HORIZON (lower WIS / RMSE = better) ===\n")
print(by_horizon)

cat("\n=== BY LOCATION (sorted by WIS) ===\n")
print(by_location)

# ---- Diagnostic: WIS decomposition plot -----------------------------------

# Stacked bar of mean WIS components by horizon
wis_decomp_long <- by_horizon |>
  select(horizon, sharpness = WIS_sharp,
         over_prediction = WIS_over,
         under_prediction = WIS_under) |>
  pivot_longer(-horizon, names_to = "component", values_to = "value")

p_decomp <- ggplot(wis_decomp_long,
                   aes(x = factor(horizon), y = value, fill = component)) +
  geom_bar(stat = "identity") +
  labs(title = "WIS decomposition by horizon",
       x = "horizon (weeks ahead)", y = "mean WIS contribution",
       fill = NULL) +
  scale_fill_manual(values = c(sharpness = "grey60",
                               over_prediction = "tomato",
                               under_prediction = "steelblue")) +
  theme_minimal()
ggsave(file.path(PLOTS_DIR, "wis_decomposition_by_horizon_SARIMA.png"),
       p_decomp, width = 7, height = 4, dpi = 120)

# WIS by location (bar chart, sorted)
p_loc <- by_location |>
  mutate(location = factor(location, levels = location)) |>
  ggplot(aes(location, mean_WIS)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = sprintf("%.2f", mean_WIS)),
            hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title = "Mean WIS by location", x = NULL, y = "mean WIS") +
  theme_minimal()
ggsave(file.path(PLOTS_DIR, "wis_by_location_SARIMA.png"),
       p_loc, width = 7, height = 4, dpi = 120)

# RMSE vs WIS by location (do they tell the same story?)
p_rmse_wis <- by_location |>
  ggplot(aes(mean_WIS, RMSE, label = location)) +
  geom_point(color = "steelblue", size = 3) +
  geom_text(hjust = -0.1, size = 3) +
  labs(title = "RMSE vs mean WIS across locations of SARIMA",
       x = "mean WIS", y = "RMSE") +
  theme_minimal()
ggsave(file.path(PLOTS_DIR, "rmse_vs_wis_by_location_SARIMA.png"),
       p_rmse_wis, width = 7, height = 5, dpi = 120)

cat(sprintf("\nTables saved to: %s\n", TABLES_DIR))
cat(sprintf("Plots  saved to: %s\n", PLOTS_DIR))
cat("Done.\n")
