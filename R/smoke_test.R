# Quick smoke test: 2 origin_dates, all locations, small ensemble.
# Writes 2 files to model-output/smoke/; each should have 1012 rows.
# Usage (from the repository root):  Rscript R/smoke_test.R

source("R/run_forecast.R")

test_dates <- c("2015-12-26", "2016-02-06")

test_config <- modifyList(CONFIG, list(
  N_ensemble = 100,
  output_dir = "model-output/smoke",
  model_name = "SIRS-EnKF-SMOKE"
))

run_all_forecasts(
  data_url     = "https://raw.githubusercontent.com/sjfox/ATSF2026/refs/heads/main/target-data/time-series.csv",
  origin_dates = test_dates,
  locations    = NULL,
  config       = test_config
)

cat("\n=== Verification ===\n")
files <- list.files(test_config$output_dir, pattern = "SMOKE\\.csv$", full.names = TRUE)
cat(sprintf("Files written: %d (expected 2)\n", length(files)))

for (f in files) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  cat(sprintf("\n%s\n  rows: %d (expected 1012)\n  columns: %s\n",
              basename(f), nrow(df), paste(colnames(df), collapse = ", ")))
  cat("  first 3 rows:\n")
  print(head(df, 3))
  cat(sprintf("  unique locations: %d  horizons: %s  quantiles: %d\n",
              length(unique(df$location)),
              paste(sort(unique(df$horizon)), collapse = ","),
              length(unique(df$output_type_id))))
  stopifnot(nrow(df) == 1012)
}
