# Quick smoke test: 2 origin_dates, all locations, small ensemble.
# Writes 2 files to forecasts_smoke/; each should have 1012 rows.

source("C:/Users/syed0/INF599/Final_project/files/run_forecast.R")

test_dates <- c("2015-12-26", "2016-02-06")

test_config <- modifyList(CONFIG, list(
  N_ensemble = 100,
  output_dir = "C:/Users/syed0/INF599/Final_project/forecasts_smoke",
  model_name = "SIHR-EnKF-SMOKE"
))

run_all_forecasts(
  data_url     = "https://raw.githubusercontent.com/sjfox/ATSF2026/refs/heads/main/target-data/time-series.csv",
  origin_dates = test_dates,
  locations    = NULL,
  config       = test_config
)

cat("\n=== Verification ===\n")
files <- list.files("Final_project/forecasts_smoke", full.names = TRUE)
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
}





sarima_fc <- list.files("C:/Users/syed0/INF599/Final_project/forecasts_SARIMA",
                        pattern = "\\.csv$", full.names = TRUE) |>
  lapply(read_csv, show_col_types = FALSE) |> bind_rows()
sirs_fc <- list.files("C:/Users/syed0/INF599/Final_project/forecasts_SIRS",
                      pattern = "\\.csv$", full.names = TRUE) |>
  lapply(read_csv, show_col_types = FALSE) |> bind_rows()

cat("SARIMA forecast rows:", nrow(sarima_fc), "\n")
cat("SIRS forecast rows:", nrow(sirs_fc), "\n")
cat("Common (origin, location, horizon):\n")
sarima_keys <- sarima_fc |> distinct(origin_date, location, horizon)
sirs_keys   <- sirs_fc   |> distinct(origin_date, location, horizon)
cat("  in both:", nrow(inner_join(sarima_keys, sirs_keys, by = c("origin_date","location","horizon"))), "\n")
cat("  SIRS only:", nrow(anti_join(sirs_keys, sarima_keys, by = c("origin_date","location","horizon"))), "\n")








# ---------------------------------------------------------------------------
# Rename forecast files to add "syed-" prefix before model name.
#   2015-10-24-SIRS-EnKF.csv  -> 2015-10-24-syed-SIRS-EnKF.csv
#   2015-10-24-SARIMA.csv     -> 2015-10-24-syed-SARIMA.csv
# Run this once per directory.
# ---------------------------------------------------------------------------


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

ROOT <- "C:/Users/syed0/INF599/Final_project"

SIRS_OLD <- "syed-SIRS-EnKF"   # current folder/metadata name
SIRS_NEW <- "syed-SIRSEnKF"    # target name (no hyphens)

# ---------- Step 1: Rename SIRS CSV files inside the OLD folder -------------

old_sirs_dir <- file.path(ROOT, SIRS_OLD)

if (dir.exists(old_sirs_dir)) {
  files <- list.files(old_sirs_dir, pattern = "\\.csv$", full.names = TRUE)
  cat(sprintf("\nProcessing %d files in %s\n", length(files), old_sirs_dir))
  
  n_renamed <- 0; n_already_ok <- 0; n_unknown <- 0
  for (f in files) {
    base <- basename(f)
    
    # Already in the target form?
    if (grepl("-syed-SIRSEnKF\\.csv$", base)) {
      n_already_ok <- n_already_ok + 1
      next
    }
    
    # Original form: 2015-10-24-SIRS-EnKF.csv
    if (grepl("-SIRS-EnKF\\.csv$", base)) {
      new_base <- sub("-SIRS-EnKF\\.csv$", "-syed-SIRSEnKF.csv", base)
      file.rename(f, file.path(old_sirs_dir, new_base))
      n_renamed <- n_renamed + 1
      next
    }
    
    # Already had syed but wrong form (e.g. -syed-SIRS-EnKF.csv)
    if (grepl("-syed-SIRS-EnKF\\.csv$", base)) {
      new_base <- sub("-syed-SIRS-EnKF\\.csv$", "-syed-SIRSEnKF.csv", base)
      file.rename(f, file.path(old_sirs_dir, new_base))
      n_renamed <- n_renamed + 1
      next
    }
    
    cat(sprintf("  unrecognized filename: %s\n", base))
    n_unknown <- n_unknown + 1
  }
  cat(sprintf("  renamed: %d   already-ok: %d   unrecognized: %d\n",
              n_renamed, n_already_ok, n_unknown))
}

# ---------- Step 2: Rename the SIRS folder itself ---------------------------

new_sirs_dir <- file.path(ROOT, SIRS_NEW)

if (dir.exists(old_sirs_dir) && !dir.exists(new_sirs_dir)) {
  cat(sprintf("\nRenaming folder:\n  %s\n  -> %s\n", old_sirs_dir, new_sirs_dir))
  file.rename(old_sirs_dir, new_sirs_dir)
}




























