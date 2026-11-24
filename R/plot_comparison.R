# ---------------------------------------------------------------------------
# Head-to-head comparison plots of SIRS+EnKF vs SARIMA, built from the score
# tables written by evaluate_forecasts.R.
#
# Usage (from the repository root):  Rscript R/plot_comparison.R
# Writes results/comparison/{wis_by_horizon,wis_by_location}.png
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

OUT_DIR <- file.path("results", "comparison")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MODEL_COLORS <- c("SIRS + EnKF" = "#2b6cb0", "SARIMA" = "#dd6b20")

read_scores <- function(table) {
  bind_rows(
    read_csv(file.path("results/SIRSEnKF/tables", table), show_col_types = FALSE) |>
      mutate(model = "SIRS + EnKF"),
    read_csv(file.path("results/SARIMA/tables", table), show_col_types = FALSE) |>
      mutate(model = "SARIMA")
  ) |>
    mutate(model = factor(model, levels = names(MODEL_COLORS)))
}

# ---- WIS by horizon --------------------------------------------------------

by_h <- read_scores("scores_by_horizon.csv") |>
  group_by(horizon) |>
  mutate(label_vjust = ifelse(mean_WIS == max(mean_WIS), -1, 2)) |>  # higher label above, lower below
  ungroup()

p_h <- ggplot(by_h, aes(horizon, mean_WIS, color = model)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.2f", mean_WIS), vjust = label_vjust),
            size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = MODEL_COLORS) +
  scale_x_continuous(breaks = 1:4, labels = paste(1:4, "wk")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Forecast error grows faster for the mechanistic model",
       subtitle = "Mean weighted interval score (lower is better), 5,896 forecasts",
       x = "Forecast horizon", y = "Mean WIS", color = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "wis_by_horizon.png"), p_h,
       width = 7, height = 4.5, dpi = 150, bg = "white")

# ---- WIS by location (dumbbell) -------------------------------------------

by_loc <- read_scores("scores_by_location.csv")
loc_order <- by_loc |> filter(model == "SIRS + EnKF") |> arrange(mean_WIS) |> pull(location)
by_loc <- by_loc |> mutate(location = factor(location, levels = loc_order))

p_loc <- ggplot(by_loc, aes(mean_WIS, location)) +
  geom_line(aes(group = location), color = "grey75", linewidth = 1.5) +
  geom_point(aes(color = model), size = 3.5) +
  scale_color_manual(values = MODEL_COLORS) +
  scale_x_continuous(limits = c(0, NA)) +
  labs(title = "SARIMA is better in 10 of 11 locations",
       subtitle = "Mean WIS by location (lower is better); Region 6 is hardest for both",
       x = "Mean WIS", y = NULL, color = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "wis_by_location.png"), p_loc,
       width = 7, height = 5, dpi = 150, bg = "white")

cat("Plots saved to:", OUT_DIR, "\n")
