# =============================================================================
# Predicted probability of Echinococcus infection in relation to open and
# agricultural land cover
#
# plot_echinococcus_open_agriculture.R
#
# Purpose:
#   Plot predicted probabilities of Echinococcus infection from the best-fit
#   binomial GLM (see echinococcus_glm.R):
#     echinococcus ~ prop_open_scl + prop_agriculture_scl
#
#   The figure has two panels, one per predictor (a: open land cover;
#   b: agriculture), each showing:
#     - the predicted probability of infection across the observed range of
#       that predictor, with the other predictor held at its mean (0 on the
#       scaled axis)
#     - a 95% confidence band around the prediction
#     - observed infection status of each sample (0/1), jittered slightly
#       vertically for visibility
#
#   The model is fit on z-scaled predictors (suffix "_scl"; see
#   extract_covariates.R). For interpretability, the x-axis is
#   back-transformed to the original land-cover proportion, using the mean
#   and standard deviation of each unscaled covariate.
#
# Input (path set in CONFIGURATION below):
#   Output of extract_covariates.R, with the columns:
#     echinococcus          - infection status (0 = negative, 1 = positive)
#     prop_open, prop_agriculture         - land-cover proportions (0-1)
#     prop_open_scl, prop_agriculture_scl - the same, z-scaled
#
# Output (written to OUTPUT_DIR):
#   echinococcus_open_agriculture_predicted_probability.png
#
# Requirements:
#   R (>= 4.1), ggeffects, ggplot2, patchwork, scales
#
# Author:  Laurel Serieys
# Contact: lekserieys@gmail.com
# Note:    Code was refactored and documented with assistance from
#          Claude (Anthropic), an AI assistant. All analytical decisions,
#          parameter values, and validation are the author's.
# License: see LICENSE file in this repository
# =============================================================================


# ---- 1. SETUP ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggeffects)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# Paths are relative to the repository root. Edit to match your local layout.
echino_data <- file.path("output", "echinococcus_data_covariates.csv")
OUTPUT_DIR  <- file.path("output", "figures")
OUTPUT_PNG  <- file.path(OUTPUT_DIR, "echinococcus_open_agriculture_predicted_probability.png")

# Panel colors.
COLOR_OPEN <- "#2b7a78"
COLOR_AG   <- "#c1440e"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ---- 2. READ DATA -----------------------------------------------------------

if (!file.exists(echino_data)) stop("Input file not found: ", echino_data)

df <- read.csv(echino_data, stringsAsFactors = FALSE)

required_cols <- c("echinococcus", "prop_open", "prop_agriculture",
                   "prop_open_scl", "prop_agriculture_scl")
missing_cols  <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0)
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))


# ---- 3. BACK-TRANSFORMATION TO ORIGINAL PROPORTIONS -------------------------
# The scaled predictors were calculated as (x - mean) / SD, so the original
# values are recovered as scaled * SD + mean. The check confirms that these
# values reproduce the scaled columns in the data.

open_center <- mean(df$prop_open, na.rm = TRUE)
open_scale  <- sd(df$prop_open,   na.rm = TRUE)
ag_center   <- mean(df$prop_agriculture, na.rm = TRUE)
ag_scale    <- sd(df$prop_agriculture,   na.rm = TRUE)

stopifnot(
  isTRUE(all.equal((df$prop_open        - open_center) / open_scale, df$prop_open_scl)),
  isTRUE(all.equal((df$prop_agriculture - ag_center)   / ag_scale,   df$prop_agriculture_scl))
)


# ---- 4. FIT MODEL AND PREDICT -----------------------------------------------

best_mod <- glm(
  echinococcus ~ prop_open_scl + prop_agriculture_scl,
  family    = binomial,
  na.action = na.fail,
  data      = df
)

summary(best_mod)

# Predicted probabilities across each predictor's observed values, with the
# other predictor held at its mean.
pred_open <- ggpredict(best_mod, terms = "prop_open_scl [all]")
pred_open$x_natural <- pred_open$x * open_scale + open_center

pred_ag <- ggpredict(best_mod, terms = "prop_agriculture_scl [all]")
pred_ag$x_natural <- pred_ag$x * ag_scale + ag_center


# ---- 5. BUILD FIGURE --------------------------------------------------------

#' One panel: predicted probability curve, 95% CI band, and observed data.
make_panel <- function(pred_df, raw_x, raw_y, x_label, color) {
  set.seed(42)  # reproducible vertical jitter of observed points
  raw_df <- data.frame(x = raw_x, y = raw_y + runif(length(raw_y), -0.03, 0.03))
  
  ggplot() +
    geom_ribbon(data = pred_df, aes(x = x_natural, ymin = conf.low, ymax = conf.high),
                fill = color, alpha = 0.18) +
    geom_line(data = pred_df, aes(x = x_natural, y = predicted),
              color = color, linewidth = 1.1) +
    geom_point(data = raw_df, aes(x = x, y = y),
               color = "grey30", alpha = 0.55, size = 1.8) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = x_label, y = "Predicted probability of infection") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
}

p_open <- make_panel(pred_open, df$prop_open, df$echinococcus,
                     "Proportion altered open land cover", COLOR_OPEN)
p_ag   <- make_panel(pred_ag, df$prop_agriculture, df$echinococcus,
                     "Proportion agriculture land cover", COLOR_AG)

combined <- (p_open | p_ag) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(face = "bold", size = 14))
  )

print(combined)


# ---- 6. SAVE ----------------------------------------------------------------

ggsave(OUTPUT_PNG, combined, width = 11, height = 5, dpi = 300, bg = "white")
message("Saved: ", OUTPUT_PNG)
