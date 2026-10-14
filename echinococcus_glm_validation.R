# =============================================================================
# Echinococcus GLM
#
# Purpose:
#   Model Echinococcus infection status in coyotes as a function of
#   surrounding land cover, age class, and sex using a binomial GLM.
#   Candidate models are compared with AICc using all-subsets selection
#   (MuMIn::dredge), and the best-fit model is checked for collinearity
#   and goodness of fit (DHARMa).
#
# Analysis decisions:
#   - Land-cover predictors: proportion developed, open, shrub, agriculture,
#     and grassland. Forest, wetland, ndvi, and ndmi were
#     excluded because they are strongly correlated with proportion developed.
#   - No random effect for individual: only 2 of 87 individuals were sampled
#     twice, so a standard GLM is used rather than a GLMM.
#   - Continuous predictors are centered and scaled (mean 0, SD 1; columns
#     suffixed "_scl"), so odds ratios are per 1-SD increase in the predictor.
#     Categorical predictors (age_class, sex) are not scaled.
#   - Average soil value (avg_soil_scl) has missing values, so models that
#     include it are fit in a separate selection run on complete cases only.
#
# Input (path set in CONFIGURATION below):
#   CSV with one row per sample and the columns:
#     echinococcus          - infection status (0 = negative, 1 = positive)
#     age_class, sex        - categorical covariates
#     prop_developed_scl, prop_open_scl, prop_shrub_scl,
#     prop_agriculture_scl, prop_grassland_scl, avg_soil_scl
#                           - scaled land-cover and soil predictors
#
# Output (written to OUTPUT_DIR):
#   best_model_coefficients.csv - estimates, SEs, Wald 95% CIs, p-values
#
# Requirements:
#   R (>= 4.1), dplyr, MuMIn, broom, performance, usdm, DHARMa
#
# Author:  Laurel Serieys
# Note:    Code was refactored and documented with assistance from
#          Claude (Anthropic), an AI assistant. All analytical decisions,
#          parameter values, and validation are the author's.
# License: see LICENSE file in this repository
# =============================================================================


# ---- 1. SETUP ---------------------------------------------------------------

  library(dplyr)
  library(MuMIn)
  library(broom)
  library(performance)
  library(usdm)
  library(DHARMa)


# Paths are relative to the repository root. Edit to match your local layout.
echino_data <- file.path("data", "echinococcus_data.csv")
OUTPUT_DIR  <- "output"

# Seed for DHARMa residual simulations, so diagnostics are reproducible.
SEED <- 2026

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ---- 2. READ AND PREPARE DATA -----------------------------------------------

if (!file.exists(echino_data)) stop("Input file not found: ", echino_data)

df <- read.csv(echino_data, stringsAsFactors = FALSE)

landcover_vars <- c("prop_developed_scl", "prop_open_scl", "prop_shrub_scl",
                    "prop_agriculture_scl", "prop_grassland_scl")

required_cols <- c("echinococcus", "age_class", "sex", landcover_vars, "avg_soil_scl")
missing_cols  <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

df <- df %>%
  mutate(
    age_class = factor(age_class),
    sex       = factor(sex)
  )

stopifnot(all(df$echinococcus %in% c(0, 1)))


# ---- 3. MODEL SELECTION: LAND COVER, AGE CLASS, AND SEX ---------------------
# Fit the global model, then compare all subsets by AICc.
# na.action = na.fail is required by dredge() so that every candidate model
# is fit to exactly the same rows.

global_mod <- glm(
  echinococcus ~ prop_developed_scl + prop_open_scl + prop_shrub_scl +
    prop_agriculture_scl + prop_grassland_scl +
    age_class + sex,
  family    = binomial,
  na.action = na.fail,
  data      = df
)

summary(global_mod)

global_dredge <- dredge(global_mod, rank = "AICc", trace = TRUE)

# Candidate models within 2 AICc units of the top model.
top_models <- get.models(global_dredge, subset = delta < 2)
print(subset(global_dredge, delta < 2))
print(top_models)


# ---- 4. MODEL SELECTION INCLUDING SOIL (COMPLETE CASES ONLY) ----------------
# avg_soil_scl has 5 missing values, so this selection run uses only rows with
# complete data for every candidate predictor.

soil_predictors <- c(landcover_vars, "avg_soil_scl", "age_class", "sex")

df_complete <- df[complete.cases(df[, c("echinococcus", soil_predictors)]), ]

message(nrow(df) - nrow(df_complete), " of ", nrow(df),
        " rows dropped due to NA in a candidate predictor.")

global_soil_mod <- glm(
  echinococcus ~ prop_developed_scl + prop_open_scl + prop_shrub_scl +
    prop_agriculture_scl + prop_grassland_scl + avg_soil_scl +
    age_class + sex,
  family    = binomial,
  na.action = na.fail,
  data      = df_complete
)

soil_dredge <- dredge(global_soil_mod, rank = "AICc", trace = TRUE)

top_soil_models <- get.models(soil_dredge, subset = delta < 2)
print(subset(soil_dredge, delta < 2))
print(top_soil_models)


# ---- 5. BEST-FIT MODEL ------------------------------------------------------
# Best-fit model identified by the AICc model selection above.

best_mod <- glm(
  echinococcus ~ prop_open_scl + prop_agriculture_scl,
  family    = binomial,
  na.action = na.fail,
  data      = df
)

summary(best_mod)


# ---- 6. COLLINEARITY CHECK --------------------------------------------------
# Variance inflation factors for the predictors in the best-fit model.

usdm::vif(as.data.frame(df[, c("prop_open_scl", "prop_agriculture_scl")]))

print(check_collinearity(best_mod))


# ---- 7. EFFECT SIZES --------------------------------------------------------
# Odds ratios with profile-likelihood 95% CIs (per 1-SD increase in predictor).

odds_ratios <- tidy(best_mod, conf.int = TRUE) %>%
  mutate(
    odds_ratio = exp(estimate),
    or_low     = exp(conf.low),
    or_high    = exp(conf.high)
  ) %>%
  select(term, odds_ratio, or_low, or_high, p.value)

print(odds_ratios)

# Coefficient table (log-odds scale) with Wald 95% CIs, exported for reporting.
coef_table <- summary(best_mod)$coefficients

best_mod_results <- data.frame(
  Covariate = rownames(coef_table),
  Estimate  = coef_table[, "Estimate"],
  Std_Error = coef_table[, "Std. Error"],
  Lower_CI  = coef_table[, "Estimate"] - 1.96 * coef_table[, "Std. Error"],
  Upper_CI  = coef_table[, "Estimate"] + 1.96 * coef_table[, "Std. Error"],
  p_value   = coef_table[, "Pr(>|z|)"],
  row.names = NULL
)

best_mod_results[, 2:5] <- round(best_mod_results[, 2:5], 2)

write.csv(best_mod_results,
          file.path(OUTPUT_DIR, "best_model_coefficients.csv"),
          row.names = FALSE)


# ---- 8. GOODNESS OF FIT (DHARMa) --------------------------------------------
# Simulation-based scaled residuals. For an adequately fitting model:
#   - QQ plot points fall along the 1:1 line, and the residual vs. predicted
#     plot shows no trend, funnel, or other pattern;
#   - uniformity (KS), dispersion, and outlier tests are non-significant
#     (p > 0.05).

set.seed(SEED)
sim_res <- simulateResiduals(fittedModel = best_mod, n = 1000, plot = FALSE)

plot(sim_res)

# Uniformity, dispersion, and outlier tests.
testResiduals(sim_res)

# Residuals against each predictor, to check for nonlinearity or missed structure.
plotResiduals(sim_res, form = df$prop_open_scl)
plotResiduals(sim_res, form = df$prop_agriculture_scl)

