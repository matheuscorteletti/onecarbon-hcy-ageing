# =============================================================================
# 01_models.R — shared analytic dataset and hierarchical models for H1
# Sourced by Figures 4–5 and Table 3 (after 00_setup.R)
# =============================================================================

# Levodopa users (n = 5 in the model set) are excluded in a sensitivity analysis rather than adjusted for
# ---- Complete-case analytic set for the H1 models (hospital, 3 strata) -------
MODEL_VARS <- c("log_hcy", "hcy", "stratum", "sex", "ancestry_pc1", "ancestry_pc2",
                "time_to_centrifugation_min", "fasting_hours",
                "log_egfr", "egfr", "log_holotc", "holotc", "log_rbcf", "rbc_fol", "log_plp", "egrac",
                "glim_any", "phase_angle_deg", "log_crp", "diabetes_t2",
                "hcy_genetic_risk_score", "log_mma", "lcms_batch")
mdat <- hosp %>%
  select(participant_id, all_of(MODEL_VARS), ffq_reliable, supp_vitamin_b12, sens_exclude_b12_antifolate,
         egfr_ckdepi2021, glim3, mthfr677, grs_tert, holotc_low, hhcy, acute_kidney_injury,
         infection_at_admission, med_levodopa) %>%
  drop_na(all_of(MODEL_VARS)) %>%
  mutate(s76 = as.integer(stratum == STRATA[2]), s85 = as.integer(stratum == STRATA[3]),
         male = as.integer(sex == "Male"))

# Blocks entered hierarchically (Table 3 / attenuation plot)
BLOCKS <- list(
  "M1 Age + design"     = c("s76", "s85", "male", "ancestry_pc1", "ancestry_pc2",
                            "time_to_centrifugation_min", "fasting_hours"),
  "M2 + Kidney"         = "log_egfr",
  "M3 + Active B12"     = "log_holotc",
  "M4 + Folate, B6, B2" = c("log_rbcf", "log_plp", "egrac"),
  "M5 + Nutrition, inflammation, diabetes, GRS" = c("glim_any", "phase_angle_deg", "log_crp",
                                                    "diabetes_t2", "hcy_genetic_risk_score"))
block_formula <- function(k) reformulate(unlist(BLOCKS[seq_len(k)]), response = "log_hcy")
hier_models <- map(seq_along(BLOCKS), ~ lm(block_formula(.x), data = mdat)) %>% set_names(names(BLOCKS))
FULL_FORMULA <- block_formula(length(BLOCKS))

# Determinant groups for decomposition / relative importance
DET_GROUPS <- list(
  "Kidney (eGFR-cys)"   = "log_egfr",
  "Active B12 (holoTC)" = "log_holotc",
  "Folate (RBC)"        = "log_rbcf",
  "B6 + B2"             = c("log_plp", "egrac"),
  "Malnutrition (GLIM, PhA)" = c("glim_any", "phase_angle_deg"),
  "Inflammation (CRP)"  = "log_crp",
  "Genetics (GRS)"      = "hcy_genetic_risk_score",
  "Diabetes"            = "diabetes_t2",
  "Sex, ancestry, pre-analytics" = c("male", "ancestry_pc1", "ancestry_pc2", "time_to_centrifugation_min", "fasting_hours"))
PAL_DET <- c("Kidney (eGFR-cys)" = "#2B6CB0", "Active B12 (holoTC)" = "#9E2A2B", "Folate (RBC)" = "#7FB7A4",
             "B6 + B2" = "#B8A04A", "Malnutrition (GLIM, PhA)" = "#E09F3E", "Inflammation (CRP)" = "#6C5B7B",
             "Genetics (GRS)" = "#4C956C", "Diabetes" = "#A0A0A0", "Sex, ancestry, pre-analytics" = "#D0D0D0",
             "Age stratum (unexplained)" = "#3A3A3A")

# ---- g-computation model (flexible) ------------------------------------------
# Natural splines for kidney and B12; stratum interactions allowed for kidney/B12
suppressPackageStartupMessages(library(splines))
G_FORMULA <- log_hcy ~ stratum * (ns(log_egfr, 3) + ns(log_holotc, 3)) + log_rbcf + log_plp + egrac +
  glim_any + phase_angle_deg + log_crp + male + ancestry_pc1 + ancestry_pc2 + time_to_centrifugation_min +
  fasting_hours + diabetes_t2 + hcy_genetic_risk_score

# Rank-preserving transfer of the reference (younger) distribution to target rows:
# each target individual receives the reference quantile matching their own rank
quantile_map <- function(x_target, x_ref, frac = 1) {
  u <- (rank(x_target, ties.method = "average") - 0.5) / length(x_target)
  q <- unname(quantile(x_ref, u, type = 8))
  x_target + frac * (q - x_target)          # frac in [0,1]: share of the gap closed
}

# Scenario engine: returns mean Hcy and % >15 for the target stratum under interventions
gcomp_scenarios <- function(d, scenarios, target = STRATA[3], ref = STRATA[1], formula = G_FORMULA) {
  fit <- lm(formula, data = d)
  smear <- mean(exp(resid(fit)))                       # Duan smearing for back-transform
  tg <- d %>% filter(stratum == target); rf <- d %>% filter(stratum == ref)
  sig <- sigma(fit)
  map_dfr(names(scenarios), function(nm) {
    nd <- scenarios[[nm]](tg, rf)
    mu <- predict(fit, newdata = nd)
    tibble(scenario = nm, mean_hcy = mean(exp(mu) * smear),
           p_hhcy = mean(pnorm((mu - log(15)) / sig)))  # P(Hcy > 15) under log-normal error
  })
}
