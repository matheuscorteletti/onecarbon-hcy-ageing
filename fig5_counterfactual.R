# =============================================================================
# Figure 5 — Counterfactual: what if the old had the kidney and B12 of the younger? (Q4; H1a, H1e)
#   P1  g-computation: mean Hcy under rank-preserving interventions (76–84 y and >=85 y; bootstrap CI)
#   P2  Prevalence of hyperhomocysteinaemia (>15 umol/L) under each scenario
#   P3  Individual predicted Hcy distribution: observed vs "younger kidney + B12" (>=85 y)
#   P4  Dose–response: share of the eGFR / holoTC gap closed -> mean Hcy (>=85 y)
#   P5  Interventional mediation (waterfall): age -> kidney / B12 / nutrition -> Hcy
#   P6  Robustness of the main counterfactual across sensitivity analyses
# Interventions: each old participant receives the younger-stratum quantile matching their own rank
# (rank-preserving distributional shift; Duan smearing for back-transformation)
# =============================================================================
source("R/00_setup.R"); source("R/01_models.R")
suppressPackageStartupMessages(library(ggdist))

ref_dist <- function(rf, v) rf[[v]]
set_young <- function(vars) function(tg, rf) {
  for (v in vars) tg[[v]] <- quantile_map(tg[[v]], rf[[v]])
  if ("glim_any" %in% vars) {                       # binary: match younger prevalence, worst PhA first spared
    k <- round(mean(rf$glim_any) * nrow(tg)); o <- order(tg$phase_angle_deg)
    tg$glim_any <- 0L; tg$glim_any[o[seq_len(k)]] <- 1L
  }
  tg
}
SCEN <- list(
  "Observed (model)"            = function(tg, rf) tg,
  "Younger RBC folate"          = set_young("log_rbcf"),
  "Younger eGFR"                = set_young("log_egfr"),
  "Younger holoTC"              = set_young("log_holotc"),
  "Younger eGFR + holoTC"       = set_young(c("log_egfr", "log_holotc")),
  "+ younger nutrition"         = set_young(c("log_egfr", "log_holotc", "phase_angle_deg", "glim_any")))
SCEN_LEV <- names(SCEN)

run_gc <- function(d, target) gcomp_scenarios(d, SCEN, target = target)
gc_boot <- function(target, B = N_BOOT) {
  est <- run_gc(mdat, target)
  bs <- map_dfr(seq_len(B), function(b) {
    i <- unlist(lapply(split(seq_len(nrow(mdat)), mdat$stratum), function(ix) sample(ix, replace = TRUE)))
    run_gc(mdat[i, ], target) %>% mutate(b = b)
  })
  est %>% left_join(bs %>% group_by(scenario) %>%
                      summarise(m_lo = quantile(mean_hcy, 0.025), m_hi = quantile(mean_hcy, 0.975),
                                p_lo = quantile(p_hhcy, 0.025), p_hi = quantile(p_hhcy, 0.975)), by = "scenario") %>%
    mutate(target = target)
}
gc <- bind_rows(gc_boot(STRATA[2]), gc_boot(STRATA[3])) %>%
  mutate(scenario = factor(scenario, rev(SCEN_LEV)), target = factor(target, STRATA[2:3]))
ref_y <- mdat %>% group_by(stratum) %>% summarise(m = mean(hcy), p = mean(hcy > 15))
write.csv(gc, file.path(OUT_TAB, "fig5_gcomputation.csv"), row.names = FALSE)

PAL_SCEN <- c("Observed (model)" = "#3A3A3A", "Younger RBC folate" = "#7FB7A4", "Younger eGFR" = "#2B6CB0",
              "Younger holoTC" = "#9E2A2B", "Younger eGFR + holoTC" = "#6C5B7B", "+ younger nutrition" = "#E09F3E")

# ---- P1: mean Hcy ----------------------------------------------------------------------
p1 <- ggplot(gc, aes(mean_hcy, scenario, colour = scenario)) +
  geom_vline(data = ref_y %>% filter(stratum == STRATA[1]), aes(xintercept = m), linetype = "22",
             linewidth = 0.35, colour = PAL_STR[1]) +
  geom_vline(data = gc %>% filter(scenario == "Observed (model)"), aes(xintercept = mean_hcy),
             linewidth = 0.25, colour = "grey60") +
  geom_pointrange(aes(xmin = m_lo, xmax = m_hi), size = 0.25, linewidth = 0.55) +
  geom_text(aes(x = m_hi, label = sprintf("%.1f", mean_hcy)), hjust = -0.25, size = 1.65, family = FONT, colour = "black") +
  facet_wrap(~target, nrow = 1, labeller = labeller(target = function(x) paste("Target:", x))) +
  scale_colour_manual(values = PAL_SCEN, guide = "none") +
  scale_x_continuous("Mean homocysteine (µmol/L, 95% CI)", expand = expansion(mult = c(0.05, 0.15))) +
  scale_y_discrete(NULL) +
  labs(subtitle = "Dashed: observed mean of the 56–75 y stratum")

# ---- P2: % HHcy -----------------------------------------------------------------------------
p2 <- ggplot(gc, aes(100 * p_hhcy, scenario, fill = scenario)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_errorbar(aes(xmin = 100 * p_lo, xmax = 100 * p_hi), width = 0.25, linewidth = 0.3, orientation = "y") +
  geom_vline(data = ref_y %>% filter(stratum == STRATA[1]), aes(xintercept = 100 * p), linetype = "22",
             linewidth = 0.35, colour = PAL_STR[1]) +
  geom_text(aes(x = 100 * p_hi, label = sprintf("%.0f%%", 100 * p_hhcy)), hjust = -0.25, size = 1.65, family = FONT) +
  facet_wrap(~target, nrow = 1, labeller = labeller(target = function(x) paste("Target:", x))) +
  scale_fill_manual(values = PAL_SCEN, guide = "none") +
  scale_x_continuous("Hcy > 15 µmol/L (%, 95% CI)", expand = expansion(mult = c(0, 0.2))) +
  scale_y_discrete(NULL)

# ---- P3: distribution shift (individual predictions, >=85 y) ----------------------------------
fit_g <- lm(G_FORMULA, data = mdat); smear <- mean(exp(resid(fit_g)))
tg <- mdat %>% filter(stratum == STRATA[3]); rf <- mdat %>% filter(stratum == STRATA[1])
ind <- bind_rows(
  tibble(set = "≥85 y, observed", hcy = exp(predict(fit_g, tg)) * smear),
  tibble(set = "≥85 y, younger eGFR + holoTC", hcy = exp(predict(fit_g, SCEN[["Younger eGFR + holoTC"]](tg, rf))) * smear),
  tibble(set = "56–75 y, observed", hcy = exp(predict(fit_g, rf)) * smear)) %>%
  mutate(set = factor(set, c("56–75 y, observed", "≥85 y, younger eGFR + holoTC", "≥85 y, observed")))
PAL_SET <- setNames(c(PAL_STR[[1]], "#6C5B7B", PAL_STR[[3]]), levels(ind$set))
p3 <- ggplot(ind, aes(hcy, fill = set, colour = set)) +
  geom_vline(xintercept = 15, linetype = "22", linewidth = 0.3, colour = COL_REF) +
  geom_density(alpha = 0.25, linewidth = 0.5, adjust = 1.1) +
  scale_fill_manual(values = PAL_SET, name = NULL) + scale_colour_manual(values = PAL_SET, name = NULL) +
  scale_x_log10("Model-predicted Hcy (µmol/L)", breaks = c(8, 10, 15, 20, 30)) +
  scale_y_continuous("Density", expand = expansion(mult = c(0, 0.05))) +
  theme(legend.position = "inside", legend.position.inside = c(1, 1), legend.justification = c(1, 1),
        legend.key.size = unit(5, "pt"), legend.text = element_text(size = BASE_PT - 1.5))

# ---- P4: dose–response of the gap closed ------------------------------------------------------
fracs <- seq(0, 1, by = 0.1)
dr_once <- function(d) {
  f <- lm(G_FORMULA, data = d); sm <- mean(exp(resid(f)))
  tg <- d %>% filter(stratum == STRATA[3]); rf <- d %>% filter(stratum == STRATA[1])
  map_dfr(fracs, function(fr) {
    shift <- function(vars) { x <- tg; for (v in vars) x[[v]] <- quantile_map(tg[[v]], rf[[v]], frac = fr); x }
    tibble(frac = fr,
           `eGFR only` = mean(exp(predict(f, shift("log_egfr"))) * sm),
           `holoTC only` = mean(exp(predict(f, shift("log_holotc"))) * sm),
           `eGFR + holoTC` = mean(exp(predict(f, shift(c("log_egfr", "log_holotc")))) * sm))
  })
}
dr_est <- dr_once(mdat) %>% pivot_longer(-frac, names_to = "lever", values_to = "m")
dr_bs <- map_dfr(seq_len(max(100, N_BOOT %/% 5)), function(b) {
  i <- unlist(lapply(split(seq_len(nrow(mdat)), mdat$stratum), function(ix) sample(ix, replace = TRUE)))
  dr_once(mdat[i, ]) %>% mutate(b = b)
}) %>% pivot_longer(-c(frac, b), names_to = "lever", values_to = "m") %>%
  group_by(frac, lever) %>% summarise(lo = quantile(m, 0.025), hi = quantile(m, 0.975))
dr <- dr_est %>% left_join(dr_bs, by = c("frac", "lever")) %>%
  mutate(lever = factor(lever, c("holoTC only", "eGFR only", "eGFR + holoTC")))
PAL_LEV <- c("eGFR only" = "#2B6CB0", "holoTC only" = "#9E2A2B", "eGFR + holoTC" = "#6C5B7B")
p4 <- ggplot(dr, aes(100 * frac, m, colour = lever, fill = lever)) +
  geom_hline(yintercept = ref_y$m[ref_y$stratum == STRATA[1]], linetype = "22", linewidth = 0.35, colour = PAL_STR[1]) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = PAL_LEV, name = NULL) + scale_fill_manual(values = PAL_LEV, name = NULL) +
  scale_x_continuous("Gap to the 56–75 y distribution closed (%)", breaks = seq(0, 100, 25)) +
  scale_y_continuous("Mean Hcy in ≥85 y (µmol/L)") +
  theme(legend.position = "top", legend.key.size = unit(5, "pt"))

# ---- P5: interventional mediation (waterfall) -----------------------------------------------------
MED <- list(Kidney = "log_egfr", `Active B12` = "log_holotc", Nutrition = c("phase_angle_deg", "glim_any"))
med_once <- function(d, target) {
  f <- lm(G_FORMULA, data = d); sm <- mean(exp(resid(f)))
  tg <- d %>% filter(stratum == target); rf <- d %>% filter(stratum == STRATA[1])
  ey <- function(x) mean(exp(predict(f, x)) * sm)
  y_obs <- ey(tg)                                             # Y(a = old, M ~ old)
  y_allm <- ey(set_young(unlist(MED))(tg, rf))                 # Y(a = old, all M ~ young)
  y_ref <- { x <- set_young(unlist(MED))(tg, rf); x$stratum <- factor(STRATA[1], STRATA); ey(x) }  # Y(a = young, M ~ young)
  ie <- map_dbl(MED, ~ y_obs - ey(set_young(.x)(tg, rf)))
  c(ie, `Mediator overlap` = (y_obs - y_allm) - sum(ie), `Direct (age)` = y_allm - y_ref, Total = y_obs - y_ref)
}
med_res <- map_dfr(STRATA[2:3], function(tgt) {
  est <- med_once(mdat, tgt)
  bs <- replicate(max(200, N_BOOT %/% 2), {
    i <- unlist(lapply(split(seq_len(nrow(mdat)), mdat$stratum), function(ix) sample(ix, replace = TRUE)))
    med_once(mdat[i, ], tgt) })
  tibble(comp = names(est), est = est, lo = apply(bs, 1, quantile, 0.025), hi = apply(bs, 1, quantile, 0.975), target = tgt)
}) %>% group_by(target) %>% mutate(pm = est / est[comp == "Total"]) %>% ungroup() %>%
  mutate(target = factor(target, STRATA[2:3]))
write.csv(med_res, file.path(OUT_TAB, "fig5_interventional_mediation.csv"), row.names = FALSE)
wf <- med_res %>% filter(comp != "Total") %>%
  mutate(comp = factor(comp, c("Kidney", "Active B12", "Nutrition", "Mediator overlap", "Direct (age)"))) %>%
  arrange(target, comp) %>% group_by(target) %>%
  mutate(end = cumsum(est), start = lag(end, default = 0), id = as.integer(comp)) %>% ungroup()
tot <- med_res %>% filter(comp == "Total") %>% mutate(id = 6L)
PAL_MED <- c("Kidney" = "#2B6CB0", "Active B12" = "#9E2A2B", "Nutrition" = "#E09F3E",
             "Mediator overlap" = "#A8A8A8", "Direct (age)" = "#3A3A3A", "Total" = "black")
p5 <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_rect(data = wf, aes(xmin = id - 0.36, xmax = id + 0.36, ymin = start, ymax = end, fill = comp)) +
  geom_segment(data = wf %>% group_by(target) %>% filter(id < max(id)),
               aes(x = id + 0.36, xend = id + 0.64, y = end, yend = end), linewidth = 0.2, colour = "grey50") +
  geom_rect(data = tot, aes(xmin = id - 0.36, xmax = id + 0.36, ymin = 0, ymax = est), fill = "white",
            colour = "black", linewidth = 0.35) +
  geom_errorbar(data = tot, aes(x = id, ymin = lo, ymax = hi), width = 0.15, linewidth = 0.3) +
  geom_text(data = wf, aes(id, pmax(start, end), label = sprintf("%.0f%%", 100 * pm)), vjust = -0.4,
            size = 1.6, family = FONT) +
  geom_text(data = tot, aes(id, hi, label = sprintf("+%.1f", est)), vjust = -0.5, size = 1.65, family = FONT) +
  facet_wrap(~target, nrow = 1, labeller = labeller(target = function(x) paste(x, "vs 56–75 y"))) +
  scale_fill_manual(values = PAL_MED, guide = "none") +
  scale_x_continuous(NULL, breaks = 1:6, labels = c("Kidney", "Active\nB12", "Nutrition", "Overlap", "Direct\n(age)", "Total")) +
  scale_y_continuous("Hcy difference (µmol/L)", expand = expansion(mult = c(0.02, 0.12))) +
  theme(axis.text.x = element_text(size = BASE_PT - 1.5, lineheight = 0.85))

# ---- P6: robustness across sensitivity analyses --------------------------------------------------------
sens_sets <- list(
  "Main analysis"                           = mdat,
  "Excl. AKI at admission"                  = mdat %>% filter(acute_kidney_injury == 0),
  "Excl. B12 supplement / antifolate drugs" = mdat %>% filter(sens_exclude_b12_antifolate == 0),
  "Excl. levodopa users"                    = mdat %>% filter(med_levodopa == 0),
  "Excl. infection at admission"            = mdat %>% filter(infection_at_admission == 0),
  "Reliable FFQ only"                       = mdat %>% filter(ffq_reliable == 1),
  "eGFR by creatinine (CKD-EPI 2021)"       = mdat %>% mutate(log_egfr = log(egfr_ckdepi2021)),
  "Women"                                   = mdat %>% filter(male == 0),
  "Men"                                     = mdat %>% filter(male == 1))
# Main-effects (linear) g-formula for the sensitivity sets: stable in small subsets, where the
# stratum x spline terms of G_FORMULA extrapolate poorly (e.g. only ~35 participants >=85 y with reliable FFQ)
G_LIN <- update(G_FORMULA, . ~ stratum + log_egfr + log_holotc + log_rbcf + log_plp + egrac + glim_any + phase_angle_deg +
                  log_crp + male + ancestry_pc1 + ancestry_pc2 + time_to_centrifugation_min + fasting_hours +
                  diabetes_t2 + hcy_genetic_risk_score)
gap_closed <- function(d) {
  f <- if (n_distinct(d$male) == 1) update(G_LIN, . ~ . - male) else G_LIN
  r <- gcomp_scenarios(d, SCEN[c("Observed (model)", "Younger eGFR + holoTC")], formula = f)
  # (ratio "share of gap closed" is unstable in small subsets; absolute reduction is reported)
  r$mean_hcy[1] - r$mean_hcy[2]                       # absolute reduction in mean Hcy (umol/L) in >=85 y
}
sens <- imap_dfr(sens_sets, function(d, nm) {
  est <- gap_closed(d)
  bs <- replicate(max(100, N_BOOT %/% 5), {
    i <- unlist(lapply(split(seq_len(nrow(d)), d$stratum), function(ix) sample(ix, replace = TRUE)))
    gap_closed(d[i, ]) })
  tibble(set = nm, n = nrow(d), est = est, lo = quantile(bs, 0.025, na.rm = TRUE), hi = quantile(bs, 0.975, na.rm = TRUE))
}) %>% mutate(set = factor(set, rev(names(sens_sets))))
write.csv(sens, file.path(OUT_TAB, "fig5_sensitivity_gap_closed.csv"), row.names = FALSE)
p6 <- ggplot(sens, aes(est, set)) +
  geom_vline(xintercept = sens$est[1], linetype = "22", linewidth = 0.3, colour = COL_REF) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi, colour = set == "Main analysis"), size = 0.22, linewidth = 0.45) +
  geom_text(aes(x = Inf, label = paste0("n = ", n)), hjust = 1.05, size = 1.55, family = FONT, colour = "grey35") +
  scale_colour_manual(values = c(`TRUE` = "#6C5B7B", `FALSE` = "black"), guide = "none") +
  scale_x_continuous("Reduction in mean Hcy of ≥85 y with younger eGFR + holoTC\n(µmol/L, 95% CI; main-effects g-formula)",
                     expand = expansion(mult = c(0.05, 0.25))) +
  coord_cartesian() +
  scale_y_discrete(NULL)

# ---- Assemble ---------------------------------------------------------------------------------------------
r1 <- row_of(p1, p2, widths = c(1.15, 1))
r2 <- row_of(p4, p5, widths = c(0.8, 1.35))
r3 <- row_of(p3, p6, widths = c(0.8, 1.2))
fig5 <- wrap_elements(full = r1) / wrap_elements(full = r2) / wrap_elements(full = r3) +
  plot_layout(heights = c(0.85, 1, 0.85))
save_fig(fig5, "Figure5_counterfactual", width = 180, height = 175)
save_panels(list(gcomp_mean = p1, gcomp_prev = p2, distribution = p3, dose_response = p4,
                 mediation = p5, robustness = p6), "Figure5")
