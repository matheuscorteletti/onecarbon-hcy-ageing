# =============================================================================
# Figure 4 — What explains the higher homocysteine of the oldest old? (Q3; H1a, H1b->folate paradox)
#   P1  Homocysteine by group (raincloud) with % > 15 umol/L
#   P2  Age gradient of the determinants (median, % of the 56–75 y hospital stratum)
#   P3  Attenuation: % difference in Hcy vs 56–75 y across hierarchical models M1–M5
#   P4  Oaxaca–Blinder decomposition of the >=85 vs 56–75 y gap (bootstrap 95% CI)
#   P5  Relative importance (LMG, R2 partition) with bootstrap 95% CI
#   P6–P8  Adjusted GAM curves by stratum: eGFR-cys, holoTC, RBC folate (+ crude folate curve)
#   P9  Folate paradox: crude vs adjusted (partial) Spearman rho with Hcy
#   P10 Joint surface: predicted Hcy over eGFR-cys x holoTC, with stratum centroids
#   P11 MMA vs eGFR-cys: why holoTC (not MMA) indexes B12 in the old
#   P12 Quantile regression: >=85 y effect across the Hcy distribution, crude vs adjusted
# =============================================================================
source("R/00_setup.R"); source("R/01_models.R")
suppressPackageStartupMessages({ library(ggdist); library(mgcv); library(relaimpo); library(quantreg) })
select <- dplyr::select

# ---- P1: Hcy by group -------------------------------------------------------------
d1 <- dat %>% filter(!is.na(hcy))
lab1 <- d1 %>% group_by(grp4) %>% summarise(p = mean(hhcy), y = max(hcy))
kw1 <- kruskal.test(hcy ~ grp4, data = d1)$p.value
p1 <- ggplot(d1, aes(grp4, hcy, fill = grp4, colour = grp4)) +
  geom_hline(yintercept = 15, linetype = "22", linewidth = 0.3, colour = COL_REF) +
  stat_halfeye(adjust = 0.9, width = 0.6, justification = -0.3, .width = 0, point_colour = NA,
               alpha = 0.75, slab_linewidth = 0) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_point(position = position_jitter(width = 0.06, height = 0, seed = 1), size = 0.25, alpha = 0.35, stroke = 0) +
  geom_text(data = lab1, aes(grp4, 42, label = sprintf("%.0f%%\n>15", 100 * p)), colour = "black",
            size = 1.7, family = FONT, lineheight = 0.85) +
  scale_y_log10(lab_umol, breaks = c(5, 10, 15, 20, 30)) +
  coord_cartesian(ylim = c(4.5, 48)) +
  scale_fill_manual(values = PAL_GRP4, guide = "none") +
  scale_colour_manual(values = colorspace::darken(PAL_GRP4, 0.25), guide = "none") +
  scale_x_discrete(NULL, labels = lab_grp4) + rot_x +
  labs(subtitle = paste("Kruskal–Wallis", fmt_p(kw1)))

# ---- P2: determinant gradient -----------------------------------------------------
grad_vars <- c("eGFR-cys" = "egfr", "holoTC" = "holotc", "RBC folate" = "rbc_fol",
               "Phase angle" = "phase_angle_deg", "MMA" = "mma_umol_l", "Homocysteine" = "hcy")
grad <- hosp %>% select(stratum, all_of(unname(grad_vars))) %>% pivot_longer(-stratum) %>%
  group_by(name, stratum) %>%
  summarise(m = median(value, na.rm = TRUE), lo = quantile(value, 0.25, na.rm = TRUE), hi = quantile(value, 0.75, na.rm = TRUE)) %>%
  group_by(name) %>% mutate(across(c(m, lo, hi), ~ 100 * .x / m[stratum == STRATA[1]])) %>% ungroup() %>%
  mutate(var = factor(names(grad_vars)[match(name, grad_vars)], names(grad_vars)))
PAL_GRAD <- c("eGFR-cys" = "#2B6CB0", "holoTC" = "#9E2A2B", "RBC folate" = "#4C956C",
              "Phase angle" = "#E09F3E", "MMA" = "#6C5B7B", "Homocysteine" = "black")
p2 <- ggplot(grad, aes(stratum, m, colour = var, group = var)) +
  geom_hline(yintercept = 100, linewidth = 0.3, colour = COL_REF) +
  geom_line(aes(linewidth = var == "Homocysteine"), position = position_dodge(0.25)) +
  geom_point(size = 1.1, position = position_dodge(0.25)) +
  ggrepel::geom_text_repel(data = filter(grad, stratum == STRATA[3]), aes(label = var), hjust = 0,
                           nudge_x = 0.25, direction = "y", size = 1.75, family = FONT, segment.size = 0.15,
                           min.segment.length = 0, box.padding = 0.08, seed = 1) +
  scale_linewidth_manual(values = c(`FALSE` = 0.45, `TRUE` = 0.9), guide = "none") +
  scale_colour_manual(values = PAL_GRAD, guide = "none") +
  scale_x_discrete("Hospital stratum", expand = expansion(add = c(0.25, 1.05))) +
  scale_y_continuous("Median (% of 56–75 y)", breaks = seq(40, 180, 20))

# ---- P3: attenuation ------------------------------------------------------------------
att <- imap_dfr(hier_models, function(m, nm) {
  ci <- confint(m)
  tibble(model = nm, term = c("76–84 y", "≥85 y"), b = coef(m)[c("s76", "s85")],
         lo = ci[c("s76", "s85"), 1], hi = ci[c("s76", "s85"), 2], r2 = summary(m)$r.squared)
}) %>% mutate(model = factor(model, names(BLOCKS)), across(c(b, lo, hi), pct_change))
att_lab <- att %>% filter(term == "≥85 y") %>%
  mutate(red = 100 * (1 - b / b[1]), lab = sprintf("R² %.2f", r2))
p3 <- ggplot(att, aes(model, b, colour = term, group = term)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = COL_REF) +
  geom_line(position = position_dodge(0.35), linewidth = 0.4) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.35), size = 0.2, linewidth = 0.45) +
  geom_text(data = att_lab, aes(model, -14, label = lab), inherit.aes = FALSE, size = 1.6, family = FONT, colour = "grey30") +
  scale_colour_manual(values = PAL_STR[2:3], name = "vs 56–75 y") +
  scale_x_discrete(NULL, labels = c("M1\nage, design", "M2\n+ eGFR", "M3\n+ holoTC", "M4\n+ folate,\nB6, B2",
                                    "M5\n+ nutrition,\nCRP, DM, GRS")) +
  scale_y_continuous("Difference in Hcy (%, 95% CI)") +
  theme(legend.position = "inside", legend.position.inside = c(1, 1), legend.justification = c(1, 1),
        axis.text.x = element_text(size = BASE_PT - 1.5, lineheight = 0.85))

# ---- P4: Oaxaca–Blinder (pooled twofold), bootstrap by group ------------------------------
od <- mdat %>% filter(stratum != STRATA[2]) %>% mutate(g = as.integer(stratum == STRATA[3]))
xvars <- unlist(DET_GROUPS)
ob_once <- function(d) {
  fit <- lm(reformulate(c(xvars, "g"), "log_hcy"), data = d)          # pooled with group dummy (Jann 2008)
  dx <- colMeans(d[d$g == 1, xvars]) - colMeans(d[d$g == 0, xvars])
  contrib <- dx * replace_na(coef(fit)[xvars], 0)
  grp <- map_dbl(DET_GROUPS, ~ sum(contrib[.x]))
  tot <- mean(d$log_hcy[d$g == 1]) - mean(d$log_hcy[d$g == 0])
  c(grp, "Age stratum (unexplained)" = tot - sum(grp), total = tot)
}
ob_est <- ob_once(od)
ob_bs  <- replicate(N_BOOT, ob_once(od[sample.int(nrow(od), replace = TRUE), ]))
ob <- tibble(comp = names(ob_est), est = ob_est, lo = apply(ob_bs, 1, quantile, 0.025),
             hi = apply(ob_bs, 1, quantile, 0.975)) %>%
  filter(comp != "total") %>% mutate(share = est / ob_est["total"]) %>%
  arrange(desc(est)) %>% mutate(comp = factor(comp, rev(comp)))
p4 <- ggplot(ob, aes(100 * est, comp, fill = comp)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_col(width = 0.68) +
  geom_errorbar(aes(xmin = 100 * lo, xmax = 100 * hi), width = 0.25, linewidth = 0.3, orientation = "y") +
  geom_text(aes(x = 100 * pmax(hi, 0) + 0.8, label = sprintf("%.0f%%", 100 * share)), hjust = 0, size = 1.65, family = FONT) +
  scale_fill_manual(values = PAL_DET, guide = "none") +
  scale_x_continuous("Contribution (log-Hcy × 100)",
                     expand = expansion(mult = c(0.05, 0.18))) +
  scale_y_discrete(NULL) +
  labs(subtitle = sprintf("Total gap: +%.0f%% Hcy", pct_change(ob_est["total"])))

# ---- P5: LMG relative importance ------------------------------------------------------------
X <- c("s76", "s85", xvars)
lmf <- lm(reformulate(X, "log_hcy"), data = mdat)
grp_ri <- c(list("Age stratum" = c("s76", "s85")), DET_GROUPS)
ri_once <- function(d) {
  f <- lm(reformulate(X, "log_hcy"), data = d)
  r <- calc.relimp(f, type = "lmg", groups = unname(grp_ri), groupnames = names(grp_ri))@lmg
  r[names(grp_ri)[lengths(grp_ri) == 1]] <- r[unlist(grp_ri[lengths(grp_ri) == 1])]
  r[names(grp_ri)]
}
ri_est <- ri_once(mdat)
ri_bs  <- replicate(max(200, N_BOOT %/% 5), ri_once(mdat[sample.int(nrow(mdat), replace = TRUE), ]))
ri <- tibble(g = names(ri_est), est = ri_est, lo = apply(ri_bs, 1, quantile, 0.025), hi = apply(ri_bs, 1, quantile, 0.975)) %>%
  arrange(est) %>% mutate(g = factor(g, g))
PAL_RI <- c(PAL_DET, "Age stratum" = "#3A3A3A")
p5 <- ggplot(ri, aes(100 * est, g, colour = g)) +
  geom_segment(aes(x = 0, xend = 100 * est, yend = g), linewidth = 0.9, alpha = 0.35) +
  geom_pointrange(aes(xmin = 100 * lo, xmax = 100 * hi), size = 0.22, linewidth = 0.45) +
  scale_colour_manual(values = PAL_RI, guide = "none") +
  scale_x_continuous("Share of R² (%)", expand = expansion(mult = c(0, 0.05))) +
  scale_y_discrete(NULL) +
  labs(subtitle = sprintf("LMG; model R² = %.2f", summary(lmf)$r.squared))

# ---- P6–P8: GAM curves by stratum ---------------------------------------------------------------
gam_base <- "s(log_rbcf, k = 5) + s(log_holotc, k = 5) + s(log_egfr, k = 5)"
gam_curve <- function(var, xlab, xbreaks, crude = FALSE) {
  others <- setdiff(c("log_egfr", "log_holotc", "log_rbcf"), var)
  f <- as.formula(paste("log_hcy ~ stratum + s(", var, ", by = stratum, k = 5) +",
                        paste0("s(", others, ", k = 5)", collapse = " + "),
                        "+ log_plp + egrac + glim_any + phase_angle_deg + log_crp + male + ancestry_pc1 + ancestry_pc2 + time_to_centrifugation_min + fasting_hours"))
  g <- gam(f, data = mdat, method = "REML")
  nd <- mdat %>% group_by(stratum) %>%
    reframe(x = seq(quantile(.data[[var]], 0.03), quantile(.data[[var]], 0.97), length.out = 80)) %>%
    rename(!!var := x)
  meds <- mdat %>% summarise(across(c(log_egfr, log_holotc, log_rbcf, log_plp, egrac, phase_angle_deg, log_crp,
                                      ancestry_pc1, ancestry_pc2, time_to_centrifugation_min, fasting_hours), median))
  nd <- bind_cols(nd, meds %>% select(-all_of(var)) %>% slice(rep(1, nrow(nd))))
  nd$glim_any <- 0; nd$male <- 0.5
  pr <- predict(g, nd, se.fit = TRUE)
  nd <- nd %>% mutate(fit = exp(pr$fit), lo = exp(pr$fit - 1.96 * pr$se.fit), hi = exp(pr$fit + 1.96 * pr$se.fit),
                      xx = exp(.data[[var]]))
  st <- summary(g)$s.table
  ptxt <- paste(sprintf("%s edf %.1f", str_remove(STRATA, " y$"), st[grep(var, rownames(st))[1:3], "edf"]), collapse = "; ")
  p <- ggplot(nd, aes(xx, fit, colour = stratum, fill = stratum)) +
    geom_hline(yintercept = 15, linetype = "22", linewidth = 0.3, colour = COL_REF) +
    geom_rug(data = mdat, aes(x = exp(.data[[var]]), colour = stratum), inherit.aes = FALSE, sides = "b",
             alpha = 0.35, linewidth = 0.15, length = unit(2, "pt")) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = PAL_STR, name = NULL) + scale_fill_manual(values = PAL_STR, name = NULL) +
    scale_x_log10(xlab, breaks = xbreaks) +
    scale_y_continuous("Adjusted Hcy (µmol/L)") +
    theme(legend.position = "none")
  if (crude) {
    cr <- gam(log_hcy ~ s(log_rbcf, k = 5), data = mdat)
    cnd <- tibble(log_rbcf = seq(quantile(mdat$log_rbcf, 0.03), quantile(mdat$log_rbcf, 0.97), length.out = 80))
    cnd$fit <- exp(predict(cr, cnd))
    p <- p + geom_line(data = cnd, aes(exp(log_rbcf), fit), inherit.aes = FALSE, linetype = "42",
                       linewidth = 0.5, colour = "black") +
      annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.2, size = 1.65, family = FONT,
               label = "dashed = crude (unadjusted)")
  }
  p
}
p6 <- gam_curve("log_egfr",   "eGFR-cys (mL/min/1.73 m²)", c(15, 30, 60, 90, 120)) +
  theme(legend.position = "inside", legend.position.inside = c(1, 1), legend.justification = c(1, 1),
        legend.key.size = unit(5, "pt"))
p7 <- gam_curve("log_holotc", "holoTC (pmol/L)", c(15, 25, 35, 50, 75, 100, 150))
p8 <- gam_curve("log_rbcf",   "RBC folate (ng/mL)", c(300, 400, 600, 800, 1000), crude = TRUE)

# ---- P9: crude vs adjusted correlation (folate paradox) -------------------------------------------
pc_vars <- c("eGFR-cys" = "log_egfr", "holoTC" = "log_holotc", "RBC folate" = "log_rbcf", "Serum folate" = "folate_serum_ng_ml",
             "PLP" = "log_plp", "EGRAC" = "egrac", "Phase angle" = "phase_angle_deg", "MNA-SF" = "mna_sf",
             "GLIM (any)" = "glim_any", "CRP" = "log_crp", "GRS" = "hcy_genetic_risk_score")
pdx <- mdat %>% left_join(hosp %>% select(participant_id, folate_serum_ng_ml, mna_sf), by = "participant_id")
adj_set <- c("s76", "s85", "male", "log_egfr", "log_holotc")
pcor <- map_dfr(names(pc_vars), function(nm) {
  v <- pc_vars[[nm]]
  Z <- pdx[, setdiff(adj_set, v)]
  ra <- cor(rank_resid(pdx[[v]], Z), rank_resid(pdx$log_hcy, Z), use = "complete.obs")
  tibble(var = nm, crude = cor(pdx[[v]], pdx$log_hcy, method = "spearman"), adjusted = ra)
}) %>% mutate(var = fct_reorder(var, adjusted))
n_pc <- nrow(pdx); crit <- 1.96 / sqrt(n_pc - 3)
p9 <- ggplot(pcor, aes(y = var)) +
  annotate("rect", xmin = -crit, xmax = crit, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_segment(aes(x = crude, xend = adjusted, yend = var), linewidth = 0.4, colour = "grey55",
               arrow = arrow(length = unit(2.2, "pt"), type = "closed")) +
  geom_point(aes(x = crude), shape = 21, fill = "white", colour = "black", size = 1.4, stroke = 0.35) +
  geom_point(aes(x = adjusted, fill = adjusted), shape = 21, colour = "black", size = 1.6, stroke = 0.35) +
  scale_fill_gradient2(low = PAL_DIV["low"], mid = "white", high = PAL_DIV["high"], limits = c(-0.7, 0.7), guide = "none") +
  scale_x_continuous("Spearman ρ with Hcy",
                     limits = c(-0.75, 0.6)) +
  scale_y_discrete(NULL) +
  labs(subtitle = "○ crude  ● partial (age, sex, eGFR, holoTC)")

# ---- P10: joint surface eGFR x holoTC ------------------------------------------------------------
gs <- gam(log_hcy ~ stratum + te(log_egfr, log_holotc, k = c(5, 5)) + log_rbcf + log_plp + egrac + glim_any +
            phase_angle_deg + log_crp + male + ancestry_pc1 + ancestry_pc2 + time_to_centrifugation_min + fasting_hours,
          data = mdat, method = "REML")
grid <- expand_grid(log_egfr = seq(log(12), log(130), length.out = 90), log_holotc = seq(log(12), log(180), length.out = 90)) %>%
  bind_cols(mdat %>% summarise(across(c(log_rbcf, log_plp, egrac, phase_angle_deg, log_crp, ancestry_pc1, ancestry_pc2,
                                        time_to_centrifugation_min, fasting_hours), median))) %>%
  mutate(stratum = factor(STRATA[1], STRATA), glim_any = 0, male = 0.5)
grid$fit <- exp(predict(gs, grid))
# mask grid cells far from data (convex-hull-like using 2D density)
kd <- MASS::kde2d(mdat$log_egfr, mdat$log_holotc, n = 90, lims = c(range(grid$log_egfr), range(grid$log_holotc)))
grid$dens <- as.vector(kd$z)
grid <- grid %>% mutate(fit = if_else(dens < max(dens) * 0.02, NA_real_, fit))
cent <- mdat %>% group_by(stratum) %>% summarise(e = exp(median(log_egfr)), h = exp(median(log_holotc)))
p10 <- ggplot(grid, aes(exp(log_egfr), exp(log_holotc))) +
  geom_raster(aes(fill = fit), interpolate = TRUE) +
  geom_contour(aes(z = fit), breaks = c(10, 12.5, 15, 20, 25), colour = "white", linewidth = 0.3, na.rm = TRUE) +
  geom_point(data = mdat, aes(egfr, holotc), size = 0.2, alpha = 0.25, stroke = 0, colour = "black") +
  geom_path(data = cent, aes(e, h), arrow = arrow(length = unit(3, "pt"), type = "closed"), linewidth = 0.5) +
  geom_point(data = cent, aes(e, h, colour = stratum), size = 1.8) +
  geom_point(data = cent, aes(e, h), size = 1.8, shape = 21, colour = "black", stroke = 0.4) +
  scale_colour_manual(values = PAL_STR, guide = "none") +
  scale_fill_gradientn(colours = PAL_SEQ, limits = c(7, 28), oob = squish, na.value = "white",
                       name = "Predicted\nHcy", breaks = c(10, 15, 20, 25)) +
  scale_x_log10("eGFR-cys (mL/min/1.73 m²)", breaks = c(15, 30, 60, 120), expand = c(0, 0)) +
  scale_y_log10("holoTC (pmol/L)", breaks = c(15, 35, 70, 150), expand = c(0, 0)) +
  theme(legend.key.height = unit(12, "pt"), legend.key.width = unit(4, "pt"))

# ---- P11: MMA vs eGFR -------------------------------------------------------------------------------
p11df <- hosp %>% mutate(b12 = if_else(holotc < 35, "holoTC < 35", "holoTC ≥ 35"))
p11 <- ggplot(p11df, aes(egfr, mma_umol_l, colour = b12)) +
  geom_hline(yintercept = 0.37, linetype = "22", linewidth = 0.3, colour = COL_REF) +
  geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.9, se = FALSE, linewidth = 0.6) +
  scale_colour_manual(values = c("holoTC < 35" = "#9E2A2B", "holoTC ≥ 35" = "#8FA9C4"), name = NULL) +
  scale_x_log10("eGFR-cys (mL/min/1.73 m²)", breaks = c(15, 30, 60, 120)) +
  scale_y_log10("MMA (µmol/L)") +
  annotate("text", x = Inf, y = 0, hjust = 1.02, vjust = -0.5, size = 1.65, family = FONT,
           label = spearman_lab(p11df$egfr, p11df$mma_umol_l)) +
  theme(legend.position = "inside", legend.position.inside = c(1, 1), legend.justification = c(1, 1),
        legend.key.size = unit(5, "pt"))

# ---- P12: quantile regression ----------------------------------------------------------------------
taus <- seq(0.1, 0.9, by = 0.1)
qr_fit <- function(form, lab) map_dfr(taus, function(t) {
  s <- summary(rq(form, tau = t, data = mdat), se = "boot", R = max(100, N_BOOT %/% 5))$coefficients
  tibble(tau = t, b = s["s85", 1], se = s["s85", 2], model = lab)
})
qr <- bind_rows(qr_fit(block_formula(1), "M1: age + design"),
                qr_fit(FULL_FORMULA, "M5: fully adjusted")) %>%
  mutate(across(c(b), list(lo = ~ .x - 1.96 * se, hi = ~ .x + 1.96 * se)), across(c(b, b_lo, b_hi), pct_change))
p12 <- ggplot(qr, aes(tau, b, colour = model, fill = model)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_ribbon(aes(ymin = b_lo, ymax = b_hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.55) + geom_point(size = 0.8) +
  scale_colour_manual(values = c("#9E2A2B", "#3A3A3A"), name = NULL) +
  scale_fill_manual(values = c("#9E2A2B", "#3A3A3A"), name = NULL) +
  scale_x_continuous("Quantile of Hcy distribution (τ)", breaks = taus[c(TRUE, FALSE)]) +
  scale_y_continuous("≥85 vs 56–75 y (% difference)") +
  theme(legend.position = "top", legend.key.size = unit(5, "pt"))

# ---- Assemble --------------------------------------------------------------------------------------------
r1 <- row_of(p1, p2, p3, widths = c(0.9, 1, 1.25))
r2 <- row_of(p4, p5, p9, widths = c(1.05, 0.95, 1))
r3 <- row_of(p6, p7, p8)
r4 <- row_of(p10, p11, p12, widths = c(1.1, 0.95, 0.95))
fig4 <- wrap_elements(full = r1) / wrap_elements(full = r2) / wrap_elements(full = r3) / wrap_elements(full = r4) +
  plot_layout(heights = c(1, 1, 0.82, 0.9))
save_fig(fig4, "Figure4_determinants_hcy", width = 180, height = 235)
save_panels(list(hcy_groups = p1, gradient = p2, attenuation = p3, oaxaca = p4, lmg = p5,
                 gam_egfr = p6, gam_holotc = p7, gam_folate = p8, folate_paradox = p9,
                 surface = p10, mma_egfr = p11, quantreg = p12), "Figure4")

# Export numeric results for the text
write.csv(ob, file.path(OUT_TAB, "fig4_oaxaca_decomposition.csv"), row.names = FALSE)
write.csv(ri, file.path(OUT_TAB, "fig4_lmg_relative_importance.csv"), row.names = FALSE)
write.csv(pcor, file.path(OUT_TAB, "fig4_crude_vs_partial_correlations.csv"), row.names = FALSE)
write.csv(qr, file.path(OUT_TAB, "fig4_quantile_regression.csv"), row.names = FALSE)
