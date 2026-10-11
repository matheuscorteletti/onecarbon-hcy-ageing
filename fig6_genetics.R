# =============================================================================
# Figure 6 — One-carbon genetics vs modifiable determinants (Q5; H1d) — gene x nutrient set
#   P1  Variance partition of log-Hcy by block (LMG) within each group
#   P2  Per-SD effect on Hcy: genetic score vs kidney, B12, folate, phase angle, by group
#   P3  Forest: per-allele effect of the 24 one-carbon variants on Hcy (FDR); SHMT1 = null control
#   P4  Variant x biomarker heatmap (per-allele SD change; FDR stars) — functional validation
#   P5  MTHFR 677 genotype x RBC folate (slopes by genotype)
#   P6  MTHFR 677 genotype x riboflavin status (EGRAC)
#   P7  Joint categories: genetic risk score tertile x low holoTC
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages(library(relaimpo))
select <- dplyr::select

g <- gxn %>% mutate(male = as.integer(sex == "Male"), community = as.integer(!hospital),
                    tt = as.integer(mthfr677 == "TT")) %>%
  filter(!is.na(log_hcy))
COV_G <- c("age_c10", "male", "community", "ancestry_pc1", "ancestry_pc2", "time_to_centrifugation_min")

# ---- P1: variance partition by group ------------------------------------------------------
BLK <- list("Genetics (GRS + 677TT)" = c("hcy_genetic_risk_score", "tt"),
            "Kidney" = "log_egfr", "Active B12" = "log_holotc", "Folate" = "log_rbcf", "B6 + B2" = c("log_plp", "egrac"),
            "Nutrition (GLIM, PhA)" = c("glim_any", "phase_angle_deg"), "Inflammation" = "log_crp",
            "Age, sex, pre-analytics" = c("age_years", "male", "time_to_centrifugation_min"))
PAL_BLK <- c("Genetics (GRS + 677TT)" = "#4C956C", "Kidney" = "#2B6CB0", "Active B12" = "#9E2A2B", "Folate" = "#7FB7A4",
             "B6 + B2" = "#B8A04A", "Nutrition (GLIM, PhA)" = "#E09F3E", "Inflammation" = "#6C5B7B",
             "Age, sex, pre-analytics" = "#C8C8C8", "Unexplained" = "#F1F1F1")
vp <- g %>% drop_na(all_of(unname(unlist(BLK)))) %>% group_by(grp4) %>% group_modify(function(d, k) {
  f <- lm(reformulate(unname(unlist(BLK)), "log_hcy"), data = d)
  r <- calc.relimp(f, type = "lmg", groups = unname(BLK), groupnames = names(BLK))@lmg
  r[names(BLK)[lengths(BLK) == 1]] <- r[unlist(BLK[lengths(BLK) == 1])]
  tibble(block = c(names(BLK), "Unexplained"), v = c(r[names(BLK)], 1 - summary(f)$r.squared), n = nrow(d))
}) %>% ungroup() %>% mutate(block = factor(block, rev(names(PAL_BLK))))
gen_share <- vp %>% filter(block == names(BLK)[1]) %>% select(grp4, gv = v)
p1 <- ggplot(vp, aes(grp4, 100 * v, fill = block)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  geom_text(aes(label = if_else(v > 0.045 & block != "Unexplained", sprintf("%.0f", 100 * v), "")),
            position = position_stack(vjust = 0.5), size = 1.6, family = FONT, colour = "white") +
  geom_text(data = vp %>% distinct(grp4, n), aes(grp4, 101, label = paste0("n=", n)), inherit.aes = FALSE,
            vjust = 0, size = 1.55, family = FONT, colour = "grey30") +
  scale_fill_manual(values = PAL_BLK, name = NULL, breaks = names(PAL_BLK)) +
  scale_y_continuous("Variance in log-Hcy (%)", expand = expansion(mult = c(0, 0.06)), breaks = seq(0, 100, 25)) +
  scale_x_discrete(NULL, labels = lab_grp4) + rot_x +
  theme(legend.position = "right", legend.key.size = unit(5, "pt"), legend.text = element_text(size = BASE_PT - 1.5))

# ---- P2: per-SD standardized effects by group -------------------------------------------------
SDV <- c("Genetic risk score" = "hcy_genetic_risk_score", "eGFR-cys (log)" = "log_egfr", "holoTC (log)" = "log_holotc",
         "RBC folate (log)" = "log_rbcf", "Phase angle" = "phase_angle_deg")
sd_eff <- g %>% drop_na(all_of(unname(SDV)), age_years) %>% group_by(grp4) %>% group_modify(function(d, k) {
  dz <- d %>% mutate(across(all_of(unname(SDV)), ~ as.numeric(scale(.x))))
  f <- lm(reformulate(c(unname(SDV), "age_years", "male", "time_to_centrifugation_min"), "log_hcy"), data = dz)
  ci <- confint(f)
  tibble(var = names(SDV), b = coef(f)[SDV], lo = ci[SDV, 1], hi = ci[SDV, 2])
}) %>% ungroup() %>% mutate(across(c(b, lo, hi), pct_change), var = factor(var, rev(names(SDV))))
p2 <- ggplot(sd_eff, aes(b, var, colour = grp4)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.7), size = 0.15, linewidth = 0.4) +
  scale_colour_manual(values = PAL_GRP4, name = NULL, labels = lab_grp4) +
  labs(x = "% change in Hcy per SD (95% CI)", y = NULL) +
  theme(legend.position = "bottom", legend.key.size = unit(4, "pt"), legend.text = element_text(size = BASE_PT - 1.5)) +
  guides(colour = guide_legend(nrow = 2))

# ---- P3: forest of the 24 variants ---------------------------------------------------------------
snp_assoc <- function(outcome, d = g, scale_y = FALSE) {
  map_dfr(seq_len(nrow(snp_tab)), function(i) {
    v <- snp_tab$var[i]
    dd <- d %>% filter(!is.na(.data[[v]]), !is.na(.data[[outcome]]))
    y <- dd[[outcome]]; if (scale_y) y <- as.numeric(scale(y))
    dd$yy <- y
    f <- lm(reformulate(c(v, COV_G), "yy"), data = dd)
    s <- coef(summary(f))[v, ]
    tibble(var = v, outcome = outcome, b = s[1], se = s[2], p = s[4], n = nrow(dd),
           eaf = mean(dd[[v]]) / 2)
  }) %>% left_join(snp_tab, by = "var")
}
fz <- snp_assoc("log_hcy") %>% mutate(q = p.adjust(p, "BH"), lo = b - 1.96 * se, hi = b + 1.96 * se,
                                     across(c(b, lo, hi), pct_change))
write.csv(fz, file.path(OUT_TAB, "fig6_variant_hcy_associations.csv"), row.names = FALSE)
PAL_ROLE <- c("Core one-carbon" = "#2B6CB0", "Hcy GWAS (GRS)" = "#4C956C", "Nutrient axis" = "#E09F3E",
              "Negative control" = "#9A9A9A")
fz <- fz %>% mutate(lab = sprintf("%s (%s)", label, allele), lab = fct_reorder(lab, b))
p3 <- ggplot(fz, aes(b, lab, colour = role)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.4, orientation = "y") +
  geom_point(aes(shape = q < 0.05, fill = role), size = 1.2, stroke = 0.35) +
  geom_text(aes(x = Inf, label = if_else(q < 0.05, sprintf("q=%.1e", q), "")), hjust = 1, size = 1.45,
            family = FONT, colour = "grey30") +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_colour_manual(values = PAL_ROLE, name = NULL) + scale_fill_manual(values = PAL_ROLE, guide = "none") +
  scale_x_continuous("% change in Hcy per effect allele (95% CI)", expand = expansion(mult = c(0.04, 0.3))) +
  labs(y = NULL) +
  theme(axis.text.y = element_text(size = BASE_PT - 1.8), legend.position = "bottom",
        legend.key.size = unit(4, "pt"), legend.text = element_text(size = BASE_PT - 1.5)) +
  guides(colour = guide_legend(nrow = 2, order = 1),
         shape = guide_legend(nrow = 2, order = 2, override.aes = list(fill = "grey30", colour = "grey30")))

# ---- P4: variant x biomarker heatmap -----------------------------------------------------------------
BIOM <- c("Hcy" = "log_hcy", "holoTC" = "log_holotc", "Total B12" = "vitamin_b12_pg_ml", "MMA" = "log_mma",
          "RBC folate" = "log_rbcf", "Serum folate" = "folate_serum_ng_ml", "5-MTHF" = "mthf5_plasma_nmol_l",
          "UMFA" = "log_umfa", "PLP" = "log_plp", "Betaine" = "betaine_umol_l", "Choline" = "choline_plasma_umol_l",
          "SAM/SAH" = "sam_sah")
gb <- g %>% mutate(vitamin_b12_pg_ml = log(vitamin_b12_pg_ml))
hm <- map_dfr(unname(BIOM), ~ snp_assoc(.x, d = gb, scale_y = TRUE)) %>%
  mutate(q = p.adjust(p, "BH"), bio = factor(names(BIOM)[match(outcome, BIOM)], names(BIOM)),
         lab = factor(sprintf("%s (%s)", label, allele), levels(fz$lab)))
write.csv(hm, file.path(OUT_TAB, "fig6_variant_biomarker_heatmap.csv"), row.names = FALSE)
p4 <- ggplot(hm, aes(bio, lab, fill = b)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = case_when(q < 0.001 ~ "***", q < 0.01 ~ "**", q < 0.05 ~ "*", TRUE ~ "")),
            size = 1.8, family = FONT, vjust = 0.75) +
  scale_fill_gradient2(low = PAL_DIV["low"], mid = "white", high = PAL_DIV["high"], limits = c(-0.5, 0.5),
                       oob = squish, name = "SD per\nallele") +
  scale_x_discrete(NULL, position = "top") + scale_y_discrete(NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(angle = 50, hjust = 0, vjust = 0, size = BASE_PT - 1.5),
        axis.text.y = element_text(size = BASE_PT - 1.8),
        legend.key.height = unit(12, "pt"), legend.key.width = unit(4, "pt"))

# ---- P5: MTHFR 677 x RBC folate ------------------------------------------------------------------------
g5 <- g %>% filter(!is.na(mthfr677))
f5 <- lm(log_hcy ~ mthfr677 * log_rbcf + log_egfr + log_holotc + age_c10 + male + community, data = g5)
p_i5 <- anova(update(f5, . ~ . - mthfr677:log_rbcf), f5)$`Pr(>F)`[2]
nd5 <- expand_grid(mthfr677 = factor(c("CC", "CT", "TT"), c("CC", "CT", "TT")),
                   log_rbcf = seq(quantile(g5$log_rbcf, 0.03), quantile(g5$log_rbcf, 0.97), length.out = 60)) %>%
  mutate(log_egfr = median(g5$log_egfr), log_holotc = median(g5$log_holotc), age_c10 = median(g5$age_c10),
         male = 0.5, community = 0)
pr5 <- predict(f5, nd5, se.fit = TRUE)
nd5 <- nd5 %>% mutate(fit = exp(pr5$fit), lo = exp(pr5$fit - 1.96 * pr5$se.fit), hi = exp(pr5$fit + 1.96 * pr5$se.fit))
PAL_MTHFR <- c(CC = "#9DB9D8", CT = "#4F7CAC", TT = "#9E2A2B")
p5 <- ggplot(nd5, aes(exp(log_rbcf), fit, colour = mthfr677, fill = mthfr677)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = PAL_MTHFR, name = "MTHFR 677") + scale_fill_manual(values = PAL_MTHFR, name = "MTHFR 677") +
  scale_x_log10("RBC folate (ng/mL)", breaks = c(300, 400, 600, 800, 1000)) +
  scale_y_continuous("Adjusted Hcy (µmol/L)") +
  labs(subtitle = sprintf("Genotype × folate %s", fmt_p(p_i5))) +
  theme(legend.position = "inside", legend.position.inside = c(1, 1), legend.justification = c(1, 1),
        legend.key.size = unit(5, "pt"))

# ---- P6: MTHFR 677 x riboflavin (EGRAC) --------------------------------------------------------------------
g6 <- g %>% filter(!is.na(mthfr677), !is.na(egrac)) %>%
  mutate(b2 = factor(if_else(egrac >= 1.3, "EGRAC ≥ 1.30\n(low B2)", "EGRAC < 1.30\n(adequate B2)"),
                     c("EGRAC < 1.30\n(adequate B2)", "EGRAC ≥ 1.30\n(low B2)")))
f6 <- lm(log_hcy ~ mthfr677 * egrac + log_egfr + log_holotc + log_rbcf + age_c10 + male + community, data = g6)
p_i6 <- anova(update(f6, . ~ . - mthfr677:egrac), f6)$`Pr(>F)`[2]
p6 <- ggplot(g6, aes(mthfr677, hcy, colour = mthfr677)) +
  geom_hline(yintercept = 15, linetype = "22", linewidth = 0.3, colour = COL_REF) +
  geom_point(position = position_jitter(width = 0.18, height = 0, seed = 2), size = 0.35, alpha = 0.35, stroke = 0) +
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange", size = 0.25, linewidth = 0.5, colour = "black",
               fun.args = list(B = 500)) +
  facet_wrap(~b2, nrow = 1) +
  scale_colour_manual(values = PAL_MTHFR, guide = "none") +
  scale_y_log10(lab_umol, breaks = c(5, 10, 15, 20, 30)) +
  scale_x_discrete("MTHFR 677 genotype") +
  labs(subtitle = sprintf("Genotype × EGRAC %s", fmt_p(p_i6)))

# ---- P7: GRS tertile x low holoTC ---------------------------------------------------------------------------
g7 <- g %>% mutate(joint = interaction(grs_tert, if_else(holotc < 35, "holoTC < 35", "holoTC ≥ 35"), sep = " | "))
j7 <- g7 %>% group_by(grs_tert, b12 = if_else(holotc < 35, "holoTC < 35", "holoTC ≥ 35")) %>%
  summarise(n = n(), p = mean(hhcy), lo = binom.test(sum(hhcy), n())$conf.int[1], hi = binom.test(sum(hhcy), n())$conf.int[2])
f7 <- glm(hhcy ~ as.numeric(grs_tert) * holotc_low + age_c10 + male + community + log_egfr, data = g7, family = binomial)
p_i7 <- coef(summary(f7))["as.numeric(grs_tert):holotc_low", 4]
p7 <- ggplot(j7, aes(grs_tert, 100 * p, fill = b12)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), position = position_dodge(0.8), width = 0.2, linewidth = 0.3) +
  geom_text(aes(y = 100 * hi, label = paste0("n=", n)), position = position_dodge(0.8), vjust = -0.5,
            size = 1.5, family = FONT, colour = "grey30") +
  scale_fill_manual(values = c("holoTC < 35" = "#9E2A2B", "holoTC ≥ 35" = "#D9C3B0"), name = NULL) +
  scale_y_continuous("Hcy > 15 µmol/L (%, 95% CI)", expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete("Hcy genetic risk score tertile") +
  labs(subtitle = sprintf("GRS × low holoTC (logistic) %s", fmt_p(p_i7))) +
  theme(legend.position = "inside", legend.position.inside = c(0, 1), legend.justification = c(0, 1),
        legend.key.size = unit(5, "pt"))

# ---- Assemble -------------------------------------------------------------------------------------------------
left  <- wrap_elements(full = p3)
right <- wrap_elements(full = p4)
r1 <- row_of(p1, p2, widths = c(1.15, 1))
r2 <- (left | right) + plot_layout(widths = c(1, 1.05))
r3 <- row_of(p5, p6, p7, widths = c(0.9, 1, 1))
fig6 <- wrap_elements(full = r1) / wrap_elements(full = r2) / wrap_elements(full = r3) +
  plot_layout(heights = c(0.8, 1.35, 0.75))
save_fig(fig6, "Figure6_genetics", width = 180, height = 225)
save_panels(list(variance_partition = p1, per_sd = p2, forest = p3, heatmap = p4, mthfr_folate = p5,
                 mthfr_riboflavin = p6, grs_b12 = p7), "Figure6")
