# =============================================================================
# Supplementary figures and tables
#   Fig S1  H1b: is there a folate "ceiling"? Piecewise slopes of Hcy on RBC folate below/above cut-points
#   Fig S2  H1c: folate x B12 interaction for Hcy and MMA (expected null)
#   Fig S3  All gene x nutrient interactions (24 variants x 6 nutrients; -log10 P, FDR)
#   Fig S4  Genotype QC: call rate, HWE and allele frequencies
#   Fig S5  Exploratory mortality: Kaplan–Meier by Hcy > 15 and Cox HR per 5 umol/L by stratum
#   Table S1 Sensitivity analyses of the fully adjusted model (M5)
#   Table S2 Gene x nutrient interactions (all)
#   Table S3 Tobit vs OLS for serum folate (right-censored at 24 ng/mL)
# =============================================================================
source("R/00_setup.R"); source("R/01_models.R")
suppressPackageStartupMessages({ library(survival); library(flextable); library(officer); library(AER) })
select <- dplyr::select

# ---- Fig S1: folate ceiling ---------------------------------------------------------------
cuts <- c(400, 500, 600, 700, 800)
adj_m5 <- setdiff(unlist(BLOCKS), "log_rbcf")
s1 <- map_dfr(cuts, function(k) {
  d <- mdat %>% mutate(lo = pmin(log_rbcf, log(k)), hi = pmax(log_rbcf - log(k), 0))
  f <- lm(reformulate(c("lo", "hi", adj_m5), "log_hcy"), data = d)
  ci <- confint(f)
  lh <- car::linearHypothesis(f, "lo = hi")
  tibble(cut = k, seg = c("Below cut-point", "Above cut-point"), b = coef(f)[c("lo", "hi")],
         lo = ci[c("lo", "hi"), 1], hi = ci[c("lo", "hi"), 2], p_diff = lh$`Pr(>F)`[2])
}) %>% mutate(across(c(b, lo, hi), ~ pct_change(.x * log(1.1))))   # % change in Hcy per 10% higher RBC folate
pS1 <- ggplot(s1, aes(factor(cut), b, colour = seg)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.5), size = 0.25, linewidth = 0.5) +
  geom_text(data = distinct(s1, cut, p_diff), aes(factor(cut), max(s1$hi) + 0.25, label = paste("P diff", fmt_p_short(p_diff))),
            inherit.aes = FALSE, size = 1.8, family = FONT) +
  scale_colour_manual(values = c("Below cut-point" = "#7FB7A4", "Above cut-point" = "#2F6B55"), name = NULL) +
  labs(x = "RBC folate cut-point (ng/mL)", y = "% change in Hcy per 10% higher RBC folate (95% CI)") +
  theme(legend.position = "top")
save_fig(pS1, "FigS1_folate_ceiling", width = 120, height = 80, dir = OUT_SUPP)

# ---- Fig S2: folate x B12 -------------------------------------------------------------------
f_h <- lm(log_hcy ~ log_rbcf * log_holotc + stratum + log_egfr + male + glim_any, data = mdat)
f_m <- lm(log_mma ~ log_rbcf * log_holotc + stratum + log_egfr + male + glim_any, data = mdat)
p_h <- coef(summary(f_h))["log_rbcf:log_holotc", 4]; p_m <- coef(summary(f_m))["log_rbcf:log_holotc", 4]
qf <- quantile(mdat$log_rbcf, c(0.1, 0.5, 0.9))
nd <- expand_grid(log_holotc = seq(log(15), log(150), length.out = 60), log_rbcf = qf) %>%
  mutate(stratum = factor(STRATA[1], STRATA), log_egfr = median(mdat$log_egfr), male = 0.5, glim_any = 0,
         fol = factor(round(exp(log_rbcf)), labels = paste0(c("P10", "P50", "P90"), ": ", round(exp(qf)), " ng/mL")))
s2 <- bind_rows(nd %>% mutate(out = "Homocysteine (µmol/L)", fit = exp(predict(f_h, nd))),
                nd %>% mutate(out = "MMA (µmol/L)", fit = exp(predict(f_m, nd))))
pS2 <- ggplot(s2, aes(exp(log_holotc), fit, colour = fol)) +
  geom_vline(xintercept = 35, linetype = "22", linewidth = 0.3, colour = COL_REF) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~out, scales = "free_y",
             labeller = labeller(out = function(x) paste0(x, "\nfolate × holoTC ", c(fmt_p(p_h), fmt_p(p_m))[match(x, unique(s2$out))]))) +
  scale_colour_manual(values = c("#C9E2D6", "#7FB7A4", "#2F6B55"), name = "RBC folate") +
  scale_x_log10("holoTC (pmol/L)", breaks = c(15, 35, 70, 150)) + labs(y = "Adjusted value") +
  theme(legend.position = "top")
save_fig(pS2, "FigS2_folate_x_B12", width = 140, height = 80, dir = OUT_SUPP)

# ---- Fig S3 / Table S2: all gene x nutrient interactions ----------------------------------------
g <- gxn %>% mutate(male = as.integer(sex == "Male"), community = as.integer(!hospital)) %>% filter(!is.na(log_hcy))
NUTR <- c("RBC folate" = "log_rbcf", "holoTC" = "log_holotc", "PLP" = "log_plp", "EGRAC" = "egrac",
          "Folate intake (DFE)" = "total_folate_dfe_ug_eadj", "Choline intake" = "choline_mg_eadj")
gxn_int <- expand_grid(v = snp_tab$var, n = unname(NUTR)) %>% pmap_dfr(function(v, n) {
  d <- g %>% filter(!is.na(.data[[v]]), !is.na(.data[[n]])) %>% mutate(nz = as.numeric(scale(log(.data[[n]] + 1e-6 * (.data[[n]] <= 0)))))
  if (n %in% c("log_rbcf", "log_holotc", "log_plp", "egrac")) d$nz <- as.numeric(scale(d[[n]]))
  f <- lm(as.formula(paste("log_hcy ~", v, "* nz + age_c10 + male + community + log_egfr + ancestry_pc1 + ancestry_pc2")), data = d)
  s <- coef(summary(f))[paste0(v, ":nz"), ]
  tibble(var = v, nutrient = names(NUTR)[match(n, NUTR)], b_int = s[1], se = s[2], p = s[4], n = nrow(d))
}) %>% mutate(q = p.adjust(p, "BH")) %>% left_join(snp_tab, by = "var")
pS3 <- ggplot(gxn_int, aes(nutrient, fct_rev(factor(label, snp_tab$label)), fill = sign(b_int) * -log10(p))) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = case_when(q < 0.05 ~ "★", p < 0.05 ~ "•", TRUE ~ "")), size = 2, family = FONT) +
  scale_fill_gradient2(low = PAL_DIV["low"], mid = "white", high = PAL_DIV["high"], limits = c(-4, 4), oob = squish,
                       name = "signed\n−log10 P") +
  scale_x_discrete(NULL, position = "top") + scale_y_discrete(NULL) +
  labs(caption = "★ FDR < 0.05; • nominal P < 0.05.\nInteraction per SD of the nutrient marker; adjusted for age, sex, setting, eGFR and ancestry.") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(angle = 35, hjust = 0, vjust = 0), plot.caption = element_text(size = BASE_PT - 1.5, hjust = 0))
save_fig(pS3, "FigS3_gene_x_nutrient", width = 120, height = 140, dir = OUT_SUPP)

# ---- Fig S4: genotype QC ----------------------------------------------------------------------------
qc <- geno_qc %>% janitor::clean_names() %>% mutate(label = paste(gene, rs_id)) %>%
  semi_join(snp_tab, by = "label") %>%
  mutate(label = factor(label, rev(snp_tab$label)))
pS4a <- ggplot(qc, aes(-log10(hwe_p_ctrl), label)) +
  geom_vline(xintercept = -log10(0.05 / nrow(qc)), linetype = "22", linewidth = 0.3, colour = "#9E2A2B") +
  geom_col(width = 0.6, fill = "#8FA9C4") + labs(x = "HWE −log10 P (controls)", y = NULL)
pS4b <- ggplot(qc, aes(100 * call_rate, label)) + geom_point(size = 1.2, colour = "#2B6CB0") +
  scale_x_continuous("Call rate (%)", limits = c(99, 100)) + labs(y = NULL) + theme(axis.text.y = element_blank())
pS4c <- ggplot(qc, aes(eaf_controls, eaf_longevos)) + geom_abline(linetype = "22", linewidth = 0.3) +
  geom_point(size = 1.2, colour = "#9E2A2B") +
  ggrepel::geom_text_repel(aes(label = gene), size = 1.6, family = FONT, segment.size = 0.15, max.overlaps = 30) +
  labs(x = "Effect-allele frequency, controls", y = "Effect-allele frequency, ≥85 y")
save_fig(pS4a | pS4b | pS4c, "FigS4_genotype_QC", width = 180, height = 90, dir = OUT_SUPP)

# ---- Fig S5: exploratory mortality --------------------------------------------------------------------
sv <- hosp %>% filter(vital_status != "Lost to follow-up" | TRUE) %>%
  mutate(time = followup_months, ev = death_event, hgrp = if_else(hhcy == 1, "Hcy > 15", "Hcy ≤ 15"),
         hcy5 = hcy / 5, male = as.integer(sex == "Male"))
km <- survfit(Surv(time, ev) ~ stratum + hgrp, data = sv)
kmd <- broom::tidy(km) %>% separate(strata, c("stratum", "hgrp"), sep = ", ") %>%
  mutate(stratum = str_remove(stratum, "stratum="), hgrp = str_remove(hgrp, "hgrp="),
         stratum = factor(stratum, STRATA))
pS5a <- ggplot(kmd, aes(time, estimate, colour = stratum, linetype = hgrp)) +
  geom_step(linewidth = 0.5) +
  scale_colour_manual(values = PAL_STR, name = NULL) + scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  scale_y_continuous("Survival probability", limits = c(0, 1)) + labs(x = "Months since admission") +
  theme(legend.position = "right")
cox <- map_dfr(STRATA, function(s) {
  f <- coxph(Surv(time, ev) ~ hcy5 + age_entry_exact + male + log_egfr + glim_any + charlson_age_adj,
             data = sv %>% filter(stratum == s))
  ci <- exp(confint(f))["hcy5", ]
  tibble(stratum = s, hr = exp(coef(f)["hcy5"]), lo = ci[1], hi = ci[2], p = coef(summary(f))["hcy5", 5])
}) %>% mutate(stratum = factor(stratum, STRATA))
pS5b <- ggplot(cox, aes(hr, fct_rev(stratum), colour = stratum)) +
  geom_vline(xintercept = 1, linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.3, linewidth = 0.5) +
  geom_text(aes(x = hi, label = fmt_p(p)), hjust = -0.15, size = 1.8, family = FONT, colour = "black") +
  scale_colour_manual(values = PAL_STR, guide = "none") + scale_x_log10(expand = expansion(mult = c(0.05, 0.4))) +
  labs(x = "HR per 5 µmol/L Hcy (adj. age, sex, eGFR, GLIM, Charlson)", y = NULL)
save_fig(pS5a | pS5b, "FigS5_mortality_exploratory", width = 180, height = 75, dir = OUT_SUPP)

# ---- Table S1: sensitivity analyses of M5 ---------------------------------------------------------------
sets <- list("Main (complete cases)" = mdat,
             "Excl. AKI at admission" = mdat %>% filter(acute_kidney_injury == 0),
             "Excl. B12 supplement / antifolate drugs" = mdat %>% filter(sens_exclude_b12_antifolate == 0),
             "Excl. B12 supplement users" = mdat %>% filter(supp_vitamin_b12 == 0),
             "Excl. levodopa users" = mdat %>% filter(med_levodopa == 0),
             "Reliable FFQ only" = mdat %>% filter(ffq_reliable == 1),
             "eGFR by creatinine (CKD-EPI 2021)" = mdat %>% mutate(log_egfr = log(egfr_ckdepi2021)),
             "Original design (without 76–84 y)" = mdat %>% filter(stratum != STRATA[2]))
keyv <- c(s85 = "≥85 y", log_egfr = "eGFR-cys (per SD)", log_holotc = "holoTC (per SD)",
          log_rbcf = "RBC folate (per SD)", glim_any = "GLIM (any)", phase_angle_deg = "Phase angle (per SD)")
TS1 <- imap_dfr(sets, function(d, nm) {
  cv <- c("log_egfr", "log_holotc", "log_rbcf", "phase_angle_deg")
  d <- d %>% mutate(across(all_of(cv), ~ as.numeric(scale(.x))))
  f <- FULL_FORMULA; if (nm == "Original design (without 76–84 y)") f <- update(f, . ~ . - s76)
  m <- lm(f, data = d); ci <- confint(m)
  tibble(Analysis = nm, n = nobs(m),
         !!!setNames(map(names(keyv), ~ sprintf("%.1f (%.1f, %.1f)", pct_change(coef(m)[.x]), pct_change(ci[.x, 1]), pct_change(ci[.x, 2]))),
                     unname(keyv)))
})
# Multiple imputation sensitivity (hospital set, chained equations)
suppressPackageStartupMessages(library(mice))
imp_vars <- c("log_hcy", "s76", "s85", "male", "ancestry_pc1", "ancestry_pc2", "time_to_centrifugation_min", "fasting_hours",
              "log_egfr", "log_holotc", "log_rbcf", "log_plp", "egrac", "glim_any", "phase_angle_deg", "log_crp",
              "diabetes_t2", "hcy_genetic_risk_score")
hi <- hosp %>% mutate(s76 = as.integer(stratum == STRATA[2]), s85 = as.integer(stratum == STRATA[3]),
                      male = as.integer(sex == "Male")) %>% select(all_of(imp_vars)) %>%
  mutate(across(c(log_egfr, log_holotc, log_rbcf, phase_angle_deg), ~ (.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE)))
imp <- mice(hi, m = 20, printFlag = FALSE, seed = 1)
pool_s <- summary(pool(as.mira(lapply(seq_len(imp$m), function(i) lm(FULL_FORMULA, data = complete(imp, i))))), conf.int = TRUE) %>% as_tibble() %>% filter(term %in% names(keyv))
TS1 <- bind_rows(TS1, tibble(Analysis = "Multiple imputation (m = 20)", n = nrow(hi),
                             !!!setNames(map(names(keyv), function(v) { r <- pool_s %>% filter(term == v)
                               sprintf("%.1f (%.1f, %.1f)", pct_change(r$estimate), pct_change(r$`2.5 %`), pct_change(r$`97.5 %`)) }),
                               unname(keyv))))

# ---- Table S3: Tobit vs OLS for serum folate --------------------------------------------------------------
sf <- dat %>% filter(hospital) %>% mutate(male = as.integer(sex == "Male"), lfa = log1p(folic_acid_total_ug))
ols <- lm(folate_serum_ng_ml ~ stratum + lfa + male + log_egfr, data = sf)
tob <- tobit(folate_serum_ng_ml ~ stratum + lfa + male + log_egfr, right = 24, left = -Inf, data = sf)
tidy_m <- function(m, lab) { s <- coef(summary(m)); ci <- confint(m)
  rn <- intersect(rownames(s), rownames(ci)); rn <- rn[rn != "(Intercept)"]
  tibble(term = rn, !!lab := sprintf("%.2f (%.2f, %.2f)", s[rn, 1], ci[rn, 1], ci[rn, 2])) }
TS3 <- tidy_m(ols, "OLS") %>% left_join(tidy_m(tob, "Tobit (right-censored at 24 ng/mL)"), by = "term") %>%
  mutate(term = dplyr::recode(term, "stratum76–84 y" = "76–84 y (vs 56–75 y)", "stratum≥85 y" = "≥85 y (vs 56–75 y)",
                       lfa = "Total folic acid intake (log1p µg/day)", male = "Male", log_egfr = "eGFR-cys (log)"))

# ---- Write supplementary tables ------------------------------------------------------------------------------
ftx <- function(tb, cap) flextable(tb) %>% padding(padding.top = 1, padding.bottom = 1, part = "all") %>%
  font(fontname = "Arial", part = "all") %>% fontsize(size = 8, part = "all") %>% bold(part = "header") %>%
  border_remove() %>% hline_top(border = fp_border(width = 1), part = "all") %>% hline_bottom(border = fp_border(width = 1), part = "all") %>%
  set_caption(caption = cap) %>% autofit()
TS2 <- gxn_int %>% transmute(Variant = label, Role = role, Nutrient = nutrient, n,
                             `Interaction β (SE)` = sprintf("%.3f (%.3f)", b_int, se), P = fmt_p_short(p), `FDR q` = fmt_p_short(q))
TS3qc <- geno_qc %>% janitor::clean_names() %>% mutate(label = paste(gene, rs_id)) %>% semi_join(snp_tab, by = "label") %>%
  transmute(Variant = label, `Effect allele` = effect_allele, `Call rate` = sprintf("%.3f", call_rate),
            `EAF controls` = sprintf("%.3f", eaf_controls), `EAF ≥85 y` = sprintf("%.3f", eaf_longevos),
            `HWE P (controls)` = fmt_p_short(hwe_p_ctrl))
doc <- read_docx() %>%
  body_add_flextable(fit_to_width(ftx(TS1, "Table S1. Sensitivity analyses of the fully adjusted model (M5): % change in homocysteine (95% CI)."), 9.5)) %>%
  body_add_break() %>%
  body_add_flextable(fit_to_width(ftx(TS3, "Table S2. Determinants of serum folate: ordinary least squares vs Tobit model accounting for the assay upper limit (24 ng/mL)."), 9.5)) %>%
  body_add_break() %>%
  body_add_flextable(fit_to_width(ftx(TS3qc, "Table S3. Genotyping quality control of the 24 one-carbon variants."), 9.5)) %>%
  body_add_break() %>%
  body_add_flextable(fit_to_width(ftx(TS2, "Table S4. Gene × nutrient interactions on log-homocysteine (all tests; BH-FDR)."), 9.5)) %>%
  body_end_section_landscape()
print(doc, target = file.path(OUT_TAB, "Supplementary_Tables.docx"))
write.csv(gxn_int, file.path(OUT_TAB, "TableS4_gene_x_nutrient.csv"), row.names = FALSE)
write.csv(TS1, file.path(OUT_TAB, "TableS1_sensitivity.csv"), row.names = FALSE)
message("saved: supplementary tables")
