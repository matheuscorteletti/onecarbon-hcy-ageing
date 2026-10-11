# =============================================================================
# Figure 3 — From diet to biochemical status under mandatory fortification (Q2)
#   P1  Energy-adjusted intake relative to the EAR/AI, by group (raincloud)
#   P2  Mirror plot: inadequate intake (% < EAR) vs biochemical insufficiency, by group
#   P3  Sources of folate: natural food vs fortified flour vs supplements (+ % > UL)
#   P4  Total folic acid intake -> plasma UMFA, by DHFR 19-bp deletion genotype
#   P5–P8  Intake -> biomarker: DFE -> RBC folate; B12 -> holoTC; B6 -> PLP; B2 -> EGRAC
#   P9  Diet–biomarker concordance (Spearman rho, bootstrap 95% CI) by group
#   P10 RBC folate x holoTC quadrants, coloured by homocysteine
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages(library(ggdist))

# ---- Reference values (IOM/NASEM, adults >=51 y) -------------------------------
dat <- dat %>% mutate(
  ear_folate = 320, ear_b12 = 2.0,
  ear_b6 = if_else(sex == "Female", 1.3, 1.4), ear_b2 = if_else(sex == "Female", 0.9, 1.1),
  ai_choline = if_else(sex == "Female", 425, 550),
  r_folate = total_folate_dfe_ug_eadj / ear_folate, r_b12 = vitamin_b12_ug_eadj / ear_b12,
  r_b6 = vitamin_b6_mg_eadj / ear_b6, r_b2 = riboflavin_mg_eadj / ear_b2, r_chol = choline_mg_eadj / ai_choline)
NUTS <- c("Folate (DFE)" = "r_folate", "Vitamin B12" = "r_b12", "Vitamin B6" = "r_b6",
          "Riboflavin" = "r_b2", "Choline (AI)" = "r_chol")

# ---- P1: intake relative to EAR ---------------------------------------------------
p1df <- dat %>% select(grp4, all_of(unname(NUTS))) %>% pivot_longer(-grp4) %>%
  mutate(nut = factor(names(NUTS)[match(name, NUTS)], names(NUTS)))
p1 <- ggplot(p1df, aes(value, fct_rev(grp4), fill = grp4, colour = grp4)) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.3, colour = COL_REF) +
  stat_slab(height = 0.85, justification = -0.1, alpha = 0.75, slab_linewidth = 0, normalize = "groups") +
  stat_pointinterval(.width = c(0.5, 0.95), point_size = 0.8, interval_size_range = c(0.3, 0.9),
                     position = position_nudge(y = -0.08)) +
  facet_wrap(~nut, nrow = 1) +
  scale_x_log10("Energy-adjusted intake / EAR (log scale)", breaks = c(0.25, 0.5, 1, 2, 4),
                labels = c("0.25", "0.5", "1", "2", "4")) +
  scale_y_discrete(NULL, labels = lab_grp4) +
  scale_fill_manual(values = PAL_GRP4, guide = "none") +
  scale_colour_manual(values = colorspace::darken(PAL_GRP4, 0.3), guide = "none") +
  coord_cartesian(xlim = c(0.2, 6))

# ---- P2: mirror — intake inadequacy vs biochemical insufficiency --------------------
mir <- dat %>% group_by(grp4) %>% summarise(
  `Folate_intake` = mean(r_folate < 1), `Folate_status` = mean(rbc_fol < 140 | folate_serum_ng_ml < 4),
  `B12_intake` = mean(r_b12 < 1),       `B12_status` = mean(holotc < 35),
  `B6_intake` = mean(r_b6 < 1),         `B6_status` = mean(plp_nmol_l < 30),
  `Riboflavin_intake` = mean(r_b2 < 1), `Riboflavin_status` = mean(egrac >= 1.3, na.rm = TRUE)) %>%
  pivot_longer(-grp4, names_to = c("nut", "type"), names_sep = "_") %>%
  mutate(nut = factor(nut, c("Folate", "B12", "B6", "Riboflavin"),
                      c("Folate", "Vitamin B12", "Vitamin B6", "Riboflavin")),
         x = if_else(type == "intake", -value, value))
p2 <- ggplot(mir, aes(x, fct_rev(grp4), fill = grp4, alpha = type)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_text(aes(label = round(100 * value), hjust = if_else(type == "intake", 1.15, -0.15)),
            size = 1.65, family = FONT, alpha = 1) +
  facet_wrap(~nut, nrow = 1) +
  scale_alpha_manual(values = c(intake = 0.45, status = 1), guide = "none") +
  scale_fill_manual(values = PAL_GRP4, guide = "none") +
  scale_x_continuous("Participants (%)", labels = function(x) paste0(abs(100 * x)), limits = c(-0.75, 0.75),
                     breaks = c(-0.5, -0.25, 0, 0.25, 0.5)) +
  scale_y_discrete(NULL, labels = lab_grp4) +
  labs(subtitle = "\u2190 Intake < EAR (%)     Biochemical insufficiency (%) \u2192") +
  coord_cartesian(clip = "off") +
  theme(plot.subtitle = element_text(hjust = 0.5))

# ---- P3: folate sources -----------------------------------------------------------
src <- dat %>% group_by(grp4) %>% summarise(
  `Natural food folate` = mean(food_folate_ug),
  `Fortified flour (folic acid × 1.7)` = mean(1.7 * folic_acid_fortified_ug),
  `Supplements (folic acid × 1.7)` = mean(1.7 * folic_acid_supplement_ug),
  above_ul = mean(folic_acid_above_ul == 1)) %>%
  pivot_longer(2:4, names_to = "src") %>%
  mutate(src = factor(src, c("Supplements (folic acid × 1.7)", "Fortified flour (folic acid × 1.7)", "Natural food folate"))) %>%
  group_by(grp4) %>% mutate(share = value / sum(value), tot = sum(value)) %>% ungroup()
p3 <- ggplot(src, aes(grp4, value, fill = src)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = if_else(share > 0.08, paste0(round(100 * share), "%"), "")),
            position = position_stack(vjust = 0.5), size = 1.65, family = FONT, colour = "white") +
  geom_text(data = distinct(src, grp4, tot, above_ul), aes(grp4, tot, label = sprintf(">UL %.0f%%", 100 * above_ul)),
            inherit.aes = FALSE, vjust = -0.5, size = 1.65, family = FONT) +
  geom_hline(yintercept = 320, linetype = "22", linewidth = 0.3, colour = COL_REF) +
  annotate("text", x = 4.45, y = 320, label = "EAR", vjust = -0.4, hjust = 1, size = 1.65, family = FONT, colour = COL_REF) +
  scale_fill_manual(values = c("#C9A227", "#7A9E7E", "#3D5A40"), name = NULL) +
  scale_y_continuous("Mean folate intake (µg DFE/day)", expand = expansion(mult = c(0, 0.12))) +
  scale_x_discrete(NULL, labels = lab_grp4) + rot_x +
  theme(legend.position = "top", legend.direction = "vertical", legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = BASE_PT - 1.5))

# ---- P4: folic acid -> UMFA by DHFR ---------------------------------------------------
p4df <- dat %>% filter(!is.na(dhfr))
fit4 <- lm(log_umfa ~ log1p(folic_acid_total_ug) * as.numeric(dhfr) + sex + age_c10, data = p4df)
p_int4 <- coef(summary(fit4))["log1p(folic_acid_total_ug):as.numeric(dhfr)", 4]
PAL_DHFR <- c("ins/ins" = "#9DB9D8", "ins/del" = "#4F7CAC", "del/del" = "#1B2F4B")
p4 <- ggplot(p4df, aes(folic_acid_total_ug + 1, umfa_nmol_l, colour = dhfr)) +
  geom_hline(yintercept = 0.2, linetype = "22", linewidth = 0.3, colour = COL_REF) +
  geom_point(size = 0.45, alpha = 0.45, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.55, alpha = 0.15, aes(fill = dhfr)) +
  scale_x_log10("Total folic acid intake (µg/day + 1)", breaks = c(10, 100, 1000, 5000),
                labels = label_number(big.mark = ",")) +
  scale_y_log10("Plasma UMFA (nmol/L)") +
  scale_colour_manual(values = PAL_DHFR, name = "DHFR") + scale_fill_manual(values = PAL_DHFR, name = "DHFR") +
  annotate("text", x = 10, y = Inf, hjust = 0, vjust = 1.2, size = 1.75, family = FONT,
           label = sprintf("Intake × genotype %s", fmt_p(p_int4))) +
  theme(legend.position = "inside", legend.position.inside = c(1, 0), legend.justification = c(1, 0))

# ---- P5–P8: intake -> biomarker ----------------------------------------------------------
pairs <- tribble(
  ~x,                         ~y,                 ~xl,                                   ~yl,                         ~cut,  ~rev,
  "total_folate_dfe_ug_eadj", "rbc_fol",          "Total folate (µg DFE/day, e-adj.)", "RBC folate (ng/mL)",       NA,   FALSE,
  "vitamin_b12_ug_eadj",      "holotc",           "Vitamin B12 (µg/day, e-adj.)",      "holoTC (pmol/L)",          35,    FALSE,
  "vitamin_b6_mg_eadj",       "plp_nmol_l",       "Vitamin B6 (mg/day, e-adj.)",            "PLP (nmol/L)",             30,    FALSE,
  "riboflavin_mg_eadj",       "egrac",            "Riboflavin (mg/day, e-adj.)",            "EGRAC (higher = worse)",   1.3,   TRUE)
scatter_pair <- function(i) {
  pr <- pairs[i, ]
  df <- dat %>% filter(!is.na(.data[[pr$y]]))
  rho_g <- df %>% group_by(grp4) %>%
    summarise(r = cor(.data[[pr$x]], .data[[pr$y]], method = "spearman", use = "complete.obs"))
  lab_g <- paste0(sprintf("<span style='color:%s'>●</span> %.2f", PAL_GRP4[as.character(rho_g$grp4)], rho_g$r), collapse = "  ")
  ggplot(df, aes(.data[[pr$x]], .data[[pr$y]])) +
    { if (!is.na(pr$cut)) geom_hline(yintercept = pr$cut, linetype = "22", linewidth = 0.3, colour = COL_REF) } +
    geom_point(aes(colour = grp4), size = 0.45, alpha = 0.5, stroke = 0) +
    geom_smooth(method = "loess", formula = y ~ x, span = 0.9, colour = "black", fill = "grey70",
                linewidth = 0.5, alpha = 0.35) +
    scale_x_log10(pr$xl) +
    { if (pr$rev) scale_y_continuous(pr$yl) else scale_y_log10(pr$yl) } +
    scale_colour_manual(values = PAL_GRP4, guide = "none") +
    labs(subtitle = paste0(spearman_lab(df[[pr$x]], df[[pr$y]]), "<br>", lab_g)) +
    theme(plot.subtitle = ggtext::element_markdown(size = BASE_PT - 1.5, lineheight = 1.15))
}
p5to8 <- map(1:4, scatter_pair)

# ---- P9: concordance by group (bootstrap CI) ----------------------------------------------
boot_rho <- function(x, y, B = N_BOOT) {
  ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]
  r <- cor(x, y, method = "spearman")
  bs <- replicate(B, { i <- sample.int(length(x), replace = TRUE); cor(x[i], y[i], method = "spearman") })
  c(r = r, lo = unname(quantile(bs, 0.025)), hi = unname(quantile(bs, 0.975)))
}
conc <- pmap_dfr(pairs, function(x, y, xl, yl, cut, rev) {
  dat %>% group_by(grp4) %>% group_modify(~ as_tibble_row(boot_rho(.x[[x]], .x[[y]]) * if (rev) -1 else 1)) %>%
    mutate(pair = paste0(str_extract(xl, "^[^(]+") %>% str_trim(), " → ", str_extract(yl, "^[^(]+") %>% str_trim()))
}) %>% mutate(lo2 = pmin(lo, hi), hi2 = pmax(lo, hi), pair = factor(pair, unique(pair)))
p9 <- ggplot(conc, aes(r, fct_rev(grp4), colour = grp4)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = COL_REF) +
  geom_pointrange(aes(xmin = lo2, xmax = hi2), size = 0.2, linewidth = 0.5) +
  scale_colour_manual(values = PAL_GRP4, name = NULL, labels = lab_grp4) +
  labs(x = "Diet–biomarker ρ (95% CI)\n(EGRAC sign inverted)", y = NULL) +
  facet_wrap(~pair, nrow = 1, scales = "free_y") +
  theme(legend.position = "right", legend.key.size = unit(5, "pt"), legend.text = element_text(size = BASE_PT - 1.5),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

# ---- P10: quadrants RBC folate x holoTC ---------------------------------------------------
q <- dat %>% filter(!is.na(hcy)) %>%
  mutate(q = paste(if_else(rbc_fol >= 400, "High folate", "Lower folate"), if_else(holotc >= 35, "adequate B12", "low B12"), sep = " / "))
xr <- range(q$rbc_fol)
q_lab <- q %>% group_by(q) %>% summarise(n = n(), med = median(hcy), phh = mean(hhcy)) %>%
  mutate(x = if_else(str_detect(q, "^High"), xr[2] * 0.98, xr[1] * 1.02), hj = if_else(str_detect(q, "^High"), 1.02, -0.02), y = if_else(str_detect(q, "adequate"), 200, 9))
p10 <- ggplot(q %>% arrange(hcy), aes(rbc_fol, holotc, colour = hcy)) +
  geom_vline(xintercept = 400, linetype = "22", linewidth = 0.3, colour = COL_REF) +
  geom_hline(yintercept = 35, linetype = "22", linewidth = 0.3, colour = COL_REF) +
  geom_point(size = 0.6, alpha = 0.85, stroke = 0) +
  geom_label(data = q_lab, aes(x, y, hjust = hj, label = sprintf("n = %d\nHcy %.1f\n>15: %.0f%%", n, med, 100 * phh)),
             inherit.aes = FALSE, size = 1.6, family = FONT, lineheight = 0.85, linewidth = 0,
             fill = alpha("white", 0.85), label.padding = unit(1.5, "pt")) +
  scale_colour_gradientn(colours = PAL_SEQ, trans = "log10", limits = c(6, 35), oob = squish,
                         name = "Hcy\n(µmol/L)", breaks = c(8, 15, 30)) +
  scale_x_log10("RBC folate (ng/mL)") + scale_y_log10("holoTC (pmol/L)") +
  theme(legend.key.height = unit(12, "pt"), legend.key.width = unit(4, "pt"))

# ---- Assemble ---------------------------------------------------------------------------
row3 <- (p3 | wrap_elements(full = p4) | wrap_elements(full = p10)) + plot_layout(widths = c(0.8, 1, 1.05))
row4 <- wrap_plots(p5to8, nrow = 1)
fig3 <- wrap_elements(full = p1) / wrap_elements(full = p2) / wrap_elements(full = row3) /
  wrap_elements(full = row4) / wrap_elements(full = p9) +
  plot_layout(heights = c(0.62, 0.5, 0.95, 0.72, 0.5))
save_fig(fig3, "Figure3_diet_to_status", width = 180, height = 235)
save_panels(c(list(intake_ear = p1, mirror = p2, folate_sources = p3, umfa_dhfr = p4),
              setNames(p5to8, c("dfe_rbc", "b12_holotc", "b6_plp", "b2_egrac")),
              list(concordance = p9, quadrants = p10)), "Figure3")
