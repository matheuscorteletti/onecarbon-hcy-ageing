# =============================================================================
# Figure 2 — Nutritional phenotype across age strata (Q1)
#   P1  Prevalence of malnutrition (GLIM), risk (MNA-SF) and sarcopenia by group
#   P2  Alluvial: age stratum -> MNA-SF -> GLIM -> EWGSOP2 (hospital set)
#   P3  Phase angle by group (raincloud)
#   P4  Handgrip by group and sex, EWGSOP2 cut-offs
#   P5  Phase angle x ASMI, coloured by GLIM, 68% ellipses per stratum
#   P6  Protein and energy intake vs ESPEN geriatric targets
#   P7  PCA of the nutritional phenotype: variable loadings (biplot arrows)
#   P8  PCA individuals by group with centroids
#   P9  Gaussian-mixture phenotype clusters: standardized profile heatmap
#   P10 Cluster composition by group
#   P11 Homocysteine by phenotype cluster within hospital strata (bridge to H1e)
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages({ library(ggdist); library(ggalluvial); library(FactoMineR); library(mclust) })
map <- purrr::map; select <- dplyr::select

# ---- P1: prevalence bars --------------------------------------------------------
sev_lab <- c("None", "At risk / moderate / probable", "Confirmed", "Severe / malnourished")
prev <- bind_rows(
  dat %>% transmute(grp4, tool = "GLIM malnutrition", sev = as.integer(glim3) - 1L) %>%
    mutate(sev = if_else(sev == 2L, 3L, sev)),
  dat %>% transmute(grp4, tool = "MNA-SF", sev = as.integer(mna3) - 1L) %>%
    mutate(sev = if_else(sev == 2L, 3L, sev)),
  dat %>% filter(!is.na(sarc4)) %>% transmute(grp4, tool = "Sarcopenia (EWGSOP2)", sev = as.integer(sarc4) - 1L)) %>%
  count(tool, grp4, sev) %>% group_by(tool, grp4) %>% mutate(p = n / sum(n)) %>% ungroup() %>%
  mutate(sev = factor(sev, 0:3, sev_lab),
         tool = factor(tool, c("GLIM malnutrition", "MNA-SF", "Sarcopenia (EWGSOP2)")))
SEV_COL <- setNames(c("#E3ECE8", "#EBC06A", "#DB7B3F", "#9E2A2B"), sev_lab)
p1 <- ggplot(prev, aes(p, fct_rev(grp4), fill = fct_rev(sev))) +
  geom_col(width = 0.74, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = if_else(sev != "None" & p >= 0.06, as.character(round(100 * p)), ""),
                colour = sev == "Severe / malnourished", group = fct_rev(sev)),
            position = position_stack(vjust = 0.5), size = 1.75, family = FONT) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  facet_wrap(~tool, ncol = 1) +
  scale_fill_manual(values = SEV_COL, breaks = sev_lab[-1], name = NULL) +
  scale_x_continuous("Participants (%)", labels = function(x) 100 * x, expand = c(0, 0)) +
  scale_y_discrete(NULL, labels = lab_grp4) +
  theme(legend.position = "bottom", legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = BASE_PT - 1.5)) +
  guides(fill = guide_legend(ncol = 1))

# ---- P2: alluvial -----------------------------------------------------------------
all_df <- dat %>% filter(hospital, !is.na(sarc4)) %>%
  transmute(stratum,
            MNA = fct_recode(mna3, "Risk" = "At risk", "Maln." = "Malnourished"),
            GLIM = fct_recode(glim3, "None " = "No malnutrition", "Mod." = "Moderate", "Sev." = "Severe"),
            Sarcopenia = fct_recode(sarc4, "None  " = "No sarcopenia", "Prob." = "Probable",
                                    "Conf." = "Confirmed", "Sev. " = "Severe")) %>%
  count(stratum, MNA, GLIM, Sarcopenia)
p2 <- ggplot(all_df, aes(axis1 = stratum, axis2 = MNA, axis3 = GLIM, axis4 = Sarcopenia, y = n)) +
  geom_alluvium(aes(fill = stratum), alpha = 0.55, width = 0.3, knot.pos = 0.35, curve_type = "sigmoid") +
  geom_stratum(width = 0.3, fill = "grey97", colour = "grey35", linewidth = 0.25) +
  geom_text(stat = "stratum", aes(label = str_remove(after_stat(stratum), " y$")), size = 1.6, family = FONT, lineheight = 0.8) +
  scale_x_discrete(limits = c("Age", "MNA-SF", "GLIM", "EWGSOP2"), expand = c(0.08, 0.08)) +
  scale_y_continuous("Hospital participants (n)", expand = c(0, 0)) +
  scale_fill_manual(values = PAL_STR, guide = "none") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())

# ---- P3: phase angle raincloud -------------------------------------------------
rain <- function(df, y, ylab, ref = NULL) {
  ggplot(df, aes(grp4, .data[[y]], fill = grp4, colour = grp4)) +
    { if (!is.null(ref)) geom_hline(yintercept = ref, linetype = "22", linewidth = 0.3, colour = COL_REF) } +
    stat_halfeye(adjust = 0.9, width = 0.55, justification = -0.32, .width = 0, point_colour = NA,
                 alpha = 0.75, slab_linewidth = 0, na.rm = TRUE) +
    geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.3, na.rm = TRUE) +
    geom_point(position = position_jitter(width = 0.06, height = 0, seed = 1), size = 0.25,
               alpha = 0.35, stroke = 0, na.rm = TRUE) +
    scale_fill_manual(values = PAL_GRP4, guide = "none") +
    scale_colour_manual(values = colorspace::darken(PAL_GRP4, 0.25), guide = "none") +
    scale_x_discrete(NULL, labels = lab_grp4) + rot_x +
    labs(y = ylab)
}
kw_lab <- function(y) {
  kw <- kruskal.test(dat[[y]] ~ dat$grp4)
  tr <- cor.test(as.numeric(hosp$stratum), hosp[[y]], method = "spearman", exact = FALSE)
  sprintf("Kruskal–Wallis %s\nTrend (hospital) %s", fmt_p(kw$p.value), fmt_p(tr$p.value))
}
p3 <- rain(dat, "phase_angle_deg", "Phase angle (°)") +
  annotate("text", x = 0.55, y = Inf, label = kw_lab("phase_angle_deg"), hjust = 0, vjust = 1.1,
           size = 1.75, family = FONT, lineheight = 0.9)

# ---- P4: handgrip by sex ----------------------------------------------------------
cut_hg <- tibble(sex = factor(c("Female", "Male")), cut = c(16, 27))
p4 <- rain(dat %>% filter(!is.na(handgrip_kg)), "handgrip_kg", "Handgrip strength (kg)") +
  geom_hline(data = cut_hg, aes(yintercept = cut), linetype = "22", linewidth = 0.3, colour = COL_REF) +
  facet_wrap(~sex, nrow = 1) +
  theme(axis.text.x = element_text(size = BASE_PT - 1.5))

# ---- P5: phase angle x ASMI --------------------------------------------------------
p5df <- dat %>% filter(hospital, !is.na(asmi_kg_m2), !is.na(phase_angle_deg)) %>%
  mutate(asmi_z = ave(asmi_kg_m2, sex, FUN = function(v) (v - mean(v)) / sd(v)))
p5 <- ggplot(p5df, aes(asmi_z, phase_angle_deg)) +
  geom_point(aes(colour = glim3), size = 0.6, alpha = 0.75, stroke = 0) +
  { map(STRATA, function(st) stat_ellipse(data = filter(p5df, stratum == st), level = 0.68,
                                          linewidth = 0.55, colour = PAL_STR[[st]])) } +
  geom_label(data = p5df %>% group_by(stratum) %>% summarise(x = mean(asmi_z), y = mean(phase_angle_deg)),
             aes(x, y, label = stratum), colour = PAL_STR, fill = alpha("white", 0.8), linewidth = 0,
             size = 1.75, family = FONT, fontface = "bold", label.padding = unit(1, "pt")) +
  scale_colour_manual(values = PAL_GLIM, name = "GLIM") +
  labs(x = "ASMI (sex-standardized z)", y = "Phase angle (°)") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, size = 1.75, family = FONT,
           label = spearman_lab(p5df$asmi_z, p5df$phase_angle_deg)) +
  theme(legend.position = "inside", legend.position.inside = c(1, 0), legend.justification = c(1, 0)) +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1)))

# ---- P6: intake vs ESPEN targets ---------------------------------------------------
espen <- dat %>% group_by(grp4) %>%
  summarise(across(c(protein_g_per_kg, energy_kcal_per_kg), list(m = median, lo = ~quantile(.x, 0.25), hi = ~quantile(.x, 0.75))),
            p_prot = mean(protein_g_per_kg >= 1), p_en = mean(energy_kcal_per_kg >= 30))
p6 <- ggplot(dat, aes(energy_kcal_per_kg, protein_g_per_kg)) +
  annotate("rect", xmin = 30, xmax = Inf, ymin = 1, ymax = Inf, fill = "#DCE9E3", alpha = 0.6) +
  geom_point(aes(colour = grp4), size = 0.35, alpha = 0.35, stroke = 0) +
  geom_errorbar(data = espen, aes(x = energy_kcal_per_kg_m, ymin = protein_g_per_kg_lo, ymax = protein_g_per_kg_hi,
                                  colour = grp4), inherit.aes = FALSE, width = 0, linewidth = 0.45) +
  geom_errorbar(data = espen, aes(y = protein_g_per_kg_m, xmin = energy_kcal_per_kg_lo, xmax = energy_kcal_per_kg_hi,
                                  colour = grp4), inherit.aes = FALSE, width = 0, linewidth = 0.45, orientation = "y") +
  geom_point(data = espen, aes(energy_kcal_per_kg_m, protein_g_per_kg_m, fill = grp4), shape = 21,
             size = 1.8, colour = "black", stroke = 0.3) +
  geom_vline(xintercept = 30, linetype = "22", linewidth = 0.3, colour = COL_REF) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = 0.3, colour = COL_REF) +
  annotate("text", x = 51, y = 2.25, label = "Meets ESPEN\nenergy + protein", hjust = 1, vjust = 1,
           size = 1.7, family = FONT, colour = "#2F6B55", lineheight = 0.85) +
  scale_colour_manual(values = PAL_GRP4, guide = "none") +
  scale_fill_manual(values = PAL_GRP4, guide = "none") +
  coord_cartesian(xlim = c(10, 52), ylim = c(0.3, 2.3)) +
  labs(x = "Energy intake (kcal/kg/day)", y = "Protein intake (g/kg/day)")

# ---- P7–P8: PCA ------------------------------------------------------------------------
pca_vars <- c("BMI" = "bmi", "Calf circ." = "calf_circ_cm", "Arm circ." = "arm_circ_cm",
              "Handgrip" = "handgrip_kg", "ASMI" = "asmi_kg_m2", "FFMI" = "ffmi_kg_m2",
              "Phase angle" = "phase_angle_deg", "Gait speed" = "gait_speed_m_s", "MNA-SF" = "mna_sf",
              "Albumin" = "albumin_g_dl", "Prealbumin" = "prealbumin_mg_dl", "Weight loss" = "weight_loss_6m_pct",
              "Protein/kg" = "protein_g_per_kg", "Energy/kg" = "energy_kcal_per_kg")
pdat <- dat %>% select(participant_id, grp4, stratum, hospital, sex, hcy, all_of(unname(pca_vars))) %>% drop_na(-hcy)
for (k in c("handgrip_kg", "asmi_kg_m2", "ffmi_kg_m2")) pdat[[k]] <- resid(lm(pdat[[k]] ~ pdat$sex))  # sex-standardize
pca <- PCA(pdat[, unname(pca_vars)], scale.unit = TRUE, graph = FALSE, ncp = 5)
eig <- pca$eig[, 2]
loads <- as_tibble(pca$var$coord[, 1:2], rownames = "var") %>%
  mutate(lab = names(pca_vars)[match(var, pca_vars)],
         dom = case_when(var %in% c("protein_g_per_kg", "energy_kcal_per_kg") ~ "Intake",
                         var %in% c("albumin_g_dl", "prealbumin_mg_dl") ~ "Biochemistry",
                         var %in% c("mna_sf", "weight_loss_6m_pct") ~ "Screening",
                         var %in% c("handgrip_kg", "gait_speed_m_s") ~ "Function",
                         TRUE ~ "Body composition"))
PAL_DOM <- c("Body composition" = "#2B6CB0", "Function" = "#9E2A2B", "Screening" = "#E09F3E",
             "Biochemistry" = "#6C5B7B", "Intake" = "#4C956C")
circ <- tibble(t = seq(0, 2 * pi, length.out = 200), x = cos(t), y = sin(t))
p7 <- ggplot() +
  geom_path(data = circ, aes(x, y), linewidth = 0.25, colour = "grey70") +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey80") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey80") +
  geom_segment(data = loads, aes(0, 0, xend = Dim.1, yend = Dim.2, colour = dom), linewidth = 0.4,
               arrow = arrow(length = unit(2.2, "pt"), type = "closed")) +
  ggrepel::geom_text_repel(data = loads, aes(Dim.1, Dim.2, label = lab, colour = dom), size = 1.75,
                           family = FONT, segment.size = 0.15, min.segment.length = 0.2, box.padding = 0.15,
                           max.overlaps = Inf, seed = 2) +
  scale_colour_manual(values = PAL_DOM, name = NULL) +
  coord_equal(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05)) +
  labs(x = sprintf("PC1 (%.1f%%)", eig[1]), y = sprintf("PC2 (%.1f%%)", eig[2])) +
  theme(legend.position = "bottom", legend.key.size = unit(5, "pt"), legend.text = element_text(size = BASE_PT - 1.5)) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(label = "■", size = 2)))

scores <- as_tibble(pca$ind$coord[, 1:4]) %>% bind_cols(pdat %>% select(participant_id, grp4, stratum, hospital, hcy))
cent <- scores %>% group_by(grp4) %>% summarise(across(Dim.1:Dim.2, mean))
p8 <- ggplot(scores, aes(Dim.1, Dim.2, colour = grp4)) +
  geom_point(size = 0.4, alpha = 0.4, stroke = 0) +
  stat_ellipse(level = 0.68, linewidth = 0.45) +
  geom_path(data = cent, aes(group = 1), colour = "black", linewidth = 0.35,
            arrow = arrow(length = unit(2.5, "pt"), type = "closed")) +
  geom_point(data = cent, aes(fill = grp4), shape = 21, size = 2, colour = "black", stroke = 0.35) +
  scale_colour_manual(values = PAL_GRP4, name = NULL) + scale_fill_manual(values = PAL_GRP4, guide = "none") +
  labs(x = sprintf("PC1 (%.1f%%)", eig[1]), y = sprintf("PC2 (%.1f%%)", eig[2])) +
  theme(legend.position = "none")

# ---- P9–P10: Gaussian mixture phenotypes ---------------------------------------------
gmm <- Mclust(scores %>% select(Dim.1:Dim.4), G = 2:6, verbose = FALSE)
scores$cl <- factor(gmm$classification)
prof <- pdat %>% select(all_of(unname(pca_vars))) %>% mutate(across(everything(), ~ as.numeric(scale(.x)))) %>%
  mutate(cl = scores$cl) %>% group_by(cl) %>% summarise(across(everything(), mean), n = n())
# Data-driven naming: muscle (ASMI) x recent weight loss
nm <- prof %>% transmute(cl, muscle = asmi_kg_m2 > 0, wl = weight_loss_6m_pct > 0.25,
                         name = case_when(muscle & !wl ~ "Robust", muscle & wl ~ "Acute weight loss",
                                          !muscle & !wl ~ "Sarcopenic", TRUE ~ "Cachectic"))
if (any(duplicated(nm$name))) nm$name <- paste0("Phenotype ", nm$cl)
PH_LEVELS <- intersect(c("Robust", "Acute weight loss", "Sarcopenic", "Cachectic"), nm$name)
if (length(PH_LEVELS) == 0) PH_LEVELS <- nm$name
PAL_PH <- setNames(c("#4C956C", "#E9C46A", "#E76F51", "#6D2E46", "#577590", "#999999")[seq_along(PH_LEVELS)], PH_LEVELS)
scores$phenotype <- factor(nm$name[match(scores$cl, nm$cl)], PH_LEVELS)
write.csv(scores %>% select(participant_id, phenotype), file.path(OUT_TAB, "nutritional_phenotype_clusters.csv"), row.names = FALSE)

prof_long <- prof %>% mutate(phenotype = factor(nm$name[match(cl, nm$cl)], PH_LEVELS)) %>%
  pivot_longer(all_of(unname(pca_vars)), names_to = "var", values_to = "z") %>%
  mutate(lab = factor(names(pca_vars)[match(var, pca_vars)], rev(names(pca_vars))))
p9 <- ggplot(prof_long, aes(phenotype, lab, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", z), colour = abs(z) > 0.8), size = 1.6, family = FONT) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_gradient2(low = PAL_DIV["low"], mid = PAL_DIV["mid"], high = PAL_DIV["high"], limits = c(-1.5, 1.5),
                       oob = squish, name = "z") +
  scale_x_discrete(NULL, labels = function(x) {
    nn <- prof$n[match(x, nm$name[match(prof$cl, nm$cl)])]; sprintf("%s (%d)", x, nn) }) +
  scale_y_discrete(NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        legend.key.height = unit(12, "pt"), legend.key.width = unit(4, "pt"))

comp <- scores %>% count(grp4, phenotype) %>% group_by(grp4) %>% mutate(p = n / sum(n))
p10 <- ggplot(comp, aes(grp4, p, fill = fct_rev(phenotype))) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = if_else(p > 0.07, as.character(round(100 * p)), "")), position = position_stack(vjust = 0.5),
            size = 1.75, family = FONT, colour = "white") +
  scale_fill_manual(values = PAL_PH, name = NULL, breaks = PH_LEVELS) +
  scale_y_continuous("Participants (%)", labels = function(x) 100 * x, expand = c(0, 0)) +
  scale_x_discrete(NULL, labels = lab_grp4) + rot_x +
  theme(legend.position = "right", legend.key.size = unit(6, "pt"))

# ---- P11: Hcy by phenotype within hospital strata ---------------------------------------
p11df <- scores %>% filter(hospital, !is.na(hcy))
fit_ph <- lm(log(hcy) ~ phenotype + stratum, data = p11df)
p_ph <- anova(lm(log(hcy) ~ stratum, data = p11df), fit_ph)$`Pr(>F)`[2]
p11 <- ggplot(p11df, aes(stratum, hcy, colour = phenotype, fill = phenotype)) +
  geom_hline(yintercept = 15, linetype = "22", linewidth = 0.3, colour = COL_REF) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.78, seed = 3),
             size = 0.45, alpha = 0.45, stroke = 0) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), position = position_dodge(0.78),
               geom = "pointrange", size = 0.18, linewidth = 0.4, shape = 21, colour = "black") +
  scale_colour_manual(values = PAL_PH, guide = "none") + scale_fill_manual(values = PAL_PH, guide = "none") +
  scale_y_log10(lab_umol, breaks = c(5, 10, 15, 20, 30, 40)) +
  scale_x_discrete("Hospital stratum", labels = function(x) str_remove(x, " y$")) +
  annotate("text", x = 0.5, y = Inf, hjust = 0, vjust = 1.2, size = 1.75, family = FONT,
           label = sprintf("Phenotype (adj. stratum)\n%s", fmt_p(p_ph)), lineheight = 0.9)

# ---- Assemble -------------------------------------------------------------------------
top <- (wrap_elements(full = p1) | p2 | (p3 / p4)) + plot_layout(widths = c(1, 1.15, 1))
mid <- (p5 | p6 | p7) + plot_layout(widths = c(1.25, 0.9, 1))
bot <- (wrap_elements(full = p8) | p9 | p10 | wrap_elements(full = p11)) + plot_layout(widths = c(0.95, 0.95, 0.8, 1.05))
fig2 <- wrap_elements(full = top) / wrap_elements(full = mid) / wrap_elements(full = bot) +
  plot_layout(heights = c(1.25, 1, 0.95))
save_fig(fig2, "Figure2_nutritional_phenotype", width = 180, height = 215)
save_panels(list(prevalence = p1, alluvial = p2, phase_angle = p3, handgrip = p4, pa_asmi = p5,
                 espen = p6, pca_loadings = p7, pca_scores = p8, cluster_profile = p9,
                 cluster_composition = p10, hcy_by_cluster = p11), "Figure2")
