# =============================================================================
# Figure 7 — Integration: nutrition, one-carbon metabolism, methylation and ageing markers;
#            clinical translation (Q6; H1e)
#   P1  Partial-correlation network (EBIC graphical lasso; edge stability by bootstrap)
#   P2  Clustered heatmap: one-carbon metabolites x nutritional markers (partial Spearman, adj. age, sex, eGFR)
#   P3  Structural equation model: malnutrition -> B12 -> Hcy -> SAM/SAH -> 5-mC; Hcy -> ADMA / GDF-15
#   P4  Hcy vs downstream markers (SAM/SAH, 5-mC, ADMA, GDF-15)
#   P5  Conditional inference tree for Hcy > 15 umol/L (clinical triage)
#   P6  "Etiological" profile of hyperhomocysteinaemia by stratum (renal / B12 / malnutrition)
#   P7  Cross-validated discrimination (AUC) of nested determinant sets for Hcy > 15
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages({ library(qgraph); library(igraph); library(ggraph); library(tidygraph)
  library(lavaan); library(partykit); library(pROC) })
select <- dplyr::select; map <- purrr::map

h <- hosp %>% mutate(male = as.integer(sex == "Male"),
                     asmi_z = ave(asmi_kg_m2, sex, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)),
                     grip_z = ave(handgrip_kg, sex, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)))

# ---- P1: network ----------------------------------------------------------------------------
NODES <- tribble(~var, ~lab, ~dom,
  "hcy", "Hcy", "One-carbon", "sam_sah", "SAM/SAH", "One-carbon", "cystathionine_nmol_l", "Cysta", "One-carbon",
  "betaine_umol_l", "Betaine", "One-carbon", "dmg_umol_l", "DMG", "One-carbon", "methionine_umol_l", "Met", "One-carbon",
  "holotc", "holoTC", "B vitamins", "rbc_fol", "RBC fol", "B vitamins", "plp_nmol_l", "PLP", "B vitamins",
  "egrac", "EGRAC", "B vitamins", "mma_umol_l", "MMA", "B vitamins",
  "egfr", "eGFR", "Kidney",
  "phase_angle_deg", "PhA", "Nutrition", "asmi_z", "ASMI", "Nutrition", "grip_z", "Grip", "Nutrition",
  "mna_sf", "MNA-SF", "Nutrition", "albumin_g_dl", "Alb", "Nutrition",
  "crp_mg_l", "CRP", "Inflammation", "il6_pg_ml", "IL-6", "Inflammation",
  "mc5", "5-mC", "Ageing / vascular", "adma", "ADMA", "Ageing / vascular", "gdf15", "GDF-15", "Ageing / vascular")
PAL_NDOM <- c("One-carbon" = "#9E2A2B", "B vitamins" = "#E09F3E", "Kidney" = "#2B6CB0", "Nutrition" = "#4C956C",
              "Inflammation" = "#6C5B7B", "Ageing / vascular" = "#3A3A3A")
nd <- h %>% select(age_years, male, all_of(NODES$var)) %>% drop_na()
Zc <- nd %>% select(age_years, male)
R <- sapply(NODES$var, function(v) rank_resid(nd[[v]], Zc))               # rank-based, adjusted for age and sex
glasso_fit <- function(Rm) EBICglasso(cor(Rm), n = nrow(Rm), gamma = 0.5, verbose = FALSE)
W <- glasso_fit(R)
B_NET <- max(100, N_BOOT %/% 5)
incl <- Reduce(`+`, map(seq_len(B_NET), function(b) {
  i <- sample.int(nrow(R), replace = TRUE); (abs(glasso_fit(R[i, ])) > 0) * 1 })) / B_NET
dimnames(W) <- dimnames(incl) <- list(NODES$lab, NODES$lab)
ed <- which(upper.tri(W) & W != 0, arr.ind = TRUE)
edges <- tibble(from = rownames(W)[ed[, 1]], to = colnames(W)[ed[, 2]], w = W[ed], stab = incl[ed]) %>%
  filter(abs(w) >= 0.05)
ord <- NODES %>% mutate(dom = factor(dom, names(PAL_NDOM))) %>% arrange(dom)
gr <- tbl_graph(nodes = ord %>% transmute(name = lab, dom = as.character(dom)), edges = edges, directed = FALSE) %>%
  mutate(strength = centrality_degree(weights = abs(w)))
lay <- create_layout(gr, layout = "linear", circular = TRUE)
labdf <- tibble(name = lay$name, dom = lay$dom, ang = atan2(lay$y, lay$x)) %>%
  mutate(lx = 1.13 * cos(ang), ly = 1.13 * sin(ang), hj = if_else(cos(ang) >= -1e-6, 0, 1))
p1 <- ggraph(lay) +
  geom_edge_arc(aes(edge_width = abs(w), edge_colour = w > 0, edge_alpha = stab), strength = 0.15) +
  geom_node_point(aes(fill = dom, size = strength), shape = 21, colour = "white", stroke = 0.4) +
  geom_text(data = labdf, aes(lx, ly, label = name, hjust = hj, colour = dom), size = 1.75, family = FONT, fontface = "bold") +
  scale_colour_manual(values = PAL_NDOM, guide = "none") +
  scale_edge_width(range = c(0.25, 2), guide = "none") +
  scale_edge_alpha(range = c(0.2, 1), limits = c(0.5, 1), name = "Bootstrap\ninclusion", oob = squish,
                   breaks = c(0.5, 0.75, 1)) +
  scale_edge_colour_manual(values = c(`TRUE` = "#9E2A2B", `FALSE` = "#2B6CB0"), labels = c(`TRUE` = "Positive", `FALSE` = "Negative"),
                           name = NULL) +
  scale_fill_manual(values = PAL_NDOM, name = NULL) +
  scale_size(range = c(1.6, 4), guide = "none") +
  coord_equal(xlim = c(-1.45, 1.45), ylim = c(-1.25, 1.25), clip = "off") +
  theme_void(base_family = FONT, base_size = BASE_PT) +
  theme(legend.position = "right", legend.key.size = unit(5, "pt"), legend.text = element_text(size = BASE_PT - 1.5),
        legend.title = element_text(size = BASE_PT - 1.5, face = "bold")) +
  guides(fill = guide_legend(override.aes = list(size = 2.2), order = 1), edge_colour = guide_legend(order = 2))

# ---- P2: clustered partial-Spearman heatmap ------------------------------------------------------
ROWS <- c("Hcy" = "hcy", "SAM" = "sam_nmol_l", "SAH" = "sah_nmol_l", "SAM/SAH" = "sam_sah", "Methionine" = "methionine_umol_l",
          "Cystathionine" = "cystathionine_nmol_l", "Cysteine" = "total_cysteine_umol_l", "Betaine" = "betaine_umol_l",
          "Choline" = "choline_plasma_umol_l", "DMG" = "dmg_umol_l", "5-MTHF" = "mthf5_plasma_nmol_l", "MMA" = "mma_umol_l",
          "5-mC" = "mc5", "ADMA" = "adma", "GDF-15" = "gdf15")
COLS <- c("Phase angle" = "phase_angle_deg", "ASMI (z)" = "asmi_z", "FFMI" = "ffmi_kg_m2", "Handgrip (z)" = "grip_z",
          "Gait speed" = "gait_speed_m_s", "Calf circ." = "calf_circ_cm", "BMI" = "bmi", "MNA-SF" = "mna_sf",
          "Weight loss" = "weight_loss_6m_pct", "Albumin" = "albumin_g_dl", "Prealbumin" = "prealbumin_mg_dl",
          "Protein/kg" = "protein_g_per_kg", "CRP" = "crp_mg_l")
Zh <- h %>% select(age_years, male, log_egfr)
pc <- expand_grid(r = names(ROWS), c = names(COLS)) %>% rowwise() %>% mutate({
  x <- rank_resid(h[[ROWS[[r]]]], Zh); y <- rank_resid(h[[COLS[[c]]]], Zh)
  ok <- complete.cases(x, y); ct <- cor.test(x[ok], y[ok])
  tibble(rho = unname(ct$estimate), p = ct$p.value) }) %>% ungroup() %>% mutate(q = p.adjust(p, "BH"))
M <- pc %>% select(r, c, rho) %>% pivot_wider(names_from = c, values_from = rho) %>% tibble::column_to_rownames("r") %>% as.matrix()
ro <- hclust(dist(M), "average")$order; co <- hclust(dist(t(M)), "average")$order
pc <- pc %>% mutate(r = factor(r, rownames(M)[ro]), c = factor(c, colnames(M)[co]))
p2 <- ggplot(pc, aes(c, r, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = case_when(q < 0.001 ~ "***", q < 0.01 ~ "**", q < 0.05 ~ "*", TRUE ~ "")),
            size = 1.7, vjust = 0.75, family = FONT) +
  scale_fill_gradient2(low = PAL_DIV["low"], mid = "white", high = PAL_DIV["high"], limits = c(-0.4, 0.4),
                       oob = squish, name = "Partial\nρ") +
  scale_x_discrete(NULL, position = "top") + scale_y_discrete(NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(angle = 50, hjust = 0, vjust = 0, size = BASE_PT - 1.5),
        axis.text.y = element_text(size = BASE_PT - 1.5),
        legend.key.height = unit(12, "pt"), legend.key.width = unit(4, "pt"))

# ---- P3: SEM -------------------------------------------------------------------------------------------
sd <- h %>% transmute(age = age_c10, male, egfr = log_egfr, pha = phase_angle_deg, mal = glim_any,
                      b12 = log_holotc, fol = log_rbcf, hcy = log_hcy, ssah = log(sam_sah), mc5,
                      adma, gdf = log(gdf15)) %>% drop_na() %>% mutate(across(everything(), ~ as.numeric(scale(.x))))
sem_model <- '
  pha  ~ age + male
  mal  ~ age + pha
  egfr ~ age
  b12  ~ age + mal + egfr
  fol  ~ age
  hcy  ~ egfr + b12 + fol + pha + mal + age + male
  ssah ~ hcy + egfr
  mc5  ~ ssah + age
  adma ~ hcy + egfr
  gdf  ~ hcy + egfr + age + pha
'
sem_fit <- sem(sem_model, data = sd, estimator = "MLR")
fm <- fitMeasures(sem_fit, c("cfi.robust", "rmsea.robust", "srmr"))
pe <- standardizedSolution(sem_fit) %>% filter(op == "~") %>% as_tibble()
write.csv(pe, file.path(OUT_TAB, "fig7_sem_standardized_paths.csv"), row.names = FALSE)
SN <- tribble(~id, ~lab, ~x, ~y, ~dom,
  "pha", "Phase\nangle", 0, 2.7, "Nutrition", "mal", "GLIM", 0, 1.55, "Nutrition",
  "fol", "RBC\nfolate", 0, 0.4, "B vitamins", "egfr", "eGFR", 1.75, -0.45, "Kidney",
  "b12", "holoTC", 1.5, 2.25, "B vitamins", "hcy", "Hcy", 3.0, 1.2, "One-carbon",
  "ssah", "SAM/\nSAH", 4.5, 2.55, "One-carbon", "mc5", "5-mC", 6.0, 2.55, "Ageing / vascular",
  "adma", "ADMA", 4.5, 0.6, "Ageing / vascular", "gdf", "GDF-15", 6.0, 0.6, "Ageing / vascular")
se <- pe %>% filter(lhs %in% SN$id, rhs %in% SN$id) %>%
  left_join(SN %>% select(id, x0 = x, y0 = y), by = c("rhs" = "id")) %>%
  left_join(SN %>% select(id, x1 = x, y1 = y), by = c("lhs" = "id")) %>%
  mutate(sig = pvalue < 0.05, dx = x1 - x0, dy = (y1 - y0) * 1.6, L = sqrt(dx^2 + dy^2),
         xs = x0 + dx / L * 0.36, ys = y0 + dy / L * 0.36 / 1.6, xe = x1 - dx / L * 0.40, ye = y1 - dy / L * 0.40 / 1.6,
         xm = xs + (xe - xs) * 0.5, ym = ys + (ye - ys) * 0.5) %>%
  filter(sig)
age_e <- pe %>% filter(rhs == "age", pvalue < 0.05, lhs %in% SN$id) %>%
  left_join(SN %>% select(id, x, y), by = c("lhs" = "id"))
p3 <- ggplot() +
  geom_segment(data = se, aes(xs, ys, xend = xe, yend = ye, linewidth = abs(est.std), colour = est.std > 0),
               arrow = arrow(length = unit(2.6, "pt"), type = "closed"), lineend = "butt") +
  geom_label(data = se, aes(xm, ym, label = sprintf("%.2f", est.std), colour = est.std > 0), size = 1.55, family = FONT,
             linewidth = 0, label.padding = unit(0.7, "pt"), fill = alpha("white", 0.9)) +
  geom_point(data = SN, aes(x, y, fill = dom), shape = 21, size = 9.5, colour = "white", stroke = 0.4) +
  geom_text(data = SN, aes(x, y, label = lab), size = 1.6, family = FONT, lineheight = 0.8, colour = "white",
            fontface = "bold") +
  geom_text(data = age_e, aes(x, y + 0.42, label = sprintf("age %+.2f", est.std)), size = 1.4, family = FONT,
            colour = "grey35") +
  scale_linewidth(range = c(0.25, 1.6), guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "#9E2A2B", `FALSE` = "#2B6CB0"), guide = "none") +
  scale_fill_manual(values = PAL_NDOM, guide = "none") +
  coord_cartesian(xlim = c(-0.45, 6.45), ylim = c(-0.8, 3.25)) +
  labs(subtitle = sprintf("Standardized paths with P < 0.05 (age, sex adjusted; grey = age paths)   CFI %.2f • RMSEA %.3f • SRMR %.3f",
                          fm[1], fm[2], fm[3])) +
  theme_void(base_family = FONT, base_size = BASE_PT) +
  theme(plot.subtitle = element_text(size = BASE_PT - 1.5, colour = "grey25"))

# ---- P4: Hcy vs downstream markers ---------------------------------------------------------------------------
DS <- c("SAM/SAH" = "sam_sah", "5-mC (%)" = "mc5", "ADMA (µmol/L)" = "adma", "GDF-15 (pg/mL)" = "gdf15")
p4df <- h %>% select(stratum, hcy, age_years, male, log_egfr, all_of(unname(DS))) %>%
  pivot_longer(all_of(unname(DS)), names_to = "m", values_to = "v") %>% filter(!is.na(v)) %>%
  mutate(m = factor(names(DS)[match(m, DS)], names(DS)))
p4lab <- p4df %>% group_by(m) %>% group_modify(function(d, k) {
  Z <- d %>% select(age_years, male, log_egfr)
  r <- cor.test(rank_resid(d$hcy, Z), rank_resid(d$v, Z))
  tibble(lab = sprintf("partial ρ = %.2f\n%s", r$estimate, fmt_p(r$p.value)))
})
p4 <- ggplot(p4df, aes(hcy, v)) +
  geom_point(aes(colour = stratum), size = 0.4, alpha = 0.5, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.5, fill = "grey70", alpha = 0.35) +
  geom_text(data = p4lab, aes(x = Inf, y = Inf, label = lab), hjust = 1.02, vjust = 1.15, size = 1.55,
            family = FONT, lineheight = 0.85) +
  facet_wrap(~m, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = PAL_STR, guide = "none") +
  scale_x_log10(lab_umol, breaks = c(8, 15, 30)) + scale_y_continuous(NULL)

# ---- P5: conditional inference tree --------------------------------------------------------------------------
td <- h %>% filter(!is.na(hcy)) %>%
  transmute(HHcy = factor(if_else(hcy > 15, "Hcy > 15", "Hcy ≤ 15"), c("Hcy ≤ 15", "Hcy > 15")),
            `eGFR-cys` = egfr, holoTC = holotc, `RBC folate` = rbc_fol, PLP = plp_nmol_l, EGRAC = egrac,
            `Phase angle` = phase_angle_deg, `MNA-SF` = mna_sf, GLIM = factor(glim_any, 0:1, c("No", "Yes")),
            CRP = crp_mg_l, Stratum = stratum, `MTHFR TT` = factor(as.integer(mthfr677 == "TT"), 0:1, c("No", "Yes")))
ct <- ctree(HHcy ~ ., data = td, control = ctree_control(maxdepth = 3, minbucket = 25, alpha = 0.01))
saveRDS(ct, file.path(OUT_TAB, "fig7_ctree.rds"))
# Draw the tree in ggplot (layout: leaves evenly spaced, parents centred over children)
tree_df <- function(ct) {
  dd <- data_party(ct); resp <- dd[["(response)"]]; fit <- dd[["(fitted)"]]
  vnames <- names(dd)
  nodes <- list(); edges <- list(); leaf_i <- 0
  walk <- function(nd, depth) {
    id <- id_node(nd)
    if (is.terminal(nd)) {
      leaf_i <<- leaf_i + 1
      y <- resp[fit == id]
      nodes[[length(nodes) + 1]] <<- tibble(id = id, depth = depth, x = leaf_i, terminal = TRUE,
                                            lab = NA_character_, n = length(y), p = mean(y == levels(resp)[2]))
      return(leaf_i)
    }
    sp <- split_node(nd); v <- vnames[varid_split(sp)]
    pv <- info_node(nd)$p.value
    kids <- kids_node(nd)
    xs <- map_dbl(kids, ~ walk(.x, depth + 1))
    br <- breaks_split(sp)
    elab <- if (!is.null(br)) c(paste("≤", signif(br, 3)), paste(">", signif(br, 3))) else {
      ix <- index_split(sp); lv <- levels(dd[[v]]); map_chr(seq_along(kids), ~ paste(lv[ix == .x], collapse = ", "))
    }
    x <- mean(xs)
    nodes[[length(nodes) + 1]] <<- tibble(id = id, depth = depth, x = x, terminal = FALSE,
                                          lab = sprintf("%s\n%s", v, fmt_p(pv)), n = NA, p = NA)
    for (k in seq_along(kids)) edges[[length(edges) + 1]] <<- tibble(x = x, y = depth, xend = xs[k], yend = depth + 1, lab = elab[k])
    x
  }
  walk(node_party(ct), 0)
  list(nodes = bind_rows(nodes), edges = bind_rows(edges))
}
tdf <- tree_df(ct)
leaves <- tdf$nodes %>% filter(terminal) %>% mutate(dmax = max(tdf$nodes$depth))
leaves <- leaves %>% mutate(depth = dmax)       # align leaves at the bottom
tdf$edges <- tdf$edges %>% left_join(leaves %>% select(xend = x, dl = depth), by = "xend") %>%
  mutate(yend = coalesce(dl, yend))
BW <- 0.32
p5 <- ggplot() +
  geom_segment(data = tdf$edges, aes(x, -y, xend = xend, yend = -yend + if_else(is.na(dl), 0, 0.32)), linewidth = 0.3, colour = "grey45") +
  geom_label(data = tdf$edges, aes((x + xend) / 2, -(y + yend) / 2 + 0.1, label = lab), size = 1.55, family = FONT,
             linewidth = 0, fill = "white", label.padding = unit(0.6, "pt")) +
  geom_label(data = filter(tdf$nodes, !terminal), aes(x, -depth, label = lab), size = 1.65, family = FONT,
             fill = "#F4F1EA", linewidth = 0.25, label.r = unit(3, "pt"), lineheight = 0.85) +
  geom_rect(data = leaves, aes(xmin = x - BW, xmax = x + BW, ymin = -depth - 0.25, ymax = -depth + 0.3),
            fill = "#E3ECE8", colour = "grey40", linewidth = 0.25) +
  geom_rect(data = leaves, aes(xmin = x - BW, xmax = x - BW + 2 * BW * p, ymin = -depth - 0.25, ymax = -depth + 0.3),
            fill = "#9E2A2B") +
  geom_text(data = leaves, aes(x, -depth - 0.42, label = sprintf("%.0f%% > 15\nn = %d", 100 * p, n)),
            size = 1.55, family = FONT, lineheight = 0.85, vjust = 1) +
  scale_x_continuous(expand = expansion(add = 0.45)) +
  scale_y_continuous(expand = expansion(add = c(0.45, 0.25))) +
  theme_void(base_family = FONT, base_size = BASE_PT)

# ---- P6: etiological profile of HHcy -------------------------------------------------------------------------
et <- h %>% filter(hhcy == 1) %>%
  mutate(renal = egfr < 60, b12 = holotc < 35, mal = glim_any == 1,
         prof = case_when(renal & !b12 & !mal ~ "Renal only", !renal & b12 & !mal ~ "B12 only",
                          !renal & !b12 & mal ~ "Malnutrition only", renal & b12 & mal ~ "All three",
                          renal & mal ~ "Renal + malnutrition", renal & b12 ~ "Renal + B12",
                          b12 & mal ~ "B12 + malnutrition", TRUE ~ "None of the three")) %>%
  count(stratum, prof) %>% group_by(stratum) %>% mutate(p = n / sum(n), N = sum(n)) %>% ungroup()
PROF_LEV <- c("Renal only", "Renal + malnutrition", "Renal + B12", "All three", "B12 + malnutrition", "B12 only",
              "Malnutrition only", "None of the three")
PAL_PROF <- setNames(c("#2B6CB0", "#5B8FC7", "#6C5B7B", "#3A2E39", "#B85C5C", "#9E2A2B", "#E09F3E", "#D6D6D6"), PROF_LEV)
et <- et %>% mutate(prof = factor(prof, rev(PROF_LEV)))
p6 <- ggplot(et, aes(stratum, p, fill = prof)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = if_else(p >= 0.06, sprintf("%.0f", 100 * p), "")), position = position_stack(vjust = 0.5),
            size = 1.6, family = FONT, colour = "white") +
  geom_text(data = distinct(et, stratum, N), aes(stratum, 1.02, label = paste0("n=", N)), inherit.aes = FALSE,
            vjust = 0, size = 1.55, family = FONT, colour = "grey30") +
  scale_fill_manual(values = PAL_PROF, breaks = PROF_LEV, name = "Among Hcy > 15:") +
  scale_y_continuous("Patients with Hcy > 15 (%)", labels = function(x) 100 * x, expand = expansion(mult = c(0, 0.07))) +
  scale_x_discrete("Hospital stratum", labels = function(x) str_remove(x, " y$")) +
  theme(legend.key.size = unit(5, "pt"), legend.text = element_text(size = BASE_PT - 1.5),
        legend.position = "bottom", legend.title.position = "top") +
  guides(fill = guide_legend(ncol = 2))

# ---- P7: cross-validated AUC of nested sets --------------------------------------------------------------------
cvd <- h %>% select(hhcy, stratum, male, log_egfr, log_holotc, log_rbcf, log_plp, egrac, glim_any, phase_angle_deg,
                    mna_sf, log_crp, hcy_genetic_risk_score, tt = mthfr677) %>% drop_na() %>%
  mutate(tt = as.integer(tt == "TT"))
SETS <- list("Age + sex" = c("stratum", "male"),
             "+ eGFR" = c("stratum", "male", "log_egfr"),
             "+ holoTC" = c("stratum", "male", "log_egfr", "log_holotc"),
             "+ nutrition (GLIM, PhA, MNA)" = c("stratum", "male", "log_egfr", "log_holotc", "glim_any", "phase_angle_deg", "mna_sf"),
             "+ folate, B6, B2" = c("stratum", "male", "log_egfr", "log_holotc", "glim_any", "phase_angle_deg", "mna_sf", "log_rbcf", "log_plp", "egrac"),
             "+ genetics (GRS, 677TT)" = c("stratum", "male", "log_egfr", "log_holotc", "glim_any", "phase_angle_deg", "mna_sf", "log_rbcf", "log_plp", "egrac", "hcy_genetic_risk_score", "tt"))
K <- 10; REP <- 20
cv_auc <- map_dfr(names(SETS), function(nm) {
  map_dfr(seq_len(REP), function(r) {
    fold <- sample(rep(1:K, length.out = nrow(cvd)))
    pr <- numeric(nrow(cvd))
    for (k in 1:K) {
      f <- glm(reformulate(SETS[[nm]], "hhcy"), data = cvd[fold != k, ], family = binomial)
      pr[fold == k] <- predict(f, cvd[fold == k, ], type = "response")
    }
    tibble(set = nm, rep = r, auc = as.numeric(auc(roc(cvd$hhcy, pr, quiet = TRUE))))
  })
}) %>% group_by(set) %>% summarise(m = mean(auc), lo = quantile(auc, 0.025), hi = quantile(auc, 0.975)) %>%
  mutate(set = factor(set, rev(names(SETS))))
write.csv(cv_auc, file.path(OUT_TAB, "fig7_cv_auc.csv"), row.names = FALSE)
p7 <- ggplot(cv_auc, aes(m, set)) +
  geom_vline(xintercept = 0.5, linewidth = 0.3, colour = COL_REF) +
  geom_segment(aes(x = 0.5, xend = m, yend = set), linewidth = 2.2, colour = "grey88") +
  geom_pointrange(aes(xmin = lo, xmax = hi, colour = m), size = 0.25, linewidth = 0.5) +
  geom_text(aes(x = hi, label = sprintf("%.2f", m)), hjust = -0.3, size = 1.65, family = FONT) +
  scale_colour_gradientn(colours = c("#9A9A9A", "#6C5B7B", "#9E2A2B"), guide = "none") +
  scale_x_continuous("CV AUC (Hcy > 15)", breaks = c(0.5, 0.7, 0.9)) +
  coord_cartesian(xlim = c(0.5, 1)) +
  scale_y_discrete(NULL)

# ---- Assemble ----------------------------------------------------------------------------------------------------
r1 <- row_of(p1, p2, widths = c(1.1, 1))
r2 <- row_of(p3, widths = 1)
r3 <- wrap_elements(full = p4)
r4 <- row_of(p5, p6, p7, widths = c(1.2, 0.9, 0.9))
fig7 <- wrap_elements(full = r1) / wrap_elements(full = p3) / r3 / wrap_elements(full = r4) +
  plot_layout(heights = c(1.1, 0.75, 0.48, 0.95))
save_fig(fig7, "Figure7_integration_translation", width = 180, height = 235)
save_panels(list(network = p1, heatmap = p2, sem = p3, downstream = p4, ctree = p5, etiology = p6, auc = p7), "Figure7")
