# =============================================================================
# Figure 1 — Study design, analytic sets and causal framework
#   P1  Participant flow: 4 groups -> analytic sets (hospital H1 set, GxN set)
#   P2  Directed acyclic graph (DAG) encoding H1 and the adjustment set
#   P3  Age distribution by group, with analytic strata cut-points
#   P4  Data availability matrix (assessment domain x group)
#   P5  Lifetime exposure to mandatory folic-acid fortification by group
# =============================================================================
source("R/00_setup.R")
suppressPackageStartupMessages(library(ggridges))

# ---- P1: flow diagram ---------------------------------------------------------
n_grp  <- dat %>% count(grp4)
n_hcy  <- dat %>% filter(!is.na(hcy)) %>% count(grp4, name = "n_hcy")
n_gxn  <- dat %>% filter(analysis_set_gxn == 1) %>% count(grp4, name = "n_gxn")
flow_n <- n_grp %>% left_join(n_hcy, by = "grp4") %>% left_join(n_gxn, by = "grp4") %>%
  mutate(x = 1:4, excl = n - n_hcy)

box <- function(x, y, w, h, fill, col = NA) annotate("rect", xmin = x - w/2, xmax = x + w/2,
                                                    ymin = y - h/2, ymax = y + h/2,
                                                    fill = fill, colour = col, linewidth = 0.3)
txt <- function(x, y, lab, size = 2.1, ...) annotate("text", x = x, y = y, label = lab,
                                                      size = size, family = FONT, lineheight = 0.9, ...)
n_tot <- nrow(dat)
p1 <- ggplot() +
  # Recruitment sources
  box(1, 6.1, 0.9, 0.62, "grey95", "grey40") +
  txt(1, 6.1, "ELSA-Brasil cohort\n(community, visit 3)") +
  box(3, 6.1, 2.9, 0.62, "grey95", "grey40") +
  txt(3, 6.1, sprintf("Hospital admissions, 3 centres (SP, RJ, MG, ES)\nconsecutive adults ≥56 y, 2016–2026")) +
  annotate("segment", x = 1, xend = 1, y = 5.79, yend = 5.18, linewidth = 0.3,
           arrow = arrow(length = unit(2, "pt"), type = "closed")) +
  annotate("segment", x = c(2, 3, 4), xend = c(2, 3, 4), y = 5.79, yend = 5.18, linewidth = 0.3,
           arrow = arrow(length = unit(2, "pt"), type = "closed")) +
  # Four groups
  { map(1:4, function(i) box(i, 4.85, 0.9, 0.62, alpha(PAL_GRP4[i], 0.18), PAL_GRP4[i])) } +
  { map(1:4, function(i) txt(i, 4.85, sprintf("%s\nn = %d", str_replace(GRP4_LEVELS[i], " ", "\n"), flow_n$n[i]),
                             fontface = "bold", colour = colorspace::darken(PAL_GRP4[i], 0.35))) } +
  annotate("segment", x = 1:4, xend = 1:4, y = 4.54, yend = 3.88, linewidth = 0.3,
           arrow = arrow(length = unit(2, "pt"), type = "closed")) +
  # Exclusions
  { map(1:4, function(i) txt(i + 0.05, 4.21, sprintf("Hcy not valid\n(>120 min to\ncentrifugation): %d", flow_n$excl[i]),
                             size = 1.6, colour = "grey35", hjust = 0)) } +
  # Measurements
  { map(1:4, function(i) box(i, 3.55, 0.9, 0.62, "white", "grey40")) } +
  { map(1:4, function(i) txt(i, 3.55, sprintf("Homocysteine\nn = %d", flow_n$n_hcy[i]))) } +
  # Analytic sets
  annotate("rect", xmin = 1.52, xmax = 4.48, ymin = 2.25, ymax = 2.95, fill = "#FBF3E4", colour = "#E09F3E", linewidth = 0.35) +
  txt(3, 2.6, sprintf("Primary H1 set (hospital, 3 age strata)\nn = %d", sum(flow_n$n_hcy[2:4])), fontface = "bold") +
  annotate("rect", xmin = 0.52, xmax = 4.48, ymin = 1.35, ymax = 2.05, fill = "#EAF1F8", colour = "#2B6CB0", linewidth = 0.35) +
  txt(2.5, 1.7, sprintf("Gene × nutrient set (all groups; Hcy + 10 core genotypes)\nn = %d", sum(flow_n$n_gxn)), fontface = "bold") +
  annotate("segment", x = 2:4, xend = 2:4, y = 3.24, yend = 2.95, linewidth = 0.3,
           arrow = arrow(length = unit(2, "pt"), type = "closed")) +
  annotate("segment", x = 1, xend = 1, y = 3.24, yend = 2.05, linewidth = 0.3,
           arrow = arrow(length = unit(2, "pt"), type = "closed")) +
  # Assessments strip
  txt(2.5, 0.9, "FFQ (\u226472 h)  \u2022  Anthropometry, BIA, handgrip  \u2022  MNA-SF, GLIM, EWGSOP2\nFasting blood (HPLC, LC-MS/MS, ELISA)  \u2022  DNA (TaqMan, 24 one-carbon SNPs)  \u2022  Vital status",
      size = 1.8, colour = "grey20") +
  coord_cartesian(xlim = c(0.5, 4.6), ylim = c(0.5, 6.45), expand = FALSE) +
  theme_void(base_family = FONT)

# ---- P2: DAG ------------------------------------------------------------------
nodes <- tribble(
  ~id,    ~label,                          ~x,   ~y,   ~role,
  "conf", "Sex, ancestry,\npre-analytics", 0.0,  3.85, "Confounder",
  "age",  "Age stratum",                   0.0,  1.55, "Exposure",
  "inf",  "Inflammation\n(CRP)",           0.0,  0.2,  "Confounder",
  "diet", "Diet + fortification\n(FFQ)",   1.95, 3.85, "Upstream",
  "fol",  "Folate\n(RBC)",                 1.95, 2.9,  "Mediator",
  "b12",  "Active B12\n(holoTC)",          1.95, 2.0,  "Mediator",
  "kid",  "Kidney\n(eGFR-cys)",            1.95, 1.1,  "Mediator",
  "mal",  "Malnutrition\n(GLIM, PhA)",     1.95, 0.2,  "Mediator",
  "gen",  "One-carbon\ngenotype",          3.55, 3.85, "Upstream",
  "hcy",  "Homocysteine",                  3.55, 1.55, "Outcome",
  "meth", "SAM/SAH,\n5-mC",                5.0,  2.45, "Downstream",
  "vasc", "ADMA,\nGDF-15",                 5.0,  0.65, "Downstream")
HW <- 0.55; HH <- 0.30   # node half-width / half-height (data units) for arrow clipping
edges <- tribble(~from, ~to, ~key,
  "age","kid",1, "age","b12",1, "age","fol",0, "age","mal",1,
  "conf","age",0, "diet","fol",0, "gen","hcy",0,
  "kid","hcy",1, "b12","hcy",1, "fol","hcy",0, "mal","hcy",1,
  "inf","mal",0, "hcy","meth",0, "hcy","vasc",0) %>%
  left_join(nodes %>% select(id, x0 = x, y0 = y), by = c("from" = "id")) %>%
  left_join(nodes %>% select(id, x1 = x, y1 = y), by = c("to" = "id")) %>%
  mutate(dx = x1 - x0, dy = y1 - y0,
         t  = pmin(HW / pmax(abs(dx), 1e-9), HH / pmax(abs(dy), 1e-9)),
         x = x0 + t * dx, y = y0 + t * dy, xend = x1 - t * dx * 1.05, yend = y1 - t * dy * 1.05)
PAL_ROLE <- c(Exposure = "#9E2A2B", Mediator = "#E09F3E", Outcome = "#2B6CB0",
              Confounder = "#9A9A9A", Upstream = "#7FB7A4", Downstream = "#6C5B7B")
p2 <- ggplot() +
  geom_segment(data = filter(edges, key == 0), aes(x, y, xend = xend, yend = yend),
               colour = "grey60", linewidth = 0.25, arrow = arrow(length = unit(2.4, "pt"), type = "closed")) +
  geom_segment(data = filter(edges, key == 1), aes(x, y, xend = xend, yend = yend),
               colour = "black", linewidth = 0.55, arrow = arrow(length = unit(2.8, "pt"), type = "closed")) +
  annotate("curve", x = 1.38, y = 3.85, xend = 1.36, yend = 2.08, curvature = 0.45, colour = "grey60",
           linewidth = 0.25, arrow = arrow(length = unit(2.4, "pt"), type = "closed")) +
  geom_label(data = nodes, aes(x, y, label = label, fill = role), family = FONT, size = 1.95,
             lineheight = 0.85, linewidth = 0, label.r = unit(2.5, "pt"),
             label.padding = unit(2.2, "pt"), colour = "black") +
  annotate("text", x = 2.5, y = -0.45, size = 1.75, family = FONT, colour = "grey30",
           label = "Bold arrows: paths tested in H1 (age \u2192 kidney / B12 / malnutrition \u2192 Hcy)") +
  scale_fill_manual(values = alpha(PAL_ROLE, 0.35), name = NULL) +
  coord_cartesian(xlim = c(-0.6, 5.6), ylim = c(-0.55, 4.2)) +
  theme_void(base_family = FONT) +
  theme(legend.position = "bottom", legend.text = element_text(size = BASE_PT - 1),
        legend.key.size = unit(6, "pt")) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(label = "")))

# ---- P3: age distribution -----------------------------------------------------
p3 <- ggplot(dat, aes(age_entry_exact, fct_rev(grp4), fill = grp4)) +
  geom_vline(xintercept = c(75.5, 84.5), linetype = "22", linewidth = 0.3, colour = COL_REF) +
  geom_density_ridges(scale = 1.6, linewidth = 0.25, colour = "white", alpha = 0.9,
                      quantile_lines = TRUE, quantiles = 2, vline_colour = "white", vline_width = 0.4,
                      rel_min_height = 0.01) +
  scale_fill_manual(values = PAL_GRP4, guide = "none") +
  scale_x_continuous("Age at entry (years)", breaks = seq(55, 105, 10), expand = c(0, 0)) +
  scale_y_discrete(NULL, expand = expansion(add = c(0.2, 1.4)))

# ---- P4: data availability matrix ---------------------------------------------
avail_vars <- c(
  "FFQ (reliable)" = "ffq_reliable", "BIA phase angle" = "phase_angle_deg", "Handgrip" = "handgrip_kg",
  "Homocysteine" = "hcy", "RBC folate / holoTC" = "holotc", "MMA" = "mma_umol_l",
  "PLP (B6)" = "plp_nmol_l", "EGRAC (B2)" = "egrac", "SAM / SAH" = "sam_nmol_l",
  "5-mC / ADMA / GDF-15" = "gdf15", "24 SNP panel" = "MTHFR_rs1801133", "Vital status" = "death_event")
avail <- map_dfr(names(avail_vars), function(nm) {
  v <- avail_vars[[nm]]
  dat %>% group_by(grp4) %>%
    summarise(p = if (v == "ffq_reliable") mean(.data[[v]] == 1) else
                   if (v == "death_event") mean(vital_status != "Lost to follow-up") else mean(!is.na(.data[[v]]))) %>%
    mutate(domain = nm)
}) %>% mutate(domain = factor(domain, rev(names(avail_vars))))
p4 <- ggplot(avail, aes(grp4, domain, fill = p)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = round(100 * p), colour = p > 0.8), size = 1.9, family = FONT) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_gradientn(colours = c("#F2F2F2", "#9DB9D8", "#2B6CB0", "#173B61"), limits = c(0.4, 1),
                       oob = squish, labels = percent, name = "Available") +
  scale_x_discrete(NULL, labels = function(x) str_replace(x, " ", "\n") %>% str_remove(" y$"), position = "top") +
  scale_y_discrete(NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x.top = element_text(size = BASE_PT - 1.5, lineheight = 0.85),
        legend.position = "right", legend.key.height = unit(14, "pt"), legend.key.width = unit(5, "pt"))

# ---- P5: fortification exposure ----------------------------------------------
p5 <- ggplot(dat, aes(age_entry_exact, pct_life_fortification, colour = grp4)) +
  geom_point(size = 0.55, alpha = 0.7, stroke = 0) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, colour = "black", linewidth = 0.4, span = 0.6) +
  scale_colour_manual(values = PAL_GRP4, name = NULL) +
  scale_x_continuous("Age at entry (years)", breaks = seq(55, 105, 10)) +
  scale_y_continuous("Life lived under mandatory\nfolic-acid fortification (%)", limits = c(0, NA)) +
  theme(legend.position = "none")

# ---- Assemble -------------------------------------------------------------------
fig1 <-wrap_elements(full = (p1 | p2) + plot_layout(widths = c(1.15, 1))) /
        wrap_elements(full = (p3 | p4 | p5) + plot_layout(widths = c(0.85, 1.25, 0.9))) +
        plot_layout(heights = c(1.3, 1))

save_fig(fig1, "Figure1_study_design", width = 180, height = 150)
save_panels(list(flow = p1, dag = p2, age = p3, availability = p4, fortification = p5), "Figure1")
