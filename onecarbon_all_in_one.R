# =============================================================================
#  ONE-CARBON / HOMOCYSTEINE — ALL-IN-ONE PIPELINE (Clinical Nutrition)
#  Figures 1–7, Supplementary Figures S1–S5, Tables 1–3 and S1–S4
#
#  HOW TO USE
#    1. Open this file in RStudio and click "Source" (or run: source("onecarbon_all_in_one.R")).
#    2. A window opens: select the Excel workbook (sheets Data, Data_Dictionary, Genotype_QC).
#    3. Everything is generated in a folder "onecarbon_output" next to the workbook.
#  The workbook is read ONCE; every figure and table uses the same data in memory.
#
#  Options (edit below if you want):
#    N_BOOT       bootstrap replicates (1000 = final; 200 = quick draft)
#    SAVE_PANELS  TRUE = also export each panel as a separate TIFF
# =============================================================================

N_BOOT_USER      <- 1000
SAVE_PANELS_USER <- FALSE

# ---- 0. Packages (installs anything missing) -----------------------------------
.pk <- c("readxl", "dplyr", "tidyr", "purrr", "stringr", "forcats", "tibble", "ggplot2", "patchwork", "scales",
         "ragg", "systemfonts", "colorspace", "ggridges", "ggdist", "ggalluvial", "ggrepel", "ggtext", "FactoMineR",
         "mclust", "mgcv", "relaimpo", "quantreg", "MASS", "qgraph", "igraph", "ggraph", "tidygraph",
         "lavaan", "partykit", "pROC", "survival", "flextable", "officer", "car", "AER", "mice", "janitor", "broom",
         "Hmisc")
.new <- setdiff(.pk, rownames(installed.packages()))
if (length(.new)) install.packages(.new, repos = "https://cloud.r-project.org")

# ---- 1. Choose the workbook ------------------------------------------------------
pick_workbook <- function() {
  cap <- "Select the one-carbon workbook (.xlsx)"
  # Non-interactive use (Rscript): pass the path in the ONECARBON_DATA environment variable
  if (!interactive() && file.exists(Sys.getenv("ONECARBON_DATA"))) return(normalizePath(Sys.getenv("ONECARBON_DATA")))
  f <- character(0)
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    f <- rstudioapi::selectFile(caption = cap, filter = "Excel files (*.xlsx)")
  } else if (.Platform$OS.type == "windows") {
    f <- utils::choose.files(caption = cap, multi = FALSE, filters = matrix(c("Excel (*.xlsx)", "*.xlsx"), 1))
  } else if (capabilities("tcltk")) {
    f <- tcltk::tk_choose.files(caption = cap, multi = FALSE, filters = matrix(c("Excel", ".xlsx"), 1))
  } else {
    message(cap); f <- file.choose()
  }
  if (!length(f) || is.null(f) || is.na(f) || !nzchar(f)) stop("No file selected.")
  normalizePath(f)
}
DATA_FILE <- pick_workbook()
stopifnot(all(c("Data", "Data_Dictionary", "Genotype_QC") %in% readxl::excel_sheets(DATA_FILE)))
OUT_ROOT <- file.path(dirname(DATA_FILE), "onecarbon_output")
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)
.old_wd <- setwd(OUT_ROOT)
Sys.setenv(ONECARBON_DATA = DATA_FILE, ONECARBON_NBOOT = N_BOOT_USER, ONECARBON_PANELS = SAVE_PANELS_USER)
message("Workbook: ", DATA_FILE, "\nOutputs:  ", OUT_ROOT)

# Runs one block in its own environment (objects of one figure do not leak into the next);
# an error in one block is reported and the pipeline continues.
.status <- list()
run_block <- function(name, expr) {
  message("\n==== ", name, " ====")
  t0 <- Sys.time()
  env <- new.env(parent = globalenv())
  assign("map", purrr::map, env); assign("select", dplyr::select, env); assign("recode", dplyr::recode, env)
  ok <- tryCatch({ eval(expr, envir = env); TRUE },
                 error = function(e) { message("  ERROR in ", name, ": ", conditionMessage(e)); FALSE })
  .status[[name]] <<- ok
  message(sprintf("  %s in %.0f s", if (ok) "done" else "FAILED", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# =============================================================================
# 2. SETUP (data read once)
# =============================================================================
# =============================================================================
# One-carbon metabolism, kidney, B12 and malnutrition across age strata
# 00_setup.R — packages, data preparation, palette, theme and export helpers
# Sourced by every figure/table script. Run from the project root.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(forcats)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ragg)
})

options(dplyr.summarise.inform = FALSE, warn = 1)
set.seed(20261010)

# ---- Paths ------------------------------------------------------------------
DATA_PATH <- Sys.getenv("ONECARBON_DATA", "banco_de_dados_onecarbon.xlsx")
OUT_FIG   <- "figures"
OUT_SUPP  <- "supplementary"
OUT_TAB   <- "tables"
for (d in c(OUT_FIG, OUT_SUPP, OUT_TAB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Number of bootstrap replicates (use 1000–2000 for the final manuscript)
N_BOOT <- as.integer(Sys.getenv("ONECARBON_NBOOT", "1000"))

# Export individual panels as well as the assembled figure (FALSE = composite only)
SAVE_PANELS <- as.logical(Sys.getenv("ONECARBON_PANELS", "FALSE"))

# ---- Font (Arial if installed; Liberation Sans is metric-identical) ---------
.fams <- unique(systemfonts::system_fonts()$family)
FONT  <- if ("Arial" %in% .fams) "Arial" else if ("Liberation Sans" %in% .fams) "Liberation Sans" else "sans"

# ---- Palette: one fixed colour per stratum in every figure ------------------
# Colour-blind safe (distinct in hue AND lightness)
GRP4_LEVELS <- c("Community 56–75 y", "Hospital 56–75 y",
                 "Hospital 76–84 y", "Hospital ≥85 y")
STRATA      <- c("56–75 y", "76–84 y", "≥85 y")

GRP4_SHORT <- setNames(c("Comm. 56\u201375", "Hosp. 56\u201375", "Hosp. 76\u201384", "Hosp. \u226585"), GRP4_LEVELS)
lab_grp4 <- function(x) unname(GRP4_SHORT[x])
rot_x <- theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
PAL_GRP4 <- setNames(c("#9A9A9A", "#2B6CB0", "#E09F3E", "#9E2A2B"), GRP4_LEVELS)
PAL_STR  <- setNames(c("#2B6CB0", "#E09F3E", "#9E2A2B"), STRATA)
PAL_GLIM <- c("No malnutrition" = "#7FB7A4", "Moderate" = "#E8B04B", "Severe" = "#B23A48")
PAL_DIV  <- c(low = "#2B6CB0", mid = "#F7F7F7", high = "#9E2A2B")     # diverging (beta/rho)
PAL_SEQ  <- c("#FFF7EC", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000")  # sequential (Hcy)
COL_REF  <- "grey35"                                                    # reference/cut-off lines

# ---- Theme ------------------------------------------------------------------
BASE_PT <- 7
theme_paper <- function(base = BASE_PT, grid = FALSE) {
  th <- theme_classic(base_size = base, base_family = FONT) +
    theme(
      text              = element_text(colour = "black"),
      axis.text         = element_text(size = base - 0.5, colour = "black"),
      axis.title        = element_text(size = base),
      axis.line         = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(1.2, "pt"),
      plot.title        = element_text(size = base + 0.5, face = "bold", hjust = 0,
                                       margin = margin(b = 3)),
      plot.subtitle     = element_text(size = base - 0.5, colour = "grey25", hjust = 0),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base, face = "bold", margin = margin(2, 2, 2, 2)),
      legend.title      = element_text(size = base - 0.5, face = "bold"),
      legend.text       = element_text(size = base - 0.5),
      legend.key.size   = unit(7, "pt"),
      legend.background = element_blank(),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt"),
      plot.margin       = margin(4, 5, 4, 4),
      panel.spacing     = unit(5, "pt")
    )
  if (grid) th <- th + theme(panel.grid.major = element_line(linewidth = 0.15, colour = "grey88"))
  th
}
theme_set(theme_paper())
update_geom_defaults("point", list(size = 0.8))
update_geom_defaults("text",  list(family = FONT, size = BASE_PT / .pt * 0.9))

# ---- Export helpers ---------------------------------------------------------
# TIFF, 600 dpi, LZW (Elsevier artwork standard). 180 mm = double-column width.
save_fig <- function(plot, name, width = 180, height = 200, dir = OUT_FIG, dpi = 600) {
  f <- file.path(dir, paste0(name, ".tiff"))
  ragg::agg_tiff(f, width = width, height = height, units = "mm", res = dpi,
                 compression = "lzw", background = "white")
  print(plot)
  invisible(dev.off())
  message("saved: ", f)
  invisible(f)
}
save_panels <- function(panels, prefix, width = 85, height = 70, dir = OUT_FIG) {
  if (!SAVE_PANELS) return(invisible(NULL))
  pd <- file.path(dir, paste0(prefix, "_panels")); dir.create(pd, showWarnings = FALSE)
  iwalk(panels, ~ save_fig(.x, paste0(prefix, "_", .y), width, height, dir = pd))
}

# Row of panels, each wrapped so axis titles do not get pushed by neighbours' tick labels
row_of <- function(..., widths = NULL) wrap_plots(map(list(...), ~ wrap_elements(full = .x)), nrow = 1, widths = widths)

# ---- Formatting helpers -----------------------------------------------------
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, 3, format = "f"))))
fmt_p_short <- function(p) ifelse(p < 0.001, "<0.001", formatC(p, 3, format = "f"))
stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
pct_change <- function(b) 100 * (exp(b) - 1)       # log-scale beta -> % change
lab_umol <- "Homocysteine (µmol/L)"

spearman_lab <- function(x, y) {
  ok <- complete.cases(x, y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  sprintf("ρ = %.2f, %s", ct$estimate, fmt_p(ct$p.value))
}

# Partial Spearman: rank-transform, residualise on covariates, Pearson of residuals
rank_resid <- function(v, Z) {
  r <- rank(v, na.last = "keep")
  if (is.null(Z)) return(r)
  ok <- complete.cases(r, Z)
  out <- rep(NA_real_, length(r))
  out[ok] <- resid(lm(r[ok] ~ ., data = as.data.frame(Z[ok, , drop = FALSE])))
  out
}

# ---- Data -------------------------------------------------------------------
raw <- read_excel(DATA_PATH, sheet = "Data", guess_max = 2000)
dict <- read_excel(DATA_PATH, sheet = "Data_Dictionary")
geno_qc <- read_excel(DATA_PATH, sheet = "Genotype_QC")

dat <- raw %>%
  mutate(
    grp4 = case_when(
      group == "Community control (ELSA)" ~ GRP4_LEVELS[1],
      group == "Hospital control"          ~ GRP4_LEVELS[2],
      str_detect(group, "76")              ~ GRP4_LEVELS[3],
      str_detect(group, "Longevo")         ~ GRP4_LEVELS[4]) %>% factor(GRP4_LEVELS),
    stratum = case_when(
      str_detect(age_group_analytic, "Younger") ~ STRATA[1],
      str_detect(age_group_analytic, "^Old \\(") ~ STRATA[2],
      str_detect(age_group_analytic, "Oldest")   ~ STRATA[3]) %>% factor(STRATA),
    hospital  = setting == "Hospital",
    sex       = factor(sex, c("F", "M"), c("Female", "Male")),
    hcy       = homocysteine_umol_l,
    log_hcy   = log(hcy),
    hhcy      = as.integer(hcy > 15),
    egfr      = egfr_cys_ckdepi2012,
    log_egfr  = log(egfr),
    holotc    = holotc_pmol_l,
    log_holotc = log(holotc),
    rbc_fol   = folate_rbc_ng_ml,
    log_rbcf  = log(rbc_fol),
    log_plp   = log(plp_nmol_l),
    log_mma   = log(mma_umol_l),
    log_crp   = log(crp_mg_l),
    log_umfa  = log(umfa_nmol_l),
    egrac     = egrac_riboflavin,
    glim3     = factor(case_when(
      str_detect(glim_malnutrition, "^No")       ~ "No malnutrition",
      str_detect(glim_malnutrition, "Moderate")  ~ "Moderate",
      str_detect(glim_malnutrition, "Severe")    ~ "Severe"), names(PAL_GLIM)),
    glim_any  = as.integer(glim3 != "No malnutrition"),
    sarc4     = factor(case_when(
      str_detect(sarcopenia_ewgsop2, "^No")       ~ "No sarcopenia",
      str_detect(sarcopenia_ewgsop2, "Probable")  ~ "Probable",
      str_detect(sarcopenia_ewgsop2, "Confirmed") ~ "Confirmed",
      str_detect(sarcopenia_ewgsop2, "Severe")    ~ "Severe"),
      c("No sarcopenia", "Probable", "Confirmed", "Severe")),
    mna3      = factor(mna_category, c("Normal", "At risk", "Malnourished")),
    ckd_cys   = factor(ckd_stage_cys, c("G1", "G2", "G3a", "G3b", "G4", "G5")),
    ckd_cys4  = fct_collapse(ckd_cys, "G1–G2" = c("G1", "G2"), "G4–G5" = c("G4", "G5")),
    mthfr677  = factor(MTHFR_rs1801133_dosage_A, 0:2, c("CC", "CT", "TT")),
    dhfr      = factor(DHFR_rs70991108_dosage_del, 0:2, c("ins/ins", "ins/del", "del/del")),
    grs_tert  = factor(hcy_genetic_risk_tertile, c("Low", "Intermediate", "High")),
    holotc_low = as.integer(holotc < 35),
    sam_sah   = sam_sah_ratio,
    gdf15     = gdf15_pg_ml,
    adma      = adma_umol_l,
    mc5       = global_dna_methylation_5mc_pct
  )

# Hospital analytic set (primary H1 analyses): 3 age strata, homocysteine available
hosp <- dat %>% filter(hospital, !is.na(hcy)) %>% mutate(stratum = droplevels(stratum))
# Gene x nutrient analytic set (all four groups)
gxn  <- dat %>% filter(analysis_set_gxn == 1)

# Core adjustment set used across models (pre-analytics + design)
COV_BASE <- c("sex", "ancestry_pc1", "ancestry_pc2", "time_to_centrifugation_min", "fasting_hours")

# 24 one-carbon variants (core + GWAS GRS components + nutrient-axis + negative control)
snp_tab <- dict %>%
  filter(str_detect(Variable, "_dosage_"),
         str_detect(`Role in one-carbon paper`, "Effect modifier|Negative")) %>%
  transmute(var = Variable,
            gene = str_extract(var, "^[A-Z0-9]+"),
            rsid = str_extract(var, "rs[0-9]+"),
            allele = str_remove(var, ".*_dosage_"),
            role = case_when(
              str_detect(`Role in one-carbon paper`, "core")     ~ "Core one-carbon",
              str_detect(`Role in one-carbon paper`, "GWAS")     ~ "Hcy GWAS (GRS)",
              str_detect(`Role in one-carbon paper`, "nutrient") ~ "Nutrient axis",
              str_detect(`Role in one-carbon paper`, "Negative") ~ "Negative control")) %>%
  mutate(label = paste0(gene, " ", rsid))

message(sprintf("Data loaded: %d participants; hospital set n = %d; GxN set n = %d; font = %s",
                nrow(dat), nrow(hosp), nrow(gxn), FONT))

# Guard against masking by packages attached in other scripts (mclust::map, MASS::select, car::recode)
map <- purrr::map; select <- dplyr::select; recode <- dplyr::recode

# =============================================================================
# 3. SHARED MODELS
# =============================================================================
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

# =============================================================================
# 4. FIGURE 1  (fig1_study_design.R)
# =============================================================================
run_block("Figure 1", quote({
# =============================================================================
# Figure 1 — Study design, analytic sets and causal framework
#   P1  Participant flow: 4 groups -> analytic sets (hospital H1 set, GxN set)
#   P2  Directed acyclic graph (DAG) encoding H1 and the adjustment set
#   P3  Age distribution by group, with analytic strata cut-points
#   P4  Data availability matrix (assessment domain x group)
#   P5  Lifetime exposure to mandatory folic-acid fortification by group
# =============================================================================
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
}))

# =============================================================================
# 5. FIGURE 2  (fig2_nutritional_phenotype.R)
# =============================================================================
run_block("Figure 2", quote({
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
}))

# =============================================================================
# 6. FIGURE 3  (fig3_diet_to_status.R)
# =============================================================================
run_block("Figure 3", quote({
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
}))

# =============================================================================
# 7. FIGURE 4  (fig4_determinants_hcy.R)
# =============================================================================
run_block("Figure 4", quote({
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
}))

# =============================================================================
# 8. FIGURE 5  (fig5_counterfactual.R)
# =============================================================================
run_block("Figure 5", quote({
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
}))

# =============================================================================
# 9. FIGURE 6  (fig6_genetics.R)
# =============================================================================
run_block("Figure 6", quote({
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
}))

# =============================================================================
# 10. FIGURE 7  (fig7_integration.R)
# =============================================================================
run_block("Figure 7", quote({
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
}))

# =============================================================================
# 11. TABLES 1-3  (tables_main.R)
# =============================================================================
run_block("Tables 1-3", quote({
# =============================================================================
# Tables 1–3 (Word, via flextable/officer)
#   Table 1  Characteristics by group (biomarkers excluded: they are shown in the figures)
#   Table 2  Dietary intake, adequacy vs EAR/UL, diet quality and supplements
#   Table 3  Hierarchical models of log-homocysteine (M1–M5): % change per SD (95% CI), R2, VIF
# =============================================================================
suppressPackageStartupMessages({ library(flextable); library(officer); library(car) })

# ---- Generic descriptive helpers ------------------------------------------------
med_iqr <- function(x, d = 1) {
  x <- x[!is.na(x)]
  sprintf(paste0("%.", d, "f [%.", d, "f–%.", d, "f]"), median(x), quantile(x, 0.25), quantile(x, 0.75))
}
n_pct <- function(x) { x <- x[!is.na(x)]; sprintf("%d (%.1f)", sum(x), 100 * mean(x)) }

# p across hospital strata (Kruskal / chi-square) and p-trend (Spearman / Cochran–Armitage via logistic)
p_cont <- function(v) {
  hh <- dat %>% filter(hospital, !is.na(.data[[v]]))
  c(kruskal.test(hh[[v]] ~ hh$stratum)$p.value,
    cor.test(as.numeric(hh$stratum), hh[[v]], method = "spearman", exact = FALSE)$p.value)
}
p_bin <- function(x) {
  hh <- dat %>% mutate(.x = x) %>% filter(hospital, !is.na(.x))
  tb <- table(hh$stratum, hh$.x)
  pc <- if (any(suppressWarnings(chisq.test(tb))$expected < 5)) fisher.test(tb, simulate.p.value = TRUE, B = 1e4)$p.value else chisq.test(tb)$p.value
  pt <- suppressWarnings(coef(summary(glm(.x ~ as.numeric(stratum), data = hh, family = binomial)))[2, 4])
  c(pc, pt)
}
row_cont <- function(label, v, d = 1, section = NA) {
  ps <- p_cont(v)
  tibble(section = section, Characteristic = label,
         !!!setNames(map(GRP4_LEVELS, ~ med_iqr(dat[[v]][dat$grp4 == .x], d)), GRP4_LEVELS),
         P = fmt_p_short(ps[1]), `P trend` = fmt_p_short(ps[2]))
}
row_bin <- function(label, x, section = NA) {
  ps <- p_bin(x)
  tibble(section = section, Characteristic = label,
         !!!setNames(map(GRP4_LEVELS, ~ n_pct(x[dat$grp4 == .x])), GRP4_LEVELS),
         P = fmt_p_short(ps[1]), `P trend` = fmt_p_short(ps[2]))
}
row_head <- function(label) tibble(section = label, Characteristic = label)

make_ft <- function(tab, title, note) {
  tab <- tab %>% mutate(across(everything(), ~ replace_na(as.character(.x), "")))
  is_head <- tab$Characteristic == tab$section & tab$Characteristic != ""
  tab2 <- tab %>% select(-section)
  ns <- dat %>% count(grp4)
  hdr <- c("Characteristic", sprintf("%s\n(n = %d)", GRP4_LEVELS, ns$n), "P†", "P trend†")
  ft <- flextable(tab2) %>%
    set_header_labels(values = setNames(as.list(hdr), names(tab2))) %>%
    padding(padding.top = 1, padding.bottom = 1, part = "all") %>% line_spacing(space = 1, part = "all") %>%
    bold(i = which(is_head), j = 1) %>%
    bg(i = which(is_head), bg = "#F2F2F2") %>%
    padding(i = which(!is_head), j = 1, padding.left = 10) %>%
    font(fontname = "Arial", part = "all") %>% fontsize(size = 8, part = "all") %>% fontsize(size = 8, part = "header") %>%
    bold(part = "header") %>% align(j = 2:ncol(tab2), align = "center", part = "all") %>%
    valign(valign = "center", part = "header") %>%
    border_remove() %>% hline_top(border = fp_border(width = 1), part = "all") %>%
    hline_bottom(border = fp_border(width = 1), part = "all") %>%
    hline(i = 1, border = fp_border(width = 0.5), part = "header") %>%
    set_caption(caption = title) %>%
    add_footer_lines(note) %>% fontsize(size = 7, part = "footer") %>%
    autofit() %>% width(j = 1, width = 2.3)
  ft
}

# ---- Table 1 --------------------------------------------------------------------------
d <- dat
T1 <- bind_rows(
  row_head("Sociodemographic"),
  row_cont("Age, years", "age_entry_exact"),
  row_bin("Female sex", d$sex == "Female"),
  row_bin("Self-reported White", d$race_color_ibge == "White"),
  row_bin("Self-reported Brown (Parda)", str_detect(d$race_color_ibge, "Brown")),
  row_bin("Self-reported Black", d$race_color_ibge == "Black"),
  row_cont("Schooling, years", "education_years"),
  row_cont("Household income, minimum wages", "household_income_mw"),
  row_bin("Widowed", d$marital_status == "Widowed"),
  row_bin("Lives alone", d$lives_alone == 1),
  row_bin("Long-term care facility", d$institutionalized == 1),
  row_head("Lifestyle"),
  row_bin("Current smoker", d$smoking_status == "Current"),
  row_bin("Heavy alcohol use (≥30 g/day)", str_detect(d$alcohol_use, "Heavy")),
  row_bin("Low physical activity (IPAQ <600)", str_detect(d$ipaq_category, "Low")),
  row_head("Comorbidities and medications"),
  row_bin("Hypertension", d$hypertension == 1),
  row_bin("Type 2 diabetes", d$diabetes_t2 == 1),
  row_bin("Chronic kidney disease (G3+)", d$chronic_kidney_disease == 1),
  row_bin("Heart failure", d$heart_failure == 1),
  row_bin("Coronary artery disease", d$coronary_artery_disease == 1),
  row_bin("Previous stroke/TIA", d$stroke_history == 1),
  row_bin("Dementia", d$dementia == 1),
  row_bin("Cancer history", d$cancer_history == 1),
  row_cont("Age-adjusted Charlson index", "charlson_age_adj", 0),
  row_bin("Polypharmacy (≥5 drugs)", d$polypharmacy == 1),
  row_bin("Proton-pump inhibitor", d$med_ppi == 1),
  row_bin("Metformin", d$med_metformin == 1),
  row_bin("Diuretic", d$med_diuretic == 1),
  row_head("Geriatric assessment"),
  row_cont("Katz ADL (0–6)", "katz_adl", 0),
  row_cont("Lawton IADL (0–8)", "lawton_iadl", 0),
  row_cont("MMSE (0–30)", "mmse", 0),
  row_cont("Clinical Frailty Scale (1–9)", "clinical_frailty_scale", 0),
  row_bin("SARC-F ≥4", d$sarc_f >= 4),
  row_head("Nutritional status"),
  row_cont("BMI, kg/m²", "bmi"),
  row_cont("Weight loss in 6 months, %", "weight_loss_6m_pct"),
  row_cont("Calf circumference, cm", "calf_circ_cm"),
  row_cont("MNA-SF score (0–14)", "mna_sf", 0),
  row_bin("MNA-SF at risk (8–11)", d$mna3 == "At risk"),
  row_bin("MNA-SF malnourished (0–7)", d$mna3 == "Malnourished"),
  row_bin("NRS-2002 ≥3 (hospital only)", d$nrs2002_at_risk == 1),
  row_bin("GLIM moderate malnutrition", d$glim3 == "Moderate"),
  row_bin("GLIM severe malnutrition", d$glim3 == "Severe"),
  row_bin("Sarcopenia, confirmed or severe (EWGSOP2)", d$sarc4 %in% c("Confirmed", "Severe")),
  row_head("Body composition and function"),
  row_cont("Phase angle, °", "phase_angle_deg"),
  row_cont("FFMI, kg/m²", "ffmi_kg_m2"),
  row_cont("ASMI, kg/m²", "asmi_kg_m2", 2),
  row_cont("Handgrip, kg", "handgrip_kg"),
  row_cont("Gait speed, m/s", "gait_speed_m_s", 2),
  row_head("Index admission (hospital only)"),
  row_bin("Infection at admission", d$infection_at_admission == 1),
  row_bin("Acute kidney injury at admission", d$acute_kidney_injury == 1),
  row_cont("Length of stay, days", "length_of_stay_days", 0)
)
note1 <- paste("Values are median [IQR] or n (%). Biomarkers are presented in Figs. 3–4.",
               "†P across the three hospital strata (Kruskal–Wallis or χ²/Fisher) and P for linear trend",
               "(Spearman or logistic regression on stratum rank); the community group is shown for reference only.",
               "ADL, activities of daily living; ASMI, appendicular skeletal muscle mass index; FFMI, fat-free mass index;",
               "GLIM, Global Leadership Initiative on Malnutrition; IADL, instrumental ADL; MMSE, Mini-Mental State Examination;",
               "MNA-SF, Mini Nutritional Assessment–Short Form; NRS, Nutritional Risk Screening.")
ft1 <- make_ft(T1, "Table 1. Characteristics of the participants by age group and setting.", note1)

# ---- Table 2 --------------------------------------------------------------------------
d <- dat
T2 <- bind_rows(
  row_head("Energy and protein"),
  row_cont("Energy, kcal/day", "energy_kcal", 0),
  row_cont("Energy, kcal/kg/day", "energy_kcal_per_kg"),
  row_cont("Protein, g/kg/day", "protein_g_per_kg", 2),
  row_bin("Energy ≥30 kcal/kg/day (ESPEN)", d$energy_adequate_espen == 1),
  row_bin("Protein ≥1.0 g/kg/day (ESPEN)", d$protein_adequate_espen == 1),
  row_head("One-carbon nutrients (energy-adjusted)"),
  row_cont("Total folate, µg DFE/day", "total_folate_dfe_ug_eadj", 0),
  row_cont("Natural food folate, µg/day", "food_folate_ug_eadj", 0),
  row_cont("Folic acid from fortified flour, µg/day", "folic_acid_fortified_ug", 0),
  row_cont("Total folic acid (fortified + supplements), µg/day", "folic_acid_total_ug", 0),
  row_cont("Vitamin B12, µg/day", "vitamin_b12_ug_eadj"),
  row_cont("Vitamin B6, mg/day", "vitamin_b6_mg_eadj", 2),
  row_cont("Riboflavin, mg/day", "riboflavin_mg_eadj", 2),
  row_cont("Choline, mg/day", "choline_mg_eadj", 0),
  row_cont("Betaine, mg/day", "betaine_mg_eadj", 0),
  row_cont("Methionine, g/day", "methionine_g_eadj", 2),
  row_head("Adequacy"),
  row_bin("Total folate < EAR (320 µg DFE)", d$folate_dfe_below_ear == 1),
  row_bin("Folic acid > UL (1000 µg)", d$folic_acid_above_ul == 1),
  row_bin("Vitamin B12 < EAR (2.0 µg)", d$b12_intake_below_ear == 1),
  row_bin("Vitamin B6 < EAR", d$b6_intake_below_ear == 1),
  row_bin("Riboflavin < EAR", d$riboflavin_intake_below_ear == 1),
  row_head("Diet quality and pattern"),
  row_cont("BHEI-R (0–100)", "bhei_r_score"),
  row_cont("MEDAS (0–14)", "medas_score", 0),
  row_cont("Dietary Inflammatory Index", "dii_score", 2),
  row_bin("Prudent pattern", d$dietary_pattern_apriori == "Prudent"),
  row_bin("Western/ultra-processed pattern", str_detect(d$dietary_pattern_apriori, "Western")),
  row_bin("Traditional (rice & beans) pattern", str_detect(d$dietary_pattern_apriori, "Traditional")),
  row_cont("Ultra-processed foods, % energy", "upf_pct_energy"),
  row_head("Supplements"),
  row_bin("Folic acid", d$supp_folic_acid == 1),
  row_bin("Vitamin B12 (oral or injection)", d$supp_vitamin_b12 == 1),
  row_bin("Multivitamin", d$supp_multivitamin == 1),
  row_bin("Oral nutritional supplement before admission", d$supp_oral_nutritional == 1),
  row_head("FFQ quality"),
  row_bin("Proxy informant", d$diet_informant == "Proxy"),
  row_bin("Major diet change in previous year", d$diet_change_last_year == 1),
  row_bin("FFQ considered reliable", d$ffq_reliable == 1))
note2 <- paste("Values are median [IQR] or n (%). Nutrients are energy-adjusted by the residual method (log scale).",
               "EAR/UL: Institute of Medicine values for adults ≥51 y. †As in Table 1.",
               "BHEI-R, Brazilian Healthy Eating Index–Revised; DFE, dietary folate equivalents; EAR, estimated average requirement;",
               "MEDAS, Mediterranean Diet Adherence Screener; UL, tolerable upper intake level.")
ft2 <- make_ft(T2, "Table 2. Dietary intake of one-carbon nutrients, adequacy and diet quality by age group.", note2)

# ---- Table 3: hierarchical models ------------------------------------------------------
cont <- c("ancestry_pc1", "ancestry_pc2", "time_to_centrifugation_min", "fasting_hours", "log_egfr", "log_holotc",
          "log_rbcf", "log_plp", "egrac", "phase_angle_deg", "log_crp", "hcy_genetic_risk_score")
mz <- mdat %>% mutate(across(all_of(cont), ~ as.numeric(scale(.x))))
hz <- map(seq_along(BLOCKS), ~ lm(block_formula(.x), data = mz))
LAB3 <- c(s76 = "Age 76–84 y (vs 56–75 y)", s85 = "Age ≥85 y (vs 56–75 y)", male = "Male sex",
          ancestry_pc1 = "Ancestry PC1", ancestry_pc2 = "Ancestry PC2",
          time_to_centrifugation_min = "Time to centrifugation", fasting_hours = "Fasting hours",
          log_egfr = "eGFR-cys (log)", log_holotc = "holoTC (log)", log_rbcf = "RBC folate (log)", log_plp = "PLP (log)",
          egrac = "EGRAC", glim_any = "GLIM malnutrition (any)", phase_angle_deg = "Phase angle",
          log_crp = "CRP (log)", diabetes_t2 = "Type 2 diabetes", hcy_genetic_risk_score = "Hcy genetic risk score")
cell <- function(m, v) {
  if (!v %in% names(coef(m))) return("")
  ci <- confint(m)[v, ]; p <- coef(summary(m))[v, 4]
  sprintf("%.1f (%.1f, %.1f)%s", pct_change(coef(m)[v]), pct_change(ci[1]), pct_change(ci[2]), stars(p))
}
T3 <- tibble(Predictor = unname(LAB3), !!!setNames(map(hz, function(m) map_chr(names(LAB3), ~ cell(m, .x))), paste0("M", 1:5)))
fitrow <- function(lab, f) tibble(Predictor = lab, !!!setNames(map(hz, f), paste0("M", 1:5)))
r2s <- map_dbl(hz, ~ summary(.x)$r.squared)
T3 <- bind_rows(T3,
  fitrow("R²", function(m) sprintf("%.3f", summary(m)$r.squared)),
  fitrow("Adjusted R²", function(m) sprintf("%.3f", summary(m)$adj.r.squared)),
  tibble(Predictor = "ΔR² (vs previous model)", !!!setNames(as.list(c("—", sprintf("%.3f", diff(r2s)))), paste0("M", 1:5))),
  fitrow("P for ΔR²", function(m) "") ,
  fitrow("Maximum VIF", function(m) sprintf("%.2f", max(car::vif(m)))),
  fitrow("n", function(m) as.character(nobs(m))))
pdr <- c("—", map_chr(2:5, ~ fmt_p_short(anova(hz[[.x - 1]], hz[[.x]])$`Pr(>F)`[2])))
T3[T3$Predictor == "P for ΔR²", -1] <- as.list(pdr)
note3 <- paste("Linear regression of log-homocysteine in the hospital analytic set (complete cases).",
               "Estimates are % change in homocysteine (95% CI) per SD of continuous predictors, or vs the reference category,",
               "computed as 100×(e^β − 1). M1: age stratum, sex, ancestry and pre-analytical variables; M2: + kidney function",
               "(eGFR from cystatin C); M3: + active B12 (holotranscobalamin); M4: + RBC folate, PLP and EGRAC;",
               "M5: + GLIM malnutrition, phase angle, CRP, type 2 diabetes and the homocysteine genetic risk score.",
               "*P < 0.05; **P < 0.01; ***P < 0.001. VIF, variance inflation factor (generalized).")
ft3 <- flextable(T3) %>%
  padding(padding.top = 1.5, padding.bottom = 1.5, part = "all") %>% line_spacing(space = 1, part = "all") %>%
  font(fontname = "Arial", part = "all") %>% fontsize(size = 8, part = "all") %>% bold(part = "header") %>%
  align(j = 2:6, align = "center", part = "all") %>%
  bold(i = 1:2, j = 1) %>% bg(i = 1:2, bg = "#FBF3E4") %>%
  border_remove() %>% hline_top(border = fp_border(width = 1), part = "all") %>%
  hline_bottom(border = fp_border(width = 1), part = "all") %>%
  hline(i = length(LAB3), border = fp_border(width = 0.5)) %>%
  set_caption(caption = "Table 3. Hierarchical models of homocysteine: attenuation of the age gradient by kidney function, active B12 and nutritional status.") %>%
  add_footer_lines(note3) %>% fontsize(size = 7, part = "footer") %>% autofit() %>% width(j = 1, width = 2.2)

# ---- Write -------------------------------------------------------------------------------
doc <- read_docx() %>%
  body_add_flextable(fit_to_width(ft1, max_width = 7)) %>% body_add_break() %>%
  body_add_flextable(fit_to_width(ft2, max_width = 7)) %>%
  body_end_section_portrait() %>%
  body_add_flextable(fit_to_width(ft3, max_width = 9.5)) %>%
  body_end_section_landscape()
print(doc, target = file.path(OUT_TAB, "Tables_1-3.docx"))
write.csv(T1, file.path(OUT_TAB, "Table1.csv"), row.names = FALSE)
write.csv(T2, file.path(OUT_TAB, "Table2.csv"), row.names = FALSE)
write.csv(T3, file.path(OUT_TAB, "Table3.csv"), row.names = FALSE)
message("saved: ", file.path(OUT_TAB, "Tables_1-3.docx"))
}))

# =============================================================================
# 12. SUPPLEMENTARY  (supplementary.R)
# =============================================================================
run_block("Supplementary", quote({
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
}))

# =============================================================================
# DONE
# =============================================================================
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
setwd(.old_wd)
message("\n==== Summary ====")
for (n in names(.status)) message(sprintf("  %-14s %s", n, if (.status[[n]]) "OK" else "FAILED"))
message("Outputs in: ", OUT_ROOT)
