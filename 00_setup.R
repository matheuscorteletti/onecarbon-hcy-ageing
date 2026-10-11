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
OUT_FIG   <- "output/figures"
OUT_SUPP  <- "output/supplementary"
OUT_TAB   <- "output/tables"
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
