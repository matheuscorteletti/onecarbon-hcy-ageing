# =============================================================================
# Tables 1–3 (Word, via flextable/officer)
#   Table 1  Characteristics by group (biomarkers excluded: they are shown in the figures)
#   Table 2  Dietary intake, adequacy vs EAR/UL, diet quality and supplements
#   Table 3  Hierarchical models of log-homocysteine (M1–M5): % change per SD (95% CI), R2, VIF
# =============================================================================
source("R/00_setup.R"); source("R/01_models.R")
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
