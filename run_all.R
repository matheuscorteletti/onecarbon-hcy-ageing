# =============================================================================
# Run the full pipeline (figures + tables). From the project root:
#   Rscript run_all.R
# Options (environment variables):
#   ONECARBON_DATA   path to the .xlsx   (default: banco_de_dados_onecarbon.xlsx)
#   ONECARBON_NBOOT  bootstrap replicates (default: 1000; use 200 for a quick draft)
#   ONECARBON_PANELS TRUE to also export every panel as its own TIFF
# =============================================================================
scripts <- c("R/fig1_study_design.R", "R/fig2_nutritional_phenotype.R", "R/fig3_diet_to_status.R",
             "R/fig4_determinants_hcy.R", "R/fig5_counterfactual.R", "R/fig6_genetics.R",
             "R/fig7_integration.R", "R/tables_main.R", "R/supplementary.R")
for (s in scripts) {
  message("\n==== ", s, " ====")
  t0 <- Sys.time()
  local(source(s, local = TRUE))
  message(sprintf("done in %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
