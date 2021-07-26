# Run MR, GIP and GIPR combined

library(data.table)
library(MendelianRandomization)
library(mrpipe)
library(ieugwasr)

d_colocres <- fread("results/all_results_coloc_T2DMMahajan.txt", check.names = T)

outcomes <- d_colocres[PP.H4.abf_GIP > 0.8 & PP.H4.abf_GIPR > 0.8,][["outcome"]]
exposures <- c("GIP", "GIPR")

res <- lapply(outcomes, function(outcome){
 # outcome <- outcomes[1]
 d_combined <- do.call("rbind", lapply(exposures, function(exposure){
 # exposure <- exposures[1]

  filename <- paste0("results/MR_", exposure, "_T2DM_", outcome, ".RData")
  load(filename)

  out <- d_input
  return(out)
 }))

 outname <- paste0(c("MR", outcome, "combined"), collapse = "_")
 perform_mr(d_input = d_combined, output_name = outname)

 mr_input_data <- MendelianRandomization::mr_input(bx = d_combined$BETA_exposure,
  bxse = d_combined$SE_exposure, by = d_combined$BETA_outcome,
  byse = d_combined$SE_outcome, snps = d_combined$SNP)

 (res_ivwr <- MendelianRandomization::mr_ivw(mr_input_data,
  model = "random"))

 return(res_ivwr)
})

Sys.Date()

sessionInfo()
