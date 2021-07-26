# Variant-outcome estimates (ESM Table 3)

# library(MendelianRandomization)
library(mrpipe)
library(data.table)

outcomes <- c("HF", "BMI", "CRPUKBB", "HDL", "TG")
outcomes_gip <- c("CAD", "ALT", "SBP")
outcomes_gipr <- c("LDL")

d1 <- lapply(outcomes, function(outcome){
# outcome <- outcomes[1]

 file1 <- paste0("results/MR_", outcome, "_combined.RData")
 file2 <- paste0("results/MR_missense_", outcome, ".RData")

 load(file1)
 d_out1 <- d_input

 load(file2)
 d_out <- rbind(d_out1, d_input)


 d_out <- d_out[,c("CHR", "POS", "rsid",
  "EA_exposure", "NEA_exposure", "BETA_exposure", "SE_exposure", "pval",
  "EAF_exposure", "EA_outcome", "NEA_outcome",
  "BETA_outcome", "SE_outcome", "P_outcome", "EAF_outcome", "bxg", "byg")]

 d_out$outcome <- outcome

 return(d_out)

})

d2 <- lapply(outcomes_gip, function(outcome){
# outcome <- outcomes[1]

 file <- paste0("results/MR_GIP_T2DM_", outcome, ".RData")
 load(file)
 d_out <- d_input

 d_out <- d_out[,c("CHR", "POS", "rsid",
  "EA_exposure", "NEA_exposure", "BETA_exposure", "SE_exposure", "pval",
  "EAF_exposure", "EA_outcome", "NEA_outcome",
  "BETA_outcome", "SE_outcome", "P_outcome", "EAF_outcome", "bxg", "byg")]

 d_out$outcome <- outcome

 return(d_out)

})

d3 <- lapply(outcomes_gipr, function(outcome){
# outcome <- outcomes[1]

 file <- paste0("results/MR_GIPR_T2DM_", outcome, ".RData")
 load(file)
 d_out <- d_input

 d_out <- d_out[,c("CHR", "POS", "rsid",
  "EA_exposure", "NEA_exposure", "BETA_exposure", "SE_exposure", "pval",
  "EAF_exposure", "EA_outcome", "NEA_outcome",
  "BETA_outcome", "SE_outcome", "P_outcome", "EAF_outcome", "bxg", "byg")]

 d_out$outcome <- outcome

 return(d_out)

})

d <- do.call("rbind", lapply(list(d1, d2, d3), function(x) do.call("rbind", x)))

d[,c("rsid", "EA_exposure", "NEA_exposure", "BETA_exposure", "EA_outcome", "BETA_outcome", "SE_outcome", "outcome")]

d1
d2
d3

Sys.Date()

sessionInfo()




