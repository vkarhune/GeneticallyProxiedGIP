# Collect missense and eQTL MR results

library(MendelianRandomization)
library(mrpipe)
library(data.table)
library(ggplot2)
library(ggforce)
library(patchwork)

d_colocres <- fread("results/all_results_coloc_T2DMMahajan.txt", check.names = T)
outcomes <- d_colocres[PP.H4.abf_GIP > 0.8 & PP.H4.abf_GIPR > 0.8,][["outcome"]]
outcomes2 <- setdiff(d_colocres[PP.H4.abf_GIP > 0.8 | PP.H4.abf_GIPR > 0.8,][["outcome"]], outcomes)

exposure <- "T2DM"

loci <- c("GIP", "GIPR", "gw")

d_res <- do.call("rbind", lapply(c(outcomes, outcomes2), function(outcome){

 file <- paste0("results/MR_missense_", outcome, ".RData")

 res_outcome <- extract_results(file, exposure = "T2DM", outcome = outcome, methods = c("res_ivwr"),
   heterogeneity = TRUE)
  res_outcome$exposure <- "T2DM"
  res_outcome$outcome <- outcome
  res_outcome$locus <- "missense"
  res_outcome$snps <- ""

 return(res_outcome)
}))

d_res2 <- do.call("rbind", lapply(c(outcomes, outcomes2), function(outcome){

 file <- paste0("results/MR_eQTL_", outcome, ".RData")

 res_outcome <- extract_results(file, exposure = "T2DM", outcome = outcome, methods = c("res_ivwr"),
   heterogeneity = TRUE)
  res_outcome$exposure <- "T2DM"
  res_outcome$outcome <- outcome
  res_outcome$locus <- "eQTL"
  res_outcome$snps <- ""

 return(res_outcome)
}))

d_res <- rbind(d_res, d_res2)

setDT(d_res)

# switch to per halving the risk of T2DM
d_res <- within(d_res, {
 beta <- -log(2)*beta
 se <- log(2)*se
 CIU <- -log(2)*cil # sign changed
 CIL <- -log(2)*ciu # sign changed

 or <- exp(beta)
 or_cil <- exp(CIL) # sign changed
 or_ciu <- exp(CIU) # sign changed
})

d_res <- within(d_res, {
 cil <- CIL
 ciu <- CIU
})

d_res$group <- ifelse(d_res$outcome %in% c("CAD", "CKD", "HF", "strokeIS"),
 "Binary", "Continuous")

d_res$Locus <- factor(d_res$locus, levels = c("meta", "GIP", "GIPR", "gw", "missense", "eQTL"),
 labels = c("GIP and GIPR\npooled", "GIP gene", "GIPR gene", "Glycemic control\ngenerally", "Missense variant", "eQTL"))

d_res[d_res$outcome %in% "HDL","outcome"] <- "HDL-C"
d_res[d_res$outcome %in% "LDL","outcome"] <- "LDL-C"

d_res$Outcome <- plyr::mapvalues(d_res$outcome, c("BMI", "CRPUKBB", "HDL-C", "HF", "TG"),
 c("Body mass index", "C-reactive protein", "High-density lipoprotein cholesterol", "Heart failure", "Triglycerides")
)

saveRDS(d_res, "results/results_mr_t2dm_missense_eqtl.Rds")

fwrite(d_res, "results/results_mr_t2dm_missense_eqtl.csv", sep = ";")

Sys.Date()

sessionInfo()
