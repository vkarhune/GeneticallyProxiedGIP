# xxx --no-save --no-restore --verbose code/mr_plots_colocalising_only.R > mr_plots_colocalising_only.out 2>&1

# Figures for only colocalising outcomes in each locus
# rm(list=ls())

.libPaths("/projappl/minmanni/project_rpackages/")

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

d_res <- do.call("rbind", lapply(loci, function(locus){
# locus <- loci[1]

 res_outcomes <- do.call("rbind", lapply(outcomes, function(outcome){
  # outcome <- outcomes[1]
  file <- paste0("results/MR_", locus, "_", exposure, "_", outcome, ".RData")

  if(!(file.exists(file))){
   cat(sprintf("File for %s and %s not found!\n", exposure, outcome))
   return(NULL)
  }

  res_outcome <- extract_results(file, exposure = exposure, outcome = outcome, methods = c("res_ivwr"),
   heterogeneity = TRUE)

  res_outcome$exposure <- exposure
  res_outcome$outcome <- outcome
  res_outcome$locus <- locus
  res_outcome$snps <- ifelse(locus %in% "gw", "", lapply(1, function(x){
   load(file)
   out <- paste0(d_input$rsid, collapse = ",")
   return(out)
  })[[1]])

  return(res_outcome)

 }))

 return(res_outcomes)

}))

setDT(d_res)

d_meta <- do.call("rbind", lapply(unique(d_res$outcome), function(outcome){

 file <- paste0("results/MR_", outcome, "_combined.RData")

 res_outcome <- extract_results(file, exposure = "T2DM", outcome = outcome, methods = c("res_ivwr"),
   heterogeneity = TRUE)
  res_outcome$exposure <- "T2DM"
  res_outcome$outcome <- outcome
  res_outcome$locus <- "meta"
  res_outcome$snps <- ""

 return(res_outcome)
}))

setDT(d_meta)

d_res <- rbindlist(list(d_res, d_meta), fill = T)

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

d_res$Locus <- factor(d_res$locus, levels = c("meta", "GIP", "GIPR", "gw"),
 labels = c("GIP and GIPR\npooled", "GIP gene", "GIPR gene", "Glycemic control\ngenerally"))

d_res[d_res$outcome %in% "HDL","outcome"] <- "HDL-C"
d_res[d_res$outcome %in% "LDL","outcome"] <- "LDL-C"

d_res$Outcome <- plyr::mapvalues(d_res$outcome, c("BMI", "CRPUKBB", "HDL-C", "HF", "TG"),
 c("Body mass index", "C-reactive protein", "High-density lipoprotein cholesterol", "Heart failure", "Triglycerides")
)

saveRDS(d_res, "results/results_mr_t2dm_colocalising_outcomes.Rds")



outcomes2
d_res2 <- do.call("rbind", lapply(outcomes2, function(outc){
 # outc <- outcomes2[1]
 loci <- c("gw", sub(".*_", "",
  names(which(as.matrix(d_colocres[outcome %in% outc, c("PP.H4.abf_GIP", "PP.H4.abf_GIPR")])[1,] > 0.8))
 ))

 res_loci <- do.call("rbind", lapply(loci, function(locus){
  file <- paste0("results/MR_", locus, "_", exposure, "_", outc, ".RData")

  if(!(file.exists(file))){
   cat(sprintf("File for %s and %s not found!\n", exposure, outc))
   return(NULL)
  }

  res_outcome <- extract_results(file, exposure = exposure, outcome = outc, methods = c("res_ivwr"),
   heterogeneity = TRUE)

  res_outcome$exposure <- exposure
  res_outcome$outcome <- outc
  res_outcome$locus <- locus
  res_outcome$snps <- ifelse(locus %in% "gw", "", lapply(1, function(x){
   load(file)
   out <- paste0(d_input$rsid, collapse = ",")
   return(out)
  })[[1]])

 return(res_outcome)
 }))

 return(res_loci)
}))

setDT(d_res2)

# switch to per halving the risk of T2DM
d_res2 <- within(d_res2, {
 beta <- -log(2)*beta
 se <- log(2)*se
 CIU <- -log(2)*cil # sign changed
 CIL <- -log(2)*ciu # sign changed

 or <- exp(beta)
 or_cil <- exp(CIL) # sign changed
 or_ciu <- exp(CIU) # sign changed
})

d_res2 <- within(d_res2, {
 cil <- CIL
 ciu <- CIU
})

d_res2$group <- ifelse(d_res2$outcome %in% c("CAD", "CKD", "HF", "strokeIS"),
 "Binary", "Continuous")

d_res2$Locus <- factor(d_res2$locus, levels = c("meta", "GIP", "GIPR", "gw"),
 labels = c("GIP and GIPR\npooled", "GIP gene", "GIPR gene", "Glycemic control\ngenerally"))

d_res2[d_res2$outcome %in% "HDL","outcome"] <- "HDL-C"
d_res2[d_res2$outcome %in% "LDL","outcome"] <- "LDL-C"

d_res2$Outcome <- plyr::mapvalues(d_res2$outcome, c("ALT", "CAD", "LDL-C", "SBP"),
 c("Alanine aminotransferase", "Coronary artery disease", "Low-density lipoprotein cholesterol", "Systolic blood pressure")
)

saveRDS(d_res2, "results/results_mr_t2dm_outcomes2.Rds")



Sys.Date()

sessionInfo()
