# rm(list=ls())

library(coloc)
library(data.table)

outcomes <- c("CKD", "CAD", "HF",
 "strokeIS",
 "ALT", "BMI", "CRPUKBB", "SBP",
 "HDL", "LDL", "TG"
)

exposures <- c("GIP", "GIPR")

# original results
d_res <- do.call("rbind", lapply(outcomes, function(outcome){
# outcome <- outcomes[1]
 out <- do.call("rbind", lapply(exposures, function(exposure){
 # exposure <- exposures[1]
  file <- paste0("results/T2DMMahajan_", outcome, "_", exposure, "_coloc.RData")
  load(file)
  out1 <- res_coloc[[1]]
  out1[,c("exposure", "outcome") := list(exposure, outcome)]
  return(out1)
 }))

 return(out)
}))

d_res[,c("exposure", "outcome", "PP.H3.abf", "PP.H4.abf")]

###





# exclude trait-specific peak
d_res_exclH3peak <- do.call("rbind", lapply(outcomes, function(outcome){
# outcome <- "BMI"
 out <- do.call("rbind", lapply(exposures, function(exposure){
 # exposure <- "GIPR"
  file <- paste0("results/T2DMMahajan_", outcome, "_", exposure, "_coloc.RData")
  load(file)

 res_coloc <- coloc.signals(D1, D2, # LD = LDmat,
  method="single",
  p12=1e-4, pthr=1e-6,
  r2thr=0.1
 )

 if(with(res_coloc[[1]], PP.H3.abf > PP.H4.abf)){

 best2 <- res_coloc[[1]][["best2"]]
 excl_indices <- which(abs(LDmat[best2,])^2 > 0.2)

 cat(sprintf("Excluding variants in LD with %s for %s in %s\n", best2, outcome, exposure))

 D1 <- lapply(D1, function(x){
  if(length(x) > 1) x <- x[-excl_indices]
  return(x)
 })

 D2 <- lapply(D2, function(x){
  if(length(x) > 1) x <- x[-excl_indices]
  return(x)
 })

 }

 res_coloc <- coloc.signals(D1, D2, # LD = LDmat,
  method="single",
  p12=1e-4, pthr=1e-6,
  r2thr=0.1
 )

  out1 <- res_coloc[[1]]
  excluded <- ifelse(exists("excl_indices"), best2, NA)
  out1[,c("exposure", "outcome", "excluded") := list(exposure, outcome, excluded)]
  return(out1)
 }))

 return(out)
}))

d_res_exclH3peak[,c("exposure", "outcome", "PP.H3.abf", "PP.H4.abf", "excluded")]

d_out <- dcast(d_res_exclH3peak[,c("exposure", "outcome", "PP.H4.abf", "excluded")], outcome ~ exposure,
 value.var = c("PP.H4.abf", "excluded"))
d_out[,"PP.H4.abf_GIP"] <- sprintf("%.2f", d_out[["PP.H4.abf_GIP"]])
d_out[,"PP.H4.abf_GIPR"] <- sprintf("%.2f", d_out[["PP.H4.abf_GIPR"]])

fwrite(d_out, "results/all_results_coloc_T2DMMahajan.txt", sep = ",")
fwrite(d_out, "results/all_results_coloc_T2DMMahajan.csv", sep = ";")

Sys.Date()

sessionInfo()
