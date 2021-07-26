# Run MR using HbA1c exposure associations

library(data.table)
library(MendelianRandomization)
library(mrpipe)

outcomes <- c("BMI", "CRPUKBB", "HDL", "HF", "TG")
exposures <- c("GIP", "GIPR", "gw")

d_raw <- fread("data/hba1c_raw.txt", check.names = T)
d_raw_gw <- fread("data/hba1c_raw_gw.txt", check.names = T)

d_res <- do.call("rbind", lapply(outcomes, function(outcome){
 out <- do.call("rbind", lapply(exposures, function(exposure){
  file <- paste0("results/MR_", exposure, "_T2DM_", outcome, ".RData")
  load(file)

  setDT(d_input)

  if(exposure %in% "gw"){
   d_input <- d_input[d_raw_gw[,c("rsid", "beta", "se")], on = c("rsid" = "rsid"), nomatch = NULL]
  } else {
   d_input <- d_input[d_raw[,c("rsid", "beta", "se")], on = c("rsid" = "rsid"), nomatch = NULL]
  }

  all.equal(sign(d_input$BETA_2), sign(d_input$beta))

  mr_input_data <- MendelianRandomization::mr_input(bx = d_input$beta, 
   bxse = d_input$se, by = d_input$BETA_outcome, 
   byse = d_input$SE_outcome, snps = d_input$SNP)
  res_ivwr <- MendelianRandomization::mr_ivw(mr_input_data, 
    model = "random")
  d_out <- data.table(exposure = exposure, outcome = outcome,
   beta = res_ivwr@Estimate, se = res_ivwr@StdError,
   cil = res_ivwr@CILower, ciu = res_ivwr@CIUpper,
   p = res_ivwr@Pvalue)
  return(d_out)
 }))

 d_combined <- do.call("rbind", lapply(c("GIP", "GIPR"), function(exposure){

  filename <- paste0("results/MR_", exposure, "_T2DM_", outcome, ".RData")
  load(filename)

  out <- d_input
  return(out)
 }))

 setDT(d_combined)

 d_combined <- d_combined[d_raw[,c("rsid", "beta", "se")], on = c("rsid" = "rsid"), nomatch = NULL]

 all.equal(sign(d_combined$BETA_2), sign(d_combined$beta))

 mr_input_data <- MendelianRandomization::mr_input(bx = d_combined$beta,
  bxse = d_combined$se, by = d_combined$BETA_outcome,
  byse = d_combined$SE_outcome, snps = d_combined$SNP)

 (res_ivwr <- MendelianRandomization::mr_ivw(mr_input_data,
  model = "random"))

 out2 <- data.table(exposure = "meta", outcome = outcome,
  beta = res_ivwr@Estimate, se = res_ivwr@StdError,
  cil = res_ivwr@CILower, ciu = res_ivwr@CIUpper,
  p = res_ivwr@Pvalue)

 out <- rbind(out, out2)

 return(out)

}))

# switch to per decrease in HbA1c
d_res <- within(d_res, {
 beta <- -beta
 se <- se
 CIU <- -cil # sign changed
 CIL <- -ciu # sign changed

 or <- exp(beta)
 or_cil <- exp(CIL) # sign changed
 or_ciu <- exp(CIU) # sign changed
})

d_res <- within(d_res, {
 cil <- CIL
 ciu <- CIU
})

d_res

saveRDS(d_res, "results/results_mr_hba1c_colocalising_outcomes.Rds")

d_m <- d_res[exposure %in% "meta",c("outcome", "beta", "se")]
d_g <- d_res[exposure %in% "gw",c("outcome", "beta", "se")]
setnames(d_m, c("outcome", "beta_meta", "se_meta"))
setnames(d_g, c("outcome", "beta_gw", "se_gw"))

d_wide <- d_m[d_g, on = c("outcome" = "outcome"), nomatch = NULL]

d_wide[,"z" := (beta_meta - beta_gw)/sqrt(se_meta^2 + se_gw^2)]
d_wide[,"p" := 2*pnorm(abs(z), lower.tail = F)]
d_wide[,"p_onesided" := pnorm(abs(z)*(2*(sign(beta_meta) == sign(beta_gw))-1), lower.tail = F)]

d_wide



Sys.Date()

sessionInfo()

