# rm(list=ls())

args <- commandArgs(trailingOnly = T)

exposure <- args[1]
outcome <- args[2]
pthresh_exposure <- as.numeric(args[3])
clump_kb <- as.numeric(args[4])
# clump_r2 <- as.numeric(args[5])
clump_r2 <- args[5]
proxy_r2 <- as.numeric(args[6])

library(data.table)
library(MendelianRandomization)
library(mrpipe)
library(ieugwasr)

if(0){ # for testing
 exposure <- "missense"
 outcome <- "CRPUKBB"
 pthresh_exposure <- 5e-8
 clump_kb <- 10000
 clump_r2 <- "varying"
 proxy_r2 <- -1
 args <- c(exposure, outcome, pthresh_exposure, clump_kb, clump_r2, proxy_r2, "test")
}

cat(sprintf("Applying the following parameters:\n"))
cat(sprintf("Exposure: %s\n", exposure))
cat(sprintf("Outcome: %s\n", outcome))
cat(sprintf("P-value threshold for instruments: %s\n", pthresh_exposure))
cat(sprintf("Clumping window: %s\n", clump_kb))
cat(sprintf("Clumping r2: %s\n", clump_r2))
cat(sprintf("Proxy r2: %s\n", proxy_r2))

d_colocres <- fread("results/all_results_coloc_T2DMMahajan.txt", check.names = T)

d_input <- d_colocres[PP.H4.abf_GIP > 0.8 | PP.H4.abf_GIPR > 0.8,]

d_exposures <- readRDS("data/instruments_T2DM_5e-6_0.05_0_all.Rds")
d_exposures[[1]]$id <- "GIP_T2DM"
d_exposures[[2]]$id <- "GIPR_T2DM"
d_exposures[[3]]$id <- "gw_T2DM"

### read outcome dataset
source("code/switch_params_gip.R")

d_outcome <- read_summarystats(
 phenotype = outcome, type = "outcome",
 file = outcome_file, cols = outcome_cols,
 is_nealelab = switch(outcome, "CRPUKBB" = TRUE, NULL),
 n = outcome_n)

if(outcome %in% "SBP"){
 d_outcome[,"BETA_outcome" := BETA_outcome/19]
 d_outcome[,"SE_outcome" := SE_outcome/19]
}



de <- d_exposures[[1]]
de <- de[rsid %in% "rs2291725",]

 expo <- unique(de$id)

 cat(sprintf("Exposure: %s\n", expo))

 if(expo == outcome) return(NULL)

 prc <- proc.time()

 d <- de[d_outcome, nomatch = NULL, on = c(rsid = "rsid")]

 if(nrow(d) == 0){
  cat(sprintf("No variants available for %s in %s dataset\n", expo, outcome))
  cat(sprintf("%s done\n", expo))
  return(NULL)
 }

 dups <- d[duplicated(d$rsid)][["rsid"]]
 if(length(dups) > 0){
  warning(sprintf("Duplicate rsids in data for %s\n", expo))
  # d <- d[!(rsid %in% dups & nomatch),]
  cat(sprintf("Removing duplicates...\n"))
  dups <- d[duplicated(d$rsid)][["rsid"]]
  d[rsid %in% dups,]
 }

 d <- d[order(d[,"pval"]),]
 d <- d[!(duplicated(rsid)),]

 expotrait <- sub(".*._", "", expo)
 expogene <- sub("_.*", "", expo)


if(!(expogene %in% c("gw"))){

 ind <- which(d_colocres$outcome %in% outcome)
 d_c <- d_colocres[ind,]
 snpexclude <- d_c[[paste0("excluded_", expogene)]]

 if(snpexclude %in% "") snpexclude <- NA

 if(!(is.na(snpexclude))){
 # exclude
 d1 <- data.table(
  id  = d$id,
  rsid = d$rsid,
  pval = d$pval,
  CHR = unique(d$CHR)
 )

 d2 <- data.table(
  id  = expo,
  rsid = snpexclude, # d_exclude[which.max(PPA_2),][["id"]],
  pval = 1e-300,
  CHR = unique(d$CHR)
 )

 dd1 <- rbindlist(list(d1, d2))
 dd1 <- dd1[order(pval),]
 dd1 <- dd1[!(duplicated(rsid)),]

 d_clumped <- ld_clump(dd1, clump_kb = clump_kb, clump_r2 = 0.2,
  bfile = "data/EUR",
  plink_bin = genetics.binaRies::get_plink_binary()
 )

 d <- d[rsid %in% setdiff(d_clumped[["rsid"]], snpexclude),]

 }
}

 clump_r2 <- ifelse(grepl("gw_", expo), 0.01, 0.1)
 d <- ld_clump(d, clump_r2 = clump_r2, clump_kb = clump_kb,
  bfile = "data/EUR",
  plink_bin = genetics.binaRies::get_plink_binary()
 )

 ### check proxies ###
 cat(sprintf("Proxysearch not yet implemented!\n"))

 # back to "data.frame" object
 class(d) <- "data.frame"

 ### align effect sizes ###
 d <- align_eas(d, expos = expo)

 expo <- "missense"
 outname <- paste0(c("MR", expo, outcome), collapse = "_")
 perform_mr(d_input = d, output_name = outname)

 timer <- proc.time() - prc

 cat(sprintf("%s on %s done in %.1f seconds\n", expo, outcome, timer[[3]]))

#}))

Sys.time()

sessionInfo()
