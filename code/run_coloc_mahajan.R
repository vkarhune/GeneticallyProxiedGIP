# run coloc

# rm(list=ls())

args <- commandArgs(trailingOnly = T)

pheno1 <- args[1]
pheno2 <- args[2]
chr <- as.numeric(args[3])
start <- as.numeric(args[4])
stop <- as.numeric(args[5])
window <- as.numeric(args[6])
gene <- args[7]
write <- as.logical(args[8])

# chr17:47,035,916-47,045,958

if(0){ # for testing
 pheno1 <- "T2DMMahajan"
 pheno2 <- "BMI"
 chr <- 17
 start <- 47035916
 stop <- 47045958
 window <- 0
 gene <- "GIP"
 write <- FALSE
}

library(mrpipe)
library(data.table)
library(coloc)
library(ieugwasr)

pheno1_file <- switch(pheno1,
 "T2DMMahajan" = "data/Mahajan.NatGenet2018b.T2D.European.zip",
 NULL
)

pheno1_cols <- switch(pheno1,
 "T2DMMahajan" = c("Chr", "Pos", "rsid", "EA", "NEA", "Beta", "SE", "Pvalue", "EAF"),
 NULL
)

d1 <- read_summarystats(
 phenotype = pheno1, type = "pheno1",
 file = pheno1_file,
 cols = pheno1_cols,
 no_rsid = switch(pheno1, "T2DMMahajan" = TRUE, NULL),
 chrpos_column = switch(pheno1, "T2DMMahajan" = "SNP", NULL),
 custom = switch(pheno1, "T2DMahajan" = list("'Pvalue' := as.numeric(Pvalue)"), NULL),
 n = switch(pheno1, "T2DMMahajan" = "Neff", NULL),
 keyfile = "../chrpos/"
)



# genomic position
d1 <- d1[CHR %in% chr & POS > start - 1 - window*1000 & POS < stop + 1 + window*1000,]

# exclude palindromic if MAF > 0.4
d1 <- d1[!(paste0(EA1, NEA1) %in% c("AT", "TA", "CG", "GC") & (-abs(EAF1-0.5)+0.5 > 0.4)),]

outcome <- pheno2
source("code/switch_params_gip.R")
rm(outcome)

d2 <- read_summarystats(
 phenotype = pheno2, type = "pheno2",
 file = outcome_file, cols = outcome_cols,
 is_nealelab = switch(pheno2, "CRPUKBB" = TRUE, NULL),
 chrpos_column = switch(pheno2, "CRPUKBB" = "variant", NULL),
 custom = switch(pheno2,
  "SBP" = list("'Effect' := Effect/19", "'StdErr' := StdErr/19"),
 NULL),
 n = outcome_n
 )



if(pheno2 %in% "CKD"){
 d2[,c("Ncases_pheno2", "Ncontrols_pheno2") := list(Npheno2*(64164/(64164+561055)), Npheno2*(561055/(64164+561055)))]
}

d <- d1[d2, on = c(SNPID = "SNPID"), nomatch = NULL]
d[,c("EA2", "NEA2") := list(toupper(EA2), toupper(NEA2))]
d <- d[(EA1 == EA2) | (EA1 == NEA2 & EA2 == NEA1),]

if(pheno2 %in% "BMI"){
 d <- d[!(SNPID %in% "rs6504588" & EA2 %in% "A"),]
}

if(write){

LDmat <- ld_matrix(d$SNPID,
 plink_bin = genetics.binaRies::get_plink_binary(),
 bfile = "data/EUR"
)
rownames(LDmat) <- sub("_.*", "", rownames(LDmat))
colnames(LDmat) <- sub("_.*", "", colnames(LDmat))

d <- d[SNPID %in% rownames(LDmat),]

D1 <- list(
 type = "cc",
 beta = d$BETA1,
 varbeta = d$SE1^2,
 pvalues = d$P1,
 N = 74124+824006,
 s = 74124/(74124+824006),
 MAF = -(abs(d$EAF1 - 0.5)) + 0.5,
 snp = d$SNPID
)

outcometype <- ifelse(pheno2 %in%
 c("CAD", "CKD", "HF", "strokeIS"),
 "cc", "quant")

N_outcome <- switch(pheno2,
 "ALT" = 344136, "BMI" = 484680, "CRPUKBB" = 343524,
 "SBP" = 745820, 
 "HDL" = 315133, "TG" = 343992, "LDL" = 343621,
 "CAD" = 60801+123504, "HF" = 47309+930014,
 "strokeIS" = 34217+406111,
 "CKD" = 64164+561055,
NULL
)

s_outcome <- switch(pheno2,
 "T2DM" = 0.16,
 "CKD" = 64164/(64164+561055),
 "CAD" = 60801/(60801+123504),
 "HF" = 47309/(47309+930014),
 "strokeIS" = 34217/(34217+406111),
NULL)

D2 <- list(
 type = outcometype,
 beta = d$BETA2,
 varbeta = d$SE2^2,
 pvalues = d$P2,
 N = N_outcome,
 s = s_outcome,
 MAF = -(abs(d$EAF2 - 0.5)) + 0.5,
 snp = d$SNPID
)

if(is.null(D2$s)){
 D2 <- D2[c("type", "beta", "varbeta", "pvalues", "N", "MAF", "snp")]
}

res_coloc <- coloc.signals(D1, D2, LD = LDmat,
 method="single",
 p12=1e-4, pthr=1e-6,
 r2thr=0.1
)

print(res_coloc[[1]])

d <- d[order(CHR, POS),]

cat(sprintf("Write data for colocalisation of %s and %s within %s gene (%i kb window):\n", pheno1, pheno2, gene, window))

objects_to_save <- ls()[ls() %in% c("d", "D1", "D2", "LDmat", "res_coloc")]

save(
 list = objects_to_save,
 file = paste0("results/", pheno1, "_", pheno2, "_", gene, "_coloc.RData")
)

cat(sprintf("DONE\n"))
}

Sys.Date()

sessionInfo()