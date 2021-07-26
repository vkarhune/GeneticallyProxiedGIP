# Create instruments

args <- commandArgs(trailingOnly = T)

exposure <- args[1]
pthresh <- as.numeric(args[2])
pthresh2 <- as.numeric(args[3])
window <- as.numeric(args[4])

if(0){
 pthresh <- 5e-6
 pthresh2 <- 0.05
 window <- 0
}

library(mrpipe)
library(data.table)
library(ieugwasr)
library(MendelianRandomization)


# READ EXPOSURE DATA
d_t2dm <- read_summarystats(
 phenotype = "T2DM", type = "exposure",
 file = "data/Vujkovic2020Discovery.txt.gz",
 cols = c("CHR", "POS", "RSID", "EA", "NEA", "BETA", "SE", "P", "EAF"),
 is_nealelab = NULL,
 chrpos_column = NULL,
 custom = "'P' := as.numeric(P)",
 n = "N"
)

str(d_t2dm)

d_hba1cukbb <- read_summarystats(
 phenotype = "HbA1cUKBB", type = "exposure",
 file = "data/30750_irnt.gwas.imputed_v3.both_sexes.tsv.bgz",
 cols = c("CHR", "POS", "rsid", "EA", "NEA", "beta", "se", "pval", "minor_AF"),
 is_nealelab = TRUE,
 chrpos_column = "variant",
 custom = NULL,
 n = "n_complete_samples"
)

str(d_hba1cukbb)
### END READ EXPOSURE

# exposures <- c("T2DM", "HbA1cUKBB")
genes <- c("GIP", "GIPR", "genome-wide")

 d1 <- switch(exposure, "T2DM" = d_t2dm, "HbA1cUKBB" = d_hba1cukbb, NULL)
 d2 <- switch(exposure, "T2DM" = d_hba1cukbb, "HbA1cUKBB" = d_t2dm, NULL)

 d2cols <- c("rsid", "EA_2", "NEA_2", "BETA_2", "SE_2", "P_2")

 setnames(d2, c("rsid", "EA_exposure", "NEA_exposure", "BETA_exposure", "SE_exposure", "pval"),
  d2cols)

 d <- d1[d2[,..d2cols], on = c("rsid" = "rsid"), nomatch = NULL]


 d_genes <- lapply(genes, function(gene){

 cat(sprintf("Running %s, %s...\n", exposure, gene))

# GIP locus:
# Genecards: chr17:47,035,916-47,045,958 (GRCh37/hg19)

# GIPR locus:
# Genecards: chr19:46,171,502-46,186,982

source <- "Genecards"
if(gene %in% "GIP"){
 chr <- 17
 start <- switch(source, 'Genecards' = 47035916, NULL)
 end <- switch(source, 'Genecards' = 47045958, NULL)
}

if(gene %in% "GIPR"){
 chr <- 19
 start <- switch(source, 'Genecards' = 46171502, NULL)
 end <- switch(source, 'Genecards' = 46186982, NULL)
}

# filter by p-value, location
if(gene %in% c("GIP", "GIPR")){
 d_gene <- d[CHR %in% chr & POS >= start - window & POS <= end + window & pval < pthresh,]
} else {
 d_gene <- d[!(CHR %in% 17 & POS >= 47035916 - window & POS <= 47045958 + window) &
  !(CHR %in% 19 & POS >= 46171502 - window & POS <= 46186982 + window),]
 d_gene <- d_gene[pval < pthresh,]
}

cat(sprintf("Filtered out %i SNPs based on genomic location and p-value threshold\n", nrow(d) - nrow(d_gene)))

# filter by the supporting phenotype
cat(sprintf("Filter %s in %s by the supporting phenotype:\n", exposure, gene))
 d_out <- d_gene[sign(BETA_exposure) == sign(BETA_2)*(2*(EA_exposure == EA_2) - 1) & P_2 < pthresh2,]

cat(sprintf("Filtered out %i SNPs based on associations with the secondary exposure\n", nrow(d_gene) - nrow(d_out)))

# output
return(d_out)

})

outfile <- paste0("results/instruments_", paste0(args, collapse = "_"), "_all.Rds")

saveRDS(d_genes, outfile)

Sys.Date()

sessionInfo()
