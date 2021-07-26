# an accessory file for GIP analyses

outcome_file <- switch(outcome,
 "CAD" = "data/cad.add.160614.website.txt",
 "strokeIS" = "data/MEGASTROKE.2.IS.EUR.GC_filtered_X_nocases_het.TBL.gz",
 "HF" = "data/HERMES_Jan2019_HeartFailure_summary_data.txt.gz",
 "BMI" = "data/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz",
 "SBP" = "data/sbp_gwas.csv",
 "CRPUKBB" = "data/30710_irnt.gwas.imputed_v3.both_sexes.tsv.bgz",
 "HDL" = "data/hdl_gwas_neale.csv.gz",
 "LDL" = "data/ldl_gwas_neale.csv.gz",
 "TG" = "data/tg_gwas_neale.csv",
 "ALT" = "data/alt_ast_neale.zip",
 "CKD" = "data/CKD_overall_ALL_JW_20180223_nstud30.dbgap.txt.gz"
)

outcome_cols <- switch(outcome,
 "CAD" = c("markername", "effect_allele", "noneffect_allele", "beta", "se_dgc", "p_dgc", "effect_allele_freq"),
 "strokeIS" = c("MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "Freq1"),
 "HF" = c("SNP", "A1", "A2", "b", "se", "p", "freq"),
 "BMI" = c("SNP", "Tested_Allele", "Other_Allele", "BETA", "SE", "P", "Freq_Tested_Allele"),
 "SBP" = c("SNP", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "Freq1"),
 "CRPUKBB" = c("rsid", "EA", "NEA", "beta", "se", "pval", "minor_AF"),
 "HDL" = c("SNP", "Allele1", "Allele2", "beta", "se", "pval", "minor_AF"),
 "LDL" = c("SNP", "Allele1", "Allele2", "beta", "se", "pval", "minor_AF"),
 "TG" = c("SNP", "Allele1", "Allele2", "beta", "se", "pval", "minor_AF"),
 "ALT" = c("SNP", "A1", "A2", "beta_alt", "se_alt", "pval_alt", "minor_AF"),
 "CKD" = c("RSID", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "Freq1")
)

outcome_n <- switch(outcome,
 "CAD" = c(60801, 123504),
 "CKD" = "n_total_sum",
 "BMI" = "N",
 "HDL" = "n_complete_samples",
 "LDL" = "n_complete_samples",
 "TG" = "n_complete_samples",
 "ALT" = "n",
 "SBP" = "TotalSampleSize",
 "HF" = c(47309, 930014),
 "CRPUKBB" = "n_complete_samples",
 NULL
)

prevalence <- switch(outcome,
 "CAD" = 0.05,
 "HF" = 0.025,
 NULL
)
