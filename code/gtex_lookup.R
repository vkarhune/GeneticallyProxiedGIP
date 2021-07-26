# GTEx lookup for proxies

library(data.table)

d <- fread("data/GTEx_Analysis_v8.metasoft.txt.gz", check.names = T)

# Ensmbl ids for GIP and GIPR from Genecards
dd <- d[grepl("ENSG00000159224|ENSG00000010310", RSID),]

d_key <- fread("data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz", check.names = T)

dd[,"variant_id" := sub(",.*", "", RSID)]

dd <- dd[d_key, on = c("variant_id" = "variant_id"), nomatch = NULL]

fwrite(dd, file = "data/gip_gipr_gtex.txt", quote = F, sep = ";")

dd[rs_id_dbSNP151_GRCh38p7 %in% c("rs55936433", "rs58274617", "rs9749185", "rs55669001", "rs12709891") & PVALUE_FE < 5e-8,]

Sys.Date()

sessionInfo()
