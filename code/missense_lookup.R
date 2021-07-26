# phenoscanner look-up

library(phenoscanner)
library(data.table)

res <- phenoscanner(genequery="GIP", catalogue = "GWAS", pvalue = 1)
rese <- phenoscanner(genequery="GIP", catalogue = "eQTL", pvalue = 5e-8)

table(res$results$consequence)
table(rese$results$consequence)

unique(res$results$rsid[which(res$results$consequence %in% "missense")])
unique(rese$results$rsid[which(rese$results$consequence %in% "missense")])

# check association with type 2 diabetes liability
#  RSID CHR      POS EA NEA    EAF   BETA     SE         P       N    OR
# 1: rs2291725  17 47039132  C   T 0.5376 0.0341 0.0039 4.259e-18 1335930 1.035

res_gipr <- phenoscanner(genequery="GIPR", pvalue = 1)
rese_gipr <- phenoscanner(genequery="GIPR", catalogue = "eQTL", pvalue = 5e-8)

table(res_gipr$results$consequence)
table(rese_gipr$results$consequence)

unique(res_gipr$results$rsid[
  which(res_gipr$results$consequence %in% "missense")])
unique(rese_gipr$results$rsid[
  which(rese_gipr$results$consequence %in% "missense")])

# rs143430880
# rs1800437

# associations with type 2 diabetes liability
# RSID CHR      POS EA NEA   EAF   BETA     SE         P       N    OR
# 1: rs143430880  19 46180976  G   A 0.001 0.0805 0.0613    0.1891 1023280 1.084
# 2:   rs1800437  19 46181392  C   G 0.205 0.0228 0.0047 1.209e-06 1333220 1.023

# rs143430880 has MAF < 0.01
# rs1800437 showed evidence for pleiotropy in colocalization

Sys.Date()

sessionInfo()
