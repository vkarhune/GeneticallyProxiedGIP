# Calculate F-statistics for ESM Table 2

# library(MendelianRandomization)
library(mrpipe)
library(data.table)

dlist <- readRDS("data/instruments_T2DM_5e-6_0.05_0_all.Rds")

d <- rbind(
 dlist[[1]][rsid %in% c("rs55936433", "rs58274617", "rs2291725"),],
 dlist[[2]][rsid %in% c("rs9749185", "rs55669001", "rs12709891"),]
)

# d[,"R2" := R2(BETA_exposure, EAF_exposure, ncases = 228499, ncontrols = 1178783)]
d[,"pseudoR2" := 2*EAF_exposure*(1-EAF_exposure)*BETA_exposure^2]
d[,"R2" := pseudoR2/(pseudoR2+pi^2/3)]
d[,"F" := Fstat(R2, n = 4/(1/228499 + 1/1178783))]

d

d_hba1c <- fread(cmd = "zcat data/30750_raw.gwas.imputed_v3.both_sexes.tsv.bgz", check.names = T)
d_hba1c[,c("CHR", "POS") := lapply(1:2, function(x) as.numeric(sapply(strsplit(variant, ":"), "[", x)))]

dd <- d_hba1c[
 (CHR %in% 17 & POS %in% c(47036042, 47036409, 47039132)) |
 (CHR %in% 19 & POS %in% c(46175416, 46177235, 46185217)),
]

dd

d_mahajan <- fread(cmd = "zcat data/Mahajan.NatGenet2018b.T2D.European.zip", check.names = T)

dd2 <- d_mahajan[
 (Chr %in% 17 & Pos %in% c(47036042, 47036409, 47039132)) |
 (Chr %in% 19 & Pos %in% c(46175416, 46177235, 46185217)),
]

dd2

#            SNP Chr      Pos EA NEA  EAF   Beta     SE  Pvalue   Neff
# 1: 17:47036042  17 47036042  A   C 0.28 -0.051 0.0071 6.4e-13 231420
# 2: 17:47036409  17 47036409  C   G 0.14  0.041 0.0093 1.2e-05 231420
# 3: 17:47039132  17 47039132  T   C 0.47 -0.043 0.0064 2.2e-11 231420
# 4: 19:46175416  19 46175416  A   G 0.11  0.069 0.0100 2.3e-11 231420
# 5: 19:46177235  19 46177235  T   C 0.89 -0.066 0.0100 6.9e-11 231420
# 6: 19:46185217  19 46185217  A   C 0.25  0.042 0.0074 2.5e-08 231420



Sys.Date()

sessionInfo()




