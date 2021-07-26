# HbA1c associations

library(data.table)

d1 <- fread(cmd = "zcat data/30750_raw.gwas.imputed_v3.both_sexes.tsv.bgz", check.names = T)
d1 <- d1[!(low_confidence_variant),]

d2 <- fread(cmd = "zcat data/30750_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", check.names = T)
d2 <- d2[!(low_confidence_variant),]

all.equal(d1$variant, d2$variant)
all.equal(d1$minor_AF, d2$minor_AF)

d1 <- d1[minor_AF > 0.01,]
d2 <- d2[minor_AF > 0.01,]

all.equal(d1$variant, d2$variant)

betaratio <- d1$beta/d2$beta
summary(betaratio)
cor(d1$tstat, d2$tstat)

# extract those used as instruments:
# rs55936433, rs58274617, rs9749185, rs12709891, rs55669001

system("zcat data/variants.tsv.bgz | grep -P 'rs55936433|rs58274617|rs9749185|rs12709891|rs55669001'")
if(0){"
17:47036042:C:A 17      47036042        C       A       rs55936433      17:47036042_C_A 3_prime_UTR_variant     non_coding      9.85705e-01     1.00000e+00     210394  2.91248e-01     A       2.91248e-01     3.38464e-01     361194  0       181557  148880  30757   179637  4.12188e-01     4.84052e+00    4.12846e-01
17:47036409:G:C 17      47036409        G       C       rs58274617      17:47036409_G_C intron_variant  non_coding      9.72164e-01     9.99997e-01    100889   1.39661e-01     C       1.39661e-01     3.41075e-01     361193  1       267280  86937   6976    93913   2.40694e-01     1.24623e+01     2.40312e-01
19:46175416:G:A 19      46175416        G       A       rs9749185       19:46175416_G_A intron_variant  non_coding      9.97562e-01     1.00000e+00    71952    9.96030e-02     A       9.96030e-02     5.25139e-01     361194  0       292859  64718   3617    68335   1.79178e-01     1.78927e+01     1.79365e-01
19:46177235:T:C 19      46177235        T       C       rs55669001      rs55669001      intron_variant  non_coding      1.00000e+00     1.00000e+00    72214    9.99657e-02     C       9.99657e-02     7.01174e-01     361194  0       292610  64954   3630    68584   1.79831e-01     1.78937e+01     1.79945e-01
19:46185217:C:A 19      46185217        C       A       rs12709891      19:46185217_C_A 3_prime_UTR_variant     non_coding      9.98098e-01     1.00000e+00     172850  2.39276e-01     A       2.39276e-01     9.07217e-01     361194  0       209036  131466  20692   152158  3.63976e-01     6.35347e+00    3.64046e-01
"}

d_out <- d1[variant %in% c("17:47036042:C:A", "17:47036409:G:C", "19:46175416:G:A", "19:46177235:T:C", "19:46185217:C:A"),]
d_out$rsid <- c("rs55936433", "rs58274617", "rs9749185", "rs55669001", "rs12709891")

d22 <- d2[variant %in% c("17:47036042:C:A", "17:47036409:G:C", "19:46175416:G:A", "19:46177235:T:C", "19:46185217:C:A"),c("variant", "beta", "se", "tstat")]
setnames(d22, c("variant", paste0(c("beta", "se", "tstat"), "_inrt")))

d_out <- d_out[d22, on = c("variant" = "variant"), nomatch = NULL]

with(d_out, summary(beta/beta_inrt))

d_out

fwrite(d_out, "data/hba1c_raw.txt", sep = " ")

# also gw-estimates
d_gw <- readRDS("data/instruments_T2DM_5e-6_0.05_0_all.Rds")[[3]]

d_variants <- fread(cmd = "zcat data/variants.tsv.bgz", check.names = T)

d_gw <- d_gw[d_variants[,c("rsid", "variant")], on = c("rsid" = "rsid"), nomatch = NULL]
d_gw[,c("check1", "check2") := lapply(3:4, function(x) sapply(strsplit(variant, ":"), "[", x))]

d_gw <- d_gw[(EA_exposure == check1 & NEA_exposure == check2) | (EA_exposure == check2 & NEA_exposure == check1),]

d_out2 <- d_gw[d1[,c("variant", "se", "beta")], on = c("variant" = "variant"), nomatch = NULL]

fwrite(d_out2, "data/hba1c_raw_gw.txt", sep = " ")

Sys.Date()

sessionInfo()
