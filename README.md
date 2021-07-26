# GeneticallyProxiedGIP

Scripts for manuscript: Karhunen V et al. (2021).
Leveraging human genetic data to investigate the cardiometabolic effects of glucose-dependent insulinotropic polypeptide signaling. _Diabetologia_.
In press.

# Preparation

The documentation of R package 'mrpipe' is under construction, but the package can be downloaded via:

```
devtools::install_github("vkarhune/mrpipe")
```

For rsid-chrpos correspondence, the file for chromosome ${i} is "../chrpos/chr${i}.Rds".

The outcome summary statistics sources are given in the manuscript, and are assumed to be located in folder "data/".

The 1000G EUR correspondence was downloaded from http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz, and the population-specific files are assumed to be located in folder "data/".

The GTEx eQTL datasets are downloaded from:
https://storage.googleapis.com/gtex_analysis_v8/multi_tissue_qtl_data/GTEx_Analysis_v8.metasoft.txt.gz and https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz.


# Workflow

## Colocalization

1. Colocalization

```
for outcome in ALT BMI CAD CKD CRPUKBB HDL HF LDL SBP TG strokeIS
do
Rscript code/run_coloc_mahajan.R T2DMMahajan ${outcome} 17 47035916 47045958 0 GIP TRUE
Rscript code/run_coloc_mahajan.R T2DMMahajan ${outcome} 19 46171502 46186982 0 GIPR TRUE
done
```

2. Process colocalization results and remove pleiotropic variants

```
Rscript code/process_coloc_mahajan.R
```

## Mendelian randomization

3. Extract instruments for MR

```
Rscript code/create_instruments_all.R T2DM 5e-6 0.05 0
```


4. Run MR in GIP and GIPR separately
```
Rscript code/mr_gip_gipr_separate.R
```

5. Run MR, GIP and GIPR combined
```
Rscript code/mr_gip_gipr_combined.R
```

6. Collect MR results

```
code/collect_mr_results.R
```

## Sensitivity analysis I: HbA1c associations

7. Extract HbA1c associations
```
Rscript code/hba1c_associations.R
```

8. Run MR using HbA1c associations
```
Rscript code/mr_gip_gipr_hba1c.R
```

## Sensitivity analysis II: functionally relevant variants

9. Look-up for missense variants:
```
Rscript code/missense_lookup.R
```

10. GTEx look-up for eQTL:
```
Rscript code/gtex_lookup.R
```

11. MR using missense variant:
```
for outcome in BMI CRPUKBB HDL HF TG
do
Rscript code/mr_missense.R missense ${outcome} 5e-8 10000 varying -1
done
```
12. MR using eQTL variant:
```
for outcome in BMI CRPUKBB HDL HF TG
do
Rscript code/mr_eqtl.R eQTL ${outcome} 5e-8 10000 varying -1
done
```

13. Collect missense and eQTL MR results:
```
Rscript code/collect_mr_missense_eqtl_results.R
```

14. Calculate F-statistics:
```
Rscript code/fstat.R
```

15. Variant-outcome estimates:
```
Rscript code/variant_outcome_estimates.R
```
