# UK Biobank GxE

This project explores GxE effects in the UK biobank.

## Phenotype-covariate table preparation

The phenotype-covariate table is pre-processed in the UKB RAP with the following filters:

- Outliers for heterozygosity or missing rate (field 22027) IS NOT "Yes"
- Genetic ethnic grouping (field 22006) IS "Caucasian"
- Smoking status | Instance 0 (field 20116.0) IS NOT "-3" (Prefer not to answer)
- Body mass index (BMI) | Instance 0 (field 21001.0) IS NOT NULL

I always chose instance 0 of the measurements (first assessment at recruitment) because is the one where most measures are available.
The filtered table is converted to csv using the `Table exporter` application (https://ukbiobank.dnanexus.com/app/app-GYXYYzQ9j5bFJ8x618yyb2jf).

## Other sample filters

I drop close relatives (up to 3rd degree). The list of dropped IDs is generated with the notebook `./filter_close_relatives.ipynp`.

## Genotype preparation

I use the WTCHG imputed genotypes (field 22828). I convert them to `pgen` format using `./filter_bgen.sh`.
During the conversion I also:

- Drop variants with MAF smaller than 0.01
- Drop samples dropped from the phenotype-covariate table or because they have close relatives (notebook `./get_final_id_list.ipynb`)
- Sort samples in the same order as the phenotype-covariate table
