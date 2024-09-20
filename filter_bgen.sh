#!/bin/bash
#
# Author: Saul Pierotti, European Bioinformatics Institute
#
# Filter and reorder samples in a UKB bgen file given a phenotype file from the
# eids present in such file, so that they match. Output in pgen format. Applies
# custom variant filters. It takes in input the basename of the bgen file to
# operate on, without extension, and the phenotype/covariate dataset with the
# eid column to be used. The eid column is assumed to be the first column. An
# header is assumed to be present and removed.
# The output is written with the same basename. Assumes plink2 is installed.

cat $2 | cut -d, -f1 | tail -n +2 >eids.tmp
plink2 \
    --bgen $1.bgen ref-first \
    --sample $1.sample \
    --maf 0.01 \
    --keep $EIDS \
    --indiv-sort file eids.tmp \
    --make-pgen \
    --out $1
