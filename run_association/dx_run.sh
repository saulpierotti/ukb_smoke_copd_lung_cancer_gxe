#!/bin/bash
#
# Run the ./run_association.R script on a batch of inputs (basenames).
# A file with one input per row is expected as first argument.

echo "Using ${1}"
dx rm -a association/run_association.R
dx upload run_association.R --destination  association/run_association.R

while read id; do
    echo "Submitting ${id}"
    dx run swiss-army-knife \
        -icmd="Rscript run_association.R ${id} smoke_lung_cancer_copd_sample_list.txt" \
        -iin="pheno_cov/smoke_lung_cancer_copd_analysis_ready.csv" \
        -iin="genotypes/${id}.psam" \
        -iin="genotypes/${id}.pvar.zst" \
        -iin="genotypes/${id}.pgen" \
        --name "${id}" \
        --destination "association/${id}" \
        --instance-type mem1_ssd2_v2_x2
done < $1