#!/bin/bash
#
# Run the ./filter_bgen.sh script on a batch of inputs (basenames).
# A file with one input per row is expected as first argument.

echo "Using ${1}"
dx rm -a genotypes/filter_bgen.sh
dx upload filter_bgen.sh --destination genotypes/filter_bgen.sh

while read id; do
    echo "Submitting ${id}"
    dx run swiss-army-knife \
        -icmd="bash filter_bgen.sh ${id} smoke_lung_cancer_copd_sample_list.txt" \
        -iin="genotypes/filter_bgen.sh" \
        -iin="Bulk/Imputation/UKB\ imputation\ from\ genotype/${id}.bgen" \
        -iin="Bulk/Imputation/UKB\ imputation\ from\ genotype/${id}.bgen.bgi" \
        -iin="Bulk/Imputation/UKB\ imputation\ from\ genotype/${id}.sample" \
        -iin="Bulk/Imputation/UKB\ imputation\ from\ genotype/${id}.mfi.txt" \
        -iin="pheno_cov/smoke_lung_cancer_copd_sample_list.txt" \
        --name "${id}" \
        --destination "genotypes/${id}" \
        --instance-type mem1_hdd1_x8
done < $1