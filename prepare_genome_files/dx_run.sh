#!/bin/bash
#
# Run the ./filter_bgen.sh script on a batch of inputs (basenames).
# A file with one input per row is expected as first argument.
echo "Using ${1}"

while read id; do
    echo "Submitting ${id}"
    dx run swiss-army-knife \
        -icmd="wget https://raw.githubusercontent.com/saulpierotti/ukb_smoke_copd_lung_cancer_gxe/refs/heads/main/prepare_genome_files/filter_bgen.sh && chmod +x filter_bgen.sh && ./filter_bgen.sh ${id} smoke_lung_cancer_sample_list.txt" \
        -iin="Bulk/Imputation/UKB\ imputation\ from\ genotype/${id}.bgen" \
        -iin="Bulk/Imputation/UKB\ imputation\ from\ genotype/${id}.bgen.bgi" \
        -iin="Bulk/Imputation/UKB\ imputation\ from\ genotype/${id}.sample" \
        -iin="Bulk/Imputation/UKB\ imputation\ from\ genotype/${id}.mfi.txt" \
        -iin="pheno_cov/smoke_lung_cancer_copd_sample_list.txt" \
        --name "${id}" \
        --destination "genotypes/${id}" \
        --instance-type mem1_hdd1_x8
done < $1