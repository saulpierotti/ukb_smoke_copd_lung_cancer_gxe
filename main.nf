#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
UKB smoke/COPD/lung cancer GxE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Author: Saul Pierotti
Mail: saul@ebi.ac.uk
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

process filter_bgen {
    label "plink2"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(bgen),
            path(bgi),
            path(sample),
            path(pheno_cov_samples)
        )
    
    output:
        tuple(
            val(meta),
            path("${meta.id}.pgen"),
            path("${meta.id}.psam"),
            path("${meta.id}.pvar"),
            eval("cat ${meta.id}.pvar | tail -n +2 | wc -l")
        )

    script:
        """
        plink2 \\
            --bgen ${meta.id}.bgen ref-first \\
            --sample ${meta.id}.sample \\
            --maf ${params.maf_filter} \\
            --keep ${pheno_cov_samples} \\
            --indiv-sort file ${pheno_cov_samples} \\
            --freq 'zs' \\
            --make-pgen \\
            --out ${meta.id}
        """
}

process run_association {
    label "pgenlibr_datatable"
    tag "${meta.id}"
    publishDir "${params.outdir}/", mode: 'copy'

    input:
        tuple(
            val(meta),
            path(pgen),
            path(psam),
            path(pvar),
            path(pheno_cov)
        )
    
    output:
        tuple(
            val(meta),
            path("${meta.id}.tsv.gwas.gz"
        )
    
    script:
        """
        run_association.R ${meta.chr} ${pheno_cov} ${meta.var_range_low} ${meta.var_range_high}
        """
}


workflow {
    Channel.fromPath ( params.input )
    .splitCsv ( header: true )
    .map { [ [ id: it.id, chr: it.id ], it.bgen, it.bgi, it.sample ] }
    .combine ( [ params.pheno_cov_id_list ] )
    .filter { it[0].id == "ukb22828_c1_b0_v3" }
    .view()
    .set { filter_bgen_ch }
    //filter_bgen ( filter_bgen_ch )

    //filter_bgen.out
    //.map {
    //    meta, pgen, psam, pvar, nvars ->
    //    def newmeta = meta.clone()
    //    newmeta.nvars = nvars as Integer
    //    return ( [ newmeta, pgen, psam, pvar ] )
    //}
    //.set { pgen_ch }

    //pgen_ch.flatMap {
    //    meta, pgen, psam, pvar ->
    //    def chunk_breakpoints = ( 1..meta.nvar ).by ( params.chunk_nvars )
    //    def tuple_list = (
    //        1..( chunk_breakpoints.size() - 1 )
    //    ).collect { [ chunk_breakpoints[ it - 1 ], chunk_breakpoints[ it ] ] }
    //    def endgroup = [ [ chunk_breakpoints[-1], n + 1 ] ]
    //    tuple_list += endgroup
    //    def res = [ [ meta ], tuple_list ].combinations()
    //    return ( res )
    //}
    //.set { ranges_ch }

    //pgen_ch.combine ( ranges_ch, by: 1 )
    //.map {
    //    meta, pgen, psam, pvar, ranges ->
    //    def newmeta = meta.clone()
    //    def low = ranges[0]
    //    def high = ranges[1]
    //    newmeta.id = meta.id + "_" + low + "-" + high
    //    newmeta.var_range_low = low
    //    newmeta.var_range_high = high
    //    return( [ newmeta, pgen, psam, pvar ] )
    //}
    //.combine ( [ params.pheno_cov ] )
    //.set { run_association_ch }
    //run_association ( run_association_ch )
}