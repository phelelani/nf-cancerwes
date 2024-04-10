#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS & INPUT
ref                       = file(params.ref, type: 'file')
fun_data_sources          = file(params.fun_data_sources, type: 'dir')
target_regions            = file(params.target_regions, type: 'file')
gnomad_af_only            = file(params.gnomad_af_only, type: 'file')
known_sites               = file(params.known_sites, type: 'file')
vep_cache                 = file(params.vep_cache, type: 'dir')
outdir                    = file(params.outdir, type: 'dir')
outdir.mkdir()

process run_Funcotator {
    tag { "Funcotator:${ind_id}" }
    label 'gatk'
    memory { 12.GB * task.attempt }
    cpus { 4 }
    publishDir "${outdir}/variant_annotation/${ind_id}", mode: 'copy', overwrite: true    
    
    input:
    tuple val(ind_id), path(vcf), path(index)
    
    output:
	tuple val(ind_id), path("${ind_id}.tumor-normal.oncefiltered.funcotated.vcf"),
    path("${ind_id}.tumor-normal.oncefiltered.funcotated.vcf.idx"), emit: annotated_variants
    
    script:
    mem = task.memory.toGiga() - 4
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
        Funcotator \
        --reference ${ref} \
        --ref-version hg38 \
        --variant ${vcf} \
        --data-sources-path ${fun_data_sources} \
        --output ${ind_id}.tumor-normal.oncefiltered.funcotated.vcf \
        --output-file-format VCF \
        --create-output-variant-index
    """
}

process run_vcf2maf {
    tag { "vcf2maf:${ind_id}" }
    label 'vep'
    publishDir "${outdir}/variant_annotation/${ind_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(ind_id), path(vcf), path(index)

    output:
    tuple val(ind_id), path("${ind_id}.tumor-normal.oncefiltered.funcotated.vep.maf"), emit: converted_mafs
    path("${ind_id}.tumor-normal.oncefiltered.funcotated.vep.vcf_summary.html"), emit: vep_vcf_summary
    
    """
    normal=`bcftools query -l ${vcf} | tr '\\n' '\\t' | awk '{ print \$1 }'`
    tumor=`bcftools query -l ${vcf} | tr '\\n' '\\t' | awk '{ print \$2 }'`

    perl /home/phelelani/applications/vcf2maf/vcf2maf.pl \
        --ncbi-build GRCh38 \
        --vep-path /opt/vep/src/ensembl-vep \
        --vep-data ${vep_cache} \
        --input-vcf ${vcf} \
        --output-maf ${ind_id}.tumor-normal.oncefiltered.funcotated.vep.maf \
        --ref-fasta ${ref} \
        --vcf-tumor-id \$tumor \
        --vcf-normal-id \$normal \
        --tumor-id \$tumor \
        --normal-id \$normal \
        --vep-forks ${task.cpus} --vep-overwrite --cache-version 106 --tmp-dir .
    """
}

process run_merge_maf {
    tag { "merge_maf" }
    publishDir "${outdir}/variant_annotation/MAF", mode: 'copy', overwrite: true
    
    input:
    path(maf)

    output:
    path("COHORT.somatic.maf"), emit: combined_mafs
    path("*.pdf"), emit: plots

    """
    maftools.R
    """
}

process run_multiqc_annotation {
    label 'multiqc'
    tag { 'multiqc:annotation' }
    publishDir "${outdir}/variant_annotation", mode: 'copy', overwrite: false

    input:
    path(dir)
    
    output:
    path("qc_annotation"), emit: multiqc_report
    
    """
    multiqc . --force --outdir qc_annotation
    """
}
