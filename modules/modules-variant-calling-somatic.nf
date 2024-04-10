#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS & INPUT
ref                       = file(params.ref, type: 'file')
target_regions            = file(params.target_regions, type: 'file')
gnomad_af_only            = file(params.gnomad_af_only, type: 'file')
known_sites               = file(params.known_sites, type: 'file')
outdir                    = file(params.outdir, type: 'dir')
outdir.mkdir()

process run_mutect2_on_pair {
    tag { "Mutect2:${ind_id}" }
    label 'gatk'
    memory { 12.GB * task.attempt }
    cpus { 4 }
    publishDir "${outdir}/variant_calling/${ind_id}", mode: 'copy', overwrite: true    
    
    input:
    tuple val(ind_id), val(sample_id), val(tumor_normal), path(bam), path(index)
    tuple path(panel_of_normals), path(index)
    
    output:
	tuple val(ind_id), path("${ind_id}.tumor-normal.unfiltered.vcf.gz"), path("${ind_id}.tumor-normal.unfiltered.vcf.gz.tbi"), emit: unfiltered_calls
    tuple val(ind_id), path("${ind_id}.tumor-normal.unfiltered.vcf.gz.stats"), emit: unfiltered_stats
    tuple val(ind_id), path("${ind_id}.tumor-normal.unfiltered.bam"), emit: unfiltered_bams
    tuple val(ind_id), path("${ind_id}.tumor-normal.unfiltered.f1r2.tar.gz"), emit: unfiltered_f1r2
    
    script:
    mem = task.memory.toGiga() - 4
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
        Mutect2 \
        -R ${ref} \
        --input ${bam.find { it =~ '_tumor.md.recal.bam$' } } \
        --tumor-sample ${sample_id.find { it =~ '_tumor$'} } \
        --input ${bam.find { it =~ '_normal.md.recal.bam$' } } \
        --normal-sample ${sample_id.find { it =~ '_normal$'} } \
        --intervals ${target_regions} \
        --panel-of-normals ${panel_of_normals} \
        --germline-resource ${gnomad_af_only} \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --f1r2-tar-gz ${ind_id}.tumor-normal.unfiltered.f1r2.tar.gz \
        --af-of-alleles-not-in-resource 0.000001 \
        --output ${ind_id}.tumor-normal.unfiltered.vcf.gz \
        --create-output-variant-index \
        --bam-output ${ind_id}.tumor-normal.unfiltered.bam \
        --create-output-bam-index
    """
}

process run_LearnReadOrientationModel {
    tag { "LearnReadOrientationModel:${ind_id}" }
    label 'gatk'
    memory { 12.GB * task.attempt }
    cpus { 4 }
    publishDir "${outdir}/variant_calling/${ind_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(ind_id), path(unfiltered_f1r2)
    
    output:
    tuple val(ind_id), path("${ind_id}.read-orientation-model.tar.gz"), emit: read_orientation_model
    
    script:
    mem = task.memory.toGiga() - 4
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
        LearnReadOrientationModel \
        --input ${unfiltered_f1r2} \
        --output ${ind_id}.read-orientation-model.tar.gz
    """
}

// RUN ON EACH TUMOR AND NORMAL BAMS SEPARATELY
process run_GetPileupSummaries {
    tag { "GetPileupSummaries:${sample_id}" }
    label 'gatk'
    memory { 12.GB * task.attempt }
    cpus { 4 }
    publishDir "${outdir}/variant_calling/${ind_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(ind_id), val(sample_id), val(tumor_normal), path(bam), path(index)
    
    output:
    tuple val(ind_id), path("${sample_id}.pileupsummaries.table"), emit: pileup_summary
    
    script:
    mem = task.memory.toGiga() - 4
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
        GetPileupSummaries \
        --input ${bam} \
        --intervals ${known_sites} \
        --variant ${known_sites} \
        --output ${sample_id}.pileupsummaries.table
    """
}

process run_CalculateContamination {
    tag { "CalculateContamination:${ind_id}" }
    label 'gatk'
    memory { 12.GB * task.attempt }
    cpus { 4 }
    publishDir "${outdir}/variant_calling/${ind_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(ind_id), path(pileup_summary)
    
    output:
    tuple val(ind_id), path("${ind_id}.tumor-normal.contamination.table"), path("${ind_id}.tumor-normal.segments.table"), emit: filtering_tables
    
    script:
    mem = task.memory.toGiga() - 4
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
        CalculateContamination \
        --input ${pileup_summary.find { it =~ '_tumor' }} \
        --matched-normal ${pileup_summary.find { it =~ '_normal' }} \
        --tumor-segmentation ${ind_id}.tumor-normal.segments.table \
        --output ${ind_id}.tumor-normal.contamination.table
    """
}

process run_FilterMutectCalls {
    tag { "FilterMutectCalls:${ind_id}" }
    label 'gatk'
    memory { 12.GB * task.attempt }
    cpus { 4 }
    publishDir "${outdir}/variant_calling/${ind_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(ind_id), path(unfiltered_vcf), path(index), path(unfiltered_stats), path(read_orientation_model), path(contamination_table), path(segments_table)
    
    output:
    tuple val(ind_id),
    path("${ind_id}.tumor-normal.oncefiltered.vcf.gz"),
    path("${ind_id}.tumor-normal.oncefiltered.vcf.gz.tbi"),
    path("${ind_id}.tumor-normal.oncefiltered.vcf.gz.filteringStats.tsv"), emit: filtered_pairs
    
    script:
    mem = task.memory.toGiga() - 4
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
        FilterMutectCalls \
        --reference ${ref} \
        --variant ${unfiltered_vcf} \
        --stats ${unfiltered_stats} \
        --ob-priors ${read_orientation_model} \
        --tumor-segmentation ${segments_table} \
        --contamination-table ${contamination_table} \
        --output ${ind_id}.tumor-normal.oncefiltered.vcf.gz \
        --create-output-variant-index
    """

}

process run_get_vcfs_samplesheet {
    tag { "write.vcfs.samplesheet" }
    memory { 4.GB * task.attempt }
    publishDir "${outdir}/variant_calling/", mode: 'copy', overwrite: false

    input:
    path(samplesheet)

    output:
    path("tumor-normal.vcfs.samplesheet.tsv"), emit: vcfs_samplesheet

    """
    echo -e "IndID\tVCF\tIndex" > tumor-normal.vcfs.samplesheet.tsv
    cat ${samplesheet} >> tumor-normal.vcfs.samplesheet.tsv
    """
}
