#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS & INPUT
ref            = file(params.ref, type: 'file')
known_indels_1 = file(params.known_indels_1, type: 'file')
known_indels_2 = file(params.known_indels_2, type: 'file')
dbsnp          = file(params.dbsnp, type: 'file')
outdir         = file(params.outdir, type: 'dir')
outdir.mkdir()

process run_bwa {
    tag { "bwa:${sample_id}" }
    label 'bwa_samtools'
    memory { 64.GB * task.attempt }
    cpus { "${params.bwa_threads}" }
    
    input:
    tuple val(ind_id), val(sample_id), val(tumor_normal), path(fastq_r1), path(fastq_r2)

    output:
    tuple val(ind_id), val(sample_id), val(tumor_normal), path("${sample_id}.*.bam"), emit: raw_bam
    
    script:
    nr_threads = task.cpus - 1

    """
    flowcell=`zcat ${fastq_r1} | head -n 1 | awk -F':' '{ print \$3 }'`
    lane=`zcat ${fastq_r1} | head -n 1 | awk -F':' '{ print \$4 }'`
    readgroup_info="@RG\\tID:\$flowcell.\$lane\\tLB:LIBA\\tSM:${sample_id}\\tPL:Illumina"
    bwa mem \
        -R \"\$readgroup_info\" \
        -t ${nr_threads}  \
        -K 100000000 \
        -Y \
        ${ref} ${fastq_r1} ${fastq_r2} | \
        samtools sort \
        -@ ${nr_threads} \
        -m ${params.memory_per_thread} \
        - > ${sample_id}.\"\$flowcell\".\"\$lane\".bam
    """
}

process run_mark_duplicates {
    tag { "MarkDuplicates:${sample_id}" }
    label 'gatk'
    memory { 16.GB * task.attempt }

    input:
    tuple val(ind_id), val(sample_id), val(tumor_normal), path(bam)

    output:
    tuple val(ind_id), val(sample_id), val(tumor_normal), path("${sample_id}.md.bam"), path("${sample_id}.md.bai"), emit: md_bam
    path("${sample_id}.metrics"), emit: md_qc
    
    script:
    mem = task.memory.toGiga() - 4
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms4g -Xmx${mem}g" \
        MarkDuplicates \
        --MAX_RECORDS_IN_RAM 5000 \
        --INPUT ${bam.findAll { it =~ '.bam$' }.join(' --INPUT ')} \
        --METRICS_FILE ${sample_id}.metrics \
        --TMP_DIR . \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --OUTPUT ${sample_id}.md.bam
    """
}

process run_create_recalibration_table {
    tag { "BaseRecalibrator:${sample_id}" }
    label 'gatk'
    memory { 16.GB * task.attempt }

    input:
    tuple val(ind_id), val(sample_id), val(tumor_normal), path(bam), path(index)

    output:
    tuple val(ind_id), val(sample_id), val(tumor_normal), path("${sample_id}.md.bam"), path("${sample_id}.md.bai"), path("${sample_id}.recal.table"), emit: recal_table

    script:
    mem = task.memory.toGiga() - 4

    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
        BaseRecalibrator \
        --input ${bam} \
        --output ${sample_id}.recal.table \
        --tmp-dir . \
        -R ${ref} \
        --known-sites ${dbsnp} \
        --known-sites ${known_indels_1} \
        --known-sites ${known_indels_2}
    """
}

process run_recalibrate_bam {
    tag { "ApplyBQSR:${sample_id}" }
    label 'gatk'
    memory { 16.GB * task.attempt }
    cpus { 2 }
    publishDir "${outdir}/alignment/${sample_id}", mode: 'copy', overwrite: false

    input:
    tuple val(ind_id), val(sample_id), val(tumor_normal), path(bam), path(index), file(recal_table)

    output:
    tuple val(ind_id), val(sample_id), val(tumor_normal), path("${sample_id}.md.recal.bam"), path("${sample_id}.md.recal.bai"), emit: recal_bam

    script:
    mem = task.memory.toGiga() - 4

    """
    gatk --java-options  "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
        ApplyBQSR \
        --input ${bam} \
        --output ${sample_id}.md.recal.bam \
        --tmp-dir . \
        -R ${ref} \
        --create-output-bam-index true \
        --bqsr-recal-file ${recal_table}
    """
}

process run_get_bams_samplesheet {
    tag { "write_samplesheet" }
    memory { 4.GB * task.attempt }
    publishDir "${outdir}/alignment/", mode: 'copy', overwrite: false

    input:
    path(samplesheet)

    output:
    path("tumor-normal-bams_samplesheet.tsv"), emit: bams_samplesheet

    """
    echo -e "IndID\tSampleID\tTumorNormal\tBAM\tIndex" > tumor-normal-bams_samplesheet.tsv
    cat ${samplesheet} >> tumor-normal-bams_samplesheet.tsv
    """
}

process run_multiqc_alignment {
    tag { 'multiqc:alignment' }    
    label 'multiqc'
    publishDir "${outdir}/alignment/", mode: 'copy', overwrite: false

    input:
    path(dir) 
    
    output:
    path("qc_alignment"), emit: multiqc_report
    
    """
    multiqc . --force --outdir qc_alignment
    """
}
