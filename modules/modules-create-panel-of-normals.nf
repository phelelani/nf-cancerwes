#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS & INPUT
ref                       = file(params.ref, type: 'file')
target_regions            = file(params.target_regions, type: 'file')
gnomad_af_only            = file(params.gnomad_af_only, type: 'file')
outdir                    = file(params.outdir, type: 'dir')
// genomicdbimport_dir       = file(params.genomicdbimport_dir, type: 'dir')
outdir.mkdir()
// genomicdbimport_dir.mkdir()

process run_mutect2_on_normal {
    tag { "Mutect2:${sample_id}.${chr}.normal" }
    label 'gatk'
    memory { 12.GB * task.attempt }
    cpus { 4 }

    input:
    tuple val(ind_id), val(sample_id), val(tumor_normal), path(bam), path(index)
	each chr
    
    output:
	tuple val(ind_id), val(sample_id), val(tumor_normal), path("${bam.baseName}.${chr}.vcf.gz"), path("${bam.baseName}.${chr}.vcf.gz.*"),  emit: normal_calls

    script:
    mem = task.memory.toGiga() - 4
    
    """
    grep "^${chr}	" ${target_regions} > ${chr}_target_regions.bed
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
        Mutect2 \
        --reference ${ref} \
        --intervals ${chr}_target_regions.bed \
        --input ${bam} \
        --max-mnp-distance 0 \
        --output ${bam.baseName}.${chr}.vcf.gz \
        --create-output-variant-index
    """
}

process run_genomics_db_import {
    tag { "GenomicsDBImport:${chr}" }
    label 'gatk'
    memory { 16.GB * task.attempt }
    time '120h'
    publishDir "${outdir}/variant_calling/panel_of_normals/pon_db/", mode: 'copy', overwrite: true
    
    input:
    path(vcf_list)
    each chr

    output:
    tuple val("${chr}"), path("${chr}_gdb"), emit: chr_genomics_db

    script:
    mem = task.memory.toGiga() - 2

    """
    grep "^${chr}	" ${target_regions} > ${chr}_target_regions.bed
    grep ".${chr}.vcf.gz" ${vcf_list} > ${chr}.vcf.list
    gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
        GenomicsDBImport \
        --reference ${ref} \
        --intervals ${chr}_target_regions.bed \
        --merge-input-intervals \
        --variant ${chr}.vcf.list \
        --batch-size 50 \
        --reader-threads 5 \
        --genomicsdb-workspace-path ${chr}_gdb \
        --tmp-dir /tmp/
    """
}

// #${genomicdbimport_dir}/${interval[0]} \
    // grep ".${interval[0]}.vcf.gz" ${vcf_list} > ${interval[0]}.vcf.list
    // gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
    //     GenomicsDBImport \
    //     --intervals ${interval[1]} \
    //     --variant ${interval[0]}.vcf.list \
    //     --batch-size 50 \
    //     --reader-threads 5 \
    //     --genomicsdb-workspace-path . \
    //     --tmp-dir /tmp/
    // echo "Processing ${interval[1]} interval was successful!" > ${interval[1].replaceAll(':', '_')}.log
    // echo "Processing ${chr} interval was successful!" > ${chr}.log


process run_create_somatic_pon {
    tag { "CreateSomaticPoNs:${chr}" }
    label 'gatk'
    memory { 16.GB * task.attempt }
    time '120h'
    
    input:
    tuple val(chr), path(chr_genomic_db)

    output:
    tuple path("panel_of_normals.${chr}.vcf.gz"), path("panel_of_normals.${chr}.vcf.gz.tbi"), emit: chr_panel_of_normals

    script:
    mem = task.memory.toGiga() - 2
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
        CreateSomaticPanelOfNormals \
        --reference ${ref} \
        --germline-resource ${gnomad_af_only} \
        --variant gendb://${chr_genomic_db} \
        --output panel_of_normals.${chr}.vcf.gz \
        --create-output-variant-index
    """
}

process run_concat_pons {
    tag { "GatherVcfs:PoN" }
    label 'gatk'
    memory { 16.GB * task.attempt }  
    publishDir "${outdir}/variant_calling/panel_of_normals", mode: 'copy', overwrite: true

    input:
    path(chr_panel_of_normals)
  
    output:
  	tuple path("panel_of_normals.vcf.gz"), path("panel_of_normals.vcf.gz.tbi"), emit: panel_of_normals
  
    script:
    mem = task.memory.toGiga() - 4
    
    """
    find . -iname "panel_of_normals.chr1.vcf.gz" > panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr2.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr3.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr4.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr5.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr6.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr7.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr8.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr9.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr10.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr11.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr12.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr13.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr14.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr15.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr16.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr17.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr18.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr19.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr20.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr21.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chr22.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chrX.vcf.gz" >> panel_of_normals.vcf.list
    find . -iname "panel_of_normals.chrY.vcf.gz" >> panel_of_normals.vcf.list

    gatk --java-options "-XX:+UseSerialGC -Xms4g -Xmx${mem}g" \
        GatherVcfs \
        --INPUT panel_of_normals.vcf.list \
        --OUTPUT panel_of_normals.vcf.gz 
    
    tabix -p vcf panel_of_normals.vcf.gz
    """
}
