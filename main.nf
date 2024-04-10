#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS & INPUT NEEDED ON MAIN
workflow         = params.workflow
chroms_all       = (1..22).toList().collect { 'chr' + "${it}" } + ["chrX", "chrY"]

// MODULES
include { run_bwa; run_mark_duplicates; run_create_recalibration_table; run_recalibrate_bam; run_get_bams_samplesheet; run_multiqc_alignment } from './modules/modules-alignment.nf'
include { run_mutect2_on_normal; run_genomics_db_import; run_create_somatic_pon; run_concat_pons } from './modules/modules-create-panel-of-normals.nf'
include { run_mutect2_on_pair; run_LearnReadOrientationModel; run_GetPileupSummaries; run_CalculateContamination;
         run_FilterMutectCalls; run_get_vcfs_samplesheet } from './modules/modules-variant-calling-somatic.nf'
include { run_Funcotator; run_vcf2maf; run_merge_maf; run_multiqc_annotation } from './modules/modules-variant-annotation-somatic.nf'

// // WORKFLOWS
// ALIGN WORKFLOW
workflow ALIGNMENT {
    take:
    samples

    main:
    run_bwa(samples)
    run_mark_duplicates(run_bwa.out.raw_bam.groupTuple(by:[0,1,2], sort:true))
    run_create_recalibration_table(run_mark_duplicates.out.md_bam)
    run_recalibrate_bam(run_create_recalibration_table.out.recal_table)
    run_recalibrate_bam.out.recal_bam.collectFile() { it ->
        [ 'bams_samplesheet.tsv', "${it[0]}" + "\t" + "${it[1]}" + "\t" + "${it[2]}" + "\t" + "${it[3]}" + "\t" +  "${it[4]}" + '\n' ]
    }.set { bams_samplesheet }
    run_get_bams_samplesheet(bams_samplesheet)
    run_multiqc_alignment(run_mark_duplicates.out.md_qc.collect())
}

// CREATE PANEL OF NORMALS
workflow CREATE_PANEL_OF_NORMALS {
    take:
    samples
    chroms_all
    
    main:
    run_mutect2_on_normal(samples.filter { it -> it[2] =~ "normal" }, chroms_all)
    run_mutect2_on_normal.out.normal_calls
        .collectFile() { it -> [ 'vcf.list', "${it[3]}" + '\n' ] }
        .set { vcf_list }
    run_genomics_db_import(vcf_list, chroms_all)
    run_create_somatic_pon(run_genomics_db_import.out.chr_genomics_db)
    run_concat_pons(run_create_somatic_pon.out.chr_panel_of_normals.collect())
}

// VARIANT CALLING WORKFLOW
workflow VARIANT_CALLING_SOMATIC {
    take:
    samples
    panel_of_normals
    chroms_all
    
    main:
    run_mutect2_on_pair(samples.groupTuple(by:0, sort:true), panel_of_normals)
    run_LearnReadOrientationModel(run_mutect2_on_pair.out.unfiltered_f1r2)
    run_GetPileupSummaries(samples)
    run_CalculateContamination(run_GetPileupSummaries.out.pileup_summary.groupTuple(by:0, sort:true))
    run_mutect2_on_pair.out.unfiltered_calls
        .join(run_mutect2_on_pair.out.unfiltered_stats)
        .join(run_LearnReadOrientationModel.out.read_orientation_model)    
        .join(run_CalculateContamination.out.filtering_tables)
        .set { filtering_input }
    run_FilterMutectCalls(filtering_input)
    run_FilterMutectCalls.out.filtered_pairs.collectFile() { it ->
        [ 'vcf_samplesheet.tsv', "${it[0]}" + "\t" + "${it[1]}" + "\t" + "${it[2]}" + '\n' ]
    }.set { vcfs_samplesheet }
    run_get_vcfs_samplesheet(vcfs_samplesheet)
}

// VARIANT ANNOTATION WORKFLOW
workflow VARIANT_ANNOTATION_SOMATIC {
    take:
    samples

    main:
    run_Funcotator(samples)
    run_vcf2maf(run_Funcotator.out.annotated_variants)
    run_merge_maf(run_vcf2maf.out.converted_mafs.map { it -> it[1] }.collect())
    run_multiqc_annotation(run_vcf2maf.out.vep_vcf_summary.collect())
}

// PICK & CHOOSE WORKFLOW
workflow {
    switch (workflow) {
        case['read-alignment']:
            Channel.fromPath(params.sample_sheet_fastq)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.IndID}",
                               "${row.SampleID}",
                               "${row.TumorNormal}",
                               "${row.FastqR1}",
                               "${row.FastqR2}" ] }
                .set { samples }
            ALIGNMENT(samples)          
            break
            // =====
        case['create-panel-of-normals']:
            //
            switch (params.sample_sheet_bams) {
                case[null]:
                    sample_sheet_bams = params.outdir + "/alignment/tumor-normal-bams_samplesheet.tsv"
                    break
                default:
                    break
            }            
            Channel.fromPath(sample_sheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.IndID}",
                               "${row.SampleID}",
                               "${row.TumorNormal}",
                               "${row.BAM}",
                               "${row.Index}" ] }
                .set { samples }
            CREATE_PANEL_OF_NORMALS(samples, chroms_all)
            break
            // =====
        case['variant-calling-somatic']:
            //
            switch (params.sample_sheet_bams) {
                case[null]:
                    sample_sheet_bams = params.outdir + "/alignment/tumor-normal-bams_samplesheet.tsv"
                    break
                default:
                   break
            }            
            Channel.fromPath(sample_sheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.IndID}",
                               "${row.SampleID}",
                               "${row.TumorNormal}",
                               "${row.BAM}",
                               "${row.Index}" ] }
                .set { samples }

            //
            switch (params.panel_of_normals) {
                case [null]:
                    panel_of_normals = params.outdir + "/variant_calling/panel_of_normals/panel_of_normals.vcf.gz"
                    break
                default:
                    break
            }
            Channel.fromPath([panel_of_normals, panel_of_normals + ".tbi"])
                .toList()
                .set { panel_of_normals_vcf }
            VARIANT_CALLING_SOMATIC(samples, panel_of_normals_vcf, chroms_all)
            break
            // =====
        case['variant-annotation-somatic']:
            //
            switch (params.sample_sheet_vcf) {
                case[null]:
                    sample_sheet_vcf = params.outdir + "/variant_calling/tumor-normal.vcfs.samplesheet.tsv"
                    break
                default:
                    break
            }            
            Channel.fromPath(sample_sheet_vcf)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.IndID}",
                               "${row.VCF}",
                               "${row.Index}" ] }
                .set { samples }
            VARIANT_ANNOTATION_SOMATIC(samples)
            break 
            // =====
        // case['complete']:
        //     break
        //     // =====            
        default:
            exit 1, "NO WORKFLOW GIVEN!"
            break
            // =====
    }
}
