#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS & INPUT


// CHANNELS


// MODULES
include {} from './modules/modules-alignment.nf'
include {} from './modules/modules-annotation.nf'
include {} from './modules/modules-qc.nf'
include {} from './modules/modules-variant-calling.nf'

// WORKFLOWS


// PICK & CHOOSE WORKFLOW
workflow {
    switch (mode) {
        case['qc']:
            break
            // =====
        case['alignment']:
            break
            // =====
        case['variant-calling']:
            break
            // =====
        case['annotation']:
            break
            // =====
        case['complete']:
            break
            // =====            
        default:
            exit 1, "NO MODE GIVEN!"
            break
            // =====
    }
}




    
