#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process proc_CP_output {
    container "${params.container__general}"

    publishDir "${params.out}", 
        mode: 'copy', 
        pattern: "processed_{data,images}"
    
    input:
    path cp_output

    output:
        path "processed_data", emit: proc_dat
        path "processed_images", emit: proc_img

    script:
    template "proc_output_${params.pipeline}.sh"
}