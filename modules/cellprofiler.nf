#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process runCP {
    container "${params.container__cellprofiler}"
    label "cellpro"

    input:
    tuple val(group), 
        val(pipeline),
        path(input_data),
        path("input_data"),
        path("metadata.csv")

    output:
        stdout emit: cp_output
        path "${group}", emit: output_files
    //tuple file("*.csv"), file("*.png"), emit: cp_output

    // publishDir "${params.out}/${params.project}/Analysis-${date}", 
    //     mode: 'copy', 
    //     pattern: "CP_output/${group}"

    script:
    template 'runCP.sh'
}