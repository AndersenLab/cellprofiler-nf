#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process runCP {

    label "cellpro"

    input:
    tuple val(group), 
        file(pipeline), 
        file(output)

    output:
        stdout emit: cp_output 
    //tuple file("*.csv"), file("*.png"), emit: cp_output

    """
        # Run cellprofiler headless
        cellprofiler -c -r -p ${pipeline} \
        -g ${group} \
        -o ${output}

    """
}