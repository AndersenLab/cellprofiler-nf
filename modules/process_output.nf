#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process proc_CP_output_dauer {

    //publishDir "${params.out}/processed_data", mode: 'copy', pattern: "*.RData"
    
    input:
    tuple val(cp_output), 
        val(out_dir), 
        val(model_name1), 
        val(model_name2), 
        file(proc_CP_out_script)

    output:
        //path "*.RData", emit: cp_out_dat

    script:
    template 'proc_CP_output_dauer.sh'
}

process proc_CP_output_toxin {

    //publishDir "${params.out}/processed_data", mode: 'copy', pattern: "*.RData"
    
    input:
    tuple val(cp_output), 
        val(out_dir), 
        val(model_name1), 
        val(model_name2),
        val(model_name3), 
        val(model_name4), 
        file(proc_CP_out_script)

    output:
        //path "*.RData", emit: cp_out_dat

    script:
    template 'proc_CP_output_toxin.sh'
}
}