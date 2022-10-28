#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process config_CP_input_dauer {
    publishDir "${params.out}/pipeline", 
        mode: 'copy', 
        pattern: "*.cppipe"
    publishDir "${params.out}/metadata", 
        mode: 'copy', 
        pattern: "metadata.csv"
    publishDir "${params.out}/groups", 
        mode: 'copy', 
        pattern: "groups.tsv"

    input:
    tuple file(raw_pipe), 
        val(meta_dir), 
        val(meta), 
        val(model_dir), 
        val(project), 
        val(mask), 
        val(group), 
        val(edited_pipe), 
        val(out),
        val(model1), 
        val(model2),
        file(config_script)
        

    output:
    path "*.cppipe", emit: cp_pipeline_file
    path "metadata.csv", emit: metadata_file
    path "groups.tsv", emit: groups_file
        

    script:
    template 'config_CP_input_dauer.sh'
}

process config_CP_input_toxin {
    publishDir "${params.out}/pipeline", 
        mode: 'copy', 
        pattern: "*.cppipe"
    publishDir "${params.out}/metadata", 
        mode: 'copy', 
        pattern: "metadata.csv"
    publishDir "${params.out}/groups", 
        mode: 'copy', 
        pattern: "groups.tsv"

    input:
    tuple file(raw_pipe), 
        val(meta_dir), 
        val(meta), 
        val(model_dir), 
        val(project), 
        val(mask), 
        val(group), 
        val(edited_pipe), 
        val(out), 
        val(model1), 
        val(model2), 
        val(model3), 
        val(model4),
        file(config_script)

    output:
    path "*.cppipe", emit: cp_pipeline_file
    path "metadata.csv", emit: metadata_file
    path "groups.tsv", emit: groups_file
        

    script:
    template 'config_CP_input_toxin.sh'
}