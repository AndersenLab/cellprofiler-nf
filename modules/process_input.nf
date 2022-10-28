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
        val(model1), 
        val(model2),
        file(config_script), 
        val(project), 
        val(mask), 
        val(group), 
        val(edited_pipe), 
        val(out)

    output:
    path "*.cppipe", emit: cp_pipeline_file
    path "metadata.csv", emit: metadata_file
    path "groups.tsv", emit: groups_file
        

    """
        # Configure the raw pipeline for CellProfiler
        awk '{gsub(/METADATA_DIR/,"${meta_dir}"); print}' ${raw_pipe} | \\
        awk '{gsub(/METADATA_CSV_FILE/,"${meta}"); print}' | \\
        awk '{gsub(/WORM_MODEL_DIR/,"${model_dir}"); print}' | \\
        awk '{gsub(/MODEL1_XML_FILE/,"${model1}"); print}' | \\
        awk '{gsub(/MODEL2_XML_FILE/,"${model2}"); print}' > pipeline.cppipe

        # Configure metadata and groups for CellProfiller with config_CP_input.R
        Rscript --vanilla ${config_script} ${project} ${mask} ${group} ${edited_pipe} ${out}

    """
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
        val(model1), 
        val(model2), 
        val(model3), 
        val(model4),
        file(config_script), 
        val(project), 
        val(mask), 
        val(group), 
        val(edited_pipe), 
        val(out)

    output:
    path "*.cppipe", emit: cp_pipeline_file
    path "metadata.csv", emit: metadata_file
    path "groups.tsv", emit: groups_file
        

    """
        # Configure the raw pipeline for CellProfiler
        awk '{gsub(/METADATA_DIR/,"${meta_dir}"); print}' ${raw_pipe} | \\
        awk '{gsub(/METADATA_CSV_FILE/,"${meta}"); print}' | \\
        awk '{gsub(/WORM_MODEL_DIR/,"${model_dir}"); print}' | \\
        awk '{gsub(/MODEL1_XML_FILE/,"${model1}"); print}' | \\
        awk '{gsub(/MODEL2_XML_FILE/,"${model2}"); print}' | \\
        awk '{gsub(/MODEL3_XML_FILE/,"${model3}"); print}' | \\
        awk '{gsub(/MODEL4_XML_FILE/,"${model4}"); print}' > pipeline.cppipe

        # Configure metadata and groups for CellProfiller with config_CP_input.R
        Rscript --vanilla ${config_script} ${project} ${mask} ${group} ${edited_pipe} ${out}

    """
}