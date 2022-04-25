#!/usr/bin/env nextflow

// Use DSL2
nextflow.preview.dsl=2

// QUEST nextflow version message
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

/*
~ ~ ~ > * PARAMETERS SETUP
*/

// Variables
date = new Date().format('yyyyMMdd')
model_name = "NonOverlappingWorms" // JUST FOR NOW

// Parameters
params.help = null
params.debug = null
params.pipe = "dauer-nf"
params.project = null
params.groups = "plate,well"
params.data_dir = "${workflow.projectDir}/input_data" // this is different for gcp
params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp
params.well_mask = "${params.data_dir}/well_masks/wellmask_98.png"
params.out = "${params.project}/Analysis-${date}"
params.raw_pipe_dir = "${params.data_dir}/CP_pipelines"
params.raw_pipe = "${params.raw_pipe_dir}/${params.pipe}.cppipe"
params.edited_pipe = "${params.out}/pipeline/pipeline.cppipe"
params.metadata_dir = "${params.out}/metadata"
params.metadata = "metadata.csv"
params.worm_model_dir = "${params.data_dir}/worm_models"
params.worm_model = "dauer.xml"

/*
~ ~ ~ > * LOG AND HELP MESSAGE SETUP
*/

if (!params.help) {
log.info '''
C E L L P R O F I L E R - N F   P I P E L I N E
===============================================
'''
    log.info ""
    log.info "Project           = ${params.project}"
    log.info "CP Pipeline       = ${params.pipe}"
    log.info "Groups            = ${params.groups}"
    log.info "Output            = ${params.out}"
    log.info ""
    } else {
log.info '''
C E L L P R O F I L E R - N F   P I P E L I N E
===============================================
'''
    log.info "Usage:"
    log.info "The typical command for running the pipeline is as follows:"
    log.info "nextflow run main.nf --project <full path to your project directory>"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--project                      The full path to your project directory"
    log.info ""
    log.info "Optional arguments:"
    log.info "--pipe                         CellProfiler pipeline name, dauer-nf is the default now."
    log.info "--groups                       comma separated metadata groupings for CellProfiler, default is plate,well"
    log.info "--outdir                       Output directory to place files, default is project/Analysis-date"
    log.info "--help                         This usage statement."
        exit 1
    }

/*
~ ~ ~ > * WORKFLOW
*/

workflow {
    // configure inputs for CellProfiler
    config_cp = Channel.fromPath("${params.raw_pipe}")
        .combine(Channel.from("${params.metadata_dir}"))
        .combine(Channel.from("${params.metadata}"))
        .combine(Channel.from("${params.worm_model_dir}"))
        .combine(Channel.from("${params.worm_model}"))
        .combine(Channel.fromPath("${params.bin_dir}/config_CP_input.R"))
        .combine(Channel.from("${params.project}"))
        .combine(Channel.from("${params.well_mask}"))
        .combine(Channel.from("${params.groups}"))
        .combine(Channel.from("${params.edited_pipe}"))
        .combine(Channel.from("${params.out}"))
        //.view()

    config_CP_input(config_cp)

    // Run CellProfiler
    groups = config_CP_input.out.groups_file
        .splitCsv(header:true, sep: "\t")
        .map { row ->
                [row.group, file("${row.pipeline}"), file("${row.output}")]
            }
        //.view()
    
    runCP(groups)
    
    // Preprocess CellProfiler output files
    proc_cp = runCP.out.cp_output
        .last() // This ensures that all items are emitted from runCP
        .combine(Channel.from("${params.out}"))
        .combine(Channel.from(model_name)) // HARDCODE VARIABLE NOW MAKE DEPENDENT ON PROFILE
        .combine(Channel.fromPath("${params.bin_dir}/proc_CP_output.R"))
        //.view()

    proc_CP_output(proc_cp)
}

/*
~ ~ ~ > * CONFIGURE FILES FOR CELLPROFILER
*/

process config_CP_input {
    publishDir "${params.out}/pipeline", mode: 'copy', pattern: "*.cppipe"
    publishDir "${params.out}/metadata", mode: 'copy', pattern: "metadata.csv"
    publishDir "${params.out}/groups", mode: 'copy', pattern: "groups.tsv"

    input:
        tuple file(raw_pipe), val(meta_dir), val(meta), val(model_dir), val(model),
        file(config_script), val(project), val(mask), val(group), val(edited_pipe), val(out)

    output:
        path "*.cppipe", emit: cp_pipeline_file
        path "metadata.csv", emit: metadata_file
        path "groups.tsv", emit: groups_file
        

    """
        # Configure the raw pipeline for CellProfiler
        awk '{gsub(/METADATA_DIR/,"${meta_dir}"); print}' ${raw_pipe} | \\
        awk '{gsub(/METADATA_CSV_FILE/,"${meta}"); print}' | \\
        awk '{gsub(/WORM_MODEL_DIR/,"${model_dir}"); print}' | \\
        awk '{gsub(/MODEL1_XML_FILE/,"${model}"); print}' > pipeline.cppipe

        # Configure metadata and groups for CellProfiller with config_CP_input.R
        Rscript --vanilla ${config_script} ${project} ${mask} ${group} ${edited_pipe} ${out}

    """
}

/*
~ ~ ~ > * RUN CELLPROFILER
*/

process runCP {

    label "cellpro"

    input:
        tuple val(group), file(pipeline), file(output)

    output:
        stdout emit: cp_output //tuple file("*.csv"), file("*.png"), emit: cp_output

    """
        # Run cellprofiler headless
        cellprofiler -c -r -p ${pipeline} \
        -g ${group} \
        -o ${output}

    """
}

/*
~ ~ ~ > * PROCESS CELLPROFILER OUTPUTS
*/

process proc_CP_output {

    //publishDir "${params.out}/processed_data", mode: 'copy', pattern: "*.RData"
    
    input:
        tuple val(cp_output), val(out_dir), val(model_name), file(proc_CP_out_script)

    output:
        //path "*.RData", emit: cp_out_dat

    """
        # remove exisitng directories if present and make fresh
        if [ -d ${out_dir}/processed_data ]; then rm -Rf ${out_dir}/processed_data; fi
        mkdir ${out_dir}/processed_data

        if [ -d ${out_dir}/processed_images ]; then rm -Rf ${out_dir}/processed_images; fi
        mkdir ${out_dir}/processed_images

        # find .csv files, concatenate them, and write new file
        find ${out_dir}/CP_output -type f -name '${model_name}.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/${model_name}.csv
        
        # move all the output images to process_images directory END WITH /?
        find ${out_dir}/CP_output -name '*.png' -exec mv {} ${out_dir}/processed_images \\;

        # Process the CellProfiler output with proc_CP_output.R
        Rscript --vanilla ${proc_CP_out_script} ${out_dir}

    """
}

/*
~ ~ ~ > * GENERATE REPORT
*/
workflow.onComplete {

    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
    { Parameters }
    ---------------------------
    Project                                 = ${params.project}
    Result Directory                        = ${params.out}
    """

    println summary

}