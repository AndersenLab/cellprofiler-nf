#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// QUEST nextflow version message
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

// INCLUDE modules
// include {
//     dauer_workflow
// } from './modules/dauerWorkflow.nf'
// include {
//     toxin_workflow
// } from './modules/toxinWorkflow.nf'

/*
~ ~ ~ > * PARAMETERS SETUP
*/

// Variables
date = new Date().format('yyyyMMdd')

// Setup pipeline parameter
params.pipeline = null
if(!params.pipeline) {
    println """
            Error: pipeline parameter not specified. Please enter --pipeline dauer or --pipeline toxin in command.
            """
            System.exit(1)
} else if(params.pipeline != "toxin" && params.pipeline != "dauer" ) {
    println """
            Error: pipeline (${params.pipeline}) does not match expected value. Please enter either dauer or toxin.
            """
            System.exit(1)
}

date = new Date().format('yyyyMMdd')

// Help:
params.help = null
// project directory
params.project = null
// Groups to use
params.groups = "plate,well"
// directory with input data
params.data_dir = "${workflow.projectDir}/input_data" // this is different for gcp
// mask for the well
params.well_mask = "input_data/well_masks/wellmask_98.png"
// location to put output files
params.out = "${params.project}/Analysis-${date}"

params.project_name = params.project.split("/").last()
params.project_tag = "Analysis-${date}"

params.container__general = "docker://andersenlab/nemascan:20220411181933701519"
params.container__cellprofiler = "cellprofiler/cellprofiler:4.2.1"

include {
    listFiles
    makeMetadata
} from './modules/process_input.nf'
include {
    runCP
} from './modules/cellprofiler.nf'
include {
    proc_CP_output
} from './modules/process_output.nf'

if (!params.help) {
log.info '''
C E L L P R O F I L E R - N F   P I P E L I N E
===============================================
'''
    log.info ""
    log.info "Projcect          = ${params.project_name}"
    log.info "Project Dir       = ${params.project}"
    log.info "CP pipeline       = ${params.pipeline}"
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
    log.info "nextflow run main.nf --pipeline <CellProfiler pipeline to use> --project <path to your project directory>"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--project                      The path to your project directory"
    log.info "--pipeline                     The CP pipeline to use: toxin, dauer"
    log.info ""
    log.info "Optional arguments:"
    log.info "--project_name                 A name for the project.  Default detects from the project directory."
    log.info "--groups                       Comma separated metadata groupings for CellProfiler, default is plate,well"
    log.info "--out                          Output directory to place files, default is project/Analysis-{current date}"
    log.info "--help                         This usage statement."
        exit 1
    }

workflow {
    in_data_dir = Channel.fromPath("${params.project}")

    listFiles(in_data_dir)
    makeMetadata(listFiles.out)
    groups = makeMetadata.out
        .splitCsv(header: true)
        .map {
            row -> ["${row.Metadata_Group}", 
                "input_data/CP_pipelines/${params.pipeline}-nf.cppipe"]
        }
    runCP_input = groups
        .combine(in_data_dir)
        .combine(
            Channel.fromPath(
                "${workflow.projectDir}/input_data"
            )
        )
        .combine(makeMetadata.out)
    runCP(runCP_input)
    concat_outputs = runCP.out.output_files.toList()
    proc_CP_output(concat_outputs)
}

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
    Project                                 = ${params.project_name}
    Project Dir                             = ${params.project}
    Pipeline Used                           = ${params.pipeline}
    Result Directory                        = ${params.out}
    """

    println summary

}