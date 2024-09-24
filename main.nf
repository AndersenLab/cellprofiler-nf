#!/usr/bin/env nextflow

// Use DSL2
nextflow.enable.dsl=2

// QUEST nextflow version message
if( !nextflow.version.matches('>=23.0') ) {
    println "This workflow requires Nextflow version 23.0 or greater -- You are running version $nextflow.version"
    println "On ${params.platform}, you can use `module load ${params.anaconda}; source activate ${params.softwareDir}/conda_envs/nf23_env`"
    exit 1
}

/*
~ ~ ~ > * PARAMETERS SETUP
*/

// Variables
date = new Date().format('yyyyMMdd')
// model_name = "NonOverlappingWorms" // JUST FOR NOW

// Setup pipeline parameter
if("${params.pipeline}" == "dauer") {
    pipe = "dauer-nf"
    worm_model1 = "dauerMod.xml"
    worm_model2 = "nondauerMod.xml"
    model_name1 = "dauerMod_NonOverlappingWorms"
    model_name2 = "nondauerMod_NonOverlappingWorms"
} else if( "${params.pipeline}" == "toxin" ) {
    pipe = "toxin-nf"
    worm_model1 = "L4_N2_HB101_100w.xml"
    worm_model2 = "L2L3_N2_HB101_100w.xml"
    worm_model3 = "L1_N2_HB101_100w.xml"
    worm_model4 = "MDHD.xml"
    model_name1 = "L4_N2_HB101_100w_NonOverlappingWorms"
    model_name2 = "L2L3_N2_HB101_100w_NonOverlappingWorms"
    model_name3 = "L1_N2_HB101_100w_NonOverlappingWorms"
    model_name4 = "MDHD_NonOverlappingWorms"
} else if(!params.pipeline) {
    println """
            Error: pipeline parameter not specified. Please enter --pipeline dauer or --pipeline toxin in command.
            """
    System.exit(1)
} else if("${params.pipeline}" != "toxin" || "${params.pipeline}" != "dauer" ) {
    println """
            Error: pipeline (${params.pipeline}) does not match expected value. Please enter either dauer or toxin.
            """
    System.exit(1)
}
if(!params.project) {
    println """
            Error: project parameter not specified. Please enter --project <valid_path_to_project_dir>.
            """
    System.exit(1)
}

// Configure other parameters
params.data_dir = "${workflow.projectDir}/input_data" // this is different for gcp
params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp
params.well_mask = "${params.data_dir}/well_masks/wellmask_98.png"
params.analysisDir = "Analysis-${date}"
params.outdir = params.project
params.out = "${params.outdir}/${params.analysisDir}"
params.raw_pipe_dir = "${params.data_dir}/CP_pipelines"
params.raw_pipe = "${params.raw_pipe_dir}/${pipe}.cppipe"
params.edited_pipe = "${params.out}/pipeline/pipeline.cppipe"
params.metadata_dir = "${params.out}/metadata"
params.worm_model_dir = "${params.data_dir}/worm_models"

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
    log.info "--groups                       comma separated metadata groupings for CellProfiler, default is plate,well"
    log.info "--outdir                       Output directory to place files, default is project/Analysis-{current date}"
    log.info "--help                         This usage statement."
        exit 1
    }

/*
~ ~ ~ > * WORKFLOW
*/

workflow {

    if (params.debug) {
        Channel.fromPath("${workflow.projectDir}")
            .combine(Channel.of("${params.pipeline}")) | prep_debug
        raw_name = prep_debug.out.flatten().last()
    } else {
        raw_name = Channel.fromPath("${workflow.projectDir}/raw_images/*").flatten().last()
    }

    Channel.fromPath("${params.project}")
        .combine(raw_name) | sanitize_names
    
    if("${params.pipeline}" == "dauer") {
        // configure inputs for CellProfiler
        config_cp = Channel.fromPath("${params.raw_pipe}")
            .combine(Channel.of("${params.metadata}"))
            .combine(Channel.of(worm_model1)) // edit here
            .combine(Channel.of(worm_model2)) // edit here
            .combine(Channel.fromPath("${params.bin_dir}/config_CP_input_dauer.R"))
            .combine(Channel.fromPath("${params.project}"))
            .combine(sanitize_names.out.toSortedList().flatten().last())
            .combine(Channel.of("${params.well_mask}"))
            .combine(Channel.of("${params.groups}"))
            .combine(Channel.of("${params.edited_pipe}"))
            .combine(Channel.of("${params.out}")) | config_CP_input_dauer

        // Run CellProfiler
        CP_in = config_CP_input_dauer.out.groups_file
            .splitCsv(header:true, sep: "\t")
            .map { row ->
                    [row.group, file("${row.pipeline}")]
                }
        CP_in.combine(Channel.fromPath("${params.out}"))
            .combine(Channel.fromPath("${params.project}"))
            .combine(config_CP_input_dauer.out.metadata_file)
            .combine(Channel.fromPath("${params.worm_model_dir}/${worm_model1}"))
            .combine(Channel.fromPath("${params.worm_model_dir}/${worm_model2}"))
            .combine(Channel.of(null))
            .combine(Channel.of(null)) | runCP
        
        // Compile csv files by model
        split_csv =  { item -> [ item.baseName + ".csv", item ] }
        runCP.out.cp_csv.flatten() | dos2unix
        nonoverlapping = dos2unix.out
            .collect()
            .flatten()
        
        nonoverlapping_joined = nonoverlapping
            .collectFile(split_csv, skip:1, keepHeader:true, storeDir:"${params.out}/processed_data")
            .collate(3)

        csv_files = nonoverlapping_joined
            .buffer(1)
            .last()
            .combine(Channel.of("${params.analysisDir}"))
            .combine(Channel.fromPath("${params.outdir}"))
            .combine(Channel.fromPath("${params.bin_dir}/proc_CP_output_dauer.R"))

        // Preprocess CellProfiler output files
        proc_CP_output(csv_files)
    }

    if("${params.pipeline}" == "toxin") {
        // configure inputs for CellProfiler
        config_cp = Channel.fromPath("${params.raw_pipe}")
            .combine(Channel.of("${params.metadata}"))
            .combine(Channel.of(worm_model1)) // edit here
            .combine(Channel.of(worm_model2)) // edit here
            .combine(Channel.of(worm_model3)) // edit here
            .combine(Channel.of(worm_model4)) // edit here
            .combine(Channel.fromPath("${params.bin_dir}/config_CP_input_toxin.R"))
            .combine(Channel.of("${params.project}"))
            .combine(sanitize_names.out.toSortedList().flatten().last())
            .combine(Channel.of("${params.well_mask}"))
            .combine(Channel.of("${params.groups}"))
            .combine(Channel.of("${params.edited_pipe}"))
            .combine(Channel.of("${params.out}")) | config_CP_input_toxin

        // Run CellProfiler
        CP_in = config_CP_input_toxin.out.groups_file
            .splitCsv(header:true, sep: "\t")
            .map { row ->
                    [row.group, file("${row.pipeline}")]
                }
            .combine(Channel.fromPath("${params.out}"))
            .combine(Channel.fromPath("${params.project}"))
            .combine(config_CP_input_toxin.out.metadata_file)
            .combine(Channel.fromPath("${params.worm_model_dir}/${worm_model1}"))
            .combine(Channel.fromPath("${params.worm_model_dir}/${worm_model2}"))
            .combine(Channel.fromPath("${params.worm_model_dir}/${worm_model3}"))
            .combine(Channel.fromPath("${params.worm_model_dir}/${worm_model4}"))
            
        runCP(CP_in)

        // Compile csv files by model
        split_csv =  { item -> [ item.baseName + ".csv", item ] }
        runCP.out.cp_csv.flatten() | dos2unix
        nonoverlapping = dos2unix.out
            .collect()
            .flatten()

        nonoverlapping_joined = nonoverlapping
            .collectFile(split_csv, skip:1, keepHeader:true, storeDir:"${params.out}/processed_data")
            .collate(5)

        csv_files = nonoverlapping_joined
            .buffer(1)
            .last()
            .combine(Channel.of("${params.analysisDir}"))
            .combine(Channel.fromPath("${params.outdir}"))
            .combine(Channel.fromPath("${params.bin_dir}/proc_CP_output_toxin.R"))

        // Preprocess CellProfiler output files
        proc_CP_output(csv_files)
    }
}

/*
~ ~ ~ > * CONFIGURE FILES FOR CELLPROFILER
*/

process prep_debug {

    executor "local"
    container null

    publishDir "${params.out}/..", mode: 'copy', pattern: "raw_images/*.TIF"

    input:
        tuple path(source_dir), val(pipeline)
    output:
        path "raw_images/*.TIF"

    """
    cp -r ${source_dir}/debug/*${pipeline}*/raw_images ./
    """
}

process sanitize_names {
    
    executor "local"
    container null

    publishDir "${params.out}/../renamed_images", mode: 'copyNoFollow', pattern: "*.TIF", overwrite: false

    input:
        tuple path(source_dir), file("dummyfile")
    output:
        path "*.TIF"

    """
    for I in ${source_dir}/raw_images/*; do
        if [[ \${I} =~ ^([a-z|A-Z|0-9|_|/|\\.|-]*/)?([0-9]+-[a-z|A-Z|0-9]+-p[0-9]+-m[0-9]+[X|x]_[A-Z][0-9]{2})(_w[0-9])?([A-Z|0-9|-]{36})?(\\.tif|\\.TIF)\$ ]]; then
            ln -s \${PWD}/\${I} \${BASH_REMATCH[2]}\${BASH_REMATCH[3]}.TIF
        fi
    done
    """
}

process config_CP_input_dauer {
    publishDir "${params.out}/pipeline", mode: 'copy', pattern: "*.cppipe"
    publishDir "${params.out}/metadata", mode: 'copy', pattern: "metadata.csv"
    publishDir "${params.out}/groups", mode: 'copy', pattern: "groups.tsv"

    label "xs"
    label "R"

    input:
        tuple file(raw_pipe), val(meta), val(model1), val(model2), \
              file(config_script), path(project), file(tif), val(mask), \
              val(group), val(edited_pipe), val(out)

    output:
        path "*.cppipe", emit: cp_pipeline_file
        path "metadata.csv", emit: metadata_file
        path "groups.tsv", emit: groups_file
        

    """
    # Configure the raw pipeline for CellProfiler
    awk '{gsub(/METADATA_DIR/,"."); print}' ${raw_pipe} | \\
    awk '{gsub(/METADATA_CSV_FILE/,"${meta}"); print}' | \\
    awk '{gsub(/WORM_MODEL_DIR/,"."); print}' | \\
    awk '{gsub(/MODEL1_XML_FILE/,"${model1}"); print}' | \\
    awk '{gsub(/MODEL2_XML_FILE/,"${model2}"); print}' > pipeline.cppipe

    # Configure metadata and groups for CellProfiller with config_CP_input.R
    Rscript --vanilla ${config_script} ${project} ${mask} ${group} ${edited_pipe} ${out}
    """
}

process config_CP_input_toxin {
    publishDir "${params.out}/pipeline", mode: 'copy', pattern: "*.cppipe"
    publishDir "${params.out}/metadata", mode: 'copy', pattern: "metadata.csv"
    publishDir "${params.out}/groups", mode: 'copy', pattern: "groups.tsv"

    label "xs"
    label "R"

    input:
        tuple file(raw_pipe), val(meta), val(model1), val(model2), val(model3), \
              val(model4), file(config_script), path(project), file(tif), val(mask), \
              val(group), val(edited_pipe), val(out)

    output:
        path "*.cppipe", emit: cp_pipeline_file
        path "metadata.csv", emit: metadata_file
        path "groups.tsv", emit: groups_file
        

    """
        # Configure the raw pipeline for CellProfiler
        awk '{gsub(/METADATA_DIR/,"."); print}' ${raw_pipe} | \\
        awk '{gsub(/METADATA_CSV_FILE/,"${meta}"); print}' | \\
        awk '{gsub(/WORM_MODEL_DIR/,"."); print}' | \\
        awk '{gsub(/MODEL1_XML_FILE/,"${model1}"); print}' | \\
        awk '{gsub(/MODEL2_XML_FILE/,"${model2}"); print}' | \\
        awk '{gsub(/MODEL3_XML_FILE/,"${model3}"); print}' | \\
        awk '{gsub(/MODEL4_XML_FILE/,"${model4}"); print}' > pipeline.cppipe

        # Configure metadata and groups for CellProfiller with config_CP_input.R
        Rscript --vanilla ${config_script} ${project} ${mask} ${group} ${edited_pipe} ${out}

    """
}

/*
~ ~ ~ > * RUN CELLPROFILER
*/

process runCP {

    label "cellpro"
    label "sm"

    publishDir "${params.out}/processed_images", mode: 'copy', pattern: "*.png"

    input:
        tuple val(group), file(pipeline), path(output), path(project), \
              file(metadata), file(model1), file(model2), file(model3), \
              file(model4)


    output:
        tuple path("*NonOverlappingWorms.csv"), path("WormObjects.csv"), emit: cp_csv
        path("*.png"), emit: cp_png

    """
    JAVA_OPTIONS=\$( echo "${task.memory}\" | awk '{ MEM=tolower(\$0); \
                                                    if ( MEM ~ /[0-9]* g[b]?\$/ ){ \
                                                      SUFFIX="g";  gsub(" gb","",MEM); gsub(" g","",MEM); \
                                                    } else { \
                                                      SUFFIX="m"; gsub(" mb","",MEM); gsub(" m","",MEM); \
                                                    }; \
                                                    if ( MEM/2 < 1 ) MINMEM=1; else MINMEM=MEM/2;
                                                    printf "-XX:ParallelGCThreads=1 -Xms%i%s -Xmx%i%s -Djava.io.tmpdir=${params.tmpDir}", \
                                                    MINMEM, SUFFIX, MEM, SUFFIX, MEM, SUFFIX;}' )
    export _JAVA_OPTIONS="\$JAVA_OPTIONS"
    echo "\$_JAVA_OPTIONS"
    # Run cellprofiler headless
    cellprofiler -c -r -p ${pipeline} \
    -g ${group} \
    -i ./ -o ./ \
    -t ${params.tmpDir}
    """
}

/*
~ ~ ~ > * Convert csv files from dos to unix
*/

process dos2unix {
    
    executor "local"
    container null

    input:
        path(dataDir)
    output:
       path "converted/*.csv"

    """
    mkdir -p converted
    for I in *.csv; do
        dos2unix -n \$I converted/\$I
    done
    """
}


/*
~ ~ ~ > * PROCESS CELLPROFILER OUTPUTS
*/

process proc_CP_output {

    publishDir "${params.out}/processed_data", mode: 'copy', pattern: "*.RData"

    label "md"
    label "R"
 
    input:
        tuple file(worms_csv), val(analysis_dir), path(project_dir),
              file(proc_CP_out_script)

    output:
        path "*.RData"

    """
    # Process the CellProfiler output with proc_CP_output.R
    Rscript --vanilla ${proc_CP_out_script} ${project_dir}/${analysis_dir}
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
    Pipeline Used                           = ${params.pipeline}
    Result Directory                        = ${params.out}
    """

    println summary

}