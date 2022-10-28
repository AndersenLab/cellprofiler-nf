#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// INCLUDE modules
include {
    config_CP_input_dauer
} from './process_input.nf'
include {
    proc_CP_output_dauer
} from './process_output.nf'
include {
    runCP
} from './cellprofiler.nf'

/*
~ ~ ~ > * PARAMETERS SETUP
*/

pipe = "dauer-nf"
worm_model1 = "dauerMod.xml"
worm_model2 = "nondauerMod.xml"
model_name1 = "dauerMod_NonOverlappingWorms"
model_name2 = "nondauerMod_NonOverlappingWorms"

workflow dauer_workflow {
    take:
    input_data_channel

    main:
    config_cp = input_data_channel
        .combine(Channel.from(worm_model1)) // edit here
        .combine(Channel.from(worm_model2)) // edit here
        .combine(Channel.fromPath("${params.bin_dir}/config_CP_input_dauer.R"))
        .combine(Channel.from("${params.project}"))
        .combine(Channel.from("${params.well_mask}"))
        .combine(Channel.from("${params.groups}"))
        .combine(Channel.from("${params.edited_pipe}"))
        .combine(Channel.from("${params.out}")) | config_CP_input_dauer
        //.view()

     // Run CellProfiler
    groups = config_CP_input_dauer.out.groups_file
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
        .combine(Channel.from(model_name1)) // HARDCODE VARIABLE NOW MAKE DEPENDENT ON PROFILE
        .combine(Channel.from(model_name2)) // HARDCODE VARIABLE NOW MAKE DEPENDENT ON PROFILE
        .combine(Channel.fromPath("${params.bin_dir}/proc_CP_output_dauer.R"))
        //.view()

    proc_CP_output_dauer(proc_cp)
}