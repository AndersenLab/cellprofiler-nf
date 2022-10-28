#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// INCLUDE modules
include {
    config_CP_input_toxin
} from './process_input.nf'
include {
    proc_CP_output_toxin
} from './process_output.nf'
include {
    runCP
} from './cellprofiler.nf'

/*
~ ~ ~ > * PARAMETERS SETUP
*/

pipe = "toxin-nf"
worm_model1 = "L4_N2_HB101_100w.xml"
worm_model2 = "L2L3_N2_HB101_100w.xml"
worm_model3 = "L1_N2_HB101_100w.xml"
worm_model4 = "MDHD.xml"
model_name1 = "L4_N2_HB101_100w_NonOverlappingWorms"
model_name2 = "L2L3_N2_HB101_100w_NonOverlappingWorms"
model_name3 = "L1_N2_HB101_100w_NonOverlappingWorms"
model_name4 = "MDHD_NonOverlappingWorms"

workflow toxin_workflow {
    take:
    input_data_channel

    main:
    config_cp = input_data_channel
        .combine(Channel.from(worm_model1)) // edit here
        .combine(Channel.from(worm_model2)) // edit here
        .combine(Channel.from(worm_model3)) // edit here
        .combine(Channel.from(worm_model4)) // edit here
        .combine(
            Channel.fromPath("${params.bin_dir}/config_CP_input_toxin.R")
            ) | config_CP_input_toxin
        //.view()

    // Run CellProfiler
    groups = config_CP_input_toxin.out.groups_file
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
        .combine(Channel.from(model_name3)) // HARDCODE VARIABLE NOW MAKE DEPENDENT ON PROFILE
        .combine(Channel.from(model_name4)) // HARDCODE VARIABLE NOW MAKE DEPENDENT ON PROFILE
        .combine(Channel.fromPath("${params.bin_dir}/proc_CP_output_toxin.R")) | proc_CP_output_toxin
        //.view()
}
