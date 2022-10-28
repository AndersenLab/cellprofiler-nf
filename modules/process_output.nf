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

    """
        # remove exisitng directories if present and make fresh
        if [ -d ${out_dir}/processed_data ]; then rm -Rf ${out_dir}/processed_data; fi
        mkdir ${out_dir}/processed_data
        if [ -d ${out_dir}/processed_images ]; then rm -Rf ${out_dir}/processed_images; fi
        mkdir ${out_dir}/processed_images
        # find .csv files, concatenate them, and write new file
        find ${out_dir}/CP_output -type f -name '${model_name1}.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/${model_name1}.csv
        
        # find .csv files, concatenate them, and write new file
        find ${out_dir}/CP_output -type f -name '${model_name2}.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/${model_name2}.csv
        
        # move all the output images to process_images directory END WITH /?
        find ${out_dir}/CP_output -name '*.png' -exec mv {} ${out_dir}/processed_images \\;
        # Process the CellProfiler output with proc_CP_output.R
        Rscript --vanilla ${proc_CP_out_script} ${out_dir}

    """
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

    """
        # remove exisitng directories if present and make fresh
        if [ -d ${out_dir}/processed_data ]; then rm -Rf ${out_dir}/processed_data; fi
        mkdir ${out_dir}/processed_data
        if [ -d ${out_dir}/processed_images ]; then rm -Rf ${out_dir}/processed_images; fi
        mkdir ${out_dir}/processed_images
        # find .csv files, concatenate them, and write new file
        find ${out_dir}/CP_output -type f -name '${model_name1}.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/${model_name1}.csv
        
        # find .csv files, concatenate them, and write new file
        find ${out_dir}/CP_output -type f -name '${model_name2}.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/${model_name2}.csv
        # find .csv files, concatenate them, and write new file
        find ${out_dir}/CP_output -type f -name '${model_name3}.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/${model_name3}.csv
        # find .csv files, concatenate them, and write new file
        find ${out_dir}/CP_output -type f -name '${model_name4}.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/${model_name4}.csv
        
        # move all the output images to process_images directory END WITH /?
        find ${out_dir}/CP_output -name '*.png' -exec mv {} ${out_dir}/processed_images \\;
        # Process the CellProfiler output with proc_CP_output.R
        Rscript --vanilla ${proc_CP_out_script} ${out_dir}

    """
}