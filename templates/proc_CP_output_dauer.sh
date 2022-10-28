#!/bin/bash

set -e
set -o pipefail
set -x

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