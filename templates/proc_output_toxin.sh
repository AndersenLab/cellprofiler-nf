#!/bin/bash

set -e
set -o pipefail
set -x

# remove exisitng directories if present and make fresh
mkdir processed_data

# find .csv files, concatenate them, and write new file
for model_name in L4_N2_HB101_100w_NonOverlappingWorms L2L3_N2_HB101_100w_NonOverlappingWorms L1_N2_HB101_100w_NonOverlappingWorms MDHD_NonOverlappingWorms; do
first_line=TRUE
for fIter in \$(ls -1 */\${model_name}.csv); do
if [ "\${first_line}" = "TRUE" ]; then 
cat \${fIter} > processed_data/\${model_name}.csv
else 
tail -n +2 \${fIter} >> processed_data/\${model_name}.csv
fi
done
done

mkdir processed_images

# move all the output images to process_images directory END WITH /?
cp */*.png processed_images/

# Process the CellProfiler output with proc_CP_output.R
proc_CP_output.R \\
${params.project_name} \\
${params.project_tag}