#!/bin/bash

set -e
set -o pipefail
set -x

# Configure the raw pipeline for CellProfiler
awk '{gsub(/METADATA_DIR/,"${meta_dir}"); print}' ${raw_pipe} | \\
awk '{gsub(/METADATA_CSV_FILE/,"${meta}"); print}' | \\
awk '{gsub(/WORM_MODEL_DIR/,"${model_dir}"); print}' | \\
awk '{gsub(/MODEL1_XML_FILE/,"${model1}"); print}' | \\
awk '{gsub(/MODEL2_XML_FILE/,"${model2}"); print}' > pipeline.cppipe

# Configure metadata and groups for CellProfiller with config_CP_input.R
Rscript --vanilla ${config_script} ${project} ${mask} ${group} ${edited_pipe} ${out}
