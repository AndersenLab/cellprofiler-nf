#!/bin/bash

set -e
set -o pipefail
set -x

pwd
ls -lh

chmod a+x ${pipeline}
export MPLCONFIGDIR=cellProfiler_tmp

mkdir cellProfiler_tmp
mkdir -p ${group}

# Run cellprofiler headless
cellprofiler -c -r \\
-i \$(pwd) \\
-p ${pipeline} \\
-g Metadata_Group=${group} \\
-o ${group} \\
-t cellProfiler_tmp