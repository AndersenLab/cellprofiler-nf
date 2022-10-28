#!/bin/bash

set -e
set -o pipefail
set -x

# Run cellprofiler headless
cellprofiler -c -r -p ${pipeline} \
-g ${group} \
-o ${output}