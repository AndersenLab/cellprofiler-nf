#!/bin/bash

set -e
set -o pipefail
set -x

makeMetadata_${params.pipeline}.R \\
${in_fileList}  \\
${params.well_mask} \\
${params.groups}
echo "done"