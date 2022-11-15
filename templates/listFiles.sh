#!/bin/bash

set -e
set -o pipefail
set -x

ls -1 ${input_dir}/raw_images/* > fileList.txt