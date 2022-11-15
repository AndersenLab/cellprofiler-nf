#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process listFiles {
    container "${params.container__general}"
    input: 
    path input_dir

    output:
    path "fileList.txt"

    script:
    template 'listFiles.sh'
}

process makeMetadata {
    container "${params.container__general}"

    input:
    path in_fileList

    output:
    path "metadata.csv"

    publishDir "${params.out}/metadata", 
        mode: 'copy', 
        pattern: "metadata.csv"

    script:
    template 'makeMetadata.sh'
}
