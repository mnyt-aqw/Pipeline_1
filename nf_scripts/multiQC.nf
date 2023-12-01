#!/usr/bin/env nextflow

nextflow.enable.dsl=2
workflow {

    MultiQC(params.trimgalore_path) 

}

process MultiQC {
    cache false
    cpus 1
    time { 40.minutes * Math.pow(1.5, task.attempt - 1) }
    errorStrategy 'terminate'
    maxRetries 3
    queue 'short'
    publishDir "${params.directory_out}MultiQC/"
    storeDir "${params.storeDir}MultiQC/"

	input:
	path files

	output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc -f ${params.trimgalore_path}*fastqc.zip
    """
}