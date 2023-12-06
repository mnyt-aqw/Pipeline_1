#!/usr/bin/env nextflow

nextflow.enable.dsl=2
workflow {

    MultiQC(params.trimgalore_path) 

}

process MultiQC {

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