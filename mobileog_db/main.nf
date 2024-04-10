#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
############################################
Thi pipeline downloads all bacterial proteins from Resfeq.
It then maps the mobile-og db against it.
The idea of this is to get a good reference id for each protein in the MGE db.
This so  we can get its proper DNA sequence.
############################################
*/

MGE_db = "${PWD}/../../databases/mobile-og/mobileOG-db_beatrix-1.6.MC.faa"
workflow {
    Download_refseq()
    // Create index and map
    Diamond_index(MGE_db)
    Diamond_mapping(Download_refseq.out.refseq_faa, Diamond_index.out.diamond_index)
}

process Download_refseq {
    cpus 1
    publishDir "${PWD}/../../Data/build_db/mobile-og/refseq_bacteria/", mode: 'symlink'
    container 'docker://mnytaqw/diamond:latest'

    output:
    path "refseq_bacteria.faa.gz", emit: refseq_faa

    script:
    """
    wget -r -l1 --no-parent -A 'bacteria.wp_protein.*.protein.faa.gz' ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/
    gunzip ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/* --stdout > refseq_bacteria.faa
    gzip refseq_bacteria.faa
    rm ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/* 
    """
}

process Diamond_index {
    cpus 1
    publishDir "${PWD}/../../Data/build_db/mobile-og/Diamond_index/", mode: 'symlink'
    container 'docker://mnytaqw/diamond:latest'

    input:
    path fasta

    output:
    path 'refseq_bacteria.dmnd', emit: diamond_index

    script:
    """
    diamond makedb --in ${fasta} -d refseq_bacteria
    """
}

process Diamond_mapping {
    cpus 20
    container 'docker://mnytaqw/diamond:latest' 
    publishDir "${PWD}/../../Data/build_db/mobile-og/Diamond_out/", mode: 'symlink'

    input:
    tuple val(sample_id), path(fasta)
    path index

    output:
    path "${sample_id}_${db[1]}*" , emit: diamond_out

    shell:
    """
    diamond blastp --db ${index} --query ${fasta} --out ${sample_id}.tsv --threads ${task.cpus} --id 100 --query-cover 100 --subject-cover 100 --outfmt 6 -k0
    """
}