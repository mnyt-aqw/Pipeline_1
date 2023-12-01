#!/usr/bin/env nextflow
/*
############################################

This pipeline takes raw untrimmed DNA, trims it, 
does quality control, annotates ARGs, MGEs and taxonomy
and finally, normalizes gene content.

############################################
*/

nextflow.enable.dsl=2
workflow {
    // Process build_ResFinder_db
    Build_ResFinder_db()

    // Process build_mobileOG_db
    Build_mobileOG_db(params.MGE_db_raw)
}

process Build_ResFinder_db {
    cache true
    cpus 1
    time { 10.minutes * Math.pow(1.5, task.attempt - 1) }
    errorStrategy 'terminate'
    maxRetries 3
    queue 'short'
    publishDir "${params.directory_out}ResFinder/"
    storeDir "${params.storeDir}ResFinder/"

    output:
    path "ResFinderARGs.dmnd"
    path "Combined_Resfinder.translations.faa"
    path "phenotypes.txt"
    path "*"

    script: // This script downloads the latest version of the ResFInder db, and translates it to protein.
    """
    # Download ResFinder
    git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
   
    # Combines  individual ARG files 
    cat resfinder_db/*.fsa > Combined_Resfinder.fsa

    # Remove unnecessary dir
    mv resfinder_db/* .
    rm -r resfinder_db/
    
    # Translates DNA to protein
    prodigal -i Combined_Resfinder.fsa -a Combined_Resfinder.translations.faa -p meta
    
    # Removes characters added by prodigal during translation
    sed -i 's/\\*//g' Combined_Resfinder.translations.faa
    sed -i 's/ .*\$//' Combined_Resfinder.translations.faa

    # Cluster db at 95% identity
    cd-hit -i Combined_Resfinder.translations.faa -o Combined_Resfinder.translations_cluster_95.fasta -c 0.95 -n 5 -M 0 -d 1000 -T 1

    # Created Diamond db
    diamond makedb --in Combined_Resfinder.translations_cluster_95.fasta -d ResFinderARGs
    """
}   

process Build_mobileOG_db {
    cache true
    cpus 1
    time { 10.minutes * Math.pow(1.5, task.attempt - 1) }
    errorStrategy 'terminate'
    maxRetries 3
    queue 'short'
    publishDir "${params.directory_out}mobileOG/"
    storeDir "${params.storeDir}mobileOG/"

    input:
    path mobile_OG

    output:
    path "mobileOG_db.dmnd"
    path '*'

    script:
    """
    cd-hit -i ${mobile_OG} -o mobile_OG_cluster_95.fasta -c 0.95 -n 5 -M 0 -d 1000 -T ${task.cpus} 
    diamond makedb --in mobile_OG_cluster_95.fasta -d mobileOG_db
    """
}  




