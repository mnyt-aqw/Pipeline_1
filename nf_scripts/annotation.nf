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
    // Collect all read and trim them
    read_pairs_ch = Channel.fromFilePairs(params.input_reads) // Collect all pair ended reads 
    trim_ch = TrimGalore(read_pairs_ch) // Trim using trim galore

    def resFinder_db_ch = tuple('ResFinder', file("${params.ARG_db_path}/ResFinderARGs.dmnd"))
    def mobileOG_db_ch = tuple('mobileOG', file("${params.MGE_db_path}/mobileOG_db.dmnd"))
    def phenotype = file("${params.ARG_db_path}/phenotypes.txt")

    // Combining the modified outputs for downstream processes
    db_ch = Channel.from(resFinder_db_ch, mobileOG_db_ch)
    //db_ch = resFinder_db_ch.concat(mobileOG_db_ch)
    diamond_ch = Diamond(trim_ch[0], db_ch) // Run Diamond

    // Taxonomic annotation
    metaxaqr_ch = MetaxaQR(trim_ch[0]) // Run MetaxaQR
    metaxaqr_ttt_ch = MetaxaQR_ttt(metaxaqr_ch[0]) // Run MetaxaQR_ttt
    metaxaQR_dc_ch = MetaxaQR_dc(metaxaqr_ttt_ch[1].collect()) // Run MetaxaQR_dc
    
    // Normalizing gene abundance
    Gene_Normalization(phenotype, diamond_ch[0].collect(), metaxaqr_ttt_ch[0].collect()) // Normalize Diamond results and output one singe file per database
    
}

process TrimGalore {

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}*R{1,2}*val*.fq.gz")
    tuple val(sample_id), path("*${sample_id}*fastqc.zip")
    path "*${sample_id}*trimming_report.txt"
  
    script:
    """
    trim_galore --paired -e 0.1 -j ${task.cpus} --phred33 -q 28 --fastqc ${reads[0]} ${reads[1]}
    """
}

process Diamond {

    input:
    tuple val(sample_id), path(trimmed_reads)
    each mode

    output:
    path "${sample_id}_${params.run_id}_${mode[0]}.tsv"

    script:
    """
    diamond blastx --query ${trimmed_reads[0]} --db ${mode[1]}  -o ${sample_id}_${params.run_id}_${mode[0]}.tsv  --id 90 --query-cover 80  --threads ${task.cpus} --outfmt 6 qseqid sseqid pident evalue length slen
    """
}

process MetaxaQR {

    input:
    tuple val(sample_id), path(trimmed_reads) 

    output:
    tuple val(sample_id), path("${sample_id}_${params.run_id}.taxonomy.txt") 
    path "*" 

    script:
    """
    # Unzip metagenomes
    zcat ${trimmed_reads[0]} > read_R1.fastq &
    zcat ${trimmed_reads[1]} > read_R2.fastq &

    # Wait for the background processes to finish
    wait

    # Run MetaxaQR
    metaxaQR -1 read_R1.fastq -2 read_R2.fastq -o ${sample_id}_${params.run_id} --cpu ${task.cpus} -g SSU -d /MetaxaQR/metaxaQR_db/SSU/mqr

    # Remove unziped metagenomes
    rm read_R*.fastq 

    # Zip fasta output
    pigz -p ${task.cpus} *.fasta 
    """
}

process MetaxaQR_ttt {

    input:
    tuple val(sample_id), path(genus_files) 

    output:
    path "${sample_id}_${params.run_id}.level_1.txt"
    path "${sample_id}_${params.run_id}.level_6.txt"
    path "*.txt"

    script:
    """
    metaxaQR_ttt -i ${genus_files} -o ${sample_id}_${params.run_id}  -m 6
    """
//  -m 6 to only go down to genus level
}

process MetaxaQR_dc {

    input:
    path level_6

    output:
    path "collected_data_${params.run_id}.txt"

    script:
    """
    metaxaQR_dc *.level_6.txt 
    mv collected_data.txt collected_data_${params.run_id}.txt
    """
}

process Gene_Normalization {

    input:
    path phenotype
    path diamond_out
    path metaxaQR_ttt
       
    output:
    path "${params.run_id}_ARG_norm.tsv"
    path "${params.run_id}_MGE_norm.tsv"

    script:
    """
    python3 ${PWD}/scripts/normalize_genes.py 
    mv ARG_norm.tsv ${params.run_id}_ARG_norm.tsv
    mv MGE_norm.tsv ${params.run_id}_MGE_norm.tsv
    """
}
