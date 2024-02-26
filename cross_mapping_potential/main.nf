#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
############################################
The purpose of this pipeline is to investigate cross mapping within the CARd database.

Processes:
1. Download_CARD - Downloads the CARD database for ARG reference.
2. Remake_CARD - Processes the CARD database to simulate genomic context around ARGs.
3. InSilicoSeq - Generates simulated metagenomic reads.
4. Generate_Reads_Python - An alternative method to generate simulated reads using Python.
5. BLASTN & BLASTN_python - Performs BLASTN analysis to map reads against the CARD database, facilitating the comparison of detection methods.
############################################
*/

workflow {

Download_CARD()
Remake_CARD(Download_CARD.out.combinedFasta)
InSilicoSeq(Remake_CARD.out.fasta, Remake_CARD.out.coverage)
Generate_Reads_Python(Download_CARD.out.combinedFasta)
BLASTN(Download_CARD.out.combinedFasta, InSilicoSeq.out.reads)
BLASTN_python(Download_CARD.out.combinedFasta, Generate_Reads_Python.out.reads)

}

process Download_CARD {
    publishDir "${PWD}/../../Data/cross_mapping_potential/CARD/", mode: 'copy'
    container 'docker://biocontainers/blast:v2.2.31_cv2' 

    output:
    path 'CARD_combined.fasta', emit: combinedFasta
    path '*'

    script:
    """
    # Download the CARD database files
    wget --quiet -O card_data.tar.bz2 --no-check-certificate https://card.mcmaster.ca/download/0/broadstreet-v3.1.4.tar.bz2
    tar -xvf card_data.tar.bz2 
    
    # Combine the necessary CARD database files into one
    cat nucleotide_fasta_protein_homolog_model.fasta nucleotide_fasta_protein_variant_model.fasta nucleotide_fasta_protein_overexpression_model.fasta nucleotide_fasta_protein_knockout_model.fasta > CARD_combined.fasta
    """
}

process Remake_CARD {
    publishDir "${PWD}/../../Data/cross_mapping_potential/CARD/", mode: 'copy'
    container 'docker://mnytaqw/insilicoseq:latest' 
    input:
    path CARD_database_1

    output:
    path 'insilicoseq_coverage.txt', emit: coverage
    path 'all_sequences.fasta',      emit: fasta

    script:
    """
    #!/usr/bin/env python
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    import glob
    import csv
    import subprocess


    # Initialize a list to store sequences and a dictionary to hold ARG identifiers and their simulated coverage.
    all_sequences = []
    dict = {'ARG': [], 'Coverage': []}

    # Loop through each CARD database file.
    with open("CARD_combined.fasta", "r") as handle:
        # Parse each file as FASTA format and process each sequence record.
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq)
            # Add flanking sequences to each ARG sequence to simulate genomic context.
            flanking_sequence = "A" * 100 + sequence + "A" * 100
            # Create a new SeqRecord with the modified sequence.
            new_record = SeqRecord(Seq(flanking_sequence), id=record.id, description="")
            all_sequences.append(new_record)
            dict['ARG'].append(record.id)
            dict['Coverage'].append(750)

    # Specify the output path for the combined sequences in FASTA format.
    FASTA_out = "all_sequences.fasta"
    # Write all modified sequences to a single FASTA file.
    SeqIO.write(all_sequences, FASTA_out, "fasta")

    # Specify the output path for the coverage file, which will be used by InSilicoSeq.
    COV_file = 'insilicoseq_coverage.txt'
    # Write the ARG identifiers and their coverage to a tab-delimited file.
    with open(COV_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for i in range(len(dict['ARG'])):
            writer.writerow([dict['ARG'][i], dict['Coverage'][i]])
    """
}

process InSilicoSeq {

    container 'docker://mnytaqw/insilicoseq:latest' 
    publishDir "${PWD}/../../Data/cross_mapping_potential/insilicoseq/", mode: 'copy'

    input:
    path FASTA_IN
    path COVERAGE

    output:
    path 'all_sequences_out_R{1,2}.fastq.gz', emit: reads

    shell:
    """
    iss generate --model novaseq \
        --cpus 8 \
        --seed 42 \
        --compress \
        --coverage_file ${COVERAGE} \
        --genomes ${FASTA_IN} \
        --output all_sequences_out
    """
}

process Generate_Reads_Python {

    container 'docker://mnytaqw/insilicoseq:latest' 
    publishDir "${PWD}/../../Data/cross_mapping_potential/insilicoseq/", mode: 'copy'

    input:
    path FASTA_IN

    output:
    path 'all_sequences_out_python_R1.fastq.gz', emit: reads

    shell:
    """
    #!/usr/bin/env python
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    def generate_reads(fasta_input, fasta_output, read_length=150):
        output_records = []  # List to hold generated reads

        # Iterate over each record in the input FASTA file
        for record in SeqIO.parse(fasta_input, "fasta"):
            seq_length = len(record.seq)

            # Generate reads of specified length
            for start in range(seq_length - read_length + 1):
                end = start + read_length
                read_seq = record.seq[start:end]
                read_id = f"{record.id}_{start+1}-{end}"
                read_record = SeqRecord(read_seq, id=read_id, description="")
                output_records.append(read_record)

        # Write all generated reads to the output FASTA file
        with open(fasta_output, "w") as output_handle:
            SeqIO.write(output_records, output_handle, "fasta")

    # Example usage
    fasta_input = "${FASTA_IN}"
    fasta_output = "all_sequences_out_python_R1.fastq.gz"
    generate_reads(fasta_input, fasta_output)
    """
}


process BLASTN_python {

    container 'docker://biocontainers/blast:v2.2.31_cv2' 

    publishDir "${PWD}/../../Data/cross_mapping_potential/RGI/", mode: 'copy'

    input:
    path CARD_IN
    path READS_IN

    output:
    path 'blast_out.txt' , emit: blast_out

    shell:
    """
    # Create a BLAST database from the combined CARD database file
    makeblastdb -in ${CARD_IN} -dbtype nucl -out CARD_DB


    # BLAST
    blastn -query ${READS_IN} -db CARD_DB -outfmt 6 -out blast_out_python.txt -qcov_hsp_perc 90 -perc_identity 90 -num_threads 10
    
    # remove the reads file
    rm reads.fastq
    """
}

process BLASTN {

    container 'docker://biocontainers/blast:v2.2.31_cv2' 

    publishDir "${PWD}/../../Data/cross_mapping_potential/RGI/", mode: 'copy'

    input:
    path CARD_IN
    path READS_IN

    output:
    path 'blast_out.txt' , emit: blast_out

    shell:
    """
    # Create a BLAST database from the combined CARD database file
    makeblastdb -in ${CARD_IN} -dbtype nucl -out CARD_DB

    #unzip file
    zcat -d ${READS_IN[0]} > reads.fastq

    # Convert fastq t o fasta
    awk 'NR%4==1 {print ">"substr(\$0,2)} NR%4==2 {print}' reads.fastq > reads.fasta

    # BLAST
    blastn -query reads.fasta -db CARD_DB -outfmt 6 -out blast_out.txt -qcov_hsp_perc 90 -perc_identity 90 -num_threads 10
    
    # remove the reads file
    rm reads.fastq
    """
}