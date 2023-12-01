#!/bin/bash

# User-configurable variables
CONTAINER="/path/to/container.sif"  # Path to the container image file
INPUT_READS_PATH="/path/to/input/reads/*.fq.gz"  # Path to input reads for annotation.nf (supports wildcard patterns)
TRIMGALORE_PATH="/path/to/trimGalore/output/"  # Path to TrimGalore output for multiQC.nf
MGE_DB_RAW_PATH="/path/to/raw/mobileOG-db.faa"  # Path to raw mobileOG database file
MGE_DB_PATH="/path/to/processed/mobileOG/db/"  # Path to processed MGE database directory
ARG_DB_PATH="/path/to/processed/ARG/db/"  # Path to processed ARG database directory
OUTPUT_DIR_PATH="/path/to/output/directory/"  # Path to output directory
STOREDIR="/path/to/storeDir/"  # Path to directory for storing long-term cache
EXECUTOR="server"  # Executor type (e.g., 'server' or 'cluster')
CLUSTEROPTIONS="cluster_options"  # Cluster options (e.g., "-A C3SE2021-2-3 -p vera")
ID="unique_run_id"  # Unique ID for this run

##########################################################
# No need to modify anything below this line

# Default value for SCRIPT_TO_RUN
SCRIPT_TO_RUN=""

# Function to display help
usage() {
    echo "Usage: $0 -s <script_name>"
    echo "  -s  Specify the Nextflow script to run. In this order(e.g. 'create_db.nf', 'annotation.nf', 'multiQC.nf')"
    echo "  -h  Display this help and exit"
}

# Parse command-line options
while getopts 's:h' flag; do
    case "${flag}" in
        s) SCRIPT_TO_RUN="${OPTARG}" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

# Check if script name is provided
if [ -z "$SCRIPT_TO_RUN" ]; then
    echo "Error: Script name is required."
    usage
    exit 1
fi

# Launch Nextflow pipeline
nextflow run nf_scripts/$SCRIPT_TO_RUN \
    --input_reads $INPUT_READS_PATH \
    --trimgalore_path $TRIMGALORE_PATH \
    --directory_out $OUTPUT_DIR_PATH \
    -profile $EXECUTOR \
    --MGE_db_raw $MGE_DB_RAW_PATH \
    --MGE_db_path $MGE_DB_PATH \
    --ARG_db_path $ARG_DB_PATH \
    --storeDir $STOREDIR \
    --cluster_options "$CLUSTEROPTIONS" \
    --container $CONTAINER \
    --run_id $ID \
    -c nf_scripts/nextflow.config