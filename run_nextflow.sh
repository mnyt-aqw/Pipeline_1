#!/bin/bash

# User-configurable variables
CONTAINER="/storage/marwe/pipelines/Pipeline_1/container/Pipeline_1.sif"  # Path to the container
INPUT_READS_PATH="/storage/marwe/test/reads/sub_DEER_subsampled*R{1,2}.fq.gz"  # Path to input reads for annotation.nf
TRIMGALORE_PATH="/storage/marwe/test/Data_pipeline_test/storeDir/trimGalore/"  # Path to TrimGalore output for multiQC.nf
MGE_DB_RAW_PATH="/storage/marwe/test/pipeline/mobileOG-db_beatrix-1.6.MC.faa"  # Path to mobileOG-db_beatrix-1.6.MC.faa file
MGE_DB_PATH="/storage/marwe/test/Data_pipeline_test/storeDir/mobileOG/"  # Path to processed MGE database
ARG_DB_PATH="/storage/marwe/test/Data_pipeline_test/storeDir/ResFinder/"  # Path to processed ARG database
OUTPUT_DIR_PATH="/storage/marwe/test/Data_pipeline_test/data_out/"  # Path to output directory
STOREDIR="/storage/marwe/test/Data_pipeline_test/storeDir/"  # Path to where long term cache should be stored
EXECUTOR="server"  # Executor type (e.g., 'server' or 'cluster')
CLUSTEROPTIONS="-A C3SE2023-1-21 -p vera"  # type "-A {account name associated with the job submission} -p {partision} ". eg. "-A C3SE2021-2-3 -p vera"
ID="test_1"  # Unique ID for this run

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