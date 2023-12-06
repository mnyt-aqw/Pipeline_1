#!/usr/bin/env bash
#SBATCH -A C3SE2023-1-21 -p vera
#SBATCH -J riksmaten_test
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH --error={PATH}/job.%J.err 
#SBATCH --output={PATH}/job.%J.out 

# Unload unwanted packages and laod Nextflow
module purge
module load Nextflow

# User-configurable variables
CONTAINER="{PATH}/Pipeline_1/container/Pipeline_1.sif"  # Path to the container image file
CLUSTEROPTIONS="-A C3SE2023-1-21 -p vera"  # Cluster options
ID="unique_run_id"  # Unique ID for this run

# Input/output paths
OUTPUT_DIR_PATH="{PATH/Data/"  # Path to output directory
STOREDIR="{PATH}/storeDir/"  # Path to directory for storing long-term cache
INPUT_READS_PATH="{PATH}/*R{1,2}*.fastq.gz"  # Path to input reads for annotation.nf (supports wildcard patterns)
TRIMGALORE_PATH="{PATH}/TrimGalore/"  # Path to TrimGalore output for multiQC.nf
MGE_DB_RAW_PATH="{PATH}/mobileOG-db_beatrix-1.6.MC.faa"  # Path to raw mobileOG database file
MGE_DB_PATH="{PATH}/mobileOG/"  # Path to processed MGE database directory
ARG_DB_PATH="{PATH}/ResFinder/"  # Path to processed ARG database directory

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
    -profile "cluster" \
    --MGE_db_raw $MGE_DB_RAW_PATH \
    --MGE_db_path $MGE_DB_PATH \
    --ARG_db_path $ARG_DB_PATH \
    --storeDir $STOREDIR \
    --container $CONTAINER \
    --run_id $ID \
    -c nf_scripts/nextflow.config \
    --scratch_dir $TMPDIR