#!/bin/bash
#$ -l h_rt=00:50:00
#$ -j y
#$ -l mem_free=1G
#$ -l scratch=2G
#$ -r y
#$ -t 1 #-2303  # Iterating through batch_X where X is from 1 to 3
#$ -cwd

module load Sali
module load CBI miniconda3/23.5.2-0-py311 blast
conda init bash
conda activate dali
export OMPI_MCA_osc=ucx
export OMPI_MCA_pml=ucx

# Output the date and hostname for logging
date
hostname

# Adjust this varible prior to running the qsub
# Don't use trailing dash here
PROJECT_ROOT_PATH=""

# Set the base directories and query name
QUERY_NAME="DRA2A"
TARGET_DB_BASE=$PROJECT_ROOT_PATH/"batch_"
QUERY_PATH=$PROJECT_ROOT_PATH/"query_DAT"
OUTPUT_BASE=$PROJECT_ROOT_PATH/"search_output"

OUTPUT_DIR="${OUTPUT_BASE}/${QUERY_NAME}"
echo output_dir: $OUTPUT_DIR

# Use SGE_TASK_ID to iterate over the batches
BATCH_INDEX=$((SGE_TASK_ID-1))

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define the target database path for the current task
TARGET_DB_PATH="${TARGET_DB_BASE}${BATCH_INDEX}"


# Save the current directory
CURRENT_DIR=$(pwd)

# Change to the target database directory
cd $TARGET_DB_PATH


ln -s $QUERY_PATH query_dat

# Create a list of .dat files in the target directory
LIST_FILE="list_batch_${BATCH_INDEX}.txt"
ls *.dat | sed 's/\.dat$//' > "$LIST_FILE"

echo START NIDALI.PM RUN WITH: query: $QUERY_NAME list of files: $LIST_FILE query path: $QUERY_PATH target_db: $TARGET_DB_PATH

# Run the DaliLite analysis
/wynton/home/doudna/bellieny-rabelo/projects/nidali/DaliLite.v5/bin/dali.pl --cd1 "$QUERY_NAME" --db "$LIST_FILE" --dat1 query_dat --dat2 $TARGET_DB_PATH --oneway --outfmt "summary"
