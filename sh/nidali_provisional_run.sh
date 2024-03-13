#!/bin/bash
#$ -l h_rt=04:00:00
#$ -j y
#$ -l mem_free=2G
#$ -l scratch=4G
#$ -r y
#$ -t 1-2303  # Iterating through batch_X where X is from 1 to 3
#$ -cwd

# Output the date and hostname for logging
date
hostname

# Use SGE_TASK_ID to iterate over the batches
BATCH_INDEX=$((SGE_TASK_ID-1))


# Adjust this varible prior to running the qsub
# Don't use trailing dash here
PROJECT_ROOT_PATH="/wynton/home/doudna/bellieny-rabelo/nidali_db"

# Set the base directories and query name
QUERY_NAME=$1
TARGET_DB_BASE="$PROJECT_ROOT_PATH/pdb_files_DAT/batch_"
QUERY_PATH="$PROJECT_ROOT_PATH/query_DAT"
OUTPUT_BASE="/wynton/home/doudna/bellieny-rabelo/projects/nidali/pyoon_bug/dali_search_output"

OUTPUT_DIR="${OUTPUT_BASE}/${QUERY_NAME}"
echo output_dir: $OUTPUT_DIR

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
for file in *.dat; do
    echo "${file%.dat}" >> "$LIST_FILE"
done

echo START NIDALI.PM RUN WITH: query: $QUERY_NAME list of files: $LIST_FILE query path: $QUERY_PATH target_db: $TARGET_DB_PATH

rm "$TARGET_DB_PATH"/dali.lock

# Run the DaliLite analysis
apptainer exec /wynton/home/doudna/bellieny-rabelo/projects/nidali/nidali.sif dali.pl --cd1 "$QUERY_NAME" --db "$LIST_FILE" --dat1 $QUERY_PATH --dat2 $TARGET_DB_PATH --oneway --outfmt "summary,alignments" --clean

#rename and move file to ouptut
mv "${QUERY_NAME}.txt" "${OUTPUT_DIR}/${QUERY_NAME}_batch_${BATCH_INDEX}.txt"

# Return to the original directory
cd "$CURRENT_DIR"
