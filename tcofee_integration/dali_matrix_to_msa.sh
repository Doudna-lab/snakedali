#!/bin/bash

# This script processes DALI alignment text files into FASTA format, organizes them into directories, and performs multiple sequence alignments (MSA) using T-Coffee in two separate steps by first generating multiple query centric MSAs that are merged into one final MSA at the end. 
# It assumes the existence of a Python script that converts DALI output to FASTA and that T-Coffee is installed and accessible in the system's PATH.

# Assumptions:
# - A Python script named 'dali_out_to_fasta.py' exists at the specified SCRIPT_PATH.
# - T-Coffee is installed and available in the system's PATH.
# - The script expects the input directory containing DALI --matrix output files (DALI alignment .txt and newick_unrooted) to be provided as the first argument.
# - The script generates a directory structure within the input directory to store the intermediate and final output files.
# - The script handles directories and files according to a specific naming convention, appending '_pairwise' to output directories.



# Path to the Python script that processes dali alignment text file to fasta format
SCRIPT_PATH="/Users/scripts/dali_out_to_fasta.py" # Ensure this is the correct path to the Python script

# Directory containing the .txt files
INPUT_DIR="$1" # Your directory with the DALI .txt file

# Base directory where the output folders will be created
BASE_OUTPUT_DIR="$INPUT_DIR/tcoffee" # Adjust as needed for output

# Creates tcoffee directory
mkdir -p "$BASE_OUTPUT_DIR"

# Postfix to be added to each output folder
POSTFIX="_pairwise"

# Loop through all .txt files in the input directory
for file in "$INPUT_DIR"/*.txt; do
    # Extract the filename without the extension
    filename=$(basename -- "$file")
    filename="${filename%.txt}"

    # Create a new folder for this file in the base output directory with the postfix
    OUTPUT_DIR="$BASE_OUTPUT_DIR/${filename}${POSTFIX}"
    mkdir -p "$OUTPUT_DIR"

    echo "Processing $file into $OUTPUT_DIR"
    python "$SCRIPT_PATH" -i "$file" -o "$OUTPUT_DIR"
done

echo "All files processed."

# Mother directory containing the _pairwise directories
MOTHER_DIR="$BASE_OUTPUT_DIR"

# Get one directory level up from MOTHER_DIR
PARENT_DIR=$(dirname "$MOTHER_DIR")

# Define TREE_PATH by appending the filename to the parent directory
TREE_PATH="${PARENT_DIR}/newick_unrooted"

# Directory to store all .aln files
ALN_OUTPUT_DIR="$MOTHER_DIR/pairwise_aln"
mkdir -p "$ALN_OUTPUT_DIR"  # Create it if it doesn't exist

# Iterate over directories ending with "_pairwise" within the output directory
for DIR_PATH in "$MOTHER_DIR"/*_pairwise; do
    if [ -d "$DIR_PATH" ]; then  # Check if it's a directory
        echo "Processing directory: $DIR_PATH"

        # Extract the directory name without the _pairwise suffix
        DIR_NAME=$(basename "$DIR_PATH" "_pairwise")

        # Initialize ALN_FILES as an empty string
        ALN_FILES=""

        # Find all .fasta files and process them
        for fasta_file in "$DIR_PATH"/*.fasta; do
            fasta_basename=$(basename "$fasta_file")

            # Skip the file if it matches the pattern A_A.fasta, where A is the directory name
            if [[ "$fasta_basename" != "${DIR_NAME}_${DIR_NAME}.fasta" ]]; then
                # If ALN_FILES is not empty, append a comma
                if [[ -n "$ALN_FILES" ]]; then
                    ALN_FILES+=","
                fi
                # Add the basename of the file to ALN_FILES
                ALN_FILES+="$fasta_basename"
            fi
        done

        # Proceed only if ALN_FILES is not empty
        if [[ -n "$ALN_FILES" ]]; then
            # Save the list to a temporary file
            echo "$ALN_FILES" > "$DIR_PATH/temp_fasta_list.txt"

            # Run t_coffee using the list of files
            cd "$DIR_PATH"
            t_coffee -aln $(cat temp_fasta_list.txt) -output fasta_aln -matrix=blosum30mt -usetree= $TREE_PATH

            # Clean up the temporary file
            rm "$DIR_PATH/temp_fasta_list.txt"

            # Find and copy .aln files to ALN_OUTPUT_DIR
            find "$DIR_PATH" -type f -name '*.fasta_aln' -exec cp {} "$ALN_OUTPUT_DIR" \;

            cd -  # Go back to the previous directory
        fi
    fi
done

echo "All _pairwise directories processed."


# Final step: Run t_coffee on all generated .aln files in the pairwise_aln directory

# Directory to store the final .aln file
FINAL_MSA_DIR="$ALN_OUTPUT_DIR/final_msa"
mkdir -p "$FINAL_MSA_DIR"  # Create it if it doesn't exist

# Initialize a variable to hold all .aln file paths
ALN_FILES_FOR_FINAL_MSA=""

# Find all .aln files in the pairwise_aln directory and add them to the list
for aln_file in "$ALN_OUTPUT_DIR"/*.fasta_aln; do
    # If ALN_FILES_FOR_FINAL_MSA is not empty, append a space
    if [[ -n "$ALN_FILES_FOR_FINAL_MSA" ]]; then
        ALN_FILES_FOR_FINAL_MSA+=" "
    fi
    # Add the path of the .aln file to ALN_FILES_FOR_FINAL_MSA
    ALN_FILES_FOR_FINAL_MSA+="$aln_file"
done

# Proceed only if ALN_FILES_FOR_FINAL_MSA is not empty
if [[ -n "$ALN_FILES_FOR_FINAL_MSA" ]]; then
    # Run t_coffee on all .aln files using the same guide tree
    # The output is directed to the final_msa directory
    cd "$FINAL_MSA_DIR"
    t_coffee -aln $ALN_FILES_FOR_FINAL_MSA -matrix=blosum30mt -output fasta_aln -outfile=final_msa.aln -usetree=$TREE_PATH

    cd -  # Go back to the previous directory
fi

echo "Final MSA generated in $FINAL_MSA_DIR."





