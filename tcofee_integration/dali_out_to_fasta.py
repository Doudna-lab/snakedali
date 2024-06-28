import argparse
import os

def extract_block_info_and_aa_sequences(alignment_block):
    """Extracts block information and amino acid sequences for Query and Subject from an alignment block."""
    lines = alignment_block.split('\n')
    block_info = ""
    aa_query = []
    aa_subject = []

    for line in lines:
        if line.startswith('No '):
            block_info = line
        elif line.startswith('Query '):
            sequence_line = line.split()[1]  # Extract the sequence part after 'Query'
            aa_query.append(sequence_line)
        elif line.startswith('Sbjct '):
            sequence_line = line.split()[1]  # Extract the sequence part after 'Sbjct'
            aa_subject.append(sequence_line)

    return block_info, ''.join(aa_query), ''.join(aa_subject)

def process_dali_file_with_block_info(file_content, output_dir):
    """Processes the entire DALI file, extracting block information and amino acid sequences for each alignment block, and outputs to separate files without the header line."""
    alignments = file_content.split('\n\n\n')

    for alignment in alignments:
        if alignment.strip():
            block_info, aa_query, aa_subject = extract_block_info_and_aa_sequences(alignment)
            query_id = block_info.split('Query=')[1].split()[0]  # Extract Query ID
            subject_id = block_info.split('Sbjct=')[1].split()[0]  # Extract Subject ID
            aln_filename = "{}_{}.fasta".format(query_id, subject_id)  # Create filename from Query and Sbjct IDs

            # Write to a separate file for each block, excluding the header line
            output_path = os.path.join(output_dir, aln_filename)
            with open(output_path, 'w') as output_file:
                output_file.write(">{0}\n{1}\n\n>{2}\n{3}\n".format(query_id, aa_query, subject_id, aa_subject))

def main():
    parser = argparse.ArgumentParser(description='Process DALI alignment file and output each block\'s alignment information into separate files, excluding the header line.')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the DALI alignment file')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory path to save the output files')

    args = parser.parse_args()

    # Ensure output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.input_file, 'r') as file:
        file_content = file.read()

    process_dali_file_with_block_info(file_content, args.output_dir)

if __name__ == "__main__":
    main()

