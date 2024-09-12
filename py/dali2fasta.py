# == Native Modules
import sys
import os
import datetime
# == Installed Modules
from mpi4py import MPI
# == Project Modules


def distribute_work(items, num_workers):
	avg = len(items) / float(num_workers)
	out = []
	last = 0.0

	while last < (len(items) - 1):
		out.append(items[int(last):int(last + avg)])
		last += avg

	return out


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


def process_dali_file_with_block_info(file_content, query_id, output_dir):
	"""Processes the entire DALI file, extracting block information and amino acid sequences for each alignment block, and outputs to separate files without the header line."""
	alignments = file_content.split('\n\n\n')

	for alignment in alignments:
		if alignment.strip():
			block_info, aa_query, aa_subject = extract_block_info_and_aa_sequences(alignment)
			# query_id = block_info.split('Query=')[1].split()[0]  # Extract Query ID
			subject_id = block_info.split('Sbjct=')[1].split()[0]  # Extract Subject ID
			aln_filename = "{}_{}.fasta".format(query_id, subject_id)  # Create filename from Query and Sbjct IDs

			# Write to a separate file for each block, excluding the header line
			output_path = os.path.join(output_dir, aln_filename)
			with open(output_path, 'w') as output_file:
				output_file.write(">{0}\n{1}\n\n>{2}\n{3}\n".format(query_id, aa_query, subject_id, aa_subject))
			return aln_filename


def process_indices(local_indices, items_list, input_prefix, output_directory):
	list_of_processed_files_per_rank = []
	for index in local_indices:
		with open(items_list[index], 'r') as file:
			file_content = file.read()
			duo_fasta_path = process_dali_file_with_block_info(file_content, input_prefix, output_directory)
			list_of_processed_files_per_rank.append()
	return list_of_processed_files_per_rank


def main():
	# Snakemake I/O
	# === Inputs
	# dali_alignment_list = list(snakemake.input.alignment_list)
	dali_alignment_list = list(sys.argv[1])
	# === Params
	list_of_processed_files = str(sys.argv[3])
	# === Outputs
	# output_directory = str(snakemake.output)
	output_directory = str(sys.argv[2])
	# === Wildcards
	# input_prefix = str(snakemake.wildcards.query_name)
	input_prefix = str(sys.argv[4])

	#Setup output filename
	list_of_files_path = f"{list_of_processed_files}_{input_prefix}.txt"

	# Initialize MPI
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()

	# Scatter the indices to processes
	if rank == 0:
		divided_indices = distribute_work(dali_alignment_list, size)
	else:
		divided_indices = None

	local_indices = comm.scatter(divided_indices, root=0)

	current_time = datetime.datetime.now()
	print(f"PROCESSING QUERY {input_prefix} on rank {rank} | Starting at {current_time}")
	processed_files_per_rank = process_indices(local_indices, dali_alignment_list, input_prefix, output_directory)

	# Gather results from all processes
	all_processed_files = comm.gather(processed_files_per_rank, root=0)

	if rank == 0:
		# Combine results on the root process
		combined_processed_files = []
		for d in all_processed_files:
			combined_processed_files.append(d)

		print(f"EXPORTING LIST OF PROCESSED FILES TO {list_of_files_path}")
		combined_processed_files_str = " ".join(combined_processed_files)
		with open(list_of_files_path, "w") as f:
			f.write(combined_processed_files_str)


if __name__ == "__main__":
	main()
