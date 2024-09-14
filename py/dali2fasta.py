# == Native Modules
import re
import sys
import argparse
import os
import datetime
import subprocess
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


def process_dali_file_with_block_info(file_content, tcoffee_params, tcoffee_bin, input_prefix, output_dir):
	# Processes the entire DALI file, extracting block information and amino acid sequences
	# for each alignment block, and outputs to separate files without the header line
	alignments = file_content.split('\n\n\n')
	processed_hits = []
	temp_fasta_list_path = os.path.join(output_dir, "temp_fasta_list.txt")

	# Collect the filenames to write to temp_fasta_list.txt
	with open(temp_fasta_list_path, 'w') as temp_fasta_list_file:
		for alignment in alignments:
			if alignment.strip():
				# Set internal Variables
				block_info, aa_query, aa_subject = extract_block_info_and_aa_sequences(alignment)
				try:
					subject_id = block_info.split('Sbjct=')[1].split()[0]  # Extract Subject ID
				except IndexError:
					subject_id = 'no_hits'
				aln_filename = "{}_{}.fasta".format(input_prefix, subject_id)  # Create filename from Query and Sbjct IDs

				# Write to a separate file for each block, excluding the header line
				output_path = os.path.join(output_dir, aln_filename)
				with open(output_path, 'w') as output_file:
					output_file.write(">{0}\n{1}\n\n>{2}\n{3}\n".format(input_prefix, aa_query, subject_id, aa_subject))
					processed_hits.append(aln_filename)
					# Write each output filename to temp_fasta_list.txt
					temp_fasta_list_file.write(f"{output_path}\n")

				# Format the t-coffee command
				subprocess_cmd = re.sub(r"\|input\|", output_path, tcoffee_params)
				formatted_cmd = f"{tcoffee_bin} {subprocess_cmd}"
				# Execute the command using subprocess
				try:
					print(f"T-COFFEE COMMAND: {formatted_cmd}")
					subprocess.run(formatted_cmd, shell=True, check=True)
				except subprocess.CalledProcessError as e:
					print(f"Error running t_coffee: {e}")
					exit(0)

	return processed_hits


def process_indices(local_indices, items_list, tcoffee_params, tcoffee_bin, input_prefix, output_directory):
	list_of_processed_files_per_rank = []
	for index in range(len(local_indices)):
		with open(items_list[index], 'r') as file:
			file_content = file.read()
			duo_fasta_path = process_dali_file_with_block_info(file_content, tcoffee_params, tcoffee_bin,
															   input_prefix, output_directory)
			list_of_processed_files_per_rank.extend(duo_fasta_path)
	return list_of_processed_files_per_rank


def main():
	# Parse the command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--input_prefix", dest="input_prefix", required=True)
	parser.add_argument("--output_dir", dest="output_dir", required=True)
	parser.add_argument("--files_list", nargs='+', dest="files_list", required=True)
	parser.add_argument("--tcoffee_params", dest="tcoffee_params", required=True)
	parser.add_argument("--tcoffee_bin", dest="tcoffee_bin", required=True)
	parser.add_argument("--manifest_out", dest="manifest_out", required=True)

	args = parser.parse_args()

	# Snakemake I/O
	# === Inputs
	# dali_alignment_list = list(snakemake.input.alignment_list)
	# dali_alignment_list = list(sys.argv[5:])
	dali_alignment_list = list(args.files_list)
	# === Params
	# list_of_files_path = str(sys.argv[2])
	list_of_files_path = args.manifest_out
	# tcoffee_params = str(sys.argv[4])
	tcoffee_params = args.tcoffee_params
	tcoffee_bin = args.tcoffee_bin

	# === Outputs
	# output_directory = str(snakemake.output)
	# output_directory = str(sys.argv[1])
	output_directory = args.output_dir
	# === Wildcards
	# input_prefix = str(snakemake.wildcards.query_name)
	# input_prefix = str(sys.argv[3])
	input_prefix = args.input_prefix

	print(f"Input alignment captured in list: {dali_alignment_list}")

	# Initialize MPI
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	# size = comm.Get_size()

	# Scatter the indices to processes
	if rank == 0:
		divided_indices = [dali_alignment_list[start::comm.size] for start in range(comm.size)]
		print(f"Comm size: {comm.size}, Data size: {len(divided_indices)}")
		# divided_indices = distribute_work(dali_alignment_list, size)
	else:
		divided_indices = None

	try:
		local_indices = comm.scatter(divided_indices, root=0)
	except ValueError:
		print("An error occurred when distributing the work across the cores. "
			  "Make sure to scale the number of cores according to the number of taks.")
		exit(0)

	current_time = datetime.datetime.now()
	print(f"PROCESSING QUERY {input_prefix} on rank {rank} | Starting at {current_time}")
	processed_files_per_rank = process_indices(local_indices,
											   dali_alignment_list,
											   tcoffee_params,
											   tcoffee_bin,
											   input_prefix,
											   output_directory)

	# Gather results from all processes
	all_processed_files = comm.gather(processed_files_per_rank, root=0)

	if rank == 0:
		# Combine results on the root process
		print(f"Combining results from {all_processed_files}")
		combined_processed_files = []
		for d in all_processed_files:
			print(f"Processing {d}")
			combined_processed_files.extend(d)

		print(f"EXPORTING LIST OF PROCESSED FILES TO {list_of_files_path}")
		combined_processed_files_str = "\n".join(set(combined_processed_files))
		with open(list_of_files_path, "w") as f:
			f.write(combined_processed_files_str)


if __name__ == "__main__":
	main()
