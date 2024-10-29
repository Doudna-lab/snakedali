# == Native Modules
import re
import sys
import argparse
import os
import pickle
import datetime
import subprocess
from pathlib import Path
from collections import defaultdict
# == Installed Modules
from Bio import SeqIO
import pandas as pd
# from mpi4py import MPI
# # == Project Modules
#
#
# # DEBUG
# ROOT_DIR = "/Users/bellieny/projects/snakedali/dump"
# # === Inputs
# dali_alignment_list = [f"{ROOT_DIR}/3800A.txt"]
# # === Params
# list_of_files_path = ""
# input_prefix = '3800'
# fseek_clusters = f"{ROOT_DIR}/toy_fseek_clustered_afdb_representatives.txt"
# afdb_pkl = f"{ROOT_DIR}/afdb_pkl/sequences.pkl"
# id_converstion_path = f"{ROOT_DIR}/clustered_AFDB_structure_key.txt"
# # === Outputs
# output_directory = '/Users/bellieny/projects/snakedali/dump'
# multi_fasta_out = args.multi_fasta_out
# === Wildcards
# input_prefix = args.input_prefix
#


def format_conversion_dict(cluster_reps_tbl):
	cluster_reps_pre = pd.read_csv(cluster_reps_tbl, sep="\t", header=None, low_memory=False)
	cluster_reps_pre.index = cluster_reps_pre.iloc[:,1]
	cluster_reps_dict = defaultdict(list)
	for row in cluster_reps_pre.itertuples():
		cluster_reps_dict[row.Index].append(row._1)
	return cluster_reps_dict


def ingest_conversion_table(id_converstion_table_path):
	id_converstion_table = pd.read_csv(id_converstion_table_path, sep="\t", low_memory=False, names=['0', '1'])
	converted_id_dict = dict(zip(id_converstion_table.iloc[:, 1], id_converstion_table.iloc[:, 0]))
	return converted_id_dict


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


def remove_duplicates_seqrecords(seq_records):
	seen = set()
	unique_records = []

	for record in seq_records:
		key = (record.id, str(record.seq))  # Use a tuple of (ID, sequence) as the key
		if key not in seen:
			unique_records.append(record)
			seen.add(key)

	return unique_records


def process_dali_file_with_block_info(file_content, input_prefix, output_dir):
	# Processes the entire DALI file, extracting block information and amino acid sequences
	# for each alignment block, and outputs to separate files without the header line
	alignments = file_content.split('\n\n\n')
	multi_id_list = []
	processed_hits = []
	id_zscores = []
	temp_fasta_list_path = os.path.join(output_dir, "temp_fasta_list.txt")

	# Collect the filenames to write to temp_fasta_list.txt
	with open(temp_fasta_list_path, 'w') as temp_fasta_list_file:
		for alignment in alignments:
			zscore = 0
			if alignment.strip():
				# Set internal Variables
				block_info, aa_query, aa_subject = extract_block_info_and_aa_sequences(alignment)
				try:
					subject_id = block_info.split('Sbjct=')[1].split()[0]  # Extract Subject ID
					zscore = block_info.split('Z-score=')[1]
				except IndexError:
					subject_id = 'no_hits'
				aln_filename = "{}_{}.fasta".format(input_prefix, subject_id)  # Create filename from Query and Sbjct IDs

				# Establish a hard cut at DaliLIte's z-score of 6
				if float(zscore) < 6:
					continue
				# Write to a separate file for each block, excluding the header line
				output_path = os.path.join(output_dir, aln_filename)
				with open(output_path, 'w') as output_file:
					output_file.write(">{0}\n{1}\n\n>{2}\n{3}\n".format(input_prefix, aa_query, subject_id, aa_subject))
					processed_hits.append(aln_filename)
					# Write each output filename to temp_fasta_list.txt
					temp_fasta_list_file.write(f"{output_path}\n")

				# Set up aggregated multi-fasta
				# aa_subject_nogaps = re.sub('-', '', aa_subject)
				# seq_record = SeqIO.SeqRecord(id=subject_id, seq=aa_subject_nogaps)
				multi_id_list.append(subject_id)
				id_zscores.append((subject_id, float(zscore)))

	return processed_hits, multi_id_list, id_zscores


def extract_fasta_ids(input_prefix, fseek_clusters, dali_alignment_list, id_conversion_dict, output_directory):
	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> GATHERING ALL FASTA IDs ASSOCIATED WITH QUERY {input_prefix}")
	list_of_processed_files = []
	list_of_multi_fasta_out = []
	cluster_reps_dict = format_conversion_dict(fseek_clusters)
	# === Generate a list of FASTA IDs associated with the query for further processing
	for index in range(len(dali_alignment_list)):
		with open(dali_alignment_list[index], 'r') as file:
			file_content = file.read()
			(duo_fasta_path, multi_fasta_id_list, zscore_list) = process_dali_file_with_block_info(file_content,
																					  input_prefix, output_directory)

			# This will export a list of FASTA files for downstream processing by T-COFFEE
			list_of_processed_files.extend(duo_fasta_path)
			current_time = datetime.datetime.now()
			print(f"{current_time}: Start FASTA scanning")
			# This will prepare an export with a list of FASTA files for downstream processing by MMSEQS2
			converted_id_list = []
			for fasta_id in multi_fasta_id_list:
				if len(fasta_id) != 5:
					print(f"*** WARNING: The FASTA ID in the AFDB database has an invalid label: {fasta_id}")
					continue
				if fasta_id == 'no_hits':
					continue
				for item in cluster_reps_dict[id_conversion_dict[re.sub(r"A$", '', fasta_id)]]:
					converted_id_list.append(item)
			current_time = datetime.datetime.now()
			print(f"{current_time}: FASTA scanning finished")
			list_of_multi_fasta_out.extend(converted_id_list)
	return list_of_multi_fasta_out, list_of_processed_files, zscore_list


def gather_afdb_sequences(input_prefix, list_of_multi_fasta_out, afdb_pkl_dict):
	current_time = datetime.datetime.now()
	print(
		f"[{current_time}] --> GATHERING PROTEIN SEQUENCES OF {len(list_of_multi_fasta_out)} FASTA IDs FROM QUERY {input_prefix}")
	seqrecord_list = []
	track_process = 0
	threshold_increment = 100
	threshold_report_process = 100
	for fasta_id in list_of_multi_fasta_out:
		try:
			seqrecord_list.append(afdb_pkl_dict[fasta_id])
		except KeyError:
			continue
		track_process += 1
		if track_process == threshold_report_process:
			current_time = datetime.datetime.now()
			print(f"[{current_time}] --> PROCESSED {track_process} FASTA IDs FROM: {input_prefix}")
			threshold_report_process += threshold_increment
	return seqrecord_list

# def process_indices(local_indices,
# 					items_list,
# 					input_prefix,
# 					conversion_dict_id,
# 					reference_fasta,
# 					cluster_reps_tbl,
# 					output_directory):
	# list_of_processed_files_per_rank = []
	# list_of_multi_fasta_out = []
	# cluster_reps_dict = format_conversion_dict(cluster_reps_tbl)
	# for index in range(len(local_indices)):
	# 	with open(items_list[index], 'r') as file:
	# 		file_content = file.read()
	# 		(duo_fasta_path, multi_fasta_id_list) = process_dali_file_with_block_info(file_content,
	# 																				  input_prefix, output_directory)
	# 		# This will export a list of FASTA files for downstream processing by T-COFFEE
	# 		list_of_processed_files_per_rank.extend(duo_fasta_path)
	# 		# This will prepare an export with a list of FASTA files for downstream processing by MMSEQS2
	# 		converted_id_list = [item for fasta_id in multi_fasta_id_list for item in cluster_reps_dict[conversion_dict_id[fasta_id.strip('A')]]]
	# 		list_of_multi_fasta_out.extend(converted_id_list)
	#
	# # Gather FASTA records of cluster members from
	# # 	the representatives present at the DALI alignment
	# seqrecord_list = []
	# for fasta_id in list_of_multi_fasta_out:
	# 	pattern = re.compile(fasta_id)
	# 	# print(f"Compiling Pattern {fasta_id}")
	# 	# Parse the file iteratively
	# 	for record in SeqIO.parse(reference_fasta, "fasta"):  # Adjust the format if not GenBank
	# 		# Search using regex in the desired field of the record
	# 		clean_id = re.sub(r"\w+:\w+\-(\w+)\-\w+", r"\1", record.id)
	# 		# print(f"Searching {clean_id} == Original id: {record.id}")
	# 		if not pattern.search(str(clean_id)):
	# 			continue
	# 		seqrecord_list.append(record)
	# 		break
	#
	# return list_of_processed_files_per_rank, seqrecord_list


def main():
	# Parse the command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--input_prefix", dest="input_prefix", required=True)
	parser.add_argument("--afdb_pkl", dest="afdb_pkl", required=True)
	parser.add_argument("--fseek_clusters", dest="fseek_clusters", required=True)
	parser.add_argument("--id_converstion", dest="id_converstion", required=True)
	parser.add_argument("--output_dir", dest="output_dir", required=True)
	parser.add_argument("--files_list", nargs='+', dest="files_list", required=True)
	parser.add_argument("--manifest_out", dest="manifest_out", required=True)
	parser.add_argument("--multi_fasta_out", dest="multi_fasta_out", required=True)
	args = parser.parse_args()

	# Snakemake I/O
	# === Inputs
	dali_alignment_list = list(args.files_list)
	# === Params
	list_of_files_path = args.manifest_out
	fseek_clusters = args.fseek_clusters
	afdb_pkl = args.afdb_pkl
	id_converstion_path = args.id_converstion
	# === Outputs
	# output_directory = str(snakemake.output)
	# output_directory = str(sys.argv[1])
	output_directory = args.output_dir
	multi_fasta_out = args.multi_fasta_out
	# === Wildcards
	# input_prefix = str(snakemake.wildcards.query_name)
	# input_prefix = str(sys.argv[3])
	input_prefix = args.input_prefix

	print(f"Input alignment captured in list: {dali_alignment_list}")

	# === Ingest conversion table
	id_conversion_dict = ingest_conversion_table(id_converstion_path)

	# === Ingest AFDB pickled Dictionary
	seqrecord_list = []
	combined_processed_files = []
	for pkl_dict_index in range(10):
		# === Account for the 10 chunks of serialized dictionaries containing AFDB FASTA sequences
		filename = Path(f"{afdb_pkl}{pkl_dict_index + 1}")
		if not filename.exists():
			print(f"Couldn't find AFDB PICKLED DICTIONARY in {filename}")
			exit(0)

		current_time = datetime.datetime.now()
		print(f"[{current_time}] --> IMPORTING AFDB PICKLED DICTIONARY: {filename}")
		with open(filename, 'rb') as afdb_file:
			afdb_pkl_dict = pickle.load(afdb_file)
			current_time = datetime.datetime.now()
			print(f"[{current_time}] --> GATHER ASSOCIATED FASTA IDS FOLDSEEK CLUSTERS")
			list_of_multi_fasta_out, list_of_processed_files, zscore_list = extract_fasta_ids(input_prefix,
																				 fseek_clusters,
																				 dali_alignment_list,
																				 id_conversion_dict,
																				 output_directory)
			current_time = datetime.datetime.now()
			print(f"[{current_time}] --> GATHER FASTA SEQUENCES FROM SERIALIZED AFDB FASTA RECORDS")
			loop_seqrecord_list = gather_afdb_sequences(input_prefix, list_of_multi_fasta_out, afdb_pkl_dict)
		seqrecord_list.extend(loop_seqrecord_list)
		combined_processed_files.extend(list_of_processed_files)

	# # Initialize MPI
	# comm = MPI.COMM_WORLD
	# rank = comm.Get_rank()
	# # size = comm.Get_size()
	#
	# # Scatter the indices to processes
	# if rank == 0:
	# 	divided_indices = [dali_alignment_list[start::comm.size] for start in range(comm.size)]
	# 	print(f"Comm size: {comm.size}, Data size: {len(divided_indices)}")
	# 	# divided_indices = distribute_work(dali_alignment_list, size)
	# else:
	# 	divided_indices = None
	#
	# try:
	# 	local_indices = comm.scatter(divided_indices, root=0)
	# except ValueError:
	# 	print("An error occurred when distributing the work across the cores. "
	# 		  "Make sure to scale the number of cores according to the number of taks.")
	# 	exit(0)


	##### ====== TURNED INTO FUNCTION: extract_fasta_ids  ======
	# current_time = datetime.datetime.now()
	# print(f"[{current_time}] --> GATHERING ALL FASTA IDs ASSOCIATED WITH QUERY {input_prefix}")
	# list_of_processed_files = []
	# list_of_multi_fasta_out = []
	# cluster_reps_dict = format_conversion_dict(fseek_clusters)
	# # === Generate a list of FASTA IDs associated with the query for further processing
	# for index in range(len(dali_alignment_list)):
	# 	with open(dali_alignment_list[index], 'r') as file:
	# 		file_content = file.read()
	# 		(duo_fasta_path, multi_fasta_id_list) = process_dali_file_with_block_info(file_content,
	# 																				  input_prefix, output_directory)
	# 		current_time = datetime.datetime.now()
	# 		print(f"{current_time}: Gathered alignment {index} of {len(dali_alignment_list)} "
	# 			  f"\n list of multi-fasta with {len(multi_fasta_id_list)} items")
	# 		# This will export a list of FASTA files for downstream processing by T-COFFEE
	# 		list_of_processed_files.extend(duo_fasta_path)
	# 		print(f"{current_time}: Start FASTA scanning at ")
	# 		# This will prepare an export with a list of FASTA files for downstream processing by MMSEQS2
	# 		converted_id_list = []
	# 		for fasta_id in multi_fasta_id_list:
	# 			if fasta_id == 'no_hits':
	# 				continue
	# 			for item in cluster_reps_dict[id_conversion_dict[fasta_id.strip('A')]]:
	# 				converted_id_list.append(item)
	# 		list_of_multi_fasta_out.extend(converted_id_list)
	##### ====== extract_fasta_ids ======


	# Gather FASTA records of cluster members from
	# 	the representatives present at the DALI alignment
	##### ====== TURNED INTO FUNCTION: gather_afdb_sequences  ======
	# current_time = datetime.datetime.now()
	# print(f"[{current_time}] --> GATHERING PROTEIN SEQUENCES OF {len(list_of_multi_fasta_out)} FASTA IDs from: {input_prefix}")
	# seqrecord_list = []
	# track_process = 0
	# threshold_increment = 100
	# threshold_report_process = 100
	# for fasta_id in list_of_multi_fasta_out:
	# 	seqrecord_list.append(afdb_pkl_dict[fasta_id])
	# 	track_process += 1
	# 	if track_process == threshold_report_process:
	# 		current_time = datetime.datetime.now()
	# 		print(f"[{current_time}] --> PROCESSED {track_process} FASTA IDs FROM: {input_prefix}")
	# 		threshold_report_process += threshold_increment
	##### ======  gather_afdb_sequences  ======
		# pattern = re.compile(fasta_id)
		# # print(f"Compiling Pattern {fasta_id}")
		# # Parse the file iteratively
		# for record in SeqIO.parse(afdb_pkl, "fasta"):  # Adjust the format if not GenBank
		# 	# Search using regex in the desired field of the record
		# 	clean_id = re.sub(r"\w+:\w+\-(\w+)\-\w+", r"\1", record.id)
		# 	# print(f"Searching {clean_id} == Original id: {record.id}")
		# 	if not pattern.search(str(clean_id)):
		# 		continue
		# 	seqrecord_list.append(afdb_pkl_dict[fasta_id])
		# 	break
	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> ALL FASTA IDs PROCESSED. "
		  f"EXPORTING FASTA FILE TO {multi_fasta_out}")
	multifasta_unique = remove_duplicates_seqrecords(seqrecord_list)
	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> EXPORTING FASTA FILE TO {multi_fasta_out}")
	SeqIO.write(multifasta_unique, multi_fasta_out, format='fasta')

	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> EXPORTING LIST OF PROCESSED FILES TO {list_of_files_path}")
	combined_processed_files_str = "\n".join(set(combined_processed_files))
	with open(list_of_files_path, "w") as f:
		f.write(combined_processed_files_str)
	# return list_of_processed_files_per_rank, seqrecord_list

	# (processed_files_per_rank, multi_fasta_per_rank) = process_indices(local_indices,
	# 																   dali_alignment_list,
	# 																   input_prefix,
	# 																   id_conversion_dict,
	# 																   afdb_pkl,
	# 																   fseek_clusters,
	# 																   output_directory)

	# # Gather results from all processes
	# all_processed_files = comm.gather(processed_files_per_rank, root=0)
	# all_multi_fasta = comm.gather(multi_fasta_per_rank, root=0)
	#
	# if rank == 0:
	# 	# Combine results on the root process
	# 	print(f"Combining results from {all_multi_fasta}")
	# 	combined_multifasta = []
	# 	for d in all_multi_fasta:
	# 		print(f"Processing {d}")
	# 		combined_multifasta.extend(d)
	# 	# Try to export directly from comm.gather
	# 	if combined_multifasta is not None:
	# 		combined_multifasta_unique = remove_duplicates_seqrecords(combined_multifasta)
	# 		SeqIO.write(combined_multifasta_unique, multi_fasta_out, format='fasta')
	#
	# 	# Combine results on the root process
	# 	print(f"Combining results from {all_processed_files}")
	# 	combined_processed_files = []
	# 	for d in all_processed_files:
	# 		print(f"Processing {d}")
	# 		combined_processed_files.extend(d)
	#
	# 	print(f"EXPORTING LIST OF PROCESSED FILES TO {list_of_files_path}")
	# 	combined_processed_files_str = "\n".join(set(combined_processed_files))
	# 	with open(list_of_files_path, "w") as f:
	# 		f.write(combined_processed_files_str)


if __name__ == "__main__":
	main()
