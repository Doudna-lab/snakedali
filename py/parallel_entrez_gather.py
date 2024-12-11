# == Native Modules
import datetime
import re
import time
import copy
import os
import sys
import http.client
import socket
import logging
# == Installed Modules
import urllib.request
import urllib.error
import pandas as pd
from mpi4py import MPI
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
# == Project Modules


def divide_work(N, n_rank, world_size):
	# determine the workload of each rank
	workloads = [N // world_size for i in range(world_size)]
	for i in range(N % world_size):
		workloads[i] += 1
	my_start = 0
	for i in range(n_rank):
		my_start += workloads[i]
	my_end = my_start + workloads[n_rank]
	return my_start, my_end


def recover_process(recover_buf, mode):
	combined_content = []
	for item in recover_buf:
		for value in item.values():
			combined_content.extend(value)
	if mode == 'list':
		combined_list = []
		for item in combined_content:
			combined_list.extend(item)
		return combined_list
	if mode == 'dict':
		combined_dict = {}
		for item in combined_content:
			combined_dict.update(item)
		return combined_dict


def ncbi_fetch(fasta_recs_dict, coord_dict, ncbi_db, file_format, rank):
	Entrez.email = "thedoudnalab@gmail.com"
	record_list = []
	report_threshold = 100
	threshold_step = 100
	count_processed = 0

	for acc in fasta_recs_dict:
		if count_processed == report_threshold:
			logging.debug(f"*** Processed {count_processed} entries on rank {rank}***")
			report_threshold += threshold_step

		def efetch_process():
			handle = Entrez.efetch(db=ncbi_db, id=f"{acc}", rettype=file_format, retmode="text")
			return SeqIO.read(handle, file_format)

		record = handle_api_errors(efetch_process)

		# Synchronize accessions orthography
		try:
			sync_acc = re.search(acc, record.id)
			record.id = sync_acc.group()
			print(sync_acc, record.id)
			# Put it back on the record
			record_list.append(record)
		except (AttributeError, NameError):
			pass

		count_processed += 1

	return record_list


def ipg_routine(db, hit_uid):
	link_accession = None
	link_coords = None

	def ipg_handler():
		nucleotide_accession = False
		coords = False
		# Fetch the IPG table for the given protein ID
		handle = Entrez.efetch(db=f"{db}", id=f"{hit_uid}", rettype="ipg", retmode="text")
		ipg_data = handle.read().decode('utf-8')  # Decode the bytes data
		handle.close()
		# Parse the IPG table
		lines = ipg_data.split('\n')
		for line in lines:
			columns = line.split('\t')
			if len(columns) > 6 and columns[6] == hit_uid:
				try:
					nucleotide_accession = columns[2]
					start = int(columns[3])
					stop = int(columns[4])
					coords = [start, stop]
					return nucleotide_accession, coords
				except ValueError:
					continue
		return nucleotide_accession, coords

	try:
		link_accession, link_coords = handle_api_errors(ipg_handler)
	except TypeError:
		pass

	return link_accession, link_coords


def handle_api_errors(func, retries=20, delay=1):
	"""Helper function to handle retries for API calls."""
	for _ in range(retries):
		try:
			return func()
		except urllib.error.HTTPError as e:
			if e.code == 429:
				print(f"Received HTTP 429 error. Retrying in {delay} seconds...")
				time.sleep(delay)
			elif e.code == 500:
				print(f"Received HTTP 500 error. Retrying in {delay} seconds...")
				time.sleep(delay)
			elif e.code == 502:
				print(f"Received HTTP 502 error. Retrying in {delay} seconds...")
				time.sleep(delay)
			else:
				print(f"Unhandled HTTP error: {e.code}. Retrying in {delay} seconds")
				time.sleep(delay)
		except (urllib.error.URLError, RuntimeError, http.client.IncompleteRead, socket.timeout, http.client.HTTPException) as e:
			print(f"API error occurred: {e}. Retrying in {delay} seconds...")
			time.sleep(delay)
	print("Max retries reached. Skipping.")
	return None


def ptn_to_nuc(id_list):
	nucleotide_uid_list = []
	source2target = {}
	nuc_coords_dict = {}
	not_found_list = []

	for seq_id in id_list:
		# Standardize protein identifiers to NCBI UIDs through ESearch
		def search_protein():
			with Entrez.esearch(db="protein", rettype='ipg', term=f"{seq_id}", idtype="acc") as handle:
				return Entrez.read(handle)

		search_record = handle_api_errors(search_protein)
		if not search_record or not search_record.get("IdList"):
			logging.debug(f"Entrez.esearch could not find results for {seq_id} on NCBI")
			print(f"Entrez.esearch could not find results for {seq_id} on NCBI")
			not_found_list.append(seq_id)
			continue

		uid = search_record["IdList"][0]

		# Grab Nuccore UIDs from the source database
		loop_nuc_gi, nuc_coords = ipg_routine("protein", uid)
		if loop_nuc_gi is None:
			not_found_list.append(uid)
			continue
		if loop_nuc_gi:
			source2target[loop_nuc_gi] = seq_id
			nucleotide_uid_list.append(loop_nuc_gi)
			nuc_coords_dict.setdefault(loop_nuc_gi, nuc_coords)
			logging.debug(f"Nuccore GI and protein accession: {loop_nuc_gi} {seq_id}")

	# Outputs nuccore UIDs and the protein->nuc UID links
	logging.debug(f"Total not found hits: {len(not_found_list)}")
	logging.debug(f"Size of found: {len(nucleotide_uid_list)}")
	return nucleotide_uid_list, nuc_coords_dict, source2target, list(set(not_found_list))


def nuc_to_gb(uid_list):
	progress = 0
	threshold_progress_report = 20
	threshold_progress_step = 20
	# Get Genbank records for each Nuccore UID
	gb_records = {}
	fasta_records = {}
	logging.debug(f" --> GATHER REGIONS FOR HITS. PROGRESS")

	def efetch_gb():
		handle = Entrez.efetch(db="nucleotide", id=f"{uid}", rettype="gb", retmode="text")
		return SeqIO.read(handle, "genbank")

	def efecth_fasta():
		handle = Entrez.efetch(db="nucleotide", id=f"{uid}", rettype="fasta", retmode="text")
		return SeqIO.read(handle, "fasta")

	for uid in uid_list:
		progress += 1
		if progress == threshold_progress_report:
			current_time = datetime.datetime.now()
			print(f"[{current_time}] FINISHED {progress} ENTRIES OUT OF {len(uid_list)}")
			threshold_progress_report += threshold_progress_step

		# handle = Entrez.efetch(db="nucleotide", id=f"{uid}", rettype="gb", retmode="text")
		# record = SeqIO.read(handle, "genbank")
		record = handle_api_errors(efetch_gb)
		fasta_record = handle_api_errors(efecth_fasta)

		gb_records.setdefault(uid, record)
		fasta_records.setdefault(uid, fasta_record)
	# Returns a list  of Genbank SeqRecords objects
	return gb_records, fasta_records


def gb_plier(query_to_gb_dict, query_to_fasta_dict, nuc_coords_dict, uid_to_acc, win_size):
	gbk_target = {}
	prot_dict = {}
	start_coord = 0
	end_coord = 0
	highlight_feature = None
	for hit_uid in query_to_gb_dict:
		gbk = query_to_gb_dict[hit_uid]
		try:
			gbk_features = gbk.features
		except AttributeError:
			logging.warning(f"No GenBank features found for {hit_uid}")
			continue
		for seq_feature in gbk_features:
			# Avoid blank feature that may occur in GenBank entries
			try:
				qualifiers = seq_feature.qualifiers
				highlight_feature = copy.deepcopy(seq_feature)
				highlight_feature.type = "highlight"
			except AttributeError:
				continue
			# Restrict search to protein-containing features
			if "protein_id" in qualifiers:
				prot_id = qualifiers["protein_id"][0]
				# Search for the protein-ids of interest
				if re.search(prot_id, uid_to_acc[hit_uid]):
					# Process feature information for future ref
					f_start = seq_feature.location.start.real
					f_end = seq_feature.location.end.real
					f_strand = seq_feature.location.strand
					try:
						f_seq = qualifiers["translation"][0]
					except KeyError:
						f_seq = ''
					f_len = len(f_seq)
					# Set start/end coords using window size
					start_coord = max(int(min([f_start, f_end])) - win_size, 0)
					end_coord = min(int(max([f_start, f_end])) + win_size + 1, len(gbk.seq))
					# Gather protein data for reference
					prep_prot_dict = {
									  "nuccore_acc": gbk.id,
									  # "region_seq": gbk.seq[start:end + 1],
									  "window_start": start_coord,
									  "window_end": end_coord,
									  "feature_start": f_start,
									  "feature_end": f_end,
									  "strand": f_strand,
									  "feature_len": f_len,
									  "mmseqs_hit": prot_id,
									  "ptn_sequence": f_seq
									  }
					prot_dict.setdefault(uid_to_acc[hit_uid], prep_prot_dict)
			else:
				f_start, f_end = nuc_coords_dict[hit_uid]
				# Set start/end coords using window size
				start_coord = max(int(min([f_start, f_end])) - win_size, 0)
				end_coord = min(int(max([f_start, f_end])) + win_size + 1, len(gbk.seq))

		# Process FASTA record
		region_record = query_to_fasta_dict[hit_uid]
		try:
			region_sequence = region_record.seq[start_coord:end_coord]
		except AttributeError:
			logging.debug(f"WARNING: Could not find sequence for {hit_uid}. Skipping")
			continue
		# Create a SeqRecord object with the feature of interest
		gbk_focused = SeqRecord(
			id=gbk.id,
			annotations=gbk.annotations,
			dbxrefs=gbk.dbxrefs,
			seq=region_sequence,
			description=gbk.description
		)
		gbk_focused.features.append(highlight_feature)
		gbk_target.setdefault(f"{gbk.id}_{start_coord}-{end_coord}", gbk_focused)
	return gbk_target, prot_dict


def export_records(rec_list, output_path, file_format,):
	with open(f"{output_path}", "w") as handle:
		SeqIO.write(rec_list, handle, file_format)
	return


def export_entrez(query_to_gb_dict, parent_path, file_format):
	if not os.path.exists(parent_path):
		os.mkdir(parent_path)
	for hit in query_to_gb_dict:
		record = query_to_gb_dict[hit]
		filename = f"{hit}.{file_format}"
		with open(f"{parent_path}{os.sep}{filename}", "w") as handle:
			SeqIO.write(record, handle, str(file_format))


def main():
	#	Inputs
	mmseqs_search_path = str(sys.argv[1])
	#   Outputs
	retrieval_report = str(sys.argv[2])
	#   Params
	col_names = str(sys.argv[3])
	gbk_dir = str(sys.argv[4])
	fna_region_dir = str(sys.argv[5])
	window_size = int(sys.argv[6])
	id_links_table = str(sys.argv[7])
	log_file = str(sys.argv[8])

	# === get basic information about the MPI communicator
	comm = MPI.COMM_WORLD
	world_size = comm.Get_size()
	my_rank = comm.Get_rank()
	divided_work = []
	process_target_gb = {}
	process_target_hit = {}

	# Entrez authentication
	print("Entrez login")
	Entrez.email = "thedoudnalab@gmail.com"

	# Configure the logging system
	logging.basicConfig(
		level=logging.DEBUG,  # Set the minimum log level (DEBUG logs everything)
		format="%(asctime)s %(message)s",  # Define log format
		handlers=[
			logging.FileHandler(log_file),  # Log to a file
		]
	)

	# NCBI databases to search
	efetch_db = 'nuccore'

	# Import blastout table
	col_names_list = col_names.split(" ")
	blast_df = pd.read_csv(mmseqs_search_path,
						   names=col_names_list,
						   index_col=False).convert_dtypes().infer_objects()

	# Process pd Dataframe and take the unique set of hit IDs
	full_unique_mms_hits = list(set(blast_df[col_names].dropna().tolist()))
	bcast_full_unique_mms_hits = MPI.COMM_WORLD.bcast(full_unique_mms_hits, root=0)

	# === DELIMIT THE LENGTH MAIN TASK TO BE SPLIT
	N = len(bcast_full_unique_mms_hits)

	# === CREATE LIST OF TUPLES CONTAINING START AND END INDICES TO BE PROCESSED BY EACH RANK
	if my_rank == 0:
		# == Divide work across ranks
		for rank in range(world_size):
			divided_work.append(divide_work(N, rank, world_size))
		logging.debug(f"*** Divided workload: {divided_work}")

	# === SCATTER WORKLOAD ACROSS THE RANKS
	v = comm.scatter(divided_work, root=0, )

	# === SLICE MAIN LIST ACCORDING TO THE DIVIDED WORKLOAD
	unique_mms_hits = bcast_full_unique_mms_hits[v[0]:v[1]]

	# === PERFORM THE CORE CODE WITHIN THE LOOP
	logging.debug(f"*** Rank: {my_rank}: {v[0]} -> {v[1]}")

	# Query NCBI to get nuccore UIDs associated with the protein hits using ESearch/ELink
	logging.debug(f"RANK {my_rank} --> LINKING PROTEIN HIT IDS TO NUCCORE ENTRIES")
	nuc_uid_per_query, nucleotide_coords_dict, hit_to_link, hits_not_found = ptn_to_nuc(unique_mms_hits)

	# Get GenBank entries through EFetch
	logging.debug(f"RANK {my_rank} --> RETRIEVE GENBANK OBJECTS")
	gb_seqrec_per_query, fasta_seqrec_per_query = nuc_to_gb(nuc_uid_per_query)

	# Narrow down Genbank files based on hit UIDs
	logging.debug(f"RANK {my_rank} --> COMPILE GENOMIC REGIONS IN GENBANK FILES")
	target_gb_dict, target_hit_dict = gb_plier(gb_seqrec_per_query, fasta_seqrec_per_query, nucleotide_coords_dict, hit_to_link, window_size)
	logging.debug(f" --> GENOMIC REGIONS SUCCESSFULLY COMPILED IN GENBANK FILES")

	# Parse hit sequences to FASTA file
	# hit_seq_record_list = ncbi_fetch(fasta_seqrec_per_query, nucleotide_coords_dict, efetch_db, 'fasta', my_rank)

	# == APPEND THE PROCESSED ITEMS TO LIST WITHIN RANK DICTIONARY
	process_target_gb.setdefault(my_rank, []).append(target_gb_dict)
	process_target_hit.setdefault(my_rank, []).append(target_hit_dict)

	# === RECOVER PROCESSED CONTENTS FROM RANKS
	recvbuf_target_gb = comm.gather(process_target_gb, root=0)
	recvbuf_target_hit = comm.gather(process_target_hit, root=0)

	# logging.debug(f"Recovery target_gb: {recvbuf_target_gb}\n\n")
	# logging.debug(f"Recovery target_hit: {recvbuf_target_hit}\n\n")
	# logging.debug(f"Recovery hit_seq_record: {recvbuf_hit_seq_record}\n\n")

	# === COMBINE CONTENT ON ROOT RANK
	if comm.rank == 0:
		combined_target_gb = recover_process(recvbuf_target_gb, 'dict')
		combined_target_hit = recover_process(recvbuf_target_hit, 'dict')
		# combined_hit_seq_record = recover_process(recvbuf_hit_seq_record, 'list')

		# Generate summary dataframe
		current_time = datetime.datetime.now()
		logging.debug(f" --> GENERATE REPORTS")
		print(f"[{current_time}] --> GENERATE REPORTS")
		df = pd.DataFrame.from_dict(combined_target_hit, orient='index')
		# Pull information from Fseek cluster members and MMSEQS Hits
		# == Import ID links table
		id_links_df = pd.read_csv(id_links_table, sep="\t", header=None)[[0, 1]]
		seed_report_df = df.merge(id_links_df, left_on='mmseqs_hit', right_on=1).drop(columns=1)

		# Export outputs
		current_time = datetime.datetime.now()
		logging.debug(f" --> EXPORT FILES")
		print(f"[{current_time}] --> EXPORT FILES")

		# == Export GenBank Files
		export_entrez(combined_target_gb, gbk_dir, "gb")
		# == Export FASTA nucleotide of target regions + window
		export_entrez(combined_target_gb, fna_region_dir, "fasta")

		seed_report_df.to_csv(retrieval_report, index=False)

		# hit_taxid_seq_record_list = attach_label_to_fasta(hit_seq_record_list, id_to_taxid)
		# export_records(combined_hit_seq_record, output_fasta_region, 'fasta')

		logging.debug(f" --> ALL PROCESSES FINISHED")


if __name__ == "__main__":
	main()
