# Native modules
import datetime
import re
import time
import copy
import os
import sys
# Installed modules
import urllib.request
import urllib.error
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from bioservices import UniProt


# # DEBUG INPUTS
# mmseqs_search_path = "/Users/bellieny/projects/snakedali/dump/toy_3800_vs_mmseqsDB_result-mms_hits.tsv"
#
# col_names = 'hit_id'
# window_size = 2000
# parent_dir = "/Users/bellieny/projects/snakedali/dump/gbk"
# config_path = "/groups/doudna/team_resources/toolbox/phyrec_screening/config/phyrec_processing.yaml"
# # Load config files
# import yaml
# with open(config_path, "r") as f:
# 	config = yaml.load(f, Loader=yaml.FullLoader)
# blast_col_names = config["blast_custom_cols"]


def ncbi_fetch(acc_list, ncbi_db, file_format):
	Entrez.email = "thedoudnalab@gmail.com"
	record_list = []
	max_retries = 50
	for acc in acc_list:

		# Attempt the search with retries
		for _ in range(max_retries):
			try:
				handle = Entrez.efetch(db=ncbi_db, id=f"{acc}", rettype=file_format, retmode="text")
				record = SeqIO.read(handle, file_format)
				handle.close()
				break  # Break the loop if successful
			except urllib.error.HTTPError as e:
				if e.code == 429:  # HTTP 429: Too Many Requests
					print(f"Received HTTP 429 error. Retrying in 10 seconds...")
					time.sleep(10)
				else:
					continue  # Re-raise other HTTP errors
			except urllib.error.URLError:
				continue

		# Synchronize accessions orthography
		try:
			sync_acc = re.search(acc, record.id)
			record.id = sync_acc.group()
			print(sync_acc, record.id)
			# Put it back on the record
			record_list.append(record)
		except (AttributeError, NameError):
			pass
	return record_list


# def ukb2ncbi(uid):
# 	u = UniProt(verbose=False)
# 	# gbk_id = u.mapping("UniProtKB_AC-ID", "EMBL-GenBank-DDBJ", query=uid, polling_interval_seconds=3, max_waiting_time=100)["results"][0]["to"]
# 	try:
# 		prot_id = u.mapping("UniProtKB_AC-ID", "EMBL-GenBank-DDBJ_CDS", query=uid, polling_interval_seconds=3, max_waiting_time=100)["results"][0]["to"]
# 	except (TypeError, IndexError):
# 		prot_id = ''
# 	return prot_id


# def elink_routine(db, hit_uid):
# 	dup_check = []
# 	not_found = ""
# 	linked = ""
# 	link_record = ""
# 	server_attempts = 0
# 	try:
# 		handle = Entrez.elink(dbfrom="protein", db=db, id=f"{hit_uid}")
# 	except urllib.error.HTTPError as err:
# 		if err.code == 500:
# 			print(f'An internal server error occurred while handling the accession {hit_uid}')
# 			not_found = hit_uid
# 			return linked, hit_uid, not_found
# 	try:
# 		link_record = Entrez.read(handle)
# 	except RuntimeError:
# 		not_found = hit_uid
# 	if link_record:
# 		try:
# 			linked = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
# 			if linked not in dup_check:
# 				dup_check.append(linked)
# 		except (IndexError, KeyError):
# 			not_found = hit_uid
# 	handle.close()
# 	return linked, hit_uid, not_found


def elink_routine(db, dbfrom, hit_uid):
	dup_check = []
	not_found = ""
	linked = ""
	link_record = ""
	max_retries = 20
	handle = ""
	for _ in range(max_retries):
		try:
			handle = Entrez.elink(dbfrom=f"{dbfrom}", db=db, id=f"{hit_uid}")
		except urllib.error.HTTPError as err:
			if err.code == 500:
				print(f'An internal server error occurred while handling the accession {hit_uid}')
				not_found = hit_uid
				return linked, hit_uid, not_found
		except urllib.error.URLError:
			print(f'An internal server error occurred while handling the accession {hit_uid}')
			print(f"Received HTTP 403 error. Retrying in 5 seconds...")
			time.sleep(5)
			continue
	try:
		link_record = Entrez.read(handle)
	except (RuntimeError, AttributeError, ValueError, FileNotFoundError):
		not_found = hit_uid
	if link_record:
		try:
			linked = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
			if linked not in dup_check:
				dup_check.append(linked)
		except (IndexError, KeyError):
			not_found = hit_uid
	try:
		handle.close()
	except AttributeError:
		pass
	return linked, link_record, not_found

# def ptn_to_nuc(hit_list: list, db_list):
# 	progress = 0
# 	threshold_progress_report = 20
# 	threshold_progress_step = 20
# 	nucleotide_uid_list = []
# 	source2target = {}
# 	not_found_list = []
# 	dup_check = []
# 	current_time = datetime.datetime.now()
# 	print(f"[{current_time}] --> GATHER NUCLEOTIDE ID FOR HITS. PROGRESS:")
# 	# Loop through MMSEQShits
# 	for hit in hit_list:
# 		progress += 1
# 		if progress == threshold_progress_report:
# 			current_time = datetime.datetime.now()
# 			print(f"[{current_time}] FINISHED {progress} ENTRIES OUT OF {len(hit_list)}")
# 			threshold_progress_report += threshold_progress_step
# 		# Standardize protein identifiers to NCBI UIDs through ESearch
# 		handle = Entrez.esearch(db="protein", term=f"{hit}", idtype="acc")
# 		current_time = datetime.datetime.now()
# 		print(f"[{current_time}] --> SUCCESS IN PERFOMING ESEARCH ON HIT {hit} -- "
# 			  f"NOW READ HANDLE")
# 		search_record = Entrez.read(handle)
# 		try:
# 			uid = search_record['IdList'][0]
# 		except IndexError:
# 			continue
# 		handle.close()
#
# 		# Loop through databases (found in config) and grab Nuccore UIDs
# 		for db_name in db_list:
# 			if uid in set(dup_check):
# 				continue
# 			loop_nuc_gi, loop_nuc_acc, not_found_hit = elink_routine(db_name, uid)
# 			if not_found_hit:
# 				not_found_list.append(not_found_hit)
# 				continue
# 			if loop_nuc_gi:
# 				dup_check.append(uid)
# 				source2target.setdefault(loop_nuc_gi, (loop_nuc_acc, hit))
# 				nucleotide_uid_list.append(loop_nuc_gi)
#
# 	# Ouputs nuccore uids and the ptn->nuc uid links
# 	return nucleotide_uid_list, source2target, list(set(not_found_list))


def ptn_to_nuc(id_list, db_name_list):
	progress = 0
	nucleotide_uid_list = []
	source2target = {}
	not_found_list = []

	for seq_id in id_list:
		progress += 1
		max_retries = 20
		search_record = {}
		# Standardize protein identifiers to NCBI UIDs through ESearch
		# Introduce a delay of 1 second before making the request
		# Attempt the search with retries
		for _ in range(max_retries):
			try:
				handle = Entrez.esearch(db="protein", term=f"{seq_id}", idtype="acc")
				search_record = Entrez.read(handle)
				handle.close()
				break  # Break the loop if successful
			except urllib.error.HTTPError as e:
				if e.code == 429:  # HTTP 429: Too Many Requests
					print(f"Received HTTP 429 error. Retrying in 10 seconds...")
					time.sleep(10)
				else:
					continue  # Re-raise other HTTP errors
			except RuntimeError:
				continue
			except urllib.error.URLError:
				print(f'An internal server error occurred while handling the accession {hit_uid}')
				print(f"Received HTTP 403 error. Retrying in 5 seconds...")
				time.sleep(5)
				continue
		try:
			uid = search_record['IdList'][0]
		except (IndexError, AttributeError, KeyError):
			print(f"Entrez.esearch could not find results for {seq_id} on NCBI")
			not_found_list.append(seq_id)
			continue
		# Grab Nuccore UIDs from the source database
		for db_name in db_name_list:
			loop_nuc_gi, loop_nuc_acc, not_found_hit = elink_routine(db_name, 'protein', uid)
			if not_found_hit:
				not_found_list.append(not_found_hit)
				continue
			if loop_nuc_gi:
				source2target.setdefault(loop_nuc_gi, seq_id)
				nucleotide_uid_list.append(loop_nuc_gi)
				print(f"Nuccore GI and ptn accession: {loop_nuc_gi} {seq_id}")
				break

	# Ouputs nuccore uids and the ptn->nuc uid links
	print(f"Total n of not found hits: {len(not_found_list)}")
	return nucleotide_uid_list, source2target, list(set(not_found_list))


def nuc_to_gb(uid_list):
	progress = 0
	threshold_progress_report = 20
	threshold_progress_step = 20
	# Get Genbank records for each Nuccore UID
	gb_records = {}
	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> GATHER REGIONS FOR HITS. PROGRESS:")
	for uid in uid_list:
		progress += 1
		if progress == threshold_progress_report:
			current_time = datetime.datetime.now()
			print(f"[{current_time}] FINISHED {progress} ENTRIES OUT OF {len(uid_list)}")
			threshold_progress_report += threshold_progress_step
		handle = Entrez.efetch(db="nucleotide", id=f"{uid}", rettype="gb", retmode="text")
		record = SeqIO.read(handle, "genbank")
		gb_records.setdefault(uid, record)
	# Returns a list  of Genbank SeqRecords objects
	return gb_records


def gb_plier(query_to_gb_dict, uid_to_acc, win_size):
	gbk_target = {}
	prot_dict = {}
	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> COMPILE GENOMIC REGIONS IN GENBANK FILES")
	for hit_uid in query_to_gb_dict:
		gbk = query_to_gb_dict[hit_uid]
		for seq_feature in gbk.features:
			# Avoid blank feature that may occur in GenBank entries
			try:
				qualifiers = seq_feature.qualifiers
			except AttributeError:
				continue
			# Restrict search to protein-containing features
			if "protein_id" in qualifiers:
				prot_id = qualifiers["protein_id"][0]
				# Search for the protein-ids of interest
				if re.search(prot_id, uid_to_acc[hit_uid][0]):
					# Process feature information for future ref
					f_start = seq_feature.location.start.real
					f_end = seq_feature.location.end.real
					f_strand = seq_feature.strand
					f_seq = qualifiers["translation"][0]
					f_len = len(f_seq)
					highlight_feature = copy.deepcopy(seq_feature)
					highlight_feature.type = "highlight"
					# Set start/end coords using window size
					start = max(int(min([f_start, f_end])) - win_size, 0)
					end = min(int(max([f_start, f_end])) + win_size + 1, len(gbk.seq))

					# Create a SeqRecord object with the feature of interest
					gbk_focused = SeqRecord(
						id=gbk.id,
						annotations=gbk.annotations,
						dbxrefs=gbk.dbxrefs,
						seq=gbk.seq[start:end + 1],
						description=gbk.description
					)
					gbk_focused.features.append(highlight_feature)
					# Gather protein data for reference
					prep_prot_dict = {
									  "nuccore_acc": gbk.id,
									  # "region_seq": gbk.seq[start:end + 1],
									  "window_start": start,
									  "window_end": end,
									  "feature_start": f_start,
									  "feature_end": f_end,
									  "strand": f_strand,
									  "feature_len": f_len,
									  "mmseqs_hit": prot_id,
									  "ptn_sequence": f_seq
									  }
					prot_dict.setdefault(uid_to_acc[hit_uid][1], prep_prot_dict)
					gbk_target.setdefault(f"{gbk.id}_{start}-{end}", gbk_focused)
	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> GENOMIC REGIONS SUCCESSFULLY COMPILED IN GENBANK FILES")
	return gbk_target, prot_dict


def export_records(rec_list, output_path, file_format,):
	with open(f"{output_path}", "w") as handle:
		SeqIO.write(rec_list, handle, file_format)


def export_gbs(query_to_gb_dict, parent_path):
	if not os.path.exists(parent_path):
		os.mkdir(parent_path)
	for hit in query_to_gb_dict:
		gbk = query_to_gb_dict[hit]
		filename = f"{hit}.gb"
		with open(f"{parent_path}{os.sep}{filename}", "w") as gb_handle:
			SeqIO.write(gbk, gb_handle, "genbank")


def attach_label_to_fasta(seqrecords_list, id_to_label_df):
	new_seqrecords_list = []
	match_dict = {}
	id_to_label_df.apply(lambda row: match_dict.setdefault(row["sacc"], row["staxid"]), axis=1).to_dict()
	for record in seqrecords_list:
		try:
			record.id = f"{record.id}|{match_dict[record.id]}"
		except KeyError:
			pass
		new_seqrecords_list.append(record)
	return new_seqrecords_list


# # Load config file
# config_path = "../config/phylogeny_processing.yaml"
# with open(config_path, "r") as f:
# 	config = yaml.load(f, Loader=yaml.FullLoader)


def main():
	# # Snakemake Imports
	# #   SMK Inputs
	# mmseqs_search_path = str(snakemake.input.mmseqs_search_result)
	# #   SMK Outputs
	# retrieval_report = str(snakemake.output.retrieval_report)
	# output_fasta_hits = str(snakemake.output.hits_fasta)
	# #   SMK Params
	# col_names = str(snakemake.params.col_names)
	# parent_dir = str(snakemake.params.parent_dir)
	# window_size = int(snakemake.params.window_size)
	# #	SMK Wildcards
	# query_name = str(snakemake.wildcards.query_name)
	# db_prefix = str(snakemake.wildcards.db_prefix)

	#	Inputs
	mmseqs_search_path = str(sys.argv[1])
	#   Outputs
	retrieval_report = str(sys.argv[2])
	output_fasta_region = str(sys.argv[3])
	#   Params
	col_names = str(sys.argv[4])
	parent_dir = str(sys.argv[5])
	window_size = int(sys.argv[6])
	#	Wildcards
	# query_name = str(sys.argv[7])
	# db_prefix = str(sys.argv[8])

	# Entrez authentication
	print("Entrez login")
	Entrez.email = "thedoudnalab@gmail.com"

	# NCBI databases to search
	efecth_db_list = ['nuccore', 'nucleotide']

	# Import blastout table
	col_names_list = col_names.split(" ")
	blast_df = pd.read_csv(mmseqs_search_path,
	                       names=col_names_list,
	                       index_col=False).convert_dtypes().infer_objects()

	# Process pd Dataframe and take the unique set of hit IDs
	unique_mms_hits = list(set(blast_df[col_names].dropna().tolist()))

	# Query NCBI to get nuccore UIDs associated with the protein hits using ESearch/ELink
	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> Linking protein hit ids to Nuccore entries")
	nuc_uid_per_query, hit_to_link, hits_not_found = ptn_to_nuc(unique_mms_hits, efecth_db_list)

	# Get GenBank entries through EFetch
	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> Retrieving GenBank objects")
	gb_seqrec_per_query = nuc_to_gb(nuc_uid_per_query)

	# Narrow down Genbank files based on hit UIDs
	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> Extracting relevant features from GenBank objects")
	targeg_gb_dict, target_hit_dict = gb_plier(gb_seqrec_per_query, hit_to_link, window_size)

	# Parse hit sequences to FASTA file
	hit_seq_record_list = ncbi_fetch(unique_mms_hits, 'protein', 'fasta')

	# Generate summary dataframe
	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> GENERATE REPORTS")
	df = pd.DataFrame.from_dict(target_hit_dict, orient='index')

	# Export outputs
	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> EXPORT FILES")
	export_gbs(targeg_gb_dict, parent_dir)
	df.to_csv(retrieval_report, index=False)

	# hit_taxid_seq_record_list = attach_label_to_fasta(hit_seq_record_list, id_to_taxid)
	export_records(hit_seq_record_list, output_fasta_region, 'fasta')


if __name__ == "__main__":
	main()
