# == Native Modules
# == Installed Modules
import argparse
import re
import pandas as pd
# == Project Modules


def fetch_seqlength_from_dat(dat_filepath):
	with open(dat_filepath, 'r') as dat_file_handle:
		for line in dat_file_handle:
			if re.match(r"^-sequence", line):
				sequence_string = line.split('\"')[1].strip()
				return len(sequence_string)


def parse_summary(line, query_id, query_length):
	target_dict = {}
	target_id = "target_id_not_found"

	try:
		_, chain, z, rmsd, alnlen, tlen, pident, *_ = line
		# Removing the chain identifier (-A for chain A, etc)
		target_id = chain[:-2]

		z = float(z)
		rmsd = float(rmsd)
		alnlen = int(alnlen)
		tlen = int(tlen)
		pident = float(pident)
		# print(f"FOUND: {query_id} vs {target_id}: {z} {rmsd} {alnlen} {tlen} {pident}")

		# Generate a "cov" column that is alnlen/(max(query_length, tlen)). If there is no
		# labeled query_length, it will just end up alnlen/tlen.
		# print(f"Compare {type(query_length)}_{query_length} to {type(tlen)}_{tlen}")
		query_cov_perc = alnlen * 100 / query_length
		target_dict = {"Z": z,
					   "RMSD": rmsd,
					   "Alignment_Length": alnlen,
					   "Hit_Length": tlen,
					   "Query_Length": query_length,
					   "Perc_Ident": pident,
					   "Query_Cov_Perc": round(query_cov_perc, 2)}

	except ValueError:
		pass

	return target_id, target_dict


def parse_dali_out(file_path, query_id, query_length):
	processing_mode = ""
	alignment_string = ""
	current_hit = ""
	summary_dict = {}

	with open(file_path) as infile:
		for line in infile:
			if line == "":
				continue

			if re.search("# No:", line):
				processing_mode = "summary"
				continue
			elif re.search("# Pairwise alignments", line):
				processing_mode = "alignments"
				continue

			if processing_mode == "summary":
				line = line.split()
				target_id, target_dict = parse_summary(line, query_id, query_length)
				if len(target_dict) >= 1:
					summary_dict.setdefault(query_id, {}).setdefault(target_id, target_dict)

			elif processing_mode == "alignments":
				match = re.search(r"No (\d+): Query=\S+ Sbjct=(\S+) Z-score=\S+", line)
				if match:
					entry_number = match.group(1)
					hit_aligned = match.group(2)[:-1]
					if current_hit != hit_aligned:
						if not int(entry_number) == int(1):
							try:
								summary_dict[query_id][current_hit].setdefault('Full_Alignment', alignment_string)
								alignment_string = ""
								current_hit = hit_aligned
							except KeyError:
								continue
						if int(entry_number) == int(1):
							current_hit = hit_aligned
							continue

				elif not match:
					alignment_string += line.strip('\n') + "\n"
	return summary_dict


def third_level_dict_to_df(third_level_dictionary):
	df = pd.DataFrame.from_dict({(i, j): third_level_dictionary[i][j]
								 for i in third_level_dictionary.keys()
								 for j in third_level_dictionary[i].keys()},
								orient='index')
	return df


def merge_dictionaries(dict1, dict2):
	merged_dict = dict1.copy()  # Make a copy of the first dictionary
	for key, value in dict2.items():
		if key in merged_dict:
			if isinstance(value, dict) and isinstance(merged_dict[key], dict):
				merged_dict[key].update(value)  # Merge inner dictionaries if both values are dictionaries
			else:
				merged_dict[key] = [merged_dict[key], value]  # Convert to list if values are not dictionaries
		else:
			merged_dict[key] = value
	return merged_dict


def main():
	# # Snakemake I/O
	# # === Inputs
	# query_name = str(snakemake.wildcards.query_name)
	# output_dali_list = list(snakemake.input.alignment_list)
	# query_dat_path = str(snakemake.input.query_dat)
	# # === Outputs
	# aggregated_report = str(snakemake.output.aggregated_report)
	# # === Params
	# id_converstion_table_path = str(snakemake.params.id_converstion_table)

	# Parse the command-line arguments
	parser = argparse.ArgumentParser()
	# parser.add_argument("--crispr_calls", dest="crispr_call_manifest", required=True)
	parser.add_argument("--query_name", dest="query_name", required=True)
	parser.add_argument("--output_dali_list", nargs='+', dest="output_dali_list", required=True)
	parser.add_argument("--query_dat_path", dest="query_dat_path", required=True)
	parser.add_argument("--aggregated_report", dest="aggregated_report", required=True)
	parser.add_argument("--id_converstion_table_path", dest="id_converstion_table_path", required=True)
	args = parser.parse_args()

	# Snakemake I/O
	# === Inputs
	query_name = args.query_name
	output_dali_list = list(args.output_dali_list)
	query_dat_path = args.query_dat_path
	# === Outputs
	aggregated_report = args.aggregated_report
	# === Params
	id_converstion_table_path = args.id_converstion_table_path

	# === Set key variables
	merged_dali_parsed_dict = {}
	query = query_name[:-1]
	# 	== Calculate sequence length from .dat file
	query_length = fetch_seqlength_from_dat(query_dat_path)

	# === Process dalilite outputs
	for aln_file_path in output_dali_list:
		# Parse Dali's alignment output to dictionary
		parsed_dali_dict = parse_dali_out(aln_file_path, query, query_length)
		# print(f"Length of intermediate dict: {len(parsed_dali_dict[query])}")
		if len(merged_dali_parsed_dict) == 0:
			merged_dali_parsed_dict = parsed_dali_dict
			continue
		merged_dali_parsed_dict = merge_dictionaries(parsed_dali_dict, merged_dali_parsed_dict)

	if id_converstion_table_path:
		id_converstion_table = pd.read_csv(id_converstion_table_path, sep="\t", low_memory=False, names=['0', '1'])
		converted_id_dict = dict(zip(id_converstion_table.iloc[:, 1], id_converstion_table.iloc[:, 0]))
		final_converted_id_dict = {}
		try:
			data_valid = True
			query_dict = merged_dali_parsed_dict[query]
		except KeyError:
			data_valid = False
			pd.DataFrame().to_excel(aggregated_report)
		if data_valid:
			for key, value in merged_dali_parsed_dict[query].items():
				try:
					final_converted_id_dict.setdefault(query, {}).setdefault(converted_id_dict[key], value)
				except KeyError:
					continue
			# Convert the dictionary to a DataFrame with two levels of indices
			parsed_dali_df = third_level_dict_to_df(final_converted_id_dict)
			# Add in Alphafold Links
			parsed_dali_df['Alphafold_link'] = "https://alphafold.ebi.ac.uk/entry/" + parsed_dali_df.index.get_level_values(1)
			# Export Parsed Dataframe
			parsed_dali_df.to_csv(aggregated_report)


if __name__ == "__main__":
	main()
