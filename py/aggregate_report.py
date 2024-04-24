# == Native Modules
# == Installed Modules
import re
import pandas as pd


# == Project Modules


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
	query_id = query_id[:-1]
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
							summary_dict[query_id][current_hit].setdefault('Full_Alignment', alignment_string)
							alignment_string = ""
							current_hit = hit_aligned
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


def main():
	# Snakemake I/O
	# === Inputs
	query_name = str(snakemake.wildcards.query_name)
	output_dali = str(snakemake.params.output_dali)
	# === Outputs
	aggregate_report = str(snakemake.output.aggregate_report)

	# DEBUG
	aln_file_path = "/Users/bellieny/projects/nidali/dump/Y156A_batch_994.txt" # "/wynton/home/doudna/bellieny-rabelo/nidali_output/alshimary/Y156A_batch_994.txt"
	query = 'Y156A'
	query_length = 104

	# Parse Dali's alignment output to dictionary
	parsed_dali_dict = parse_dali_out(aln_file_path, query, query_length)



	# Convert the dictionary to a DataFrame with two levels of indices
	parsed_dali_df = third_level_dict_to_df(parsed_dali_dict)


if __name__ == "__main__":
	main()
