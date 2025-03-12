# == Native Modules
import argparse
import gc
import os
from pathlib import Path
import re
# == Installed Modules
import pandas as pd
# == Project Modules


# # DEBUG INPUTS:
# retrieval_report = "/groups/doudna/projects/daniel_projects/yoon_ph/type3_run_alt2/alignments/6ifb/mmseqs/results/6ifb_vs_mmseqsDB_summary_report.csv_backup"
# crispr_call_manifest = "/groups/doudna/projects/daniel_projects/yoon_ph/type3_run_alt2/alignments/6ifb/minced/6ifb_vs_mmseqsDB_array.manifest"
# fseek_clusters = "/groups/doudna/team_resources/shared_databases/afdb/fseek_clusters/fseek_clustered_afdb_representatives.txt"
# mmseqs_search_result_m8 = "/groups/doudna/projects/daniel_projects/yoon_ph/type3_run_alt2/alignments/6ifb/mmseqs/results/6ifb_vs_mmseqsDB_result-mms_m8.tsv"
# crispr_call_dir = "/groups/doudna/projects/daniel_projects/yoon_ph/type3_run_alt2/alignments/6ifb/minced"
# dali_summary_report = "/groups/doudna/projects/daniel_projects/yoon_ph/type3_run_alt2/results/6ifb/6ifb_daliout.xlsx"
# window_size = 20001
#

def format_conversion_dict(cluster_reps_tbl):
	cluster_reps_pre = pd.read_csv(cluster_reps_tbl, sep="\t", header=None, low_memory=False)
	cluster_reps_pre.index = cluster_reps_pre.iloc[:,1]
	# cluster_reps_dict = defaultdict(str)
	cluster_reps_dict = {}
	for row in cluster_reps_pre.itertuples():
		cluster_reps_dict.setdefault(row._1, row.Index)
	return cluster_reps_dict


def validate_gff(filepath: Path):
	with open(filepath, "r") as f:
		for line in f:
			if re.search(r"^#", line):
				continue
			return 1
	return 0


def compile_region_files(dir):
	regions_dict = {}
	region_bin = {}
	for root, dirs, files in os.walk(dir, topdown=False):
		for filename in files:
			filepath = Path(root, filename)
			# print(filepath.suffix)
			# return
			if filepath.suffix == ".gff":
				gff_validation = validate_gff(filepath)
				if gff_validation == 1:
					regions_dict.setdefault(filepath.stem, filepath)
					region_bin.setdefault(filepath.stem, 1)
					continue
				region_bin.setdefault(filepath.stem, 0)
	df_region_bin = pd.DataFrame.from_dict(region_bin, orient="index")
	return regions_dict, df_region_bin


def main():
	# Parse the command-line arguments
	parser = argparse.ArgumentParser()
	# parser.add_argument("--crispr_calls", dest="crispr_call_manifest", required=True)
	parser.add_argument("--entrez_report", dest="retrieval_report", required=True)
	parser.add_argument("--dali_summary", dest="dali_summary_report", required=True)
	parser.add_argument("--fseek_clusters", dest="fseek_clusters", required=True)
	parser.add_argument("--crispr_call_dir", dest="crispr_call_dir", required=True)
	parser.add_argument("--window_size", dest="window_size", required=True)
	parser.add_argument("--report_output", dest="report_output", required=True)
	args = parser.parse_args()

	# Snakemake I/O
	# === Inputs
	retrieval_report = str(args.retrieval_report)
	dali_summary_report = str(args.dali_summary_report)
	# === Outputs
	report_output = str(args.report_output)
	# === Params
	fseek_clusters = str(args.fseek_clusters)
	crispr_call_dir = str(args.crispr_call_dir)
	window_size = int(args.window_size) * 2

	# === Import Dalilite's summary
	df_dali = pd.read_excel(dali_summary_report)
	df_dali = df_dali.rename(columns={"Unnamed: 1": "cluster_rep"})
	df_dali = df_dali[["cluster_rep", "Z"]]

	# === Import Entrez Garthering Report
	df_retrieval_report = pd.read_csv(retrieval_report, low_memory=False)
	#  - Generate clean AFDB ID
	df_retrieval_report["clean_afdb_id"] = df_retrieval_report["0"].apply(
		lambda x: re.sub(r"^AFDB:AF-(\w+)-\w+$", r"\1", x)
	)
	df_retrieval_report = df_retrieval_report.rename(columns={"0": "afdb_id"})
	#  - Reinstate the region ID
	df_retrieval_report["region_id"] = df_retrieval_report["nuccore_acc"] + "_" + df_retrieval_report["window_start"].astype(str) + "-" + df_retrieval_report["window_end"].astype(str)
	#  - Calculate the effective window size of the regions
	df_retrieval_report["effective_window_size"] = abs(df_retrieval_report["window_start"].astype(int) - df_retrieval_report["window_end"].astype(int)) - (abs(df_retrieval_report["feature_end"].astype(int) - df_retrieval_report["feature_start"].astype(int)))

	# === Imports minced calls and compile the relevant information into a binary score
	region_files, df_region_bin = compile_region_files(crispr_call_dir)
	df_consolidated = df_retrieval_report.merge(df_region_bin, left_on='region_id', right_index=True, how='left')
	df_consolidated = df_consolidated.rename(columns={0: "minced_call"})

	# === Imports foldseek clusters
	#  - Format -> { member: cluster_rep }
	fseek_clusters_dict = format_conversion_dict(fseek_clusters)
	df_consolidated['cluster_rep'] = df_consolidated['clean_afdb_id'].map(fseek_clusters_dict)
	#  - Thrash large dict release memory
	del fseek_clusters_dict
	gc.collect()

	# === Integrate dalilite's Z to the final report
	df_consolidated = df_consolidated.merge(df_dali, on='cluster_rep', how='left')

	#  - Calculate the crispr score
	df_consolidated['score'] = df_consolidated['effective_window_size'] / window_size * df_consolidated['minced_call'] * (df_consolidated['Z'] / df_consolidated['Z'].max())

	# === Export Report
	df_consolidated.to_csv(report_output, index=False)


if __name__ == "__main__":
	main()
