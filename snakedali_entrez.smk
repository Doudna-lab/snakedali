# **** Imports ****
import glob

# DEFINE RULES THAT WON'T BE SUBMITTED TO COMPUTING NODES

# **** Variables ****
# Set up batch index range
batch_index = list(range(config["batch_range"][0],config["batch_range"][1] + 1))

# Prep run on UCSF's Wynthon HPC
# module load Sali CBI miniforge3/24.3.0-0 mpi/openmpi-x86_64
# conda activate snake

# Cluster run template
# nohup snakemake --snakefile snakedali.smk --configfile config/dali_bmarking.yaml --profile_sge profile_sge &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# Run MMESEQS2 clustering to reduce the initial number of hits to Representative sequences only
		expand("{run}/alignments/{query_name}/mmseqs/results/{query_name}_clustered_rep_seq.fasta",
			run=config["run"],query_name=config["query_name"]),
		# Recover genomic context from hits found through MMSEQS2
		expand("{run}/alignments/{query_name}/mmseqs/results/{query_name}_vs_mmseqsDB_summary_report.csv",
			run=config["run"],query_name=config["query_name"]),
		# Predict CRISPR arrays on the recovered sequences
		expand("{run}/alignments/{query_name}/minced/{query_name}_vs_mmseqsDB_array.manifest",
			run=config["run"],query_name=config["query_name"])

# noinspection SmkAvoidTabWhitespace
rule cluster_hits:
	input:
		mmseqs_hits_fasta = "{run}/alignments/{query_name}/mmseqs/results/{query_name}_vs_mmseqsDB_result-mms_hits.fasta",
	output:
		cluster_representatives = "{run}/alignments/{query_name}/mmseqs/results/{query_name}_clustered_rep_seq.fasta"
	params:
		tmpdir = config['tmpdir'],
		output_directory = "{run}/alignments/{query_name}/mmseqs/results",
		mmseqs_clust_params_string=config["mmseqs_clust_params"],
	conda:
		'envs/mmseqs.yaml'
	threads:
		config["threads"]
	shell:
		"""
		#module load CBI miniforge3/24.3.0-0 || True
		eval "$(conda shell.bash hook)"
		conda activate mmseqs
		mkdir -p {params.tmpdir}/{wildcards.query_name}/mmseqs_cluster
		cd {params.tmpdir}/{wildcards.query_name}/mmseqs_cluster
		mmseqs easy-cluster {input.mmseqs_hits_fasta} {wildcards.query_name}_clustered tmp_hits {params.mmseqs_clust_params_string}
		mv {wildcards.query_name}_clustered* {params.output_directory}
		"""

# noinspection SmkAvoidTabWhitespace
rule retrieve_genomic_context:
	input:
		cluster_representatives = "{run}/alignments/{query_name}/mmseqs/results/{query_name}_clustered_rep_seq.fasta",
		mmseqs_search_result_m8= "{run}/alignments/{query_name}/mmseqs/results/{query_name}_vs_mmseqsDB_result-mms_m8.tsv"
	output:
		retrieval_report = "{run}/alignments/{query_name}/mmseqs/results/{query_name}_vs_mmseqsDB_summary_report.csv"
	params:
		window_size = config["window_size"],
		col_names = "hit_id",
		log_file = "{run}/alignments/{query_name}/mmseqs/results/entrez.log",
		gbk_dir = "{run}/alignments/{query_name}/mmseqs/results/gbk",
		fasta_region_dir= "{run}/alignments/{query_name}/mmseqs/results/fna_region",
	message:
		"""
Recover genomic regions from MMSEQS2 hits:\n {input.cluster_representatives}
Use window size: {params.window_size}
Export summarized report to:\n {output.retrieval_report}
Export regions in FASTA format:\n {params.fasta_region_dir}
		"""
	shell:
		"""
		#module load CBI miniforge3/24.3.0-0 || True		
		eval "$(conda shell.bash hook)"
		conda activate biopympi
		grep ">" {input.cluster_representatives} | cut -d ' ' -f 1 | cut -d '>' -f 2 > tmp_{wildcards.query_name}_ids.tsv
		mpiexec python ../py/parallel_entrez_gather.py \
		tmp_{wildcards.query_name}_ids.tsv \
		{output.retrieval_report} \
		{params.col_names} \
		{params.gbk_dir} \
		{params.fasta_region_dir} \
		{params.window_size} \
		{input.mmseqs_search_result_m8} \
		{params.log_file}
		rm tmp_{wildcards.query_name}_ids.tsv
		"""

# noinspection SmkAvoidTabWhitespace
rule call_crispr_regions:
	input:
		retrieval_report = "{run}/alignments/{query_name}/mmseqs/results/{query_name}_vs_mmseqsDB_summary_report.csv",
	output:
		crispr_call_manifest = "{run}/alignments/{query_name}/minced/{query_name}_vs_mmseqsDB_array.manifest",
	params:
		fasta_region_dir = "{run}/alignments/{query_name}/mmseqs/results/fna_region",
		crispr_call_dir = "{run}/alignments/{query_name}/minced"
	message:
		"""
Predict CRISPR regions from genomic regions listed in the entrez report:\n {input.retrieval_report}
Export predictions to:\n {params.crispr_call_dir}
		"""
	shell:
		"""
		#module load CBI miniforge3/24.3.0-0 || True		
		eval "$(conda shell.bash hook)"
		conda activate minced
		files_list={params.fasta_region_dir}/*.fasta
		rm -f {output.crispr_call_manifest}
        for fasta_file in $files_list
        do
            base_name=$(basename "$fasta_file" .fasta)
            minced "$fasta_file" {params.crispr_call_dir}/"$base_name".fasta {params.crispr_call_dir}/"$base_name".gff
            echo $base_name >> {output.crispr_call_manifest}
        done;
		"""

# noinspection SmkAvoidTabWhitespace
# rule summarize_crispr_call:
# 	input:
# 		crispr_call_manifest = "{run}/alignments/{query_name}/minced/{query_name}_vs_mmseqsDB_array.manifest",
# 		retrieval_report = "{run}/alignments/{query_name}/mmseqs/results/{query_name}_vs_mmseqsDB_summary_report.csv",
# 	output:
# 		crispr_call_report = "{run}/alignments/{query_name}/minced/{query_name}_crispr_summary.txt",
