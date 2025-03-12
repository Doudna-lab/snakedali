# **** Imports ****
import glob

# **** Variables ****
# Set up batch index range
batch_index = list(range(config["batch_range"][0],config["batch_range"][1] + 1))

# Prep run on UCSF's Wynthon HPC
# module load Sali CBI miniforge3/24.3.0-0 mpi/openmpi-x86_64
# conda activate snake

# Cluster run template
# nohup snakemake --snakefile snakedali_align.smk --configfile config/dali_bmarking.yaml --profile profile_sge &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# Exports PDB Formatted Files to DAT files
		expand("{run}/query/{query_name}/{query_name}A.dat",
			run=config["run"],query_name=config["query_name"]),
		# Creates a text list containing all database entry IDs in preparation for dalilite
		expand("{run}/dali_run_elements/batches/batch_{batch_index}/list_batch_{batch_index}.txt",
			run=config["run"],batch_index=batch_index),
		# Runs the structural alignment and allocate the final outputs according to the config file variables
		expand("{run}/alignments/{query_name}/batches/batch_{batch_index}/{query_name}A.txt",
			run=config["run"],query_name=config["query_name"], batch_index=batch_index),
		# Aggregate structural alignment results in a single xlsx spreadsheet
		expand("{run}/results/{query_name}/{query_name}_daliout.xlsx",
			run=config["run"],query_name=config["query_name"]),
		# Generate FASTA alignments from Dalilite results
		expand("{run}/alignments/{query_name}/fasta/pairwise_alignments/fasta_pairwise_manifest.txt",
			run=config["run"],query_name=config["query_name"]),
		# Expand list of sequence hits by running MMSEQS2 against different databases
		expand("{run}/alignments/{query_name}/mmseqs/results/{query_name}_vs_mmseqsDB_result-mms_hits.tsv",
			run=config["run"],query_name=config["query_name"]),
		# Recover FASTA-formatted sequences of the MMSEQS2 hits
		expand("{run}/alignments/{query_name}/mmseqs/results/{query_name}_vs_mmseqsDB_result-mms_hits.fasta",
			run=config["run"],query_name=config["query_name"]),


# noinspection SmkAvoidTabWhitespace
rule dali_import:
	input:
		query_pdb = lambda wildcards: glob.glob("{input_dir}/{query_name}.pdb".format(
			input_dir=config["input_dir"],query_name=wildcards.query_name
		))
	output:
		query_dat = "{run}/query/{query_name}/{query_name}A.dat",
	params:
		singularity_path = config['singularity_path'],
		dali_path=config['dali_path'],
		dat_directory_parent="{run}/query",
		tmp_dat_directory="{run}/query/tmp",
		dat_directory="{run}/query/{query_name}"
	message:
		"""
Perform PDB Data Import with DaliLIte:
PDB Query path:
	--> {input.query_pdb}
DAT Formatted Files:
	--> {output.query_dat}
		"""
	shell:
		"""
		mkdir -p {params.dat_directory}
		mkdir -p {params.tmp_dat_directory}
		cd {params.tmp_dat_directory}
		apptainer exec --bind {params.dat_directory_parent} {params.singularity_path} import.pl --pdbfile {input.query_pdb} --pdbid {wildcards.query_name} --dat .		
		cd {params.dat_directory}
		mv {params.tmp_dat_directory}/{wildcards.query_name}*.dat {output.query_dat}
		"""

# noinspection SmkAvoidTabWhitespace
rule create_list_input:
	input:
		dat_database = lambda wildcards: glob.glob("{db_dir}/pdb_files_DAT/batch_{batch_index}".format(
			db_dir=config["db_dir"], batch_index=wildcards.batch_index
		))
	output:
		target_entries_list = "{run}/dali_run_elements/batches/batch_{batch_index}/list_batch_{batch_index}.txt"
	params:
		output_directory = "{run}/dali_run_elements/batches/batch_{batch_index}"
	message:
		"""
Create list of inputs on: 
	--> {output.target_entries_list} 
Based on DAT database: 
	--> {input.dat_database}
		"""
	script:
		"py/create_list_file.py"

# noinspection SmkAvoidTabWhitespace
rule dali_run:
	input:
		query_dat = "{run}/query/{query_name}/{query_name}A.dat",
		target_entries_list = "{run}/dali_run_elements/batches/batch_{batch_index}/list_batch_{batch_index}.txt"
	output:
		output_dali ="{run}/alignments/{query_name}/batches/batch_{batch_index}/{query_name}A.txt"
	params:
		tmpdir = config['tmpdir'],
		dali_path = config['dali_path'],
		run_time_bechmark = "dali_run_time",
		input_dir = "{run}/query/{query_name}",
		dat_database = lambda wildcards: glob.glob("{db_dir}/pdb_files_DAT/batch_{batch_index}".format(
			db_dir=config["db_dir"], batch_index=wildcards.batch_index)),
	threads:
		config["threads"]
	message:
		"""
Perform sequence alignments with DaliLIte:
Query path:
	--> {input.query_dat}
Target DB:
	--> {input.target_entries_list}
DAT database:
	--> {params.dat_database}
Threads:
	--> {threads}
TMPDIR:
	--> {params.tmpdir}
OUTPUT:
	--> {output.output_dali}
		"""
	shell:
		"""
		# LOAD MODULE
		module load CBI || True
		
		# CREATE SCRATCH FOLDER
		mkdir -p {wildcards.run}/alignments/{wildcards.query_name}/batches/batch_{wildcards.batch_index}/
		mkdir -p {params.tmpdir}/{wildcards.query_name}/batches/batch_{wildcards.batch_index}/
		cd {params.tmpdir}/{wildcards.query_name}/batches/batch_{wildcards.batch_index}/
		
		# CREATE SYMLINKS
		ln -s {params.input_dir} query_dat || true
		ln -s {params.dat_database} db_dat || true
		ln -s {input.target_entries_list} db_entry_list || true		
		start_time=`date +%s`
		{params.dali_path}/dali.pl \
		--cd1 {wildcards.query_name}A \
		--db db_entry_list \
		--dat1 query_dat \
		--dat2 db_dat \
		--np {threads} \
		--clean \
		--oneway \
		--outfmt "summary,alignments"
		
		echo $(expr `date +%s` - $start_time) >> {params.run_time_bechmark}_{threads}p.txt
		mv {params.run_time_bechmark}_{threads}p.txt {wildcards.run}/alignments/{wildcards.query_name}/batches/batch_{wildcards.batch_index}/{params.run_time_bechmark}_{threads}p_{wildcards.query_name}A.txt
		
		# MOVE RESULTS TO PERMANENT PATH
		mv {params.tmpdir}/{wildcards.query_name}/batches/batch_{wildcards.batch_index}/{wildcards.query_name}A.txt {wildcards.run}/alignments/{wildcards.query_name}/batches/batch_{wildcards.batch_index}/{wildcards.query_name}A.txt
		"""

# noinspection SmkAvoidTabWhitespace
rule consolidate_reports:
	input:
		alignment_list = expand("{run}/alignments/{{query_name}}/batches/batch_{batch_index}/{{query_name}}A.txt",
			run=config["run"],batch_index=batch_index),
		query_dat = "{run}/query/{query_name}/{query_name}A.dat"
	output:
		aggregated_report = "{run}/results/{query_name}/{query_name}_daliout.csv"
	params:
		id_converstion_table = config["id_convert"]
	message:
		"""
Aggregate dali outputs for query {wildcards.query_name}: 
	--> Pull Query dat from: {input.query_dat}
	--> Generate Report on: {output.aggregated_report}
	
	Wildcads: {wildcards} 
		"""
	shell:
		"""
		module load CBI miniforge3/24.3.0-0 || True		
		eval "$(conda shell.bash hook)"
		conda activate biopympi
		python -u ../py/aggregate_report.py \
		--query_name {wildcards.query_name} \
		--output_dali_list {input.alignment_list} \
		--query_dat_path {input.query_dat} \
		--aggregated_report {output.aggregated_report} \
		--id_converstion_table_path {params.id_converstion_table}		
		"""

# noinspection SmkAvoidTabWhitespace
rule dali_to_fasta:
	input:
		alignment_list = expand("{run}/alignments/{{query_name}}/batches/batch_{batch_index}/{{query_name}}A.txt",
			run=config["run"],batch_index=batch_index)
	output:
		fasta_manifest = "{run}/alignments/{query_name}/fasta/pairwise_alignments/fasta_pairwise_manifest.txt",
		aggregated_multi_fasta = "{run}/alignments/{query_name}/fasta/{query_name}_hits.fasta",
	params:
		pairwise_dir = "{run}/alignments/{query_name}/fasta/pairwise_alignments",
		fseek_clustered_representatives = config["fseek_clustered_representatives"],
		id_converstion_table = config["id_convert"],
		afdb_pkl = config["afdb_pkl"]
		# tree_unrooted = "{run}/alignments/{query_name}/tcoffee/newick_unrooted",
		# tcoffee_bin = config["tcoffee_bin"],
		# tcoffee_params = config["tcoffee_1st_params"]
	threads:
		config["threads"]
	resources:
		mem_mb=340000
	shell:
		"""	
		module load CBI miniforge3/24.3.0-0 || True		
		eval "$(conda shell.bash hook)"
		conda activate biopympi
		python -u ../py/dali2fasta.py \
		 --output_dir {params.pairwise_dir} \
		 --multi_fasta_out {output.aggregated_multi_fasta} \
		 --manifest_out {output.fasta_manifest} \
		 --input_prefix {wildcards.query_name} \
		 --afdb_pkl {params.afdb_pkl} \
		 --fseek_clusters {params.fseek_clustered_representatives} \
		 --id_converstion {params.id_converstion_table} \
		 --files_list {input.alignment_list}	
		"""

# noinspection SmkAvoidTabWhitespace
rule broad_seq_search_aa:
	input:
		aggregated_multi_fasta = "{run}/alignments/{query_name}/fasta/{query_name}_hits.fasta",
	output:
		mmseqs_search_result = "{run}/alignments/{query_name}/mmseqs/results/{query_name}_vs_mmseqsDB_result-mms_hits.tsv",
		mmseqs_search_result_m8= "{run}/alignments/{query_name}/mmseqs/results/{query_name}_vs_mmseqsDB_result-mms_m8.tsv"
	params:
		tmpdir = config['tmpdir'],
		mms_n_iterations = config['mms_n_iterations'],
		mmseqs_source_db=config["mmseqs_db_path"],
		mmseqs_params_string=config["mmseqs_search_params"]
	message:
		"""
Create MMSEQS query databse from FASTA file:\n {input.aggregated_multi_fasta}
Pull Source Database: {params.mmseqs_source_db}
Writes single column hit results on: {output.mmseqs_search_result}
		"""
	shell:
		"""
		module load CBI miniforge3/24.3.0-0 || True
		eval "$(conda shell.bash hook)"
		conda activate mmseqs
		mkdir -p {params.tmpdir}/{wildcards.query_name}/mmseqs
		cd {params.tmpdir}/{wildcards.query_name}/mmseqs
		
		mmseqs createdb {input.aggregated_multi_fasta} {wildcards.query_name}_queryDB
		mmseqs createindex {wildcards.query_name}_queryDB tmp
		
		mmseqs search {wildcards.query_name}_queryDB {params.mmseqs_source_db} {wildcards.query_name}_resultDB tmp --num-iterations {params.mms_n_iterations} {params.mmseqs_params_string}
		mmseqs convertalis {wildcards.query_name}_queryDB {params.mmseqs_source_db} {wildcards.query_name}_resultDB {wildcards.query_name}_resultDB.m8
		mmseqs convertalis {wildcards.query_name}_queryDB {params.mmseqs_source_db} {wildcards.query_name}_resultDB --format-output "target" {wildcards.query_name}_resultDB_hits
		mv {wildcards.query_name}_resultDB_hits {output.mmseqs_search_result}
		mv {wildcards.query_name}_resultDB.m8 {output.mmseqs_search_result_m8}
		"""

# noinspection SmkAvoidTabWhitespace
rule gather_hit_sequences:
	input:
		mmseqs_search_result = "{run}/alignments/{query_name}/mmseqs/results/{query_name}_vs_mmseqsDB_result-mms_hits.tsv",
	output:
		mmseqs_hits_fasta = "{run}/alignments/{query_name}/mmseqs/results/{query_name}_vs_mmseqsDB_result-mms_hits.fasta",
	params:
		mmseqs_fasta_db = lambda wildcards: glob.glob("{mmseqs_db}.faa".format(mmseqs_db=config["mmseqs_db_path"])),
		tmpdir = config['tmpdir'],
	message:
		"""
Recover FASTA-formatted sequences from mmseqs search hits:
{input.mmseqs_search_result}
Use the mmseqs database reference:
{params.mmseqs_fasta_db}
Output hits FASTA to:
{output.mmseqs_hits_fasta}			
		"""
	shell:
		"""
		module load CBI miniforge3/24.3.0-0 || True
		eval "$(conda shell.bash hook)"
		conda activate seqtk
		mkdir -p {params.tmpdir}/{wildcards.query_name}/mmseqs_gather
		cd {params.tmpdir}/{wildcards.query_name}/mmseqs_gather
		seqtk subseq {params.mmseqs_fasta_db} {input.mmseqs_search_result} > tmp_hits.fasta
		
		conda activate seqkit
		seqkit rmdup -s < tmp_hits.fasta > tmp_dedup_hits.fasta
		mv tmp_dedup_hits.fasta {output.mmseqs_hits_fasta}
		"""
