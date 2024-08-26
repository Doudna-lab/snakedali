# **** Imports ****
import glob

# **** Variables ****
# Set up batch index range
batch_index = list(range(config["batch_range"][0],config["batch_range"][1] + 1))

# Prep run on UCSF's Wynthon HPC
# module load Sali conda-forge/py310-24.3.0 mpi/openmpi-x86_64 CBI
# conda activate snake

# Cluster run template
# nohup snakemake --snakefile snakedali.smk --configfile config/dali_template.yaml --profile profile/ &

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
			run=config["run"],query_name=config["query_name"])

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
		apptainer exec --bind {params.dat_directory_parent} {params.singularity_path} import.pl --pdbfile {input.query_pdb} --pdbid {wildcards.query_name} --dat {params.tmp_dat_directory}		
		cd {params.dat_directory}
		mv {params.tmp_dat_directory}/{wildcards.query_name}* {output.query_dat}
		"""

# noinspection SmkAvoidTabWhitespace
rule create_list_input:
	input:
		dat_database = lambda wildcards: glob.glob("{db_dir}/pdb_files_DAT/batch_{batch_index}".format(
			db_dir=config["db_dir"], batch_index=wildcards.batch_index
		))
	output:
		target_entries_list = "{run}/dali_run_elements/batches/batch_{batch_index}/list_batch_{batch_index}.txt"
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
		echo LOAD MODULE
		module load CBI		
		
		echo CREATE SCRATCH FOLDER
		echo {params.tmpdir}/{wildcards.query_name}/batches/batch_{wildcards.batch_index}/
		mkdir -p {wildcards.run}/alignments/{wildcards.query_name}/batches/batch_{wildcards.batch_index}/
		mkdir -p {params.tmpdir}/{wildcards.query_name}/batches/batch_{wildcards.batch_index}/
		cd {params.tmpdir}/{wildcards.query_name}/batches/batch_{wildcards.batch_index}/
		
		echo CREATE SYMLINKS
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
		
		echo MOVE RESULTS TO PERMANENT PATH
		mv {params.tmpdir}/{wildcards.query_name}/batches/batch_{wildcards.batch_index}/{wildcards.query_name}A.txt {wildcards.run}/alignments/{wildcards.query_name}/batches/batch_{wildcards.batch_index}/{wildcards.query_name}A.txt
		"""

# noinspection SmkAvoidTabWhitespace
rule consolidate_reports:
	input:
		alignment_list = expand("{run}/alignments/{{query_name}}/batches/batch_{batch_index}/{{query_name}}A.txt",
			run=config["run"],batch_index=batch_index),
		query_dat = "{run}/query/{query_name}/{query_name}A.dat"
	output:
		aggregated_report = "{run}/results/{query_name}/{query_name}_daliout.xlsx"
	params:
		id_converstion_table = config["id_convert"]
	message:
		"""
Aggregate dali outputs for query {wildcards.query_name}: 
	--> Pull Query dat from: {input.query_dat}
	--> Generate Report on: {output.aggregated_report}
	
	Wildcads: {wildcards} 
		"""
	script:
		"py/aggregate_report.py"
