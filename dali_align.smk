# **** Variables ****
# configfile: "config/dali_align.yaml"

# Run on wynthon
# module load Sali anaconda/py311-2024.02 mpi/openmpi-x86_64 CBI

# Set up batch index range
batch_index = list(range(config["batch_range"][0],config["batch_range"][1] + 1))

# **** Imports ****
import glob

# Cluster run template
# nohup snakemake --snakefile dali_align_bmarking.smk --configfile config/dali_template.yaml --profile profile/ &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		#
		expand("{run}/dali_run_elements/batches/batch_{batch_index}/list_batch_{batch_index}.txt",
			run=config["run"],batch_index=batch_index),
		#
		expand("{run}/alignments/{query_name}/batches/batch_{batch_index}/{query_name}.txt",
			run=config["run"],query_name=config["query_name"], batch_index=batch_index),
		#
		expand("{run}/results/{query_name}/{query_name}_daliout.xlsx",
			run=config["run"],query_name=config["query_name"])

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

rule dali_run:
	input:
		query_dat = lambda wildcards: glob.glob("{input_dir}/{query_name}.dat".format(
			input_dir=config["input_dir"],query_name=wildcards.query_name
		)),
		target_entries_list = "{run}/dali_run_elements/batches/batch_{batch_index}/list_batch_{batch_index}.txt"
	output:
		output_dali ="{run}/alignments/{query_name}/batches/batch_{batch_index}/{query_name}.txt"
	params:
		dali_path = config['dali_path'],
		run_time_bechmark = "dali_run_time",
		input_dir = config["input_dir"],
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
OUTPUT:
	--> {output.output_dali}
		"""
	shell:
		"""
		module load CBI
		echo CREATE TEMP RESULT DIRECTORY
		cd {wildcards.run}/alignments/{wildcards.query_name}/batches/batch_{wildcards.batch_index}/
		echo CREATE SYMLINKS
		ln -s {params.input_dir} query_dat || true
		ln -s {params.dat_database} db_dat || true
		ln -s {input.target_entries_list} db_entry_list || true		
		start_time=`date +%s`
		{params.dali_path}/dali.pl \
		--cd1 {wildcards.query_name} \
		--db db_entry_list \
		--dat1 query_dat \
		--dat2 db_dat \
		--np {threads} \
		--clean \
		--oneway \
		--outfmt "summary,alignments"
		echo $(expr `date +%s` - $start_time) >> {params.run_time_bechmark}_{threads}p.txt
		"""

rule consolidate_reports:
	input:
		expand("{run}/alignments/{{query_name}}/batches/batch_{batch_index}/{{query_name}}.txt",
			run=config["run"],batch_index=batch_index)
	output:
		aggregated_report = "{run}/results/{query_name}/{query_name}_daliout.xlsx"
	message:
		"""
Aggregate dali outputs for query {wildcards.query_name}: 
	--> {output.aggregated_report} 
		"""
	script:
		"py/aggregate_report.py"
