# **** Variables ****
configfile: "config/dali_align.yaml"

# Set up batch index range
batch_index = list(range(config["batch_range"][0],config["batch_range"][1]))

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile dali_align.smk -j 5 --cluster "qsub -l h_rt={cluster.time} -j y -cwd" --cluster-config config/cluster.yaml --latency-wait 120 --use-singularity &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		#
		expand("{db_dir}/pdb_files_DAT/batch_{batch_index}/list_batch_{batch_index}.txt",
			db_dir=config["db_dir"],batch_index=batch_index),
		#
		expand("{db_dir}/results_{batch_index}",
			db_dir=config["db_dir"],batch_index=batch_index),
		#
		expand("{db_dir}/pdb_files_DAT/batch_{batch_index}/{query_name}.txt",
			db_dir=config["db_dir"],query_name=config["query_name"], batch_index=batch_index),

rule create_list_input:
	input:
		dat_database = "{db_dir}/pdb_files_DAT/batch_{batch_index}"
	output:
		target_entries_list = "{db_dir}/pdb_files_DAT/batch_{batch_index}/list_batch_{batch_index}.txt"
	message:
		"""
Create list of inputs on: 
	--> {output.target_entries_list} 
Based on DAT database: 
	--> {input.dat_database}
		"""
	script:
		"py/create_list_file.py"

rule setup_outdir:
	input:
		target_entries_list = "{db_dir}/pdb_files_DAT/batch_{batch_index}/list_batch_{batch_index}.txt",
	output:
		directory("{db_dir}/results_{batch_index}"),
	message:
		"""
Anchor on: 
	--> {input.target_entries_list} 
Set up output directory: 
	--> {output}	
		"""
	shell:
		"""
		mkdir -p {output}
		"""

rule dali_run:
	input:
		query_dat = "{db_dir}/query_DAT/{query_name}.dat",
		tracked_outdir="{db_dir}/results_{batch_index}",
		target_entries_list = "{db_dir}/pdb_files_DAT/batch_{batch_index}/list_batch_{batch_index}.txt"
	output:
		temp_output_dali ="{db_dir}/pdb_files_DAT/batch_{batch_index}/{query_name}.txt"
	params:
		query_dir="{db_dir}/query_DAT",
		dat_database = "{db_dir}/pdb_files_DAT/batch_{batch_index}",
		output_dali = lambda wildcards: glob.glob("{run}".format(
			run=config["run"]
		))
	singularity:
		"nidali_mpi.sif"
	threads:
		config["threads"]
	message:
		"""
Perform sequence alignments with DaliLIte:
Query path:
	--> {params.query_dir}
Query ID:
	--> {wildcards.query_name}
Target DB:
	--> {input.target_entries_list}
DAT database:
	--> {params.dat_database}
OUTPUT:
	--> {params.output_dali}
		"""
	shell:
		"""
		cd {wildcards.db_dir}/pdb_files_DAT/batch_{wildcards.batch_index}
		dali.pl \
		--cd1 {wildcards.query_name} \
		--db {input.target_entries_list} \
		--dat1 {params.query_dir} \
		--dat2 {params.dat_database} \
		--oneway \
		--outfmt "summary,alignments" \
		--title {output.temp_output_dali} \
		--clean
		mv {output.temp_output_dali} {params.output_dali}
        """
