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
		expand("{run}/pdb_files_DAT/batch_{batch_index}/list_batch_{batch_index}.txt",
			run=config["run"],batch_index=batch_index),
		#
		expand("{run}/nidali_search_output/{query_name}/{query_name}_batch_{batch_index}",
			run=config["run"],query_name=config["query_name"], batch_index=batch_index),

rule create_list_input:
	input:
		dat_database= lambda wildcards: glob.glob("{db_dir}/pdb_files_DAT/batch_{batch_index}".format(
			db_dir=config['db_dir'],batch_index=wildcards.batch_index))
	output:
		output_list = "{run}/pdb_files_DAT/batch_{batch_index}/list_batch_{batch_index}.txt"
	message:
		"""
Perform sequence alignments with DaliLIte:
Create list of inputs on: 
	--> {output.output_list} 
DAT database: 
	--> {input.dat_database}
		"""
	script:
		"py/create_list_file.py"

rule dali_run:
	input:
		query_dat = lambda wildcards: glob.glob("{db_dir}/query_DAT/{query_name}.dat".format(
			db_dir=config['db_dir'],query_name=wildcards.query_name)),
		dat_database=lambda wildcards: glob.glob("{db_dir}/pdb_files_DAT/batch_{batch_index}".format(
			db_dir=config['db_dir'],batch_index=wildcards.batch_index)),
	output:
		output_dir = "{run}/nidali_search_output/{query_name}/{query_name}_batch_{batch_index}"
	params:
		query_dir = "{run}/query_DAT"
	singularity:
		"nidali.sif"
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
	--> {wildcards.batch_index}
DAT database:
	--> {input.dat_database}
		"""
	shell:
		"""
		dali.pl \
		--cd1 {wildcards.query_name} \
		--db batch_{wildcards.batch_index} \
		--dat1 {params.query_dir} \
		--dat2 {input.dat_database} \
		--oneway \
		--outfmt "summary,alignments" \
        """
