# **** Variables ****
configfile: "config/cluster.yaml"
configfile: "config/dali_align.yaml"

# Set up batch index range
batch_index = list(range(config["batch_range"][0],config["batch_range"][0]))

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile expanded_search.smk -j 5 --cluster "sbatch -t {cluster.time}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# Generate BLAST database
		# expand("{run}/custom_db/{db_prefix}.phr",
		# 	run=config["run"],db_prefix=config["db_prefix"]),
		# PSI-Blast the MSA input using a list of sequence databases
		expand("{run}/nidali_search_output/",
			run=config["run"],db_prefix=config["db_prefix"], input_prefix=config["input_prefix"], cluster=config["cluster"]),

rule dali_run:
	input:
		query_dat = "{run}/query_DAT/{query_name}.dat",
		dat_database = "{run}/pdb_files_DAT/batch_{batch_index}",
	output:
		output_dir = "{run}/nidali_search_out/{query_name}_batch_{batch_index}"
	params:
		dali_bin = config["dali_path"],
		dat_id_list_prefix = config["dat_list"],
		query_dir = "{run}/query_DAT"
	singularity:
		"docker://"
	threads:
		config["threads"]
	message:
		"""
Perform sequence alignments with DaliLIte:
Create list of inputs on: 
	--> {input.dat_database}/{params.dat_id_list_prefix}_{wildcards.batch_index} 
Bin path: 
	--> {params.dali_bin}
Query path: 
	--> {input.query_dat}
DAT database: 
	--> {input.dat_database}
		"""
	shell:
		"""
		python py/create_list_file.py {input.dat_database} {wildcards.batch_index} {params.dat_id_list_prefix}		
		{params.dali_bin}/dali.pl \
		--cd1 {wildcards.query_name} \
		--db {params.dat_id_list_prefix}_{wildcards.batch_index} \
		--dat1 {params.query_dir} \
		--dat2 {input.dat_database} \
		--oneway \
		--outfmt "summary,alignments" \
		--clean
        """