# **** Variables ****
# configfile: "config/dali_align.yaml"

# Run on wynthon
# module load Sali anaconda/py311-2024.02 mpi/openmpi-x86_64 CBI

# Set up batch index range
batch_index = list(range(config["batch_range"][0],config["batch_range"][1]))

# **** Imports ****
import glob

# Cluster run template
# nohup snakemake --snakefile dali_align.smk --configfile config/dali_align_cas12.yaml -j 110 --cluster "qsub -l h_rt={cluster.time} -j y -pe smp 4 -cwd" --cluster-config config/cluster.yaml --rerun-triggers 'mtime' --latency-wait 120 --use-singularity --singularity-args '--bind /scratch --bind /wynton/home/doudna/bellieny-rabelo/nidali_output --bind /wynton/home/doudna/bellieny-rabelo/nidali_db' --rerun-incomplete &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		#
		expand("{db_dir}/pdb_files_DAT/batch_{batch_index}/list_batch_{batch_index}.txt",
			db_dir=config["db_dir"],batch_index=batch_index),
		# #
		# expand("{db_dir}/results_{batch_index}",
		# 	db_dir=config["db_dir"],batch_index=batch_index),
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

# rule setup_outdir:
# 	input:
# 		target_entries_list = "{db_dir}/pdb_files_DAT/batch_{batch_index}/list_batch_{batch_index}.txt",
# 	output:
# 		tracked_outdir = "{db_dir}/results_{batch_index}"
# 	message:
# 		"""
# Anchor on:
# 	--> {input.target_entries_list}
# Set up output directory:
# 	--> {output}
# 		"""
# 	shell:
# 		"""
# 		touch {output}
# 		"""

rule dali_run:
	input:
		query_dat = "{db_dir}/query_DAT/{query_name}.dat",
		target_entries_list = "{db_dir}/pdb_files_DAT/batch_{batch_index}/list_batch_{batch_index}.txt",
	output:
		temp_output_dali ="{db_dir}/pdb_files_DAT/batch_{batch_index}/{query_name}.txt"
	params:
		dali_path = config['dali_path'],
		run_time_bechmark = "dali_run_time",
		target_entries_list = "{db_dir}/pdb_files_DAT/batch_{batch_index}/list_batch_{batch_index}.txt",
		query_dir="{db_dir}/query_DAT",
		dat_database = "{db_dir}/pdb_files_DAT/batch_{batch_index}",
		output_dali = lambda wildcards: glob.glob("{run}".format(
			run=config["run"]
		))
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
	--> {params.target_entries_list}
DAT database:
	--> {params.dat_database}
Threads:
	--> {threads}
OUTPUT:
	--> {params.output_dali}
		"""
	shell:
		"""
		module load CBI
		echo CREATE TEMP RESULT DIRECTORY
		mkdir -p {wildcards.db_dir}/pdb_files_DAT/batch_{wildcards.batch_index}/tmp_{wildcards.query_name}
		cd {wildcards.db_dir}/pdb_files_DAT/batch_{wildcards.batch_index}/tmp_{wildcards.query_name}
		echo CREATE SYMLINKS
		ln -s {params.query_dir} query_dat || true
		ln -s {params.dat_database} db_dat || true
		ln -s {params.target_entries_list} db_entry_list || true		
		start_time=`date +%s`
		{params.dali_path}/dali.pl \
		--cd1 {wildcards.query_name} \
		--db db_entry_list \
		--dat1 query_dat \
		--dat2 db_dat \
		--np {threads} \
		--oneway \
		--outfmt "summary,alignments" \
		--title {output.temp_output_dali}
		echo $(expr `date +%s` - $start_time) >> {params.run_time_bechmark}_{threads}p.txt
		touch {output.temp_output_dali}
		cp {wildcards.db_dir}/pdb_files_DAT/batch_{wildcards.batch_index}/tmp_{wildcards.query_name}/{wildcards.query_name}.txt {params.output_dali}/{wildcards.query_name}_batch_{wildcards.batch_index}.txt
        """

rule consolidate_reports:
	input:
		expand("{run}/{query_name}_batch_{batch_index}.txt", run='',query_name='',batch_index='')
	output:
		aggregated_report = ""
