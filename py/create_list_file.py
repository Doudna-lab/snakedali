# == Native Modules ==
import os
from sys import argv
# == Installed Modules ==
# == Project Modules ==


def mkdir_p(directory):
	"""
	Creates a directory if it does not exist, similar to 'mkdir -p' in shell.

	Parameters:
		directory (str): The path of the directory to create.
	"""
	try:
		os.makedirs(directory, exist_ok=True)
		print(f"Directory '{directory}' is ready.")
	except Exception as e:
		print(f"An error occurred while creating the directory '{directory}': {e}")


def main():
	# Snakemake I/O
	# === Inputs
	source_dat_path = str(snakemake.input.dat_database)
	# === Outputs
	dat_list_file = str(snakemake.output.target_entries_list)
	# === Params
	output_directory = str(snakemake.params.output_directory)

	# Set Variables and Output dir
	mkdir_p(output_directory)
	names_list = []

	for root, dirs, files in os.walk(source_dat_path, topdown=False):
		for name in files:
			if name.endswith(".dat"):
				names_list.append(name)

	with open(dat_list_file, 'w') as dat_file_handle:
		for dat_name in names_list:
			dat_file_handle.write(f"{dat_name}\n")


if __name__ == "__main__":
	main()
