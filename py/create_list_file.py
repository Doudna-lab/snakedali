# == Native Modules ==
import os
from sys import argv
# == Installed Modules ==
# == Project Modules ==


def main():
	# Snakemake I/O
	# === Inputs
	source_dat_path = str(snakemake.input.dat_database)
	# === Outputs
	dat_list_file = str(snakemake.output.target_entries_list)

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
