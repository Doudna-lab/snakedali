# == Native Modules
import argparse
import subprocess
# == Installed Modules

# == Project Modules


def main():
	# THIS IS A SNIPPET OF CODE TAKEN FROM dali2fasta.py
	# ALl the inputs for the process carried out by this scripts are generated upstream by dali2fasta
	# Parse the command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--tcoffee_params", dest="tcoffee_params", required=True)
	parser.add_argument("--tcoffee_bin", dest="tcoffee_bin", required=True)
	args = parser.parse_args()

	tcoffee_params = args.tcoffee_params
	tcoffee_bin = args.tcoffee_bin

	# Format the t-coffee command
	subprocess_cmd = re.sub(r"\|input\|", output_path, tcoffee_params)
	formatted_cmd = f"{tcoffee_bin} {subprocess_cmd}"
	# Execute the command using subprocess
	try:
		print(f"T-COFFEE COMMAND: {formatted_cmd}")
		subprocess.run(formatted_cmd, shell=True, check=True)
	except subprocess.CalledProcessError as e:
		print(f"Error running t_coffee: {e}")
		exit(0)


if __name__ == "__main__":
	main()
