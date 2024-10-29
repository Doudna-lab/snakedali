# == Native Modules
import pickle
import datetime
import re
import sys
from itertools import batched
# == Installed Modules
from Bio import SeqIO
# == Project Modules

#DEBUG
# ROOT_DIR = "/Users/bellieny/projects/snakedali/dump"
# fasta_file = f"{ROOT_DIR}/toy_sequences.fasta"
# breakdown_parts = 11


def process_breakdown(total_size, n_of_parts):
	break_increment = (int(total_size) // int(n_of_parts)) + (int(total_size) % int(n_of_parts))
	return break_increment


def chunks(data, SIZE=10000):
	return map(dict, batched(data.items(), SIZE))


def main():
	fasta_file = sys.argv[1]
	output_pickle = sys.argv[2]
	breakdown_parts = int(sys.argv[3])

	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> IMPORTING AFDB FASTA FILE: {fasta_file}")
	indexed_dict = {}
	for record in SeqIO.parse(fasta_file, "fasta"):
		clean_id = re.sub(r"\w+:\w+\-(\w+)\-\w+", r"\1", record.id)
		indexed_dict.setdefault(clean_id, record)

	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> IMPORTING COMPLETE")

	breakdown_size = process_breakdown(len(indexed_dict), breakdown_parts)
	track_file_numner = 0
	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> BREAKING DOWN FASTA INTO CHUNKS")
	for item in chunks(indexed_dict, breakdown_size):
		track_file_numner += 1
		filename = f"{output_pickle}{track_file_numner}"
		current_time = datetime.datetime.now()
		print(f"[{current_time}] --> EXPORTING CHUNK {track_file_numner}")
		with open(filename, "wb") as pickle_handle:
			pickle.dump(item, pickle_handle)
		current_time = datetime.datetime.now()
		print(f"[{current_time}] --> EXPORT COMPLETE")

	current_time = datetime.datetime.now()
	print(f"[{current_time}] --> FINISHED SERIALIZING AFDB FASTA")

if __name__ == "__main__":
	main()
