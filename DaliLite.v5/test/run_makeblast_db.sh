#! /usr/bin/env bash
#$ -S /bin/bash  # run job as a Bash shell [IMPORTANT]
#$ -cwd          # run job in the current working directory

module load CBI
module load blast
makeblastdb  -in pdb.seq -out pdb.blast -dbtype prot
