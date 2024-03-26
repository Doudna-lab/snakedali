#!/bin/bash
#$ -l h_rt=04:00:00
#$ -j y
#$ -l mem_free=2G
#$ -l scratch=4G
#$ -r y
#$ -cwd

module load mpi/openmpi-x86_64
rm dali.lock
mpirun -n 3 apptainer exec /wynton/home/doudna/bellieny-rabelo/projects/nidali/nidali.sif dali.pl --cd1 13AFA --cd2 G8EQA --dat1 /wynton/home/doudna/bellieny-rabelo/nidali_db/query_DAT --dat2 /wynton/home/doudna/bellieny-rabelo/nidali_db/pdb_files_DAT/batch_1000 --oneway --outfmt "summary,alignments" --np 3 --clean
