# select serial (NP=1)  or parallel version (NP>1)
NP=1
MPI="/usr/lib64/openmpi/bin/mpirun" # /usr/lib64/openmpi/bin/mpirun

# run tests in subdirectory
cd ./test
# clean up after previous trial
rm -f *
# import PDB entries
echo "> > > Importing PDB structures ..."
echo "cd ./test; ../bin/import.pl --pdbfile ../toy_PDB/pdb1ppt.ent --pdbid 1ppt --dat ./"
../bin/import.pl --pdbfile ../toy_PDB/pdb1ppt.ent --pdbid 1ppt --dat ./ > import.stdout 2> import.stderr
echo
echo "* * * Result of data import: ./test/1pptA.dat * * *"
echo
cat 1pptA.dat
echo

# create BLAST database of structures imported to DAT/
# (We use here a small mock database containing only a few structures.)
echo "> > > Creating Blast database ..."
echo "cd ./test; ls ../DAT | perl -pe 's/\.dat//' > pdb.list; ../bin/dat2fasta.pl ../DAT < pdb.list  | awk -v RS=\">\" -v FS=\"\n\" -v ORS=\"\" ' { if (\$2) print \">\"\$0 } ' > pdb.fasta"
echo "makeblastdb  -in pdb.fasta -out pdb.blast -dbtype prot"
ls ../DAT | perl -pe 's/\.dat//' > pdb.list
../bin/dat2fasta.pl ../DAT < pdb.list  | awk -v RS=">" -v FS="\n" -v ORS="" ' { if ($2) print ">"$0 } ' > pdb.fasta
# assuming blastp and makeblastp are in your PATH
# (Install the BLAST executables blastp and makeblastdb from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
# create BLAST database
makeblastdb  -in pdb.fasta -out pdb.blast -dbtype prot > makeblastdb.stdout 2> makeblastdb.stderr
echo "* * * Blast database done * * *"
echo
# pairwise comparison
echo "> > > Starting pairwise comparison ..."
echo "cd ./test; ../bin/dali.pl --pdbfile1 ../toy_PDB/pdb1ppt.ent --pdbfile2 ../toy_PDB/pdb1bba.ent.gz --dat1 ./ --dat2 ./ --outfmt summary,alignments" 
../bin/dali.pl --np $NP --MPIRUN_EXE $MPI --pdbfile1 ../toy_PDB/pdb1ppt.ent --pdbfile2 ../toy_PDB/pdb1bba.ent.gz --dat1 ./ --dat2 ./ --outfmt "summary,alignments" > pairwise.stdout 2> pairwise.stderr
echo
echo "* * * Result file of pairwise comparison: ./test/mol1A.txt * * *"
echo
cat mol1A.txt

echo "> > > Starting systematic search ..."
echo "cd ./test; ../bin/dali.pl --cd1 2nrmA --db pdb.list --TITLE systematic --dat1 ../DAT --dat2 ../DAT"
# pairwise comparison - result is 2nrmA.txt
../bin/dali.pl --np $NP --MPIRUN_EXE $MPI --cd1 2nrmA --db pdb.list --TITLE systematic --dat1 ../DAT --dat2 ../DAT > systematic.stdout 2> systematic.stderr
echo
echo "* * * Result file of systematic search: ./test/2nrmA.txt * * *"
echo
cat 2nrmA.txt

echo "> > > Starting hierarchical search ..."
echo "cd ./test; ../bin/dali.pl --cd1 5oheA --db pdb.list --TITLE hierarchical --hierarchical -repset ../matrix_example.list --dat1 ../DAT --dat2 ../DAT"
# pairwise comparison - result is 2nrmA.txt
../bin/dali.pl --np $NP --MPIRUN_EXE $MPI --cd1 5oheA --db pdb.list --TITLE hierarchical --hierarchical -repset ../matrix_example.list --dat1 ../DAT --dat2 ../DAT > hierarchical.stdout 2> hierarchical.stderr
echo
echo "* * * Result file of hierarchical search: ./test/5oheA.txt * * *"
echo
cat 5oheA.txt


# knowledge-based search of our small database - result is 5cnbA.txt
echo "> > > Starting knowledge-based search ..."
echo "cd ./test; ../bin/dali.pl --cd1 5cnbA --db pdb.list --walk --repset pdb.list --BLAST_DB pdb.blast --TITLE walk --dat1 ../DAT --dat2 ../DAT"
../bin/dali.pl --np $NP --MPIRUN_EXE $MPI --cd1 5cnbA --db pdb.list --walk --repset pdb.list --BLAST_DB pdb.blast --TITLE walk --dat1 ../DAT --dat2 ../DAT > walk.stdout 2> walk.stderr
echo
echo "* * * Result file of knowledge-based search: ./test/5cnbA.txt * * *"
echo
cat 5cnbA.txt

# all-aainst-all comparison -- results include newick and ordered
echo "> > > Starting all-against-all comparison ..."
echo "cd ./test; ../bin/dali.pl --query ../matrix_example.list --matrix --dat1 ../DAT --dat2 ../DAT --TITLE matrix"
../bin/dali.pl --np $NP --MPIRUN_EXE $MPI --query ../matrix_example.list --matrix --dat1 ../DAT --dat2 ../DAT --TITLE matrix > matrix.stdout 2> matrix.stderr
echo
echo "* * * Result files of all-against-all comparison: ./test/ordered ./test/newick ./test/newick_unrooted ./test/101mA.txt ./test/1a00A.txt ./test/1a07A.txt ./test/1allA.txt ./test/1binA.txt * * *"
echo
cat ordered
echo
cat newick
echo
cat newick_unrooted
echo
cat 101mA.txt
echo
cat 1a00A.txt
echo
cat 1a87A.txt
echo
cat 1allA.txt
echo
cat 1binA.txt
echo
