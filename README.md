# nidali
New Improved DALI

DaliLite.v5 is a standalone program for protein structural alignment and database search.

### Replace /home/you by the parent directory where DaliLite is installed

<details>
<summary>download</summary>
<br>
	cd /home/you
	wget http://ekhidna2.biocenter.helsinki.fi/dali/DaliLite.v5.tar.gz
	tar -zxvf DaliLite.v5.tar.gz
</details>
<details>
<summary>compile</summary> 
	cd /home/you/DaliLite.v5/bin
        make # ignore Warnings
        # if using openmpi (check OPENMPI_PATH in Makefile)
	# make parallel
</details>
<details>
<summary>test</summary> 
        cd /home/you/DaliLite.v5
	# the script assumes that blastp and makeblastdb are in your PATH
	# if not, get them from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
	./test.csh
	# compare output to ./test_output
	# if using openmpi, set NP=2 in test.csh and rerun
</details>
<details>
<summary>usage</summary>
/home/you/DaliLite.v5/bin/import.pl
/home/you/DaliLite.v5/bin/dali.pl [-h]
