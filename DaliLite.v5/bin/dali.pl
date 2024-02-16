#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin"; # modules and scripts in same directory as this script
use mpidali;
use FSSP;
use walk;
use auxiliary_lib;
use warnings;
use 5.010;
use Getopt::Long qw(GetOptions);


# BLAST config
# Install BLAST from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
my $BLAST_DB="pdb.blast"; # "/home/luholm/dalitest/pdb.blast"; # makeblastdb -in /data/uniprot/pdb.fasta -out pdb.blast -dbtype prot
my $BLASTP_EXE="blastp"; #"/home/luholm/ncbi-blast-2.8.1+/bin\/blastp";
my $BLAST_NUM_THREADS=32;
# temporary BLAST input/output
my $tmpblastin="$$.fasta";
my $tmpblastout="$$.blast";
# dali programs
my $MPIRUN_EXE= auxiliary_lib::which('mpirun'); #"/usr/lib64/openmpi/bin/mpirun "; #-output-filename x  ";
my $MPIDALI_BIN=$FindBin::Bin;
my $MPICOMPARE_EXE="$MPIDALI_BIN/mpicompare";
my $SERIALCOMPARE_EXE="$MPIDALI_BIN/serialcompare";
# dali params
my $DALIDATDIR_1="./DAT/";
my $DALIDATDIR_2="./DAT/";
my $JOBID="test";
my $TITLE="test";
my $NPARA=1;
my $DOWOLF=1;
my $DOPARSI=1;
my $matrix;
my $hierarchical;
my $walk;
my $list1file;
my $list2file;
my $cd1;
my $cd2;
my $PDBID1='mol1';
my $PDBID2='mol2';
my $PDBFILE1;
my $PDBFILE2;
my $ONEWAY=0;
my $H=100;
my $HMAX=200;
my $KMAX=2000;
my $MAX_HITS=10000;
my $MAX_DALICON=10000;
my $PDB_list;
my $PDB25_list;
my $PDB90_list;
my $OUTFMT='summary';
my $OUTHTML;
my $help;
my $clean;
my $MAXNPARA=99;
my @tmpfiles=('input2.list','soaplist1', 'wolflist1','list0', 'list1', 'list2', 'wolf_output', 'list12', 'filter_input', 'dalicon_input', 'dp_input', 'ww');

my $USAGE= <<"EOB";
Dali is a program for pairwise structural alignment and database searching. Structural alignments with a Z-score above 2 are reported. Outputs are xxxxX.txt for each query structure xxxxX. 

* Pairwise alignment: 
	$0 --query <list-of-query-structure-identifiers> --db <list-of-target-structure-identifiers> [OPTIONS]

* All-against-all comparison: 
	$0 --query <list-of-query-structure-identifiers> --matrix [OPTIONS]

* Database search:
	$0 --query <list-of-query-structure-identifiers> --db <list-of-pdb-structure-identifiers> ( --hierarchical | --walk ) --repset <list-of-pdb-subset-structure-identifiers> [OPTIONS] 

	OPTIONS:
        --dat1 <path>             path to directory containing query data [default: $DALIDATDIR_1]
        --dat2 <path>             path to directory containing target data [default: $DALIDATDIR_2]
        --np <integer>            number of processes if using openmpi [default: $NPARA]
	-h, --help		  long help text, more options

All structures must be imported beforehand using import.pl. 

Citation: Holm L (2019) Benchmarking fold detection by DaliLite.v5 (unpublished)

EOB

my $LONGHELP= <<"EOB";
Dali is a program for pairwise structural alignment and database searching. Structural alignments with a Z-score above 2 are reported. Outputs are xxxxX.txt for each query structure xxxxX. 

USAGE: $0 [ BASIC-OPTIONS] [MPI-OPTIONS] \\
	( --cd1 <xxxxX> |  --pdbfile1 <first.pdb> [ --pdbid1 <mol1> ] | --query <query.list> ) \\
	( --matrix | --cd2 <yyyyY> | --pdbfile2 <second.pdb> [ --pdbid2 <mol2> ] | --db <target.list> \\
	  [ ( --hierarchical | --walk [WALK-OPTIONS] ) --repset <pdb25.list> [BLAST-OPTIONS] ] )

        --cd1 <xxxxX>             query structure identifier
        --pdbfile1 <filename>     query structure in PDB format
        --pdbid1 <xxxx>           four-letter query structure identifier (chain identifier will be appended automatically) [default: $PDBID1]
        --query <filename>        list of query structure identifiers
        --matrix                  all-against-all comparison. Generates additional outputs called 'ordered' (similarity matrix) and 'newick_unrooted' (dendrogram).
        --cd2 <xxxxX>             target structure identifier
        --pdbfile2 <filename>     target structure in PDB format
        --pdbid2 <yyyyy>          four-letter target structure identifier (chain identifier will be appended automatically) [default: $PDBID2]
	--db <filename>	  	  list of target structure identifiers 
	--hierarchical		  hierarchical structure database search 
	--walk			  knowledge-based structure database search
        --repset <filename>       list of structure identifiers of non-redundant subset of PDB

	BASIC-OPTIONS:
        --dat1 <path>             path to directory containing query data [default: $DALIDATDIR_1]
        --dat2 <path>             path to directory containing target data [default: $DALIDATDIR_2]
        --oneway                  asymmetric structure comparison (A,B) only [default: symmetric (A,B) and (B,A)]
	--title <string>          written to output [default: $TITLE]
	--outfmt <string>	  result blocks to output: summary,alignments,equivalences,transrot [default: $OUTFMT]
	--clean			  remove temporary files

	MPI-OPTIONS:
        --np <integer>            number of processes if using openmpi (between 1 and $MAXNPARA) [default: $NPARA]
`	--MPIRUN_EXE <string>	  location of mpirun executable [default: $MPIRUN_EXE]

	BLAST-OPTIONS:
	--HMAX <integer>          number of top scoring representatives to send to final BLAST [default: $HMAX]
	--KMAX <integer>          number of final BLAST hits to align structurally [default: $KMAX]
	--BLAST_DB <string>	  location of Blast database [default: $BLAST_DB]
	--BLASTP_EXE <string>	  location of Blast executable [default: $BLASTP_EXE]
	--BLAST_NUM_THREADS <integer>	number of threads when runnign Blast [default: $BLAST_NUM_THREADS]

	WALK-OPTIONS:
        --targetset <pdb25.list>  used with H to limit the radius of the walk [default: same as --repset]
        --H <integer>             walk radius is Z-score of Hth hit in the target set [default: $H]
        --MAX_HITS <integer>      number of hits returned from walk [default: $MAX_HITS]
        --MAX_DALICON <integer>   max number of comparisons performed during walk [default: $MAX_DALICON]


* Pairwise alignment: 
	$0 ( --cd1 <xxxxX> | --query <query.list> ) ( --cd2 <yyyyY> | --db <target.list> ) [BASIC-OPTIONS]
	bin/dali.pl --pdbfile1 first.pdb [ --pdbid1 mol1 ] --pdbfile2 second.pdb [ --pdbid2 mol2 ] [BASIC-OPTIONS]

* All-against-all comparison: 
	$0 --matrix --query <query.list> [BASIC-OPTIONS]

* Database searches:
	$0 --hierarchical --repset <pdb25.list> ( --cd1 <xxxxX> | --query <query.list> ) --db <pdb.list> [BASIC-OPTIONS] [BLAST-OPTIONS]
	$0 --walk --repset <pdb25.list> ( --cd1 <xxxxX> | --query <query.list> ) --db <pdb.list> [BASIC-OPTIONS] [BLAST-OPTIONS] [WALK-OPTIONS]

All structures must be imported beforehand using import.pl (unless using --pdbfile1 and --pdbfile2 options).  For database searches, PDB-BLAST must be installed locally. The knowledge base is accessed remotely over the Internet.  Non-redundant subsets of the PDB can be generated using, for example, CD-HIT. 

Citation: Holm L (2019) Benchmarking fold detection by DaliLite.v5 (unpublished)

EOB

GetOptions(
	'h|help' => \$help,
        'cd1=s' => \$cd1,
	'cd2=s' => \$cd2,
        'query=s' => \$list1file,
        'db=s' => \$list2file,
	'PDBFILE1=s' => \$PDBFILE1,
	'PDBFILE2=s' => \$PDBFILE2,
	'PDBID1=s' => \$PDBID1,
	'PDBID2=s' => \$PDBID2,
	'outfmt=s' => \$OUTFMT,
	'html' => \$OUTHTML,
        'matrix' => \$matrix,
	'hierarchical' => \$hierarchical,
	'walk' => \$walk,
        'repset=s' => \$PDB25_list,
        'targetset=s' => \$PDB90_list,
        'dat1=s' => \$DALIDATDIR_1,
        'dat2=s' => \$DALIDATDIR_2,
        'np=i' => \$NPARA,
	'MPIRUN_EXE=s' => \$MPIRUN_EXE,
	'BLASTP_EXE=s' => \$BLASTP_EXE,
	'BLAST_DB=s' => \$BLAST_DB,
	'BLAST_NUM_THREADS=i' => \$BLAST_NUM_THREADS,
        'jobid=s' => \$JOBID,
        'title=s' => \$TITLE,
        'wolf=i' => \$DOWOLF,
        'parsi=i' => \$DOPARSI,
        'oneway' => \$ONEWAY,
	'clean' => \$clean,
        'KMAX=i' => \$KMAX,
        'MAX_HITS=i' => \$MAX_HITS,
        'MAX_DALICON=i' => \$MAX_DALICON) or die $USAGE;

# print long help
if($help) { die $LONGHELP; }
my $BLAST_EXE="$BLASTP_EXE -db $BLAST_DB -num_threads $BLAST_NUM_THREADS -outfmt '6 qseqid sseqid pident length qlen slen bitscore evalue' -max_target_seqs 2000 -evalue 0.01 ";

# check mandatory parameters
if (!$list1file && !$cd1 && !$PDBFILE1) { die $USAGE; }
if (!$list2file && !$cd2 && !$PDBFILE2 && !$matrix) { die $USAGE; }
if ($matrix && !$list1file) { die $USAGE; }
if ( ($hierarchical || $walk) && !($PDB25_list && $list2file) ) { die $USAGE; } 
if(!$PDB90_list) { $PDB90_list=$PDB25_list; }
if($NPARA > $MAXNPARA) { $NPARA=$MAXNPARA; warn "# NPARA reduced to $MAXNPARA\n"; }
# add / to dalidatpaths if missing
if($DALIDATDIR_1 !~ /\/$/) { $DALIDATDIR_1.='/'; }
if($DALIDATDIR_2 !~ /\/$/) { $DALIDATDIR_2.='/'; }
# select serialcompare or mpicompare
my $COMPARE_EXE=$SERIALCOMPARE_EXE;
if($NPARA>1) { $COMPARE_EXE="$MPIRUN_EXE --np $NPARA $MPICOMPARE_EXE"; }

# hierarchical parameters
my $qcovcutoff=0.6;
my $nbest=1;
my $MINSSE=3;
my $MAXGAP=30; # merge short gaps in tuples2fasta
my %pdb90; # count top hits


# querylist
my @querylist;
if($cd1) { push(@querylist,$cd1); }
if($list1file) {
        open(IN,"<$list1file");
        while(<IN>) {
                my($cd1)=/^\s*(\w{5})/;
                push(@querylist,$cd1);
        }
        close(IN);
}
if($PDBFILE1) {
	my $tmpfile="$$.tmp";
        foreach my $cd1 (&import_one($PDBFILE1,$PDBID1,$DALIDATDIR_1,$tmpfile,1)) {
        	push(@querylist,$cd1);
		my $dccpfile="$cd1\.dccp";
		if($PDBID1 eq 'mol1' && -e $dccpfile) { die "Alignment file $dccpfile already exists - delete or rename it!\n"; }
	}
}

# target list
if($cd2) { system("echo $cd2 > input2.list"); $list2file="input2.list"; }
if($PDBFILE2) {
	my $tmpfile="$$.tmp";
	$list2file="input2.list";
	open(LIST,">$list2file");
        foreach my $cd2 (&import_one($PDBFILE2,$PDBID2,$DALIDATDIR_2,$tmpfile,1)) {
                print LIST "$cd2\n";
        }
	close(LIST);
}

# memorize targetset
if($PDB90_list) {
	open(IN,"<$PDB90_list");
	while(<IN>) {
	        my($cd1)=/^(\w+)/;
	        $pdb90{$cd1}=1;
	}
	close(IN);
}


# lock directory
&lock();

# main loop
if($matrix) { # all against all structure comparisons: generates extra outputs
	warn "# matrix option detected\n";
        my $soaplist1='soaplist1';
        my $wolflist1='wolflist1';
        my(@list)=&split_list1($list1file,$DALIDATDIR_1,$soaplist1,$wolflist1);
        &matrix($soaplist1,$wolflist1,$list1file,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_1,$JOBID,$TITLE,$OUTFMT,$OUTHTML);
} else {
	foreach $cd1 (@querylist) {
		if($hierarchical) {
	        	&do_one_hierarchical($cd1);
		} elsif($walk) {
			&do_one_walk($cd1);
		} else {
			&do_one_pairwise($cd1);
		}
	        # output HTML x PDB
		warn "# $cd1 generate_FSSP\n";
	        &generate_FSSP($cd1,$list2file,"$cd1\.html","$cd1\.txt",$DALIDATDIR_1,$DALIDATDIR_2,$JOBID,$TITLE,$OUTFMT,$OUTHTML);
	}
}

# clean up
if($clean) { 
	my $cmd=join(" ","rm -f ",@tmpfiles);
	system($cmd);
	# remove *.dccp and *.html
	foreach my $cd1 (@querylist) {
		my $cmd="rm -f $cd1\.dccp $cd1\.html";
		system($cmd);
	}
}

# unlock directory
&unlock();

exit();

###############################################################################

sub do_one_pairwise {
	my($cd1)=@_;
	my($nres,$nsse)=&get_size_nsse($cd1,$DALIDATDIR_1);
	system("echo $cd1 > list1 ; cp $list2file list2"); 
	# run soap if needed
	if($nsse<$MINSSE) {
	        &wolf('SOAP',$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
        } else {
	        if($DOWOLF > 0 ) { &wolf('WOLF',$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0); }
                if($DOPARSI > 0 ) { &parsi($COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0); }
                if($ONEWAY) {
                        # do nothing
                } else { # swap list12 back and forth
                        &swaplist12();
                        if($DOWOLF > 0 ) { &wolf('WOLF',$COMPARE_EXE,$DALIDATDIR_2,$DALIDATDIR_1,2); }
                        if($DOPARSI > 0 ) { &parsi($COMPARE_EXE,$DALIDATDIR_2,$DALIDATDIR_1,2); }
                        &swaplist12();
                 }
	}
}

###############################################################################

sub do_one_hierarchical {
        my($cd1)=@_;
        my $dccpfile="$cd1\.dccp";
        my %blasted;

        # self alignment, initializes cd1.dccp
        system("echo $cd1 > list1");
        my($nres,$nsse)=&get_size_nsse($cd1,$DALIDATDIR_1);
        &self_alignment($cd1,$nres,$DALIDATDIR_1,$COMPARE_EXE);
        # set walk cutoff to sqrt(Zself)
        my($zcutoff)=&set_zcutoff($cd1); # sqrt(Zself)
        if($zcutoff<2) { $zcutoff=2.0; }
        warn "Initial $cd1 zcutoff is $zcutoff\n";

        # always run WOLF x PDB25
        system("echo $cd1 > list1 ; cp $PDB25_list list2");
        if($nsse>=$MINSSE) { &wolf('WOLF',$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0); }
        &generate_PARSI_entrypoints($cd1,$nsse,$PDB25_list,$COMPARE_EXE);

        # final BLAST to fill holes, new chains only
        my ($m)=&generate_BLAST_list2($dccpfile,$zcutoff-2.0, $HMAX,$KMAX,$DALIDATDIR_2,1); # exclude hits already in dccpfile
        # remove already blasted from list2
        my @tmp;
        my $n=0;
        open(IN,"< list2"); while (<IN>) { next if(!/\w/); my($cd2)=/(\w+)/; next if($blasted{$cd2}); push(@tmp,$cd2); } close(IN);
        open(OUT,"> list2"); foreach my $cd2 (@tmp) { print OUT "$cd2\n"; $n++; } close(OUT);
        warn "# got $n / $m hits from pdb-BLAST\n";
        &compare_bothways($cd1,$nsse,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2);
}

sub get_sequence {
        my($cd1,$datfile)=@_;
        $_=`grep sequ $datfile`;
        my($seq)=/\"(\w+)/;
        return($seq);
}

###############################################################################

sub do_one_walk {
        my($cd1)=@_;
        my $dccpfile="$cd1\.dccp";
        my %blasted;

        # self alignment, initializes cd1.dccp
        system("echo $cd1 > list1");
        my($nres,$nsse)=&get_size_nsse($cd1,$DALIDATDIR_1);
        &self_alignment($cd1,$nres,$DALIDATDIR_1,$COMPARE_EXE);
        # set walk cutoff to sqrt(Zself)
        my($zcutoff)=&set_zcutoff($cd1); # sqrt(Zself)
        if($zcutoff<2) { $zcutoff=2.0; }
        warn "Initial $cd1 zcutoff is $zcutoff\n";

        my(@list2)=&test_BLAST($cd1,$nsse,$nres,$zcutoff);
        my $n=$#list2;
        warn "# got $n hits from query-BLAST\n";
        foreach my $cd2 (@list2) { $blasted{$cd2}=1; }

        # always run WOLF x PDB25
        system("echo $cd1 > list1 ; cp $PDB25_list list2");
        if($nsse>=$MINSSE) { &wolf('WOLF',$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0); }
        # update cutoff
        my(@x)=&adjust_zcutoff($cd1,$zcutoff,\%pdb90,$H); # returns rank of last PDB25 hit
        $n=$x[0];
        $zcutoff=$x[1];
        # PARSI/SOAP, if no strong hits
        if($n<2) {
                warn "# Running PARSI/SOAP\n";
                &generate_PARSI_entrypoints($cd1,$nsse,$PDB25_list,$COMPARE_EXE);
                # update cutoff
                my(@x)=&adjust_zcutoff($cd1,$zcutoff,\%pdb90,$H); # adjusts global zcutoff
                $zcutoff=$x[1];
        }

        #  walk from entrypoints with transitive alignments
        &walk($cd1,$dccpfile,$zcutoff,$DALIDATDIR_1,$DALIDATDIR_2,$COMPARE_EXE,$MAX_HITS,$MAX_DALICON,\%pdb90,$H);

        # update cutoff
        (@x)=&adjust_zcutoff($cd1,$zcutoff,\%pdb90,$H); # returns rank of last PDB25 hit
        $zcutoff=$x[1];

        # final BLAST to fill holes, new chains only
        my ($m)=&generate_BLAST_list2_repset($dccpfile,$zcutoff-2.0, $HMAX,$KMAX,$DALIDATDIR_2,1,\%pdb90); # exclude hits already in dccpfile, BLAST reps only
        # remove already blasted from list2
        my @tmp;
        $n=0;
        open(IN,"< list2"); while (<IN>) { next if(!/\w/); my($cd2)=/(\w+)/; next if($blasted{$cd2}); push(@tmp,$cd2); } close(IN);
        open(OUT,"> list2"); foreach my $cd2 (@tmp) { print OUT "$cd2\n"; $n++; } close(OUT);
        warn "# got $n / $m hits from pdb-BLAST\n";
        if($nsse<$MINSSE) {
                &wolf("SOAP",$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
        } else {
                system "echo $cd1 > list1";
                &wolf("WOLF",$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
                &parsi($COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
                if($ONEWAY) {
                        # skip
                } else {
                        &swaplist12();
                        &wolf("WOLF",$COMPARE_EXE,$DALIDATDIR_2,$DALIDATDIR_1,2);
                        &parsi($COMPARE_EXE,$DALIDATDIR_2,$DALIDATDIR_1,2);
                        &swaplist12();
                }
        }

        # what zcutoff was used?
        warn "Final $cd1 zcutoff is $zcutoff\n";
}

sub test_BLAST {
        my($cd1,$nsse,$nres,$zcutoff)=@_;
        # case: BLAST against PDB
        my $method="WOLF"; if($nsse<$MINSSE) { $method="SOAP"; }
        my $datfile="$DALIDATDIR_1\/$cd1\.dat";
        my $seq=&get_sequence($cd1,$datfile);
        my @list2;
        foreach(&run_BLAST(&dat2fasta($DALIDATDIR_1,$cd1))) {
                my($x,$y)=split(/\s+/);
                push(@list2,$y);
        }
        my $nlist2=1+$#list2;
        warn "# run_BLAST returned $nlist2:\n@list2\n";
        if($nlist2<1) { return; } # no hits/
        # wolf/parsi against BLAST set
        open(OUT,">list2");
        foreach my $cd2 (@list2) { print OUT "$cd2\n"; }
        close(OUT);
        # cd1 x list2
        &wolf($method,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
        &parsi($COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
        # test ALL BLAST hits in PDB
        if(&check_strong_match($cd1,$nres,$zcutoff,$qcovcutoff)) { # require zscore > sqrt(self)+2, qcov>0#
                warn "# $cd1 has homolog in PDB\n";
                return(@list2);
        }
        # list2 x cd1
        &swaplist12();
        &wolf($method,$COMPARE_EXE,$DALIDATDIR_2,$DALIDATDIR_1,2);
        &parsi($COMPARE_EXE,$DALIDATDIR_2,$DALIDATDIR_1,2);
        &swaplist12();
        # test ALL BLAST Shits in PDB
        if(&check_strong_match($cd1,$nres,$zcutoff,$qcovcutoff)) { # require zscore > sqrt(self)+2, qcov>0#
                warn "# $cd1 has homolog in PDB\n";
                return(@list2);
        }
        return(@list2);
}

sub generate_PARSI_entrypoints {
        my($cd1,$nsse,$PDB90_list,$COMPARE_EXE)=@_;
        # copy cd1 as list1, PDB90 as list2
        my $cmd="echo $cd1 > list1 ; cp $PDB90_list list2";
        warn "$cmd\n";
        system($cmd);
        if($nsse<$MINSSE) {
                &wolf('SOAP',$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
        } else {
                &parsi($COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
                if($ONEWAY) {
                        # skip
                } else {
                        # swap list1,list2
                        &swaplist12(); # swap reverse
                        &parsi($COMPARE_EXE,$DALIDATDIR_2,$DALIDATDIR_1,2);
                        &swaplist12(); # swap back
                }
        }
}

sub check_strong_match { # return TRUE if z>zcutoff and qcov>qcovcutoff
        my($cd1,$nres,$zcut,$qcovcutoff)=@_;
        my $dccpfile="$cd1\.dccp";
        open(IN,"<$dccpfile");
        while(<IN>) {
                # DCCP   1   3161.4 0.0 197    42.2         100      1                1colA 1colA
                if(/DCCP.{17}\s*(\d+)\s+(\S+).*\s+(\d+)\s+(\w+)\s+(\w+)\s*$/) {
                        next if($4 eq $5); # self hit
                        my $lali=$1;
                        my $z=$2;
                        my $pide=$3;
# old server:           if($z>20 or $z>$nres/10 or $pide>20) {close(IN); return(1); } # TRUE: strong match found
                        if($1 > $qcovcutoff*$nres and $z > $zcut) {close(IN); return(1); } # TRUE: strong match found
                }
        }
        close(IN);
        return(0); # FALSE: no strong match
}


####################################################################
# run BLAST of FASTA file against PDB                              #
####################################################################

sub run_BLAST_files {
        # call as &run_BLAST2file($infile,$outfile,@fasta)
        my($infile)=shift;
        my($outfile)=shift;
        if($#_<0) { return; }
        open(OUT,">$infile");
        foreach(@_) { print OUT $_; }
        close(OUT);
        my $cmd="$BLAST_EXE < $infile > $outfile";
        #warn "$cmd\n";
        system($cmd);
        # extract all hits, in order of  first seen
        my @list3;
        my %seen;
        open(IN,"<$outfile");
        while(<IN>) {
                my(@data)=split(/\s+/);
                my($cd3)=$data[1];
                next if($seen{$cd3});
                push(@list3,$cd3);
                $seen{$cd3}=1;
        }
        close(IN);
        return(@list3);
}

sub generate_BLAST_list2 {
        my($dccpfile,$zcutoff,$H,$K,$DALIDATDIR_2,$includedccp)=@_;
        my(@tmp)=&get_zscore_tuples($dccpfile,$zcutoff,0,$H);
        # hash already done
        my %seen;
        foreach(@tmp) {
                my($cd2)=split(/\s/);
                next if($includedccp>0); # want to include self hits (external dccpfile)
                $seen{$cd2}=1;
        }
        my(@fasta)=&tuples2fasta($DALIDATDIR_2, 50, @tmp);
        my(@list3)=&run_BLAST_files($tmpblastin,$tmpblastout,@fasta);
        open(OUT,">list2");
        # accept top K BLAST hits
        my $n=0;
	my $k=$#list3;
	if($k>$K) { $k=$K; }
        foreach my $cd3 (@list3[0..$k]) {
                next if(defined($seen{$cd3}));
                print OUT "$cd3\n";
                $n++;
        }
        close(OUT);
        return($n);
}

sub generate_BLAST_list2_repset {
        my($dccpfile,$zcutoff,$H,$K,$DALIDATDIR_2,$includedccp,$REPSET)=@_;
        my(@tmp)=&get_zscore_tuples($dccpfile,$zcutoff,0,$H);
        my %seen;
        my @tmprep; # BLAST representatives only
        foreach(@tmp) {
                my $line=$_;
                my($cd2)=split(/\s/);
                if($REPSET->{$cd2}) { push(@tmprep,$line); }
                next if($includedccp>0); # want to include self hits (external dccpfile)
                $seen{$cd2}=1;
        }
        my(@fasta)=&tuples2fasta($DALIDATDIR_2, 50, @tmprep);
        my(@list3)=&run_BLAST_files($tmpblastin,$tmpblastout,@fasta);
        open(OUT,">list2");
        my $n=0;
	my $k=$#list3; 
	if($k>$K) { $k=$K; }
        foreach my $cd3 (@list3[0..$k]) {
		next if(defined($seen{$cd3}));
		#next if($seen{$cd3}>0);
                print OUT "$cd3\n";
                $n++;
        }
        close(OUT);
        return($n);
}

sub run_BLAST {
        # call as &run_BLAST(@fastalines);
        warn "# run_BLAST $BLAST_EXE\n";
        my $pid = open3(\*CHLD_IN, \*CHLD_OUT, \*CHLD_ERR, "$BLAST_EXE") or return("open3() failed $!");
        foreach (@_) { print CHLD_IN $_; }
        close CHLD_IN;
        my(@lines) = <CHLD_OUT>; # list of cd2's
        close CHLD_OUT;
        waitpid($pid,0); # wait for child process to finish
        # error exit
        my $status=($?>> 8);
        if($status != 0) { warn "BLAST job failed\n"; warn <CHLD_ERR>,"\n"; }
        close CHLD_ERR;
        return(@lines);
}
####################################################################
# convert cd-list to FASTA                                         #
####################################################################
sub dat2fasta {
        # return FASTA for cd-list
        my($DALIDATDIR,@list)=@_;
        warn "# dat2fasta: $DALIDATDIR $#list inputs\n";
        my $seq='';
        my @lines;
        my %seen;
        foreach (@list) {
                chomp;
                my($id)=/(\w+)/;
                # get description
                my $file="$DALIDATDIR\/$id\.dat";
#       -compnd   " MOLECULE: AVIAN PANCREATIC POLYPEPTIDE;                              "
                my $cmd="egrep '^-compnd|^-sequence' $file";
                #warn "# id: $id cmd: $cmd\n";
                my(@tmp)=`$cmd`;
                foreach(@tmp) {
                        chomp;
                        if(/^\-compnd/) {
                                s/^\S+\s+//;
                                s/\"//g;
                                s/MOLECULE://;
                                s/\;\s+$//;
				if(!defined($seen{$id})) {
                                #if($seen{$id}<1) {
                                        push(@lines,"\>$id $_\n$seq\n");
                                        $seen{$id}=1;
                                }
                                $seq='';
                        } elsif(/^\-sequence/) {
                                s/^\S+\s+//;
                                s/\"//g;
                                $seq=$_;
                        }
                }

        }
        #warn "# dat2fasta returns $#lines\n";
        return(@lines);
}
