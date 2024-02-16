#use lib qw(/home/luholm/DaliLite.v5/bin);
use strict;
use FindBin qw($Bin);
use upgma;

# programs
my $MPIDALI_BIN=$FindBin::Bin;
my $FILTER95FITZ_EXE="$MPIDALI_BIN/filter95fitz";
my $PIPE96_EXE="perl $MPIDALI_BIN/pipe96.pl";
my $PIPE96FREE_EXE="perl $MPIDALI_BIN/pipe96-free.pl";
my $PIPE_EXE="$MPIDALI_BIN/pipe";
my $SORTDCCP_EXE="perl $MPIDALI_BIN/sortdccp.pl";
my $PIPEDCCP_EXE="$MPIDALI_BIN/pipedccp";
my $SELFDCCP_EXE="$MPIDALI_BIN/selfdccp";
my $FSSPFILTER_EXE="perl $MPIDALI_BIN/fsspfilter.pl";
my $FSSP_EXE="$MPIDALI_BIN/fssp";
my $FSSPSELECT_EXE="perl $MPIDALI_BIN/fsspselect.pl";
my $HTML_EXE="perl $MPIDALI_BIN/htmljs.pl";
my $DSSP_EXE="$MPIDALI_BIN/dsspcmbi";
my $PUU_EXE="$MPIDALI_BIN/puu";
my $PUUTOS_EXE="$MPIDALI_BIN/puutos";# parameters
# parameters
my $max_filename=80; # hardcoded in Fortran programs
my $MINLEN=29;
my $doubleflag=0; # override in matrix; 0 = save cd1 only; 1 = save cd1 & cd2; 2 = save cd2 only
my $POLLINTERVAL=1;
my $TRIPLET_THRESHOLD=0.9;
my $TRIPLET_ITERATIONS=10;
my $CONVERGENCE_THRESHOLD=1.05;
my $PIPENBEST=2; # 1 # send suboptimal parsi alignments to dalicon
my $zcut=2.0;
my $nbest=1;
my $MINSSE=3;
# global similarity matrix
my %sola;
my %z;
# lock work directory
my $LOCKFILE="dali.lock"; # this is generated in PWD. Needs to be removed if it exists before serverparser runs
my $TMP; # temporary file
my $DEVNULL='> /dev/null';

sub lock {
	if (-s $LOCKFILE) { die "Directory is locked by $LOCKFILE\n\tthere may be another DALI process running in this work directory\n\t or, the previous run crashed: remove the dali.lock file \n"; }
	    else { system("echo lock  > $LOCKFILE"); }
}

sub unlock {
	system("rm -f $LOCKFILE");
}

sub FSSP_output { # called from pairwise.pl, matrix.pl
	my($list2file, $DALIDATDIR_1, $DALIDATDIR_2, $JOBID, $TITLE, $OUTFMT, $OUTHTML, @list)=@_;
	my $dbsize=$#list; #list+1-$[;
	foreach my $cd1 (@list) {
	        next if($cd1!~/\w{5}/);
	        my $cmd="$SORTDCCP_EXE < $cd1\.dccp | $FSSPFILTER_EXE $cd1 $zcut $nbest 1 | $FSSP_EXE $cd1 $dbsize $DALIDATDIR_1 $DALIDATDIR_2 | $FSSPSELECT_EXE $list2file $zcut $nbest | $HTML_EXE $DALIDATDIR_1 $DALIDATDIR_2 $JOBID \"$TITLE\" \"$OUTFMT\" $OUTHTML > $cd1\.html 2> $cd1\.txt";
	        warn "$cmd\n";
	        system($cmd);
	}
}

sub generate_FSSP { # called from search.pl
        my($cd1,$pdblist,$out_html,$out_txt,$DALIDATDIR_1, $DALIDATDIR_2, $JOBID, $TITLE, $OUTFMT, $OUTHTML)=@_;
        $_=`wc -l $pdblist`;
        my $dbsize=/(\d+)/;
        my $cmd="$SORTDCCP_EXE < $cd1\.dccp | $FSSPFILTER_EXE $cd1 $zcut $nbest 1 | $FSSP_EXE $cd1 $dbsize $DALIDATDIR_1 $DALIDATDIR_2 | $FSSPSELECT_EXE $pdblist $zcut $nbest | $HTML_EXE $DALIDATDIR_1 $DALIDATDIR_2 $JOBID \"$TITLE\" \"$OUTFMT\" $OUTHTML > $out_html 2> $out_txt";
        warn "$cmd\n";
        system($cmd);
}

sub matrix {
        my($soaplist1,$wolflist1,$listfile,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$JOBID,$TITLE,$OUTFMT,$OUTHTML)=@_;
        # save symmetric DCCP
        my $doubleflag=1;
        # create list0 = list1 = list2; remove annotations
        my(@list);
        open(IN,"<$listfile") || die "Can't open $listfile\n";
        open(OUT,">list0");
        while(<IN>) {
                my($cd)=/^\s*(\S+)/;
                print OUT "$cd\n";
                push(@list,$cd);
		warn "# save $cd in list0\n";
        }
        close(IN);
        close(OUT);
        my $cmd="cp list0 list2; cp list0 list1";
	warn "$cmd\n";
        system($cmd);
        # self comparisons
        &selfdccp('list0',$DALIDATDIR_1,$DALIDATDIR_2,$COMPARE_EXE);
        print "# self zscoresum is ",zscoresum(@list),"\n";
        # do we need to run SOAP?
        if(-s $soaplist1 > 1) {
                system("cp $soaplist1 list1");
                &wolf('SOAP',$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag);
        }
        # do we need to run WOLF?
        if(-s $wolflist1 > 1) {
                system("cp $wolflist1 list1");
                &wolf('WOLF',$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag);
        }
        # run WOLF
        my $summa=zscoresum(@list);
        print "# wolf/soap zscoresum is $summa\n";
        system("rm -f fort.*");
        # run PARSI
        if(-s $wolflist1 > 1) {
                &parsi($COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag);
                my $summa=&zscoresum(@list);
                print "# parsi zscoresum is $summa\n";
                system("rm -f fort.*");
        }
        # run triplets
        my $oldsum=$summa;
        my $iter=0;
        my $lfitz="T";
        while ($iter<$TRIPLET_ITERATIONS) {
                $iter++;
                if($lfitz eq "F") { $lfitz="T"; } else { $lfitz="F"; }
                &triplets('dalicon_input',$lfitz,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag,@list);
                my $summa=&zscoresum(@list);
                print "# triplet $lfitz $iter zscoresum is $summa\n";
                last if($summa-$oldsum<$CONVERGENCE_THRESHOLD); # exit loop
                $oldsum=$summa;
                system("rm -f fort.*");
        }
        # output similarity matrix
        &output_simmatrix("ww",@list);
        # generate oredered matrix, tree and heatmap
        &upgma($listfile,"ww","ordered","tree.xml","tree.html","heatmap.html");
        # generate fssp-html files for each structure
	&FSSP_output('list0',$DALIDATDIR_1,$DALIDATDIR_2,$JOBID,$TITLE,$OUTFMT,$OUTHTML,@list);
        # clean up
        $cmd="rm -f $TMP fort.*";
        warn "$cmd\n";
        system($cmd);
}

sub swaplist12 {
        system("mv list1 e ; mv list2 list1; mv e list2");
}

sub zscoresum {
        # calculate sum of DCCP z-scores over all DCCP files in list
        my(@list)=@_;
        my $total=0;
        foreach my $cd1 (@list) {
                my $file="$cd1\.dccp";
                open(IN,"<$file") || next;
                my %best;
                while(<IN>) {
                        if(/DCCP.{21}\s+(\S+).*\s+(\S+)\s+(\S+)\s*$/){
                                my $cd2=$2;
                                if($2 eq $cd1) { $cd2=$3; }
                                if($1 > $best{$cd2}) { $best{$cd2}=$1; }
                        }
                }
                close(IN);
                foreach (keys %best) { $total+=$best{$_}; }
        }
       my $result=sprintf("%10.1f",$total);
       return($result);
}

sub selfdccp {
        my($list1,$DALIDATDIR_1,$DALIDATDIR_2,$COMPARE_EXE)=@_;
	my $TMP="$$.tmp";
	my $cmd="$SELFDCCP_EXE $DALIDATDIR_1 < $list1 | $COMPARE_EXE  $DALIDATDIR_1 $DALIDATDIR_1 DP $DEVNULL ; cat fort.* > $TMP ";
        warn "$cmd\n";
        system($cmd);
        # distribute tmp to cd1.dccp
        open(IN,"< $TMP");
        my %data;
        my $cd1;
        while(<IN>) {
                if(/DCCP.*\s(\S+)$/) { $cd1=$1; }
                $data{$cd1}.=$_;
        }
        close(IN);
        foreach $cd1 (keys %data) {
                open(OUT,">>$cd1\.dccp") || next;
                print OUT $data{$cd1};
                close(OUT);
        }
}

sub parsi {
        my($COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag)=@_;
        # mpicompare outputs to unit 100+rank
        my $cmd="rm -f fort.* ; $COMPARE_EXE $DALIDATDIR_1 $DALIDATDIR_2 PARSI $DEVNULL";
        warn "$cmd\n";
        my $status=system($cmd);
        if($status != 0) { warn "ERROR running PARSI\n"; return; }
	$cmd="cat fort.1[0-9][0-9] > filter_input ; rm -f fort.*; $COMPARE_EXE $DALIDATDIR_1 $DALIDATDIR_2 FILTER95 < filter_input ";
	warn "$cmd\n";
        $status=system($cmd);
        if($status != 0) { warn "ERROR running FILTER95\n"; return; }
	$cmd="sort -nr fort.1[0-9][0-9] | uniq | $PIPE96FREE_EXE 1.0 $PIPENBEST > dalicon_input ; rm -f fort.* ";
        warn "$cmd\n";
        $status=system($cmd);
        if($status != 0) { warn "ERROR running $PIPE96FREE_EXE\n"; return; }
        &dalicon('dalicon_input',"T",$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag);
}

sub parsi12 {
        my($list12,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag)=@_;
        # input is tuples (cd1,cd2) from list12
        my $cmd="$COMPARE_EXE $DALIDATDIR_1 $DALIDATDIR_2 PARSI12 < $list12 $DEVNULL";
        warn "$cmd\n";
        my $status=system($cmd);
        if($status != 0) { warn "ERROR running PARSI12\n"; return; }
        $cmd="cat fort.1[0-9][0-9] > filter_input ; rm -f fort.*; $COMPARE_EXE $DALIDATDIR_1 $DALIDATDIR_2 FILTER95 < filter_input ";
        warn "$cmd\n";
        $status=system($cmd);
        if($status != 0) { warn "ERROR running FILTER95\n"; return; }
        $cmd="sort -nr fort.1[0-9][0-9] | uniq | $PIPE96FREE_EXE 1.0 $PIPENBEST > dalicon_input ; rm -f fort.* ";
        warn "$cmd\n";
        $status=system($cmd);
        if($status != 0) { warn "ERROR running $PIPE96FREE_EXE\n"; return; }
        &dalicon('dalicon_input',"T",$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag);
}

sub compare_bothways {
        my($cd1,$nsse,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2)=@_;
	# create  list2 outside
        if($nsse<$MINSSE) {
        	&wolf("SOAP",$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
        } else {
                system "echo $cd1 > list1";
                &wolf("WOLF",$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
                &parsi($COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
		# reverse 
		&swaplist12();
                &wolf("WOLF",$COMPARE_EXE,$DALIDATDIR_2,$DALIDATDIR_1,2);
                &parsi($COMPARE_EXE,$DALIDATDIR_2,$DALIDATDIR_1,2);
		# restore list2
		&swaplist12();
        }
}

sub compare_bothways_wolfsoap {
        my($cd1,$nsse,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2)=@_;
        if($nsse<$MINSSE) {
                &wolf("SOAP",$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
        } else {
                system "echo $cd1 > list1";
                &wolf("WOLF",$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
		&swaplist12();
                &wolf("WOLF",$COMPARE_EXE,$DALIDATDIR_2,$DALIDATDIR_1,2);
        	&swaplist12();
	}

}

sub wolf {
        my($wolfsoap,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag)=@_;
	warn "# sub wolf($wolfsoap,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag)\n";
        # postprocess WOLF results
        my $cmd="$COMPARE_EXE $DALIDATDIR_1 $DALIDATDIR_2 $wolfsoap $DEVNULL ; cat fort.1[0-9][0-9] > wolf_output ; rm fort.1[0-9][0-9] ; $COMPARE_EXE $DALIDATDIR_1 $DALIDATDIR_2 DP < wolf_output $DEVNULL ; cat fort.1[0-9][0-9] | perl $MPIDALI_BIN\/dccp2dalicon.pl > dalicon_input ; rm -f fort.1[0-9][0-9] ";
	warn "wolf: $cmd\n";
	my $status=system($cmd);
        if($status != 0) { warn "ERROR running WOLF\n"; return; }
        # concatenate fort.1[0-9][0-9] for mpidalicon
        &dalicon('dalicon_input',"T",$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag);
}

sub dalicon {
        my($dalicon_input,$lfitz,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag)=@_;
 	warn "# sub dalicon($dalicon_input,$lfitz,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag)\n";
	$_=$COMPARE_EXE;
	my($NPARA)=/\-\-np (\d+)/;
        my $cmd="$COMPARE_EXE $DALIDATDIR_1 $DALIDATDIR_2 DALICON $lfitz < $dalicon_input $DEVNULL"; # writes to fort.(10+rank)
        warn "# $cmd\n";
        system($cmd);
        # format conversion to proper DCCP
	$cmd="cat fort.1[0-9][0-9] > dp_input ; rm fort.* ; $COMPARE_EXE $DALIDATDIR_1 $DALIDATDIR_2 DP < dp_input $DEVNULL";
        warn "# $cmd\n";
        system($cmd);
        # collect alignments from fort.10*rank; append to cd1.dccp and cd2.dccp
        my %data;
        my $cd1='';
        my $cd2='';
        my $block='';
	my $NPARA=1;
	if($COMPARE_EXE =~ /\-\-np (\d+) /) { $NPARA=$1; }
	my $k=100+$NPARA-1;
	if($k<101) { $k=101; } # serialcompare
        foreach my $inunit (100+1..$k) {
                my $file="fort\.$inunit"; # dp output
                open(IN,"<$file") || next;
                # collect all DCCP blocks in big hash{cd1}
                while(<IN>) {
                        if(/DCCP/) {
                                if($cd1 and $doubleflag<2) { $data{$cd1}.=$block; }
                                if($cd2 and $doubleflag>0) { $data{$cd2}.=$block; }
                                $block='';
                                ($cd1,$cd2)=/(\w+)\s+(\w+)\s*$/;
                        }
                        $block.=$_;
                }
                close(IN);
                # last block
                if($cd1 and $doubleflag<2) { $data{$cd1}.=$block; }
                if($cd2 and $doubleflag>0) { $data{$cd2}.=$block; }
        }
        foreach my $cd1 (keys %data) {
                my $file=$cd1.".dccp";
                open(OUT,">>$file") || warn "Can't open $file\n";
                print OUT $data{$cd1};
		#warn "# dalicon ($doubleflag) outputting $cd1\n$data{$cd1}\n";
                close(OUT);
        }
	# clean up
	system("rm -f fort.*");
}

sub triplets {
        # test all triplets, induce weakest link
        my($outfile,$lfitz,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag,@list)=@_;
        open(OUT,">$outfile") || return; ## "Can't open $outfile\n";
        # read DCCP files, store zscores in hash
        my %block; # alignments
        foreach my $cd1 (@list) {
                my $file="$cd1\.dccp";
                open(IN,"<$file") || next;
                my $i="";
                my $j="";
                my $z=0;
                my $lkeep=1;
                while(<IN>) {
                        if(/DCCP/) {
                                ($z,$i,$j)=/^.{26}\s+(\S+)\s+\S+\s+\S+\s+(\w+)\s+(\w+)/;
                                if($z{$i}{$j}>$z) { $lkeep=0; } else { $lkeep=1; }
                                #warn "# triplet: $lkeep $_\n";
                                if($lkeep==1) { $z{$i}{$j}=$z; $block{$i}{$j}=''; } # overwrite block
                        }
                        if($lkeep==1) { $block{$i}{$j}.=$_; }
                }
                close(IN);
        }
        # symmetrize
        foreach my $cd1 (@list) {
                 foreach my $cd2(@list) {
                        if($z{$cd1}{$cd2}>$z{$cd2}{$cd1}) {
                                $z{$cd2}{$cd1}=$z{$cd1}{$cd2};
                                $block{$cd2}{$cd1}=$block{$cd1}{$cd2};
                        } else {
                                $z{$cd1}{$cd2}=$z{$cd2}{$cd1};
                                $block{$cd1}{$cd2}=$block{$cd2}{$cd1};
                        }
                 }
        }
        # repair all pairs
        open(OUT2,">list12");
        my $n=$#list;
        foreach my $i (0..$n) {
                my $cd1=$list[$i];
                print OUT "$cd1  \n";
                foreach my $j ($i+1..$n) {
                        my $cd2=$list[$j];
                        my $z12=$z{$cd1}{$cd2};
                        my $cdx='';
                        my $bestx=0.0;
                        foreach my $k ($j+1..$n) {
                                my $cd3=$list[$k];
                                my $z13=$z{$cd1}{$cd3};
                                my $z32=$z{$cd3}{$cd2};
                                if($z13 > $bestx && $z32 > $bestx) {
                                        $cdx=$cd3;
                                        $bestx=$z13; if($z32<$z13) { $bestx=$z32; }
                                }
                        }
                        next if($bestx <= $sola{$cd1}{$cd2});
                        $sola{$cd1}{$cd2}=$bestx; # tested pass height
                        next if($bestx < $z12/$TRIPLET_THRESHOLD);
                        next if($cdx eq "");
                        print join("\t","# triplet",$cd1,$cd2,$cdx,$z12,$bestx,$z{$cd1}{$cdx},$z{$cdx}{$cd2}),"\n"; ##,$block{$cd1}{$cdx},$block{$cdx}{$cd2};
                        my(@ali13)=&dccptoali($cd1,$block{$cd1}{$cdx});
                        my(@ali32)=&dccptoali($cdx,$block{$cdx}{$cd2});
                        my @ali12;
                        # init ali12
                        my $n=0;
                        my @pr1;
                        my @pr2;
                        foreach my $i1 (1..$#ali13) {
                                $ali12[$i1]=0; # initialize
                                my $i3=$ali13[$i1];
                                next if($i3<1);
                                my $i2=$ali32[$i3];
                                next if($i2<1);
                                $ali12[$i1]=$i2;
                                $n++;
                                push(@pr1,$i1,$i1);
                                push(@pr2,$i2,$i2);
                        }
                        next if($n<$MINLEN);
                        print OUT "$cd2\* \n$n\n@pr1\n@pr2\n";
                        print OUT2 "$cd1 $cd2\n";
                }
                print OUT "END\n";
        }
        print OUT "END\n";
        close(OUT);
        close(OUT2);
        # refine alignments
        &dalicon("$outfile",$lfitz,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag);
        &parsi12("list12",$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,$doubleflag);
}

sub output_simmatrix {
        # assumes %z has been filled beforehand (in triplets)
        my($outfile,@fulllist)=@_;
	my $BASE=0.1;
	# set baseline similarity to 0.01
	foreach my $cd1 (@fulllist) {
		foreach my $cd2 (@fulllist) {
			if($z{$cd1}{$cd2}<$BASE) { $z{$cd1}{$cd2}=$BASE; }
		}
	}
	# remove zero incidence rows/columns: assumes %z is non-zero
	my(@list);
	foreach my $cd1 (@fulllist) {
		my $summa=0.0;
		foreach my $cd2 (@fulllist) {
			next if($cd2 eq $cd1);
			$summa+=$z{$cd1}{$cd2};
		}
		warn "z-summa of $cd1 is $summa\n";
		if($summa > 0.0) { push(@list,$cd1); } else { warn "# excluded $cd1 from simmatrix\n"; }
	}
        # output similarity matrix
        open(OUT,">$outfile");
        print OUT sprintf "%5i\n", $#list+1-$[,"\n";
        foreach my $cd1 (@list) {
                my(@x)=($cd1);
                foreach my $cd2 (@list) {
                        push(@x,sprintf("%8.1f",$z{$cd1}{$cd2}));
                }
                print OUT join("\t",@x),"\n";
        }
        close(OUT);
}

#______________________________________________________________________________
# Title     : dccptoali
# Usage     :
# Function  :
# Example   :
# Keywords  :
# Options   :
# Author    : holm@ebi.ac.uk, jong@ebi.ac.uk
# Category  :
# Returns   :
# Version   : 1.0
#------------------------------------------------------------------------------
sub dccptoali {
        my ($cd1,$block)=@_;
        my (@l1,@l2,@r1,@r2,@ali,$i,$i1,$i2,$shift,$l);
        my (@lines)=split("\n",$block);

        $_=$lines[$[];
        my ($nblock,$id1,$id2)=/\s+(\d+)\s+(\w+)\s+(\w+)\s{0,2}$/;

# DCCP   1    204.5 3.5  73     1.6           5     11                1amx  1kcw
# alignment
#   5     8   9    12  18    21  31    35  36    43  66    72  93    97 102   114
# 122   131 138   141 142   150
# 918   921 930   933 935   938 945   949 951   958 964   970 972   976 977   989
# 991  10001001  10041008  1016

        my (@ranges)=split(/\s+/,$block);
        if ($id1 eq $cd1) {
            foreach my $i (1..$nblock) {
                ($r2[$i])=pop(@ranges);
                if(length($r2[$i])>4) { $_=$r2[$i]; my ($x,$y)=/^(\d+)(\d{4})/; $r2[$i]=$y; push(@ranges,$x); }
                ($l2[$i])=pop(@ranges);
                if(length($l2[$i])>4) { $_=$l2[$i]; my ($x,$y)=/^(\d+)(\d{4})/; $l2[$i]=$y; push(@ranges,$x); }
            }
            foreach my $i (1..$nblock) {
                ($r1[$i])=pop(@ranges);
                if(length($r1[$i])>4) { $_=$r1[$i]; my ($x,$y)=/^(\d+)(\d{4})/; $r1[$i]=$y; push(@ranges,$x); }
                ($l1[$i])=pop(@ranges);
                if(length($l1[$i])>4) { $_=$l1[$i]; my ($x,$y)=/^(\d+)(\d{4})/; $l1[$i]=$y; push(@ranges,$x); }
            }
        } else {                                        ## swap cd2-cd1
            foreach my $i (1..$nblock) {
                ($r1[$i])=pop(@ranges);
                if(length($r1[$i])>4) { $_=$r1[$i]; my ($x,$y)=/^(\d+)(\d{4})/; $r1[$i]=$y; push(@ranges,$x); }
                ($l1[$i])=pop(@ranges);
                if(length($l1[$i])>4) { $_=$l1[$i]; my ($x,$y)=/^(\d+)(\d{4})/; $l1[$i]=$y; push(@ranges,$x); }
             }
            foreach my $i (1..$nblock) {
                ($r2[$i])=pop(@ranges);
                if(length($r2[$i])>4) { $_=$r2[$i]; my ($x,$y)=/^(\d+)(\d{4})/; $r2[$i]=$y; push(@ranges,$x); }
                ($l2[$i])=pop(@ranges);
                if(length($l2[$i])>4) { $_=$l2[$i]; my ($x,$y)=/^(\d+)(\d{4})/; $l2[$i]=$y; push(@ranges,$x); }
            }
        }
        undef @ali;
        foreach my $i (1..$nblock) {               # invert negative ranges !!!!
            #print "$i -> $l1[$i] .. $r1[$i]   $l2[$i] .. $r2[$i]\n";
            my $i1=$l1[$i];
            my $j1=$r1[$i];
            my $i2=$l2[$i];
            my $shift=$i2-$i1;
            foreach my $k ($i1..$j1) { $ali[$k]=$k+$shift; }
        }
        return (@ali);
}

#______________________________________________________________________________
# Title     : get_type
# Usage     :
# Function  :
# Example   :
# Keywords  :
# Options   :
# Author    : holm@ebi.ac.uk, jong@ebi.ac.uk
# Category  :
# Returns   :
# Version   : 1.0
#------------------------------------------------------------------------------
sub get_type {
        #
        # call as kind=get_type($set1,"argument") or    -- processing arguments
        #         kind=get_type($item,"list")           -- processing list items
        #
        # calls itself recursively until gets a terminal kind
        #
        # terminals are: code | brkfile | list | unknown -- argument context
        #               code | brkfile | unknown        -- list context
        #
        # syntax:       argument -> file | code | unknown
        #               file -> list | brkfile
        #               list -> code | brkfile
        #
        #               file is a text file
        #               code is 5 characters and no text file
        #               brkfile is a file with ATOM records
        #               list is a file that is not a brkfile
        #
        my ($item,$context)=@_;
        my $isfile;

        if (-T $item) {$isfile=1;} else {$isfile=0;};
        if ($context eq "argument") {   # expect file | code
                my $x=length($item);
                if(($x == 4 || $x==5) && (!$isfile)) { $context="code";
                } elsif ($isfile){ $context="file";
                } else { $context="unknown";
                }
        } elsif ($context eq "file") {  # expect list | brkfile
                        if(!$isfile) {$context="unknown";} else
                        {
                          open(IN1,"<$item"); my $x="list"; while(<IN1>)
                          {if (/^ATOM  /) {$x="brkfile"; last;};}; close(IN1);
                          $context=$x;
                        }
        } elsif ($context eq "list") {  # expect code | brkfile
                my $x=length($item);
                if(( $x==5 ) && (!$isfile) ) { $context="code";
                } else {
                        $x="unknown"; open(IN1,"<$item") || warn "Can't open $item\n"; 
			while(<IN1>) 
                        	{if (/^ATOM  /) {$x="brkfile"; last;};}; 
			close(IN1);
                        $context=$x;
                }
        } else { warn "get_type called in unknown context \"$context\""; return; }
        while ($context ne "list" && $context ne "code" &&
                $context ne "brkfile" && $context ne "unknown")
                {$context=&get_type($item,$context);}
        return $context;
}


#______________________________________________________________________________
# Title     : prepare_brkfile
# Usage     :
# Function  : PDB reader, output is internal data file
#             This also generates DSSP.
# Example   :
# Keywords  :
# Options   :
# Author    : holm@ebi.ac.uk, jong@ebi.ac.uk
# Category  :
# Returns   :
# Version   : 1.1
#------------------------------------------------------------------------------
sub prepare_brkfile {       # call as prepare_brkfile($query_pdbfile,$short)
	my ($query_pdbfile,$queryid,$DALIDATDIR)=@_;
	my ($brkheader,$brkcompnd,$brksource,$brkauthor);

# filename length check
if(length($query_pdbfile) > $max_filename) { warn "STOP: file names incl. path must be max $max_filename long\n"; return; }


#----------------------------------------------------------------------------
# 	It is assumed that query-pdbfile is a full-atom entry ! --> CHECK!
#		else return empty chain list + warning
#----------------------------------------------------------------------------

	open(IN,"<$query_pdbfile");
	my $nca=0;
	my $natom=0;
	my @lines;
	while(<IN>) {
		if(/^ATOM  /){ push(@lines,$_); $natom++;if(/ CA  /){$nca++;}}
		last if(/^ENDMDL/);
	} close(IN);
	my $x=$natom / ($nca+1);
	if($x<3){
		print "suspect $query_pdbfile is CA-only: $nca CAs / $natom total atoms\n";
		print @lines;
		return;
	}

#----------------------------------------------------------------------------
# 	always generate DSSP file 
#----------------------------------------------------------------------------

	my $query_dsspfile = $queryid . ".dssp";
	my $cmd="$DSSP_EXE $query_pdbfile $query_dsspfile";
	warn "$cmd\n";
	my $status=system($cmd);
	if($status) { warn "DSSP error\n"; }

#----------------------------------------------------------------------------
# 	chain names arrays
#----------------------------------------------------------------------------

	my @CHAIN;
        my @CHAINLEN;
	open(IN,$query_pdbfile); 
	my $first=0; 
	my $chain_len=0;
	my $chain_num=0;
	my $chain="?";
	while(<IN>) {
		chomp; 
		my $line= $_;
		if(/^ATOM/) {
			my ($atom,$xchain)=/^.{13}([" "\w]{4}).{4}([" "\w])/;
			if ($xchain ne $chain) {
				if($first != 0) { push(@CHAINLEN,$chain_len); }
				$first=1; $chain = $xchain;push(@CHAIN, $chain);
				$chain_num++; $chain_len=0;
			}
			if($atom eq "CA  ") { $chain_len++; }
		}
		last if (/^ENDMDL/); # use only first model of NMR entries
	}
	close(IN);
	#last line
	push(@CHAINLEN,$chain_len);

	$x=$#CHAIN+1; 
	foreach ($[..$#CHAIN){print $CHAIN[$_]," length: ", $CHAINLEN[$_],"\n";}

	# remove shor:t chains
	my @long;
	foreach($[..$#CHAIN){if($CHAINLEN[$_]>$MINLEN){push(@long,$CHAIN[$_]);}}
    	@CHAIN=@long;

        # make non-redundant list (case: Swissmodel had edited blank-chain residues in A chain)
	my %x;
	foreach (@CHAIN) { $x{$_}=1; }
	undef(@CHAIN);
	foreach (sort keys %x) { push(@CHAIN,$_); }
	
	$x=$#CHAIN+1-$[;
        print "# $x chains @CHAIN\n";

	my @shortcode;
	foreach my $i ( $[ .. $#CHAIN) {
		if($CHAIN[$i] ne " ") {
			push(@shortcode,($queryid . $CHAIN[$i]));
		} else {
			push(@shortcode,$queryid);
		}
	}

#----------------------------------------------------------------------------
# 	domain decomposition
#----------------------------------------------------------------------------

	# create puu.default file
	open(OUT,">puu.default");
	print OUT join("\n","#lhetatm","TRUE","#printdepth","-99","#lsplit","TRUE"),"\n";
	close(OUT);

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# It prints what chains are in PDB file
	#__________________________________________
	foreach (my $i=$[; $i< @CHAIN;  $i++) {
        	print "\n CHAIN $CHAIN[$i]\n";
		system("rm -f  domains.puu units.puu subunits.puu");
		my $cmd="$PUU_EXE $shortcode[$i] $query_pdbfile $query_dsspfile";
		my $status=system($cmd);
		if($status) { die "Error in puu: $cmd\n"; }
		$cmd="$PUUTOS_EXE units.puu $query_dsspfile $DALIDATDIR\/";
		$status=system($cmd);
		if($status) { die "Error in puutos: $cmd\n"; }
	}
	my @idlist;
	foreach my $chain (@CHAIN) {   
		print "\n CHAIN : $queryid$chain"; 
		my $x=$queryid;
		if($chain=~/\w/) { $x.=$chain; }
		push(@idlist,$x);
	}
	print "\n\n";

	return(@idlist);
}

#----------------------------------------------------------------------------
# end sub: prepare_brkfile
#----------------------------------------------------------------------------

sub matrix_html {
        my($file,$TITLE,@list1)=@_;
        open(OUT,">$file") || die "Can't open $file\n";
        print OUT<<EOB;
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>$TITLE</title>
  <link rel="stylesheet" href="/css/jquery-ui.css">
  <script src="/js/jquery.js"></script>
  <script src="/js/jquery-ui.js"></script>
  <script src="/js/plotly-latest.min.js"></script>
  <link rel="stylesheet" href="/css/style.css">
  <script>
        \$(function(){
                \$( "#tabs" ).tabs();
                \$("#includeContentHeatmap").load("./heatmap.html");
                \$("#includeContentScatter").load("./scatter.html");
        });
  </script>
  <script type="text/javascript" src="/js/raphael/raphael-min.js" ></script>
  <script type="text/javascript" src="/js/jsphylosvg-min.js"></script>
  <script type="text/javascript" src="/js/yui.js"></script>
  <script type="text/javascript">
        window.onload = function(){
                YUI().use('oop', 'json-stringify', 'io-base', 'event', 'event-delegate', function(Y){
                        var uri = "tree.xml";
                        function complete(id, o, args) {
                                var data = o.responseXML; // Response data.
                                var dataObject = {
                                                        xml: data,
                                                        fileSource: true
                                                };
                                phylocanvas = new Smits.PhyloCanvas(
                                        dataObject,
                                        'svgCanvas',
                                        800, 800
//                                      ,'circular'
                                );
                        };
                        Y.on('io:complete', complete, Y);
                        var request = Y.io(uri); // links to web page
                });
        };
  </script>
</head>
<body>
<h1>Results: $TITLE</h1>
<div id="tabs">
  <ul>
    <li><a href="#tabs-1">Dendrogram</a></li>
    <li><a href="#tabs-2">Heatmap</a></li>
    <li><a href="#tabs-3">Projection</a></li>
    <li><a href="#tabs-4">Summaries</a></li>
    <li><a href="#tabs-5">Download</a></li>
  </ul>
  <div id="tabs-1">
<p>Structural similarity dendrogram. Labels are linked to structural summaries.
The dendrogram is derived by average linkage clustering of the structural similarity matrix (Dali Z-scores).
</p>
    <div id="svgCanvas"> </div>
  </div>
  <div id="tabs-2">
        <div id="includeContentHeatmap"></div>
  </div>
  <div id="tabs-3">
        <div id="includeContentScatter"></div>
  </div>
  <div id="tabs-4">
    <UL>
EOB
	print OUT "<TABLE><TR><TH>Structure</TH><TH>Interactive</TH><TH>Download</TH></TR>\n";
        foreach my $cd1 (@list1) {
                #print OUT "<LI><A HREF=$cd1\.html>$cd1</A></LI>\n";
		print OUT "<TR><TD>$cd1</TD><TD><A HREF=$cd1\.html>html</A></TD><TD><A HREF=$cd1\.txt>txt</A></TD></TR>\n";
        }
	print OUT "</TABLE>\n";
        print OUT<<EOT;
  </div>
  <div id="tabs-5">
    <UL>Data
        <LI><A HREF=ordered>Similarity matrix</A></LI>
        <LI><A HREF=eigen.txt>Eigenvectors</A> from Correspondence Analysis</LI>
        <LI>Newick dendrogram: <A HREF=newick>rooted</A>, <A HREF=newick_unrooted>unrooted</A></LI>
        <LI><A HREF=tree.xml>PhyloXML</A> dendrogram</LI>
    </UL>
  </div>
</div>
<script language="javascript" type="text/javascript">
     $(window).load(function() {
     $('#loading').hide();
  });
</script>
</body>
</html>
EOT
        close(OUT);
}

sub get_size_nsse {
        my($cd1,$WORKDIR)=@_;
        my $nres=0;
        my $nsse=0;
        $_=`head -1 $WORKDIR\/$cd1\.dat` || return(0,0);
        ($nres,$nsse)=/>>>>\s+\S+\s+(\d+)\s+(\d+)/;
        return($nres,$nsse);
}

sub split_list1 {
	# create soaplist1 and wolflist1 in CWD
	my($list1file,$DALIDATDIR,$soaplist1,$wolflist1)=@_;
	open(IN,"<$list1file");
	my @list;
	while(<IN>) {
		my($cd1)=/(\w{5})/; 
		push(@list,$cd1);
	}
	close(IN);
	warn "split_list1: list is @list\n";
	open(SOAP,">$soaplist1");
	open(WOLF,">$wolflist1");
	foreach my $cd1 (@list) {
	        my($nres,$nsse)=&get_size_nsse($cd1,$DALIDATDIR);
	        warn "Structure $cd1 has $nres residues and $nsse secondary structure elements\n";
	        next if($nres<$MINLEN);
	        if($nsse<$MINSSE) { print SOAP "$cd1\n"; } else { print WOLF "$cd1\n"; }
	}
	close(SOAP);
	close(WOLF);
	return(@list);
}

sub sortdccp {
	my($dccpfile,$nkeep)=@_;
	if($nkeep<1) { $nkeep=999999; } # keep all
	my $TRUE=1;
	my $FALSE=0;
	my @array;
	my $first=$TRUE;
	my $block="";
	open(IN,"<$dccpfile") || return(); 
	while(<IN>) {
#DCCP   1   2149.4 0.0 157    34.4         100      1                1mup  1mup
	        if(/^ DCCP\s+(\d+)/) {
	                next if($1!=1); # remove redundant DCCP   2 ... lines
	                # save old block
	                if(!$first) { push(@array,$block); } else { $first=$FALSE; }
	                # initialize new block
	                $block = $_;
	        } else {
	                $block .= $_;   # append to block
	        }
	}
	push(@array,$block);
	close(IN);

	my $i;
	my @datakeys;
	foreach $i ($[..$#array) {
	# DCCP   1   8983.1 0.0 451    66.8         100      1                1j5sA 1j5sA
	        $_=$array[$i]; push(@datakeys,(/^ DCCP.{21}\s+([\-\.\d]+)/));
	}

	#sub bydatakeys { $datakeys[$b] <=> $datakeys[$a]; }
	my @sortarray=@array[sort bydatakeys $[..$#array ];      # sort array on value

	#warn "# datakeys: @datakeys\n";	
	#foreach (0..100) { last if($_>$#array); warn "#sortarray $_: $sortarray[$_]\n"; }

	my @sorted;
	my %seen;
	foreach (@sortarray) {
	        my $block=$_;
	        my(@x)=split(/\n/);
	        my(@y)=split(/\s+/,$x[0]);
	        my $cda=pop(@y);
		my $cdb=pop(@y);
        	my $key="$cda\_$cdb";
       		#warn "key=$key seen=$seen{$key}\n";
        	$seen{$key}++;
        	next if($seen{$key}>$nkeep);
        	foreach (split("\n",$block)) { push(@sorted,$_); }
	}
	return(@sorted);

	sub bydatakeys { $datakeys[$b] <=> $datakeys[$a]; }
}

sub self_alignment {
        my($cd1,$nres,$DALIDATDIR_1,$COMPARE_EXE)=@_;
        open(OUT,">dalicon_input");
        print OUT "$cd1\n$cd1\*\n1\n1 $nres\n1 $nres\nEND\nEND\n";
        close(OUT);
        &dalicon('dalicon_input',"T",$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_1,0);
}

sub get_zself {
        my($cd1)=@_;
        open(IN,"<$cd1\.dccp");
        my $zself=0.0;
        while(<IN>) {
                if(/DCCP.* (\S+)\s+\d+\s+\d+\s+(\w+)\s+(\w+)$/) {
                        next if($2 ne $3);
                        if($1 > $zself) { $zself=$1; }
                }
        }
        return($zself);
}

sub set_zcutoff {
        my($cd1)=@_;
        my $zself=&get_zself($cd1);
        my $zcutoff=sqrt($zself); # sqrt with safety margin
        warn "# cd1=$cd1 zself=$zself zcutoff=$zcutoff\n";
        return($zcutoff);
}

sub tuples2fasta {
	# input @tuples from get_zscore_tuples
	# output FASTA format
	# merge gaps shorted than maxgap
	my($DALIDATDIR,$maxgap,@tuples)=@_;
	my $halfgap=int($maxgap/2);
        my $seq='';
        my @lines;
        my %seen;
        foreach (@tuples) {
                chomp;
                my($id,$z,$tmp)=split(/\t/);
		next if($id eq "");
		next if($seen{$id});
		$seen{$id}=1;
		#warn "# tuples2fasta ($id,$z,$tmp)\n";
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
				# segments substring
				my $subseq='';
				my(@x)=split(/\s+/,$tmp);
				my $where=0;
				while($#x>0) {
					my($from)=shift(@x);
					my($to)=shift(@x);
					if($from-$where<$maxgap) { # merge short gap
						$subseq.=substr($seq,$where,$from-$where); 
					}else {  # long gap
						$subseq.=substr($seq,$where,$halfgap);
						my $x=$from-$halfgap; if($x<0) { $x=0; }
						$subseq.=substr($seq,$x,$halfgap);
					}
					$subseq.=substr($seq,$from,$to-$from); # matched segmetn
					$where=$to;
				}
				#warn "subseq\t$id\t$z\t$tmp\t$subseq\n";
                                push(@lines,"\>$id $_\n$subseq\n");
                                $seq='';
                        } elsif(/^\-sequence/) {
                                s/^\S+\s+//;
                                s/\"//g;
                                $seq=$_;
                        }
                }

        }
        #warn "# tuples2fasta returns $#lines\n";
        return(@lines);	
}

sub get_zscore_tuples {
	my($dccpfile,$zcutoff,$MINLALI,$H)=@_;
	my @tuples;
	my %z; # key = cd2
	my %ranges; # key = cd2
	open(IN,"<$dccpfile");
	while(<IN>) {
		if(/DCCP .* (\S+)\s+\S+\s+(\S+)\s+(\w+)\s+(\w+)\s*$/) {
			#warn $_;
			#warn "z=$1 / $zcutoff cda=$3 $z{$3} cdb=$4 $z{$4}\n";
			next if($1 < $zcutoff); 
			next if($1 <= $z{$3} && $1 <= $z{$4});
			my $nseg=$2;
			my($a,$b,$lali)=&get_DCCP_ranges($nseg);
			next if($lali < $MINLALI);
			#warn join("\t",$1,$nseg,$lali,$MINLALI,$3,$a,$4,$b),"\n";
			if($z{$3}<$1) { $z{$3}=$1; $ranges{$3}=$a; }
			if($z{$4}<$1) { $z{$4}=$1; $ranges{$4}=$b; }
		}
	}
	close(IN);
	my $n=0;
	foreach my $cd2 (sort { $z{$b} <=> $z{$a} } keys %z) {
		#warn "# tuple:",join("\t",$cd2,$z{$cd2},$ranges{$cd2}),"\n";
		push(@tuples,join("\t",$cd2,$z{$cd2},$ranges{$cd2}));
		$n++;
		last if($n>=$H);
	}
	return(@tuples);
}

sub get_DCCP_ranges {
	my($nseg)=@_;
	# assumes <IN> handle is open DCCP-file
	my $i=0;
	my $j=0;
	my @x; # from to from to ...
	# first block
	while(<IN>) {
		#warn "get_DCCP_ranges in1: $_\n";
		chomp;
		next if(/alignment/);
#  143   152 155   166 176   179 180   190 191   197 199   202 219   222 224   262
		s/^\s+//; # remove leading space
		my(@tmp)=split(/\s+/); # from, to, from, to, ...
		while($#tmp >= 0) {
			my($x)=shift(@tmp);
			if($x>=10000) { 
				$_=$x;
				my($a,$b)=/(\d+)(\d{4})/; 
				unshift(@tmp,$b);
				$x=$a;
			}
			push(@x,$x);
		}
		last if($#x+1>=2*$nseg);
	}
	#warn "result1 $#x/$nseg: @x\n";
	# second section
	my @y; 
	while(<IN>) {
		#warn "get_DCCP_ranges in2: $_\n";
		chomp;
		s/^\s+//; # remove leading space
                my(@tmp)=split(/\s+/); # from, to, from, to, ...
                while($#tmp >= 0) {
                        my($x)=shift(@tmp);
                        if($x>=10000) {
                                my($a,$b)=/(\d+)(\d{4})/;
                                unshift(@tmp,$b);
                                $x=$a;
                        }
                        push(@y,$x);
                }
                last if($#y+1>=2*$nseg);
	}	
	#warn "result2 $#y/$nseg: @y\ny),
	my $lali=0;
	my $i=0;
	while($i<$#x) {
		$lali+=$x[$i+1]-$x[$i]+1;
		$i+=2;
	}
	return(join(' ',@x),join(' ',@y),$lali);  
}

sub import_one {
        my($pdbfile,$short,$DALIDATDIR,$tmpfile,$verbose)=@_;
        # uncompress
        if($pdbfile =~ /\.gz$/) {
                if($verbose) { warn "Detected compressed file\n"; }
                system("gzip -cd $pdbfile > $tmpfile");
                $pdbfile=$tmpfile;
        }
        my($x)=&get_type($pdbfile,'list');
        if($x ne 'brkfile') { warn "$pdbfile is $x not a PDB file -skipped\n"; return(); }
        if($short!~/^\w{4}$/) { warn "invalid identifier $short - must be 4 characters [e.g., 1ppt] - skipped\n"; return(); }
        # return list of chain identifiers
        return(&prepare_brkfile($pdbfile,$short,$DALIDATDIR));
}


return(1);
