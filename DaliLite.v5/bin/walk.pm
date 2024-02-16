use strict;
use warnings;
use mpidali;
use FSSP;

# parameters
#my $MAX_HITS=10000; # max size of first shell; walk to twice outputted amount (500) so as to drop fewer hits with medium Z-scores
#my $MAX_DALICON=10000; # max number of comparisons performed
my $MIN_LALI=25;
my $BATCH_SIZE=100; # was 100 # small bs atch-size => largest domain dominates!
my $tmp_dccpfile='tmp.dccp';
my $TMP_DALICON='tmp.dalicon';
my $MAX_ITER=50;
my $LFITZ='T';
my $MAX_RADIUS=1000000; # max size of second shell; stop search when number of induced structures exceeds this
my %lali_induced;
my %lali_observed;
my %inducer;
my %maxtried_lali;
my %map21;
my %map23;
my %z12;
my %z23;
my %iteration; # iteration of improved z12
my %tested; # allow two tests
my $iteration;
my $ndalicon;

sub adjust_zcutoff {
        # call as my(@x)=&adjust_cutoff($cd1,$zcutoff,\%pdb25, $H);
        my($cd1,$initial_zcutoff,$REPSET,$H)=@_;
        my $dccpfile="$cd1\.dccp";
        my $zcutoff=$initial_zcutoff;
        my $i=0;
        foreach(&get_zscore_tuples($dccpfile,$zcutoff,0,100000)) {
                my($cd2,$z)=split(/\t/);
		next if($cd2 !~ /\w{5}/);
		#warn "# adjust_cutoff: i=$i cd1=$cd1 cd2=$cd2 rep=$REPSET->{$cd2} z=$z zcutoff=$zcutoff\n";
                next if(!defined($REPSET->{$cd2}));
                last if($z<$zcutoff);
                $i++;
                if($i==$H) {
                        if($z>$zcutoff) {
                                $zcutoff=$z;
                                warn "# Adjusted zcutoff to $zcutoff\n";
                        }
                        last;
                }
        }
        warn "# $i PDB25 hits above zcutoff=$zcutoff\n";
        return($i,$zcutoff);
}

#call as &walk(cd1 dccpfile zcutoff DALIDATDIR_1 DALIDATDIR_2 NPARA \%pdb25 H)
sub walk {
        my($cd1,$dccpfile,$zcutoff,$DALIDATDIR_1,$DALIDATDIR_2,$COMPARE_EXE,$MAX_HITS,$MAX_DALICON,$REPSET,$H)=@_;
        warn "# This is walk($cd1,$dccpfile,$zcutoff,$DALIDATDIR_1,$DALIDATDIR_2,$COMPARE_EXE,$MAX_HITS,$MAX_DALICON,$REPSET,$H)\n";
# global hashes; key is cd2/3
$iteration=0;
$ndalicon=0;
undef %lali_induced; undef %lali_observed; undef %inducer; undef %maxtried_lali; undef %map21; undef %map23; undef %z12; undef %z23; undef %iteration; undef %tested;

# MAIN LOOP
# terminate if more than MAXHITS
# select cd3 which has largest lali_induced{cd3}-lali_observed{cd3}
#       generate cd1-cd3 alignment
#       add cd3 as cd2
#       induce cd3
my $n_cd3=0;
while(1) {
        $iteration++;
        # pool domain lists
	my($nhit)=&read_dccp($dccpfile,$cd1,$zcutoff); # returns map21{cd2}, lali_observed{cd2}, maxtried_lali{cd2}, $z12{cd2}
        my(@y)=&sortlist_z($zcutoff); #&sortlist_laligain($zcutoff); ###&sortlist_z($zcutoff);
        my $n=$#y;
        if($n>$BATCH_SIZE) { $n=$BATCH_SIZE; }
        last if($n<0);
        my(@keys)=@y[0..$n];

        warn "# calling cd1=$cd1 zcutoff=$zcutoff run_dalicon with list: @keys\n";
        ($n_cd3)=&run_dalicon($cd1,$zcutoff,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,@keys);
        print "# $n_cd3 seed alignments sent to $cd1 dalicon\n";
        $ndalicon+=$n_cd3;
        # adjuct zcutoff
        my(@tmp)=&adjust_zcutoff($cd1,$zcutoff,$REPSET,$H);
        $zcutoff=$tmp[1];
	# test termination
	if($n_cd3<1) {
		print "# Terminating walk at $iteration iterations: converged\n";
		last;
	}
        if($ndalicon>=$MAX_DALICON) {
                print "# Terminating walk at $iteration iterations $ndalicon dalicon runs > MAX_DALICON=$MAX_DALICON\n";
                last;
        }
        if ($iteration>=$MAX_ITER) {
                print "# Terminating walk at $iteration iterations > MAX_ITER=$MAX_ITER\n";
                last;
        }
        warn "# iteration $iteration: $nhit hits returned\n";
        if($nhit>=$MAX_HITS) {
                print "# Terminating walk at $iteration iterations / $nhit hits > MAX_HITS=$MAX_HITS\n";
                last;
        }
}

#print_stack($iteration);
        my %x;
        foreach (keys %z12,keys %z23) { $x{$_}=1; }
        my $cd3;
        my $rank=0;
        foreach $cd3 (keys %z23) {
                $rank++;
                my $cd2=$inducer{$cd3};
		my $a='.'; if(defined($z23{$cd3})) { $a=$z23{$cd3}; }
		my $b='.'; if(defined($lali_induced{$cd3})) { $b=$lali_induced{$cd3};  }
		my $c='.'; if(defined($maxtried_lali{$cd3})) { $c=$maxtried_lali{$cd3}; }
		my $d='.'; if(defined($lali_observed{$cd3})) { $d=$lali_observed{$cd3}; }
		my $e='.'; if(defined($z12{$cd3})) { $d=$z12{$cd3}; }
		my $f='.'; if(defined($iteration{$cd3})) { $f=$iteration{$cd3}; }
		print join("\t",$iteration,$rank,$cd3,$f,$a,$b,$c,$d,$cd2,$e),"\n";
        }


# clean up
#system "rm -f $TMP_DALICON"

########################## SUBROUTINES ########################################

sub sortlist_z { #>> select best per each z12_domain
	my($zcutoff)=@_;
        my %x;
        foreach (keys %lali_induced) {
		#warn "# sortlist_z  key: $_ inducer: $z12{$inducer{$_}} lali_induced: $lali_induced{$_}\n";
                next if($lali_induced{$_}<$MIN_LALI);
                my $x=0.0;
		if(defined($inducer{$_})) { $x=$z12{$inducer{$_}}; }
                my $y=0.0;
		if(defined($z23{$_})) { $y=$z23{$_}; }
                if ($y<$x) { $x=$y; } # min z-score
                next if($x<$zcutoff); # poor hit
		my $a=0.0; if(defined($maxtried_lali{$_})) { $a=$maxtried_lali{$_}; }
		my $b=0.0; if(defined($lali_observed{$_})) { $b=$lali_observed{$_}; }
                if( ($lali_induced{$_} > $a) && ($lali_induced{$_} > $b) ) { $x{$_}=$x; }
        }
        my(@keys)=sort { $x{$b} <=> $x{$a} } keys %x;
        return(@keys);
}

sub sortlist_lalifraction {
	my($zcutoff)=@_;
        my %x;
        foreach (keys %lali_induced) {
                next if($lali_induced{$_}<$MIN_LALI);
                my $x=$z12{$inducer{$_}};
                my $y=$z23{$_};
                if ($y<$x) { $x=$y; } # min z-score
                next if($x<$zcutoff); # poor hit
                my $f=$lali_induced{$_}/(1+$lali_observed{$x});
                if($f > $x{$_}) { $x{$_}=$f; }
        }
        my(@keys)=sort { $x{$b} <=> $x{$a} } keys %x;
        return(@keys);
}

sub sortlist_laligain {
	my($zcutoff)=@_;
        my %x;
        foreach (keys %lali_induced) {
                #warn "# sortlist_z  key: $_ inducer: $z12_domain{$inducer{$_}}\n";
                next if($lali_induced{$_}<$MIN_LALI);
                my $x=$z12{$inducer{$_}};
                my $y=$z23{$_};
                if ($y<$x) { $x=$y; } # min z-score
                next if($x<$zcutoff); # poor hit
                if( ($lali_induced{$_} > $maxtried_lali{$_}) && ($lali_induced{$_} > $lali_observed{$_}) ) { $x{$_}=$lali_induced{$_}-$maxtried_lali{$_}-$lali_observed{$_} ; }
        }
        my(@keys)=sort { $x{$b} <=> $x{$a} } keys %x;
        return(@keys);
}

sub print_stack {
        my($iteratien)=@_;
        my %x;
        foreach (keys %z12,keys %z23) { $x{$_}=1; }
        my $cd3;
        my $rank=0;
        foreach $cd3 (keys %z23) {
                $rank++;
                my $cd2=$inducer{$cd3};
                print join("\t",$iteration,$rank,$cd3,$iteration{$cd3},$z23{$cd3},$lali_induced{$cd3},$maxtried_lali{$cd3},$lali_observed{$cd3},$cd2,$z12{$cd3}),"\n";
        }
}

sub run_dalicon {
        my($cd1,$zcutoff,$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,@cd3list)=@_;
        warn "# This is run_dalicon cd1=$cd1 zcutoff=$zcutoff\n";
        # do batch of 100 best
        my $cd3;
        open(OUT,">$TMP_DALICON");
        print OUT "$cd1\n";
        my $n_cd3=0;
        foreach $cd3 (@cd3list) {
                next if(defined($tested{$cd3}) && $tested{$cd3}>2); # case: cd3 produces no improvement because it is not in pdb.list
                #last if($n_cd3>$BATCH_SIZE);
                my $a=0.0; if(defined($lali_induced{$cd3})) { $a=$lali_induced{$cd3}; }
		my $b=0.0; if(defined($maxtried_lali{$cd3})) { $b=$maxtried_lali{$cd3}; }
		my $x=$a-$b;

                # output cd3,z13, preali to 'tmp.dalicon'
                my($l13,$m13)=&generate_map13($map21{$inducer{$cd3}},$map23{$cd3});
                next if($lali_induced{$cd3} < $MIN_LALI); # reject short transitive alignments
                #warn "# generate_map13 returned $cd3 $l13 $m13\n";
                printf OUT "%-5.5s\*\n", $cd3;
                print OUT "$l13\n$m13";
                $n_cd3++;
                print "# dalicon iteration=$iteration cd1=$cd1 cd3=$cd3 z12=$z12{$inducer{$cd3}} z23=$z23{$cd3} x=$x n_cd3=$n_cd3\n";
		$x=0.0; if(defined($maxtried_lali{$cd3})) { $x=$maxtried_lali{$cd3}; }
                if($l13>$x) { $maxtried_lali{$cd3}=$l13; }
                $tested{$cd3}++;
        }
        print OUT "END\nEND\n";
        close(OUT);
        # run dalicon
        &dalicon($TMP_DALICON,"T",$COMPARE_EXE,$DALIDATDIR_1,$DALIDATDIR_2,0);
        return($n_cd3);
}

sub generate_map13 {
        my($map21,$map23)=@_;
        # generate induced alignment map13 from map21 x map23
        my %m21;
        my ($cd2starts,$cd1starts,$lengths)=split(/ /,$map21);
        my(@cd2starts)=split(/,/,$cd2starts);
        my(@cd1starts)=split(/,/,$cd1starts);
        my(@lengths)=split(/,/,$lengths);
        my $from2;
        foreach $from2 (@cd2starts) {
                my $from1=shift(@cd1starts);
                my $l=shift(@lengths);
                foreach my $i (0..$l-1) { $m21{$from2+$i}=$from1+$i; }
        }
        my %m23;
	my $cd3starts;
        ($cd2starts,$cd3starts,$lengths)=split(/ /,$map23);
        (@cd2starts)=split(/,/,$cd2starts);
        my(@cd3starts)=split(/,/,$cd3starts);
        (@lengths)=split(/,/,$lengths);
        foreach $from2 (@cd2starts) {
                my $from3=shift(@cd3starts);
                my $l=shift(@lengths);
                foreach my $i (0..$l-1) { $m23{$from2+$i}=$from3+$i; }
        }
        my $i2;
        my $l13=0;
        my @x;
        my @y;
        foreach $i2 (sort { $a <=> $b } keys %m21) {
                my $i1=$m21{$i2};
		next if(!defined($m23{$i2}));
                my $i3=$m23{$i2};
                if($i3>0) { push(@x,$i1); push(@y,$i3); push(@x,$i1); push(@y,$i3); $l13++; } # ranges!
        }
        my $m13="@x\n@y\n"; # tmp.dalicon input!
        return($l13,$m13);
}

sub induce {
        my($cd1,$zcutoff,$cd2)=@_;
        if($z12{$cd2}<$zcutoff) { return; } # too poor quality
        #print "# this is induce: $cd2\n";
#  use dali; select 1 as orderflag, cd2 as cd3, zscore, t1.align_id, from1, from2, blocklen from dali_dccp as t1 left join dali_segments as t2 on t1.align_id=t2.align_id where zscore>2 and cd1='2anuA' order by zscore desc, t1.align_id, from1; select 2 as orderflag, cd1 as cd3, zscore, t1.align_id, from1, from2, blocklen from dali_dccp as t1 left join dali_segments as t2 on t1.align_id=t2.align_id where zscore>2 and cd2='2anuA' order by zscore desc, t1.align_id, from1;
#orderflag       cd3     zscore  align_id        from1   from2   blocklen
#1       2anuA   99.9    2488932 1       1       224
#1       2anuF   43.7    2488123 1       1       224
#1       2anuE   43.5    2488124 1       1       223
#1       2anuC   43.5    2488125 1       1       141
#1       2anuC   43.5    2488125 143     142     81
#1       2anuB   43.0    2488126 1       1       140
#1       2anuB   43.0    2488126 144     141     81
        my $oldalignid=0;
        my $map23;
        my @xstarts;
        my @ystarts;
        my @lengths;
        my $orderflag=0;
        my $cd3;
        my $z23;
        #warn "# call $cd1 get_FSSP_shell_segments($cd2,$zcutoff)\n";
        my(@tmp)=&get_FSSP_shell_segments($cd2,$zcutoff);
        foreach (@tmp) {
                chomp;
                next if(/orderflag/); # title row
                my(@x)=split(/\t/);
		if($#x<3) { warn "# get_FSSP returned too short line @x\n"; next; }
                my $alignid=$x[3];
                if($alignid != $oldalignid) {
                        if($oldalignid>0) {
                                if($orderflag==1) { $map23=join(' ',join(',',@xstarts), join(',',@ystarts), join(',',@lengths)); } else  { $map23=join(' ',join(',',@ystarts), join(',',@xstarts), join(',',@lengths)); }
                                &induce_alignment($cd2,$cd3,$z23,$map23);
                        }
                        # init new alignid
                        @xstarts=();
                        @ystarts=();
                        @lengths=();
                }
                $cd3=$x[1];
                $z23=$x[2];
                $oldalignid=$alignid;
                $orderflag=$x[0];
                my $from1=$x[4];
                my $from2=$x[5];
                my $blocklen=$x[6];
                push(@xstarts,$from1);
                push(@ystarts,$from2);
                push(@lengths,$blocklen);
        }
        # last alignment
        if($oldalignid>0) {
                if($orderflag==1) { $map23=join(' ',join(',',@xstarts), join(',',@ystarts), join(',',@lengths)); } else  { $map23=join(' ',join(',',@ystarts), join(',',@xstarts), join(',',@lengths)); }
                &induce_alignment($cd2,$cd3,$z23,$map23);
        }
}

sub induce_alignment {
        my($cd2,$cd3,$z23,$map23)=@_;
        # generate induced alignment map13 from map21 x map23
        my($l13,$m13)=&generate_map13($map21{$cd2},$map23);
        # case: z12{$cd3} is defined => only keep if induced alignment is longer than maxtried
        # case: z12{$cd3} is undefined and z23{$cd3} is defined => only keep if induced alignment is longer than maxtried and current induced
        # case: z12{$cd3} is undefined and z23{$cd3} is undefined => keep induced alignment
        my $keep=0;
        if(defined($z12{$cd3})) {
                if($l13 > $maxtried_lali{$cd3}) { $keep=1; }
        } else {
                if(!defined($z23{$cd3})) {
                        $keep=1;
                } elsif(defined($maxtried_lali{$cd3}) && defined($lali_induced{$cd3})) {
                        if(($l13 > $maxtried_lali{$cd3}) && ($l13 > $lali_induced{$cd3})) { $keep=1; }
                }
        }
        # overwrite cd3 alignment
        if($keep) {
                $inducer{$cd3}=$cd2;
                $lali_induced{$cd3}=$l13;
                $map23{$cd3}=$map23;
                $z23{$cd3}=$z23;
                #warn "# $cd2 induced alignment of $cd3 $l13 $z23 $map23\n";
        }
}

sub read_dccp {
        my($dccpfile,$cd1,$zcutoff)=@_;
        my(@lines)=&sortdccp($dccpfile);
        my $n=0;
        my $block='';
        my $cd2;
        my %order;
        my %alignment12;
        my %z;
        my $keep=0;
        #while(<IN>) {
        foreach(@lines) {
                if (/DCCP/) {
                        if($block ne '' and $keep) { $alignment12{$cd2}=$block; }
                        my(@x)=split(/\s+/);
                        my $cdb=pop(@x);
                        my $cda=pop(@x);
                        $cd2=$cdb;
                        $order{$cd2}=1;
                        if($cdb eq $cd1) { $cd2=$cda; $order{$cd2}=2; }
                        my($nfrag)=pop(@x); pop(@x);
                        my($z)=pop(@x);
			if(!defined($z{$cd2})) { $keep=1; }
                        elsif($z>$z{$cd2}) { $keep=1; } else { $keep=0; }
                        if($keep) { $z{$cd2}=$z; }
                        if($z>=2) { $n++; }
                        $block="$nfrag\n";
                        #warn "# read_dccp got cda=$cda cdb=$cdb z=$z cd2=$cd2 z12=$z12{$cd2} nhit=$n\n";
                } elsif(/alignment/) {
                        next;
                 } else {
                        # fix nres>1000 here, later split on whitespace
# dccp   1    384.7 6.1 139     0.6          78     21                mol1  1yqwr
# alignment
# 338   341 345   350 384   388 898   902 904   908 915   9191104  11081123  1128
#1129  11451150  11541155  11631166  11691170  11751206  12091210  12131304  1308
#1310  13131314  13201338  13471357  13751376  1379
#   8    11  34    39  42    46  54    58  59    63  64    68  82    86  87    92
#  95   111 112   116 118   126 163   166 168   173 198   201 265   268 290   294
# 300   303 309   315 416   425 426   444 446   449
                        s/^\s+//;
                        my $line=$_;
                        foreach (split(/\s+/,$line)) {
                                if($_>9999) {
                                        my($a,$b)=/^(\d+)(\d{4})$/;
                                        $block.="$a $b ";
                                } else {
                                        $block.="$_ ";
                                }
                        }
                }
        }
        #close(IN);
        # last alignment
        if($block ne '' and $keep) { $alignment12{$cd2}=$block; }

        # convert alignment12{cd2} and order{cd2} to map21{cd2}, lali_observed{cd2}, maxtried_lali{cd2}
        foreach $cd2 (keys %alignment12) {
                #warn "# read_dccp cd2=$cd2 z=$z{$cd2} z12=$z12{$cd2}\n";
		my $x=0.0;
		if(defined($z12{$cd2})) { $x=$z12{$cd2}; }
                if($z{$cd2}>$x) { # if z-score is better, then overwrite & induce
                        $iteration{$cd2}=$iteration;
                        # build alignment
                        my $order=$order{$cd2};
                        my @xstarts;
                        my @ystarts;
                        my @lengths;
                        my $lali=0;
                        $_=$alignment12{$cd2};
                        my(@x)=split(/\s+/);
                        my($nfrag)=$x[0];
                        my $i=1;
                        while($i<=$nfrag*2) {
                                my $from1=$x[$i];
                                my $to1=$x[$i+1];
                                my $jres=$x[$i+$nfrag*2];
                                my $length=$to1-$from1+1;
                                if ($order==2) {
                                        my $x=$from1;
                                        $from1=$jres;
                                        $jres=$from1;
                                }
                                push(@xstarts,$from1);
                                push(@ystarts,$jres);
                                push(@lengths,$length);
                                $lali+=$length;
                                $i+=2;
                        }
                        $z12{$cd2}=$z{$cd2};
                        $map21{$cd2}=join(' ',join(',', @ystarts), join(',', @xstarts), join(',', @lengths)); # @cd2starts @cd1starts @lengths
                        $lali_observed{$cd2}=$lali;
                        $lali_induced{$cd2}=$lali;
                        if(!defined($maxtried_lali{$cd2})) { $maxtried_lali{$cd2}=0; }
                        $inducer{$cd2}=$cd1; # because cd2 now has direct alignment to cd1
                        #print "# read_dccp: cd2=$cd2 lali=$lali z12=$z12{$cd2}\n";

                        # induce new neighbours using improved alignments to cd2
                        &induce($cd1,$zcutoff,$cd2);
                }
        }

        return($n);
}


} # ens sub walk

return(1); # end module

__END__


