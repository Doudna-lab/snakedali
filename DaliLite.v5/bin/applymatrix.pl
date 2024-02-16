#!/usr/bin/perl
# apply translation-rotation matrix from STDIN to PDB-file, write PDB-file with transformed coordinates to STDOUT
# grep matrix lines from dali.pl --outfmt '... transrot ...' output
#
# Translation-rotation matrices
#-matrix  "1ppt-A 1bba-A  U(1,.)   0.631906 -0.761372 -0.144939           -0.890845"
#-matrix  "1ppt-A 1bba-A  U(2,.)   0.512616  0.550832 -0.658642          -10.882093"
#-matrix  "1ppt-A 1bba-A  U(3,.)   0.581308  0.341902  0.738366            4.946664"

use strict;
if($#ARGV<0) { die "USAGE: $0 <pdbfile> < 3-lines-representation-of-rotation-matrix-and-translation-vector > transformed.pdbfile\n"; }
my $tmpfile="$$.tmp";
my($pdbfile)=@ARGV;
my $cleanup;
# uncompress gzipped
if($pdbfile =~ /\.gz$/) {
	warn "# Detected gzipped file - uncompressing\n";
	system("gzip -cd $pdbfile > $tmpfile");
	$pdbfile=$tmpfile;
	$cleanup=1;
}
# get matrix from STDIN
my $u12=my $u13=my $u21=my $u23=my $u31=my $u32=my $t1=my $t2=my $t3=0.0;
my $u11=my $u22=my $u33=1.0;
my $iline=0;
while(<STDIN>) {
	next if(/^#/); # skip comments
	# grep first three lines with four floats
	if(/(\-*\d+\.\d+)\s+(\-*\d+\.\d+)\s+(\-*\d+\.\d+)\s+(\-*\d+\.\d+)/) {
		$iline++;
		warn "# $iline data input: $1 $2 $3 $4 from $_\n";
		if($iline==1) { $u11=$1; $u12=$2; $u13=$3; $t1=$4; }
		elsif($iline==2) { $u21=$1; $u22=$2; $u23=$3; $t2=$4; }
		elsif($iline==3) { $u31=$1; $u32=$2; $u33=$3; $t3=$4; }
		else { last; }
	}
}

# apply transform to PDB file
open(IN,"<$pdbfile") || die "Can't open PDB-file $pdbfile\n";
while(<IN>) {
	if(/^ATOM  / || /^HETATM/) {
#ATOM    229  CA  MET A  17      -8.661   0.837  -3.571  1.00  0.00           C
		my($left,$x,$y,$z,$rite)=/^(.{30})\s*(\S+)\s+(\S+)\s+(\S+)(\s.*$)/;
		my $a=sprintf("%8.3f",$t1+$u11*$x+$u12*$y+$u13*$z);
		my $b=sprintf("%8.3f",$t2+$u21*$x+$u22*$y+$u23*$z);
		my $c=sprintf("%8.3f",$t3+$u31*$x+$u32*$y+$u33*$z);
		print $left.$a.$b.$c.$rite,"\n";
	} else {
		print $_;
	}
}
close(IN);

# clean up
if($cleanup) { system("rm -f $tmpfile"); }

