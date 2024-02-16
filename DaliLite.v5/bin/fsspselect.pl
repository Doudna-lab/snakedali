#!/usr/bin/perl
#
# filter fssp-data: pass through  Represetatives == Z >= 2.0 & in list0
#				 at least nbest representatives
#

($list0file,$cut2,$nbest)=@ARGV;
$TRUE=1; $FALSE=0;

open(IN,$list0file); while(<IN>) {($cd2)=/^(\w+)/; $list0{$cd2}=1; } close(IN);

$ide=$z=$cd1=$cd2='';

$n=0;
while(<STDIN>) {
	next if (!/^-/);
	# pass through search structure
	if(/^-dbsize/ || /^-query/ || /^-nres/ || /^-header/ || /^-nchain/
		|| /^-dsspin/) {
		print $_;
		if(/^-nres/) { $nres1=/(\d+)/; }
	} else { # cache match between -protein lines
		if(/^-summar/) {
			($cd1,$cd2,$z,$lali,$nres2,$ide)=
	/.*\"(\S+)\s+(\S{4,6})\s{0,4}([\d.]+)\s+\S+\s+(\d+)\s+(\d+)\s+(\d+)/;
			$cd1=~s/\-//; # cd1 is always query structure!
			$cd2=~s/\-//;
		}
		if(/^-protein/) { # print old block, init new block
			if(defined($list0{$cd2})&&$z>=$cut2) {print $block; $n++;}
			else { push(@suboptimals,$block); push(@z,$z); } # to print anyway nbest
			$block = $_;
			$ide=$z=$lali=$nres2=$cd1=$cd2='';
		} else {
			$block .= $_;
		}
	}
}
#last line special
if(defined($list0{$cd2})&&$z>=$cut2) {print $block; $n++;}
else { push(@suboptimals,$block); push(@z,$z); } # to print anyway nbest

# always print at least nbest alignments
# print "* * * * * $n $nbest in fsspselect\n";
if($n<$nbest) {
	# sort suboptimals by z-score, print nbest blocks
	foreach $k ($n+1..$nbest) {
		# find index of highest z-score; set used z-score to -1
		$zmax=-10.0; $ibest=0; $i=$[;
		foreach $z (@z) { if($z > $zmax) { $ibest=$i; $zmax=$z; } $i++; }
		# print "round $k: $zmax $ibest $z[$ibest] @z\n";
		$z[$ibest]=-1.0;
		last if($zmax<0);
		print $suboptimals[$ibest];
	}
}

