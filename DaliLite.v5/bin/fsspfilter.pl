#!/usr/bin/perl
#
# fsspfilter.pl
#
# input SORTED dccp file + cd1, passes through highest z-score per protein pair, sorted
#

$TRUE=1;
$FALSE=0;

if (@ARGV < 4) { die "Usage cat 1mup.dccp | fsspfilter.pl 1mup z-cutoff nbest [nbestpercd2]\n";}

my ($cd1,$zcut,$nbest,$nbestpercd2)=@ARGV;

if($nbestpercd2<1) { $nbestpercd2=1; }

my $align_id=0;
my %blockarray;
my %zscore;
my %cd2;
my $skip=1;
while(<STDIN>) {
#DCCP   1   2149.4 0.0 157    34.4         100      1                1mup  1mup
	if(/^ DCCP\s+(\d+)\s+(\S+).*\s+(\w+)\s+(\w+)\s*$/) {
		next if($1!=1); # remove redundant DCCP   2 ... lines
		# save old block
		if(!$skip) {
#			$zold=$zarray{$cd2};
#			if ( $z > $zold ) {  # ali to keep per cd2 selected by Z
#				$zarray{$cd2}=$z;
#				$blockarray{$cd2}=$block;
#			}
			$align_id++;
			$blockarray{$align_id}=$block;
			$zscore{$align_id}=$z;
			$cd2{$align_id}=$cd2;
		}
		# initialize new block
		if($3 ne $cd1 && $4 ne $cd1) {$skip=$TRUE} else {$skip=$FALSE};
		$z=substr($_,26,8);
 		$id1=$3;
		$id2=$4;
		if($id1 eq $cd1) { $cd2=$id2;} else  {$cd2=$id1;}
		$block = $_;
	} else {
		$block .= $_;	# append to block
	}
}
# last line
if (!$skip) {
#	$zold=$zarray{$cd2};
#	if ( $z > $zold ) {
#		$zarray{$cd2}=$z;
#		$blockarray{$cd2}=$block;
#	}
	$align_id++;
	$blockarray{$align_id}=$block;
	$zscore{$align_id}=$z;
	$cd2{$align_id}=$cd2;
}

#foreach $cd2 ( keys(%blockarray) ) {
#	$key=$zarray{$cd2} . " " . $cd2; # compose sort key from Z + cd2
#	@srt{$key}=$cd2;
#}

$n=0;
my %npercd2;
#foreach $z (reverse sort(keys(%srt))) {	# sort associative array on keys
my $m=$align_id;
#foreach $align_id (keys %zscore) {
foreach $align_id (1..$m) {
	$n++;
	$npercd2{$cd2{$align_id}}++;
	next if($npercd2{$cd2{$align_id}}>$nbestpercd2);
	#warn "# reverse sorted n=$n z=$z align_id=$align_id\n# block=$blockarray{$align_id} \n";
	#$_=$z; ($sim)=/([\d\.]+)\s+(\w+)/;
	#last if($sim<$zcut && $n>$nbest);
	last if($zscore{$align_id}<$zcut && $n>$nbest);
	#print $blockarray{$srt{$z}};	# transitive association !
	print $blockarray{$align_id};
}



