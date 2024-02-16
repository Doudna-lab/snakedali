# convert filtered DCCP format to falicon input

use strict;

my $nblock=0;
my $oldcd1='?????';
while(<STDIN>) {
	chomp;
	# DCCP   1   2027.2 0.0 195    19.6          87    220                4cpaB 3k2kA
	if(/DCCP.*\s+(\d+)\s+(\w+)\s+(\w+)\s*$/) {
		print "\n";
		if($2 ne $oldcd1) {
			$oldcd1=$2;
			print "END\n$2\n";
		}
		print "$3\*\n";
		print "$1\n";
	} elsif(/alignment/) {
		next;
	} else { # alignments
		s/(\d{4})$/ $1/;
		s/(\d)(\d{4}) /$1 $2 /g;
		print "$_\n";
	}
}
print "\nEND\n";
