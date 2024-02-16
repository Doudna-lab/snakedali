#!/usr/bin/perl
#______________________________________________________________________________
# Title     : pipe96.perl
# Usage     :
# Function  :
# Example   :
# Keywords  :
# Options   :
# Author    : jong@bio.cc,
# Category  :
# Returns   :
# Version   : 1.0
#------------------------------------------------------------------------------
#
#   pipe96.perl
#
#       selects Nbest pairs per cd1cd2 with zscore > cutoff [first argument]
#
#       use as: sort fort.96 | pipeforwarder.perl cutoff Nbest | pipe96.perl
#

use strict;

my ($cutoff,$Nbest) = @ARGV;

my %n;
my $oldcd1='?????';

while(<STDIN>) {
	chomp;
	s/^\s+//;
	s/\-99//g;
	my($z,$key,@x)=split(/\s+/);
        $n{$key}++;
	my $cd1=substr($key,0,5);
	my $cd2=substr($key,5,5);
        if( ($n{$key} <= $Nbest) && ($z >= $cutoff) ) { 
		if($cd1 ne $oldcd1) { 
			if($oldcd1 ne '?????') { print "END\n"; }
			print "$cd1\n";
		}
		$oldcd1=$cd1;
		# reformat for dalicon
		shift(@x); # idom
		shift(@x); # score
		my $n=$#x/4;
		shift(@x); # nseq with -99 placeholders
		print "$cd2\*\n$n\n";
		print join(" ",@x[0..2*$n-1]),"\n";
		print join(" ",@x[2*$n..4*$n-1]),"\n";

	}
}
print "END\nEND\n";
