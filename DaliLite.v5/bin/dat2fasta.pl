#!/usr/bin/perl
# input: pdb.fasta
# output: description added from DAT/

use strict;

if($#ARGV<0) { die "USAGE: $0 <dalidatdir>\n"; }
my($DATDIR)=@ARGV; #"/data/DaliLite/DAT";

my $seq='';
while(<STDIN>) {
	#warn "IN: $_\n";
		my($id)=/(\w+)/;
		# get description
		my $file="$DATDIR\/$1\.dat";
#	-compnd   " MOLECULE: AVIAN PANCREATIC POLYPEPTIDE;                              "
		my $cmd="egrep '^-compnd|^-sequence' $file";
		my(@lines)=`$cmd`;
		foreach(@lines) {
			chomp;
			if(/^\-compnd/) {
				s/^\S+\s+//;
				s/\"//g;
				s/MOLECULE://;
				s/\;\s+$//;
				print "\>$id $_\n$seq\n";
				$seq='';
			} elsif(/^\-sequence/) {
				s/^\S+\s+//;
                                s/\"//g;
				$seq=$_;
			}
		}
}
