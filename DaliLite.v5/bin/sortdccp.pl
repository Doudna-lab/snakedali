#!/usr/bin/perl
#______________________________________________________________________________
# Title     : sortdccp.pl
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
# sortdccp.perl: read a DCCP file from STDIN, sort by Daliscore, write to STDOUT
#

use strict;

my($nkeep)=@ARGV; # optional to filter nkeep best per cd1-cd2
if($nkeep<1) { $nkeep=999999; } # keep all

my $TRUE=1;
my $FALSE=0;
my @array;
my $first=$TRUE;
my $self=$FALSE;
my $block="";
while(<STDIN>) {
#DCCP   1   2149.4 0.0 157    34.4         100      1                1mup  1mup
        if(/^ DCCP\s+(\d+)\s.*\s(\w+)\s+(\w+)$/) {
                next if($1!=1); # remove redundant DCCP   2 ... lines
                # save old block
                #??? if(!$first) { if($self) { } else {  push(@array,$block); } } else { $first=$FALSE; }
                if($block ne "") { push(@array,$block); }
		$self=($2 eq $3);
                #warn "# got first=$first self=$self from ($2 eq $3) <= $_\n";
                # initialize new block
                $block = $_;
        } else {
                $block .= $_;   # append to block
        }
}
# last line
push(@array,$block);

my $i;
my @datakeys;
foreach $i ($[..$#array) {
# DCCP   1   8983.1 0.0 451    66.8         100      1                1j5sA 1j5sA
        $_=$array[$i]; push(@datakeys,(/^ DCCP.{21}\s+([\-\.\d]+)/));
#        $_=$array[$i]; push(@datakeys,(/^ DCCP\s+\d+\s+([\-\.\d]+)/));
}
#warn "datakeys: @datakeys\n";
sub bydatakeys { $datakeys[$b] <=> $datakeys[$a]; }
my (@tmp)=sort bydatakeys $[..$#array;
#warn "tmp: @tmp\n";
my @sortarray=@array[@tmp];      # sort array on value
#warn "sortarray: @sortarray\n";

# output
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
        print $block;
}


