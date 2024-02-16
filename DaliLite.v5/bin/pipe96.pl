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

while(<STDIN>) {
        my $z=substr($_,0,5);
        my $key=substr($_,11,5);
        $n{$key}++;
        if( ($n{$key} <= $Nbest) && ($z >= $cutoff) ) { print $_; }
}

