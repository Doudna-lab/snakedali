#!/usr/bin/perl
#
# read a pdb file
#
use FindBin qw($Bin);
use lib "$FindBin::Bin"; # modules and scripts in same directory
use strict;
use mpidali;
use warnings;
use 5.010;
use Getopt::Long qw(GetOptions);

my $DALIDATDIR="./DAT";
my $PDBDIR="/data/pdb";
my $rsync;
my $help;
my $infile;
my $listfile;
my $short;
my $verbose;
my $clean;
my $tmpfile="$$\.tmp";
my @tmpfiles=('puu.default', 'units.puu', 'subunits.puu', 'domains.puu');

my $USAGE= <<"EOB";
Convert PDB entry into internal data format used by Dali.

* Import single PDB entry to ./DAT:
	$0 --pdbfile <filename> --pdbid <xxxx> [ --dat <path> ]

* Import a subset of the PDB archive:
	$0 --pdblist <filename> [ --dat <path> ]

* Automated PDB mirroring:
	$0 --rsync [ --pdbmirrordir <path> ] [ --dat <path> ]

* Options:
	--dat <path>		directory to store imported data [default: $DALIDATDIR]
	--pdbfile <filename>	PDB formatted file, may be compressed (.gz). <filename> includes path
	--pdbid <xxxx>		four-letter PDB identifier
	--pdblist <filename>	list of PDB entries. <filename> includes path and must match the pattern pdb????.ent
	--rsync			automated PDB mirroring
	--pdbmirrordir <path>	PDB mirror directory [default: $PDBDIR]
	--clean			remove temporary files
	--verbose		verbose

EOB

GetOptions(
	'h|help' => \$help,
	'dat=s' => \$DALIDATDIR,
	'pdbfile=s' => \$infile,
	'pdblist=s' => \$listfile,
	'pdbid=s' => \$short,
	'rsync' => \$rsync,
	'verbose' => \$verbose,
	'clean' => \$clean,
	'pdbmirrordir=s' => \$PDBDIR) or die $USAGE;

if($help) { die $USAGE; }
if(! ($infile || $listfile || $rsync) ) { die $USAGE; }
# main
if(!$DALIDATDIR) { $DALIDATDIR="./"; }
&checkdir($DALIDATDIR);
my @list;
if($rsync) { 
	if(!$PDBDIR) { $PDBDIR="./"; }
	&checkdir($PDBDIR);
	(@list)=&mirror($PDBDIR,$DALIDATDIR); 
} elsif($listfile) {
	&lock();
	open(LIST,"<$listfile") || die "Can't open file: $listfile\n";
	while(<LIST>) {
		my($infile)=/^(\S+)/;
		my($short)=/pdb(\w{4})\.ent/;
		&import_one($infile,$short,$DALIDATDIR,$tmpfile,$verbose);
	}
	close(LIST);
	&unlock();
} else { 
	&lock(); # lock directory (units.puu, subunits.puu, domains.puu)
	&import_one($infile,$short,$DALIDATDIR,$tmpfile,$verbose); 
	push(@list,$short);
	&unlock();
}

# clean up
system("rm -f $tmpfile");
if($clean) {
	system("rm -f @tmpfiles");
	# remove dssp files
	foreach my $cd (@list) { 
		my $fn="$cd\.dssp";
		system("rm -f $fn");
	}
}

exit();

sub checkdir {
	my($dir)=@_;
	if(-d $dir) { return } else { warn "# Directory $dir does not exist - creating it!\n"; system('mkdir $dir'); }
}

sub mirror {
	my($PDBDIR,$DALIDATDIR)=@_;
	# call rsync
	my $mirrorlog="pdb_update.log";
	my @list;
	my $cmd="rsync -rlpt -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/pdb/ $PDBDIR > $mirrorlog";
	if($verbose) { warn "# Starting rsync\n$cmd\n"; }
	system($cmd);
	(@list)=`grep ent.gz $mirrorlog`;
	if($verbose) { warn "# rsync done, $#list new entries in $mirrorlog\n"; }
	# import new
	&lock(); # lock directory (units.puu, subunits.puu, domains.puu)
	foreach (@list) {
		chomp;
		my $pdbfile="$PDBDIR\/$_";
		my($short)=/pdb(\w{4})\.ent/;
		if($verbose) { warn "# importing $short from $_\n"; }
		&import_one($pdbfile,$short,$DALIDATDIR,$tmpfile,$verbose);
	}
	&unlock();
	return(@list);
}

