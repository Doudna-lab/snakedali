#use lib qw(/home/luholm/mpidali4/);
use strict;
use IPC::Open3;
use mpidali;

####################################################################
# CONFIGURATION SECTION --> MODIFY PATHS                           #
# ##################################################################

# We use a web server to retrieve computed structural similarities 
my $FSSP_LOCAL_MYSQL_SERVER=0; # 0 = remote; 1 = local
my $FSSP_SERVER_URL='http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/FSSP.cgi';
# below are needed if FSSP_LOCAL_MYSQL_SERVER is 1
my $DCCPTABLE='dali.dccp';
my $SEGMENTSTABLE='dali.segments';
my $LOCAL_MYSQL_EXE="mysql -uluholm -px8Jv8EtNjk";


####################################################################
# END OF CONFIGURATION SECTION                                     #
# ##################################################################


####################################################################
# retrieve data from precomputed konowledgebase FSSP_SERVER_URL    #
####################################################################
# retrieve segments
# curl --form cd1=2pnkA --form zcutoff=32.0  http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/FSSP.cgi 
#
# retrieve cdlist
# curl --form cd1=2pnkA --form zcutoff=32.0  --form dccp=1 http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/FSSP.cgi 


sub get_FSSP_shell_dccp {
	if($FSSP_LOCAL_MYSQL_SERVER > 0) {
		return(&get_FSSP_shell_dccp_local(@_));
	} else {
		return(&get_FSSP_shell_dccp_remote(@_));
	}
}

sub get_FSSP_shell_segments {
        if($FSSP_LOCAL_MYSQL_SERVER > 0) {
                return(&get_FSSP_shell_segments_local(@_));
        } else {
                return(&get_FSSP_shell_segments_remote(@_));
        }
}
sub get_FSSP_shell_dccp_remote {
	my($cd2,$zcutoff)=@_;
	if($zcutoff<2) { $zcutoff=2.0; }
	my $cmd="curl --form cd1=$cd2 --form zcutoff=$zcutoff --form dccp=1 $FSSP_SERVER_URL";
	warn "$cmd\n";
	my(@lines)=`$cmd`;
	my %x;
	foreach (@lines) {
		chomp;
		if(/^(\w{5})$/) { $x{$1}=1; }
	}
	return(keys %x);
}

sub get_FSSP_shell_segments_remote {
        my($cd2,$zcutoff)=@_;
        if($zcutoff<2) { $zcutoff=2.0; }
        my $cmd="curl --form cd1=$cd2 --form zcutoff=$zcutoff $FSSP_SERVER_URL";
	warn "$cmd\n";
        return(`$cmd`);
}

sub get_FSSP_shell_dccp_local {
	my($cd1,$zcutoff)=@_;
	if($zcutoff<2) { $zcutoff=2.0; }
	my %x;
        my $cmd1="$LOCAL_MYSQL_EXE -e \"select distinct(cd1) from $DCCPTABLE where cd2='$cd1' and zscore>$zcutoff order by zscore desc;\"";
        my $cmd2="$LOCAL_MYSQL_EXE -e \"select distinct(cd2) from $DCCPTABLE where cd1='$cd1' and zscore>$zcutoff order by zscore desc;\"";
        my(@a)=`$cmd1`;
        my(@b)=`$cmd2`;
	foreach (@a,@b) {
		chomp;
		if(/^(\w{5})$/) { $x{$1}=1; }
	}
	return(keys %x);
}

sub get_FSSP_shell_segments_local {
	my($cd1,$zcutoff)=@_;
        if($zcutoff<2) { $zcutoff=2.0; }
        my $cmd="$LOCAL_MYSQL_EXE -e \"select distinct 1 as orderflag, cd2 as cd3, zscore, t1.align_id, from1, from2, blocklen from dali.dccp as t1 left join $SEGMENTSTABLE as t2 on t1.align_id=t2.align_id where zscore> $zcutoff and cd1='$cd1'; select distinct 2 as orderflag, cd1 as cd3, zscore, t1.align_id, from1, from2, blocklen from dali.dccp as t1 left join $SEGMENTSTABLE as t2 on t1.align_id=t2.align_id where zscore> $zcutoff and cd2='$cd1';\"";
        my(@lines)=`$cmd`;
        return(@lines);
}

return(1);

