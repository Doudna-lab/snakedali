# input: distmat
# output: newick

use strict;

sub upgma {
        my($listfile,$simtablefile,$matrixfile,$xmlfile,$treefile,$heatmapfile)=@_;

# pick up annotations from listfile
my %anno;
open(IN,"<$listfile");
while(<IN>) {
        chomp;
        my($cd1)=/^(\S+)\s*(.*)$/;
        $anno{$cd1}=$2;
	$anno{$cd1}=~s/[():,;]//g; # remove characters confusing parsing of newick string
        warn "#anno of $cd1 is $anno{$cd1}\n";
}
close(IN);

# input similarity matrix
my %sim;
open(IN,"<$simtablefile") || die "Can't open $simtablefile\n";
my @lines=<IN>;
close(IN);
$_=$lines[0];
my($n)=/(\d+)/;
# parse table;
my %list;
my %matrix;
my %orig_matrix;
foreach my $i (1..$n) {
        $_=$lines[$i];
        chomp;
        s/^\s*//;
        my(@x)=split(/\s+/);
        $list{$i}=$x[0];
        foreach my $j (1..$n) {
                $matrix{$i}{$j}=$x[$j];
                $orig_matrix{$i}{$j}=$x[$j];
        }
}
# symmwtrize
foreach my $i (1..$n) {
        foreach my $j ($i+1..$n) {
                if($matrix{$i}{$j}>$matrix{$j}{$i}) { $matrix{$j}{$i}=$matrix{$i}{$j}; }
                else { $matrix{$i}{$j}=$matrix{$j}{$i}; }
        }
}
# initialize
my %newick;
my %depth;
my %size;
my %members;
foreach my $i (1..$n) {
        $newick{$i}="$list{$i} $anno{$list{$i}}"; # list key is number, anno key is cd1
        $newick{$i}=~s/\s+$//;
        warn "# init newick $i as $newick{$i}\n";
        $depth{$i}=$matrix{$i}{$i};
        $size{$i}=1;
        $members{$i}="$i";
}
# agglomerative contraction: overwrite cluster head
my $i=0;
my $j=0;
foreach (2..$n) {
        # find best group-edge
        my $best=0;
        foreach my $x (1..$n) {
                next if($size{$x}<1);
                foreach my $y ($x+1..$n) {
                        next if($size{$y}<1);
                        my $z=$matrix{$x}{$y};
                        if($z>$best) { $i=$x; $j=$y; $best=$z; }
                }
        }
        # heavier left
        if($size{$j}>$size{$i}) {
                my $e=$i;
                $i=$j;
                $j=$e;
        }
        # merge i<-j
        my $d1=$depth{$i}-$best; if($d1<0) { $d1=0; }
        my $d2=$depth{$j}-$best; if($d2<0) { $d2=0; }
        $depth{$i}=$best;
        $newick{$i}="($newick{$i}\:$d1\,$newick{$j}\:$d2)";
        warn "# internal node $i:$d1 $j:$d2 $newick{$i}\n";
        # overwrite matrix{$i}{*} with average similarities
        my $m=$size{$i}+$size{$j};
	if($m<1) { $m=1; }
        foreach my $k (1..$n) {
                my $s=$size{$i}*$matrix{$i}{$k}+$size{$j}*$matrix{$j}{$k};
                $matrix{$i}{$k}=$s/$m;
        }
        $size{$i}=$m;
        $size{$j}=0;
        $members{$i}=join(" ",$members{$i},$members{$j});
        $members{$j}="";
}
$_="$newick{$i}";
s/ /_/g;
my $newick_unrooted="$newick{$i}\;";
$newick{$i}=$_;
my $newickstring="\($newick{$i}\:$depth{$i}\)\;";
my $inputstring="\($newick{$i}\:0.1\)\;";
my $xml=&convert_newick_to_xml($inputstring,$depth{$i});
open(OUT,">newick");
print OUT "$newickstring\n";
close(OUT);
open(OUT,">newick_unrooted");
print OUT "$newick_unrooted\n";
close(OUT);
# output similarity matrix in tree order
open(OUT,">$matrixfile");
my(@x)=reverse split(/\s+/,$members{$i});
my $cdstring="[";
my $zstring="[";
print OUT $n,"\n";
foreach my $i (@x) {
        print OUT $list{$i};
        $cdstring.="\"$list{$i}\",";
        $zstring.="[";
        foreach my $j (@x) {
                print OUT "\t",$orig_matrix{$i}{$j};
                $zstring.="$orig_matrix{$i}{$j},";
        }
        $zstring.="],";
        print OUT "\n";
}
$zstring.="]";
$zstring=~s/,]/]/g;
$cdstring.="]";
$cdstring=~s/,]/]/;
close(OUT);
# output xml file
open(OUT,">$xmlfile") || die "Can't open $xmlfile\n";
print OUT "$xml\n";
close(OUT);
# output html file
open(OUT,">$treefile") || die "Can't open $treefile\n";
print OUT &html_xmltree($xmlfile);
close(OUT);
# output heatmap
my $pixels=$n*20;
if($pixels<800) { $pixels=800; }
&html_heatmap($cdstring,$zstring,$pixels,$heatmapfile);

} # end sub upgma

sub html_newicktree { # not used
        my($newickstring,$outfile)=@_;
        open(OUT,">$outfile") || die "html_tree Can't open $outfile\n";
        print OUT<<EOB;
<html>
<head>
        <script type="text/javascript" src="/js/raphael/raphael-min.js"></script>
        <script type="text/javascript" src="/js/jsphylosvg-min.js"></script>
        <script type="text/javascript">
        window.onload = function(){
                        var dataObject = { newick: \'$newickstring\' };
                        phylocanvas = new Smits.PhyloCanvas(
                                dataObject,
                                'svgCanvas',
                                800,800
//                      ,'circular'
                        );
        };
        </script>
</head>
<body>
        <div>Structural similarity dendrogram.</div>
        <div id="svgCanvas"> </div>
</body>
</html>
EOB
        close(OUT);
}

sub html_xmltree {
	my($xmlfilename)=@_;
	my $html="";
        $html.=<<EOB;
<html>
<head>
    <script type="text/javascript" src="/js/jquery/jquery-1.4.4.min.js" ></script>
    <script type="text/javascript" src="/js/raphael/raphael-min.js" ></script>
    <script type="text/javascript" src="/js/jsphylosvg-min.js"></script>
    <script type="text/javascript" src="/js/yui.js"></script>

    <script type="text/javascript">
        window.onload = function(){
                YUI().use('oop', 'json-stringify', 'io-base', 'event', 'event-delegate', function(Y){
                        var uri = "$xmlfilename";
                        function complete(id, o, args) {
                                var data = o.responseXML; // Response data.
                                var dataObject = {
                                                        xml: data,
                                                        fileSource: true
                                                };
                                phylocanvas = new Smits.PhyloCanvas(
                                        dataObject,
                                        'svgCanvas',
                                        800, 800
//                                      ,'circular'
                                );
                        };
                        Y.on('io:complete', complete, Y);
                        var request = Y.io(uri); // links to web page
                });
        };
</script>

</head>
<body>
    <div>Structural similarity dendrogram.</div>
The dendrogram is derived by average linkage clustering of the structural similarity matrix (Dali Z-scores).
    <div id="svgCanvas"> </div>
</body>
</html>
EOB
        return($html);
}

sub html_heatmap {
        # generate HTML page with heatmap of reordered similarity table
        my($cdstring,$zstring,$pixels,$outfile)=@_;
        open(OUT,">$outfile") || die "html_heatmap Can't open $outfile\n";
        print OUT<<EOB;
<head>
  <!-- Plotly.js -->
  <script src="/js/plotly-latest.min.js"></script>
</head>

<body>
        <div>Structural similarity matrix (Dali Z-scores).</div>

  <div id="heatmap"><!-- Plotly chart will be drawn inside this DIV --></div>
  <script>
   var data = [
  {
    z: $zstring,
    x: $cdstring,
    y: $cdstring,
    type: 'heatmap',
    colorscale: 'Reds'
  }
];
   var layout = { height: $pixels, width: $pixels };

Plotly.newPlot('heatmap', data=data, layout=layout);

  </script>

</body>

EOB
        close(OUT);
}

sub convert_newick_to_xml {
	my($newick,$zroot)=@_;
        $_=$newick;
        #warn "# newick string is $_\n";
        s/\(/(\n/g;
        s/\)/\n)\n/g;
        s/:/\n:/g;
        s/,/\n,\n/g;
        #warn "# parsed string is\n$_\n";

        my $xml="";

        $xml.="<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"
 xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.10/phylo
xml.xsd\" xmlns=\"http://www.phyloxml.org\">\n";
        $xml.="<phylogeny rooted=\"true\">\n";

        my $indent="";
        foreach (split(/\n/)) {
                if(/\(/) { # (
                        $xml.="$indent<clade>\n";
                        $indent.="\t";
                } elsif(/\)/) { # )
                        chop($indent);
                        $xml.="$indent</clade>\n";
                } elsif(/^:(\S+)/) { # :1.234
                        $xml.="$indent<branch_length>$1</branch_length>\n";
                } elsif(/,/) { # clade separator
                        chop($indent);
                        $xml.="$indent</clade>\n";
                        $xml.="$indent<clade>\n";
                        $indent.="\t";
                } elsif(/^;*$/) { # skip punctuation
                        next;
                } else { # node name
                        my($name,$anno)=/^([a-zA-Z0-9]+)\_*\s*(.*)/; ##/^([a-zA-Z0-9]+)[_ ]*(.*)$/;
			warn "got ($name,$anno) from $_\n";
                        $xml.="$indent<name>$name $anno</name>\n";
                        $xml.="$indent<annotation>\n";
                        $xml.="$indent\t<desc>$anno</desc>\n";
                        $xml.="$indent\t<uri>$name\.html<\/uri>\n";
                        $xml.="$indent</annotation>\n";
                }
        }
	# patch zroot
	warn "# patch before: $xml\n";
	$xml =~ s/<\/clade>\s*$//;
	$xml.="<name>Z=".sprintf("%5.1f",$zroot)."</name>\n$indent</clade>\n";
	warn "# patch after: $xml\n";
	# finish off
        $xml.="</phylogeny>\n";
        $xml.="</phyloxml>\n"; 
       return($xml);
}

#######################################################################
#
#######################################################################
return 1;

