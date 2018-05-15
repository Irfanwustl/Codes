#! /usr/bin/perl 

if (scalar(@ARGV) != 3) {
    print <<EOF;
USAGE: $0 <list file> <Result File> <output file>

EOF
exit;
}
if (-e "$ARGV[2]") {
    print "File '$ARGV[2]' existed. Delete it? (yes or no) ";
    my $input = <STDIN>;
    if ($input =~ /^yes$/i) {
	unlink $ARGV[2];
    } else {
	die "Could not establish file '$ARGV[2]'.";
    }
}
open INFILE,   "$ARGV[0]"  or die "File '$ARGV[0]' open error!";
#open INFILE2,  "$ARGV[1]"  or die "File '$ARGV[1]' open error!";
open OUTFILE,  ">$ARGV[2]" or die "File '$ARGV[2]' open error!";
open OUTFILE2, ">$ARGV[2]_else" or die "File '$ARGV[2]_else' open error!";
my $ID="";
undef my $contigIDtemp;
undef my @temp;
undef my $start;
undef my $stop;
undef my @line;
undef my @coverage;
undef my $count;
undef my %length;
#$contigIDtemp ="";
#$count =0 ;
#$size= 0;
$bool=1;
$keyword =0;
$minimum =1;
while (<INFILE>){
	chomp;
	$GENEID	 = $_;
    $name{$GENEID} = $GENEID;
}
open INFILE2,  "$ARGV[1]"  or die "File '$ARGV[1]' open error!";
while(<INFILE2>){
    chomp;
    @token = split("\t",$_);
    $GENEID = $_;
    if($name{$GENEID} eq $GENEID){
        print OUTFILE  "$_\n";
    }else{
        print OUTFILE2 "$_\n";
    }
}

close INFILE;
close INFILE2;
close OUTFILE;
close OUTFILE2;
exit;
