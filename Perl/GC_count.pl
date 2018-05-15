#! /usr/bin/perl
if (scalar(@ARGV) != 1) {
    print <<EOF;

USAGE: $0 <FASTA file>

EOF

exit;

}
open INFILE,   "$ARGV[0]"  or die "File '$ARGV[0]' open error!";
$GC_max = 0;
$start  = 1;
while (<INFILE>){
        chomp;
        if(/\>/){
              @token = split(">",$_);
              if($start==0){
              	$GC = sprintf "%.2f",(($G + $C)/$length)*100;
              	if ($GC > $GC_max){
              	        $GC_max = $GC;
              	        $ID_max = $ID;
              	}
	      }
	      $ID =$token[1];
              if($start==1){
                $ID_max = $ID;
              }
	      $start  = 0;	
	      $length = 0;
	      $A=$C=$G=$T=0;
        }else{
              $length += length $_;
              $A  +=($_=~tr/Aa//);
              $C  +=($_=~tr/Cc//);
              $G  +=($_=~tr/Gg//);
              $T  +=($_=~tr/Tt//);
        } 
}  
$GC = sprintf "%.2f",(($G + $C)/$length)*100;
        if ($GC > $GC_max){
              $GC_max = $GC;
              $ID_max = $ID;
        }
print "   [Sequence $ID_max has the highest GC%: $GC_max], \n";
close INFILE;
exit;
