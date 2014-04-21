#!/usr/bin/perl -w
use strict;

#input file assumed to be phastConse elements file 
## a wig file with 6 columns: bin, chr, start, stop, name, score

unless(scalar(@ARGV) == 7){
  print "usage: $0 <input CNC file> <overlap exons> <overlap UTRs> <output file> <chr> <beg pos> <end pos>\n";
  exit;
}

my $skipLen = 0;
my %SKIPregs = ();
open(IN,"$ARGV[1]") or die "cannot read exon file\n";
while(<IN>){
    chomp;
    my @sp = split(' ',$_);
    $SKIPregs{$sp[0]}{$sp[1]}++;
    $skipLen += $sp[1]-$sp[0]+1;
    if(($sp[1]-$sp[0]+1)%3 != 0 || $sp[0]>$sp[1]){
	die "EXON $sp[0]-$sp[1] not multple of 3!\n";
    }
}
close(IN);
open(IN,"$ARGV[2]") or die "cannot read UTR file\n";
while(<IN>){
    chomp;
    my @sp = split(' ',$_);
    foreach my $s (keys %SKIPregs){
	my $c=0;
	foreach my $e (keys %{$SKIPregs{$s}}){
	    $c++;
	    if($c>1){
		die "multiple regions starting at $s!!!\n";
	    }
	    if($sp[0]<=$e && $sp[1]>=$s){
		die "UTR overlaps CDS! $sp[0] $sp[1] and $s $e\n";
	    }
	}
    }
    $SKIPregs{$sp[0]}{$sp[1]}++;
    $skipLen += $sp[1]-$sp[0]+1;
    if($sp[0]>$sp[1]){
	die "UTR $sp[0] > $sp[1]!\n";
    }
}
close(IN);

my @SKIP = ();
my $nSKIP = 0;
foreach my $beg (sort {$a<=>$b} keys %SKIPregs){
    foreach my $end (sort {$a<=>$b} keys %{$SKIPregs{$beg}}){
	if($nSKIP == 0){
	    $SKIP[$nSKIP][0] = $beg;
	    $SKIP[$nSKIP][1] = $end;
	    $nSKIP++;
	}
	else{
	    if($beg <= $SKIP[$nSKIP-1][1]){
		die "[$beg,$end] in [$SKIP[$nSKIP-1][0],$SKIP[$nSKIP-1][1]]!\n";
	    }
	    else{
		$SKIP[$nSKIP][0] = $beg;
		$SKIP[$nSKIP][1] = $end;
		$nSKIP++;
	    }
	}
    }
}
print "read $nSKIP regions to be skipped with length $skipLen\n";

open(IN,"$ARGV[0]") or die "cannot read $ARGV[0]\n";
open(OUT,">$ARGV[3]") or die "cannot write to $ARGV[3]\n";

my $tmp = 0;
my $totLen = 0;
my $nCNCs = 0;
my $pb = 0;
my $pe = 0;
while(<IN>){
    chomp;
    my @sp = split(' ',$_);
    if($sp[1] ne "chr$ARGV[4]" || $sp[2] < $ARGV[5] || $sp[3] > $ARGV[6]){
	next;
    }
    my $beg = $sp[2];
    my $end = $sp[3];
    if($beg > $end){
	die "$beg > $end !!!\n";
    }
    elsif($beg<=$pe){
	die "overlapping elements!  $pb $pe and $beg $end\n";
    }
    $pb = $beg;
    $pe = $end;
    $tmp = scalar(@SKIP)-1 if($tmp>=scalar(@SKIP));
    while($tmp > 0 && $SKIP[$tmp][0] > $beg){
	$tmp--;
    }
    while($tmp < scalar(@SKIP) && $SKIP[$tmp][1] < $beg){
	$tmp++;
    }
    if($tmp<scalar(@SKIP) && $SKIP[$tmp][0]<=$beg && $SKIP[$tmp][1]>=$beg){
	$beg = $SKIP[$tmp][1]+1;
	$tmp++;
    }
    while($beg < $end && $tmp < scalar(@SKIP) && $SKIP[$tmp][0] <= $end){
	my $len = $SKIP[$tmp][0]-$beg;
	if($len > 0){
	    print OUT "$beg ".($SKIP[$tmp][0]-1)."\n";
	    $nCNCs++;
	    $totLen += $SKIP[$tmp][0]-$beg;
	}
	$beg = $SKIP[$tmp][1]+1;
	$tmp++;
    }
    if($beg < $end){
	print OUT "$beg $end\n";
    $nCNCs++;
    $totLen += $end-$beg+1;
  }
}

close(IN);
close(OUT);
print "$nCNCs CNCs with total length $totLen\n";
print "total selected length = ".($totLen+$skipLen)."\n";
