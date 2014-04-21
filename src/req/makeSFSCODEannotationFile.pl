#!/usr/bin/perl -w
use strict;


# Alterations made by Lawrence Uricchio on 3/14/2014
# The hash checks on %REGS were not working correctly 
# because of heterogeneity in the data stored (some array,
# some not arrays).  Should be fixed now

unless(scalar(@ARGV) == 11){
    die "usage:  $0 <MINPOS> <MAXPOS> <NEUTRAL BUFFER> <OUTFILE> <GENES FILE> <UTRs FILE> <CNCs FILE> <DENSE BEG> <DENSE END> <SEL/NEUT> <WITHSEQ>\n";
}

my $BLOCKSIZE = 5000;
my $MINPOS = $ARGV[0];
my $MAXPOS = $ARGV[1];
my $BUFF = $ARGV[2];
my $DB = $ARGV[7];
my $DE = $ARGV[8];
my $WITHSEQ=$ARGV[10];
if($WITHSEQ ne 0){
  die "sorry, please set WITHSEQ to 0, as this is not yet implemented...\n";
}
unless($WITHSEQ == 0 || $WITHSEQ == 1){
  die "arg 10 \"WITHSEQ\" must be 0 or 1\n";
}

my %REGS = ();
my $nRegs = 0;
my @CREGS = ();
my $nCREGS = 0;
open(IN,"$ARGV[4]") or die "cannot read genes file: $ARGV[4]\n";
my $simMin = 1e10;
my $simMax = 0;
my $lines = 0;
while(<IN>){
    $lines++;
    chomp;
    my @sp = split(' ',$_);
    if($sp[0] < $MINPOS || $sp[0] > $MAXPOS){
	next;
    }
    $REGS{$sp[0]}{$sp[1]-$sp[0]+1}[0] = "E";
    $REGS{$sp[0]}{$sp[1]-$sp[0]+1}[1] = $sp[2];
    $CREGS[$nCREGS][0] = $sp[0];
    $CREGS[$nCREGS][1] = $sp[1]-$sp[0]+1;
    $nCREGS++;
    if($sp[0] > $sp[1]){
	die "E $lines $sp[0] > $sp[1]\n";
    }
    $nRegs++;
}
close(IN);

open(IN,"$ARGV[5]") or die "cannot read UTR file $ARGV[5]\n";
$lines = 0;
while(<IN>){
    $lines++;
    chomp;
    my @sp = split(' ',$_);
    if($sp[0] < $MINPOS || $sp[0] > $MAXPOS){
	next;
    }
    if(exists($REGS{$sp[0]}) && exists($REGS{$sp[0]}{$sp[1]-$sp[0]+1})){
      next;
    }
    $REGS{$sp[0]}{$sp[1]-$sp[0]+1}[0] = "U";
    if(exists($sp[2])) {
        $REGS{$sp[0]}{$sp[1]-$sp[0]+1}[1] = $sp[2];
    }
    if($sp[0] > $sp[1]){
	die "U $lines $sp[0] > $sp[1]\n";
    }
    $nRegs++;
}
close(IN);
open(IN,"$ARGV[6]") or die "cannot read CNCs file $ARGV[6]\n";
$lines = 0;
while(<IN>){
  $lines++;
  chomp;
  my @sp = split(' ',$_);
  if($sp[0] < $MINPOS || $sp[0] > $MAXPOS){
    next;
  }
  if(exists($REGS{$sp[0]}) && exists($REGS{$sp[0]}{$sp[1]-$sp[0]+1})){
    next;
  }
  $REGS{$sp[0]}{$sp[1]-$sp[0]+1}[0] = "C";
  if(exists($sp[2])){
      $REGS{$sp[0]}{$sp[1]-$sp[0]+1}[1] = $sp[2];
  }
  if($sp[0] > $sp[1]){
    die "C $lines $sp[0] > $sp[1]\n";
  }
  $nRegs++;
}
close(IN);
print "read $nRegs regions\n";

my @srtREGS = ();
my $nsrtREGS = 0;
foreach my $s (sort {$a<=>$b} keys %REGS){
  foreach my $l (keys %{$REGS{$s}}){
    $srtREGS[$nsrtREGS][0] = $s;
    $srtREGS[$nsrtREGS][1] = $l;
    $srtREGS[$nsrtREGS][2] = $REGS{$s}{$l}[0];
    $srtREGS[$nsrtREGS][3] = $REGS{$s}{$l}[1];
    if($nsrtREGS>0 && $srtREGS[$nsrtREGS][0] < $srtREGS[$nsrtREGS-1][0]+$srtREGS[$nsrtREGS-1][1]){
      die "overlapping elements!!  $srtREGS[$nsrtREGS][0]+$srtREGS[$nsrtREGS][1] and $srtREGS[$nsrtREGS-1][0]+$srtREGS[$nsrtREGS-1][1]\n";
    }
    $nsrtREGS++;
  }
}

#foreach my $s (keys %REGS) {
#     foreach my $l (keys %{$REGS{$s}}) {
#
#         print "$s $l @{$REGS{$s}{$l}}\n";
#     }
#}
#exit();

#add neutral buffer around coding regions
my @NEUT = ();
my $nNEUT = 0;
my $tmp = 0;
if(1){
    while($tmp<$nCREGS && $CREGS[$tmp][0]+$CREGS[$tmp][1] < $DB){
	$tmp++;
    }
    $NEUT[0][0] = $CREGS[$tmp][0]-$BUFF;
    if($NEUT[0][0] < 0){
	die "need to go back and adjust values to make upstream buffer non-negative!\n";
    }
    $NEUT[0][1] = $BLOCKSIZE;
    $nNEUT++;
    my $tsum = $BLOCKSIZE;
    while($tsum<=$BUFF-$BLOCKSIZE){
	$NEUT[$nNEUT][0] = $NEUT[$nNEUT-1][0]+$NEUT[$nNEUT-1][1];
	$NEUT[$nNEUT][1] = $BLOCKSIZE;
	$nNEUT++;
	$tsum += $BLOCKSIZE;
    }
    if($tsum < $BUFF){
	$NEUT[$nNEUT][0] = $NEUT[$nNEUT-1][0]+$NEUT[$nNEUT-1][1];
	$NEUT[$nNEUT][1] = $BUFF-$tsum;
	$nNEUT++;
    }
    my $neutLen = $BUFF; #before first element
    my $numLess = 0;
    $tmp++;
    while($tmp<$nCREGS && $CREGS[$tmp][0]<$DE){
	my $dist = $CREGS[$tmp][0]-($CREGS[$tmp-1][0]+$CREGS[$tmp-1][1]);
	if($dist >= 2*$BUFF){
	    $neutLen += 2*$BUFF;
	    $NEUT[$nNEUT][0] = $CREGS[$tmp-1][0]+$CREGS[$tmp-1][1];
	    $NEUT[$nNEUT][1] = $BLOCKSIZE;
	    $nNEUT++;
	    my $tsum = $BLOCKSIZE;
	    while($tsum<$BUFF-$BLOCKSIZE){
		$NEUT[$nNEUT][0] = $NEUT[$nNEUT-1][0]+$NEUT[$nNEUT-1][1];
		$NEUT[$nNEUT][1] = $BLOCKSIZE;
		$nNEUT++;
		$tsum += $BLOCKSIZE;
	    }
	    if($tsum < $BUFF){
		$NEUT[$nNEUT][0] = $NEUT[$nNEUT-1][0]+$NEUT[$nNEUT-1][1];
		$NEUT[$nNEUT][1] = $BUFF-$tsum;
		$nNEUT++;
	    }
	    
	    $NEUT[$nNEUT][0] = $CREGS[$tmp][0]-$BUFF;
	    $NEUT[$nNEUT][1] = $BLOCKSIZE;
	    $nNEUT++;
	    $tsum = $BLOCKSIZE;
	    while($tsum<$BUFF-$BLOCKSIZE){
		$NEUT[$nNEUT][0] = $NEUT[$nNEUT-1][0]+$NEUT[$nNEUT-1][1];
		$NEUT[$nNEUT][1] = $BLOCKSIZE;
		$nNEUT++;
		$tsum += $BLOCKSIZE;
	    }
	    if($tsum < $BUFF){
		$NEUT[$nNEUT][0] = $NEUT[$nNEUT-1][0]+$NEUT[$nNEUT-1][1];
		$NEUT[$nNEUT][1] = $BUFF-$tsum;
		$nNEUT++;
	    }
	}
	elsif($dist > 0){
	    $neutLen += $dist;
	    $numLess++;
	    if($dist <= $BLOCKSIZE){
		$NEUT[$nNEUT][0] = $CREGS[$tmp-1][0]+$CREGS[$tmp-1][1];
		$NEUT[$nNEUT][1] = $dist;
		$nNEUT++;
	    }
	    else{
		$NEUT[$nNEUT][0] = $CREGS[$tmp-1][0]+$CREGS[$tmp-1][1];
		$NEUT[$nNEUT][1] = $BLOCKSIZE;
		$nNEUT++;
		$tsum = $BLOCKSIZE;
		while($tsum<$dist-$BLOCKSIZE){
		    $NEUT[$nNEUT][0] = $NEUT[$nNEUT-1][0]+$NEUT[$nNEUT-1][1];
		    $NEUT[$nNEUT][1] = $BLOCKSIZE;
		    $tsum += $BLOCKSIZE;
		    $nNEUT++;
		}
		if($tsum < $dist){
		    $NEUT[$nNEUT][0] = $NEUT[$nNEUT-1][0]+$NEUT[$nNEUT-1][1];
		    $NEUT[$nNEUT][1] = $dist-$tsum;
		    $nNEUT++;
		}
	    }
	}
	$tmp++;
    }
    if($tmp<$nCREGS){
	$neutLen += $BUFF;  #after last element
	$tsum = $BLOCKSIZE;
	$NEUT[$nNEUT][0] = $CREGS[$tmp-1][0]+$CREGS[$tmp-1][1];
	$NEUT[$nNEUT][1] = $BLOCKSIZE;
	$nNEUT++;
	while($tsum<$BUFF-$BLOCKSIZE){
	    $NEUT[$nNEUT][0] = $NEUT[$nNEUT-1][0]+$NEUT[$nNEUT-1][1];
	    $NEUT[$nNEUT][1] = $BLOCKSIZE;
	    $tsum += $BLOCKSIZE;
	    $nNEUT++;
	}
	if($tsum < $BUFF){
	    $NEUT[$nNEUT][0] = $NEUT[$nNEUT-1][0]+$NEUT[$nNEUT-1][1];
	    $NEUT[$nNEUT][1] = $BUFF-$tsum;
	    $nNEUT++;
	}
    }
}

my $si = 0;
my $beg = 0;
for(my $i=0; $i<$nNEUT; $i++){
    $beg = $NEUT[$i][0];

    while($si>0 && $srtREGS[$si][0]+$srtREGS[$si][1] >= $NEUT[$i][0]){
	$si--;
    }
    while($si<$nsrtREGS && $srtREGS[$si][0]+$srtREGS[$si][1] <= $NEUT[$i][0]){
	$si++;
    }
    if($si>=$nsrtREGS || $srtREGS[$si][0] >= $NEUT[$i][0]+$NEUT[$i][1]){
	my $d = ($NEUT[$i][0]+$NEUT[$i][1])-$beg;
	if($d>0){
	    $REGS{$beg}{$d}[0] = "N";
	}
	next;
    }
    else{
	while($si<$nsrtREGS && $srtREGS[$si][0]<=$NEUT[$i][0]+$NEUT[$i][1]){
	    my $d = $srtREGS[$si][0]-$beg;
	    if($d>0){
		$REGS{$beg}{$d}[0] = "N";
	    }
	    $beg = $srtREGS[$si][0]+$srtREGS[$si][1];
	    $si++;
	}
	my $d = ($NEUT[$i][0]+$NEUT[$i][1])-$beg;
	if($d>0){
	    $REGS{$beg}{$d}[0] = "N";
	}
    }
}

#now loop through REGS again, and count up non-simulated blocks.
my $cnt = 0;
my $prev = 0;
my @FIN = ();
my $nFIN = 0;
my $nLOCI = 0;
my $sumSim = 0;
my $finNEUT = 0;
my $finNEUTlen = 0;
my $TDIST = 0;
my @REM = ();
my $pend = 0;
foreach my $s (sort {$a<=>$b} keys %REGS){
    my $foo=0;
    foreach my $l (keys %{$REGS{$s}}){
	$simMin = $s if($s<$simMin);
	$simMax = $s+$l if($s+$l>$simMax);
	$foo++;
	if($s < $pend){
	    die "$s < $pend!! $REGS{$s}{$l}: $REM[$nFIN-1][0] $REM[$nFIN-1][1] $FIN[$nFIN-1]\n";
	}
	$pend = $s+$l;
	if($foo>1){
	    die "REGS{$s} has multiple lengths!!\n";
	}
	if($cnt == 0){
	    $TDIST = $s;
	}
	if($cnt>0){
	    my $d = $s-$prev;
	    if($d>0){
		$FIN[$nFIN] = "$d;";
		$REM[$nFIN][0] = $s;
		$REM[$nFIN][1] = $d;
		$nFIN++;
		$TDIST += $d;
	    }
	}
	if($REGS{$s}{$l}[0] eq "E"){
	    if($ARGV[9] eq "NEUT"){
		$FIN[$nFIN] = "$l,,C,0,1.0;";
	    }
	    elsif($ARGV[9] eq "SEL"){
		$FIN[$nFIN] = "$l,,C,2 0 0 0 0.184 0.00040244,1.0;";
	    }
	    else{
		die "$ARGV[9] must be either SEL or NEUT\n";
	    }
	    $REM[$nFIN][0] = $s;
	    $REM[$nFIN][1] = $l;
	    $nFIN++;
	    if($TDIST != $s){
		warn "0 FIN[$nFIN] = $TDIST, s=$s!!\n";
		for(my $j=$nFIN; $j>$nFIN-10; $j--){
		    warn "$j $FIN[$j]; $REM[$j][0]+$REM[$j][1]\n";
		}
		exit;
	    }
	    $TDIST += $l;
	}
	elsif($REGS{$s}{$l}[0] eq "N"){
	    $FIN[$nFIN] = "$l,,N,0,1.0;";
	    $REM[$nFIN][0] = $s;
	    $REM[$nFIN][1] = $l;
	    $nFIN++;
	    $finNEUT++;
	    $finNEUTlen += $l;
	    if($TDIST != $s){
		warn "1 FIN[$nFIN] = $TDIST, s=$s!!\n";
		for(my $j=$nFIN; $j>$nFIN-10; $j--){
		    warn "$j $FIN[$j]; $REM[$j][0]+$REM[$j][1]\n";
		}
		exit;
	    }
	    $TDIST += $l;
	}
	else{
	    if($ARGV[9] eq "NEUT"){
		$FIN[$nFIN] = "$l,,N,0,1.0;";
	    }
	    elsif($ARGV[9] eq "SEL"){
		$FIN[$nFIN] = "$l,,N,2 0 0 0 0.0415 0.0015625,1.0;";
	    }
	    else{
		die "$ARGV[9] must be either SEL or NEUT\n";
	    }
	    $REM[$nFIN][0] = $s;
	    $REM[$nFIN][1] = $l;
	    $nFIN++;
	    if($TDIST != $s){
		warn "2 FIN[$nFIN] = $TDIST, s=$s!!\n";
		for(my $j=$nFIN; $j>$nFIN-10; $j--){
		    warn "$j $FIN[$j]; $REM[$j][0]+$REM[$j][1]\n";
		}
		exit;
	    }
	    $TDIST += $l;
	}
	$sumSim += $l;
	$nLOCI++;
	$prev = $s+$l;
	$cnt++;
    }
}

open(OUT,">$ARGV[3]") or die "cannot write to $ARGV[3]\n";
print OUT "$nLOCI\n";
for(my $i=0; $i<$nFIN; $i++){
    print OUT "$FIN[$i]\n";
}
close(OUT);

print "total simulated (excluding overlap) = $sumSim\n";
print "total neutral regions = $finNEUT; len = $finNEUTlen\n";
print "simulating bases:\n$simMin $simMax\n";
