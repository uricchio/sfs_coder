#!/usr/bin/perl -w 
use strict;

unless(scalar(@ARGV) == 4){
    die "usage: $0 <input map file> <simMin> <simMax> <output file>\n";
}

my $infile = $ARGV[0];
my $beg = $ARGV[1];
my $end = $ARGV[2];

my %MAP = ();
my $minMap = -1;
my $maxMap = 0;
my $lastPoint = 0;
my $nPoints = 0;
my $less = 0;
my $bigger = 0;

open(IN, "$infile")
  or die "cannot read $infile\n";
my $line = <IN>;
my @PREV = ();
my @PREV2 = ();
my $PMAP = 0;
my $minMapPos = 0;
my $maxMapPos = 0;
while(<IN>){
  chomp;
  my @sp = split(' ',$_);
  if($sp[3] == 0 || $sp[3] <= $PMAP){ #skip first point with rate 0
    next;
  }
  if($sp[1] >= $beg && $sp[1] <= $end){ #physical position in interval
    if($minMap == -1){
      if(scalar(@PREV)>0){
	$minMap = $PREV[2]+($sp[3]-$PREV[2])*$beg/($sp[1]+$PREV[0]);
	$minMapPos = $sp[1];
      }
      else{
	$minMap = $sp[3]*$beg/$sp[1];
	$minMapPos = $sp[1];
      }
    }
    if($sp[3] > 0){
      $MAP{$sp[1]-$beg} = $sp[3];
      $lastPoint = $sp[1];
      $nPoints++;
      if($sp[3] > $maxMap){
	$maxMap = $sp[3];
	$maxMapPos = $sp[1];
      }
    }
  }
  elsif($bigger == 0 && $sp[1] >= $end){
    $nPoints++;
    $maxMap = $maxMap + ((($end-$lastPoint)/($sp[1]-$lastPoint))*
			 ($sp[3]-$maxMap));
    $maxMapPos = $sp[1];
    $bigger = 1;
    $lastPoint = $sp[1];
  }
  if(scalar(@PREV) > 0){
    $PREV2[0] = $PREV[0];
    $PREV2[1] = $PREV[1];
    $PREV2[2] = $PREV[2];
  }
  $PREV[0] = $sp[1];
  $PREV[1] = $sp[2];
  $PREV[2] = $sp[3];
  $PMAP = $sp[3];
  if($sp[1] > $end){
    last;
  }
}
if($lastPoint < $end){
  $nPoints++;
  $maxMap = $maxMap+($PREV[2]-$PREV2[2])/($PREV[0]-$PREV2[0])*($end-$PREV[0]);
  $maxMapPos = $PREV[0];
}
print "recMap = $minMap ($minMapPos) - $maxMap ($maxMapPos)\n";
close(IN);

open(OUT,">$ARGV[3]")
  or die "cannot write to $ARGV[3] file\n";
print OUT "$nPoints\n";

print "minMap=$minMap\nmaxMap=$maxMap\n";

foreach my $pos (sort {$a<=>$b} keys %MAP){
  my $val = ($MAP{$pos}-$minMap)/($maxMap-$minMap);
  print OUT "$pos $val\n";
}
print OUT ($end-$beg)," 1.0\n";
close(OUT);

