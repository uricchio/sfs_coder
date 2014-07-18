#!/usr/bin/perl -w
use strict;

unless(scalar(@ARGV) == 7){
  print "usage:  $0 <input file> <output UTR file> <output exon file> <chr> <BEG POSITION> <END POSITION> <ref genome>\n";
  exit;
}

my $CHR = $ARGV[3];
my $BEG = $ARGV[4];
my $END = $ARGV[5];

open(IN,"$ARGV[0]") or die "cannot read $ARGV[0]\n";
if($ARGV[0] =~ /\.gz$/){
  close(IN);
  open(IN,"gunzip -c $ARGV[0] |") or die "cannot gunzip $ARGV[0]\n";
}

my $line = <IN>;
my %exons = ();
my $numExons = 0;
my %UTRs = ();
my $numUTRs = 0;
my $numGenes = 0;
my $totCodingLen = 0;
my $lines = 0;
while(<IN>){
  $lines++;
  chomp;
  my @sp = split(' ', $_);
  if($sp[2] ne "chr$CHR"){
    next;
  }
  elsif($sp[5]<$BEG){
    next;
  }
  elsif($sp[4]>$END){
    next;
  }
  my $CDSb = $sp[6];
  my $CDSe = $sp[7];
  my $nx = $sp[8];
  my $begs = $sp[9];
  my $ends = $sp[10];
  my @b = split(',',$begs);
  my @e = split(',',$ends);
  if(scalar(@b) != $nx || scalar(@e)!= $nx){
    die "b and e have different number of elements! ".(scalar(@b))." ".(scalar(@e))."\n\n$_\n";
  }
  
  for(my $i=0; $i<$nx; $i++){
    if($e[$i] < $CDSb){
      if(!exists($UTRs{$b[$i]})){
	$UTRs{$b[$i]} = $e[$i];
      }
    }
    elsif($b[$i] < $CDSb){
      if(!exists($UTRs{$b[$i]})){
	$UTRs{$b[$i]} = $CDSb-1;
	if($e[$i]>$CDSb+2){
	  $exons{$CDSb} = $e[$i];
	  $totCodingLen += $e[$i]-$CDSb+1;
	}
	if($CDSb>$e[$i]){
	  die "CDSb=$CDSb > e[$i]=$e[$i]\n";
	}
	$numExons++;
      }
    }
    elsif($b[$i] > $CDSe){
      if(!exists($UTRs{$b[$i]})){
	$UTRs{$b[$i]} = $e[$i];
      }
    }
    elsif($e[$i] > $CDSe){
      if(!exists($UTRs{$CDSe+1})){
	if($CDSe > $b[$i]+2){
	  $exons{$b[$i]} = $CDSe;
	  $totCodingLen += $CDSe-$b[$i]+1;
	}
	die "b[$i]=$b[$i] > CDSe=$CDSe\n" if $b[$i]>$CDSe;
	$numExons++;
	$UTRs{$CDSe+1} = $e[$i];
      }
    }
    elsif(!exists($exons{$b[$i]}) || $e[$i]>$exons{$b[$i]}){
      if($e[$i] > $b[$i]+2){
	$exons{$b[$i]} = $e[$i];
	$totCodingLen += $e[$i]-$b[$i]+1;
      }
      die "b[$i]=$b[$i] > e[$i]=$e[$i]\n" if $b[$i]>$e[$i];
      $numExons++;
    }
  }
  $numGenes++;
}
print "read $numGenes genes and $numExons exons with total length $totCodingLen\n";
close(IN);


#NOW SORT OUT UNIQ ELEMENTS AND ROUND LENGHT OF CDS EXONS TO WHOLE CODONS
my @srtEx = ();
my $nsrt = 0;
foreach my $s (sort {$a<=>$b} keys %exons){
  $srtEx[$nsrt][0] = $s;
  $srtEx[$nsrt][1] = $exons{$s};
  $nsrt++;
}
my @uniqEx = ();
$uniqEx[0][0] = $srtEx[0][0];
$uniqEx[0][1] = $srtEx[0][1];
if((($uniqEx[0][1]-$uniqEx[0][0]+1) % 3) == 1){
  $uniqEx[0][1] += 2;
}
elsif((($uniqEx[0][1]-$uniqEx[0][0]+1) % 3) == 2){
  $uniqEx[0][1]++;
}

my $nuex = 0;
for(my $i=1; $i<$nsrt; $i++){
  if($srtEx[$i][0] > $uniqEx[$nuex][1]){
    $nuex++;
    $uniqEx[$nuex][0] = $srtEx[$i][0];
    $uniqEx[$nuex][1] = $srtEx[$i][1];
  }
  elsif($srtEx[$i][1] > $uniqEx[$nuex][1]){
    $uniqEx[$nuex][1] = $srtEx[$i][1];
  }
  if((($uniqEx[$nuex][1]-$uniqEx[$nuex][0]+1) % 3) == 1){
    $uniqEx[$nuex][1] += 2;
  }
  elsif((($uniqEx[$nuex][1]-$uniqEx[$nuex][0]+1) % 3) == 2){
    $uniqEx[$nuex][1]++;
  }
}

open(OUTEX,">$ARGV[2]") or die "cannot write to $ARGV[2]\n";
my $finalLen = 0;
for(my $i=0; $i<=$nuex; $i++){
  $finalLen += $uniqEx[$i][1]-$uniqEx[$i][0]+1;
  my $seq = "";
  print OUTEX "$uniqEx[$i][0] $uniqEx[$i][1] $seq\n";
  if($uniqEx[$i][0] > $uniqEx[$i][1]){
    die "uniqEx[$i] $uniqEx[$i][0] > $uniqEx[$i][1]\n";
  }
}
print "$nuex uniq exons; length = $finalLen\n";
close(OUTEX);

my @srtUtr = ();
$nsrt = 0;
foreach my $s (sort {$a<=>$b} keys %UTRs){
  $srtUtr[$nsrt][0] = $s;
  $srtUtr[$nsrt][1] = $UTRs{$s};
  $nsrt++;
}
print "$nsrt sorted UTRs\n";
my @uniqUtr = ();
$uniqUtr[0][0] = $srtUtr[0][0];
$uniqUtr[0][1] = $srtUtr[0][1];

my $nutr = 0;
for(my $i=1; $i<$nsrt; $i++){
  if($srtUtr[$i][0] > $uniqUtr[$nutr][1]){
    $nutr++;
    $uniqUtr[$nutr][0] = $srtUtr[$i][0];
    $uniqUtr[$nutr][1] = $srtUtr[$i][1];
  }
  elsif($srtUtr[$i][1] > $uniqUtr[$nutr][1]){
    $uniqUtr[$nutr][1] = $srtUtr[$i][1];
  }
}
print "$nutr uniq UTRs\n";

my @UTRfin = ();
my $nUTRfin = 0;
for(my $i=0; $i<scalar(@uniqUtr); $i++){
  my $KEEP = 1;
  for(my $j=0; $j<scalar(@uniqEx); $j++){
    if($uniqEx[$j][0] > $uniqUtr[$i][1]){
      last;
    }
    elsif($uniqEx[$j][1] < $uniqUtr[$i][0]){
      next;
    }
    elsif($uniqEx[$j][0] <= $uniqUtr[$i][0] &&
	  $uniqEx[$j][1] >= $uniqUtr[$i][1]){
      $KEEP = 0;
      last;
    }
    elsif($uniqEx[$j][0] >= $uniqUtr[$i][0] &&
	  $uniqEx[$j][0] <= $uniqUtr[$i][1]){
      $uniqUtr[$i][1] = $uniqEx[$j][0]-1;
      if($uniqUtr[$i][1] <= $uniqUtr[$i][0]){
	$KEEP = 0;
	last;
      }
    }
    elsif($uniqEx[$j][1] >= $uniqUtr[$i][0] &&
	  $uniqEx[$j][1] <= $uniqUtr[$i][1]){
      $uniqUtr[$i][0] = $uniqEx[$j][1]+1;
      if($uniqUtr[$i][1] <= $uniqUtr[$i][0]){
	$KEEP = 0;
	last;
      }
    }
  }
  if($KEEP == 1 && $uniqUtr[$i][1] > $uniqUtr[$i][0]){
    $UTRfin[$nUTRfin][0] = $uniqUtr[$i][0];
    $UTRfin[$nUTRfin][1] = $uniqUtr[$i][1];
    #    print "$nUTRfin $UTRfin[$nUTRfin][0] $UTRfin[$nUTRfin][1]\n";
    $nUTRfin++;
    #    if($nUTRfin>10){
    #      last;
    #    }
  }
}

open(OUTUTR,">$ARGV[1]") or die "cannot write to $ARGV[1]\n";
$finalLen = 0;
for(my $i=0; $i<$nUTRfin; $i++){
  $finalLen += $UTRfin[$i][1]-$UTRfin[$i][0]+1;
  my $seq = "";
  print OUTUTR "$UTRfin[$i][0] $UTRfin[$i][1]\n";
}
print "$nUTRfin uniq utrs; length = $finalLen\n";
close(OUTUTR);
