#!/bin/perl

use strict;
use warnings;

my $pairs = $ARGV[0];


open PAIRS, $pairs;

my %hash;
while(<PAIRS>){
  chomp;
  my @a = split(/\ /, $_);
  my $first_code = $a[0];
  my $second_code = $a[1];
  $hash{$first_code} = 1;
  $hash{$second_code} = 1;
}
close PAIRS;
my $line_number = 0;
my $print_line = 0;
while(<STDIN>) {
  $line_number++;
  chomp;
  if ($line_number % 4 == 1){
    my @a = split(/\ /, $_);
    my $code = $a[0];
    $code =~ s/@//g;
    if ($hash{$code}){
      $print_line = 0;
    }else{
      $print_line = 1;
    }
  }
  if ($print_line == 1){
    print "$_\n";
  }
}

