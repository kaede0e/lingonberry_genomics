#!/user/bin/perl
use warnings;
use strict;


#This takes a fasta file of genes with the header format:
#>STRG.960.3.p1 STRG.960   Chr01_Vaccinium_myrtillus_NK2018_v1_RagTag:8573929-8592545(+)
# It picks the longest version of each gene and prints it out.

my %gene;
my %info;
my $current_gene = "NA";
my $last_title;
while(<STDIN>){
  chomp;
  if ($_ =~ m/^>/g){
    my @a = split(/\ /, $_);
    my $gene = $a[1];
    my $subgene = $a[0];
    if ($gene eq $current_gene){
      #Add it to the current
      $info{'name'}{$gene}{$subgene} = $_;
      $info{'length'}{$gene}{$subgene} = 0;
      $last_title = $subgene;
    }else{

      #It's a new gene part of the file
      #Print out best gene from last set
      #Find the longest length $subgene
      if ($current_gene ne "NA"){
        my $hash_ref = $info{'name'}{$current_gene};

        # Initialize variables to store the longest $subgene and its length
        my $longest_subgene = '';
        my $longest_length = 0;

        # Iterate through the keys of the inner hash
        foreach my $subgene (keys %$hash_ref) {

          # Get the length of the value for the current $subgene
          my $length = $info{'length'}{$current_gene}{$subgene};
          #print "gene is $subgene with $length AA\n";
          # Update the longest $subgene and its length if necessary
          if ($length > $longest_length) {
            $longest_subgene = $subgene;
            $longest_length = $length;
          }
        }
        #print "Longest gene is $longest_subgene at $longest_length\n";
        print "$info{'name'}{$current_gene}{$longest_subgene}\n";
        print "$info{'sequence'}{$current_gene}{$longest_subgene}\n";
      }
      #Load in new file
      $current_gene = $gene;
      $info{'name'}{$gene}{$subgene} = $_;
      $info{'length'}{$gene}{$subgene} = 0;
      $last_title = $subgene;
    }
  }else{
    #if its a part of the file with AA sequence
    my $length = length($_);
    $info{'length'}{$current_gene}{$last_title} += $length;
    $info{'sequence'}{$current_gene}{$last_title} .= $_;

  }


}
#Print last genes

my $hash_ref = $info{'name'}{$current_gene};

# Initialize variables to store the longest $subgene and its length
my $longest_subgene = '';
my $longest_length = 0;

# Iterate through the keys of the inner hash
foreach my $subgene (keys %$hash_ref) {

  # Get the length of the value for the current $subgene
  my $length = $info{'length'}{$current_gene}{$subgene};
  #print "gene is $subgene with $length AA\n";
  # Update the longest $subgene and its length if necessary
  if ($length > $longest_length) {
    $longest_subgene = $subgene;
    $longest_length = $length;
  }
}
#print "Longest gene is $longest_subgene at $longest_length\n";
print "$info{'name'}{$current_gene}{$longest_subgene}\n";
print "$info{'sequence'}{$current_gene}{$longest_subgene}\n";
