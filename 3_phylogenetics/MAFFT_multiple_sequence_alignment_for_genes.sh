#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1200M

## make species tree with single-copy genes only ##

# ---------------------------------------------------------------------

#1. Align each gene using multiple sequence aligner. 
## MAFFT - Multiple sequence aligner ##

module load mafft 
for fasta in `ls Single_Copy_Orthologue_Sequences/OG*.fa`;
do
  mafft --globalpair --maxiterate 1000 $fasta > ${fasta}_mafft_aln.fa;
done

#2. Change the header to represent each species. The headers need to be the same among all the genes to run IQTREE correctly.
module load seqkit/2.3.1  StdEnv/2020
#mkdir iqtree_input_aln
for fasta in `ls iqtree_input_aln/OG*.fa_mafft_aln.fa`;
do
  cat $fasta | grep ">" | sed s/'>'//g > ${fasta}.header.txt;
  paste -d '  ' ${fasta}.header.txt input_species.txt > ${fasta}.replace_header_names.txt;
  seqkit replace -p "(.+)" -r '{kv}' -k ${fasta}.replace_header_names.txt ${fasta} > ${fasta}.header.replaced.fa;
done

#3. IQ-Tree to infer species tree. 

## IQ-TREE - species tree inference and concordance factor from gene trees ##

# ---------------------------------------------------------------------

module load iq-tree #may need to module --force purge and then reload
iqtree2 -p iqtree_input_aln/ --prefix lingonberry_single_copy_gene_concat -B 1000
iqtree2 -S iqtree_input_aln/ --prefix lingonberry_single_copy_gene_loci
iqtree2 -t lingonberry_single_copy_gene_concat.treefile --gcf lingonberry_single_copy_gene_loci.treefile --prefix lingonberry_single_copy_gene_concord


## make species tree with BUSCO genes ##

# ---------------------------------------------------------------------

#Extract BUSCO outputs
for species in `cat list_of_species.txt`;
do
  cp -r ${species}_BUSCO_output/run_eudicots_odb10/busco_sequences/single_copy_busco_sequences ${species}_single_copy_busco_genes/;
done 
#Do this for each species from the subdirectory to change the header
for species in `cat ../species.txt`;
do
  for fasta in `ls *.faa`
  do
    cat $fasta | sed "1s/.*/>${species%.faa}/" > header_changed_${fasta};
  done
done
#concatenate the same busco genes 
ls */header_changed_* | cut -f 2 -d "/" | sort | uniq > list_of_uniq_busco_genes_all.txt
for busco_gene in `cat list_of_uniq_busco_genes_all.txt`;
do
  cat */$busco_gene > single_copy_busco_genes/combined_${busco_gene}
done
#concatenate only those busco genes shared by all species in the dataset
ls */header_changed_* | cut -f 2 -d "/" | sort | uniq -c | grep "11 header" | sed s/'     11 '/''/g > list_of_uniq_shared_busco_genes_all_2.txt
for busco_gene in `cat list_of_uniq_shared_busco_genes_all_2.txt`;
do
  cat */$busco_gene > shared_single_copy_busco_genes_2/combined_${busco_gene}
done

#1. Align every busco gene using MAFFT
module load mafft 
mkdir iqtree_input_aln
for fasta in `ls *.faa`;
do
  mafft --globalpair --maxiterate 1000 $fasta > iqtree_input_aln/${fasta}.mafft.aln.faa;
done

#2. IQ-Tree to build indivdiual trees. Should probaly do this as a job.
module load iq-tree #may need to module --force purge and then reload
## the mafft aln with less than three species won't be accepted. 
## IQ-Tree only takes alignment file of min. 3 species. 
N=40 #this script allows the following loop to be run in parallel, of N threads, each loop uses one thread. 
for fasta in `ls iqtree_input_aln/*.mafft.aln.faa`;
do
  ((i=i%N)); ((i++==0)) && wait
  name=`echo $fasta | cut -d/ -f2`
  iqtree -s $fasta -pre ${name}_gene_loci -nt 1 &
done
mv single_copy_busco_genes/*.treefile single_copy_busco_genes_trees
cat single_copy_busco_genes_trees/*.treefile > single_copy_busco_genes_trees/combined_busco_genes.11sp.trees

#2.2 for use in MEGAX, concatenate alignment files into one .meg file with header:
#------------ text file shold have a header with --------
#mega
!Title busco1;
!Format DataType=Protein indel=-;
#--------------------------------------------------------
module load StdEnv/2020
module load seqkit/2.3.1
seqkit concat iqtree_input_aln/*.faa > mega_input_aln/combined_header_changed_shared_busco_genes.aln.faa
cat combined_header_changed_shared_busco_genes.aln.faa | sed s/'>'/'#'/g > combined_header_changed_shared_busco_genes.aln.meg #this file should be ready for MEGAX TimeTree/RelTime analysis. 

#3. filter gene trees 
tree_folder=/~/out_busco/single_copy_busco_genes_trees
module load python/3.7 r/4.2.1
source /~/bin/python_env/bin/activate
export PATH=$PATH:/project/ctb-grego/khirabayashi/bin/TreeShrink
run_treeshrink.py -t single_copy_busco_genes_trees/combined_busco_genes.11sp.trees \
-o single_copy_busco_genes_trees/ \
-O combined_busco_genes.11sp.trees.treeshrink

#3. Astral to build species tree from individual gene trees
module load java StdEnv/2020
export PATH=$PATH:/~/bin/Astral
java -jar -Xmx25000M /~/bin/Astral/astral.5.7.8.jar \
-i combined_busco_genes.11sp.trees.treeshrink.trees \
-o combined_busco_genes.11sp.trees.treeshrink.astral.out.tre \
2> species_tree_astral.out.log

#4. Calculate gene concordance and get individual gene trees, with Astral species tree.
iqtree2 -t combined_busco_genes.11sp.trees.treeshrink.astral.out.tre --gcf combined_busco_genes.11sp.trees.treeshrink.trees --prefix lingonberry_busco_genes_concord --cf-verbose







