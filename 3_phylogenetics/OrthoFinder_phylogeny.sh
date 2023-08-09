#!/bin/bash

#--------------------------------------#
module load python/3.8 StdEnv/2020 gcc/9.3.0
source /~/bin/python_env/bin/activate
export PATH=$PATH:/~/bin/OrthoFinder_source

folder_with_protein_seqs_in_aminoacids=/~/reference_data/proteins_faa
output_directory=/home/kaedeh/scratch/Lingonberry/output/orthofinder

orthofinder.py -t 16 -a 4 -f $folder_with_protein_seqs_in_aminoacids -o $output_directory
#preliminary 10 species dataset interactive job with ntasks=1 thread died because of memory outage
#AND make sure to NOT blow up your project directory but accidentally saving all the tmp results in /project!

#--------------------------------------#
#make species tree with single-copy genes only. 
module load iq-tree
iqtree2 -p mafft_aln_fasta/ --prefix lingonberry_single_copy_gene_concat -B 1000 -T AUTO
iqtree2 -S mafft_aln_fasta/ --prefix lingonberry_single_copy_gene_loci -T AUTO
iqtree2 -t lingonberry_single_copy_gene_concat.treefile --gcf lingonberry_single_copy_gene_loci.treefile --prefix lingonberry_single_copy_gene_concord



