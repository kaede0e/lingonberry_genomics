#!/bin/bash

#OrthoFinder for phylogenetic analysis with single-copy gene synteny
module load python
source /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/python_env/bin/activate

python /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/OrthoFinder_source/orthofinder.py -f reference_data

#--------------------------------------#
#in fossa
export PATH=$PATH:/project/ctb-grego/khirabayashi/bin/OrthoFinder
#Phylogenetic analysis using whole protein sequences
module load StdEnv/2020 python/3.7 gcc/9.3.0 blast+/2.13.0
source /project/ctb-grego/khirabayashi/bin/python_env/bin/activate
folder_with_protein_seqs_in_aminoacids=/project/ctb-grego/khirabayashi/reference_data/protein_seq/representatives/

cp /project/ctb-grego/khirabayashi/Lingonberry/out_Lingonberry_RedCandy_sup_model_duplexed_all_reads_combined/gene_annotation/Lingonberry_RedCandy_genes_anno.faa $folder_with_protein_seqs_in_aminoacids
orthofinder -t 50 -f $folder_with_protein_seqs_in_aminoacids
orthofinder -t 4 -d -f $folder_with_gene_seqs_in_nucleotides

#something is not working. 
export PATH=$PATH:/project/ctb-grego/khirabayashi/bin/OrthoFinder
module load StdEnv/2020 python/3.8 gcc/9.3.0 blast/2.2.26
virtualenv python_env_orthofinder
source python_env_orthofinder/bin/activate
pip install numpy
pip install scipy

#--------------------------------------#
#in cedar; use Orthofinder_source directory, and smcpp_python_env as virtualenv. 
module load python/3.8 StdEnv/2020 gcc/9.3.0
source /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/smcpp_python_env/bin/activate
export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/OrthoFinder_source

folder_with_protein_seqs_in_aminoacids=/home/kaedeh/projects/rrg-gowens/kaedeh/Lingonberry/reference_data/proteins_faa
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



