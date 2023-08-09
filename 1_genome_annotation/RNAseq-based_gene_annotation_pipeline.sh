#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12000M

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

## RNAseq evidence-based gene annotation pipeline for novel genome assembly ##

# ---------------------------------------------------------------------

# Description of this pipeline:
#1. Map RNAseq/cDNA seq to draft genome assembly (HISAT2)
#2. Assemble transcript, retaining hits with BLAST alignment (StringTie)
#3. Find gene structure including gene, mRNA, CDS, UTR, etc. within transcripts (TransDecoder)
#4. Annotate genes with orthologs (EggNOG-mapper)

# ---------------------------------------------------------------------

#####################################
#### Execution of programmes ########
#####################################

export PATH=$PATH:/home/kaedeh/scratch/Lingonberry/reference_data/RNAseq/
export PATH=$PATH:/home/kaedeh/projects/rrg-gowens/kaedeh/Lingonberry/raw_fastq_files/ #store my Illumina reads here
module load hisat2 stringtie samtools bedtools
export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/
export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/util/ #this is not working so you need full path to invoke scripts
module load StdEnv/2020 gcc/9.3.0 blast+/2.13.0

# ---------------------------------------------------------------------

#### 0. QC RNA reads & check all files are there ########
# Make sure to run fastqc/0.11.9 or fastp before using the data.
# Make sure to copy one final draft genome to working directory /home/kaedeh/scratch/Lingonberry/output/gene_annotation_pipeline in 'Lingonberry_RedCandy_*asm.fasta' format.
# Make sure to have all RNAseq .fastq files (paired end, stranded) in /home/kaedeh/projects/rrg-gowens/kaedeh/Lingonberry/raw_fastq_files/*_1.fastq or *_2.fastq, matching prefix for paired reads.

# ---------------------------------------------------------------------

input_refgenome=/home/kaedeh/projects/rrg-gowens/kaedeh/Lingonberry/output/scaffold/Lingonberry_RedCandy_asm_7.2.ragtag.scaffold.fasta
#### 1. Transcript alignment to genome ####
hisat2-build $input_refgenome Lingonberry_RedCandy_hisat2_index #took about 20min
ls /home/kaedeh/projects/rrg-gowens/kaedeh/Lingonberry/raw_fastq_files/*_1*.fastq.gz | tr '\n' ',' | sed 's/.$//' > RNAseq_m1_list_of_files.txt
ls /home/kaedeh/projects/rrg-gowens/kaedeh/Lingonberry/raw_fastq_files/*_2*.fastq.gz | tr '\n' ',' | sed 's/.$//' > RNAseq_m2_list_of_files.txt
hisat2 \
-q -x Lingonberry_RedCandy_hisat2_index -1 `cat RNAseq_m1_list_of_files.txt` -2 `cat RNAseq_m2_list_of_files.txt` -S RedCandy_RNAseq_on_genome_hisat2_aln.sam
samtools sort -o RedCandy_RNAseq_on_genome_hisat2_sorted.aln.bam RedCandy_RNAseq_on_genome_hisat2_aln.sam

#echo "Finished Hisat2 alignment: `date`"

# ---------------------------------------------------------------------

#### 2. Transcript assembly and extract fasta sequence into fasta format ####
stringtie -o RedCandy_stringtie_gene_structure_output.gtf RedCandy_RNAseq_on_genome_hisat2_sorted.aln.bam #assembles the aligned transcripts into a file with structural definitions
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/EVidenceModeler-v2.0.0/EvmUtils/misc/cufflinks_gtf_to_alignment_gff3.pl \
RedCandy_stringtie_gene_structure_output.gtf > RedCandy_stringtie_transcript_asm_output.gff #this changes structure file into gff

#echo "Finished StringTie transcript assembly: `date`"

# ---------------------------------------------------------------------

#### 3. Find CDS within transcripts ####
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl \
RedCandy_stringtie_gene_structure_output.gtf $input_refgenome > RedCandy_stringtie_transcripts.fasta #extracts fasta seq of assembled transcript, in StringTie transcript coordinate (i.e. STR1.1)
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl \
RedCandy_stringtie_gene_structure_output.gtf > RedCandy_stringtie_transcripts.gff3 #still in transcript coordinate
TransDecoder.LongOrfs -S -t RedCandy_stringtie_transcripts.fasta #-S adds the "stranded" seq mode so that it doesn't mess up the orientation.

#echo "Finished TransDecoder long ORF prediction: `date`"

# BlastP library is prepared from UniProt Arabidopsis and Vaccinium known proteins#
blastp -query RedCandy_stringtie_transcripts.fasta.transdecoder_dir/longest_orfs.pep \
    -db blastp_lib/uniprot_Vaccinium_Arabidopsis_known_proteins \
    -outfmt 6 -evalue 1e-5 -num_threads 8 -out RedCandy_blastp_Vaccinium_Arabidopsis_homology_based_genes.txt

#echo "Finished BlastP gene search for orthologous genes in Vaccinium/Arabidopsis: `date`"

TransDecoder.Predict -t RedCandy_stringtie_transcripts.fasta --retain_blastp_hits RedCandy_blastp_Vaccinium_Arabidopsis_homology_based_genes.txt
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/util/cdna_alignment_orf_to_genome_orf.pl \
     RedCandy_stringtie_transcripts.fasta.transdecoder.gff3 \
     RedCandy_stringtie_transcripts.gff3 \
     RedCandy_stringtie_transcripts.fasta > RedCandy_stringtie_transcripts.fasta.transdecoder.genome.gff3 #save this final .gff3 file and make sure it is accessible to you.

awk -F ";Name" '{print $1}' RedCandy_stringtie_transcripts.fasta.transdecoder.genome.gff3 > RedCandy_stringtie_transcripts.fasta.transdecoder.genome.clean.gff3 #cleans up the annotation file with long headline.

echo "Finished TransDecoder CDS prediction: `date`"

mkdir final_annotation_files
scp $input_refgenome final_annotation_files/ #genome assembly
scp RedCandy_stringtie_transcripts.fasta.transdecoder.genome.clean.gff3 final_annotation_files/Lingonberry_RedCandy_genes_anno.gff3 #gene positions and features annotated on the above genome

echo "Finished TranDecoder genome annotation with genes: `date`"

# ---------------------------------------------------------------------

#### 4. Gene functional annotation on EggNOG-mapper web ####
# get a protein.fasta from annotation (.gff3 mapped to genome)
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/util/gff3_file_to_bed.pl RedCandy_stringtie_transcripts.fasta.transdecoder.genome.clean.gff3 > RedCandy_stringtie_transcripts.fasta.transdecoder.genome.clean.bed
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/util/gff3_file_to_proteins.pl \
--gff3 final_annotation_files/Lingonberry_RedCandy_genes_anno.gff3 \
--fasta $input_refgenome > final_annotation_files/Lingonberry_RedCandy_genes_anno.faa
# take the longest gene from isoforms/splicing variants
cat final_annotation_files/Lingonberry_RedCandy_genes_anno.faa | perl /home/kaedeh/projects/rrg-gowens/kaedeh/Lingonberry/scripts/fasta2longestisoform.pl > final_annotation_files/Lingonberry_RedCandy_longest_genes_anno.faa

echo "Finished converting nt .fasta into proteins .fasta: `date`"
echo "Congrats, your final_annotation_files/*anno.faa is ready for submission to EggNOG-mapper!: `date`"

# ---------------------------------------------------------------------

#### 5. prepare files for EDTA 
cd final_annotation_files
cat Lingonberry_RedCandy_genes_anno.gff3 | grep "gene" > Lingonberry_RedCandy_genes_only_anno.gff3
cat Lingonberry_RedCandy_genes_anno.gff3 | grep "CDS" > Lingonberry_RedCandy_CDS_anno.gff3
module load bedops 
gff2bed < Lingonberry_RedCandy_CDS_anno.gff3 > Lingonberry_RedCandy_CDS_anno.bed #inputs to EDTA
module load bedtools
bedtools getfasta Lingonberry_RedCandy_asm_7.2.ragtag.scaffold.fasta -bed Lingonberry_RedCandy_CDS_anno.bed > Lingonberry_RedCandy_CDS_anno.fasta
mkdir my_data
cp Lingonberry_RedCandy_asm_7.2.ragtag.scaffold.fasta Lingonberry_RedCandy_CDS_anno.bed Lingonberry_RedCandy_CDS_anno.fasta my_data
tar -xzvf my_data.tar.gz my_data

echo "Copy the my_data.tar.gz to input for EDTA pipeline; TE annotation."

#### 6. get basic stats on your gene annotation
echo "The number of genes annotated: `cat Lingonberry_RedCandy_genes_anno.gff3 | grep "gene" | wc -l`. "
echo "The total length of genes annotated: `cat Lingonberry_RedCandy_genes_anno.gff3 | grep "gene" | awk '{print $5-$4}' | paste -s -d+ - | bc` bp. "

# ---------------------------------------------------------------------


Job ID: 2750996
Cluster: cedar
User/Group: kaedeh/kaedeh
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 8
CPU Utilized: 1-20:15:46
CPU Efficiency: 23.91% of 7-17:08:56 core-walltime
Job Wall-clock time: 23:08:37
Memory Utilized: 4.98 GB
Memory Efficiency: 5.32% of 93.75 GB

