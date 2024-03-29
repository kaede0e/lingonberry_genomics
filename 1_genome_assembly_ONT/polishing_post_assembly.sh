#!/bin/bash

# The three methods below are not used. Skip to "SmartDenovo pipeline" for the actual pipeline used in the paper. 
# --------------------------------------------------------------------

##### Flye pipeline #####
# draft assembly with Flye
# generally 1 or 2 round of long-read polishing with Racon
# then 1 or 2 round of short-read (Illumina) polishing with Pilon

#minimap2 - mapping raw Nanpore reads to assembly
/~/bin/minimap2/minimap2 \ #check if ref.fa and qry.fa is opposite. I did ref=assembly
-x map-ont flye_assembly.fasta ../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz > minimap2_flye_assembly_to_raw_reads_overlap.paf
#Racon - polishing with long-reads
/~/bin/racon/build/bin/racon \
-u \
../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz \
minimap2_flye_assembly_to_raw_reads_overlap.paf \
flye_assembly.fasta > racon_ONT_polished_flye_assembly_round1.fasta

# --------------------------------------------------------------------

##### miniasm pipeline #####
# draft assembly with miniasm
# generally 3 rounds of long-read polishing with Racon which dramatically improves BUSCO score
# then 3 rounds of short-read (Illumina) polishing with Pilon which increases BUSCO from 85 to 95+

#minimap2 round 1 - align draft assembly to raw reads
/~/bin/minimap2/minimap2 \
-x map-ont miniasm_assembly.fasta ../../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz > minimap2_miniasm_assembly_to_raw_reads_overlap_round1.paf

#Racon round 1
/~/bin/racon/build/bin/racon \
-u ../../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz \
minimap2_miniasm_assembly_to_raw_reads_overlap_round1.paf \
miniasm_assembly.fasta > racon_ONT_polished_miniasm_assembly_round1.fasta

#minimap2 round 2 - align polished assembly round 1 to raw reads
/~/bin/minimap2/minimap2 \
-x map-ont racon_ONT_polished_miniasm_assembly_round1.fasta ../../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz > minimap2_miniasm_assembly_to_raw_reads_overlap_round2.paf

#Racon round 2
/~/bin/racon/build/bin/racon \
-u ../../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz \
minimap2_miniasm_assembly_to_raw_reads_overlap_round2.paf \
racon_ONT_polished_miniasm_assembly_round1.fasta > racon_ONT_polished_miniasm_assembly_round2.fasta

#minimap2 round 3 - align polished assembly round 2 to raw reads
/~/bin/minimap2/minimap2 \
-x map-ont racon_ONT_polished_miniasm_assembly_round2.fasta ../../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz > minimap2_miniasm_assembly_to_raw_reads_overlap_round3.paf

#Racon round 3
/~/bin/racon/build/bin/racon \
-u ../../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz \
minimap2_miniasm_assembly_to_raw_reads_overlap_round3.paf \
racon_ONT_polished_miniasm_assembly_round2.fasta > racon_ONT_polished_miniasm_assembly_round3.fasta

# --------------------------------------------------------------------

##### Wtdbg2 pipeline #####
# draft assembly with wtdbg2
# polishing with racon x3

#minimap2 round 1 - align draft assembly to raw reads
/~/bin/minimap2/minimap2 \
-x map-ont wtdbg2_assembly.fasta ../../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz > ../alignments/minimap2_wtdbg2_assembly_to_raw_reads_overlap.paf
#Racon round 1
/project/ctb-grego/khirabayashi/bin/racon/build/bin/racon \
-u \
../../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz \
../alignments/minimap2_wtdbg2_assembly_to_raw_reads_overlap.paf \
wtdbg2_assembly.fasta > ../polishing/racon_ONT_polished_wtdbg2_assembly_round1.fasta

#minimap2 round 2 - align polished assembly round 1 to raw reads
/~/bin/minimap2/minimap2 \
-x map-ont ../polishing/racon_ONT_polished_wtdbg2_assembly_round1.fasta ../../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz > ../alignments/minimap2_wtdbg2_assembly_to_raw_reads_overlap_round2.paf

#Racon round 2
/~/bin/racon/build/bin/racon \
-u ../../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz \
../alignments/minimap2_wtdbg2_assembly_to_raw_reads_overlap_round2.paf \
../polishing/racon_ONT_polished_wtdbg2_assembly_round1.fasta > ../polishing/racon_ONT_polished_wtdbg2_assembly_round2.fasta

#minimap2 round 3 - align polished assembly round 2 to raw reads
/~/bin/minimap2/minimap2 \
-x map-ont ../polishing/racon_ONT_polished_wtdbg2_assembly_round2.fasta ../../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz > ../alignments/minimap2_wtdbg2_assembly_to_raw_reads_overlap_round3.paf

#Racon round 3
/project/ctb-grego/khirabayashi/bin/racon/build/bin/racon \
-u ../../Lingonberry_RedCandy_MinION_SREXS_raw_reads.sup.fastq.gz \
../alignments/minimap2_wtdbg2_assembly_to_raw_reads_overlap_round3.paf \
../polishing/racon_ONT_polished_wtdbg2_assembly_round2.fasta > ../polishing/racon_ONT_polished_wtdbg2_assembly_round3.fasta

# --------------------------------------------------------------------

##### SmartDenovo pipeline #####
# draft assembly with SmartDenovo
# polishing with NextPolish x3 - uses python so activate python_env

#NextPolish x1 (make a long-read folder file (fofn), config file where the genome = draft assembly to polish)
ls reads1.fq reads2.fa.gz > lgs.fofn #path to your raw long-reads

#run.cfg
[General]
job_type = local
job_prefix = nextPolish
task = default
rewrite = yes
rerun = 3
parallel_jobs = 2
multithread_jobs = 3
genome = /~/smartdenovo_outputs/smartdenovo_assembly.fasta
genome_size = auto
workdir = ./01_rundir
polish_options = -p {multithread_jobs}

[lgs_option]
lgs_fofn = ./lgs.fofn
lgs_options = -min_read_len 1k -max_depth 100
lgs_minimap2_options = -x map-ont

#call the command as follows (also make sure to be in python_env)
./nextPolish output/run.cfg
#rename output: genome.nextpolish.fasta accordingly

# --------------------------------------------------------------------
##### Polishing with short-read with Pilon #####
module load pilon samtools bwa

#align reads to draft genome to improve
bwa index Lingonberry_RedCandy_smartdenovo_asm_4.fasta
bwa mem -t 40 Lingonberry_RedCandy_smartdenovo_asm_4.fasta ../../Lingonberry_RedCandy_F122990_1.2_paired_trimmomatic.fastq.gz ../../Lingonberry_RedCandy_F122990_2.2_paired_trimmomatic.fastq.gz > Lingonberry_RedCandy_smartdenovo_asm_4.2_aln-pe.sam
samtools sort -o Lingonberry_RedCandy_smartdenovo_asm_4.2_sorted.aln-pe.bam Lingonberry_RedCandy_smartdenovo_asm_4.2_aln-pe.sam
samtools index Lingonberry_RedCandy_smartdenovo_asm_4.2_sorted.aln-pe.bam
java -Xmx250g -jar /~/bin/pilon-1.24.jar --genome Lingonberry_RedCandy_smartdenovo_asm_4.fasta --frags Lingonberry_RedCandy_smartdenovo_asm_4.2_sorted.aln-pe.bam --output Lingonberry_RedCandy_smartdenovo_asm_5.2 --diploid &> Lingonberry_RedCandy_smartdenovo_asm_4.2to5.2_pilon.stdout.log

bwa index Lingonberry_RedCandy_smartdenovo_asm_5.2.fasta
bwa mem -t 40 Lingonberry_RedCandy_smartdenovo_asm_5.2.fasta ../../Lingonberry_RedCandy_F122990_1.2_paired_trimmomatic.fastq.gz ../../Lingonberry_RedCandy_F122990_2.2_paired_trimmomatic.fastq.gz > Lingonberry_RedCandy_smartdenovo_asm_5.2_aln-pe.sam
samtools sort -o Lingonberry_RedCandy_smartdenovo_asm_5.2_sorted.aln-pe.bam Lingonberry_RedCandy_smartdenovo_asm_5.2_aln-pe.sam
samtools index Lingonberry_RedCandy_smartdenovo_asm_5.2_sorted.aln-pe.bam
java -Xmx250g -jar /~/bin/pilon-1.24.jar --genome Lingonberry_RedCandy_smartdenovo_asm_5.2.fasta --frags Lingonberry_RedCandy_smartdenovo_asm_5.2_sorted.aln-pe.bam --output Lingonberry_RedCandy_smartdenovo_asm_6.2 --diploid &> Lingonberry_RedCandy_smartdenovo_asm_5.2to6.2_pilon.stdout.log

bwa index Lingonberry_RedCandy_smartdenovo_asm_6.2.fasta
bwa mem -t 40 Lingonberry_RedCandy_smartdenovo_asm_6.2.fasta ../../Lingonberry_RedCandy_F122990_1.2_paired_trimmomatic.fastq.gz ../../Lingonberry_RedCandy_F122990_2.2_paired_trimmomatic.fastq.gz > Lingonberry_RedCandy_smartdenovo_asm_6.2_aln-pe.sam
samtools sort -o Lingonberry_RedCandy_smartdenovo_asm_6.2_sorted.aln-pe.bam Lingonberry_RedCandy_smartdenovo_asm_6.2_aln-pe.sam
samtools index Lingonberry_RedCandy_smartdenovo_asm_6.2_sorted.aln-pe.bam
java -Xmx250g -jar /~/bin/pilon-1.24.jar --genome Lingonberry_RedCandy_smartdenovo_asm_6.2.fasta --frags Lingonberry_RedCandy_smartdenovo_asm_6.2_sorted.aln-pe.bam --output Lingonberry_RedCandy_smartdenovo_asm_7.2 --diploid &> Lingonberry_RedCandy_smartdenovo_asm_6.2to7.2_pilon.stdout.log


#align reads to draft genome to improve
bwa index Lingonberry_minus_canu_smartdenovo_asm_4.fasta
bwa mem Lingonberry_minus_canu_smartdenovo_asm_4.fasta ../../Lingonberry_minus_F122991_1_paired_trimmomatic.fastq.gz ../../Lingonberry_minus_F122991_2_paired_trimmomatic.fastq.gz > Lingonberry_minus_canu_smartdenovo_asm_4_aln-pe.sam
samtools sort -o Lingonberry_minus_canu_smartdenovo_asm_4_sorted.aln-pe.bam Lingonberry_minus_canu_smartdenovo_asm_4_aln-pe.sam
samtools index Lingonberry_minus_canu_smartdenovo_asm_4_sorted.aln-pe.bam
java -Xmx250G -jar /~/bin/pilon-1.24.jar --genome Lingonberry_minus_canu_smartdenovo_asm_4.fasta --frags Lingonberry_minus_canu_smartdenovo_asm_4_sorted.aln-pe.bam --output Lingonberry_minus_canu_smartdenovo_asm_5.fasta --diploid &> Lingonberry_minus_canu_smartdenovo_asm_4to5_pilon.stdout.log

bwa index Lingonberry_minus_canu_smartdenovo_asm_5.fasta
bwa mem Lingonberry_minus_canu_smartdenovo_asm_5.fasta ../../Lingonberry_minus_F122991_1_paired_trimmomatic.fastq.gz ../../Lingonberry_minus_F122991_2_paired_trimmomatic.fastq.gz > Lingonberry_minus_canu_smartdenovo_asm_5_aln-pe.sam
samtools sort -o Lingonberry_minus_canu_smartdenovo_asm_5_sorted.aln-pe.bam Lingonberry_minus_canu_smartdenovo_asm_5_aln-pe.sam
samtools index Lingonberry_minus_canu_smartdenovo_asm_5_sorted.aln-pe.bam
java -Xmx250G -jar /~/bin/pilon-1.24.jar --genome Lingonberry_minus_canu_smartdenovo_asm_5.fasta --frags Lingonberry_minus_canu_smartdenovo_asm_5_sorted.aln-pe.bam --output Lingonberry_minus_canu_smartdenovo_asm_6.fasta --diploid &> Lingonberry_minus_canu_smartdenovo_asm_5to6_pilon.stdout.log

bwa index Lingonberry_minus_canu_smartdenovo_asm_6.fasta
bwa mem Lingonberry_minus_canu_smartdenovo_asm_6.fasta ../../Lingonberry_minus_F122991_1_paired_trimmomatic.fastq.gz ../../Lingonberry_minus_F122991_2_paired_trimmomatic.fastq.gz > Lingonberry_minus_canu_smartdenovo_asm_6_aln-pe.sam
samtools sort -o Lingonberry_minus_canu_smartdenovo_asm_6_sorted.aln-pe.bam Lingonberry_minus_canu_smartdenovo_asm_6_aln-pe.sam
samtools index Lingonberry_minus_canu_smartdenovo_asm_6_sorted.aln-pe.bam
java -Xmx250G -jar /~/bin/pilon-1.24.jar --genome Lingonberry_minus_canu_smartdenovo_asm_6.fasta --frags Lingonberry_minus_canu_smartdenovo_asm_6_sorted.aln-pe.bam --output Lingonberry_minus_canu_smartdenovo_asm_7 --diploid &> Lingonberry_minus_canu_smartdenovo_asm_6to7_pilon.stdout.log

bwa index Lingonberry_minus_canu_smartdenovo_asm_7.fasta
bwa mem Lingonberry_minus_canu_smartdenovo_asm_7.fasta ../../Lingonberry_minus_F122991_1_paired_trimmomatic.fastq.gz ../../Lingonberry_minus_F122991_2_paired_trimmomatic.fastq.gz > Lingonberry_minus_canu_smartdenovo_asm_7_aln-pe.sam
samtools sort -o Lingonberry_minus_canu_smartdenovo_asm_6_sorted.aln-pe.bam Lingonberry_minus_canu_smartdenovo_asm_7_aln-pe.sam
samtools index Lingonberry_minus_canu_smartdenovo_asm_7_sorted.aln-pe.bam
java -Xmx250G -jar /~/bin/pilon-1.24.jar --genome Lingonberry_minus_canu_smartdenovo_asm_7.fasta --frags Lingonberry_minus_canu_smartdenovo_asm_7_sorted.aln-pe.bam --output Lingonberry_minus_canu_smartdenovo_asm_8 --diploid &> Lingonberry_minus_canu_smartdenovo_asm_7to8_pilon.stdout.log


# --------------------------------------------------------------------
#### Analyzing polished assembly reuslts with Merqury ####
#this calculates QV score (consensus accuracy)

export PATH=$PATH:/~/bin/meryl-1.4/bin
export PATH=$PATH:/~/bin/meryl-1.4/bin/merqury/eval/
module load StdEnv/2020 r/4.2.1 #package ggplot2 and scales should be already installed
#best_k.sh suggests kmer = 19.45 so I'll go with 19; make the meryl db 
#Merqury doesn't use paired reads, so it doesn't matter which single read I use (fwrd reads only used)
meryl k=19 count output Lingonberry_minus_F122991_1_paired_read.meryl ../../../Lingonberry_minus_F122991_1_paired_trimmomatic.fastq.gz
meryl k=19 count output Lingonberry_RedCandy_F122990_1.2_paired_read.meryl ../../../Lingonberry_RedCandy_F122990_1.2_paired_trimmomatic.fastq.gz
#calculates qv based on kmer distribution of Illumina reads
qv.sh Lingonberry_minus_F122991_1_paired_read.meryl ../Lingonberry_minus_canu_smartdenovo_asm_4.fasta Lingonberry_minus_canu_smartdenovo_asm_4
for n in {5..7};
do
	qv.sh Lingonberry_RedCandy_F122990_1.2_paired_read.meryl ../Lingonberry_RedCandy_smartdenovo_asm_${n}.2.fasta Lingonberry_RedCandy_smartdenovo_asm_${n}.2
done


 
