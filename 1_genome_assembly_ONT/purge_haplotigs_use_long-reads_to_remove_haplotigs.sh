#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6000M

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

## Remove haplotigs with purge_haplotigs

#1. Map ONT raw reads to draft assemly 
module load StdEnv/2020 samtools bedtools minimap2 gcc/11.3.0 r/4.2.1 perl #load bedtools first (because it has incompatible programmes such as StdEnv2020, needs to be older StdEnv2016)
minimap2 -t 4 -ax map-ont \
/~/output/assembly/Lingonberry_RedCandy_smartdenovo_asm_4.fasta \
/~/output/sequencing_guppy/Lingonberry_RedCandy_all_reads.fastq.gz | 
samtools sort -m 6G -o Lingonberry_RedCandy_smartdenovo_asm_4_minimap2_ont.bam -T tmp.ali

export PATH=$PATH:/project/ctb-grego/khirabayashi/bin/minimap2
minimap2 -t 12 -ax map-ont \
/~/polishing/Lingonberry_RedCandy_smartdenovo_asm_7.fasta \
/~/Lingonberry_RedCandy_all_reads.fastq.gz | 
samtools sort -m 6G -o Lingonberry_RedCandy_smartdenovo_asm_7_minimap2_ont.bam -T tmp.ali

minimap2 -t 12 -ax map-ont \
/~/polishing/Lingonberry_minus_canu_smartdenovo_asm_7.fasta \
/~/Lingonberry_minus_all_reads.fastq.gz | 
samtools sort -m 6G -o Lingonberry_minus_canu_smartdenovo_asm_7_minimap2_ont.bam -T tmp.ali

minimap2 -t 20 -ax map-ont \
/~/polishing/Lingonberry_RedCandy_smartdenovo_asm_7.2.fasta \
/~/Lingonberry_RedCandy_all_reads.fastq.gz | 
samtools sort -m 6G -o Lingonberry_RedCandy_smartdenovo_asm_7.2_minimap2_ont.bam -T tmp.ali

#2. Visualize coverage plot and manually decide thresholds
export PATH=$PATH:/~/bin/purge_haplotigs/bin
export PATH=$PATH:/~/bin/purge_haplotigs/bin

purge_haplotigs  hist  -b Lingonberry_RedCandy_smartdenovo_asm_7.2_minimap2_ont.bam \
-g /~/polishing/Lingonberry_RedCandy_smartdenovo_asm_7.2.fasta \
-t 40

purge_haplotigs  hist  -b Lingonberry_minus_canu_smartdenovo_asm_7_minimap2_ont.bam \
-g /~/polishing/Lingonberry_minus_canu_smartdenovo_asm_7.fasta \
-t 40

purge_haplotigs  cov  -i Lingonberry_RedCandy_smartdenovo_asm_7.2_minimap2_ont.bam.gencov  -l 5 -m 42 -h 95 \
            -o Lingonberry_RedCandy_smartdenovo_asm_7.2_minimap2_ont_coverage_stats.csv -j 70  -s 70

purge_haplotigs  cov  -i Lingonberry_minus_canu_smartdenovo_asm_7_minimap2_ont.bam.gencov  -l 5 -m 40 -h 95 \
            -o Lingonberry_RedCandy_smartdenovo_asm_7_minimap2_ont_coverage_stats.csv -j 70  -s 70

purge_haplotigs  purge  -g /~/polishing/Lingonberry_RedCandy_smartdenovo_asm_7.2.fasta \
-c Lingonberry_RedCandy_smartdenovo_asm_7.2_minimap2_ont_coverage_stats.csv \
-o Lingonberry_RedCandy_smartdenovo_asm_7.2_purge_haplotig_curated \
-t 40 -d -b Lingonberry_RedCandy_smartdenovo_asm_7.2_minimap2_ont.bam #before running this, rename or delete the tmp_purge_haplotigs/ directory otherwise it won't run (tries to rerun from the previous outputs).

purge_haplotigs  purge  -g /~/polishing/Lingonberry_minus_canu_smartdenovo_asm_7.fasta \
-c Lingonberry_RedCandy_smartdenovo_asm_7_minimap2_ont_coverage_stats.csv \
-o Lingonberry_minus_canu_smartdenovo_asm_7_purge_haplotig_curated \
-t 40 -d -b Lingonberry_minus_canu_smartdenovo_asm_7_minimap2_ont.bam



