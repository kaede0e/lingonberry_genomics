#!/bin/bash
## MSMC - multiple sequentially Markovian coalescent m##

# ---------------------------------------------------------------------
module load StdEnv/2020 gcc/9.3.0 gsl msmc2/2.1.3
module load python scipy-stack r 
export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/msmc-tools

#0.1 Align Illumina pe reads to my ref reference assembly (I can do this in rcs.uvic.ca)
#MUST ALIGN BOTH SPECIES TO THE SAME REFERENCE GENOME (chose ssp. minus for better assembly quality)
module load pilon samtools bwa
input_refgenome=Lingonberry_minus_asm_7.ragtag.scaffold
minus_Illumina_1=/project/ctb-grego/khirabayashi/Lingonberry/Lingonberry_minus_F122991_1_paired_trimmomatic.fastq.gz
minus_Illumina_2=/project/ctb-grego/khirabayashi/Lingonberry/Lingonberry_minus_F122991_2_paired_trimmomatic.fastq.gz
RedCandy_Illumina_1=/project/ctb-grego/khirabayashi/Lingonberry/Lingonberry_RedCandy_F122990_1.2_paired_trimmomatic.fastq.gz
RedCandy_Illumina_2=/project/ctb-grego/khirabayashi/Lingonberry/Lingonberry_RedCandy_F122990_2.2_paired_trimmomatic.fastq.gz

bwa index ${input_refgenome}.fasta
bwa mem -t 40 ${input_refgenome}.fasta $minus_Illumina_1 $minus_Illumina_2 > ${input_refgenome}.aln-pe.sam
samtools sort --threads 20 -O BAM -o ${input_refgenome}.aln-pe.bam ${input_refgenome}.aln-pe.sam
samtools index ${input_refgenome}.aln-pe.bam

bwa mem -t 40 ${input_refgenome}.fasta $RedCandy_Illumina_1 $RedCandy_Illumina_2 > ${input_refgenome}.RedCandy_aln-pe.2.sam
samtools sort --threads 20 -O BAM -o ${input_refgenome}.RedCandy_aln-pe.2.bam ${input_refgenome}.RedCandy_aln-pe.2.sam
samtools index ${input_refgenome}.RedCandy_aln-pe.2.bam

#0.2 Remove PCR duplicates with Picard (it's on GATK)
module load picard
#the best practice guideline is found in https://gatk.broadinstitute.org/hc/en-us/articles/360039568932
#0.2.1-3 (bams are already sorted so and others are unnecessary so moving on to step 4)
#0.2.4 Mark duplicates using MarkDuplicates
java -Xmx150g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
INPUT=${input_refgenome}.aln-pe.bam \
OUTPUT=${input_refgenome}.aln-pe.markdup.bam \
METRICS_FILE=metrics.txt \
CREATE_INDEX=true

java -Xmx150g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
INPUT=${input_refgenome}.RedCandy_aln-pe.2.bam \
OUTPUT=${input_refgenome}.RedCandy_aln-pe.markdup.2.bam \
METRICS_FILE=metrics.txt \
CREATE_INDEX=true

#0.3 the resulting .bam file should be free of duplicates and ready to use in msmc. 
# Calculate mean sequencing depth from the .bam file 
samtools depth -a ${input_refgenome}.aln-pe.markdup.bam | awk '{c++;s+=$3}END{print s/c}'
#minus: 35.327
samtools depth -a ${input_refgenome}.RedCandy_aln-pe.markdup.2.bam | awk '{c++;s+=$3}END{print s/c}'
#RedCandy: 37.3833

#I can't do this in rcs.uvic.ca so I had to transfer files over to cedar... 
#1. Make a VCF and mask file per sample/population for each chromosome. 
#do this in Cedar
for chr in `cat Lingonberry_master_inputs_scaff/Lingonberry_minus_ragtag.scaffold.chrnames.txt`;
do
  bcftools mpileup -q 20 -Q 20 -C 50 -Ou -r $chr -f Lingonberry_master_inputs_scaff/${input_refgenome}.chr.fasta Lingonberry_master_inputs_scaff/${input_refgenome}.aln-pe.markdup.bam | bcftools call -c -V indels |
  bamCaller.py 35 Lingonberry_master_inputs_scaff/Lingonberry_minus_F122991_paired_reads_${chr}.mask.bed.gz | gzip -c > Lingonberry_master_inputs_scaff/Lingonberry_minus_F122991_paired_reads_${chr}.vcf.gz;
done

for chr in `cat Lingonberry_master_inputs_scaff/Lingonberry_minus_ragtag.scaffold.chrnames.txt`;
do
  bcftools mpileup -q 20 -Q 20 -C 50 -Ou -r $chr -f Lingonberry_master_inputs_scaff/${input_refgenome}.chr.fasta Lingonberry_master_inputs_scaff/${input_refgenome}.RedCandy_aln-pe.markdup.2.bam | bcftools call -c -V indels |
  bamCaller.py 37 Lingonberry_master_inputs_scaff/Lingonberry_minus_F122990.2_paired_reads_RedCandy_${chr}.mask.bed.gz | gzip -c > Lingonberry_master_inputs_scaff/Lingonberry_minus_F122990.2_paired_reads_RedCandy_${chr}.vcf.gz;
done

# ---------------------------------------------------------------------
Job ID: 2804115
Cluster: cedar
User/Group: kaedeh/kaedeh
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 4
CPU Utilized: 06:46:03
CPU Efficiency: 46.65% of 14:30:28 core-walltime
Job Wall-clock time: 03:37:37
Memory Utilized: 95.34 MB
Memory Efficiency: 0.40% of 23.44 GB
# ---------------------------------------------------------------------


#2. Make a mask file for each species/population's genome using GenMap. 
#do this in rcs.uvic.ca
export PATH=$PATH:/project/ctb-grego/khirabayashi/bin/genmap-build/bin

genmap index -F /project/ctb-grego/khirabayashi/Lingonberry/out_Lingonberry_minus_sup_model_duplexed_all_reads_combined/scaffolding/${input_refgenome}.chr.fasta -I genmap_index_minus
genmap index -F /project/ctb-grego/khirabayashi/Lingonberry/out_Lingonberry_RedCandy_sup_model_duplexed_all_reads_combined/scaffolding/Lingonberry_RedCandy_bilberry_smartdenovo_asm_7_ragtag.scaffold.chr.fasta -I genmap_index_RedCandy_bilberry_refgenome

genmap map -K 30 -E 2 -I genmap_index_minus/ -O genmap_map_minus_output -t -w -bg -T 30 -v &> genmap_map_minus.log
cat genmap_map_minus_output.bedgraph | perl filter_mappability_bed.pl > excluded_regions.bed #this is regions to exclude because of repetitiveness. Score zero on the genmap.
cat genmap_map_minus_output.bedgraph | perl filter_mappability_nonrepetitive_bed.pl > included_regions.bed #this is regions to include for non-repetitive. Score >0 on the genmap.
genmap map -K 30 -E 2 -I genmap_index_RedCandy_bilberry_refgenome/ -O genmap_map_RedCandy_bilberry_refgenome_output -t -w -bg -T 30 -v &> genmap_map_Redcandy_bilberry_refgenome.log
cat genmap_out | perl filter_mappability_bed.pl > excluded_regions.bed #243708242 bp

#3. Combine the two populations in to prepare inputs for MSMC2. 
#MASTERVARDIR=/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/msmc-tools/MSMC-tutorial-files/cg_data

#generate_multihetsep.py --chr 1 \
#--mask test_output/NA12878.chr1.mask.bed.gz --mask test_output/NA19240.chr1.mask.bed.gz \
#--mask hs37d5_chr1.mask.bed \
#test_output/NA12878.chr1.vcf.gz test_output/NA19240.chr1.vcf.gz > test_output/NA12878_NA19240.chr1.multihetsep.txt

# ---------------------------------------------------------------------
### this is unnecessary ###
#OUTDIR=/home/kaedeh/scratch/Lingonberry/output/msmc2/Lingonberry_minus_outputs
#INDIR=/home/kaedeh/scratch/Lingonberry/output/msmc2/Lingonberry_minus_inputs

#for chr in `cat Lingonberry_ragtag.scaffold.chrnames.txt`;
#do
#	generate_multihetsep.py --chr $chr \
#	--mask $INDIR/Lingonberry_minus_canu_smartdenovo_asm_7_ragtag_F122991_paired_reads_${chr}.mask.bed.gz \
#    --mask $INDIR/Lingonberry_minus_excluded_regions.mask.bed \
#    $INDIR/Lingonberry_minus_canu_smartdenovo_asm_7_ragtag_F122991_paired_reads_${chr}.vcf.gz > $OUTDIR/Lingonberry_minus_canu_smartdenovo_asm_7_ragtag_F122991_paired_reads_${chr}.multihetsep.txt
#done

#INDIR=/home/kaedeh/scratch/Lingonberry/output/msmc2/Lingonberry_RedCandy_inputs
#OUTDIR=/home/kaedeh/scratch/Lingonberry/output/msmc2/Lingonberry_RedCandy_outputs

#for chr in `cat Lingonberry_ragtag.scaffold.chrnames.txt`;
#do
#	generate_multihetsep.py --chr $chr \
#	--mask $INDIR/Lingonberry_RedCandy_smartdenovo_asm_7_ragtag_F122990_paired_reads_${chr}.mask.bed.gz \
#    --mask $INDIR/Lingonberry_RedCandy_excluded_regions.mask.bed \
#    $INDIR/Lingonberry_RedCandy_smartdenovo_asm_7_ragtag_F122990_paired_reads_${chr}.vcf.gz > $OUTDIR/Lingonberry_RedCandy_smartdenovo_asm_7_ragtag_F122990_paired_reads_${chr}.multihetsep.txt
#done
# ---------------------------------------------------------------------

OUTDIR=/home/kaedeh/scratch/Lingonberry/output/msmc2/Lingonberry_master_outputs_scaff
INDIR=/home/kaedeh/scratch/Lingonberry/output/msmc2/Lingonberry_master_inputs_scaff
chrnames_prefix=_Vaccinium_vitis-idaea_ssp_minus

for chr in `cat Lingonberry_minus_ragtag.scaffold.chrnames.txt`;
do
	generate_multihetsep.py --chr $chr \
	--mask $INDIR/Lingonberry_minus_F122991_paired_reads_${chr}.mask.bed.gz \
	--mask $INDIR/Lingonberry_minus_F122990.2_paired_reads_RedCandy_${chr}.mask.bed.gz \
  --negative_mask $INDIR/Lingonberry_minus_excluded_regions.mask.bed \
    $INDIR/Lingonberry_minus_F122991_paired_reads_${chr}.vcf.gz \
    $INDIR/Lingonberry_minus_F122990.2_paired_reads_RedCandy_${chr}.vcf.gz > $OUTDIR/Lingonberry_minus_RedCandy_${chr}.multihetsep.txt
done

#4. Single population - effective population size estimate with MSMC2. 
#msmc2 -t 1 -p 1*2+15*1+1*2 -o test_output/NA12878.chr1.msmc2 -I 0,1 test_output/NA12878_NA19240.chr1.multihetsep.txt
#msmc2 -t 1 -p 1*2+15*1+1*2 -o test_output/NA19240.chr1.msmc2 -I 2,3 test_output/NA12878_NA19240.chr1.multihetsep.txt

msmc2 -t 1 -p 1*2+16*1+1*2 -o $OUTDIR/Lingonberry_minus.msmc2 -I 0,1 \
$OUTDIR/Lingonberry_minus_RedCandy_Chr01${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr02${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr03${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr04${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr05${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr06${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr07${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr08${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr09${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr10${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr11${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr12${chrnames_prefix}.multihetsep.txt

msmc2 -t 1 -p 1*2+16*1+1*2 -o $OUTDIR/Lingonberry_RedCandy.msmc2 -I 2,3 \
$OUTDIR/Lingonberry_minus_RedCandy_Chr01${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr02${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr03${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr04${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr05${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr06${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr07${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr08${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr09${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr10${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr11${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr12${chrnames_prefix}.multihetsep.txt

#5. Population separation history - estimate of split between two populations
#msmc2 -t 1 -I 0-2,0-3,1-2,1-3 -s -p 1*2+15*1+1*2 -o test_output/NA12878_NA19240.chr1.msmc2 test_output/NA12878_NA19240.chr1.multihetsep.txt

msmc2 -t 1 -I 0-2,0-3,1-2,1-3 -s -p 1*2+16*1+1*2 -o $OUTDIR/Lingonberry_minus_RedCandy.msmc2 \
$OUTDIR/Lingonberry_minus_RedCandy_Chr01${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr02${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr03${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr04${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr05${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr06${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr07${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr08${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr09${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr10${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr11${chrnames_prefix}.multihetsep.txt \
$OUTDIR/Lingonberry_minus_RedCandy_Chr12${chrnames_prefix}.multihetsep.txt

#6. Combine results for two population coalescence estimate
combineCrossCoal.py $OUTDIR/Lingonberry_minus_RedCandy.msmc2.final.txt $OUTDIR/Lingonberry_minus.msmc2.final.txt \
    $OUTDIR/Lingonberry_RedCandy.msmc2.final.txt > $OUTDIR/Lingonberry_minus_RedCandy.msmc2.combined.msmc2.final.txt

# ---------------------------------------------------------------------
Job ID: 2830955
Cluster: cedar
User/Group: kaedeh/kaedeh
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 4
CPU Utilized: 00:57:13
CPU Efficiency: 36.14% of 02:38:20 core-walltime
Job Wall-clock time: 00:39:35
Memory Utilized: 3.26 GB
Memory Efficiency: 13.90% of 23.44 GB
# ---------------------------------------------------------------------

#7. Estimate of homozygosity/heterozygosity to check if what I'm seeing is real. 
#filter variant calls (.vcf) with mask files 
module load vcftools
INDIR=raw_msmc_output
genmap_directory=/project/ctb-grego/khirabayashi/Lingonberry/out_genmap_mappabilitymask/genmap_results_minus

#1.get allele freq from unfiltered vcf
#2.extract filtered allele freq (mapped coverage + repetitive regions)
for i in {1..9}; 
do
  vcftools --freq2 --gzvcf $INDIR/Lingonberry_minus_F122991_paired_reads_Chr0${i}_Vaccinium_vitis-idaea_ssp_minus.vcf.gz --max-alleles 2 --out Lingonberry_minus_Chr0${i}_vcftools_allele;
  cat Lingonberry_minus_Chr0${i}_vcftools_allele.frq | awk '{print $1, $2, $2, $5}' | sed s/' '/'  '/g | tail -n +2 > Lingonberry_minus_Chr0${i}.vcf.allele.frq.bed;
  bedtools intersect -a Lingonberry_minus_Chr0${i}.vcf.allele.frq.bed -b $INDIR/Lingonberry_minus_F122991_paired_reads_Chr0${i}_Vaccinium_vitis-idaea_ssp_minus.mask.bed.gz > Lingonberry_minus_Chr0${i}.vcf.allele.frq.mapped.bed #mappable regions included
  bedtools intersect -a Lingonberry_minus_Chr0${i}.vcf.allele.frq.bed -b $genmap_directory/excluded_regions.bed -v > Lingonberry_minus_Chr0${i}.vcf.allele.frq.filtered.bed #repetitive regions excluded
done
for i in {10..12}; 
do
  vcftools --freq2 --gzvcf $INDIR/Lingonberry_minus_F122991_paired_reads_Chr${i}_Vaccinium_vitis-idaea_ssp_minus.vcf.gz --max-alleles 2 --out Lingonberry_minus_Chr${i}_vcftools_allele;
  cat Lingonberry_minus_Chr${i}_vcftools_allele.frq | awk '{print $1, $2, $2, $5}' | sed s/' '/'  '/g | tail -n +2 > Lingonberry_minus_Chr${i}.vcf.allele.frq.bed;
  bedtools intersect -a Lingonberry_minus_Chr${i}.vcf.allele.frq.bed -b $INDIR/Lingonberry_minus_F122991_paired_reads_Chr${i}_Vaccinium_vitis-idaea_ssp_minus.mask.bed.gz > Lingonberry_minus_Chr${i}.vcf.allele.frq.mapped.bed #mappable regions included
  bedtools intersect -a Lingonberry_minus_Chr${i}.vcf.allele.frq.bed -b $genmap_directory/excluded_regions.bed -v > Lingonberry_minus_Chr${i}.vcf.allele.frq.filtered.bed #repetitive regions excluded
done

for i in {1..9}; 
do
  vcftools --freq2 --gzvcf $INDIR/Lingonberry_minus_F122990.2_paired_reads_RedCandy_Chr0${i}_Vaccinium_vitis-idaea_ssp_minus.vcf.gz --max-alleles 2 --out Lingonberry_RedCandy_Chr0${i}_vcftools_allele;
  cat Lingonberry_RedCandy_Chr0${i}_vcftools_allele.frq | awk '{print $1, $2, $2, $5}' | sed s/' '/'  '/g | tail -n +2 > Lingonberry_RedCandy_Chr0${i}_vcftools_allele.frq.bed;
  bedtools intersect -a Lingonberry_RedCandy_Chr0${i}_vcftools_allele.frq.bed -b $INDIR/Lingonberry_minus_F122990.2_paired_reads_RedCandy_Chr0${i}_Vaccinium_vitis-idaea_ssp_minus.mask.bed.gz > Lingonberry_RedCandy_Chr0${i}.vcf.allele.frq.mapped.bed #mappable regions included
  bedtools intersect -a Lingonberry_RedCandy_Chr0${i}_vcftools_allele.frq.bed -b $genmap_directory/excluded_regions.bed -v > Lingonberry_RedCandy_Chr0${i}.vcf.allele.frq.filtered.bed #repetitive regions excluded
done
for i in {10..12}; 
do
  vcftools --freq2 --gzvcf $INDIR/Lingonberry_minus_F122990.2_paired_reads_RedCandy_Chr${i}_Vaccinium_vitis-idaea_ssp_minus.vcf.gz --max-alleles 2 --out Lingonberry_RedCandy_Chr${i}_vcftools_allele;
  cat Lingonberry_RedCandy_Chr${i}_vcftools_allele.frq | awk '{print $1, $2, $2, $5}' | sed s/' '/'  '/g | tail -n +2 > Lingonberry_RedCandy_Chr${i}_vcftools_allele.frq.bed;
  bedtools intersect -a Lingonberry_RedCandy_Chr${i}_vcftools_allele.frq.bed -b $INDIR/Lingonberry_minus_F122990.2_paired_reads_RedCandy_Chr${i}_Vaccinium_vitis-idaea_ssp_minus.mask.bed.gz > Lingonberry_RedCandy_Chr${i}.vcf.allele.frq.mapped.bed #mappable regions included
  bedtools intersect -a Lingonberry_RedCandy_Chr${i}_vcftools_allele.frq.bed -b $genmap_directory/excluded_regions.bed -v > Lingonberry_RedCandy_Chr${i}.vcf.allele.frq.filtered.bed #repetitive regions excluded
done

#cat Lingonberry_minus_*vcf.allele.frq.mapped.filtered.bed > Lingonberry_minus.vcf.allele.frq.mapped.filtered.bed
#cat Lingonberry_RedCandy_*vcf.allele.frq.mapped.filtered.bed > Lingonberry_RedCandy.vcf.allele.frq.mapped.filtered.bed
cat Lingonberry_minus_*vcf.allele.frq.filtered.bed > Lingonberry_minus.vcf.allele.frq.filtered.bed
cat Lingonberry_RedCandy_*vcf.allele.frq.filtered.bed > Lingonberry_RedCandy.vcf.allele.frq.filtered.bed

#3.similarly, how many callable sites are there within the window? 
cat included_regions.bed | grep "Chr01" > Lingonberry_minus_Chr01_included_regions.mask.bed
cat included_regions.bed | grep "Chr02" > Lingonberry_minus_Chr02_included_regions.mask.bed
cat included_regions.bed | grep "Chr03" > Lingonberry_minus_Chr03_included_regions.mask.bed
cat included_regions.bed | grep "Chr04" > Lingonberry_minus_Chr04_included_regions.mask.bed
cat included_regions.bed | grep "Chr05" > Lingonberry_minus_Chr05_included_regions.mask.bed
cat included_regions.bed | grep "Chr06" > Lingonberry_minus_Chr06_included_regions.mask.bed
cat included_regions.bed | grep "Chr07" > Lingonberry_minus_Chr07_included_regions.mask.bed
cat included_regions.bed | grep "Chr08" > Lingonberry_minus_Chr08_included_regions.mask.bed
cat included_regions.bed | grep "Chr09" > Lingonberry_minus_Chr09_included_regions.mask.bed
cat included_regions.bed | grep "Chr10" > Lingonberry_minus_Chr10_included_regions.mask.bed
cat included_regions.bed | grep "Chr11" > Lingonberry_minus_Chr11_included_regions.mask.bed
cat included_regions.bed | grep "Chr12" > Lingonberry_minus_Chr12_included_regions.mask.bed

for i in {1..9}; 
do
  bedtools intersect -a Lingonberry_minus_asm_7.ragtag.scaffold.chr.bed -b $INDIR/Lingonberry_minus_F122991_paired_reads_Chr0${i}_Vaccinium_vitis-idaea_ssp_minus.mask.bed.gz > Lingonberry_minus_Chr0${i}_callable_sites.mapped.bed #mappable regions included 
  bedtools intersect -a Lingonberry_minus_Chr0${i}_callable_sites.mapped.bed -b $genmap_directory/Lingonberry_minus_Chr0${i}_included_regions.mask.bed > Lingonberry_minus_Chr0${i}_callable_sites.mapped.filtered.bed #repetitive regions excluded 
done
for i in {10..12}; 
do
  bedtools intersect -a Lingonberry_minus_asm_7.ragtag.scaffold.chr.bed -b $INDIR/Lingonberry_minus_F122991_paired_reads_Chr${i}_Vaccinium_vitis-idaea_ssp_minus.mask.bed.gz > Lingonberry_minus_Chr${i}_callable_sites.mapped.bed #mappable regions included 
  bedtools intersect -a Lingonberry_minus_Chr${i}_callable_sites.mapped.bed -b $genmap_directory/Lingonberry_minus_Chr${i}_included_regions.mask.bed > Lingonberry_minus_Chr${i}_callable_sites.mapped.filtered.bed #repetitive regions excluded 
done

for i in {1..9}; 
do
  bedtools intersect -a Lingonberry_minus_asm_7.ragtag.scaffold.chr.bed -b $INDIR/Lingonberry_minus_F122990.2_paired_reads_RedCandy_Chr0${i}_Vaccinium_vitis-idaea_ssp_minus.mask.bed.gz > Lingonberry_RedCandy_Chr0${i}_callable_sites.mapped.bed
  bedtools intersect -a Lingonberry_RedCandy_Chr0${i}_callable_sites.mapped.bed -b $genmap_directory/Lingonberry_minus_Chr0${i}_included_regions.mask.bed > Lingonberry_RedCandy_Chr0${i}_callable_sites.mapped.filtered.bed
done
for i in {10..12}; 
do
  bedtools intersect -a Lingonberry_minus_asm_7.ragtag.scaffold.chr.bed -b $INDIR/Lingonberry_minus_F122990.2_paired_reads_RedCandy_Chr${i}_Vaccinium_vitis-idaea_ssp_minus.mask.bed.gz > Lingonberry_RedCandy_Chr${i}_callable_sites.mapped.bed
  bedtools intersect -a Lingonberry_RedCandy_Chr${i}_callable_sites.mapped.bed -b $genmap_directory/Lingonberry_minus_Chr${i}_included_regions.mask.bed > Lingonberry_RedCandy_Chr${i}_callable_sites.mapped.filtered.bed
done

cat Lingonberry_minus_*callable_sites.mapped.filtered.bed > Lingonberry_minus.callable_sites.mapped.filtered.bed
cat Lingonberry_RedCandy_*callable_sites.mapped.filtered.bed > Lingonberry_RedCandy.callable_sites.mapped.filtered.bed


#Use ROHan to compute run of homozygosity and inbreeding coefficient?
export PATH=$PATH:/project/ctb-grego/khirabayashi/bin/rohan/bin/
refgenome=/project/ctb-grego/khirabayashi/Lingonberry/out_Lingonberry_minus_sup_model_duplexed_all_reads_combined/markdups/Lingonberry_minus_bilberry_canu_smartdenovo_asm_7_ragtag.scaffold.fasta
minus_BAM=/project/ctb-grego/khirabayashi/Lingonberry/out_Lingonberry_minus_sup_model_duplexed_all_reads_combined/markdups/Lingonberry_minus_bilberry_canu_smartdenovo_asm_7_ragtag.scaffold.aln-pe.markdup.bam
RedCandy_BAM=/project/ctb-grego/khirabayashi/Lingonberry/out_Lingonberry_minus_sup_model_duplexed_all_reads_combined/markdups/Lingonberry_minus_bilberry_canu_smartdenovo_asm_7_ragtag.scaffold.RedCandy_aln-pe.markdup.2.bam

rohan -t 25 --tstv 1.71 --rohmu 2e-5 -o Lingonberry_minus_ROHan $refgenome $minus_BAM
rohan -t 25 --tstv 1.69 --rohmu 2e-5 -o Lingonberry_RedCandy_ROHan $refgenome $RedCandy_BAM



