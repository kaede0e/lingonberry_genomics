#!/bin/bash
## PSMC - pairwise sequentially Markovian coalescent model ##

# ---------------------------------------------------------------------
export PATH=$PATH:/project/ctb-grego/khirabayashi/bin/psmc/utils
export PATH=$PATH:/project/ctb-grego/khirabayashi/bin/psmc
# ---------------------------------------------------------------------

module load samtools bcftools

#psmc from: https://github.com/lh3/psmc    easy to compile without permissions 
genome=/project/ctb-grego/khirabayashi/Lingonberry/out_Lingonberry_minus_sup_model_duplexed_all_reads_combined/scaffolding/Lingonberry_minus_asm_7.ragtag.scaffold.chr.fasta
bamdir=/project/ctb-grego/khirabayashi/Lingonberry/out_Lingonberry_minus_sup_model_duplexed_all_reads_combined/markdups

#bam prefix as positional argument eg for SAMPLE1.bam: sbatch -J SAMPLE1 PSMC.sh SAMPLE1  
RUN=Lingonberry_minus_asm_7.ragtag.scaffold.aln-pe.markdup
RUN2=Lingonberry_minus_asm_7.ragtag.scaffold.RedCandy_aln-pe.markdup.2
MU=3e-09 #mutation rate
GEN=5 #generation time

#ensure you have an older version of samtools: samtools0.1.19 WORKS CONFIRMED, you will also need 'gnuplot' for plotting
#set -d to 1/3 coverage and -D to 2x coverage (min and max)
#this creates consensus sequence
bcftools mpileup -q 20 -Q 20 -C 50 -Ou -f $genome $bamdir/${RUN}.bam | bcftools call -c -V indels \
| vcfutils.pl vcf2fq -d 12 -D 76 | gzip > ${RUN}_diploid.fq.gz
bcftools mpileup -q 20 -Q 20 -C 50 -Ou -f $genome $bamdir/${RUN2}.bam | bcftools call -c -V indels \
| vcfutils.pl vcf2fq -d 12 -D 74 | gzip > ${RUN2}_diploid.fq.gz

#creates psmc input 
fq2psmcfa -q20 ${RUN}_diploid.fq.gz > ${RUN}_diploid_5yr.psmcfa
fq2psmcfa -q20 ${RUN2}_diploid.fq.gz > ${RUN2}_diploid_5yr.psmcfa
#this runs psmc, largely default settings appropriate for almost everything 
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${RUN}_diploid_5yr.psmc ${RUN}_diploid_5yr.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${RUN2}_diploid_5yr.psmc ${RUN2}_diploid_5yr.psmcfa

#plot
#psmc_plot.pl -p -u $MU -g $GEN -Y 60 minus_psmc_yaxis_max60 Lingonberry_minus_bilberry_canu_smartdenovo_asm_7_ragtag.scaffold.aln-pe.markdup_diploid.psmc
psmc_plot.pl -R -p -u $MU -g $GEN -Y 60 Lingonberry_minus_5yr_psmc ${RUN}_diploid_5yr.psmc
psmc_plot.pl -R -p -u $MU -g $GEN Lingonberry_RedCandy_5yr_psmc_2 ${RUN2}_diploid_5yr.psmc
#the .txt file generated can be imported in R for plotting; 1st column = Time, 2nd column = Ne x 10^4 

###if you want 100 bootstraps... takes much longer 
splitfa ${RUN}_diploid.psmcfa > ${RUN}_diploid.split.psmcfa
seq 100 | xargs -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" \
          -o bootstrap100/round-{}.${RUN}_psmc ${RUN}_diploid.split.psmcfa | sh
cat ${RUN}_diploid.psmc bootstrap100/round-*.${RUN}_psmc > ${RUN}_combined.psmc
psmc_plot.pl -R -p -u $MU -g $GEN ${RUN} ${RUN}_combined.psmc


