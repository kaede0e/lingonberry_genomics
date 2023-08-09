#!/bin/bash
#SBATCH --time=14-00:00:00
#SBATCH --mem=192000M
#SBATCH --cpus-per-task=48

#Assembly with SmartDenovo
#Includes mapping, trimming, and bunch of cleaning and consensus finding

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

export PATH=$PATH:/~/bin
export PATH=$PATH:/~/bin/smartdenovo/
export PATH=$PATH:/~/bin/smartdenovo/smartdenovo/

#All in one-go
smartdenovo.pl \
-c 1 -t 20 -p Lingonberry_RedCandy_smartdenovo_asm_1 /~/Lingonberry_RedCandy_all_reads.fastq.gz > Lingonberry_RedCandy_smartdenovo_asm_1.mak
make -f Lingonberry_RedCandy_smartdenovo_asm_1.mak #you have to call this 'make file' from //bin/smartdenovo/

smartdenovo.pl \
-c 1 -t 50 -p Lingonberry_RedCandy_canu_corrected_smartdenovo_asm_1 /~/Lingonberry_RedCandy_canu_correction.trimmedReads.fasta.gz > Lingonberry_RedCandy_canu_corrected_smartdenovo_asm_1.mak
make -f Lingonberry_RedCandy_canu_corrected_smartdenovo_asm_1.mak #you have to call this 'make file' from //bin/smartdenovo/


# ---------------------------------------------------------------------
echo "Done assembly with smartdenovo."
# ---------------------------------------------------------------------

echo "Finished job at `date`"
