#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=125G
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --job-name=Lingonberry_RedCandy_guppy_basecall_GPU

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

## Post-sequencing process with Nanopore reads

#### Basecalling ########

/~/bin/ont-guppy/bin/guppy_basecaller \
-c dna_r9.4.1_450bps_sup.cfg \
-i ./fast5 \
-s ./output/guppy/Lingonberry_SREXS_01_May_05_2022 \
--compress_fastq \
-x "auto" \
--num_callers 14 \
--gpu_runners_per_device 8 \
--chunk_size 500 \
--chunks_per_runner 240
#--resume

# ---------------------------------------------------------------------
echo "Finished basecalling at: `date`"
# ---------------------------------------------------------------------
echo ""

##Outputs in /pass/fastq have xxxx(dataname).fastq
## /pass  contains reads with above Q7 quality read score value, /fail is below.
## so first thing you'd want to do after is to concatenate all and gzip.
cat ./fastq.gz/* > ./fastq/Lingonberry.fastq
gzip ./fastq/Lingonberry.fastq
