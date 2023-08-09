#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=125G
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --job-name=Lingonberry_RedCandy_guppy_basecall_duplex_GPU

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

## Post-sequencing process with Nanopore reads (duplex mode for R10.4)

module load gcc python/3.10 parasail
source /~/bin/python_env/bin/activate
export PATH=$PATH:/~/bin
export PATH=$PATH:/~/bin/ont-guppy/bin
module load cuda/11.7
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUDA_HOME/extras/CUPTI/lib64/

###### Duplex basecalling pipeline by guppy #########
#This will perform simplex basecall with `fast` model then use this to find duplex reads, then duplex basecall all the pairs that>
#Simplex basecalling with guppy `sup` model performed at last.

### 1. guppy runs 'fast' mode, find duplex pairs, do duplex basecalling for those.
guppy_duplex \
-i /~/raw_fast5_files/Lingonberry_RedCandy_SREXS_R10.4_3_Oct8_2022_fast5 \
-s ./output/guppy/Lingonberry_RedCandy_Takara_04_SREXS_03_Oct_8_2022_duplex_by_guppy \
--duplex_chunks_per_runner 240 #\
#--call_non_duplex_reads

gzip ./output/guppy/Lingonberry_RedCandy_Takara_04_SREXS_03_Oct_8_2022_duplex_by_guppy/final/duplex/pass/*

# ---------------------------------------------------------------------
echo "Duplex basecalled reads are found in output/guppy/*duplex_by_guppy/final/duplex/pass/"
# ---------------------------------------------------------------------

for file in `cat MinION_run_names.txt`;
do
  ### 2. remove duplex reads from 'sup' mode reads (do the singleplex basecalled fastq.gz file at once)
  zcat /~/output/guppy/${file}/pass/* | perl /~/scripts/filter_fastq.pl /~/output/guppy/${file}_duplex_by_guppy/final_pairs.txt > ${file}_duplexfiltered.sup.fastq
  gzip ${file}_duplexfiltered.sup.fastq

  # ---------------------------------------------------------------------
  echo "Duplex filtered singleplex 'sup' reads are found in output/guppy/*duplex_by_guppy/final/simplex_sup/"
  # ---------------------------------------------------------------------

  ### 3. combine duplex reads with 'sup' simplex reads into a single fastq.gz
  cat ${file}_duplexfiltered.sup.fastq.gz ./output/guppy/${file}_duplex_by_guppy/final/duplex/pass/*.fastq.gz > ${file}_duplex_basecalled.fastq.gz
  cp ${file}_duplex_basecalled.fastq.gz /~/output/sequencing_guppy/

  # ---------------------------------------------------------------------
  echo "Duplex reads + simplex 'sup' reads combined to SEQ_ID.duplex_basecalled.fastq.gz in /~/output/sequencing_guppy/"
  # ---------------------------------------------------------------------
done

# ---------------------------------------------------------------------
echo "Finished duplex basecalling with guppy at: `date`"
# ---------------------------------------------------------------------
echo ""
