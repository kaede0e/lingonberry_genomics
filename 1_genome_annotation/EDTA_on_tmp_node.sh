#!/bin/bash

## EDTA - Extensive De novo TE annotator pipeline m##

#required files: genome.fasta, CDS.fasta, CDS.bed
#done in 2 days for lingonberry genome assembly (518Mbp genome, 12 chromosomes)

# ---------------------------------------------------------------------
##### WARNING: EDTA on whole genome creates 500k+ of tmp files. This script executes it in local disk in Cedar to avoid bloating storage space ######
module load StdEnv/2020 python
module load singularity
# ---------------------------------------------------------------------

cd $SLURM_TMPDIR
mkdir work
cd work
tar -xzf ~/projects/rrg-gowens/kaedeh/Lingonberry/output/annotation/TE_annotation_pipeline_3/my_data.tar.gz
cd my_data
ls #list everything that's there

# Now do my computations here on the local disk using the contents of the extracted archive...

export SINGULARITYENV_SLURM_TMPDIR=$SLURM_TMPDIR
singularity exec -B /home -B /project -B /scratch -B /localscratch -B $SLURM_TMPDIR -W $SLURM_TMPDIR /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/EDTA/EDTA-2.0.0.sif \
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/EDTA/EDTA.pl \
--genome Lingonberry_RedCandy_asm_7.2.ragtag.scaffold.chr.fasta \
--cds Lingonberry_RedCandy_CDS_anno.fasta \
--exclude Lingonberry_RedCandy_CDS_anno.bed \
--overwrite 1 --sensitive 1 --anno 1 --threads ${SLURM_CPUS_PER_TASK}
  #You can change --overwrite 0 to allow running from where you left off - when a job timed out.

# The computations are done, so clean up the data set...
cd $SLURM_TMPDIR
ls -R #recursive list of every file in local scratch (this will be long...)
tar -cf ~/projects/rrg-gowens/kaedeh/Lingonberry/output/annotation/TE_annotation_pipeline_3/EDTA_results.tar work


