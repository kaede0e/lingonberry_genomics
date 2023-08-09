#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=Lingonberry_RedCandy_assembly_ragtag_with_Vmac_refgenome

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

##### Ragtag - reference genome based scaffolding #####

module load StdEnv/2020 minimap2/2.24
source /~/busco_env/bin/activate

ragtag.py scaffold \
/home/kaedeh/scratch/Lingonberry/reference_data/V_myrtillus_genome.fasta \
~/output/assembly/Lingonberry_RedCandy_asm_7_purge_haplotig.fasta \
-u -e /home/kaedeh/scratch/Lingonberry/reference_data/V_myrtillus_genome_scaffold_names.txt
