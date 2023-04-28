#!/bin/bash

#SBATCH --job-name=phage_assembly	## job name
#SBATCH -A katrine_lab			## account to charge
#SBATCH -p standard			## partition/queue name
#SBATCH --nodes=1               	## number of nodes to use
#SBATCH --ntasks=1              	## number of tasks to launch
#SBATCH --cpus-per-task=32		## number of cores the job needs
#SBATCH --error=slurm-%J.err		## error log file
#SBATCH --output=slurm-%J.out		##output info file

## load modules
module load unicycler/0.5.0

## define directories
input=data/processed/phage
output=data/processed/assembly

# unicylcer assembly
unicycler \
-1 $input/Phage_Car_rep3_sub_1.fastq.gz \
-2 $input/Phage_Car_rep3_sub_2.fastq.gz \
-t $SLURM_CPUS_PER_TASK \
-o $output/Phage_Car_sub \
--keep 0
