#!/bin/bash

#SBATCH --job-name=phage_assembly	## job name
#SBATCH -A <account>			## account to charge
#SBATCH -p standard			## partition/queue name
#SBATCH --nodes=1               	## number of nodes to use
#SBATCH --ntasks=1              	## number of tasks to launch
#SBATCH --cpus-per-task=16		## number of cores the job needs
#SBATCH --error=slurm_%A_%a.err		## error log file name: %A is job id, %a is array task id
#SBATCH --output=slurm_%A_%a.out	## output filename
#SBATCH --array=1-33			## number of array tasks

## load modules
module load unicycler/0.5.0

## define directories
input=data/processed/phage
output=data/processed/assembly

## create file name list
temp=$(basename -s _clean_1.fastq.gz $input/*_clean_1.fastq.gz | sort -u)

## select the file prefix corresponding to the array ID job number
prefix=`echo "$temp" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## job commands

# unicylcer assembly
unicycler \
-1 $input/${prefix}_clean_1.fastq.gz \
-2 $input/${prefix}_clean_2.fastq.gz \
-t $SLURM_CPUS_PER_TASK \
-o $output/${prefix}
