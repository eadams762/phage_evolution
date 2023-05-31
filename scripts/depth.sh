#!/bin/bash

#--------------------------SBATCH settings--------------------------

#SBATCH --job-name=depth		## job name
#SBATCH -A katrine_lab			## account to charge
#SBATCH -p standard			## partition/queue name
#SBATCH --nodes=1               	## number of nodes to use
#SBATCH --ntasks=1              	## number of tasks to launch
#SBATCH --cpus-per-task=1		## number of cores the job needs
#SBATCH --error=tmp/slurm_%A_%a.err	## error log file name: %A is job id, %a is array task id
#SBATCH --output=tmp/slurm_%A_%a.out	## output filename
#SBATCH --array=1-288			## number of array tasks

#-------------------------------------------------------------------

## load module
module load samtools/1.15.1

## define directories
dir=data/processed/coverage

## create file name list
temp=$(basename -s _sort.bam $dir/*_sort.bam | sort -u)

## select the file prefix corresponding to the array ID job number
file=`echo "$temp" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## job commands
samtools depth -a $dir/"$file"_sort.bam \
| awk -v f=$file '{sum+=$3; sumsq+=$3*$3} END { print f"\t"sum/NR"\t"sqrt(sumsq/NR - (sum/NR)**2)}' \
>> $dir/depth.txt
