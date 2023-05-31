#!/bin/bash

#--------------------------SBATCH settings--------------------------

#SBATCH --job-name=breseq_S_I		## job name
#SBATCH -A <account>			## account to charge
#SBATCH -p standard			## partition/queue name
#SBATCH --nodes=1               	## number of nodes to use
#SBATCH --ntasks=1              	## number of tasks to launch
#SBATCH --cpus-per-task=32		## number of cores the job needs
#SBATCH --error=tmp/slurm_%A_%a.err	## error log file name: %A is job id, %a is array task id
#SBATCH --output=tmp/slurm_%A_%a.out	## output filename
#SBATCH --array=1-24			## number of array tasks

#-------------------------------------------------------------------

## load module
module load anaconda/2020.07

## activate breseq environment
source activate breseq

## define directories
input=data/processed/sister_I
bill=data/ref/Bill.gbk
bob=data/ref/Bob.gbk
car=data/ref/Car.gbk
ccs4=data/ref/CCS4.gbk
sds1=data/ref/SDS1.gbk
output=data/processed/breseq/sister_I

## create file name list
temp=$(cat $input/prefix.txt)

## select the file prefix corresponding to the array ID job number
prefix=`echo "$temp" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## job commands

# run breseq
breseq -p \
-j $SLURM_CPUS_PER_TASK \
--polymorphism-reject-indel-homopolymer-length 0 \
--polymorphism-reject-surrounding-homopolymer-length 0 \
--polymorphism-bias-cutoff 0 \
--polymorphism-minimum-total-coverage-each-strand 0 \
-o $output/${prefix} \
-r $bill \
-r $bob \
-r $car \
$input/${prefix}_clean_1.fastq.gz \
$input/${prefix}_clean_2.fastq.gz

# remove unnecessary output files
rm -rf $output/${prefix}/01_sequence_conversion/
rm -rf $output/${prefix}/02_reference_alignment/
rm -rf $output/${prefix}/03_candidate_junctions/
rm -rf $output/${prefix}/04_candidate_junction_alignment/
rm -rf $output/${prefix}/05_alignment_correction/
rm -rf $output/${prefix}/06_bam/

# rename output.gd file
mv $output/${prefix}/output/output.gd $output/${prefix}/output/${prefix}.gd
