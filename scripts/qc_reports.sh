#!/bin/bash

#--------------------------SBATCH settings------

#SBATCH --job-name=reports		## job name
#SBATCH -A <account>			## account name
#SBATCH -p standard			## partition/queue name
#SBATCH --nodes=1			## number of nodes to use
#SBATCH --ntasks=1			## number of tasks to launch
#SBATCH --cpus-per-task=32		## number of cores the job needs
#SBATCH --error=tmp/slurm-%J.err	## error log file
#SBATCH --output=tmp/slurm-%J.out	## output log file

#-----------------------------------------------

# load modules
module load fastqc/0.11.9
module load anaconda/2020.07

# define directories
raw=data/raw
clean=data/processed/reads
mqc=output/reports/qc
fqc=data/fastqc

# fastqc reports
fastqc \
-o $fqc/raw/ \
--noextract \
-t $SLURM_CPUS_PER_TASK \
$raw/*.fastq.gz

fastqc \
-o $fqc/clean/ \
--noextract \
-t $SLURM_CPUS_PER_TASK \
$clean/*.fastq.gz

# multiqc reports
source activate multiqc
multiqc $fqc/raw/ -o $mqc/raw/
multiqc $fqc/clean/ -o $mqc/clean/
