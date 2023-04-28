#!/bin/bash

#--------------------------SBATCH settings------

#SBATCH --job-name=gdtools        ## job name
#SBATCH -A katrine_lab          ## account name
#SBATCH -p standard             ## partition/queue name
#SBATCH --nodes=1               ## number of nodes to use
#SBATCH --ntasks=1              ## number of tasks to launch
#SBATCH --cpus-per-task=8       ## number of cores the job needs
#SBATCH --error=slurm-%J.err    ## error log file
#SBATCH --output=slurm-%J.out   ##output info file

#-----------------------------------------------

# load modules
module load anaconda/2020.07

# activate conda environment
source activate breseq

# Directories
reference=data/ref
SI=data/processed/breseq/sister_I
SM=data/processed/breseq/sister_M
CI=data/processed/breseq/cousin_I
CM=data/processed/breseq/cousin_M

# Sister_I .tsv output
gdtools COMPARE -f TSV -o $SI/sister_I.tsv -r $reference/Bill.gbk -r $reference/Bob.gbk -r $reference/Car.gbk $SI/*/output/*.gd

# Sister_M .tsv output
gdtools COMPARE -f TSV -o $SM/sister_M.tsv -r $reference/Bill.gbk -r $reference/Bob.gbk -r $reference/Car.gbk $SM/*/output/*.gd

# Cousin_I .tsv output
gdtools COMPARE -f TSV -o $CI/cousin_I.tsv -r $reference/Bob.gbk -r $reference/CCS4.gbk -r $reference/SDS1.gbk $CI/*/output/*.gd

# Cousin_M .tsv output
gdtools COMPARE -f TSV -o $CM/cousin_M.tsv -r $reference/Bob.gbk -r $reference/CCS4.gbk -r $reference/SDS1.gbk $CM/*/output/*.gd
