#!/bin/bash

#--------------------------SBATCH settings--------------------------

#SBATCH --job-name=coverage             ## job name
#SBATCH -A katrine_lab                  ## account to charge
#SBATCH -p standard                     ## partition/queue name
#SBATCH --nodes=1                       ## number of nodes to use
#SBATCH --ntasks=1                      ## number of tasks to launch
#SBATCH --cpus-per-task=4               ## number of cores the job needs
#SBATCH --error=tmp/slurm_%A_%a.err     ## error log file name: %A is job id, %a is array task id
#SBATCH --output=tmp/slurm_%A_%a.out    ## output filename
#SBATCH --array=1-96                    ## number of array tasks

#-------------------------------------------------------------------

# load modules
module load bwa/0.7.8
module load samtools/1.15.1
module load anaconda/2020.07

# define directories
dirlist1=data/processed/dir_s.txt
dirlist2=data/processed/dir_c.txt
files1=data/processed/sis.txt
files2=data/processed/cus.txt
ref=/pub/adamsed/bwa.ref/Bob.fasta
output=data/processed/coverage

# define directories and prefix names
dir=`cat "$dirlist1" "$dirlist2" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
prefix=`cat "$files1" "$files2" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

# Run read mapping
bwa mem \
-t $SLURM_CPUS_PER_TASK \
-M $ref \
$dir/${prefix}_clean_1.fastq.gz \
$dir/${prefix}_clean_2.fastq.gz \
| samtools view -bS - \
| samtools sort - -o $output/${prefix}_bob_sort.bam

samtools index $output/${prefix}_bob_sort.bam

# activate conda environment
source activate deeptools

# extract region of interest with samtools
bamCoverage \
-p $SLURM_CPUS_PER_TASK \
-b $output/${prefix}_bob_sort.bam \
-o $output/${prefix}_bob.bw \
-bs 1
