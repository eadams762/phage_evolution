#!/bin/bash

#SBATCH --job-name=clean_reads		## job name
#SBATCH -A <account>			## account to charge
#SBATCH -p standard			## partition/queue name
#SBATCH --nodes=1               	## number of nodes to use
#SBATCH --ntasks=1              	## number of tasks to launch
#SBATCH --cpus-per-task=20		## number of cores the job needs
#SBATCH --mem-per-cpu=6G		## requested memory (6G = max)
#SBATCH --error=slurm_%A_%a.err		## error log file name: %A is job id, %a is array task id
#SBATCH --output=slurm_%A_%a.out	## output filename
#SBATCH --array=1-141			## number of array tasks

## load modules
module load bbmap/38.87
module load bowtie2/2.4.4

## define directories
input=data/raw
trim=data/raw/trim
dedupe=data/raw/dedupe
clean=data/processed/reads
stats=data/processed/stats
grch38=/dfs6/pub/adamsed/bowtie2.refs/GRCh38/grch38_1kgmaj_snvs

## create file name list
temp=$(basename -s _1.fastq.gz $input/*_1.fastq.gz | sort -u)

## select the file prefix corresponding to the array ID job number
prefix=`echo "$temp" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

## job commands

# trim adapters, phix, and low quality bases
bbduk.sh \
in=$input/${prefix}_1.fastq.gz \
in2=$input/${prefix}_2.fastq.gz \
ref=adapters,phix \
ktrim=r \
mink=11 \
hdist=1 \
qtrim=rl \
trimq=30 \
minlen=15 \
out=$trim/${prefix}_trim_1.fastq.gz \
out2=$trim/${prefix}_trim_2.fastq.gz \
stats=$stats/${prefix}_bbduk_stats.txt \
tpe \
tbo \
threads=$SLURM_CPUS_PER_TASK

#deduplicate reads
dedupe.sh \
in1=$trim/${prefix}_trim_1.fastq.gz \
in2=$trim/${prefix}_trim_2.fastq.gz \
out=$dedupe/${prefix}_dedupe.fastq.gz \
threads=$SLURM_CPUS_PER_TASK

# remove human reads
bowtie2 \
-x $grch38 \
-p $SLURM_CPUS_PER_TASK \
-S ${prefix}_human_reads.sam \
--interleaved $dedupe/${prefix}_dedupe.fastq.gz \
--un-conc-gz $clean/${prefix}_clean_%.fastq.gz

# remove intermediate files
rm ${prefix}_human_reads.sam
rm $trim/${prefix}_trim_*.fastq.gz
rm $dedupe/${prefix}_dedupe.fastq.gz
