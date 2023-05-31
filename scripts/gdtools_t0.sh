#!/bin/bash

#--------------------------SBATCH settings------

#SBATCH --job-name=gdtools      ## job name
#SBATCH -A <account>            ## account name
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
sis_1=data/processed/breseq/sister_t0/Phage_Sisters_rep1/output
sis_2=data/processed/breseq/sister_t0/Phage_Sisters_rep2/output
sis_3=data/processed/breseq/sister_t0/Phage_Sisters_rep3/output
sis_out=data/processed/breseq/sister_t0
cus_1=data/processed/breseq/cousin_t0/Phage_Cousins_rep1/output
cus_2=data/processed/breseq/cousin_t0/Phage_Cousins_rep2/output
cus_3=data/processed/breseq/cousin_t0/Phage_Cousins_rep3/output
cus_out=data/processed/breseq/cousin_t0

# Sister_t0 .html output
gdtools COMPARE -o $sis_out/sister_t0.html -r $reference/Bill.gbk -r $reference/Bob.gbk -r $reference/Car.gbk $sis_1/Phage_Sisters_rep1.gd $sis_2/Phage_Sisters_rep2.gd $sis_3/Phage_Sisters_rep3.gd

# Sister_t0 .tsv output
gdtools COMPARE -f TSV -o $sis_out/sister_t0.tsv -r $reference/Bill.gbk -r $reference/Bob.gbk -r $reference/Car.gbk $sis_1/Phage_Sisters_rep1.gd $sis_2/Phage_Sisters_rep2.gd $sis_3/Phage_Sisters_rep3.gd

# Sister_t0 .gd output
gdtools COMPARE -f GD -o $sis_out/sister_t0.gd -r $reference/Bill.gbk -r $reference/Bob.gbk -r $reference/Car.gbk $sis_1/Phage_Sisters_rep1.gd $sis_2/Phage_Sisters_rep2.gd $sis_3/Phage_Sisters_rep3.gd

# Cousin_t0 .html output
gdtools COMPARE -o $cus_out/cousin_t0.html -r $reference/Bob.gbk -r $reference/CCS4.gbk -r $reference/SDS1.gbk $cus_1/Phage_Cousins_rep1.gd $cus_2/Phage_Cousins_rep2.gd $cus_3/Phage_Cousins_rep3.gd

# Cousin_t0 .tsv output
gdtools COMPARE -f TSV -o $cus_out/cousin_t0.tsv -r $reference/Bob.gbk -r $reference/CCS4.gbk -r $reference/SDS1.gbk $cus_1/Phage_Cousins_rep1.gd $cus_2/Phage_Cousins_rep2.gd $cus_3/Phage_Cousins_rep3.gd

# Cousin_t0 .gd output
gdtools COMPARE -f GD -o $cus_out/cousin_t0.gd -r $reference/Bob.gbk -r $reference/CCS4.gbk -r $reference/SDS1.gbk $cus_1/Phage_Cousins_rep1.gd $cus_2/Phage_Cousins_rep2.gd $cus_3/Phage_Cousins_rep3.gd
