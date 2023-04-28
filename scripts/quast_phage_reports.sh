#!/bin/bash

#--------------------------SBATCH settings------

#SBATCH --job-name=quast	## job name
#SBATCH -A katrine_lab		## account name
#SBATCH -p standard		## partition/queue name
#SBATCH --nodes=1		## number of nodes to use
#SBATCH --ntasks=1		## number of tasks to launch
#SBATCH --cpus-per-task=8	## number of cores the job needs
#SBATCH --error=slurm-%J.err	## error log file
#SBATCH --output=slurm-%J.out	##output info file

#-----------------------------------------------

# load modules
module load anaconda/2020.07

# define directories
assembly=data/processed/assembly

# quast assembly report
source activate quast
quast.py \
-o $assembly \
-t $SLURM_CPUS_PER_TASK \
--no-icarus \
$assembly/Phage*.fasta
