#!/bin/bash

#--------------------------SBATCH settings------

#SBATCH --job-name=sym_lnk              ## job name
#SBATCH -A katrine_lab			## account name
#SBATCH -p standard                     ## partition/queue name
#SBATCH --cpus-per-task=1               ## number of cores the job needs
#SBATCH --error=tmp/slurm_%J.err        ## error log file
#SBATCH --output=tmp/slurm_%J.out       ## output info file

#-----------------------------------------------

DestDir=data/processed/t20
CI=data/processed/cousin_I
CM=data/processed/cousin_M/
SI=data/processed/sister_I/
SM=data/processed/sister_M/

for f in $(basename -a $CI/*_T20_*)
do
   echo "Processing $f file..."
   ln -sr $CI/$f $DestDir/$f
done

for f in $(basename -a $CM/*_T20_*)
do
   echo "Processing $f file..."
   ln -sr $CM/$f $DestDir/$f
done

for f in $(basename -a $SI/*_T20_*)
do
   echo "Processing $f file..."
   ln -sr $SI/$f $DestDir/$f
done

for f in $(basename -a $SM/*_T20_*)
do
   echo "Processing $f file..."
   ln -sr $SM/$f $DestDir/$f
done
