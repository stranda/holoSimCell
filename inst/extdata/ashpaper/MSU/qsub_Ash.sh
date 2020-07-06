root=/mnt/research/TIMBER/Ash
nrep=20
#usr=f0007250
#usr=f0007264
#usr=f0007265
#usr=f0007266
#usr=f0011312

cd $root/SHELL

echo '#!/bin/sh
#SBATCH -C [intel16|intel18]
#SBATCH -N 1 -c 1
#SBATCH -t 4:00:00
#SBATCH --mem 8G
#SBATCH -J Ash_'$USER'

newgrp TIMBER 

cd '$root'/code/

singularity run /mnt/research/TIMBER/Ash/holosim.simg Rscript hSC_Ash2.R ${SLURM_ARRAY_TASK_ID} '$nrep' '$USER' AshESA '$root'/OUT/

scontrol show job ${SLURM_JOB_ID}' > hSC_Ash_$USER.sh

cd /mnt/research/TIMBER/Ash/SHELL
sbatch -a 1-10 -o "/mnt/research/TIMBER/Ash/QSTAT/hSC_%a_%A.out" hSC_Ash_$USER.sh

