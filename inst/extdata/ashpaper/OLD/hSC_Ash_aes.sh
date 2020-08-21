#!/bin/sh                                                                                                     
#SBATCH -C [intel16|intel18]                                                                                  
#SBATCH -N 1 -c 1                                                                                             
#SBATCH -t 4:00:00                                                                                            
#SBATCH --mem 8G                                                                                              
#SBATCH -J Ash_aes                                                                                          

newgrp TIMBER

/mnt/research/TIMBER/Ash/SHELL/runscriptMSU.sh ${TMPDIR}

scontrol show job ${SLURM_JOB_ID}


