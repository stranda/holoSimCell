#!/bin/bash
#
# first command line param is the $TMPDIR assigned from slurm.  It comes in as $1 and is copied to TMPDIR
# 


TMPDIR=$1
cd $TMPDIR
mkdir csv
mkdir log
mkdir simdir

pwd

singularity exec /mnt/research/TIMBER/Ash/holosim.simg Rscript /mnt/research/TIMBER/Ash/code/hSC_Ash2.R 1 10 AES AshESA $TMPDIR/simdir $TMPDIR/csv

mv csv/*.csv /mnt/research/TIMBER/Ash/OUT
