#!/bin/bash

#SBATCH --job-name=dps_vsvi_post
#SBATCH --output=dps_vsvi_post.out
#SBATCH --error=dps_vsvi_post.err
#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --time=24:00:00
#SBATCH --export=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=sijie.xiong@rutgers.edu	

HDIR=$SLURM_SUBMIT_DIR
SDIR="/scratch/$USER/$SLURM_JOB_ID"

mkdir -p $SDIR
cp $HDIR/helper.py $SDIR
cp $HDIR/NetComponents.py $SDIR
cp $HDIR/DPS_VSVI_postproc.py  $SDIR
cd $SDIR

python DPS_VSVI_postproc.py

rm -rf  $SDIR


