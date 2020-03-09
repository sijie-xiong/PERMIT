#!/bin/bash

#SBATCH --job-name=vsci_csci_post
#SBATCH --output=vsci_csci_post.out
#SBATCH --error=vsci_csci_post.err
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
cp $HDIR/VSCI_CSCI_postproc.py  $SDIR
cd $SDIR

python VSCI_CSCI_postproc.py

rm -rf  $SDIR


