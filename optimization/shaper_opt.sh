#!/bin/bash

#SBATCH --job-name=shaper_opt
#SBATCH --output=shaper_opt.out
#SBATCH --error=shaper_opt.err
#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16000
#SBATCH --time=72:00:00
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --mail-user=sijie.xiong@rutgers.edu	

HDIR=$SLURM_SUBMIT_DIR
SDIR="/scratch/$USER/"
JDIR="$SDIR/$SLURM_JOB_ID"

mkdir -p $JDIR
cp $HDIR/NetSetting.m $JDIR
cp $HDIR/shaper_opt.m $JDIR
cp $HDIR/data/IoT/*.mat $JDIR
cd $JDIR

module load MATLAB/R2018b

export i_eps=$1 
export i_rho=$2
export shaper=$3

matlab -nodisplay -nosplash -r "shaper_opt, exit"

cp $JDIR/*.mat $HDIR/res
cd $SDIR
rm -rf  $SLURM_JOB_ID
