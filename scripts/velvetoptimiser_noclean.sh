#!/bin/bash

#SBATCH --time 48:00:00
#SBATCH --job-name=VelvetOptimiser
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal
#SBATCH --mem=180GB
#SBATCH --mail-type ALL
#SBATCH -A <user_account>
#SBATCH --mail-type ALL
#SBATCH --mail-user <myLinkBlueID@uky.edu>

echo "SLURM_NODELIST: "$SLURM_NODELIST

strainID=$1

lowK=$2

highK=$3

step=$4

mkdir $strainID 

  cp $strainID*1_paired*f*q* $strainID/

  cp $strainID*2_paired*f*q* $strainID/

  cd $strainID

# create hard-coded read names

  cp $strainID*1_paired*f*q* forward.fq

  cp $strainID*2_paired*f*q* reverse.fq

# run velvetoptimiser in singularity

  singularity run --app perlvelvetoptimiser226 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf VelvetOptimiser.pl \
  -s $lowK -e $highK -x $step -d velvet_${strainID}_${lowK}_${highK}_${step}_noclean -o ' -clean no' -f ' -shortPaired -fastq -separate forward.fq reverse.fq'

# remove read files with hard-coded names

  rm forward.fq

  rm reverse.fq 

  cd ..
