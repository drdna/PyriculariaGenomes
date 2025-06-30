#!/bin/bash

# Runs SNPcaller on files in specified directory

#SBATCH --time 1:00:00
#SBATCH --job-name=SNP_calling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal
#SBATCH --mem=32GB
#SBATCH --mail-type ALL
#SBATCH	-A cea_farman_s25abt480
#SBATCH --mail-type ALL
#SBATCH --mail-user farman@uky.edu,linkblue@uky.edu

echo "SLURM_NODELIST: "$SLURM_NODELIST
echo "PWD :" $PWD

blastfolder=$1
snpsfolder=${blastfolder/BLAST/SNPs}

export PERL5LIB=/project/farman_s25abt480/PERL_MODULES:$PERL5LIB
export PERL5LIB=/project/farman_uksr/PERL_MODULES:$PERL5LIB

perl /project/farman_s25abt480/SCRIPTs/Run_SU4.pl $blastfolder $snpsfolder


