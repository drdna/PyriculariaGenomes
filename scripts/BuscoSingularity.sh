#!/bin/bash

#SBATCH --time 8:00:00
#SBATCH --job-name=busco
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal
#SBATCH --mem=180GB
#SBATCH --mail-type ALL
#SBATCH	-A cea_farman_s25abt480
#SBATCH --mail-type ALL
#SBATCH --mail-user farman@uky.edu

echo "SLURM_NODELIST: "$SLURM_NODELIST
echo "PWD :" $PWD

in=$1
out=${in/\.fasta/}_busco

singularity run --app busco570 /share/singularity/images/ccs/conda/amd-conda14-rocky8.sinf busco \
 --in $in --out $out --mode genome --lineage_dataset ascomycota_odb10 -f

