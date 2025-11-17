#!/usr/bin/env bash

#SBATCH -p nbi-medium,nbi-long,jic-medium,jic-long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 0-4:00                           # time (D-H:m)
#SBATCH -o ./slurm_output/multiqc.%J.out    # STDOUT
#SBATCH --mem=20G                           # memory pool for ALL cores

echo "--> $(date "+%F %T") : run_multiqc.slurm start"

# load MultiQC v1.7
source package 08be3662-986c-46d5-bb9a-3505b68cbad4

multiqc --version

cd qc

multiqc --force .

echo "--> $(date "+%F %T") : run_multiqc.slurm end"
