#!/usr/bin/env bash

# enter bash "strict" mode
set -euo pipefail

fq_files=$1     # list of sequencing files

# set job array from list of sequencing files
num_files=$(awk 'END {print NR}' "${fq_files}")
array_range="1-${num_files}%25"

mkdir -p slurm_output

unset SLURM_CPU_BIND    # technical thing for submitting to HPC from interactive node

a_job_id=$(sbatch --parsable --array=${array_range} master.slurm "${fq_files}")     # run alignment script master.slurm
m_job_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:${a_job_id}             run-multiqc.slurm)      # after all alignments complete, perform multiqc

echo "Alignment job id : $a_job_id"
echo "MultiQC job id   : $m_job_id"

squeue -u $USER