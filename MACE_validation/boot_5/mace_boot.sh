#!/bin/bash
#SBATCH --partition=cpu              # or cpu_il
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48           # one socket on EPYC nodes
#SBATCH --mem=32G
#SBATCH --time=15:00:00
#SBATCH --hint=nomultithread
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
#SBATCH --job-name=mace_mg11_cpu

set -euo pipefail
eval "$(conda shell.bash hook)"; conda activate gpaw-env
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
RUN_NAME="MACE_$(basename "$PWD")"
srun --cpu-bind=cores --hint=nomultithread \
  mace_run_train --config Mg11_config.yaml --name "$RUN_NAME"
