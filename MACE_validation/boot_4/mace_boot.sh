#!/bin/bash
#SBATCH --partition=cpu              # or cpu_il
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16           # one socket on EPYC nodes
#SBATCH --mem=32G
#SBATCH --time=15:00:00
#SBATCH --hint=nomultithread
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
#SBATCH --job-name=mace_mg11_cpu

set -euo pipefail
eval "$(conda shell.bash hook)"; conda activate gpaw-env
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

srun --cpu-bind=cores --hint=nomultithread \
  mace_run_train --config Mg11_config.yaml --name MACE_Mg11
