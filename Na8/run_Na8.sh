#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sila.horozoglu@partner.kit.edu
#SBATCH --job-name=Na8_cluster
#SBATCH --output=job_Na8.out
#SBATCH --error=job_Na8.err

# === Activate Conda Environment ===
eval "$(conda shell.bash hook)"
conda activate gpaw-env

# === Load Turbomole (after activating conda) ===
module load chem/turbomole/7.9

# === Set SMP parallelization ===
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS

# Debugging: confirm Turbomole executables are found
echo "Using Turbomole define at: $(which define)"
echo "PARNODES = $PARNODES"

# === Run your ASE GA script ===
python create_db_na.py
python run_na.py
