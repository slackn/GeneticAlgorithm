#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sila.horozoglu@partner.kit.edu
#SBATCH --job-name=ge_cluster
#SBATCH --output=job.out
#SBATCH --error=job.err

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
python create_db_ger.py
python run_ger.py
