#!/bin/bash
#SBATCH --time=60:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sila.horozoglu@partner.kit.edu
#SBATCH --job-name=na8_cluster

# === Activate Conda Environment ===
eval "$(conda shell.bash hook)"
conda activate gpaw-env

# === Load Turbomole (after activating conda) ===
module load chem/turbomole/7.9

# === Set SMP parallelization ===
export PARA_ARCH=SMP
export PARNODES=$SLURM_NTASKS


# === Run your ASE GA script ===
python create_db.py
python run_ga_na.py