#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sila.horozoglu@partner.kit.edu
#SBATCH --job-name=ge_cluster


export PARNODES=$SLURM_NTASKS

eval "$(conda shell.bash hook)"
conda activate gpaw-env
python create_db_ger.py
python run_ger.py
