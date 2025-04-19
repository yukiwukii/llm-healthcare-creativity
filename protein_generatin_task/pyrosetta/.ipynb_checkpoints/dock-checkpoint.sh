#!/bin/bash
#SBATCH --partition=UGGPU-TC1
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --time=60
#SBATCH --nodes=1
#SBATCH --job-name=binding_energy
#SBATCH --output=binding_energy_%j.out
#SBATCH --error=binding_energy_%j.err

# Load anaconda module
module load anaconda

# Initialize conda
eval "$(conda shell.bash hook)"
conda activate pyrosenv

echo "Calculating binding energy"
python workaround_docking.py

echo "Analysis completed"
