#!/bin/bash
#SBATCH --partition=UGGPU-TC1
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --time=120
#SBATCH --nodes=1
#SBATCH --job-name=large_dock
#SBATCH --output=large_dock_%j.out
#SBATCH --error=large_dock_%j.err

# Load anaconda module
module load anaconda

# Initialize conda
eval "$(conda shell.bash hook)"
conda activate pyrosenv

echo "Starting large-scale docking"

python workaround_docking.py --receptor ompT_cleaned.pdb --ligand_list ligands.txt --output_dir large_docking_results

echo "Docking completed"
