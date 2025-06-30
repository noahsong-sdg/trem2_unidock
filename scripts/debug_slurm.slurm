#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --account=ucb-general
#SBATCH --job-name=debug_network
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=debug_network_%j.out
#SBATCH --error=debug_network_%j.err

echo "=== TESTING NETWORK FROM COMPUTE NODE ==="
echo "Hostname: $(hostname)"
echo "Date: $(date)"
echo ""

cd /scratch/alpine/noso3320/trem2/trem2_unidock/scripts
conda activate /projects/noso3320/hts_env

echo "Running network debug test..."
python debug_download.py

echo ""
echo "=== END DEBUG TEST ===" 
