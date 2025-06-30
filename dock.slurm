#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=normal 
#SBATCH --partition=aa100
#SBATCH --account=ucb-general
#SBATCH --job-name=unidock_hts
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --gres=gpu:a100:2
#SBATCH --output=unidock_hts_%j.out
#SBATCH --error=unidock_hts_%j.err

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=noso3320@colorado.edu

echo this is job $SLURM_JOBID
echo $SLURM_NTASKS

SCRIPT_DIR="/scratch/alpine/noso3320/trem2_unidock/scripts"
DATA_DIR="/scratch/alpine/noso3320/trem2_unidock/data"
RESULTS_DIR="/scratch/alpine/noso3320/trem2_unidock/results"

module purge
module load intel
module load impi
module load cuda/11.8
module load gcc
module load cmake          
module load anaconda

echo "=== UniDock High-Throughput Screening Workflow ==="
cd /scratch/alpine/noso3320/trem2_unidock

# Step 2: Check receptor file
RECEPTOR_FILE="$DATA_DIR/receptor/cluster1_receptor.pdbqt"
if [ ! -f "$RECEPTOR_FILE" ]; then
    echo "ERROR: Receptor file not found at $RECEPTOR_FILE"
    echo "Please ensure the TREM2 receptor is prepared and converted to PDBQT format."
    exit 1
else
    echo "Receptor file found: $RECEPTOR_FILE"
fi

# Step 3: Run UniDock molecular docking
echo "Step 6: Running UniDock molecular docking..."
echo "Using PDBQT files from: $DATA_DIR/ligands_pdbqt/"
python $SCRIPT_DIR/03_dock.py
