#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=normal 
#SBATCH --partition=cpu
#SBATCH --account=ucb-general
#SBATCH --job-name=unidock_hts
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=unidock_hts_%j.out
#SBATCH --error=unidock_hts_%j.err

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=noso3320@colorado.edu

echo this is job $SLURM_JOBID
echo $SLURM_NTASKS

SCRIPT_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock/scripts"
DATA_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock/data"
RESULTS_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock/results"

module purge
module load intel
module load impi
module load cuda/11.8
module load rocm
module load gcc
module load cmake          
module load anaconda

echo "=== UniDock High-Throughput Screening Workflow ==="
cd /scratch/alpine/noso3320/trem2/trem2_unidock
conda activate /projects/noso3320/hts_env

python $SCRIPT_DIR/getdata.py
