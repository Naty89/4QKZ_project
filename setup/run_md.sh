#!/bin/bash
#SBATCH --job-name=md_100ns
#SBATCH --partition=gpuA100x4
#SBATCH --account=bfam-delta-gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem-per-gpu=40G
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=md_%j.out
#SBATCH --error=md_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=gatnatiwos1@gmail.com

# Move to your working directory
cd /u/agebremedhin/rayca

# Load the Delta GROMACS module (GPU-enabled)
module purge
module load gromacs/2022.5.cuda

# Optional: print diagnostic info
echo "Running on host: $(hostname)"
echo "Using GPUs: $CUDA_VISIBLE_DEVICES"
echo "Starting time: $(date)"
echo ""

# Run MD (non-MPI, single GPU)
gmx mdrun -deffnm md -nb gpu -pme gpu -bonded gpu -ntomp 8

echo ""
echo "Job finished at: $(date)"

