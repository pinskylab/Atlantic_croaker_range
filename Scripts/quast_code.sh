#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=hifi_quast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=56000
#SBATCH --time=4:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /scratch/ksf63/
source ~/.bashrc 
conda activate quast
quast.py -t 32 -o hifi_filter_quast /scratch/ksf63/hifi4523.asm.bp.p_ctg.filtered.fa
