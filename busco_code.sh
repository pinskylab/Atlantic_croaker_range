#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=busco
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=56000
#SBATCH --time=12:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /home/ksf63
source ~/.bashrc
conda activate busco
busco -i hifi4523.asm.bp.p_ctg.fa -l actinopterygii_odb10 -o busco_hifi_primary_output -m genome
