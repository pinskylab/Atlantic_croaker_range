#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=hifi_filter_busco
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=56000
#SBATCH --time=12:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /scratch/ksf63
source ~/.bashrc
conda activate busco
busco -i /scratch/ksf63/hifi4523.asm.bp.p_ctg.filtered.fa -l actinopterygii_odb10 -o hifi_filter_busco -m genome -c 32  
