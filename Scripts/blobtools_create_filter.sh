#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=blobtools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=28000
#SBATCH --time=24:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /scratch/ksf63
source ~/.bashrc
conda activate blobtoolkit
blobtools create \
    --fasta /scratch/ksf63/hifi4523.asm.bp.p_ctg.filtered.fa \
    /scratch/ksf63/blobtools_outf

