#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=hifi_blast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100000
#SBATCH --time=48:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /scratch/ksf63
export SINGULARITY_BIND=/scratch/ksf63
export BLASTDB=/scratch/ksf63/nt
source ~/.bashrc
conda activate blast
blastn -db nt \
        -query hifi4523.asm.bp.p_ctg.fa \
        -outfmt "6 qseqid staxids bitscore std" \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads 40 \
        -out hififull.ncbi.blastn.out
