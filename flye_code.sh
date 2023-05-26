#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=flye
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180000
#SBATCH --time=32:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /scratch/ksf63
source ~/.bashrc
conda activate flye
flye --pacbio-hifi m64190e_230310_212131.hifi_reads.fastq.gz --out-dir flye_out --threads 32
