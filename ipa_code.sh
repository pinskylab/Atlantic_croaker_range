#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=ipagenome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180000
#SBATCH --time=18:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /scratch/ksf63
source ~/.bashrc
conda activate ipa
ipa local --nthreads 7 --njobs 4 -i m64190e_230310_212131.hifi_reads.fastq.gz --resume
