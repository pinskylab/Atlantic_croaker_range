#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=canu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180000
#SBATCH --time=24:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /scratch/ksf63
source ~/.bashrc
conda activate canu
canu \
-p canu -d canu \
genomeSize=800m \
useGrid=false \
-pacbio-hifi m64190e_230310_212131.hifi_reads.fastq.gz
