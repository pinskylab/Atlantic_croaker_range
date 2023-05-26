#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=jellyfish
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=56000
#SBATCH --time=12:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /home/ksf63/Jellyfish
./jellyfish count -C -m 21 -s 56000M -t 32 /home/ksf63/m64190e_230310_212131.hifi_reads.fasta.gz
