#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=hifiasm_assembly
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=56000
#SBATCH --time=12:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /scratch/ksf63/hifiasm
./hifiasm -o hifi4523.asm -t 32 m64190e_230310_212131.hifi_reads.fastq.gz
