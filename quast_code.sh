#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=hifiasm_quast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=56000
#SBATCH --time=4:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /home/ksf63/
source ~/.bashrc
conda activate quast
quast.py -t 32 /home/ksf63/hifi4523.asm.bp.hap1.p_ctg.fa /home/ksf63/hifi4523.asm.bp.hap2.p_ctg.fa /home/ksf63/hifi4523.asm.bp.p_ctg.fa /home/ksf63/hifi4523.asm.bp.p_utg.fa /home/ksf63/hif$
