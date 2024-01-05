#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=multiqc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=56000
#SBATCH --time=24:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /projects/f_mlp195/kyraf/croaker/lcwgs/
source ~/.bashrc
conda activate multiqc

multiqc /projects/f_mlp195/kyraf/croaker/lcwgs/qual_filtered2_q20  --filename run2_q20_report

