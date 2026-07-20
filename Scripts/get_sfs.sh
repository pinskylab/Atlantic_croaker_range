#!/bin/bash
#SBATCH --partition=p_deenr_1
#SBATCH --job-name=sfs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100000
#SBATCH --time=24:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /projects/f_mlp195/kyraf/croaker/lcwgs/angsd

for POP in NJ DE MD VA NC OB LB SC GA FL
do
        echo $POP
       /home/ksf63/angsd/misc/realSFS $POP.saf.idx -P 32 > $POP.sfs
done
