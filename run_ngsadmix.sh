#!/bin/bash
#SBATCH --partition=p_deenr_1
#SBATCH --job-name=ngsadmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100000
#SBATCH --time=72:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /projects/f_mlp195/kyraf/croaker/lcwgs/

## This script is used to run ngsAdmix. See http://www.popgen.dk/software/index.php/NgsAdmix for details.

INPUT_PATH=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/angsd/} # Path to the directory where the beagle formatted genotype likelihood file is stored. An example for the Greenland cod data is: /workdir/cod/g$
BEAGLE=${2:-bam_list_realigned_mindp132_maxdp4000_minind0.beagle.gz} # Name the beagle formatted genotype likelihood file. An example for the Greenland cod data is: bam_list_realigned_mincov_contamination_f$
MINMAF=${3:-0.05} # Minimum allele frequency filter
MINK=${4:-1} # Minimum value of K
MAXK=${5:-10} # Maximum value of K
THREADS=${6:-24} # Number of threads to use, default value is 8, but the program can use a lot more if they are made available
NGSADMIX=${7:-/home/ksf63/angsd/misc/NGSadmix} # Path to NGSAdmix, default value is '/programs/NGSadmix/NGSadmix'

PREFIX=${BEAGLE%%.*}

for ((K = $MINK; K <= $MAXK; K++)); do
        #run ngsAdmix
        echo $K
        $NGSADMIX \
        -likes $INPUT_PATH$BEAGLE \
        -K $K \
        -P $THREADS \
        -o $INPUT_PATH'/ngsadmix_'$PREFIX'_k'$K \
        -minMaf $MINMAF
done

