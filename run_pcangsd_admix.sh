#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=admix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100000
#SBATCH --time=24:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /projects/f_mlp195/kyraf/croaker/lcwgs/

source ~/.bashrc
conda activate pcangsd

## This script is used to run PCA using pcangsd. It can be used to run individual-based PCA, estimate selection, inbreeding coefficient, kinship, admixture, and others. The input is a beagle formatted genot$
## See https://github.com/Rosemeis/pcangsd for details

BASEDIR=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/} # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be $
BEAGLE=${2:-/projects/f_mlp195/kyraf/croaker/lcwgs/angsd/bam_list_realigned_mindp132_maxdp4000_minind0.beagle.gz} # Path to the beagle formatted genotype likelihood file
MINMAF=${3:-0.05} # Minimum allele frequency filter
ANALYSIS=${4:-admix} # Type of analysis with pcangsd. It can be one of the following: pca, selection, inbreedSites, kinship, admix
MINE=${5:-1} # Minimum number of eigenvectors to use in the modelling of individual allele frequencies (relevant for admix)
MAXE=${6:-5} # Maximum number of eigenvectors to use in the modelling of individual allele frequencies (relevant for admix)

PREFIX=`echo $BEAGLE | sed 's/\..*//' | sed -e 's#.*/\(\)#\1#'`

if [ $ANALYSIS = pca ]; then
        pcangsd -beagle $BEAGLE -minMaf $MINMAF -threads 16 -o $BASEDIR'angsd/pcangsd_'$PREFIX

elif [ $ANALYSIS = selection ]; then
        pcangsd -beagle $BEAGLE -selection -minMaf $MINMAF -threads 32 -o $BASEDIR'angsd/pcangsd_'$PREFIX -sites_save

elif [ $ANALYSIS = inbreedSites ]; then
        python2 /workdir/programs/pcangsd/pcangsd.py -beagle $BEAGLE -inbreedSites -minMaf $MINMAF -threads 16 -o $BASEDIR'angsd/pcangsd_'$PREFIX -sites_save

elif [ $ANALYSIS = kinship ]; then
        python2 /workdir/programs/pcangsd/pcangsd.py -beagle $BEAGLE -kinship -minMaf $MINMAF -threads 16 -o $BASEDIR'angsd/pcangsd_'$PREFIX

elif [ $ANALYSIS = admix ]; then
        #for ((E = $MINE; E <= $MAXE; E++)); do
                #echo $E
        pcangsd -beagle $BEAGLE -admix -minMaf $MINMAF -threads 32 -o $BASEDIR'angsd/pcangsd_'$PREFIX'_autoselect'
        #done


fi
