#!/bin/bash
#SBATCH --partition=p_deenr_1
#SBATCH --job-name=beagle
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100000
#SBATCH --time=72:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /projects/f_mlp195/kyraf/croaker/lcwgs/bam_overlapclipped

## This script is used to get beagle formatted genotype likelihood file

BAMLIST=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/bam_list_realigned.txt} # Path to textfile listing bamfiles to include in global SNP calling with absolute paths
BASEDIR=${2:-/projects/f_mlp195/kyraf/croaker/lcwgs/} # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be $
REFERENCE=${3:-/projects/f_mlp195/kyraf/croaker/lcwgs/hifi4523.asm.bp.p_ctg.filtered.fa} # Path to reference genome
SNPLIST=${4:-/projects/f_mlp195/kyraf/croaker/lcwgs/angsd/global_snp_list_bam_list_realigned_mindp132_maxdp4000_minind0.50_minq20.txt} # Path to the SNP list
CHRLIST=${5:-/projects/f_mlp195/kyraf/croaker/lcwgs/angsd/global_snp_list_bam_list_realigned_mindp132_maxdp4000_minind0.50_minq20.chrs} #Path to chromosome list

## Build base name of output files
OUTBASE=`echo $SNPLIST | sed 's/\..*//' | sed -e 's#.*snp_list_\(\)#\1#'`

## Get beagle formatted genotype likelihood
/home/ksf63/angsd/angsd \
-b $BAMLIST \
-anc $REFERENCE \
-out $BASEDIR'angsd/'$OUTBASE'_minmaf0.05' \
-GL 1 \
-doGlf 2 \
-doMaf 1 \
-minMaf 0.05 \
-doMajorMinor 3 \
-P 32 \
-sites $SNPLIST \
-rf $CHRLIST \
>& $BASEDIR'nohups/'$OUTBASE'_minmaf0.05_get_beagle.log'
