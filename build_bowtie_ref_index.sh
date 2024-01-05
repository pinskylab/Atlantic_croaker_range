#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=build_bowtie
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=56000
#SBATCH --time=24:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /projects/f_mlp195/kyraf/croaker/lcwgs
#source ~/.bashrc
#module use /projects/community/modulefiles
#module load java/17.0.6


## This script is used to build the bow tie reference index.
# Run this only when working with a new reference that has not been formatted for bowtie2

REFERENCE=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/hifi4523.asm.bp.p_ctg.filtered.fa}
REFNAME=${2:-hifi}

REFBASENAME="${REFERENCE%.*}"

## First create .bai and .dict files if they haven't been created
if [ ! -f $REFERENCE'.fai' ] ; then
        /home/ksf63/samtools-1.17/samtools faidx $REFERENCE
fi

#if [ ! -f $REFBASENAME'.dict' ] ; then
        #java -jar /home/ksf63/picard.jar CreateSequenceDictionary R=$REFERENCE O=$REFBASENAME'.dict'
#fi

#conda activate bowtie2

## Build the reference index
#bowtie2-build $REFERENCE $REFBASENAME
