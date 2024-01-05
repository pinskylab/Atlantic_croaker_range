#!/bin/bash
#SBATCH --partition=p_deenr_1
#SBATCH --job-name=realign_indels
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=100000
#SBATCH --time=72:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /projects/f_mlp195/kyraf/croaker/lcwgs/bam_overlapclipped


## This script is used to quality filter and trim poly g tails. It can process both paired end and single end data.
BAMLIST=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/bam_list_dedup_overlapclipped.list} # Path to a list of merged, deduplicated, and overlap clipped bam files. Full paths should be included. A$
BASEDIR=${2:-/projects/f_mlp195/kyraf/croaker/lcwgs/} # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be $
REFERENCE=${3:-/projects/f_mlp195/kyraf/croaker/lcwgs/hifi4523.asm.bp.p_ctg.filtered.fa} # Path to reference fasta file and file name, e.g /workdir/cod/reference_seqs/gadMor2.fasta
SAMTOOLS=${4:-/home/ksf63/samtools-1.17/samtools} # Path to samtools
JOBS=${5:-64} # Number of indexing jobs to run in parallel (default 1)

JOB_INDEX=0

## Loop over each sample
#for SAMPLEBAM in `cat $BAMLIST`; do

#if [ -e $SAMPLEBAM'.bai' ]; then
        #echo "the file already exists"
#else
        ## Index bam files
        #$SAMTOOLS index $SAMPLEBAM &

        #JOB_INDEX=$(( JOB_INDEX + 1 ))
        #if [ $JOB_INDEX == $JOBS ]; then
                #wait
                #JOB_INDEX=0
        #fi
#fi

#done

#wait

## Realign around in-dels
# This is done across all samples at once

## Use an older version of Java
module use /projects/community/modulefiles
module load java/1.8.0_362

## Create list of potential in-dels
if [ ! -f $BASEDIR'bam_overlapclipped/all_samples_for_indel_realigner.intervals' ]; then
        java -XX:ParallelGCThreads=64 -Xmx100g -jar /home/ksf63/GenomeAnalysisTK.jar \
           -T RealignerTargetCreator \
           -R $REFERENCE \
           -I $BAMLIST \
           -o $BASEDIR'bam_overlapclipped/all_samples_for_indel_realigner.intervals' \
           -drf BadMate
fi

## Run the indel realigner tool
java -XX:ParallelGCThreads=64 -Xmx100g -jar /home/ksf63/GenomeAnalysisTK.jar \
   -T IndelRealigner \
   -R $REFERENCE \
   -I $BAMLIST \
   -targetIntervals $BASEDIR'bam_overlapclipped/all_samples_for_indel_realigner.intervals' \
   --consensusDeterminationModel USE_READS  \
   --nWayOut _realigned.bam

