#!/bin/bash
#SBATCH --partition=p_deenr_1
#SBATCH --job-name=deduplicate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200000
#SBATCH --time=72:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /projects/f_mlp195/kyraf/croaker/lcwgs/bams_local_sort

module use /projects/community/modulefiles
module load java/17.0.6


## This script is used to deduplicate bam files and clipped overlapping read pairs for paired end data. It can process both paired end and single end data.
BAMLIST=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/bam_list_merged.tsv} # Path to a list of merged, deduplicated, and overlap clipped bam files. Full paths should be included. An example of su$
SAMPLETABLE=${2:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/sample_table_merged.tsv} # Path to a sample table where the 1st column is the prefix of the MERGED bam files. The 4th column is the sampl$
PICARD=${3:-/home/ksf63/picard.jar} # Path to picard tools
BAMUTIL=${4:-/home/ksf63/bamUtil} # Path to bamUtil


## Loop over each sample
for SAMPLEBAM in `cat $BAMLIST`; do

        ## Extract the file name prefix for this sample
        SAMPLESEQID=`echo $SAMPLEBAM | sed 's/_bt2_.*//' | sed -e 's#.*/bam/\(\)#\1#'`
        SAMPLEPREFIX=`echo ${SAMPLEBAM%.bam}`

        ## Remove duplicates and print dupstat file
        # We used to be able to just specify picard.jar on the CBSU server, but now we need to specify the path and version
        java -XX:ParallelGCThreads=64 -Xmx200g -jar $PICARD MarkDuplicates I=$SAMPLEBAM O=$SAMPLEPREFIX'_dedup.bam' M=$SAMPLEPREFIX'_dupstat.txt' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

        ## Extract data type from the merged sample table
        #DATATYPE=`grep -P "${SAMPLESEQID}\t" $SAMPLETABLE | cut -f 6`

        #if [ $DATATYPE != se ]; then
                ## Clip overlapping paired end reads (only necessary for paired end data)
                #$BAMUTIL ./bam clipOverlap --in $SAMPLEPREFIX'_dedup.bam' --out $SAMPLEPREFIX'_dedup_overlapclipped.bam' --stats
        #fi

done
