#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=mapping
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100000
#SBATCH --time=72:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /projects/f_mlp195/kyraf/croaker/lcwgs/

source ~/.bashrc

## This script is used to quality filter and trim poly g tails. It can process both paired end and single end data.
SAMPLELIST=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/Sample_List_Combined.tsv} #Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample ta$
SAMPLETABLE=${2:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/Sample_Table_Combined.tsv} #Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sampl$
FASTQDIR=${3:-/projects/f_mlp195/kyraf/croaker/lcwgs/qual_filtered_q30/} #Path to the directory where fastq file are stored. An example for the quality-filtered Greenland pe data is: /workdir/cod/greenland-$
BASEDIR=${4:-/projects/f_mlp195/kyraf/croaker/lcwgs/} #Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be w$
FASTQSUFFIX1=${5:-_adapter_clipped_qual_filtered_f_paired.fastq.gz} #Suffix to fastq files. Use forward reads with paired-end data. An example for the quality-filtered Greenland paired-end data is: _adapter$
FASTQSUFFIX2=${6:-_adapter_clipped_qual_filtered_r_paired.fastq.gz} #Suffix to fastq files. Use reverse reads with paired-end data. An example for the quality-filtered Greenland paired-end data is: _adapter$
MAPPINGPRESET=${7:-very-sensitive} #The pre-set option to use for mapping in bowtie2 (very-sensitive for end-to-end (global) mapping [typically used when we have a full genome reference], very-sensitive-loc$
REFERENCE=${8:-/projects/f_mlp195/kyraf/croaker/lcwgs/hifi4523.asm.bp.p_ctg.filtered.fa} #Path to reference fasta file and file name, e.g /workdir/cod/reference_seqs/gadMor2.fasta
REFNAME=${9:-hifi} #Reference name to add to output files, e.g. gadMor2
THREADS=${10:-32} #Number of threads to use. Default is 8
MINQ=${11:-0} #Minimum mapping quality filter. Default is 0 (no filter)
SAMTOOLS=${12:-/home/ksf63/samtools-1.17}

## Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST`; do

        ## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
        SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
        SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
        LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
        SAMPLE_SEQ_ID=$SAMPLE_ID'_'$SEQ_ID'_'$LANE_ID

        ## Extract data type from the sample table
        DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`

        ## The input and output path and file prefix
        SAMPLETOMAP=$FASTQDIR$SAMPLE_SEQ_ID
        SAMPLEBAM=$BASEDIR'bam/'$SAMPLE_SEQ_ID

        ## Define platform unit (PU), which is the lane number
        PU=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`

        ## Define reference base name
        REFBASENAME="${REFERENCE%%.*}"
conda activate bowtie2

        ## Map reads to the reference
        # Map the paired-end reads
        if [ $DATATYPE = pe ]; then
        # We ignore the reads that get orphaned during adapter clipping because that is typically a very small proportion of reads. If a large proportion of reads get orphaned (loose their mate so they beco$
                bowtie2 -q --phred33 --$MAPPINGPRESET -p $THREADS -I 0 -X 1500 --fr --rg-id $SAMPLE_SEQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x /projects/f_mlp195/kyraf/croaker/lcwgs/hifi4523.asm.bp.p_ctg.filtered -1 $SAMPLETOMAP$FASTQSUFFIX1 -2 $SAMPLETOMAP$FASTQSUFFIX2 -S $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

        # Map the single-end reads
        elif [ $DATATYPE = se ]; then
                bowtie2 -q --phred33 --$MAPPINGPRESET -p $THREADS --rg-id $SAMPLE_SEQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -U $SAMPLETOMAP$FASTQSUFFIX1 -S  $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

        fi

                ## Convert to bam file for storage
        $SAMTOOLS view -b -@ $THREADS $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam' > $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam'
        rm $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

        ## Filter the mapped reads
        # Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20) -- do we want the quality score filter??
        $SAMTOOLS sort -@ $THREADS $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam' -o $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'_minq'$MINQ'_sorted.bam'

        $SAMTOOLS index -M -b -@ $THREADS $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'_minq'$MINQ'_sorted.bam'

done

