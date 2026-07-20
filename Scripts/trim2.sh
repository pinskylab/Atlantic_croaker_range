#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=quality_filter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100000
#SBATCH --time=24:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /home/ksf63/

## This script is used to quality filter and trim poly g tails. It can process both paired end and single end data.
SAMPLELIST=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/Sample_List_Combined.tsv}
SAMPLETABLE=${2:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/Sample_Table_Combined.tsv}
BASEDIR=${3:-/projects/f_mlp195/kyraf/croaker/lcwgs/}
RAWFASTQDIR=${4:-/projects/f_mlp195/kyraf/croaker/lcwgs/qual_filtered_q30/}
RAWFASTQSUFFIX1=${5:-_adapter_clipped_qual_filtered_f_paired.fastq.gz}
RAWFASTQSUFFIX2=${6:-_adapter_clipped_qual_filtered_r_paired.fastq.gz}
FILTER=${7:-quality}
THREADS=${8:-16}
FASTP=${9:/home/ksf63/fastp}

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
	RAWFASTQ_ID=$RAWFASTQDIR$SAMPLE_SEQ_ID
	SAMPLEQUAL=$BASEDIR'qual_filtered_q30_trim/'$SAMPLE_SEQ_ID
	
	## Trim polyg tail or low quality tail with fastp.
	# --trim_poly_g forces polyg trimming, --cut_right enables cut_right quality trimming
	# -Q disables quality filter, -L disables length filter, -A disables adapter trimming
	# Go to https://github.com/OpenGene/fastp for more information
	if [ $DATATYPE = pe ]; then
		if [ $FILTER = polyg ]; then
			$FASTP --trim_poly_g --cut_right -L -A --thread $THREADS -i $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' -I $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' -o $SAMPLEQUAL'_adapter_clipped_qual_filtered_f_paired.fastq.gz' -O $SAMPLEQUAL'_adapter_clipped_qual_filtered_r_paired.fastq.gz' -h $SAMPLEQUAL'_adapter_clipped_fastp.html' -j $SAMPLEQUAL'_adapter_clipped_fastp.json'
		elif [ $FILTER = quality ]; then
			$FASTP ./fastp --trim_front1 10 --trim_front2 10 --cut_front -L -A -Q --disable_trim_poly_g --thread $THREADS -i $RAWFASTQ_ID$RAWFASTQSUFFIX1 -I $RAWFASTQ_ID$RAWFASTQSUFFIX2 -o $SAMPLEQUAL'_adapter_clipped_qual_filtered_f_paired.fastq.gz' -O $SAMPLEQUAL'_adapter_clipped_qual_filtered_r_paired.fastq.gz' -h $SAMPLEQUAL'_adapter_clipped_fastp.html' -j $SAMPLEQUAL'_adapter_clipped_fastp.json'
		elif [ $FILTER = length ]; then
			$FASTP --max_len1 $MAXLENGTH -Q -L -A --thread $THREADS -i $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' -I $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' -o $SAMPLEQUAL'_adapter_clipped_qual_filtered_f_paired.fastq.gz' -O $SAMPLEQUAL'_adapter_clipped_qual_filtered_r_paired.fastq.gz' -h $SAMPLEQUAL'_adapter_clipped_fastp.html' -j $SAMPLEQUAL'_adapter_clipped_fastp.json'
		fi
	elif [ $DATATYPE = se ]; then
		if [ $FILTER = polyg ]; then
			$FASTP --trim_poly_g --cut_right -L -A --thread $THREADS -i $SAMPLEADAPT'_adapter_clipped_se.fastq.gz' -o $SAMPLEQUAL'_adapter_clipped_qual_filtered_se.fastq.gz' -h $SAMPLEQUAL'_adapter_clipped_fastp.html' -j $SAMPLEQUAL'_adapter_clipped_fastp.json'
		elif [ $FILTER = quality ]; then
			$FASTP -L -A --thread $THREADS -i $SAMPLEADAPT'_adapter_clipped_se.fastq.gz' -o $SAMPLEQUAL'_adapter_clipped_qual_filtered_se.fastq.gz' -h $SAMPLEQUAL'_adapter_clipped_fastp.html' -j $SAMPLEQUAL'_adapter_clipped_fastp.json'
		elif [ $FILTER = length ]; then
			$FASTP --max_len1 $MAXLENGTH -Q -L -A --thread $THREADS -i $SAMPLEADAPT'_adapter_clipped_se.fastq.gz' -o $SAMPLEQUAL'_adapter_clipped_qual_filtered_se.fastq.gz' -h $SAMPLEQUAL'_adapter_clipped_fastp.html' -j $SAMPLEQUAL'_adapter_clipped_fastp.json'
		fi
	fi
	
done
