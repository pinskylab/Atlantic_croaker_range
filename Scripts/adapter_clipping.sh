#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=adapter_clip
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100000
#SBATCH --time=24:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /projects/f_mlp195/kyraf/croaker/lcwgs
source ~/.bashrc
conda activate trimmomatic

 
SAMPLELIST=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/Sample_List_Run1.tsv}
SAMPLETABLE=${2:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/Sample_Table_Run1.tsv}
RAWFASTQDIR=${3:-/projects/f_mlp195/kyraf/croaker/lcwgs/data_run1/}
BASEDIR=${4:-/projects/f_mlp195/kyraf/croaker/lcwgs/}
RAWFASTQSUFFIX1=${5:-.fastq.gz}
RAWFASTQSUFFIX2=${6:-.fastq.gz}
ADAPTERS=$7 # Path to a list of adapter/index sequences. For Nextera libraries: /workdir/cod/reference_seqs/NexteraPE_NT.fa For BEST libraries: /workdir/cod/reference_seqs/BEST.fa
TRIMMOMATIC=${8:-/programs/trimmomatic/trimmomatic-0.39.jar} ## Path to trimmomatic (default /programs/trimmomatic/trimmomatic-0.39.jar)
THREADS=${9:-24}
JOBS=${10:-1}

JOB_INDEX=0
## Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST`; do
	
	## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
	SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
	SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
	LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
	SAMPLE_SEQ_ID=$SAMPLE_ID'_'$SEQ_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely
	
	## Extract data type from the sample table
	DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`
	
	## The input and output path and file prefix
	RAWFASTQ_ID=$RAWFASTQDIR$SAMPLEFILE
	SAMPLEADAPT=$BASEDIR'adapter_clipped/'$SAMPLE_SEQ_ID
	
	## Adapter clip the reads with Trimmomatic
	# The options for ILLUMINACLIP are: ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
	# For definitions of these options, see http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

	if [ $DATATYPE = pe ]; then
		 PE -threads $THREADS -phred33 $RAWFASTQ_ID$RAWFASTQSUFFIX1 $RAWFASTQ_ID$RAWFASTQSUFFIX2 $SAMPLEADAPT'_adapter_clipped_f_paired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_f_unpaired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_r_paired.fastq.gz' $SAMPLEADAPT'_adapter_clipped_r_unpaired.fastq.gz' 'ILLUMINACLIP:'$ADAPTERS':2:30:10:1:true' &
	
	elif [ $DATATYPE = se ]; then
		SE -threads $THREADS -phred33 $RAWFASTQ_ID$RAWFASTQSUFFIX1 $SAMPLEADAPT'_adapter_clipped_se.fastq.gz' 'ILLUMINACLIP:'$ADAPTERS':2:30:10' &
	fi
	
	JOB_INDEX=$(( JOB_INDEX + 1 ))
	if [ $JOB_INDEX == $JOBS ]; then
		wait
		JOB_INDEX=0
	fi
done
