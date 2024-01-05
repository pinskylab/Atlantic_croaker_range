#!/bin/bash
#SBATCH --partition=p_deenr_1
#SBATCH --job-name=pop_freq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100000
#SBATCH --time=14:30:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /projects/f_mlp195/kyraf/croaker/lcwgs/bam_overlapclipped

BASEDIR=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/} #  Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be$
SAMPLETABLE=${2:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/sample_table_merged_9pop.txt} # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the s$
POPCOLUMN=${3:-5} # The column index of the variable that you want to group by in the sample table above. In the Greenland project, it's the fifth column,  and thus 5
BAMLISTPREFIX=${4:bam_list_realigned_} # Prefix of the bam lists. An example from the Greenland cod project is bam_list_realigned_
REFERENCE=${5:-/projects/f_mlp195/kyraf/croaker/lcwgs/hifi4523.asm.bp.p_ctg.filtered.fa} # Path to reference genome
SNPLIST=${6:-/projects/f_mlp195/kyraf/croaker/lcwgs/angsd/global_snp_list_bam_list_realigned_mindp132_maxdp4000_minind0.50_minq20.txt} # Path to the SNP list
CHRLIST=${7:-/projects/f_mlp195/kyraf/croaker/lcwgs/angsd/global_snp_list_bam_list_realigned_mindp132_maxdp4000_minind0.50_minq20.chrs}
MINDP=${8:-132} # Minimum depth filter
MAXDP=${9:-4000} # Maximum depth filter
MININD=${10:-0.50} # Minimum individual filter
MINQ=${11:-20} # Minimum quality filter
MINMAPQ=${12:-20} # Minimum mapping quality (alignment score) filter, default value is 20
ANGSD=${13:-/home/ksf63/angsd/angsd} # Path to ANGSD, default value is /workdir/programs/angsd0.931/angsd/angsd
THREADS=${14:-10} # Number of parallel threads to use, default value is 8.
EXTRA_ARG=${15:-'-remove_bads 1 -only_proper_pairs 1 -C 50'} # Extra arguments when running ANGSD, default value is '-remove_bads 1 -only_proper_pairs 1 -C 50'
SAMPLEDIR=${16:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/}

OUTBASE=`echo $SNPLIST | sed 's/\..*//' | sed -e 's#.*\/\(\)#\1#'`
OUTDIR=$BASEDIR'angsd/popminind'$MININD'/'
if [ ! -d "$OUTDIR" ]; then
        mkdir $OUTDIR
fi

for POP in DE MD VA NC OB LB SC GA FL; do
        echo $POP
        $ANGSD \
        -b $SAMPLEDIR/$POP'_bams.txt' \
        -anc $REFERENCE \
        -ref $REFERENCE \
        -out $OUTDIR$POP'_'$OUTBASE'_popminind'$MININD \
        -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
        -P $THREADS \
        -setMinDepth $MINDP -setMaxDepth $MAXDP -minInd $MININD -minQ $MINQ -minMapQ $MINMAPQ \
        -sites $SNPLIST -rf $CHRLIST \
        $EXTRA_ARG \
        >& $BASEDIR'nohups/'$POP'_'$OUTBASE'_popminind'$MININD'_maf.log'
done
