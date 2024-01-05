#!/bin/bash
#SBATCH --partition=p_deenr_1
#SBATCH --job-name=call_snps
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=100000
#SBATCH --time=72:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /projects/f_mlp195/kyraf/croaker/lcwgs/bam_overlapclipped
## This script is used to call SNPs using angsd

BAMLIST=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/bam_list_realigned.txt} # Path to textfile listing bamfiles to include in global SNP calling with absolute paths
BASEDIR=${2:-/projects/f_mlp195/kyraf/croaker/lcwgs/} # Path to the base directory where output files will be written to a subdirectory named "angsd/". An example for the Greenland cod data is: /workdir/cod$
REFERENCE=${3:-/projects/f_mlp195/kyraf/croaker/lcwgs/hifi4523.asm.bp.p_ctg.filtered.fa} # Path to reference genome
MINDP=${4:-132} # Minimum depth filter, e.g. 0.33 x number of individuals
MAXDP=${5:-4000} # Maximum depth filter, e.g. mean depth + 4 s.d.
MININD=${6:-0.50} # Minimum individual filter, minimum number of individuals a read has to be present in, e.g. 50% of individuals
MINQ=${7:-20} # Minimum quality filter
MINMAF=${8:-0.05} # Minimum minor allele frequency filter
MINMAPQ=${9:-20} # Minimum mapping quality (alignment score) filter, default value is 20
ANGSD=${10:-/home/ksf63/angsd/angsd} # Path to ANGSD, default value is /programs/angsd-0.940/angsd
THREADS=${11:-64} # Number of parallel threads to use, default value is 8.
EXTRA_ARG=${12:-'-remove_bads 1 -only_proper_pairs 1 -C 50'} # Extra arguments when running ANGSD, default value is '-remove_bads 1 -only_proper_pairs 1 -C 50'

## Extract the name of the bam list (excluding path and suffix)
BAMLISTNAME=`echo $BAMLIST | sed 's/\..*//' | sed -e 's#.*/\(\)#\1#'`

## Build base name of output files
OUTBASE=$BAMLISTNAME'_mindp'$MINDP'_maxdp'$MAXDP'_minind'$MININD'_minq'$MINQ

## Call SNPs
$ANGSD -b $BAMLIST -ref $REFERENCE -out $BASEDIR'angsd/'$OUTBASE \
-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 10000 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
-setMinDepth $MINDP -setMaxDepth $MAXDP -minInd $MININD \
-minQ $MINQ -minMapQ $MINMAPQ \
-SNP_pval 1e-6 -minMaf $MINMAF \
-P $THREADS \
$EXTRA_ARG \
>& $BASEDIR'nohups/'$OUTBASE'.log'

## Create a SNP list to use in downstream analyses
gunzip -c $BASEDIR'angsd/'$OUTBASE'.mafs.gz' | cut -f 1,2,3,4 | tail -n +2 > $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt'
$ANGSD sites index $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt'

## Also make it in regions format for downstream analyses
cut -f 1,2 $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt' | sed 's/\t/:/g' > $BASEDIR'angsd/global_snp_list_'$OUTBASE'.regions'

## Lastly, extract a list of chromosomes/LGs/scaffolds for downstream analysis
cut -f1 $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt' | sort | uniq > $BASEDIR'angsd/global_snp_list_'$OUTBASE'.chrs'


