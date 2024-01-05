#!/bin/bash
#SBATCH --partition=p_deenr_1
#SBATCH --job-name=sfs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100000
#SBATCH --time=24:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /projects/f_mlp195/kyraf/croaker/lcwgs/bam_overlapclipped

## This script is used to get beagle formatted genotype likelihood file

DATA=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/bam_overlapclipped}
#BAMLIST=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/bam_list_realigned.txt} # Path to textfile listing bamfiles to include in global$
DIR=${2:-/projects/f_mlp195/kyraf/croaker/lcwgs/} # Path to the base directory where adapter clipped fastq file are stored in a subdirectory $
REF=${3:-/projects/f_mlp195/kyraf/croaker/lcwgs/hifi4523.asm.bp.p_ctg.filtered.fa} # Path to reference genome
#ANC=${4:-/projects/f_mlp195/kyraf/croaker/lcwgs/GCA_004119915.2_ASM411991v2_genomic.fa}

for POP in NJ DE MD VA NC OB LB SC GA FL
do
        echo $POP
        /home/ksf63/angsd/angsd -bam $DIR/sample_lists/$POP'_bams.txt' -ref $REF -anc $REF -out $POP -P 32 \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
                -minMapQ 20 -minQ 20 -minInd 0.5 -setMinDepth 132 -setMaxDepth 4000  -doCounts 1 \
                -GL 1 -doSaf 1
done

