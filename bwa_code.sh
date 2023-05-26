#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=hifi_bwa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=56000
#SBATCH --time=24:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /scratch/ksf63
module use /projects/community/modulefiles
module load Minimap2/minimap2-2.14
module load samtools/1.13-gc563
export SINGULARITY_BIND=/scratch/ksf63
minimap2 -ax map-pb -t 3 -2 hifi4523.asm.bp.p_ctg.fa m64190e_230310_212131.hifi_reads.fastq.gz > hifi.aln.sam
samtools view -Sb -o hifi.aln.bam hifi.aln.sam
samtools sort -o hifi.sort.bam hifi.aln.bam
samtools view -f 0x2 -b hifi.sort.bam hifi.sort.filt.bam
