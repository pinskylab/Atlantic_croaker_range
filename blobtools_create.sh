#!/bin/bash
#SBATCH --partition=p_mlp195
#SBATCH --job-name=blobtools_create
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=28000
#SBATCH --time=4:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err
cd /scratch/ksf63
source ~/.bashrc
conda activate blobtoolkit
blobtools create \
    --fasta ~/hifi4523.asm.bp.p_ctg.fa \
    /scratch/ksf63
blobtools add \
        --hits /scratch/ksf63/hifi.ncbi.blastn.out
        --taxdump /scratch/ksf63/taxdump \
        --taxrule bestsumorder \
        /scratch/ksf63
blobtools add \
        --busco /home/ksf63/busco_hifi_primary_output/run_actinopterygii_odb10/full_table.tsv \
        /scratch/ksf63
blobtools add \
        --cov /scratch/ksf63/hifi.sort.bam \
        /scratch/ksf63
