#!/bin/bash
#SBATCH --partition=p_deenr_1
#SBATCH --job-name=call_snps
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100000
#SBATCH --time=72:00:00
#SBATCH --output=build_snpeff.%j.out
#SBATCH --error=build_snpeff.%N.%j.err

module use /projects/community/modulefiles
module load java/17.0.6

cd /home/ksf63/snpEff
java -Xmx8g -jar snpEff.jar ac45.23 /home/ksf63/snpEff/data/2320snps.vcf > 2320snps.ann.vcf
