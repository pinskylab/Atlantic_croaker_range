#!/bin/bash
#SBATCH --partition=p_deenr_1
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100000
#SBATCH --time=24:00:00
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%N.%j.err

## This script is used to get pairwise Fst estimates from angsd for each population / group pair

SAFDIR=${1:-/projects/f_mlp195/kyraf/croaker/lcwgs/angsd/} #  Path to per population saf.gz files. An example for the Greenland cod data is: /workdir/cod/greenland-cod/angsd/popminind2/
SAMPLETABLE=${2:-/projects/f_mlp195/kyraf/croaker/lcwgs/sample_lists/sample_table_merged_9pop.tsv} # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the s$
POPCOLUMN=${3:-5} # The column index of the variable that you want to group by in the sample table above. In the Greenland project, it's the fifth column, and thus 5
#BASENAME=$4 # Base name of the saf files excluding ".saf.gz". It will be used as the base name of all output files. An example from the Greenland cod project is _global_snp_list_bam_list_realigned_mindp161$
REALSFS=${4:-/home/ksf63/angsd/misc/realSFS} # Path to the realSFS module in ANGSD. Default value is /workdir/programs/angsd0.935/angsd/misc/realSFS
THREADS=${5:-32} # Number of parallel threads to use, default value is 8.
EXTRA_ARG=${6:-''} # Extra arguments for the SFS estimation step, default value is ''

cd $SAFDIR

I=1
for POP1 in `tail -n +2 $SAMPLETABLE | cut -f $POPCOLUMN | sort | uniq`; do
        J=1
        for POP2 in `tail -n +2 $SAMPLETABLE | cut -f $POPCOLUMN | sort | uniq`; do
                if [ $I -lt $J ]; then
                        echo $POP1'_'$POP2
                        if [ ! -f $POP1'.saf.idx' ] || [ ! -f $POP2'.saf.idx' ]; then
                                echo 'One or both of the saf.idx files do not exist. Will proceed to the next population pair.'
                        else
                                # Check if Fst output already exists
                                if [ ! -f $POP1'_'$POP2'.fst' ]; then
                                        # Generate the 2dSFS to be used as a prior for Fst estimation (and individual plots)
                                        $REALSFS $POP1'.saf.idx' $POP2'.saf.idx' -P $THREADS $EXTRA_ARG > $POP1'_'$POP2'.2dSFS'
                                        # Estimating Fst in angsd
                                        $REALSFS fst index  $POP1'.saf.idx' $POP2'.saf.idx' -sfs $POP1'_'$POP2'.2dSFS' -fstout $POP1'_'$POP2'.alpha_beta'
                                        $REALSFS fst print $POP1'_'$POP2'.alpha_beta.fst.idx' > $POP1'_'$POP2'.alpha_beta.txt'
                                        awk '{ print $0 "\t" $3 / $4 }' $POP1'_'$POP2'.alpha_beta.txt' > $POP1'_'$POP2'.fst'
                                fi
                                # Check if average Fst output already exists
                                if [ ! -f $POP1'_'$POP2'.average_fst.txt' ]; then
                                        # Estimating average Fst in angsd
                                        $REALSFS fst stats $POP1'_'$POP2'.alpha_beta.fst.idx' > $POP1'_'$POP2'.average_fst.txt'
                                fi
                        fi
                fi
                J=$((J+1))
        done
        I=$((I+1))
done

