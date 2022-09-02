# Fst estimates - angsd
http://www.popgen.dk/angsd/index.php/Fst \
For estimating per SNP Fst and its genome average\
*Had to install my own angsd conda here, for some reasone the hummingbird wide version did not have realSFS*

        #!/bin/bash

        #SBATCH --job-name angsd-fst
        #SBATCH --output %A_fst_angsd.out
        #SBATCH --error %A_fst_angsd.err
        #SBATCH --mail-type=ALL
        #SBATCH --mail-user=mabrober@ucsc.edu
        #SBATCH --time=15-00:00:00
        #SBATCH --partition=128x24
        #SBATCH --nodes=1
        #SBATCH --time=30-00:00:00
        #SBATCH --mem=120GB
        #SBATCH --ntasks-per-node=24

        ## This script is used to get pairwise Fst estimates from angsd for each population / group pair

        module load miniconda3.9
        conda activate angsd

        SAFDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/angsd/popminind6/ #  Path to per population saf.gz files.
        SAMPLETABLE=/hb/groups/bernardi_lab/may/DTR/population-analysis/sample_lists/clean_sample_table_merged.tsv # Path to the sample table
        POPCOLUMN=5 # The column index of the variable that you want to group by in the sample table above.
        BASENAME=_global_snp_list_DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_popminind6 # Base name of the saf files excluding ".saf.gz". It will be used as the base name of all output files.
        THREADS=24 # Number of parallel threads to use, default value is 8.
        EXTRA_ARG='' # Extra arguments for the SFS estimation step, default value is ''

        cd $SAFDIR

        I=1
        for POP1 in `tail -n +2 $SAMPLETABLE | cut -f $POPCOLUMN | sort | uniq`; do
            J=1
            for POP2 in `tail -n +2 $SAMPLETABLE | cut -f $POPCOLUMN | sort | uniq`; do
                if [ $I -lt $J ]; then
                    echo $POP1'_'$POP2
                    if [ ! -f $POP1$BASENAME'.saf.idx' ] || [ ! -f $POP2$BASENAME'.saf.idx' ]; then
                        echo 'One or both of the saf.idx files do not exist. Will proceed to the next population pair.'
                    else
                        # Check if Fst output already exists
                        if [ ! -f $POP1'_'$POP2$BASENAME'.fst' ]; then
                            # Generate the 2dSFS to be used as a prior for Fst estimation (and individual plots)
                            realSFS $POP1$BASENAME'.saf.idx' $POP2$BASENAME'.saf.idx' -P $THREADS $EXTRA_ARG > $POP1'_'$POP2$BASENAME'.2dSFS'
                            # Estimating Fst in angsd
                            realSFS fst index  $POP1$BASENAME'.saf.idx' $POP2$BASENAME'.saf.idx' -sfs $POP1'_'$POP2$BASENAME'.2dSFS' -fstout $POP1'_'$POP2$BASENAME'.alpha_beta'
                            realSFS fst print $POP1'_'$POP2$BASENAME'.alpha_beta.fst.idx' > $POP1'_'$POP2$BASENAME'.alpha_beta.txt'
                            awk '{ print $0 "\t" $3 / $4 }' $POP1'_'$POP2$BASENAME'.alpha_beta.txt' > $POP1'_'$POP2$BASENAME'.fst'
                        fi
                        # Check if average Fst output already exists
                        if [ ! -f $POP1'_'$POP2$BASENAME'.average_fst.txt' ]; then
                            # Estimating average Fst in angsd
                            realSFS fst stats $POP1'_'$POP2$BASENAME'.alpha_beta.fst.idx' > $POP1'_'$POP2$BASENAME'.average_fst.txt'
                        fi
                    fi
                fi
                J=$((J+1))
            done
            I=$((I+1))
        done

The output files for each pairwise comparison should be some thing like:\

        OMA_PHI_global_snp_list_DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_popminind6.2dSFS
        OMA_PHI_global_snp_list_DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_popminind6.alpha_beta.fst.gz
        OMA_PHI_global_snp_list_DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_popminind6.alpha_beta.fst.idx
        OMA_PHI_global_snp_list_DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_popminind6.alpha_beta.txt
        OMA_PHI_global_snp_list_DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_popminind6.average_fst.txt
        OMA_PHI_global_snp_list_DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_popminind6.fst
Some notes about some of these:\
- the headers for the `.fst` file should be:\
        `chromosome_name` `position` `alpha` `beta` `fst`
- the headers for the `average_fst.txt` should be * I THINK but and just deducing from clues here and there*\
        `unweighted_avg_fst` `weighted_avg_fst`
        - we want to use weighted_avg_fst(i. e. ratio of averages) (based on Bhatia 2013 https://genome.cshlp.org/content/23/9/1514.full.pdf)
        - 
