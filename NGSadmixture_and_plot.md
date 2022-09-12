# Admixture 

The input is the global beagle file, whatever k you want and the minimum allele frequecy filture. \
The output will be a `.fopt.gz`, a `.qopt`, and a `.log` file. 

`ngsadmix.mpi`

			#!/bin/bash

			#SBATCH --job-name 12-admix
			#SBATCH --output %A_12-admix.out
			#SBATCH --error %A_12-admix.err
			#SBATCH --mail-type=ALL
			#SBATCH --mail-user=mabrober@ucsc.edu
			#SBATCH --time=15-00:00:00
			#SBATCH --partition=128x24
			#SBATCH --nodes=1
			#SBATCH --time=30-00:00:00
			#SBATCH --mem=120GB
			#SBATCH --ntasks-per-node=24

			#http://www.popgen.dk/software/index.php/NgsAdmix

			module load miniconda3.9
			conda activate angsd

			BASEDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/angsd/
			BEAGLE=$BASEDIR/DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.beagle.gz
			MINMAF=0.05
			K=12
			PREFIX=`echo $BEAGLE | sed 's/\..*//' | sed -e 's#.*/\(\)#\1#'`

			NGSadmix -likes $BEAGLE -K $K -minMaf $MINMAF -seed 4 -P 24 -o $BASEDIR'NGSadmix'/$PREFIX'struc_k'$K'_'

If you want to run a bunch of different ks this is a loop - it takes FOREVER. 
`ngsadmix-kloop.mpi`
	
			#!/bin/bash

			#SBATCH --job-name loop-admix
			#SBATCH --output %A_loop-admix.out
			#SBATCH --error %A_loop-admix.err
			#SBATCH --mail-type=ALL
			#SBATCH --mail-user=mabrober@ucsc.edu
			#SBATCH --time=15-00:00:00
			#SBATCH --partition=128x24
			#SBATCH --nodes=1
			#SBATCH --time=30-00:00:00
			#SBATCH --mem=120GB
			#SBATCH --ntasks-per-node=24

			#http://www.popgen.dk/software/index.php/NgsAdmix

			module load miniconda3.9
			conda activate angsd

			BASEDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/angsd/
			BEAGLE=DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.beagle.gz
			MINMAF=0.05
			MINK=5 # Minimum value of K
			MAXK=14 # Maximum value of K
			THREADS=24 # Number of threads to use, default value is 8, but the program can use a lot more if they are made available
			PREFIX=${BEAGLE%%.*}

			for ((K = $MINK; K <= $MAXK; K++)); do
				#run ngsAdmix
				echo $K
				NGSadmix \
				-likes $BASEDIR$BEAGLE \
				-K $K \
				-P $THREADS \
				-o $BASEDIR'NGSadmix/out_put_admix/loop_output/ngsadmix_'$PREFIX'_k'$K \
				-minMaf $MINMAF
			done

## Plot admixture barplot
`plot_admix.mpi`

		#!/bin/bash

		#SBATCH --job-name admix_plot
		#SBATCH --output %A_admix_plot.out
		#SBATCH --error %A_admix_plot.err
		#SBATCH --mail-type=ALL
		#SBATCH --mail-user=mabrober@ucsc.edu
		#SBATCH --time=15-00:00:00
		#SBATCH --ntasks=24
		#SBATCH --mem=120GB

		#run Rscript
		module load miniconda3.9
		conda activate tidyverse

		Rscript /hb/groups/bernardi_lab/may/DTR/population-analysis/angsd/NGSadmix/plot_admix.r

The r script: 		
		
`plot_admix.r` 

		#read in .qopt file - has individual ancestry proportion (admixture)
		q<-read.table("/hb/groups/bernardi_lab/may/DTR/population-analysis/angsd/NGSadmix/out_put_admix/loop_output/ngsadmix_DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_k5.qopt")
		#read in sample_table

		#pop<-read.table("/hb/groups/bernardi_lab/may/DTR/population-analysis/sample_lists/clean_sample_table_merged_groups_headers.txt", header = TRUE)
		pop<-read.table("/hb/groups/bernardi_lab/may/DTR/population-analysis/sample_lists/beagle_sample_order.tsv", header = FALSE)
		#pop<-read.table("/hb/groups/bernardi_lab/may/DTR/population-analysis/sample_lists/clean_sample_table_merged_headers.tsv", header = TRUE)
		k<-5
		png(paste0("/hb/groups/bernardi_lab/may/DTR/population-analysis/angsd/NGSadmix/out_put_admix/DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_k",k,"all_beagle_order.png"), width = 2500, height =1000, units ="px")
		barplot(t(q),
						col=1:k,
						names=pop$V1,
						las=2,
						space=0,
						border=NA,
		#       xlab="Individuals",
						ylab="Admixture proportions",
						main=paste("k=",k))
		dev.off()
