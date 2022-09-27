- [Estimate LD](#Estimate-LD)
- 	[Prep inout](#prep-input-files-and-info)
- [PCA](#pca-by-population)
	- [Load PCAngsd](#load-pcangsd)
	- [selection](#run-selection)
		- [plot pca](#plot-pca)
	
# Estimate LD 
### *only works on head node - will run tonight * 
From: https://github.com/nt246/lcwgs-guide-tutorial/blob/main/tutorial3_ld_popstructure/markdowns/ld.md
The estimation of linkage disequilibrium (LD) has important applications, e.g. for inference of population size, demographic history, selection, and for the discovery of structural variants. In addition, since many downstream analyses make assumptions about the independence of genomic loci, LD estimates are essential for trimming the list of loci to be included in these analyses (LD pruning).

### prep input files and info
The `--geno` file can be prepped by removing the the header and first 3 columns of the `.beagle.gz` file 
	
	zcat DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.beagle.gz | sed '1d' | gzip -c > DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.noheader.beagle.gz
	zcat DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.noheader.beagle.gz | cut -f 4- |  gzip -c > DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_ngsld_formatted.beagle.gz
	
The `--pos` file can be prepped by removing the 3rd column of the `.pos.gz` file

	zcat DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.pos.gz | cut -f 1,2 | gzip -c > DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_ngsld_formatted.pos.gz
	
To calculate n_sites I 

	zcat DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.beagle.gz |wc -l	
Which gave me 22573846  

`ngsLD.sh`

		#!/bin/bash

		#SBATCH --job-name ngsld
		#SBATCH --output %A_ngsld.out
		#SBATCH --error %A_ngsld.err
		#SBATCH --mail-type=ALL
		#SBATCH --mail-user=mabrober@ucsc.edu
		#SBATCH --time=15-00:00:00
		#SBATCH --partition=128x24
		#SBATCH --nodes=1
		#SBATCH --time=12-00:00:00
		#SBATCH --mem= 120GB
		#SBATCH --ntasks-per-node=24

		#Muh input
		BASEDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/angsd/
		BASENAME=DTR_dedup_bams_mindp100_maxdp590_minind13_minq20
		BEAGLE=$BASEDIR/DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_ngsld_formatted.beagle.gz
		POS=$BASEDIR/DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_ngsld_formatted.pos.gz
		N_IND=132
		N_SITES=22573846
		
		#Turn it on!
		module load ngsLD

		#Runit
		ngsLD \
		--geno $BEAGLE \
		--posH $POS \
		--probs \
		--n_ind $N_IND \
		--n_sites $N_SITES \
		--max_kb_dist 0 \
		--outH ${BASENAME}.ld

		# Options:
		# --probs: just says the input is genotype probabilities (likelihoods or posteriors)?
		# --n_ind INT: sample size (number of individuals).
		# --n_sites INT: total number of sites.
		# --max_kb_dist DOUBLE: maximum distance between SNPs (in Kb) to calculate LD. Set to 0(zero) to disable filter. [100]
		# --max_snp_dist INT: maximum distance between SNPs (in number of SNPs) to calculate LD. Set to 0 (zero) to disable filter. [0]
		# --n_threads INT: number of threads to use. [1]
		# --out FILE: output file name. [stdout]

# PCA by population
### Load pcangsd
Had a bit of a time getting this to run, turned out I just needed to re intall pcangsd. Since there's no simple conda install command available for this Here's how to do that (info from https://github.com/Rosemeis/pcangsd and https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-from-file)

	cd /hb/home/mabrober/programs
	mkdir pcangsd
	vim environment.yml
	
From the pcangsd page, open the `environment.yml` to see what dependencies are required. Basically you need this .yml file in the directory. It looks like this:

	name: pcangsd
	channels:
	    - defaults
	dependencies:
	    - python>=3.6
	    - numpy
	    - scipy
	    - cython
	    
Then we create the environment, activate it, and git clone and build:
	
	conda env create -f environment.yml
	conda activate pcangsd
	git clone https://github.com/Rosemeis/pcangsd.git
	cd pcangsd/
	python setup.py build_ext --inplace
	pip3 install -e .
	
Now we should have a working pcangsd program\

## Run selection
Runtime: 1h18m
`selection.mpi`

		#!/bin/bash

		#SBATCH --job-name selection
		#SBATCH --output %A_selection.out
		#SBATCH --error %A_selection.err
		#SBATCH --mail-type=ALL
		#SBATCH --mail-user=mabrober@ucsc.edu
		#SBATCH --time=15-00:00:00
		#SBATCH --partition=128x24
		#SBATCH --nodes=1
		#SBATCH --time=12-00:00:00
		#SBATCH --mem= 128GB
		#SBATCH --ntasks-per-node=24

		## This script is used to run PCA using pcangsd. It can be used to run individual-based PCA, estimate selection. The input is a beagle formatted genotype likelihood file.
		## https://github.com/Rosemeis/pcangsd

		module load miniconda3.9
		conda activate pcangsd

		BASEDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/angsd/
		BEAGLE=$BASEDIR/DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.beagle.gz
		MINMAF=0.05
		ANALYSIS=selection
		PCANGSD=pcangsd

		PREFIX=`echo $BEAGLE | sed 's/\..*//' | sed -e 's#.*/\(\)#\1#'`

		if [ $ANALYSIS = pca ]; then
			$PCANGSD --beagle $BEAGLE --minMaf $MINMAF --threads 16 -o $BASEDIR'pcangsd/pcangsd_'$PREFIX

		elif [ $ANALYSIS = selection ]; then
			$PCANGSD --beagle $BEAGLE --selection --minMaf $MINMAF --threads 16 -o $BASEDIR'pcangsd/pcangsd_'$PREFIX --sites_sav
		fi

Output files from this:

	pcangsd_DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.args
	pcangsd_DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.cov
	pcangsd_DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.selection.npy
	pcangsd_DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.sites
	
## Plot PCA from selection output
Input included `PCA_DAPC_source_functions.r`, the covariance matrix (.cov file), and your sample table - make sure you have headers - you'll need those to tell R what columns to plot.\
This is run on the head node, it takes 2 secs
	
	module load miniconda3.9
	conda activate tidyverse
	R
	
	library(tidyverse)
	library(ggplot2)
	library(cowplot)
	library(miscTools)

	#directions to PCA function
	source("/hb/groups/bernardi_lab/may/DTR/population-analysis/angsd/pcangsd/PCA_DAPC_source_functions.r")

	#define variables
	genome_cov <- read_table("/hb/groups/bernardi_lab/may/DTR/population-analysis/angsd/pcangsd/pcangsd_DTR_dedup_bams_mindp100_maxdp590_minind13_minq20_contrm.cov",col_names = FALSE) #covariance matrix
	sample_table<-read_tsv("/hb/groups/bernardi_lab/may/DTR/population-analysis/sample_lists/clean_sample_table_merged_headers.txt") #sample table with headers 

	#plot pcas
	PCA(genome_cov, sample_table$sampleID, as.character(sample_table$population), 1,2, show.ellipse = F, show.line = T, show.label = F) #
	ggsave(paste0('/hb/groups/bernardi_lab/may/DTR/population-analysis/angsd/pcangsd/pca/pca_unpruned_contrm__PC1_PC2_pops.png'),height=10, width=8)
	
	## This function (PCA) takes a covariance matrix and performs PCA.
	  # cov_matrix: a square covariance matrix generated by most pca softwares
	  # ind_label: a vector in the same order and length as cov_matrix; it contains the individual labels of the individuals represented in the covariance matrix
	  # pop_label: a vector in the same order and length as cov_matrix; it contains the population labels of the individuals represented in the covariance matrix
	  # x_axis: an integer that determines which principal component to plot on the x axis
	  # y_axis: an integer that determines which principal component to plot on the y axis
	  # show.point: whether to show individual points
	  # show.label: whether to show population labels
	  # show.ellipse: whether to show population-specific ellipses
	  # show.line: whether to show lines connecting population means with each individual point
	  # alpha: the transparency of ellipses
	  # index_exclude: the indices of individuals to exclude from the analysis
	
# current error "Warning messages:
1: In cov.trob(df) : Probable convergence failure
2: In cov.trob(df) : Probable convergence failure
3: In cov.trob(df) : Probable convergence failure
Warning messages:
1: In cov.trob(df) : Probable convergence failure
2: In cov.trob(df) : Probable convergence failure
3: In cov.trob(df) : Probable convergence failure"

	
