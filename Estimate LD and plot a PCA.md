- [Estimate LD](#Estimate-LD)
- [PCA](#pca-by-population)\
	- [Load PCAngsd](#load-pcangsd)\
	- [selection](#run-selection)\
	- [pca]
	
# Estimate LD 
### *STILL HAVE NOT BEEN ABLE TO GET THIS TO RUN* 
From: https://github.com/nt246/lcwgs-guide-tutorial/blob/main/tutorial3_ld_popstructure/markdowns/ld.md
The estimation of linkage disequilibrium (LD) has important applications, e.g. for inference of population size, demographic history, selection, and for the discovery of structural variants. In addition, since many downstream analyses make assumptions about the independence of genomic loci, LD estimates are essential for trimming the list of loci to be included in these analyses (LD pruning).

To calculate n_sites I 

	zcat DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.beagle.gz |wc -l
	
Which gave me 22573847 - because my beagle file has a header I subtract 1 from that. You could do this with your .pos file too since both files have a row for each snp (+ a header). 

`ngsLD.mpi`

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
		BEAGLE=$BASEDIR/DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.beagle.gz
		POS=$BASEDIR/DTR_dedup_bams_mindp100_maxdp590_minind13_minq20.pos.gz
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
	
## Plot PCA

