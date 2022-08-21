# Stage 1: Data pre-processing
Low coverage whole genome sequencing analysis journey...all the nitty gritty 'cause we're learnin'\
All work run through slurm/university cluster\
This is version # 2.6k

- [Data Info](#data-info)
- [Data pre-process](#data-pre-process)
	- [Moving data](#Gather-sequences)
	- [Making support files of sample tables and sample lists for later use](#Create-sample-table-and-sample-lists)
	- [Fastp - QC, adapter trimming](#Fastp---QC,-adapter-trimming)
	- [Map, filter, reindex, reads](#Map-and-filter-reads)
	- [Merge replicate libraries if you got 'em](#Merge-replicate-libraries-of-samples)
	- [Deduplicate and clip overlapping portions of reads and reindex](#Deduplicate-and-clip-overlapping-portions-of-reads)
	- [Calculate coverage and determine appropriate snp calling parameters](#Calculate-coverage)

# Data Info
Low coverage data of 120+ samples of fish in 12 populations around the Pacific\
Libraries were prep'ed using Tn5 transposase (prob wouldn't rec at this point)\
Pools of Tn5 libraries named MR01, MR02, MR03, MR05, MR06, MR07, sent to Novogene in Sacramento, CA, USA

Data is stored in 3 directories: \
/hb/groups/bernardi_lab/may/data/pop_data/lcwgs/ \
WGS_09102021_NovogeneF001 libraries: (MR01,MR02) \
WGS_09102021_NovogeneF002 library: (MR03) \
WGS_05082022_NovogeneF001 libraries: (MR05, MR06, MR07)

# Data pre-process
## Get data in the right place and make tables and lists to be used throughout
### Gather sequences
To process, I copy all the `fq.gz` files in each of their individual directories to a directory called `all_raw_seq_data` using `gather_code.txt` which will look for and `find` all files within directories that end in `.fq.gz`and copy them to the destination file `GATHERED_DIR`.\
I run this from the directory that contains all the individual directories that hold the sequence files that I want (ie `WGS_05082022_NovogeneF001/raw_data/`. I did this directly in the command line (copy and paste lines below) on the head node.\

`gather_code.txt`

    DATA_DIR=/hb/groups/bernardi_lab/may/data/pop_data/lcWGS_data/pools01-03/WGS_09102021_NovogeneF001/raw_data
    GATHERED_DIR=/hb/groups/bernardi_lab/may/data/pop_data/lcWGS_data/pools01-03/all_raw_seq_data
   
    cd $DATA_DIR
    find . -type f -name *.fq.gz -exec cp {} $GATHERED_DIR \;

### Check MD5 to make sure data transfer was complete

    ls *.gz > list.txt
    list=./list.txt
    head list.txt
    for file in `cat $list`; do md5sum $file >> MD5check.txt; done
    
Here, I just copy and paste the lists of MD5s from the file of pre-transfer MD5s from the sequenceing facility and the list of MD5s produced here and paste into excel where I check for duplicate cells function (Home>Conditional Formatting>Highlight Cell Rules>Duplicate Values). Hopefully all the MD5s are duplicated showing that the file tranfers were successfully completed. 

## Create sample table and sample lists
Adapted from: Physalia Day 1 tutorial: https://github.com/nt246/physalia-lcwgs/blob/main/day_1/markdowns/data_processing.md#4-make-sure-youre-familiar-with-for-loops-and-how-to-assign-and-call-variables-in-bash

For the scripts below to work, the sample table has to be a tab deliminated table with the following six columns, strictly in this order:

1. prefix: the raw fastq file names **without** the `_1.fq.gz` ** ie. one file name per pair of reads

2. lane_number lane number; each sequencing lane or batch should be assigned a unique identifier. This is important so that if you sequence a library across multiple different sequencing lanes, you can keep track of which lane/batch a particular set of reads came from (important for accounting for sequencing error patterns or batch effects).

3. seq_id sequence ID; this variable is only relevant when different libraries were prepared out of the same sample and were run in the same lane (e.g. if you wanted to include a replicate). In this case, seq_id should be used to distinguish these separate libraries. If you only have a single library prepared from each of your samples (even if you sequence that same library across multiple lanes), you can just put 1 for all samples in this column.

4. sample_id sample ID; a unique identifier for each individual sequenced

5. population population name; the population or other relevant grouping variable that the individual belongs to

6. data_type data type; there are only two allowed entries here: pe (for paired-end data) or se (for single end data). We need this in the table because for some of our processing steps, the commands are slightly different for paired-end and single-end data.

7. *I added a column of 'uniq_ids' which is essentially: col4'_'col5'_'col3'_'col2'_' - it just simplified the scripts down the line for me.

It is important to make sure that the combination of sample_id, seq_id, and lane_number is unique for each fastq file. *ie. uniq_id

We need a second file that we call fastq list. This is simply a list of prefixes for the samples we want to analyze. Our sample table can contain data for all individuals in our study, but at any given time, we may only want to perform an operation on a subset of them. 

To make the table, first I wanted to see what all the headers were for each file. To make a list of file names and their headers:\
`get_headers.sh`

        cd /hb/groups/bernardi_lab/may/data/pop_data/lcWGS_data/new_raw_data
        DIR=/hb/groups/bernardi_lab/may/data/pop_data/lcWGS_data/sample_lists
        OUTFILE=file_headers2.txt
        CPDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/sample_lists

        for i in ./*.fq.gz;
        do
                echo $(echo $i; zcat $i| head -n 1) >> $DIR/$OUTFILE
        done

Then I pasted the info into excel to make the table and list because, call me a *dingdong*, but I like seeing WTF is going on. \
How I chose what info went into the columns for the table:
`sample_table.tsv`\
1. prefix: The whole, messy, file name without the fq.gz file ending
2. lane #: The number after L - caveat: for the most recent sequencing run, in order to make sure that those reads would be differentiated from the first batch, I changed lane 1 to 3 and 2 to 4. 
3. seq_ID: the 7th or 8th letter of the file name. If it is a number it is a replicate sequence (differently indexed individual), if it is a "C" there is no replicate. - Need to check whether I am resequencing some individuals. - I am IRI16!
4. sample_ID: First 5 letters of the file name \ ex:AMA04 \
5. population_ID: First 3 letters of the file name \ ex: AMA\
6. data type: they are all paired end reads so "pe" for all \
`sample_list`

## Fastp - QC, adapter trimming
Runtime : 9+ hours \
Can't remember if I tried making this an array - would be good here though if possible.\

`1-fastp.mpi`

        #!/bin/bash
        #SBATCH --job-name fp_all
        #SBATCH --output=%j-fastp_all.out
        #SBATCH --error=%j-fastp_all.err
        #SBATCH --mail-type=ALL
        #SBATCH --mail-user=mabrober@ucsc.edu
        #SBATCH --partition=128x24
        #SBATCH --nodes=1
        #SBATCH --time=30-00:00:00
        #SBATCH --mem=120GB
        #SBATCH --ntasks-per-node=16

        module load miniconda3.9
        conda activate /hb/scratch/mabrober/programs/fastp

        #Hybrid code with physalia (which uses fastqc - I use fastp) to mostly to have naming of downstream files be consistent.

        RAWDATA=/hb/groups/bernardi_lab/may/data/pop_data/lcWGS_data/all_raw_data #path to raw fq.gz files
        OUTDIR=/hb/groups/bernardi_lab/may/data/pop_data/lcWGS_data/
        WORKDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/
        SAMPLELIST=$WORKDIR/sample_lists/sample_list.txt # Path to a list of prefixes of the raw fastq files. It can be a subset of the the 1st column of the sample table.
        SAMPLETABLE=$WORKDIR/sample_lists/sample_table_all.tsv # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se.
        RAW_R1=_1.fq.gz # Suffix to raw fastq files. Use forward reads with paired-end data.
        RAW_R2=_2.fq.gz # Suffix to raw fastq files. Use reverse reads with paired-end data.
        HTML=$OUTDIR/reports_fastp/html_reports #where output reports will live
        JSON=$OUTDIR/reports_fastp/json_reports
        TXT=$OUTDIR/reports_fastp/txt_reports

        ## Loop over each sample

        for SAMPLEFILE in `cat $SAMPLELIST`; do

            ## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library. This is for the naming of trimmed/processed files
            SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
            POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
            SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
            LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
            SAMPLE_UNIQ_ID=$SAMPLE_ID'_'$POP_ID'_'$SEQ_ID'_'$LANE_ID  # When a sample has been sequenced in multiple lanes, we need to be able to identify the files from each run uniquely

            ## The input and output path and file prefix
                SAMPLEADAPT=$OUTDIR'/adapter_clipped/'$SAMPLE_UNIQ_ID

               fastp -i $RAWDATA/$SAMPLEFILE$RAW_R1\
                -I $RAWDATA/$SAMPLEFILE$RAW_R2\
                -o ${SAMPLEADAPT}_adapter_clipped_f_paired.fastq.gz\
                -O ${SAMPLEADAPT}_adapter_clipped_r_paired.fastq.gz\
                -l 30\
                -h $REPORT_DIR/$HTML/${SAMPLE_UNIQ_ID}_fastp.html
                -j $REPORT_DIR/$HTML/${SAMPLE_UNIQ_ID}_fastp.json
                -R $OUT_DIR/${FILES}fastp_report\

        done
        
Once `fastp` is done, go to the `reports_html` directory where all the individual .html files were output to (/hb/groups/bernardi_lab/may/data/pop_data/lcWGS_data/reports_fastp/html_reports). Then,

    module load miniconda3.9
    conda activate multiqc
    multiqc .
    
This will pull together and summarize all the reports into one pretty html report. <3

## Map and filter reads
Bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer\
### Index reference bowtie2-build
Before we map the reads to the reference genome, we index them
runtime:17min\
`1-bowtie_index.mpi`

        #!/bin/bash
        #
        #SBATCH --job-name bt-index
        #SBATCH --output=%j-bt-index.out
        #SBATCH --error=%j-bt-index.err
        #SBATCH --mail-type=ALL
        #SBATCH --mail-user=mabrober@ucsc.edu
        #SBATCH --partition=128x24
        #SBATCH --nodes=1
        #SBATCH --time=30-00:00:00
        #SBATCH --mem=120GB
        #SBATCH --ntasks-per-node=12

        module load miniconda3.9
        conda activate bowtie2
        
        DIR=/hb/groups/bernardi_lab/may/DTR/DTRgenome/assembly/assembly-versions/masurca-wkfl/5-final-hic-mappped/
        REF=$DIR/kuro_gaps_homogenized_output.fasta

        cd $DIR

        bowtie2-build $REF ./

### Map and filter - bowtie2
Using an array since I have 190 pe seqs to process.\
Run tim: 7 hours 
`2-mapnfilt_array.mpi`

        #!/bin/bash
        #SBATCH --job-name mapnfilt
        #SBATCH --output=%A_%a_mapnfilt.out
        #SBATCH --error=%A_%a_mapnfilt.err
        #SBATCH --mail-type=ALL
        #SBATCH --mail-user=mabrober@ucsc.edu
        #SBATCH --partition=128x24
        #SBATCH --time=30-00:00:00
        #SBATCH --mem=5GB
        #SBATCH --array=1-139%24 #need this to tell slurm I want an array with 139 iterations(samples) with up to 24 running at a time (calc available mem and threads)
        #SBATCH --ntasks=1

        echo job: mapping (bowtie) and filtering (samtools)
        module load miniconda3.9
        conda activate bowtie2
        module load samtools

    ## Where is everything, what are they called
        BASEDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/
        OUTDIR=$BASEDIR/1-mapnfilter/bams/
        SAMPLELIST=$BASEDIR/sample_lists/sample_list.txt # Path to a list of prefixes of the raw fastq files (looks like "PHI06_3_CKDL210014824-1a-AK13786-AK40299_HJM3WDSX2_L2"). Subset of 1st column ofsample table.
        SAMPLETABLE=$BASEDIR/sample_lists/sample_table_all_uid.tsv # Path to a sample table of info for each sample seqs
        FASTQDIR=/hb/groups/bernardi_lab/may/data/pop_data/lcWGS_data/adapter_clipped/ # Path to the directory of input fastq files
        FASTQSUFFIX1=_adapter_clipped_f_paired.fastq.gz # Suffix to forward fastq files
        FASTQSUFFIX2=_adapter_clipped_r_paired.fastq.gz # Suffix to reverse fastq files
        REF=/hb/groups/bernardi_lab/may/DTR/DTRgenome/assembly/assembly-versions/masurca-wkfl/5-final-hic-mappped/kuro_gaps_homogenized_output.fasta # Path to reference fasta file and file name
        REFNAME=kuro # Reference genome name to add to output files

    ## Set parameter
        MAPPINGPRESET=very-sensitive # The pre-set option to use for mapping in bowtie2 (very-sensitive for end-to-end (global) mapping [typically used when we have a full genome reference], very-sensitive-local for partial read mapping that allows soft-clipping [typically used when mapping genomic reads to a transcriptome]

    ## Make sure the array creates jobs for each line in sample_list.txt
        SAMPLEFILE=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $SAMPLELIST)
        echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
        echo Samplefile is $SAMPLEFILE

    ## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
        SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
        POP_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 5`
        SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
        SAMPLE_UNIQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 7 `
        PU=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
        DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`

    ## The input and output path and file prefix
        SAMPLETOMAP=$FASTQDIR$SAMPLE_UNIQ_ID
        SAMPLEBAM=$OUTDIR$SAMPLE_UNIQ_ID

    ## Define reference base name
        REFBASENAME="${REF%.*}"

    # Map reads to the reference
        echo sample uniq ID is $SAMPLE_UNIQ_ID

        bowtie2 -q --phred33 --$MAPPINGPRESET -p 1 -I 0 -X 1500 --fr --rg-id $SAMPLE_UNIQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -1 $SAMPLETOMAP$FASTQSUFFIX1 -2 $SAMPLETOMAP$FASTQSUFFIX2 -S $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

        ## Convert to bam file for storage (including all the mapped reads)
        samtools view -bS -F 4 $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam' > $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam'
        rm -f $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'

        ## Filter the mapped reads (to only retain reads with high mapping quality)
        # Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20) -- do we want the quality score filter??
        samtools view -h -q 20 $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam' | samtools view -buS - | samtools sort -o $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'_minq20_sorted.bam'
    echo Weeeeeeeeee!
    
### 
'3-checkoutput.mpi'\


## Merge replicate libraries of samples
Run time : ~2 hours to merge over ~110 files into ~30\
`4-merge_reps.mpi`

        #!/bin/bash
        #
        #SBATCH --job-name merge-reps
        #SBATCH --output=%j-merge_reps.out
        #SBATCH --error=%j-merge_reps.err
        #SBATCH --mail-type=ALL
        #SBATCH --mail-user=mabrober@ucsc.edu
        #SBATCH --partition=128x24
        #SBATCH --nodes=1
        #SBATCH --time=30-00:00:00
        #SBATCH --mem=120GB
        #SBATCH --ntasks-per-node=6

        module load samtools

        cd /hb/groups/bernardi_lab/may/DTR/population-analysis/1-mapnfilter/bams/
        OUTDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/1-mapnfilter/merge_test/
        #file=("BAK12_BAK" "BAK13_BAK" "FAL02_FAL" "FAL04_FAL" "FAL05_FAL" "FAL06_FAL" "FAL07_FAL""IRI07_IRI" "IRI08_IRI" "IRI09_IRI" "IRI13_IRI" "IRI16_IRI" "KIN05_KIN" "KIN08_KIN" "KIN06_KIN" "MOO01_MOO" "MOO04_MOO" "MOO05_MOO" "MOO06_MOO" "MOO07_MOO" "OOUR01_OKI" "OOUR01_OKI" "OZPA03_OKI" "OZPA04_OKI" "PHI05_PHI" "PHI06_PHI" "PHI07_PHI") #list file name groupings here #must do "OZPA07_OKI" still
        #file=("FAL07_FAL") #fix
        file=("MOO07_MOO")  #another fix
	
        for i in ${file[@]};   #for all items in $file
        do
        input_file_1="${i}_1_*_pe_bt2_kuro_minq20_sorted.bam"
        input_file_2="${i}_2_*_pe_bt2_kuro_minq20_sorted.bam"
        #input_file_3="${i}_3_*_pe_bt2_kuro_minq20_sorted.bam"
        out_file="${i}_M_2_pe_bt2_kuro_minq20_sorted.bam" #M for merged # luckily, all file to be merged were in lane 2

        samtools merge $OUTDIR/$out_file $input_file_1 $input_file_2 #$input_file_3
        done

## Deduplicate and clip overlapping portions of reads
Runtime: damn, deleted the email with the info I guess\
'5-dedup_clip.mpi'

        #!/bin/bash
        #
        #SBATCH --job-name dedup_clip
        #SBATCH--output=%A_%a_dedup_clip.out
        #SBATCH --error=%A_%a_dedup_clip.err
        #SBATCH --mail-type=ALL
        #SBATCH --mail-user=mabrober@ucsc.edu
        #SBATCH --partition=128x24
        ##SBATCH --nodes=1
        #SBATCH --time=30-00:00:00
        #SBATCH --array=1-2%
        #SBATCH --mem=90GB
        #SBATCH --ntasks-per-node=24
        echo `date` job: removing duplicate reads (picard) and clipping overlapping ends of the f/r paired reads (bamutil)

        module load miniconda3.9
        conda activate /hb/scratch/mabrober/programs/arima_programs/picard

        BAMDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/1-mapnfilter/bams
        BAMLIST=/hb/groups/bernardi_lab/may/DTR/population-analysis/sample_lists/dedup_fixlist2.txt
        REF=/hb/groups/bernardi_lab/may/DTR/DTRgenome/assembly/assembly-versions/masurca-wkfl/5-final-hic-mappped//kuro_gaps_homogenized_output.fasta
        REFNAME=kuro # Reference name to add to output files
        PICARD=/hb/scratch/mabrober/programs/arima_programs/picard/share/picard-2.26.2-0/picard.jar
        BAMUTIL=/hb/home/mabrober/.conda/envs/bamutil/bin/bam #bam -h for usage

        SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
        echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
        echo Sample bam is $SAMPLEBAM'pe_bt2_kuro_minq20_sorted.bam'

        ## Loop over each sample
        #for SAMPLEBAM in `cat $BAMLIST`; do

        ## Remove duplicates and print dupstat file
            java -Xmx60g -jar $PICARD MarkDuplicates I=$BAMDIR/$SAMPLEBAM'pe_bt2_kuro_minq20_sorted.bam' O=$BAMDIR/$SAMPLEBAM'bt2_kuro_minq20_sorted_dedup.bam' M=$BAMDIR/$SAMPLEBAM'bt2_kuro_minq20_sorted_dupstat.txt' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

        ## Clip overlapping paired end reads (only necessary for paired-end data)
            conda activate bamutil
            bam clipOverlap --in $BAMDIR/$SAMPLEBAM'bt2_kuro_minq20_sorted_dedup.bam' --out $BAMDIR/deduped_overlap_clipped_bams/$SAMPLEBAM'bt2_kuro_minq20_sorted_dedup_overlapclipped.bam' --stats

        echo done-zo woot! 
        
### reindex the deduplicated and overlap - clipped reads 
 - while also making a new list for next step\
Runtime: 20 min?
 '6-reindex.mpi"
 
        #!/bin/bash
        #
        #SBATCH --job-name index
        #SBATCH--output=%A_%a-index.out
        #SBATCH --error=%A_%a-index.err
        #SBATCH --mail-type=ALL
        #SBATCH --mail-user=mabrober@ucsc.edu
        #SBATCH --partition=128x24
        #SBATCH --nodes=1
        #SBATCH --time=30-00:00:00
        #SBATCH --mem=20GB
        #SBATCH --array=1-138%18
        #SBATCH --ntasks-per-node=4

        echo `date` job: re-index bam files

        BAMDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/1-mapnfilter/bams/deduped_overlap_clipped_bams/
        BAMLIST=/hb/groups/bernardi_lab/may/DTR/population-analysis/sample_lists/index_bam_list_dedup_overlapclipped_fix.list
        REF=/hb/groups/bernardi_lab/may/DTR/DTRgenome/assembly/assembly-versions/masurca-wkfl/5-final-hic-mappped/kuro_gaps_homogenized_output.fasta
        REFNAME=kuro # Reference name to add to output files

        #Set up array - assigns sample to individula job within array from BAMLIST
        SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
        echo Slurm array task ID is: $SLURM_ARRAY_TASK_ID
        echo Sample bam is $SAMPLEBAM'bt2_kuro_minq20_sorted_dedup_overlapclipped.bam'

        module load samtools

        cd $BAMDIR/
        samtools index $SAMPLEBAM'bt2_kuro_minq20_sorted_dedup_overlapclipped.bam'
      
### cursory check bams
Here, I don't know if this is legit but I just `ls -lah` the dir with all the $SAMPLEBAMbt2_kuro_minq20_sorted_dedup_overlapclipped.bams and `mv`ed any bams that were <500M in size into a`bad_bams` dir. \_"]_/

### make list 
Now we need lists of bams preceded by paths to the bam files.  
 
This one below makes a list of ALL the bams. To make population level lists, use `grep "AMA" > AMA.list` for example, to make the population lists. 
 
`make.list.sh`
        LISTDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/sample_lists
        BAMDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/bams/deduped_overlap_clipped_bams
	LIST=bam.list
	
        for SAMPLEBAM in `cat $LISTDIR/$LIST`; do
                echo $BAMDIR/$SAMPLEBAM'bt2_kuro_minq20_sorted_dedup_overlapclipped.bam' >> $LISTDIR/bam_path_list.txt
        done
       
## Calculate final coverage
This uses `samtools depth` and runs an array to count depth per position for each .bam file. The out put is a `bam.depth.gz` in the same place where your .bam files live. \
Run time: 20 min \
`1-count_depth_per_position_per_sample.mpi`

	#!/bin/bash
	#
	#SBATCH --job-name 1readdepth
	#SBATCH --output=%A_%a-1readdepth.out
	#SBATCH --error=%A_%a-1readdepth.err
	#SBATCH --mail-type=ALL
	#SBATCH --mail-user=mabrober@ucsc.edu
	#SBATCH --partition=128x24
	##SBATCH --nodes=1
	#SBATCH --time=7-00:00:00
	#SBATCH --mem=30GB
	##SBATCH --ntasks-per-node=6
	#SBATCH --array=1-132%12
	## This script is used to count per position depth for bam files. It will create one depth file per bam file.
	# Advantage of this version is that you can use R to process the depth files.

	module load samtools

	BAMLIST=/hb/groups/bernardi_lab/may/DTR/population-analysis/sample_lists/bam_path_list.txt #$1 # Path to a list of merged bam files. Full paths should be included
	MINBASEQ=20 # Minimum base quality filter
	MINMAPQ=20 # Minimum mapping quality filter

	SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID},1p" $BAMLIST)
	## Count per position depth per sample
	samtools depth -aa $SAMPLEBAM -q $MINBASEQ -Q $MINMAPQ | cut -f 3 |gzip > ${SAMPLEBAM}.depth.gz
	
Now we summarize the depth info \
For this we use R which means we need 2 scripts. 1, is the slurm job sumbission and the second is the actual R script. We also need to load all the necessary packages into the tidyverse conda environment. so just activate tidyverse and load the condas for `r-r.utils` `r-data.table` `r cowplot` \

`2-summmarize_depth_per_position.mpi` 

	#!/bin/bash
	#
	#SBATCH --job-name 2sumreaddepth
	#SBATCH --output=%A_2sumdepth.out
	#SBATCH --error=%A_2sumdepth.err
	#SBATCH --mail-type=ALL
	#SBATCH --mail-user=mabrober@ucsc.edu
	#SBATCH --partition=128x24
	#SBATCH --nodes=1
	#SBATCH --time=9-00:00:00
	#SBATCH --mem=120GB
	#SBATCH --ntasks-per-node=20


	module load miniconda3.9
	conda activate tidyverse

	Rscript summarize_depth_per_position.r
	
The Rscript:\
`summmarize_depth_per_position.r`

	## This script is used to summarize all the individual bam.depth.gz files (output of samtools depth),
	# It will come up with summary statistics for each individual as well per position depth and presence/absence summed across all individuals

	args <- commandArgs(trailingOnly = TRUE)
	BAMLIST <- "/hb/groups/bernardi_lab/may/DTR/population-analysis/sample_lists/bam_path_list.tsv" #args[1] # Path to a list of merged bam files. Full paths should be included. An example of such a bam list is /workdir/cod/greenland-cod/sample_lists/bam_list_merged.tsv
	SAMPLETABLE <- "/hb/groups/bernardi_lab/may/DTR/population-analysis/sample_lists/sample_table_merged.tsv" # args[2] # Path to a sample table where the 1st column is the prefix of the MERGED bam files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The 5th column is population name and 6th column is the data type. An example of such a sample table is: /workdir/cod/greenland-cod/sample_lists/sample_table_merged.tsv
	BASEDIR <- "/hb/groups/bernardi_lab/may/DTR/population-analysis/bams/deduped_overlap_clipped_bams" # args[3] # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be written to separate subdirectories. An example for the Greenland cod data is: /workdir/cod/greenland-cod/

	library(tidyverse)
	library(data.table)
	bam_list <- read_tsv(BAMLIST, col_names = F)$X1
	bam_list_prefix <- str_extract(BAMLIST, "[^.]+")
	sample_table <- read_tsv(SAMPLETABLE)
	for (i in 1:length(bam_list)){
	  print(i)
	  print(pryr::mem_used())
	  depth <- fread(paste0(bam_list[i], ".depth.gz"))$V1
	  mean_depth <- mean(depth)
	  sd_depth <- sd(depth)
	  presence <- as.logical(depth)
	  proportion_of_reference_covered <- mean(presence)
	  if (i==1){
	    output <- data.frame(sample_seq_id=sample_table[i,1], mean_depth, sd_depth, proportion_of_reference_covered)
	    total_depth <- depth
	    total_presence <- presence
	  } else {
	    output <- rbind(output, cbind(sample_seq_id=sample_table[i,1], mean_depth, sd_depth, proportion_of_reference_covered))
	    total_depth <- total_depth + depth
	    total_presence <- total_presence + presence
	  }
	}
	write_tsv(output, paste0(bam_list_prefix, "_depth_per_position_per_sample_summary.tsv"))
	write_lines(total_depth, paste0(bam_list_prefix, "_depth_per_position_all_samples.txt"))
	write_lines(total_presence, paste0(bam_list_prefix, "_presence_per_position_all_samples.txt"))

	# Total Depth per Site across All Individuals (on server)
	total_depth <- fread(paste0(bam_list_prefix, "_depth_per_position_all_samples.txt"))
	total_presence <- fread(paste0(bam_list_prefix, "_presence_per_position_all_samples.txt"))
	total_depth_summary <- count(total_depth, by=V1)
	total_presence_summary <- count(total_presence, by=V1)
	
Also I just want to make sure I'm understanding the purpose because it's still fuzzy in my head and it helps me to step back and see it in writing. 
So the previous step was summarize_depth_per_position.R done and now I have files like:

	dedup_bams_depth_per_position_all_samples.txt  (list of coverage per position)
	dedup_bams_depth_per_position_per_sample_summary.tsv (table of with row per sample and columns for mean depth, sd depth, and proportion of reference covered)
	dedup_bams_presence_per_position_all_samples.txt (list of number of samples present at each position)
	
Now we use `2.1-summary_stats_plot.mpi` to further summarized stats and plot a histogram (output of read_count.Rmd). We use the histogram which should be a normal distribution (though might have a long tail) to figure out what our min adn max depth should be for the next steps of snp calling and filtering.  

SO: for this part, I did a bunch of test runs of just the "## Design filters for SNP calling " on the head node to adjust the min and max depth values. Apparently this part is kind of eyeballed. We're going for chopping off the end of the distributions to make it a normal distribution looking plot.

 `2.1-summary_stats_plot.mpi` This calls on the R script 
 
	#!/bin/bash
	#
	#SBATCH --job-name sumstat
	#SBATCH --output=%A-sumstat.out
	#SBATCH --error=%A-sumstat.err
	#SBATCH --mail-type=ALL
	#SBATCH --mail-user=mabrober@ucsc.edu
	#SBATCH --partition=128x24
	##SBATCH --nodes=1
	#SBATCH --time=7-00:00:00
	#SBATCH --mem=100GB
	#SBATCH --cpus_per-task=1


	module load miniconda3.9
	conda activate tidyverse

	Rscript summary_stats.r
	
`summary_stats.r` The R script:

	#!/bin/bash
	#
	#SBATCH --job-name sumstat
	#SBATCH --output=%A-sumstat.out
	#SBATCH --error=%A-sumstat.err
	#SBATCH --mail-type=ALL
	#SBATCH --mail-user=mabrober@ucsc.edu
	#SBATCH --partition=128x24
	##SBATCH --nodes=1
	#SBATCH --time=7-00:00:00
	#SBATCH --mem=100GB
	#SBATCH --cpus_per-task=1


	module load miniconda3.9
	conda activate tidyverse

	Rscript summary_stats.r
	(tidyverse) [mabrober@hb fresh]$ cat summary_stats.r
	## Summary statistics
	library(tidyverse)
	library(cowplot)
	library(knitr)

	basedir <- '/hb/groups/bernardi_lab/may/DTR/population-analysis'

	# Read in the data

	per_position_per_ind <- read_tsv(paste0(basedir, "/sample_lists/dedup_bam_paths_depth_per_position_per_sample_summary.tsv"))
	depth_hist <- read_tsv(paste0(basedir,"/sample_lists/dedup_bam_paths_depth_per_position_all_samples_histogram.tsv"))
	presence_hist <- read_tsv(paste0(basedir,"/sample_lists/dedup_bam_paths_presence_per_position_all_samples_histogram.tsv"))

	# mean depth per position across all individuals
	mean_depth_per_position_all <- round(sum(depth_hist$by*depth_hist$n)/sum(depth_hist$n),2)

	# standard deviation
	sd <- sqrt(sum((depth_hist$by-254.58)^2*depth_hist$n)/sum(depth_hist$n))

	# total sites covered at least once
	sites_cov <- sum(filter(depth_hist, by>0)$n)
	# % of the reference genome
	percent_cov <- round((sites_cov)/sum(depth_hist$n)*100,2)

	# Put this in writing - find this in the slurm .out file
	cat("The mean depth per position across all individuals is", mean_depth_per_position_all, "and the standard deviation is", sd)
	cat(". A total of", sites_cov, "sites were covered at least once. This is", percent_cov, "% of the reference genome")

	# Let's plot shit - not totally necessary but interesting

	png(paste0(basedir,"/sample_lists/presence.png"),  width = 2000, height = 1500, units = "px")
	ggplot(presence_hist) +
	  geom_point(aes(x=by, y=n)) +
	  theme_cowplot()
	dev.off()

	png(paste0(basedir,"/sample_lists/depth.png"),  width = 2000, height = 1500, units = "px")
	ggplot(depth_hist) +
	  geom_point(aes(x=by, y=n)) +
	  theme_cowplot()
	dev.off()

	png(paste0(basedir,"/sample_lists/bases_vs_prop_ref_covered.png"),  width = 2000, height = 1500, units = "px")
	ggplot(per_position_per_ind) +
	  geom_point(aes(x=mean_depth, y=proportion_of_reference_covered)) +
	  theme_cowplot()
	dev.off()

	## Design filters for SNP calling

	## Super low and super high coverage sites are cut from this figure
	png(paste0(basedir,"/sample_lists/super_low_and_super_high_cut1.png"),  width = 2000, height = 1500, units = "px")

	filter(depth_hist, by>10, by<700) %>%
	  ggplot(aes(x=by, y=n)) +
	  geom_freqpoly(stat = "identity") +
	  theme_cowplot()

	## 176 is the mode of the second peak
	filter(depth_hist, by>10) %>%
	  arrange(by=n) %>%
	  tail(n=1)

	filter(depth_hist, by<345) %>% #176# by<x, where x = mode of the second peak
	  arrange(by=n) %>%
	  head(n=1)

	mode_line <- 345
	max_depth <- 590 #345-100+345 #350 # try and get this so that it makes the depth a normal distribution
	paste0('If MaxDepth=', max_depth, ' filter is used, ', round(sum(filter(depth_hist, by>max_depth)$n)/sum(depth_hist$n)*100,2), '% of all sites and ', round(sum(filter(depth_hist, by>max_depth)$n*filter(depth_hist, by>max_depth)$by)/sum(depth_hist$n*depth_hist$by)*100,2), '% of the final mapped data will be lost')

	min_depth <- 100#39 # Set this to be the first trough
	paste0('If minDepth=', min_depth, ' filter is used, ', round(sum(filter(depth_hist, by<min_depth)$n)/sum(depth_hist$n)*100,2), '% of all sites and ', round(sum(filter(depth_hist, by<min_depth)$n*filter(depth_hist, by<min_depth)$by)/sum(depth_hist$n*depth_hist$by*100,2)), '% of the final mapped data will be lost')

	## If these filters are used
	filter(depth_hist, by>1, by<700) %>%
	ggplot(aes(x=by, y=n)) +
	  geom_freqpoly(stat = "identity") +
	  geom_vline(xintercept = c(max_depth,min_depth), color="red") +
	  theme_cowplot()

	dev.off()

Now we call gloabal (for all samples regardless of population) snps using the snps we determined in previous step\

`1-call_global_snps_all.mpi`

	#!/bin/bash
	#
	#SBATCH --job-name glblsnpcall_angsd
	#SBATCH --output %A_glblsnpcall_angsd.out
	#SBATCH --error %A_glblsnpcall_angsd.err
	#SBATCH --mail-type=ALL
	#SBATCH --mail-user=mabrober@ucsc.edu
	#SBATCH --time=15-00:00:00
	#SBATCH --ntasks=1
	#SBATCH --mem=120GB   #

	## This script is used to call SNPs using angsd
	module load angsd

	BAMLIST=dedup_bam_paths.tsv # Path to textfile listing bamfiles to include in global SNP calling with absolute paths
	LIST_DIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/sample_lists/$BAMLIST
	#POPSLIST=populations.txt
	BASEDIR=/hb/groups/bernardi_lab/may/DTR/population-analysis/ # Path to the base directory where output files will be written to a subdirectory named "angsd/". An example for the Greenland cod data is: /workdir/cod/greenland-cod/
	REFERENCE=/hb/groups/bernardi_lab/may/DTR/DTRgenome/genome_annot/genome/kuro_filt_s500.fasta # Path to reference genome
	MINDP=100 # Minimum depth filter 
	MAXDP=590 # Maximum depth filter 
	MININD=13 #10% of individuals
	MINQ=20  #Minimum quality filter #The minimum base quality score
	MINMAF=0.5 #Minimum minor allele frequency filter
	MINMAPQ=20 ### Minimum mapping quality (alignment score) filter, default value is 20
	EXTRA_ARG='-remove_bads 1 -only_proper_pairs 1 -C 50' # Extra arguments when running ANGSD, default value is '-remove_bads 1 -only_proper_pairs 1 -C 50'
	#THREADS=${11:-8} # Number of parallel threads to use, default value is 8.
	#-SNP_pval      1.000000        (Remove sites with a pvalue larger)


	## Extract the name of the bam list (excluding path and suffix)
	BAMLISTNAME=`echo $BAMLIST | sed 's/\..*//'` #| sed -e 's#.*/\(\)#\1#'`

	## Build base name of output files
	OUTBASE=$BAMLISTNAME'_mindp'$MINDP'_maxdp'$MAXDP'_minind'$MININD'_minq'$MINQ

	## Call SNPs #GL 1 is the samtools model
	angsd -b $LIST_DIR -ref $REFERENCE -out $BASEDIR'angsd'/$OUTBASE \
		-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 10000 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
		-setMinDepth $MINDP -setMaxDepth $MAXDP -minInd $MININD \
		-minQ $MINQ -minMapQ $MINMAPQ \
		-SNP_pval 1e-6 -minMaf $MINMAF \
		$EXTRA_ARG
		>& $BASEDIR'angsd/'$OUTBASE'.log'

	## Create a SNP list to use in downstream analyses
	gunzip -c $BASEDIR'angsd/'$OUTBASE'.mafs.gz' | cut -f 1,2,3,4 | tail -n +2 > $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt'
	angsd sites index $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt'

	## Also make it in regions format for downstream analyses
	cut -f 1,2 $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt' | sed 's/\t/:/g' > $BASEDIR'angsd/global_snp_list_'$OUTBASE'.regions'

	## Lastly, extract a list of chromosomes/LGs/scaffolds for downstream analysis
	cut -f1 $BASEDIR'angsd/global_snp_list_'$OUTBASE'.txt' | sort | uniq > $BASEDIR'angsd/global_snp_list_'$OUTBASE'.chrs'
	
  
