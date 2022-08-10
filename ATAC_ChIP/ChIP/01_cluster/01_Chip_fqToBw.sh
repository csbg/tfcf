#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=Chip_fqToBw
#SBATCH --cpus-per-task=12
#SBATCH --mem=30G
#SBATCH --time=00-08:00:00
#SBATCH -p medium
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 
##SBATCH --dependency=afterany:624523


##SBATCH --mail-type=END
##SBATCH --mail-user=user@mail.es
# HOW TO RUN ME
# for i in *fastq.gz; do echo $i | sed 's/_R._001.fastq.gz//g' ; done | sort | uniq > samplesNames.txt  
# N=`cat samplesNames.txt | wc -l`
# sbatch --array=1-${N} /PATH/TO/SCRIPTS/FOLDER/01_Chip_fqToBw.sh \
#/PATH/TO/FASTQ/TREE/FOLDER \
#/PATH/TO/OUTPUT \
#/PATH/TO/referenceGenomes/mm10_reordered/mm10.reordered


# UPDATES
# it works with gziped fastq files
# not valid for atac seq, dosnt take into account the 9bp mapping sift (Tn5, hichikers paper)
# i can use nextflow or prepare time points in a file in case the code breaks
##===============================================================================
## GLOBAL VARIABLES
PROJECT_DIR=$1
# In my case, the fastq files are in project folder inside demux_fastq

RAW_FASTQ_DIR=$PROJECT_DIR"/demux_fastq"
# not used, but i will store here final files
EDITED_DIR=$PROJECT_DIR"/pipelineOut"
FASTQ_DIR=$EDITED_DIR"/fastq" 
libMetricsP="${EDITED_DIR}/libMetrics"
# set to lowercase yes or not
removeTemp="yes"

## path to the main fully processed data output
basePath=$2
# here we will stored the final filtered bam files
bamsPath="${basePath}/bamfiles"


##  REFERENCE Genome info
REFERENCE_DIR=$3
# REFERENCE_DIR includes the path and naming base of the reference genome, we will
#   add .sizes, .blacklist.bed to the base to get the rest of paths. Genome
#   indexes must follow the same base name

GenomeIndex=$REFERENCE_DIR
chr_genome_size=$REFERENCE_DIR".sizes"
BlackList=$REFERENCE_DIR".blacklist.bed"
wigToBigWig="/PATH/TO/wigToBigWig"
picardPath='/PATH/TO/picard.jar'

##===============================================================================
## Required Software
#module load Trimmomatic/0.38-Java-1.8
#module load Trim_Galore/0.6.0
# trim_galore will requilre fastqc and cutadapt
# I have trim galore and cutadapt in conda
module load FastQC/0.11.8-Java-1.8
module load Bowtie2/2.3.4.2-foss-2018b
module load SAMtools/1.12-GCC-10.2.0
module load picard/2.18.17-Java-1.8
module load R/4.0.0-foss-2018b
module load Java/1.8.0_192
module load BEDTools/2.27.1-foss-2018b
#module load MACS2/2.1.0.20151222-foss-2018b-Python-2.7.15
#module load deepTools/3.2.0-foss-2018b-Python-2.7.15
# this is for bamCoverage, but i already have it in conda and newer version
#module load Homer/4.10-foss-2018b

##===============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

##===============================================================================
if [ ! -e ${EDITED_DIR} ]; then
    mkdir -p ${EDITED_DIR}
fi


cd $EDITED_DIR

#Choose file
# If we always get in list all files from folder
#FILES=($(ls $RAW_FASTQ_DIR | grep "fastq.gz" | sed -r 's/.{16}$//' | uniq))
# If we want to do only the files in samplesNames.txt file. Perfecto for re-running 
#dead jobs
FILES=($(cat $RAW_FASTQ_DIR/samplesNames.txt))
filename=${FILES[$SLURM_ARRAY_TASK_ID - 1]}
echo -e $SLURM_JOB_NODELIST 
#https://slurm.schedmd.com/faq.html#task_prolog
echo "print =========================================="
echo "print SLURM_JOB_ID = $SLURM_JOB_ID"
echo "print SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
echo "print =========================================="
echo $filename
##===============================================================================

# get some paths
read1_path="${RAW_FASTQ_DIR}/${filename}_R1_001.fastq.gz"
read2_path="${RAW_FASTQ_DIR}/${filename}_R2_001.fastq.gz"
stepControl="${EDITED_DIR}/QC/pipelineStep_${filename}.txt"
summaryFile="${basePath}/QC/summary_${filename}.txt"

#############
# BEGIN: Create summary file
#############								
if [ ! -e ${EDITED_DIR}/QC/ ] && echo exists ; then
    mkdir -p ${EDITED_DIR}/QC/
fi

if [ ! -e ${basePath}/QC/ ] && echo exists ; then
    mkdir -p ${basePath}/QC/
fi

if [ ! -e ${stepControl} ] ; then
    touch ${stepControl}
fi

echo -e "Starting Summary file -------------------------------------- \n"

# check content of first line of step control file
linec=`sed "1q;d" ${stepControl}`
if [[ ${linec} != "Summary" ]]; then 
    echo -e "STARTING \n $(date) \n" 
    echo "SAMPLE: ${filename}" 
    echo -e "sample name\tfastq name\tread count\tmillions" >> ${summaryFile}

    # QC: read counts if file is gziped or not
    if file --mime-type ${read1_path} | grep -q gzip$; then
        Counts1="$(zcat ${read1_path} | echo $((`wc -l`/4)))"
        Counts2="$(zcat ${read2_path} | echo $((`wc -l`/4)))" 
    else
        echo "Not gzipped files"
        echo $read1_path
        #Counts1="$(cat ${RAW_FASTQ_DIR}/${filename}_R1_001.fastq | echo $((`wc -l`/4)))"
        #Counts2="$(cat ${RAW_FASTQ_DIR}/${filename}_R2_001.fastq| echo $((`wc -l`/4)))"
    fi


    echo -e "READ COUNTS \n" >> ${summaryFile}
    rc=$((Counts1/1000000))
    echo -e "${filename} \t ${filename}_R1 \t ${Counts1} \t ${rc}" >> ${summaryFile}
    rc=$((Counts2/1000000))
    echo -e "${filename} \t ${filename}_R2 \t ${Counts2} \t ${rc} \n" >> ${summaryFile}

    echo -e "Summary file - done -------------------------------------- \n"
    # store stage control info
    echo "Summary" > ${stepControl}
else
    echo -e "Summary file - already done before -------------------------------------- \n"
fi

#############
# TR: Trimming
#############
# We have to trim for repetitive seq, adapters, quality etc.
echo -e "Starting Trimming and FASTQC -------------------------------------- \n"

# check content of second line of step control file
linec=`sed "2q;d" ${stepControl}`
if [[ ${linec} != "Trim" ]]; then 
    if [ ! -e ${EDITED_DIR}/fastQC/ ]; then
        mkdir -p ${EDITED_DIR}/trimming/  
        mkdir -p ${EDITED_DIR}/fastQC/
    fi

    if [ ! -e ${EDITED_DIR}/fastQC/${filename}_R2_001.fastq_fastqc.html ]; then
        trim_galore --cores $SLURM_CPUS_PER_TASK --paired --fastqc --gzip \
        --output_dir ${EDITED_DIR}/trimming/ --fastqc_args "--outdir ${EDITED_DIR}/fastQC/" \
        ${read1_path} ${read2_path}
    fi
    echo -e "Trimming - done -------------------------------------- \n"

    # store stage control info
    echo "Trim" >> ${stepControl}

else
    echo -e "Trimming - already done before -------------------------------------- \n"
fi

# define new path variables
trimedRead1="${EDITED_DIR}/trimming/${filename}_R1_001_val_1.fq.gz"
trimedRead2="${EDITED_DIR}/trimming/${filename}_R2_001_val_2.fq.gz"

#############
#Aling # Mapping
############
echo -e "Starting Alignment -------------------------------------- \n"

if [ ! -e ${EDITED_DIR}/BAM/ ]; then
    mkdir -p ${EDITED_DIR}/BAM/
    mkdir -p ${EDITED_DIR}/unMapped/
fi

samFile="${EDITED_DIR}/BAM/${filename}.sam"

# check content of third line of step control file
extr_unmap="${EDITED_DIR}/unMapped/${filename}_R%_001_unmap.fastq.gz"
linec=`sed "3q;d" ${stepControl}`
if [[ ${linec} != "Align" ]]; then 
    bowtie2 -p $SLURM_CPUS_PER_TASK -X 1000 --no-discordant --no-mixed \
    --very-sensitive -x $GenomeIndex -1 ${trimedRead1} -2 ${trimedRead2} \
    -S ${samFile} --un-conc-gz ${extr_unmap}
        
    echo -e "Alignment - done -------------------------------------- \n"
    # store stage control info
    echo "Align" >> ${stepControl}

    # remove temporal trimmed fastq files
    if [[ $removeTemp == 'yes' ]] ; then
        rm ${trimedRead1}
        rm ${trimedRead2}
    fi
else
    echo -e "Alignment - already done before ------------------------------ \n"
fi

# SAM to BAM (remove SAM)
echo -e "Starting SAM to BAM -------------------------------------- \n"

bamFile="${EDITED_DIR}/BAM/${filename}.bam"

# check content of fourth line of step control file
linec=`sed "4q;d" ${stepControl}`
if [[ ${linec} != "SamBam" ]]; then 
    samtools view -o ${bamFile} -bhS -@ $SLURM_CPUS_PER_TASK \
    ${samFile} 

    echo -e "SAM to BAM - done -------------------------------------- \n"
    # store stage control info
    echo "SamBam" >> ${stepControl}
else
    echo -e "SAM to BAM - already done before ------------------------------ \n"
fi


# Sort BAM by location (remove original BAM)

echo -e "Starting BAM sorting --------------------------------------\n"

bamSort="${EDITED_DIR}/BAM/${filename}.sort.bam"

# check content of fifth line of step control file
linec=`sed "5q;d" ${stepControl}`
if [[ ${linec} != "BamSort" ]]; then 
    samtools sort -o ${bamSort} ${bamFile} 
    echo -e "BAM sorting - done --------------------------------------\n"
    # store stage control info
    echo "BamSort" >> ${stepControl}

    # delete intermediate files
    if [[ $removeTemp == 'yes' ]] ; then
        rm ${samFile}
        rm ${bamFile}
    fi
else
    echo -e "BAM sorting - already done before -----------------------------\n"
fi


##############################
#Dupicates marking and removal
##############################
echo -e "Starting mark and remove duplicates ------------------------------------------------ \n"
# Samtools flags to remove duplicates + unmapped reads:
#   a) -F 1028 => read unmapped // duplicate
#   b) -F 1804 => read unmapped // mate unmapped // not primary alignment // read failing platform // duplicate

if [ ! -e ${libMetricsP} ]; then
    mkdir -p ${libMetricsP}
fi

bamSortMarkDup="${EDITED_DIR}/BAM/${filename}.sort.markdup.bam"
bamSortRmDup="${EDITED_DIR}/BAM/${filename}.sort.rmdup.bam"

# check content of sixth line of step control file
linec=`sed "6q;d" ${stepControl}`
if [[ ${linec} != "RmDup" ]]; then 
    # Estimate the numbers of unique molecules in a sequencing library
    java -Xmx24G -jar ${picardPath} EstimateLibraryComplexity I=${bamSort} \
        O=${libMetricsP}/${filename}.bam.lib_comp.txt VALIDATION_STRINGENCY=SILENT
    # Check insert size distribution and read orientation
    java -Xmx24G -jar ${picardPath} CollectInsertSizeMetrics I=${bamSort} \
        O=${libMetricsP}/${filename}.bam.insert_size.txt H=${libMetricsP}/${filename}.bam.insert_size.pdf VALIDATION_STRINGENCY=SILENT
    # Identify duplicated reads by flaging them in a new ioutput BAM file
    java -Xmx24G -jar ${picardPath} MarkDuplicates I=${bamSort} O=${bamSortMarkDup} \
        METRICS_FILE=${libMetricsP}/${filename}.bam.mkdup.txt
    # remove reads marked as duplicated
    samtools view -o ${bamSortRmDup} -@ $SLURM_CPUS_PER_TASK -bh -F 1804 ${bamSortMarkDup}
    # Check insert size distribution and read orientation after removing duplicates
    java -Xmx24G -jar ${picardPath} CollectInsertSizeMetrics I=${bamSortRmDup} \
        O=${libMetricsP}/${filename}.rmdup.bam.insert_size.txt H=${libMetricsP}/${filename}.rmdup.bam.insert_size.pdf VALIDATION_STRINGENCY=SILENT

    # remove orphan reads?
    # bampe_rm_orphan.py ${bam[0]} ${prefix}.bam --only_fr_pairs

    # QC: Show duplicates in samtools flags
    echo -e "SAMTOOLS FLAGSTAT - DUPLICATES \n" >> ${summaryFile}

    samtools flagstat ${bamSortMarkDup} >> ${summaryFile}

    echo -e "\n"  >> ${summaryFile}
   
    if [[ $removeTemp == 'yes' ]] ; then
        rm ${bamSort}
        rm ${bamSortMarkDup}
    fi
    echo -e "Remove duplicates - done ------------------------------------------------ \n"
    # store stage control info
    echo "RmDup" >> ${stepControl}
else
    echo -e "Remove duplicates - already done before ------------------------------------ \n"

fi 

echo -e "DATE CHECK \n $(date) \n" 



###################################
# BlackList Filter # Remove Blacklist regions from BAM
#################################
echo -e "Starting Remove blacklist regions -----------------------------------------\n"

bamSortMarkDupBlack="${EDITED_DIR}/BAM/${filename}.sort.rmdup.rmblackls.bam"

# check content of seventh line of step control file
linec=`sed "7q;d" ${stepControl}`
if [[ ${linec} != "Blacklist" ]]; then 
    bedtools intersect -v -abam ${bamSortRmDup} -b ${BlackList} > ${bamSortMarkDupBlack}
    if [[ $removeTemp == 'yes' ]] ; then
        rm ${bamSortRmDup}
    fi
    echo -e "Remove blacklist regions - done -----------------------------------------\n"
    # store stage control info
    echo "Blacklist" >> ${stepControl}
else
    echo -e "Remove blacklist regions - already done before ------------------------------\n"
fi




# Remove chrM, chrUn, _random ... (remove sort.rmdup.rmblackls.bam)
echo -e "Starting remove chrM and useless chromosomes ------------------------------\n"

if [ ! -e ${bamsPath} ]; then
    mkdir -p ${bamsPath}
fi

bamSortMarkDupBlackChr="${bamsPath}/${filename}.sort.rmdup.rmblackls.rmchr.bam"

# check content of eigth line of step control file
linec=`sed "8q;d" ${stepControl}`
if [[ ${linec} != "Remove" ]]; then 
	samtools view -h ${bamSortMarkDupBlack} | \
    awk '(!index($3, "random")) && (!index($3, "chrUn")) && ($3 != "chrM") && ($3 != "chrEBV")' | \
    samtools view -Sb - > ${bamSortMarkDupBlackChr}

    # QC: Show final reads
    echo -e "SAMTOOLS FLAGSTAT - FINAL READS \n" >> ${summaryFile}
    samtools flagstat ${bamSortMarkDupBlackChr} >> ${summaryFile}
    echo -e "\n" >> ${summaryFile}

    echo -e "Remove chrM and useless chromosomes - done ----------------------------------\n"
    # store stage control info
    echo "Remove" >> ${stepControl}
    if [[ $removeTemp == 'yes' ]] ; then
        rm ${bamSortMarkDupBlack}
    fi
else
    echo -e "Remove chrM and useless chromosomes - already done before--------------------------\n"
fi


# Index BAM
echo -e "Starting Index BAM ------------------------------------------------------\n"

# check content of ninth line of step control file
linec=`sed "9q;d" ${stepControl}`
if [[ ${linec} != "IndexBam" ]]; then 
	samtools index ${bamSortMarkDupBlackChr} -@ $SLURM_CPUS_PER_TASK
    echo -e "Index BAM - done ------------------------------------------------------\n"
    # store stage control info
    echo "IndexBam" >> ${stepControl}
else
    echo -e "Index BAM - already done before --------------------------------------\n"
fi


###################################
# BigWigs 
#################################
# deeptools bamCoverage
# Normalize by CPM (This is the scaled bigWig we will use)
echo -e "Starting BigWigs --------------------------------------------------\n"

if [ ! -e ${basePath}/BigWig/ ]; then
	mkdir ${basePath}/BigWig/
fi

if [ ! -e ${EDITED_DIR}/BigWig/ ]; then
	mkdir ${EDITED_DIR}/BigWig/
fi

bigWigOut="${basePath}/BigWig/${filename}.sort.rmdup.rmblackls.rmchr.norm.bw"

# check content of tenth line of step control file
linec=`sed "10q;d" ${stepControl}`
if [[ ${linec} != "BigWnorm1" ]]; then 
	bamCoverage --binSize 5 --normalizeUsing CPM --exactScaling \
    -b ${bamSortMarkDupBlackChr} -of bigwig \
    -o ${bigWigOut} --numberOfProcessors $SLURM_CPUS_PER_TASK
	
    echo -e "BigWig norm 1 - done ---------------------------------------------\n"
    # store stage control info
    echo "BigWnorm1" >> ${stepControl}
else
    echo -e "BigWig norm 1 - already done before ------------------------------\n"
fi


echo -e "END --------------------------------------------------"

seff $SLURM_JOBID

exit 0