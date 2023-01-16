#!/bin/bash
# -*- ENCODING: UTF-8 -*-

#####===============================================================================
####   MECC_sgRNA_v1.sh
####   Counts for CrisprScreens
#####===============================================================================

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=CRISPR
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH -p short
#SBATCH -o /home/X/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/X/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
# for i in *_R1_001.fastq.gz; do echo $i | sed 's/_R1_001.fastq.gz//g' ; done | \
#sort | uniq > samplesNames.txt
# N=`cat samplesNames.txt | wc -l`
# sbatch --array=1-${N} /home/X/programas/PhD/CRISPR-screen/00_MECC_sgRNA_v1.1.sh \
#/PATH/TO/FOLDER/WITH/FASTQ \
#/PATH/TO/OUTPUT/FOLDER \
#PATH/TO/LIBRARY/REFERENCE/GENOMES \
#'CACCG(.{20})GT'

# the files stored in samplesNames.txt are the ones that will be analysed
# Files must have the S[0-9]_R[12]_001.fastq.gz structure, otherwise R script in
# next step will fail with rbind message

# Fastq file naming
# [Cell]_[extra]_[libraryID]_[date?]_... 
# The most important is to have the CRISPR library ID in the 3rd section
#This will tell the script against which libraries reference to align
# Allowed library IDs for now:
#As  BBr  LibP1  MGLibA  OP1234  OP3  R1Br  TF1Br  TF3Br 
#B   Br   LibP2  OP1     OP2     OP4  R2Br  TF2Br

##===============================================================================
## GLOBAL VARIABLES
# the ones that are given at script run

# -- base path were fastqs are located (inside a folder named demux_fastq/)
PARENT_DIR=$1
#PARENT_DIR="/home/X/data/2021/CRISPR/sequencedData/merge4-492"
final_dir=$2
#final_dir="/home/X/data/2021/CRISPR/allProcessed/merge4-492"

## Guide REFERENCE 
## These are the indexes for the alignment of the sgRNAs. 
## We will have just one genome index for every library and
# an index for all the libraries together allGuides/finalGuides.fa
## So here we pass the main folder with
# each of the IDs named ID/ID.fa as in the guide ID stated in every sample
GenomeIndex_all=$3
#/home/X/referenceGenomes/sgRNA_indexes/bowtie2

# bulk library design (More info about pattern in line 54)
flanquingSeq=$4
#flanquingSeq='CACCG(.{20})GT'

# if the fastq files have 001 or not, we add the string or not
extraSTR="_001"
#extraSTR=""


##===============================================================================
## FIXED USER VARIABLES

# path to the location of my git repo
gitPath="/home/X/programas/PhD"

# these ones shouldnt need to be modified if we are using the right folder structure
FASTQ_DIR=$PARENT_DIR"/demux_fastq"
EDITED_DIR=$PARENT_DIR"/pipelineOut"
#FASTQC_DIR=$EDITED_DIR"/fastQC"

# path to the python script to extract the guide sequence
# the script will rely on the first '(.{' symbol in the pattern to get 
#the number of nucleotides left-surrounding our set of nucleotides of interest
# the script will also assume that the number of nucleotides of interest 
#is enclosed between '(.{' and '})', and get accordingly the (.{n}) 
#nucleotides after the nucleotides preceding '(.{' in the matched string
# Hence, adding a filter pattern with '(.{' at the left of when we ask for 
#our n nucleotides will cause fail
# first input is GZIP-ed fastq file path, then you can optionally add a regexp pattern 
#(default is 'CACCG(.{20})GTTTTAGAGC') and 'True' if you want the output fastq file
#to be also GZIP-ed
extractScript="${gitPath}/CRISPR-screen/00_NR_CRISPR-extract.py"
outGZ="False"

# path to R script for final guide analysis
RscriptP="${gitPath}/CRISPR-screen/01a_MECC_sgRNA_pre_analysis_v2.r"

# path to the python mapping report script
reportScript="${gitPath}/CRISPR-screen/01b_mergedMappingReport.py"


##===============================================================================
## Required Software
#module load FastQC/0.11.8-Java-1.8
#module load MultiQC/1.7-foss-2018b-Python-2.7.15
#module load Bowtie/1.2.2-foss-2018b
module load Bowtie2/2.3.4.2-foss-2018b
module load SAMtools/1.9-foss-2018b
#module load Perl/5.28.0-GCCcore-7.3.0

# set python paths
export PATH="/home/X/programas/miniconda3/bin:$PATH"
export PYTHONPATH=/home/X/programas/miniconda3/bin/python3.8
##===============================================================================

export PS4='$LINENO+ '
set -x

echo Maren Calleja  
echo me.callejac@gmail.com  mcallejac@unav.es
echo updated
echo julenmendieta92@gmail.com  jmendietaes@unav.es
echo Date ; date
echo script -- aln -- counts
echo LibraryType=SE
echo genome = mm10

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

echo -e $SLURM_JOB_NODELIST 
#https://slurm.schedmd.com/faq.html#task_prolog
echo "print =========================================="
echo "print SLURM_JOB_ID = $SLURM_JOB_ID"
echo "print SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
echo "print =========================================="

##===============================================================================
## Create or move around the corresponding directories
START_TIME=$SECONDS

if [ ! -e ${EDITED_DIR} ]; then
    mkdir -p ${EDITED_DIR}
    mkdir -p ${EDITED_DIR}/fastq_extr
fi

if [ ! -e ${EDITED_DIR}/QC/ ] && echo exists ; then
    mkdir -p ${EDITED_DIR}/QC/
fi

if [ ! -e ${final_dir}/QC/ ] && echo exists ; then
    mkdir -p ${final_dir}/QC/
fi

cd $EDITED_DIR
##===============================================================================
##Choose files from samplesNames.txt file
FILES=($(cat $FASTQ_DIR/samplesNames.txt))
filename=${FILES[$SLURM_ARRAY_TASK_ID - 1]}
sampleRunName=`basename ${final_dir}`

# get library ID to which map this file
mapLib=(${filename//_/ })
mapLib=${mapLib[2]}
# get path to focus library
GenomeIndex="$GenomeIndex_all/$mapLib/$mapLib.fa"
# get path to all libraries used in the lab (contamination control)
GenomeIndex_allGuide="$GenomeIndex_all/allGuides/finalGuides.fa"

# ensure that we have this index
if [[ ! -f ${GenomeIndex} ]]
then
    echo "Inferred index file does not exist"
    echo ${GenomeIndex}
    exit 1
fi

# get some paths
read1_path="${FASTQ_DIR}/${filename}_R1${extraSTR}.fastq.gz"
stepControl="${EDITED_DIR}/QC/pipelineStep_${filename}.txt"
summaryFile="${final_dir}/QC/summary_${filename}.txt"

if [ ! -e ${stepControl} ] ; then
    touch ${stepControl}
fi

##===============================================================================
#############
## BEGIN: create summary file
#############

echo -e "Starting Summary file -------------------------------------- \n"

# check content of first line of step control file
linec=`sed "1q;d" ${stepControl}`
if [[ ${linec} != "Summary" ]]; then 
    echo -e "STARTING \n $(date) \n" 
    echo "SAMPLE: ${filename}" 
    echo -e "READ COUNTS \n" >> ${summaryFile}
    echo -e "sample name\tfastq name\tread count\tmillions" >> ${summaryFile}

    # QC: read counts if file is gziped or not
    if file --mime-type ${read1_path} | grep -q gzip$; then
        Counts1="$(zcat ${read1_path} | echo $((`wc -l`/4)))"
    else
        echo "Not gzipped files"
        echo $read1_path
        exit 1
        #Counts1="$(cat ${RAW_FASTQ_DIR}/${filename}_R1_001.fastq | echo $((`wc -l`/4)))"
        #Counts2="$(cat ${RAW_FASTQ_DIR}/${filename}_R2_001.fastq| echo $((`wc -l`/4)))"
    fi

    rc=$((Counts1/1000000))
    echo -e "${sampleRunName} \t ${filename} \t ${Counts1} \t ${rc}" >> ${summaryFile}
    
    echo -e "Summary file - done -------------------------------------- \n"
    # store stage control info
    echo "Summary" > ${stepControl}
else
    echo -e "Summary file - already done before -------------------------------------- \n"
fi


#############
# Extract: Extraction of guides sequence
#############

# Extraction
echo -e "Starting CRISPR-extract-carpools------------------- \n"

# Create output dir
if [[ ${outGZ} == "False" ]]; then
    extr_fastq="${EDITED_DIR}/fastq_extr/${filename}_R1${extraSTR}_extracted.fastq"
    extr_unmap="${EDITED_DIR}/fastq_extr/${filename}_R1${extraSTR}_extracted_unmap.fastq"
else
    extr_fastq="${EDITED_DIR}/fastq_extr/${filename}_R1${extraSTR}_extracted.fastq.gz"
    extr_unmap="${EDITED_DIR}/fastq_extr/${filename}_R1${extraSTR}_extracted_unmap.fastq.gz"
fi

# Python script arguments
#	[Regular expression used to extract 20 bp target sequence]
#	[FastqFile] compressed, must finish in .gz
#	[True to GZIP reulting fastq]


# check content of second line of step control file
linec=`sed "2q;d" ${stepControl}`
if [[ ${linec} != "Extract" ]]; then 
    echo -e "\nGUIDE EXTRACTION\n" >> ${summaryFile}
    python ${extractScript} --fastq_path ${read1_path} --pattern ${flanquingSeq} --outgz ${outGZ} >> ${summaryFile}
    mv ${FASTQ_DIR}/${filename}_R1${extraSTR}_extracted.fastq ${EDITED_DIR}/fastq_extr/

    echo -e "CRISPR-extract-carpools - done -------------------------------------- \n"
    # store stage control info
    echo "Extract" > ${stepControl}
else
    echo -e "CRISPR-extract-carpools - already done before ---------------------- \n"
fi



#############
# Alignemnt
#############

# Alignment and postprocess

echo -e "Starting Alignment -------------------------------------- \n"

## David comment: This is the alignment process with Bowtie2, usually we are above 90% aligment rate, it would be good to flag the samples that do not yield such %aligmmnet

if [ ! -e ${EDITED_DIR}/BAM/ ]; then
    mkdir -p ${EDITED_DIR}/BAM/
fi


samPath="${EDITED_DIR}/BAM/${filename}.sam"
samPath_unM="${EDITED_DIR}/BAM/${filename}.unMap.sam"
# check content of third line of step control file
linec=`sed "3q;d" ${stepControl}`
if [[ ${linec} != "Align" ]]; then 
    # For relatively short reads (e.g. less than 50 bp) Bowtie 1 is sometimes 
    #faster and/or more sensitive.
    # http://bowtie-bio.sourceforge.net/bowtie2/faq.shtml
    #bowtie -p $SLURM_CPUS_PER_TASK ${GenomeIndex} ${extr_fastq} -S ${samPath} --best

    # The mapping takes 22 secods, mapping percentaje difference goes from 96.58 
    #in bowtie2 to 97.23 in bowtie but I like that bowtie 2 states number of 
    #alignments of reads in more than one ref
    echo -e "\nALIGNMENT\n" >> ${summaryFile}
    bowtie2 -p $SLURM_CPUS_PER_TASK -x $GenomeIndex -U ${extr_fastq} \
            -S ${samPath} --un ${extr_unmap} >> ${summaryFile} 2>&1
    
    # with --un we write unmapped reads to another fastq
    # then we try to map them to the reference genome with all the guides
    echo -e "\nALIGNMENT of previously unmapped\n" >> ${summaryFile}
    bowtie2 -p $SLURM_CPUS_PER_TASK -x $GenomeIndex_allGuide -U ${extr_unmap} \
            -S ${samPath_unM} >> ${summaryFile} 2>&1

    # Default mode: search for multiple alignments, report the best one
    # Bowtie2 does not guarantee that the alignment reported is the best 
    #possible in terms of alignment score.
    # Increasing -R makes Bowtie 2 slower, but increases the likelihood that 
    #it will report the correct alignment for a read that aligns many places.
    # store stage control info
    echo "Align" >> ${stepControl}

    echo -e "Alignment - done -------------------------------------- \n"
else
    echo -e "Alignment - already done before -------------------------------------- \n"
 
fi


#############
# toBam: convert to BAM and create sorted index
#############

echo -e "Starting SAM to BAM -------------------------------------- \n"

bamPath="${EDITED_DIR}/BAM/${filename}.bam"
bamSortPath="${EDITED_DIR}/BAM/${filename}.sort.bam"

bamPath_unM="${EDITED_DIR}/BAM/${filename}.unMap.bam"
bamSortPath_unM="${EDITED_DIR}/BAM/${filename}.unMap.sort.bam"
# check content of forth line of step control file
linec=`sed "4q;d" ${stepControl}`
if [[ ${linec} != "toBam" ]]; then 
    # first we filter out reads mapping more than twice
    #samtools view -h ${samPath} | grep -v XS:i > ${samPathUniq}
    #skipped=`samtools view ${samPath} | grep XS:i | wc -l`
    skipped=0
    echo -e "\nFILTER READS MAPPING MORE THAN ONCE\n" >> ${summaryFile}
    echo -e "$skipped reads filtered out\n" >> ${summaryFile}
    echo -e "Number of reads that map more than once, and their top target\n" >> ${summaryFile}
    samtools view -h ${samPath} | grep XS:i | awk '{print $3}' | sort | uniq -c >> ${summaryFile}
    
    # Then we convert to bam
    samtools view -o ${bamPath} -bhS -@ $SLURM_CPUS_PER_TASK ${samPath}
    samtools view -o ${bamPath_unM} -bhS -@ $SLURM_CPUS_PER_TASK ${samPath_unM}

    # I will discard this and retrieve idxstats ONLY
    # samtools view ${EDITED_DIR}/BAM/${filename}.bam| cut -f3 | sort | \
    #uniq -c > ${EDITED_DIR}/BAM/${filename}.cnt 
    # cat ${EDITED_DIR}/BAM/${filename}.cnt| awk '{print $1;}' > file2 
    # cat ${EDITED_DIR}/BAM/${filename}.cnt| awk '{print $2;}' > file1
    # paste file1 file2 > ${EDITED_DIR}/BAM/${filename}.txt
    # rm file1 file2
    # head=( ${EDITED_DIR}/BAM/${filename}.txt COUNTS )
    # ( IFS=$'\t'; echo "${head[*]}"; \
    #cat ${EDITED_DIR}/BAM/${filename}.txt ) > ${EDITED_DIR}/BAM/${filename}.final.txt

    samtools sort -o ${bamSortPath} ${bamPath} 
    samtools index -b ${bamSortPath}

    samtools sort -o ${bamSortPath_unM} ${bamPath_unM} 
    samtools index -b ${bamSortPath_unM}

    # store stage control info
    echo "toBam" >> ${stepControl}

    echo -e "SAM to BAM - done -------------------------------------- \n"
else
    echo -e "SAM to BAM - already done before ---------------------------- \n"

fi


#############
# lastChecks: create flagstats and idxstats files and remove redundant data
#############

if [ ! -e ${final_dir}/idxstats/ ]; then
    mkdir -p ${final_dir}/idxstats/
fi

echo -e "Starting last checks -------------------------------------- \n"

# check content of fifth line of step control file
linec=`sed "5q;d" ${stepControl}`
if [[ ${linec} != "lastChecks" ]]; then 

    echo -e "\nSAMTOOLS FLAGSTAT\n" >> ${summaryFile}
    samtools flagstat ${bamSortPath} | head -5 | tail -1 >> ${summaryFile}
    echo -e "\nSAMTOOLS FLAGSTAT OF UNMAPPED\n" >> ${summaryFile}
    samtools flagstat ${bamSortPath_unM} | head -5 | tail -1 >> ${summaryFile}
    echo -e "\n"  >> ${summaryFile}

    # idxstats = counts in this case
    samtools idxstats ${bamSortPath} > ${final_dir}/idxstats/${filename}.idxstats
    samtools idxstats ${bamSortPath_unM} > ${final_dir}/idxstats/${filename}.unMap.idxstats

    # Remove intermediate files, pickup the right ones... as needed
    rm -rf ${samPath}
    rm -rf ${samPath_unM}
    rm -rf ${bamPath}
    rm -rf ${bamPath_unM}
    rm -rf ${extr_fastq}

    # store stage control info
    echo "lastChecks" >> ${stepControl}

    echo -e "last checks - done -------------------------------------- \n"
else
    echo -e "last checks - already done before ---------------------------- \n"

fi

#############
# Rscript: prepare code to run R script
#############

if [ ! -e ${final_dir}/RSession/ ]; then
    mkdir -p ${final_dir}/RSession/
fi

echo -e "Starting Rscript text ---------------------------------- \n"

# check content of sixth line of step control file
linec=`sed "6q;d" ${stepControl}`
if [[ ${linec} != "Rscript" ]]; then 
    echo "module load R/4.0.5-foss-2020b" > ${final_dir}/RSession/toRunR.txt
    echo "python ${reportScript} -ip ${final_dir}" >> ${final_dir}/RSession/toRunR.txt
    echo "Rscript --vanilla ${RscriptP} ${final_dir} > \
${final_dir}/RSession/${sampleRunName}.Rout.txt" >> ${final_dir}/RSession/toRunR.txt
    echo "cd ${final_dir}/RSession" >> ${final_dir}/RSession/toRunR.txt
    echo "zip ${sampleRunName}.zip *" >> ${final_dir}/RSession/toRunR.txt

    # store stage control info
    echo "Rscript" >> ${stepControl}

    echo -e "Rscript - done -------------------------------------- \n"
else
    echo -e "Rscript - already done before ---------------------------- \n"
fi


echo -e "FINISHED... ------------------------------------------------------\n"

echo -e  Elapsed time $(($SECONDS - $START_TIME)) s

seff $SLURM_JOBID
