# General
All scripts are numerated in order of execution, and are built to be run in a SLURM queueing system.

# FACS based CRISPR data analysis
## FILE NAMING FORMAT
[Cell]_[extra]_[libraryID]_[date?]_...._S[0-9]_R[12]_001.fastq.gz
Only the first and third sections between '_' will be used for file handling

## Scripts
    * 00_MECC_sgRNA_v1.1.sh should be run first and contains all the steps to go from Fastq files to filtered BAM files. Specifically, the steps to extract guide RNAs (gRNA) from the Fastq files, Align them to gRNA library references, generate indexed bam files, flagstat and idxstats of them, and output a script to run the remaining scripts in the folder. The script expects four inputs:
        *  The path to the Fastq files. It is expected to contain a demux_fastq folder with the fastq files and the samplesNames.txt file (instructions in the script). Temporal analysis folders will be created in it during the run
        * The output path
        * The base name (with full path) to the folder containing the reference genomes for each of the libraries.
        * The pattern to extract the guides from the Fastq files. E.g.: 'CACCG(.{20})GT'
    This scripts generates a file at /OUTPUT/FOLDER/RSession/toRunR.txt containing the commands needed to run the remaining scripts.
    * 01a_MECC_sgRNA_pre_analysis_v2.r contains the steps to go make count tables of the gRNAs per sample. It has been tested in R v4.1.0. The script expects one input:
        * The output path indicated in 00_MECC_sgRNA_v1.1.sh
    *  01b_mergedMappingReport.py generates an alignment report. The script expects one input:
        * The output path indicated in 00_MECC_sgRNA_v1.1.sh
    * Scripts containing NR are called by the previous scripts
