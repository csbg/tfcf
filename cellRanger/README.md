# General
All scripts are numerated in order of execution, and are built to be run in a SLURM queueing system.

# FACS based CRISPR data analysis
## FILE NAMING FORMAT
There is no specific file naming format, but each experiments ID comes from the base name of the libraries.csv file declaring FASTQ paths and library types of input libraries in cellRanger. E.g: ${id}_Library.csv

## Scripts
    * 00_cellRcount_new.sh should be run first and contains all the steps to go from Fastq files to cellRanger's output. Input variables are defined and explained in the script.
    This scripts calls 00_NR_PREP_Genome.sh to generate a merged reference genome with GFP and BFP sequences.
