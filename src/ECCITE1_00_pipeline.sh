# TO BE RUN ON THE CAME SERVER

source $CODEBASE/tfcf/setup.sh

# ### DATA:
#
# # The Original data was SLX-19954.HY75WDRXX
# # - Pattern of files with guide RNA reads: SLX-19954.HY75WDRXX.s_?.r_?.lostreads.fq.gz
# # --> These were moved to the folder gRNA/
# # - Pattern of files with mRNA reads (transcriptome): SLX-19954.HY75WDRXX.s_?.r_?.fq.gz
# # --> These were moved to the folder RNA/ and  then renamed and moved to RNA_illumina/
#
# # The guide RNAs were resequenced in SLX-20379.HYHJ2DRXX
# # - Pattern of files with guide RNA reads: SLX-20379.HYHJ2DRXX.s_?.r_?.fq.gz
# # - In this case, the transcriptome was not sequenced
# # --> These were moved to the folder gRNA/
#
# basedir=$DATA
#
# # Here we rename the lost reads to illumina names so that we can run them through cellranger (this was not used later as CITE-seq worked)
# cd $DATA/Raw_ECCITE1/gRNA/
# ls $CODE/crukci_to_illumina.py # END OF THIS NOTE - not used here but done manually
# mv SLX-19954.HY75WDRXX.s_1.r_1.lostreads.fq.gz SITTH12_S1_L001_R1_001.fastq.gz
# mv SLX-19954.HY75WDRXX.s_1.r_2.lostreads.fq.gz SITTH12_S1_L001_R2_001.fastq.gz
# mv SLX-19954.HY75WDRXX.s_1.i_1.lostreads.fq.gz SITTH12_S1_L001_I1_001.fastq.gz
# mv SLX-19954.HY75WDRXX.s_1.i_2.lostreads.fq.gz SITTH12_S1_L001_I2_001.fastq.gz
#
# # CITESEQ - first run using all cell barcodes (obtained from the cellranger code)
# # conda activate tfcf
# # cd ~/GFS/DATA_David/
# # CITE-seq-Count -R1 Raw_data_ECCITE/gRNA/SITTH12_S1_L001_R1_001.fastq.gz -R2 Raw_data_ECCITE/gRNA/SITTH12_S1_L001_R2_001.fastq.gz -trim 18 -t TAG_LIST.csv -cbf 1 -cbl 16 -umif 17 -umil 26 --expected_cells 10000 -o citeseq_all_barcodes --whitelist /home/people/nfortelny/code/cellranger-5.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt &> citeseq_all_barcodes.log
# # conda deactivate
#
# # CELLRANGER - gRNA - Not used
# # cd $HOME/omicstmp/
# # ~/code/cellranger-5.0.1/cellranger count --id=gRNA_cellranger --transcriptome=$HOME/GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/ --no-bam --expect-cells=10000 --localcores=1 --fastqs=$HOME/GFS/DATA_David/Raw_ECCITE1/gRNA/ --localmem=64 &> gRNA_cellranger.log
#
# # CELLRANGER - mRNA (transcriptome)
# cd $DATA/Raw_ECCITE1/
# cp -R RNA/ RNA_illumina/
# cd RNA_illumina/
# # Convert names to cellranger compatible names
# python3 $CODE/crukci_to_illumina.py
# cd ~/omicstmp/nf/
# # Cellranger analysis (5.0.1) - not used
# # ~/code/cellranger-5.0.1/cellranger count --id=RNA_cellranger --transcriptome=$HOME/GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/ --no-bam --expect-cells=10000 --localcores=1 --fastqs=$HOME/GFS/DATA_David/Raw_data_ECCITE/RNA_illumina/ --localmem=64 &> RNA_cellranger.log
# # Cellranger analysis (6.0.1) - USED
# $CODEBASE/cellranger-6.0.1/cellranger count --id=ECCITE1 --transcriptome=$GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/ --no-bam --expect-cells=10000 --localcores=1 --fastqs=$DATA/Raw_ECCITE1/RNA_illumina/ --localmem=64 &> ECCITE1.log

#
# # CITESEQ - second run using barcodes from the dataset processed with cellranger (not used)
# # conda activate tfcf
# # zcat RNA_cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | sed 's/\-[0-9]*$//g' > RNA_cellranger_barcodes.txt
# # cd ~/GFS/DATA_David/
# # CITE-seq-Count -R1 Raw_data_ECCITE/gRNA/SITTH12_S1_L001_R1_001.fastq.gz -R2 Raw_data_ECCITE/gRNA/SITTH12_S1_L001_R2_001.fastq.gz -trim 18 -t Raw_ECCITE1/TAG_LIST.csv -cbf 1 -cbl 16 -umif 17 -umil 26 --expected_cells 10000 -o citeseq_filtered_barcodes --whitelist RNA_cellranger_barcodes.txt &> citeseq_filtered_barcodes.log
# # conda deactivate
#
# # CITESEQ - using barcodes from the cellranger code, analyzing both sets of data together (SLX-19954.HY75WDRXX and SLX-20379.HYHJ2DRXX)
# conda activate tfcf
# cd $basedir
# CITE-seq-Count \
#  -R1 Raw_ECCITE1/gRNA/SITTH12_S1_L001_R1_001.fastq.gz,Raw_ECCITE1/gRNA/SLX-20379.HYHJ2DRXX.s_2.r_1.lostreads.fq.gz \
#  -R2 Raw_ECCITE1/gRNA/SITTH12_S1_L001_R2_001.fastq.gz,Raw_ECCITE1/gRNA/SLX-20379.HYHJ2DRXX.s_2.r_2.lostreads.fq.gz \
#  -trim 18 -t Raw_ECCITE1/TAG_LIST.csv -cbf 1 -cbl 16 -umif 17 -umil 26 --expected_cells 10000 -o citeseq_combined \
#  --whitelist $CODEBASE/cellranger-6.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt &> citeseq_combined.log
# conda deactivate




# ANALYZE WITHOUT GUIDES
cd $HOME/omicstmp

# If genome doesn't exist, then create it
if [ ! -f "$HOME/omicstmp/newGenomeExtended/fasta/genome.fa" ]; then
    bash $CODE/src/PREP_Genome.sh
fi

id=ECCITE1_RNA_cellranger_601_GFPBFP
$CODEBASE/cellranger-6.0.1/cellranger count \
	--id=$id \
	--transcriptome=newGenomeExtended/ \
	--no-bam \
	--expect-cells=10000 \
	--localcores=24 \
	--fastqs=$RAWDATA/ECCITE1/RNA_illumina/ \
	--localmem=64 &> $id.log

mkdir -p ~/GFS/PROJECTS/TfCf/Data/$id/
mv $id/outs ~/GFS/PROJECTS/TfCf/Data/$id/
