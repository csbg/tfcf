# TO BE RUN ON THE CAME SERVER

# To install cellranger
# wget --no-check-certificate -O cellranger-6.0.1.tar.gz 

basedir=$HOME/GFS/PROJECTS/TfCf/

cd $HOME/GFS/DATA_David/Raw_data_ECCITE/gRNA/
ls  $HOME/code/tfcf/crukci_to_illumina.py # END OF THIS NOTE
mv SLX-19954.HY75WDRXX.s_1.r_1.lostreads.fq.gz SITTH12_S1_L001_R1_001.fastq.gz
mv SLX-19954.HY75WDRXX.s_1.r_2.lostreads.fq.gz SITTH12_S1_L001_R2_001.fastq.gz
mv SLX-19954.HY75WDRXX.s_1.i_1.lostreads.fq.gz SITTH12_S1_L001_I1_001.fastq.gz
mv SLX-19954.HY75WDRXX.s_1.i_2.lostreads.fq.gz SITTH12_S1_L001_I2_001.fastq.gz

# CITESEQ - all barcodes
#conda create -n tfcf
conda activate tfcf
#conda install pip
#pip install CITE-seq-Count==1.4.4
cd ~/GFS/DATA_David/
CITE-seq-Count -R1 Raw_data_ECCITE/gRNA/SITTH12_S1_L001_R1_001.fastq.gz -R2 Raw_data_ECCITE/gRNA/SITTH12_S1_L001_R2_001.fastq.gz -trim 18 -t TAG_LIST.csv -cbf 1 -cbl 16 -umif 17 -umil 26 --expected_cells 10000 -o citeseq_all_barcodes --whitelist /home/people/nfortelny/code/cellranger-5.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt &> citeseq_all_barcodes.log
conda deactivate

# CELLRANGER
cd $HOME/GFS/DATA_David/Raw_data_ECCITE/
nano  crukci_to_illumina.py # END OF THIS NOTE
cp -R RNA/ RNA_illumina/
python3 crukci_to_illumina.py RNA_illumina/
~/code/cellranger-5.0.1/cellranger count --id=RNA_cellranger --transcriptome=$HOME/GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/ --no-bam --expect-cells=10000 --localcores=1 --fastqs=$HOME/GFS/DATA_David/Raw_data_ECCITE/RNA_illumina/ --localmem=64 &> RNA_cellranger.log


# CITESEQ - barcodes from cellranger
conda activate tfcf
zcat RNA_cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | sed 's/\-[0-9]*$//g' > RNA_cellranger_barcodes.txt
cd ~/GFS/DATA_David/
CITE-seq-Count -R1 Raw_data_ECCITE/gRNA/SITTH12_S1_L001_R1_001.fastq.gz -R2 Raw_data_ECCITE/gRNA/SITTH12_S1_L001_R2_001.fastq.gz -trim 18 -t Raw_data_ECCITE/TAG_LIST.csv -cbf 1 -cbl 16 -umif 17 -umil 26 --expected_cells 10000 -o citeseq_filtered_barcodes --whitelist RNA_cellranger_barcodes.txt &> citeseq_filtered_barcodes.log
conda deactivate

# CITESEQ - combined - barcodes from cellranger
conda activate tfcf
cd $basedir
CITE-seq-Count -R1 Raw_data_ECCITE/gRNA/SITTH12_S1_L001_R1_001.fastq.gz,Raw_data_ECCITE/gRNA/SLX-20379.HYHJ2DRXX.s_2.r_1.lostreads.fq.gz -R2 Raw_data_ECCITE/gRNA/SITTH12_S1_L001_R2_001.fastq.gz,Raw_data_ECCITE/gRNA/SLX-20379.HYHJ2DRXX.s_2.r_2.lostreads.fq.gz -trim 18 -t Raw_data_ECCITE/TAG_LIST.csv -cbf 1 -cbl 16 -umif 17 -umil 26 --expected_cells 10000 -o citeseq_combined --whitelist /home/people/nfortelny/code/cellranger-5.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt &> citeseq_combined.log
conda deactivate

# CELLRANGER
cd $HOME/omicstmp/
~/code/cellranger-5.0.1/cellranger count --id=gRNA_cellranger --transcriptome=$HOME/GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/ --no-bam --expect-cells=10000 --localcores=1 --fastqs=$HOME/GFS/DATA_David/Raw_data_ECCITE/gRNA/ --localmem=64 &> gRNA_cellranger.log

