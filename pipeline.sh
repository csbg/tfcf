# TO BE RUN ON THE CAME SERVER

# CITESEQ
#conda create -n tfcf
conda activate tfcf
#conda install pip
#pip install CITE-seq-Count==1.4.4
cd ~/GFS/DATA_David/Raw_data_ECCITE/
CITE-seq-Count -R1 gRNA/SLX-19954.HY75WDRXX.s_1.r_1.lostreads.fq.gz -R2 gRNA/SLX-19954.HY75WDRXX.s_1.r_2.lostreads.fq.gz -trim 18 -t TAG_LIST.csv -cbf 1 -cbl 16 -umif 17 -umil 26 --expected_cells 10000 -o citeseq --whitelist /home/people/nfortelny/code/cellran$

# CELLRANGER
cd $HOME/GFS/DATA_David/Raw_data_ECCITE/
nano  crukci_to_illumina.py # END OF THIS NOTE
cp -R RNA/ RNA_illumina/
python3 crukci_to_illumina.py RNA_illumina/
~/code/cellranger-5.0.1/cellranger count --id=RNA_cellranger --transcriptome=$HOME/GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/ --no-bam --expect-cells=10000 --localcores=1 --fastqs=$HOME/GFS/DATA_David/Raw_data_ECCITE/RNA_illumina/ --localmem=64 &> RNA_cell$






