cd $HOME/omicstmp/
mkdir nf
cd nf

~/code/cellranger-6.0.1/cellranger reanalyze \
	--id=CITESEQ1_CLEAN \
    --matrix=$GFS/PROJECTS/TfCf/Data/CITESEQ1/outs/filtered_feature_bc_matrix.h5 \
	--barcodes=$GFS/PROJECTS/TfCf/Analysis/CITESEQ1_01_InitialAnalysis/Cells_Keep.csv \
	--localcores=24 \
    --localmem=64 &> CITESEQ1_CLEAN.log

cp -R CITESEQ1_CLEAN* ~/GFS/PROJECTS/TfCf/Data/
