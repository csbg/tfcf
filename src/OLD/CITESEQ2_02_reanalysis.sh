cd $HOME/omicstmp/

id=CITESEQ2_CLEAN

~/code/cellranger-6.0.1/cellranger reanalyze \
	--id=$id \
    --matrix=$GFS/PROJECTS/TfCf/Data/CITESEQ2/outs/filtered_feature_bc_matrix.h5 \
	--barcodes=$GFS/PROJECTS/TfCf/Analysis/CITESEQ2_01_InitialAnalysis/Cells_Keep.csv \
	--localcores=24 \
    --localmem=64 &> ${id}.log

mkdir -p ~/GFS/PROJECTS/TfCf/Data/$id/
mv $id/outs ~/GFS/PROJECTS/TfCf/Data/$id/
#rm -rf ${id}*