
cd $HOME/omicstmp/



# AGGREGATE DATA

echo "sample_id,molecule_h5" > INT_00_Aggr.csv
echo "ECCITE1,$DATA/ECCITE1_RNA_cellranger_601/outs/molecule_info.h5"  >> INT_00_Aggr.csv
echo "CITESEQ1,$DATA/CITESEQ1_RNAonly/outs/molecule_info.h5" >> INT_00_Aggr.csv
echo "CITESEQ2,$DATA/CITESEQ2/outs/molecule_info.h5" >> INT_00_Aggr.csv
echo "ECCITE2,$DATA/ECCITE2/outs/molecule_info.h5" >> INT_00_Aggr.csv
cat INT_00_Aggr.csv

$HOME/code/cellranger-6.0.1/cellranger aggr --id=INT_00_Aggr --csv=INT_00_Aggr.csv --normalize=none &> INT_00_Aggr.log

rm -r INT_00_Aggr/SC_RNA_AGGREGATOR_CS/
mv INT_00_Aggr* $DATA/



# # REMOVE UNWANTED CELLS
#
# name=INT_00_Aggr_CLEAN
#
# echo "Barcode" > ${name}_barcodes.csv
# cat $DATA/ECCITE1_RNA_cellranger_601/outs/filtered_feature_bc_matrix/barcodes.tsv >> ${name}_barcodes.csv
# cut -f1 -d, $ANALYSIS/CITESEQ1_01_InitialAnalysis/Cells_Keep.csv | grep -v "Barcode" | sed "s/\-1/-2/g" >> ${name}_barcodes.csv
#
# grep -n "Barcode" ${name}_barcodes.csv
# grep "\-1" ${name}_barcodes.csv | wc -l
# grep "\-2" ${name}_barcodes.csv | wc -l
# gunzip -c $DATA/INT_00_Aggr/outs/count/filtered_feature_bc_matrix/barcodes.tsv.gz | grep "\-1" | wc -l
# gunzip -c $DATA/INT_00_Aggr/outs/count/filtered_feature_bc_matrix/barcodes.tsv.gz | grep "\-2" | wc -l
#
# $CODEBASE/cellranger-6.0.1/cellranger reanalyze \
#  --id=$name \
#  --matrix=$DATA/INT_00_Aggr/outs/count/filtered_feature_bc_matrix.h5 \
#  --barcodes=${name}_barcodes.csv \
#  --localcores=24 \
#  --localmem=64 &> $name.log
#
# rm -r $name/SC_RNA_REANALYZER_CS/
# mv ${name}* $DATA/