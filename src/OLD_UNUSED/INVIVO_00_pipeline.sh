source $CODEBASE/tfcf/setup.sh


cd $HOME/omicstmp/

id=INVIVO_00_Aggr

# AGGREGATE DATA

echo "sample_id,molecule_h5" > $id.csv
echo "ECCITE3_high_7d,$DATA/ECCITE3_high_7d_onlyRNA/outs/molecule_info.h5" >> $id.csv
echo "ECCITE3_high_14d,$DATA/ECCITE3_high_14d_onlyRNA/outs/molecule_info.h5" >> $id.csv
echo "ECCITE3_low_7d,$DATA/ECCITE3_low_7d_onlyRNA/outs/molecule_info.h5"  >> $id.csv
echo "ECCITE3_low_14d,$DATA/ECCITE3_low_14d_onlyRNA/outs/molecule_info.h5" >> $id.csv
echo "ECCITE6,$DATA/ECCITE6_onlyRNA/outs/molecule_info.h5" >> $id.csv

cat $id.csv

$HOME/code/cellranger-6.0.1/cellranger aggr --id=$id --csv=$id.csv --normalize=none &> $id.log

mkdir -p $DATA/$id/
mv $id/outs $DATA/$id/

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