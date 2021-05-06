

cd $HOME/omicstmp/nf/

echo "sample_id,molecule_h5" > INT_00_Aggr.csv
echo "ECCITE1,$DATA/ECCITE1_RNA_cellranger_601/outs/molecule_info.h5"  >> INT_00_Aggr.csv
echo "CITESEQ1,$DATA/CITESEQ1_RNAonly/outs/molecule_info.h5" >> INT_00_Aggr.csv

cat INT_00_Aggr.csv

$HOME/code/cellranger-6.0.1/cellranger aggr --id=INT_00_Aggr --csv=INT_00_Aggr.csv --normalize=none &> INT_00_Aggr.log