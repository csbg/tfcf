

cd $HOME/omicstmp/nf/

echo "sample_id,molecule_h5" > AGGR1.aggr.csv
echo "ECCITE1,$DATA/ECCITE1_RNA_cellranger/outs/molecule_info.h5"  >> AGGR1.aggr.csv
echo "CITESEQ1,$DATA/CITESEQ1/outs/molecule_info.h5" >> AGGR1.aggr.csv

cat AGGR1.aggr.csv

$HOME/code/cellranger-6.0.1/cellranger aggr --id=AGGR1 --csv=AGGR1.aggr.csv --normalize=none &> AGGR1.log
