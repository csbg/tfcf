

cd $HOME/omicstmp/

echo "sample_id\tmolecule_h5" > AGGR1.aggr.csv
echo "ECCITE1\t$DATA/ECCITE1_RNA_cellranger/outs/molecule_info.h5"  >> AGGR1.aggr.csv
echo "CITESEQ1\t$DATA/CITESEQ1_CLEAN/outs/molecule_info.h5" >> AGGR1.aggr.csv

cd $HOME/omicstmp/
~/code/cellranger-6.0.1/cellranger aggr --id=AGGR1 --csv=AGGR1.aggr.csv --normalize=none
