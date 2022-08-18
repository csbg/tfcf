source $CODEBASE/tfcf/setup.sh
to=$HOME/Desktop/share/
mkdir -p $to


cd $ANALYSIS
zip -r --exclude=*.RData $to/invivo.zip INVIVO_01_01_Integration/*
zip -r --exclude=*.RData $to/invitro.zip ECCITE4_04_Integration/*
zip -r --exclude=*.RData $to/leukemia.zip INT_03_02_SeuratIntegration_CellCycle_Leukemia/*

ds="leukemia"
cp INT_03_02_SeuratIntegration_CellCycle_Leukemia/GeneCount* $to
mv $to/GeneCount.pdf $to/GeneCount_$ds.pdf
mv $to/GeneCount.tsv $to/GeneCount_$ds.tsv

ds="invivo"
cp INVIVO_01_01_Integration/GeneCount* $to
mv $to/GeneCount.pdf $to/GeneCount_$ds.pdf
mv $to/GeneCount.tsv $to/GeneCount_$ds.tsv

ds="invitro"
cp ECCITE4_04_Integration/GeneCount* $to
mv $to/GeneCount.pdf $to/GeneCount_$ds.pdf
mv $to/GeneCount.tsv $to/GeneCount_$ds.tsv


cd $DATA
cp INT_00_Aggr/outs/count/cloupe.cloupe $to/leukemia.cloupe
cp ECCITE4_INT/outs/count/cloupe.cloupe $to/invitro.cloupe
cp INVIVO_00_Aggr/outs/count/cloupe.cloupe $to/invivo.cloupe


