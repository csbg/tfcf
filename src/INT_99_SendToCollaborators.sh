cd $ANALYSIS

to=$HOME/Desktop/share/

mkdir -p $to

zip -r --exclude=*.RData $to/invivo.zip INVIVO_01_01_Integration/*
zip -r --exclude=*.RData $to/invitro.zip ECCITE4_04_Integration/*
zip -r --exclude=*.RData $to/leukemia.zip INT_03_02_SeuratIntegration_CellCycle_Leukemia/*

cd $DATA

cp INT_00_Aggr/outs/count/cloupe.cloupe $to/leukemia.cloupe
cp ECCITE4_INT/outs/count/cloupe.cloupe $to/invitro.cloupe
