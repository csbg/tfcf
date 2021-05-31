cd $GFS/PROJECTS/TfCf/Data/Raw_ECCITE2

python $CODEBASE/tfcf/crukci_to_illumina.py


# ANALYZE WITH ANTIBODIES
mkdir -p $HOME/omicstmp/nf
cd $HOME/omicstmp/nf

~/code/cellranger-6.0.1/cellranger count --id=ECCITE2 \
 --no-bam \
 --libraries=$CODEBASE/tfcf/metadata/ECCITE2_Library.csv \
 --transcriptome=$HOME/GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/ \
 --feature-ref=$CODEBASE/tfcf/metadata/ECCITE2_Features.csv \
 --localcores=24 \
 --localmem=64 \
 --expect-cells=10000 &> ECCITE2.log

mkdir -p ~/GFS/PROJECTS/TfCf/Data/ECCITE2/
cp -R ECCITE2/outs ~/GFS/PROJECTS/TfCf/Data/ECCITE2/