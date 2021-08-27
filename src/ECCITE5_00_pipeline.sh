source $CODEBASE/tfcf/setup.sh



# If genome doesn't exist, then create it
if [ ! -f "$HOME/omicstmp/newGenomeExtended/fasta/genome.fa" ]; then
    bash $CODE/src/PREP_Genome.sh
fi




######### BASIC ANALYSIS


cd $HOME/omicstmp

id=ECCITE5

echo $id

echo "Files"
cat $CODEBASE/tfcf/metadata/${id}_Library.csv

echo "Features"
cat $CODEBASE/tfcf/metadata/ECCITE4_Features.csv

~/code/cellranger-6.0.1/cellranger count --id=$id \
 --no-bam \
 --libraries=$CODEBASE/tfcf/metadata/${id}_Library.csv \
 --transcriptome=newGenomeExtended/ \
 --feature-ref=$CODEBASE/tfcf/metadata/ECCITE4_Features.csv \
 --localcores=24 \
 --localmem=128 \
 --expect-cells=10000 &> ${id}.log

mkdir -p ~/GFS/PROJECTS/TfCf/Data/$id
mv $id/outs ~/GFS/PROJECTS/TfCf/Data/$id
