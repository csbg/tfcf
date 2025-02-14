source $CODEBASE/tfcf/setup.sh



# If genome doesn't exist, then create it
if [ ! -f "$HOME/omicstmp/newGenomeExtended/fasta/genome.fa" ]; then
    bash $CODE/src/PREP_Genome.sh
fi




# ######### BASIC ANALYSIS
#
#
# cd $HOME/omicstmp
#
# id=ECCITE5
#
# echo $id
#
# echo "Files"
# cat $CODEBASE/tfcf/metadata/${id}_Library.csv
#
# echo "Features"
# cat $CODEBASE/tfcf/metadata/ECCITE4_Features.csv
#
# ~/code/cellranger-6.0.1/cellranger count --id=$id \
#  --no-bam \
#  --libraries=$CODEBASE/tfcf/metadata/${id}_Library.csv \
#  --transcriptome=newGenomeExtended/ \
#  --feature-ref=$CODEBASE/tfcf/metadata/ECCITE4_Features.csv \
#  --localcores=24 \
#  --localmem=128 \
#  --expect-cells=10000 &> ${id}.log
#
# mkdir -p $DATA/$id
# mv $id/outs $DATA/$id



######### WITHOUT GUIDES

cd $HOME/omicstmp

id=ECCITE5_onlyRNA

grep -v "CRISPR Guide Capture" $CODE/metadata/ECCITE5_Library.csv > $CODE/metadata/${id}_Library.csv

echo $id

echo "Files"
cat $CODE/metadata/${id}_Library.csv


~/code/cellranger-6.0.1/cellranger count --id=$id \
 --no-bam \
 --libraries=$CODE/metadata/${id}_Library.csv \
 --transcriptome=newGenomeExtended/ \
 --localcores=8 \
 --localmem=24 \
 --expect-cells=10000 &> ${id}.log

mkdir -p $DATA/$id
mv $id/outs $DATA/$id
