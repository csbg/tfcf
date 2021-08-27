source $CODEBASE/tfcf/setup.sh


# If genome doesn't exist, then create it
if [ ! -f "$HOME/omicstmp/newGenomeExtended/fasta/genome.fa" ]; then
    bash $CODE/src/PREP_Genome.sh
fi



# ######### BASIC ANALYSIS
#
# cd $HOME/omicstmp
#
# id=ECCITE6
#
# echo $id
#
# echo "Files"
# cat $CODEBASE/tfcf/metadata/${id}_Library.csv
#
# echo "Features"
# cat $CODEBASE/tfcf/metadata/ECCITE6_Features.csv
#
# ~/code/cellranger-6.0.1/cellranger count --id=$id \
#  --no-bam \
#  --libraries=$CODEBASE/tfcf/metadata/${id}_Library.csv \
#  --transcriptome=newGenomeExtended/ \
#  --feature-ref=$CODEBASE/tfcf/metadata/ECCITE6_Features.csv \
#  --localcores=24 \
#  --localmem=128 \
#  --expect-cells=10000 &> ${id}.log
#
# mkdir -p ~/GFS/PROJECTS/TfCf/Data/$id
# cp -R $id/outs ~/GFS/PROJECTS/TfCf/Data/$id


######### BASIC ANALYSIS without guides

cd $HOME/omicstmp

id_original=ECCITE6

echo $id_original
id="${id_original}_onlyRNA"
echo $id

grep -v "CRISPR Guide Capture" $CODE/metadata/${id_original}_Library.csv > $CODE/metadata/${id}_Library.csv

echo "Files"
cat $CODEBASE/tfcf/metadata/${id}_Library.csv


~/code/cellranger-6.0.1/cellranger count --id=$id \
 --no-bam \
 --libraries=$CODE/metadata/${id}_Library.csv \
 --transcriptome=newGenomeExtended/ \
 --localcores=4 \
 --localmem=12 \
 --expect-cells=10000 &> ${id}.log

mkdir -p $DATA/$id
mv $id/outs $DATA/$id