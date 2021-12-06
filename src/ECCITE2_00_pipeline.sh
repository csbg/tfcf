source $CODEBASE/tfcf/setup.sh

basedir=$HOME/GFS/PROJECTS/TfCf/


# # PREPARE DATA
# cd Data/Raw_data_ECCITE1
# python $CODEBASE/tfcf/crukci_to_illumina.py
#
#
# # ANALYZE WITH ANTIBODIES
# mkdir -p $HOME/omicstmp/nf
# cd $HOME/omicstmp/nf
#
# ~/code/cellranger-6.0.1/cellranger count --id=ECCITE2 \
# 	--no-bam \
#     --libraries=$CODEBASE/tfcf/metadata/ECCITE1_Library.csv \
#     --transcriptome=$HOME/GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/ \
#     --feature-ref=$CODEBASE/tfcf/metadata/ECCITE1_Features.csv \
#     --localcores=24 \
#     --localmem=64 \
#     --expect-cells=10000 &> ECCITE1.log
#
# mkdir -p $DATA/$id/
# mv $id/outs $DATA/$id/
#
#
#
# # ANALYZE WITHOUT GUIDES
# cd $HOME/omicstmp
#
# grep -v "CRISPR Guide Capture" $CODEBASE/tfcf/metadata/ECCITE2_Library.csv > $CODEBASE/tfcf/metadata/ECCITE2_Library_onlyRNA.csv
#
# id="ECCITE2_onlyRNA"
#
# ~/code/cellranger-6.0.1/cellranger count --id=$id \
# 	--no-bam \
#     --libraries=$CODEBASE/tfcf/metadata/ECCITE2_Library_onlyRNA.csv \
#     --transcriptome=$GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/ \
#     --localcores=24 \
#     --localmem=64 \
#     --expect-cells=10000 &> $id.log
#
# mkdir -p $DATA/$id/
# mv $id/outs $DATA/$id/


# ANALYZE WITHOUT GUIDES
cd $HOME/omicstmp

# If genome doesn't exist, then create it
if [ ! -f "$HOME/omicstmp/newGenomeExtended/fasta/genome.fa" ]; then
    bash $CODE/src/PREP_Genome.sh
fi

id="ECCITE2_onlyRNA_GFPBFP"

~/code/cellranger-6.0.1/cellranger count --id=$id \
	--no-bam \
    --libraries=$CODEBASE/tfcf/metadata/ECCITE2_Library_onlyRNA.csv \
    --transcriptome=newGenomeExtended/ \
    --localcores=24 \
    --localmem=64 \
    --expect-cells=10000 &> $id.log

mkdir -p $DATA/$id/
mv $id/outs $DATA/$id/
