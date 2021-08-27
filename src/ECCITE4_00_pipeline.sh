source $CODEBASE/tfcf/setup.sh


# If genome doesn't exist, then create it
if [ ! -f "$HOME/omicstmp/newGenomeExtended/fasta/genome.fa" ]; then
    bash $CODE/src/PREP_Genome.sh
fi



######### BASIC ANALYSIS


cd $HOME/omicstmp

for id in ECCITE4_Cas9 ECCITE4_WT; do
    
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
done



######### AGGREGATE DATA
source $CODEBASE/tfcf/setup.sh
cd $HOME/omicstmp/

id=ECCITE4_INT

echo "sample_id,molecule_h5" > $id.csv
echo "Cas9,$DATA/ECCITE4_Cas9/outs/molecule_info.h5"  >> $id.csv
echo "WT,$DATA/ECCITE4_WT/outs/molecule_info.h5" >> $id.csv
cat $id.csv

$HOME/code/cellranger-6.0.1/cellranger aggr --id=$id --csv=$id.csv --normalize=none &> $id.log

mkdir -p ~/GFS/PROJECTS/TfCf/Data/$id/
mv $id/outs ~/GFS/PROJECTS/TfCf/Data/$id/