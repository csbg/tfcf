source $CODEBASE/tfcf/setup.sh

cd $HOME/omicstmp

for id in ECCITE3_low_7d ECCITE3_low_14d ECCITE3_high_7d ECCITE3_high_14d; do
        
    echo $id
    #id="ECCITE3_low_7d"
    #id="ECCITE3_low_14d"
    
    echo "Files"
    cat $CODEBASE/tfcf/metadata/${id}_Library.csv
    
    echo "Features"
    cat $CODEBASE/tfcf/metadata/ECCITE3_Features.csv
    
	~/code/cellranger-6.0.1/cellranger count --id=$id \
	 --no-bam \
	 --libraries=$CODEBASE/tfcf/metadata/${id}_Library_CRISPR.csv \
	 --transcriptome=$HOME/GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/ \
	 --feature-ref=$CODEBASE/tfcf/metadata/ECCITE3_Features.csv \
	 --localcores=24 \
	 --localmem=128 \
	 --expect-cells=10000 &> ${id}.log
    
    rm -rf $id/SC_RNA_COUNTER_CS
    
    mv ${id}* $GFS/PROJECTS/TfCf/Data/
done