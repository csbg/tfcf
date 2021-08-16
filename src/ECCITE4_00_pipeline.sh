source $CODEBASE/tfcf/setup.sh



### CREATE GENOME

# cd $HOME/omicstmp
#
# mkdir newGenome/
# oldGenomePath="$GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/"
#
# cp $oldGenomePath/fasta/genome.fa newGenome/genome.fa
# cp $oldGenomePath/genes/genes.gtf newGenome/genes.gtf
#
# x="GFP"
# for x in GFP BFP; do
#     pathFA=$CODEBASE/tfcf/metadata/${x}.fa
#
#     numberBases=$(cat $pathFA | grep -v "^>" | tr -d "\n" | wc -c)
#     gtf="$x\tunknown\texon\t1\t${numberBases}\t.\t+\t.\tgene_id xxx${x}xxx; transcript_id xxx${x}xxx; gene_name xxx${x}xxx; gene_biotype xxxprotein_codingxxx;"
#     gtf=$(echo $gtf | sed 's/xxx/"/g')
#
#     echo -e $gtf > ${x}.gtf
#
#     cat ${x}.gtf >> newGenome/genes.gtf
#     fold -w 60 $pathFA >> newGenome/genome.fa
#
#     echo -e "" >> newGenome/genome.fa
#
#     rm ${x}.gtf
# done
#
#
# grep ">" newGenome/genome.fa
# tail -5 newGenome/genes.gtf
# tail -200 newGenome/genome.fa
#
# ~/code/cellranger-6.0.1/cellranger mkref --genome=newGenomeExtended --fasta=newGenome/genome.fa --genes=newGenome/genes.gtf



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
    
    rm -rf $id/SC_RNA_COUNTER_CS
    
    mv ${id}* $GFS/PROJECTS/TfCf/Data/
done