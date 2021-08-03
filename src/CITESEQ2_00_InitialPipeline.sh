source ~/code/tfcf/setup.sh

cd $GFS/PROJECTS/TfCf/NewData/
python ~/code/tfcf/crukci_to_illumina.py


### CREATE GENOME

cd $HOME/omicstmp

mkdir newGenome/
oldGenome="$GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/"
echo "ORIGINAL GENOME: $oldGenome"

cp $oldGenome/fasta/genome.fa newGenome/genome.fa
cp $oldGenome/genes/genes.gtf newGenome/genes.gtf

x="GFP"
for x in GFP BFP; do
	pathFA=$CODEBASE/tfcf/metadata/${x}.fa

	numberBases=$(cat $pathFA | grep -v "^>" | tr -d "\n" | wc -c)
	gtf="$x\tunknown\texon\t1\t${numberBases}\t.\t+\t.\tgene_id xxx${x}xxx; transcript_id xxx${x}xxx; gene_name xxx${x}xxx; gene_biotype xxxprotein_codingxxx;"
	gtf=$(echo $gtf | sed 's/xxx/"/g')

	echo -e $gtf > ${x}.gtf
    l 
	cat ${x}.gtf >> newGenome/genes.gtf
	fold -w 60 $pathFA >> newGenome/genome.fa
	
	echo -e "" >> newGenome/genome.fa

	rm ${x}.gtf
done


grep ">" newGenome/ls -algenome.fa
tail -5 newGenome/genes.gtf
tail -200 newGenome/genome.fa


~/code/cellranger-6.0.1/cellranger mkref --genome=newGenomeExtended --fasta=newGenome/genome.fa --genes=newGenome/genes.gtf




# ANALYZE WITH ANTIBODIES
cd $HOME/omicstmp/

~/code/cellranger-6.0.1/cellranger count --id=CITESEQ2 \
	--no-bam \
    --libraries=$CODEBASE/tfcf/metadata/CITESEQ2_Library.csv \
    --transcriptome=newGenomeExtended/ \
    --feature-ref=$CODEBASE/tfcf/metadata/CITESEQ2_Features.csv \
    --localcores=24 \
    --localmem=64 \
    --expect-cells=10000 &> CITESEQ2.log

mkdir -p ~/GFS/PROJECTS/TfCf/Data/CITESEQ2/
mv CITESEQ2/outs ~/GFS/PROJECTS/TfCf/Data/CITESEQ2/
rm -rf CITESEQ2*
rm -rf newGenome**


# ANALYZE WITH ANTIBODIES but normal Genome (for later integration with cellranger aggr)
id=CITESEQ2_standardGenome

cd $HOME/omicstmp/

~/code/cellranger-6.0.1/cellranger count --id=$id \
	--no-bam \
    --libraries=$CODEBASE/tfcf/metadata/CITESEQ2_Library.csv \
    --transcriptome=$GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/ \
    --feature-ref=$CODEBASE/tfcf/metadata/CITESEQ2_Features.csv \
    --localcores=24 \
    --localmem=64 \
    --expect-cells=10000 &> $id.log

mkdir -p ~/GFS/PROJECTS/TfCf/Data/$id/
mv $id/outs ~/GFS/PROJECTS/TfCf/Data/$id/
rm -rf ${id}*
