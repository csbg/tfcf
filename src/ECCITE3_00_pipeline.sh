source $CODEBASE/tfcf/setup.sh

# # # Figuring out where the guides are
# cd $RAWDATA/Raw_ECCITE3/demux_fastq/
#
# sequences=( TGTGGCGATAGAGCTGCTGT CGGGGAACTTTCACCCATCA CTCGTTCCCTAACGGCGCGG ACAGCAGCTCTATCGCCACA TGATGGGTGAAAGTTCCCCG CCGCGCCGTTAGGGAACGAG)
# files=($(ls))
#
# for f in "${files[@]}"; do
#   echo $f
#   for s in "${sequences[@]}"; do
#     n=$(gunzip -c $f | head -400000 | grep $s | wc -l)
#     echo $s $n
#   done
# done
#
# gunzip -c T_high_sgRNA_d7_24052021_S4_R2_001.fastq.gz | grep "TGATGGGTGAAAGTTCCCCG"
# gunzip -c T_high_sgRNA_d7_24052021_S4_R2_001.fastq.gz | grep "TGATGGGTGAAAGTTCCCCG" | head -20
# gunzip -c T_low_sgRNA_d14_31052021_S7_R2_001.fastq.gz | grep "CCGCGCCGTTAGGGAACGAG" | head -20
#
# # ECCITE2:
# gunzip -c SINTF10_S1_L001_R2_001.fastq.gz | grep "GACTCCGGGTACTAAATGTC"

# # COUNT SEQUENCES
# files=($(ls *sgRNA*R2_001.fastq.gz))
# for f in "${files[@]}"; do
#   echo $f
#   gunzip -c $f | head -4000 | grep "AAGCAGTGGTATCAACGCAGAGTACAT" | wc -l >   ~/$f.spacer.counts.txt
#   gunzip -c $f | grep "AAGCAGTGGTATCAACGCAGAGTACAT" | sed "s/AAGCAGTGGTATCAACGCAGAGTACAT//g" | head -10000 | sort | uniq -c | sort -k1 | tail -50 > ~/$f.guide.counts.txt
# done

# # COPY CLOUPE
# cd ~/GFS/PROJECTS/TfCf/
#
# files=($(ls Data | grep "ECCITE3" | grep -v ".log"))
# for f in "${files[@]}"; do
#   echo $f
#   cp Data/$f/outs/cloupe.cloupe $f.cloupe
# done



### CREATE GENOME

cd $HOME/omicstmp

mkdir newGenome/
oldGenomePath="$GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/"

cp $oldGenome/fasta/genome.fa newGenome/genome.fa
cp $oldGenome/genes/genes.gtf newGenome/genes.gtf

x="GFP"
for x in GFP BFP; do
	pathFA=$CODEBASE/tfcf/metadata/${x}.fa

	numberBases=$(cat $pathFA | grep -v "^>" | tr -d "\n" | wc -c)
	gtf="$x\tunknown\texon\t1\t${numberBases}\t.\t+\t.\tgene_id xxx${x}xxx; transcript_id xxx${x}xxx; gene_name xxx${x}xxx; gene_biotype xxxprotein_codingxxx;"
	gtf=$(echo $gtf | sed 's/xxx/"/g')

	echo -e $gtf > ${x}.gtf

	cat ${x}.gtf >> newGenome/genes.gtf
	fold -w 60 $pathFA >> newGenome/genome.fa
	
	echo -e "" >> newGenome/genome.fa

	rm ${x}.gtf
done


grep ">" newGenome/genome.fa
tail -5 newGenome/genes.gtf
tail -200 newGenome/genome.fa


~/code/cellranger-6.0.1/cellranger mkref --genome=newGenomeExtended --fasta=newGenome/genome.fa --genes=newGenome/genes.gtf




######### BASIC ANALYSIS


cd $HOME/omicstmp

for id in ECCITE3_low_7d ECCITE3_low_14d ECCITE3_high_7d ECCITE3_high_14d; do
#for id in ECCITE3_low_7d; do
        
    echo $id
    #id="ECCITE3_low_7d"
    #id="ECCITE3_low_14d"
    
    echo "Files"
    cat $CODEBASE/tfcf/metadata/${id}_Library.csv
    
    echo "Features"
    cat $CODEBASE/tfcf/metadata/ECCITE3_Features.csv
    
	~/code/cellranger-6.0.1/cellranger count --id=$id \
	 --no-bam \
	 --libraries=$CODEBASE/tfcf/metadata/${id}_Library.csv \
	 --transcriptome=newGenomeExtended/ \
	 --feature-ref=$CODEBASE/tfcf/metadata/ECCITE3_Features.csv \
	 --localcores=24 \
	 --localmem=128 \
	 --expect-cells=10000 &> ${id}.log
    
    rm -rf $id/SC_RNA_COUNTER_CS
    
    mv ${id}* $GFS/PROJECTS/TfCf/Data/
done





######### LOOK FOR GUIDES IN EXPRESSION


cd $HOME/omicstmp

for id in ECCITE3_low_7d ECCITE3_low_14d ECCITE3_high_7d ECCITE3_high_14d; do
        
    echo $id
    #id="ECCITE3_low_7d"
    #id="ECCITE3_low_14d"
    
    echo "Files"
    cat $CODEBASE/tfcf/metadata/${id}_Library_CRISPR.csv
    
    echo "Features"
    cat $CODEBASE/tfcf/metadata/ECCITE3_Features.csv
    
	~/code/cellranger-6.0.1/cellranger count --id=${id}_CRISPR \
	 --no-bam \
	 --libraries=$CODEBASE/tfcf/metadata/${id}_Library_CRISPR.csv \
	 --transcriptome=newGenomeExtended/ \
	 --feature-ref=$CODEBASE/tfcf/metadata/ECCITE3_Features.csv \
	 --localcores=24 \
	 --localmem=128 \
	 --expect-cells=10000 &> ${id}_CRISPR.log
    
    rm -rf ${id}_CRISPR/SC_RNA_COUNTER_CS
    
    mv ${id}_CRISPR* $GFS/PROJECTS/TfCf/Data/
done