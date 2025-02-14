source $CODEBASE/tfcf/setup.sh

# If genome doesn't exist, then create it
if [ ! -f "$HOME/omicstmp/newGenomeExtended/fasta/genome.fa" ]; then
    bash $CODE/src/PREP_Genome.sh
fi


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



# ######### BASIC ANALYSIS
#
#
# cd $HOME/omicstmp
#
# for id in ECCITE3_low_7d ECCITE3_low_14d ECCITE3_high_7d ECCITE3_high_14d; do
# #for id in ECCITE3_low_7d; do
#
#     echo $id
#     #id="ECCITE3_low_7d"
#     #id="ECCITE3_low_14d"
#
#     echo "Files"
#     cat $CODEBASE/tfcf/metadata/${id}_Library.csv
#
#     echo "Features"
#     cat $CODEBASE/tfcf/metadata/ECCITE3_Features.csv
#
# 	~/code/cellranger-6.0.1/cellranger count --id=$id \
# 	 --no-bam \
# 	 --libraries=$CODEBASE/tfcf/metadata/${id}_Library.csv \
# 	 --transcriptome=newGenomeExtended/ \
# 	 --feature-ref=$CODEBASE/tfcf/metadata/ECCITE3_Features.csv \
# 	 --localcores=24 \
# 	 --localmem=128 \
# 	 --expect-cells=10000 &> ${id}.log
#
#     rm -rf $id/SC_RNA_COUNTER_CS
#
#     mv ${id}* $GFS/PROJECTS/TfCf/Data/
# done
#
#
#
#
#
# ######### LOOK FOR GUIDES IN EXPRESSION
#
#
# cd $HOME/omicstmp
#
# for id in ECCITE3_low_7d ECCITE3_low_14d ECCITE3_high_7d ECCITE3_high_14d; do
#
#     echo $id
#     #id="ECCITE3_low_7d"
#     #id="ECCITE3_low_14d"
#
#     echo "Files"
#     cat $CODEBASE/tfcf/metadata/${id}_Library_CRISPR.csv
#
#     echo "Features"
#     cat $CODEBASE/tfcf/metadata/ECCITE3_Features.csv
#
# 	~/code/cellranger-6.0.1/cellranger count --id=${id}_CRISPR \
# 	 --no-bam \
# 	 --libraries=$CODEBASE/tfcf/metadata/${id}_Library_CRISPR.csv \
# 	 --transcriptome=newGenomeExtended/ \
# 	 --feature-ref=$CODEBASE/tfcf/metadata/ECCITE3_Features.csv \
# 	 --localcores=24 \
# 	 --localmem=128 \
# 	 --expect-cells=10000 &> ${id}_CRISPR.log
#
#	mkdir -p $DATA/$id
#	mv $id/outs $DATA/$id
# done


######### BASIC ANALYSIS without crispr guides


cd $HOME/omicstmp

for id_original in ECCITE3_low_7d ECCITE3_low_14d ECCITE3_high_7d ECCITE3_high_14d; do
#for id in ECCITE3_low_7d; do

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
done