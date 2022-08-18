source $CODEBASE/tfcf/setup.sh


# If genome doesn't exist, then create it
if [ ! -f "$HOME/omicstmp/newGenomeExtended/fasta/genome.fa" ]; then
    bash $CODE/src/PREP_Genome.sh
fi

# Generate library files for the first set of samples
for id in ECCITE8_IRR_d14 ECCITE8_OP1_d7 ECCITE8_OP1_d9 ECCITE8_OP2_d7 ECCITE8_OP2_d9 ECCITE8_OP3_d7 ECCITE8_OP3_d9 ECCITE8_Toxin_d14; do
	file=~/code/tfcf/metadata/${id}_Library.csv
	id2=$(echo $id | sed "s/ECCITE8_//g")
	touch $file
	echo "fastqs,sample,library_type" > $file
	echo "/usr/local/AGFORTELNY//PROJECTS/TfCf/NewData/ECCITE8_1/demux_fastq/,3mRNA_$id2,Gene Expression" >> $file
    echo "/usr/local/AGFORTELNY//PROJECTS/TfCf/NewData/ECCITE8_1/demux_fastq/,CRISPR_$id2,CRISPR Guide Capture" >> $file
    
    if [ $id == "ECCITE8_IRR_d14" ]; then
        echo "/usr/local/AGFORTELNY//PROJECTS/TfCf/NewData/ECCITE8_2/demux_fastq/,3mRNA_IRR_OP1_d14,Gene Expression" >> $file
    fi
	
    echo "/usr/local/AGFORTELNY//PROJECTS/TfCf/NewData/ECCITE8_3/demux_fastq/,3mRNA_$id2,Gene Expression" >> $file

    if [[ $id =~ "OP" ]]; then
        echo "/usr/local/AGFORTELNY//PROJECTS/TfCf/NewData/ECCITE8_4/demux_fastq/,3mRNA_$id2,Gene Expression" >> $file
    fi
done



# Generate library files for the second set of samples
for id in ECCITE8_IRR_OP1_d28_rep1 ECCITE8_IRR_OP1_d28_rep2; do
	file=~/code/tfcf/metadata/${id}_Library.csv
	id2=$(echo $id | sed "s/ECCITE8_//g")
	touch $file
	echo "fastqs,sample,library_type" > $file
	echo "/usr/local/AGFORTELNY//PROJECTS/TfCf/NewData/ECCITE8_2/demux_fastq/,3mRNA_${id2}_28102021,Gene Expression" >> $file
	echo "/usr/local/AGFORTELNY//PROJECTS/TfCf/NewData/ECCITE8_2/demux_fastq/,CRISPR_${id2}_28102021,CRISPR Guide Capture" >> $file
done


########## BASIC ANALYSIS
cd $HOME/omicstmp

# ECCITE8_IRR_d14 ECCITE8_OP1_d7 ECCITE8_OP1_d9 ECCITE8_OP2_d7 ECCITE8_OP2_d9 ECCITE8_OP3_d7 ECCITE8_OP3_d9 ECCITE8_Toxin_d14 ECCITE8_IRR_OP1_d28_rep1 ECCITE8_IRR_OP1_d28_rep2
for id in ECCITE8_OP1_d7 ECCITE8_OP1_d9 ECCITE8_OP2_d7 ECCITE8_OP2_d9 ECCITE8_OP3_d7 ECCITE8_OP3_d9; do

    echo $id

    echo "Files"
    cat $CODEBASE/tfcf/metadata/${id}_Library.csv

    echo "Features"
    cat $CODEBASE/tfcf/metadata/ECCITE8_Features_ensgs.csv

    ~/code/cellranger-6.0.1/cellranger count --id=$id \
     --no-bam \
     --libraries=$CODEBASE/tfcf/metadata/${id}_Library.csv \
     --transcriptome=newGenomeExtended/ \
     --feature-ref=$CODEBASE/tfcf/metadata/ECCITE8_Features_ensgs.csv \
     --localcores=24 \
     --localmem=128 \
     --expect-cells=10000 &> ${id}.log

    mkdir -p ~/GFS/PROJECTS/TfCf/Data/$id
    mv $id/outs ~/GFS/PROJECTS/TfCf/Data/$id
done



######### BASIC ANALYSIS without RNA
# cd $HOME/omicstmp
# 
# # ECCITE8_IRR_d14 ECCITE8_OP1_d7 ECCITE8_OP1_d9 ECCITE8_OP2_d7 ECCITE8_OP2_d9 ECCITE8_OP3_d7 ECCITE8_OP3_d9 ECCITE8_Toxin_d14 ECCITE8_IRR_OP1_d28_rep1 ECCITE8_IRR_OP1_d28_rep2
# for id_original in ECCITE8_OP1_d7 ECCITE8_OP1_d9 ECCITE8_OP2_d7 ECCITE8_OP2_d9 ECCITE8_OP3_d7 ECCITE8_OP3_d9; do
# 
#     echo $id_original
# 	id="${id_original}_onlyRNA"
# 	echo $id
# 
# 	grep -v "CRISPR Guide Capture" $CODE/metadata/${id_original}_Library.csv > $CODE/metadata/${id}_Library.csv
# 	
# 	
#     echo "Files"
#     cat $CODEBASE/tfcf/metadata/${id}_Library.csv
#     
# 	~/code/cellranger-6.0.1/cellranger count --id=$id \
# 	 --no-bam \
# 	 --libraries=$CODEBASE/tfcf/metadata/${id}_Library.csv \
# 	 --transcriptome=newGenomeExtended/ \
# 		 --localcores=24 \
# 		 --localmem=128 \
# 	 --expect-cells=10000 &> ${id}.log
#     
# 	mkdir -p $DATA/$id
# 	mv $id/outs $DATA/$id
# done



