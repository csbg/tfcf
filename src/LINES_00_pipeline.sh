source $CODEBASE/tfcf/setup.sh

# If genome doesn't exist, then create it
if [ ! -f "$HOME/omicstmp/genomeLINES_DONE/fasta/genome.fa" ]; then
    bash $CODE/src/PREP_Genome_LINES.sh
fi

cd $HOME/omicstmp

id="ECCITE7"
for id in ECCITE2 ECCITE3_high_14d ECCITE3_high_7d ECCITE3_low_14d ECCITE3_low_7d ECCITE4_Cas9 ECCITE4_WT ECCITE5 ECCITE6 ECCITE7_Lib1Rep1 ECCITE7_Lib1Rep2 ECCITE7_Lib2 ECCITE8_IRR_d14 ECCITE8_OP1_d7 ECCITE8_OP1_d9 ECCITE8_OP2_d7 ECCITE8_OP2_d9 ECCITE8_OP3_d7 ECCITE8_OP3_d9 ECCITE8_Toxin_d14; do
	
	echo $id
	name=${id}_LINES
	
	lib=$CODE/metadata/${id}_Library.csv
	id_base=$(echo $id | sed "s/_.*//g")
	feat="$CODE/metadata/${id_base}_Features.csv"
	feat2="$CODE/metadata/${id_base}_Features.ensgs.csv"
	if [ -f feat2 ]; then
	    feat=feat2
	fi
    
    echo "Files"
    cat $lib

    echo "Features"
    cat $feat

    ~/code/cellranger-6.0.1/cellranger count --id=$name \
     --no-bam \
     --libraries=$lib \
     --transcriptome=genomeLINES_DONE/ \
     --feature-ref=$feat \
     --localcores=24 \
     --localmem=128 \
     --expect-cells=10000 &> ${name}.log

    mkdir -p ~/GFS/PROJECTS/TfCf/Data/$name
    mv $name/outs ~/GFS/PROJECTS/TfCf/Data/$name
	
done





