source $CODEBASE/tfcf/setup.sh

# setup
# python3 -m venv python_venvs/macs3
# source python_venvs/macs3/bin/activate
# pip install macs3
# deactivate


# Check if files are paired-end
# cd $RAWDATA/Raw_ChIP/
#
# files=($(ls *.bam))
#
# for f in "${files[@]}"; do
#     echo $f
#     paired=$(samtools view -c -f 1 $f)
#     total=$(samtools view -c $f)
#     echo "total:$total paired:$paired"
# done

cd $RAWDATA/Raw_ChIP/
inName1="DM_Men1_150521_S7.sort.rmdup.rmblackls.rmchr.bam"
inName2="DM_Men1_150521_S39.sort.rmdup.rmblackls.rmchr.bam"
outName="DM_Men1_150521.sort.rmdup.rmblackls.rmchr.bam"
samtools merge ${outName} ${inName1} ${inName2}
samtools index ${outName}


source $HOME/python_venvs/macs3/bin/activate

mkdir $DATA/ChIP_Peaks/
cd $_

ls $RAWDATA/Raw_ChIP/ | grep -v "bai$" | grep "bam$" > bam.files.txt
sampleGroups=$(sed "s/_.*$//g" bam.files.txt | sort | uniq)

for g in ${sampleGroups[@]}; do
    echo $g
    
    control=$(grep "^${g}_" bam.files.txt | grep "input")
    
    if [ -z "$control" ]; then
        control=$(grep "^${g}_" bam.files.txt | grep "IgG")
    fi
    
    if [ ! -z "$control" ]; then
        
        echo "CONTROL: $control"
        
        samples=$(grep "^${g}_" bam.files.txt | grep -v $control)
        for s in ${samples[@]}; do
            f=$(echo $s | sed "s/.sort.rmdup.rmblackls.rmchr.bam//g")
            echo "SAMPLE: $f"
            macs3 callpeak -t $RAWDATA/Raw_ChIP/$f.sort.rmdup.rmblackls.rmchr.bam -c $RAWDATA/Raw_ChIP/$control -f BAMPE -g mm -n $f &> $f.log
        done
    fi
done


deactivate