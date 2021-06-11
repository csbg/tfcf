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



source $HOME/python_venvs/macs3/bin/activate

mkdir $DATA/ChIP_Peaks/
cd $_


f="CUT_PU1-100221_S14"
macs3 callpeak -t $RAWDATA/Raw_ChIP/$f.sort.rmdup.rmblackls.rmchr.bam -c $RAWDATA/Raw_ChIP/CUT_IgG-100221_S16.sort.rmdup.rmblackls.rmchr.bam -f BAMPE -g mm -n $f &> $f.log

deactivate