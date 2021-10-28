### CREATE GENOME

source $CODEBASE/tfcf/setup.sh

cd $HOME/omicstmp

mkdir genomeLINES/
oldGenomePath="$GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/"

cp $oldGenomePath/fasta/genome.fa genomeLINES/genome.fa
cp $GFS/PROJECTS/TfCf/CollaborationData/extendedGtf_strand.gtf genomeLINES/genes.gtf


sed -E "s/chrUn_(GL[0-9]*)/\1.1/g" genomeLINES/genes.gtf > genomeLINES/genes_clean.gtf

~/code/cellranger-6.0.1/cellranger mkref --genome=genomeLINES_DONE --fasta=genomeLINES/genome.fa --genes=genomeLINES/genes_clean.gtf

