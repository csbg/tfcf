basedir=$HOME/GFS/PROJECTS/TfCf/

cd Data/Raw_data_CITESEQ1

python crukci_to_illumina.py


# ANALYZE WITH ANTIBODIES
cd $HOME/omicstmp/
mkdir nf
cd nf

~/code/cellranger-6.0.1/cellranger count --id=CITESEQ1 \
	--no-bam \
    --libraries=$CODEBASE/tfcf/metadata/CITESEQ1_Library.csv \
    --transcriptome=$HOME/GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/ \
    --feature-ref=$CODEBASE/tfcf/metadata/CITESEQ1_Features.csv \
    --localcores=24 \
    --localmem=64 \
    --expect-cells=10000 &> CITESEQ1.log

mkdir -p ~/GFS/PROJECTS/TfCf/Data/CITESEQ1/
cp -R CITESEQ1/outs ~/GFS/PROJECTS/TfCf/Data/CITESEQ1/

# ANALYZE WITHOUT ANTIBODIES
$CODEBASE/cellranger-6.0.1/cellranger count --id=CITESEQ1_RNAonly \
 --no-bam \
 --transcriptome=$GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/ \
 --libraries=$CODEBASE/tfcf/metadata/CITESEQ1_Library_RNAonly.csv \
 --localcores=24 \
 --localmem=64 \
 --expect-cells=10000 &> CITESEQ1_RNAonly.log


# gunzip -c SLX-20609.SINTG12.H5JH3DRXY.s_2.r_2.fq.gz | head -20000 | perl -ne 'print unless (0 != ($.-2) % 4)' | rev | cut -c58-84
#
#
# TCTGTC								# End of Nextera Read 1
# TCACGGCGTGCTACGGTACATTGATCC 		# 10 barcode + UMI
# AAACGATCCTGGCCGGAATTTCG				# Capture sequence
# CTTACAGATTTGATGCCGCCAAGAACGAGTAGCT	# sgRNA?
#
#
# TGCACTGACAACCAA
#
# TCTGTC
# GTAGGCAAGCAGGACGGAAGGAGCTGA
# CAACGATCCTGGCCGGAATTTCG
# AGAAGGTAC
# TGCACTGACAACCAA		#guide
# CGTATACGCG
#
# TGTTCATCAA
# CCAGAGGTGC
# CCTCAAAG
#
#   -trim START_TRIM, --start-trim START_TRIM
#                         Number of bases to discard from read2.
# ???
#   -cbf CB_FIRST, --cell_barcode_first_base CB_FIRST
#                         Postion of the first base of your cell barcodes.
# 7
#   -cbl CB_LAST, --cell_barcode_last_base CB_LAST
# 22
#                         Postion of the last base of your cell barcodes.
#   -umif UMI_FIRST, --umi_first_base UMI_FIRST
#                         Postion of the first base of your UMI.
# 23
#   -umil UMI_LAST, --umi_last_base UMI_LAST
#                         Postion of the last base of your UMI.
# 34
#
#
# the SLX-20609 but am