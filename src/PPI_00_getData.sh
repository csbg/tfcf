
outdir=$GFS/PROJECTS/TfCf/Analysis/PPI_00_getData/
echo $outdir
mkdir $outdir
cd $outdir

wget -O hu.MAP2.Drew.2021.MolSysBio.gz http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_ppis_genename_20200821.pairsWprob.gz
gunzip hu.MAP2.Drew.2021.MolSysBio.gz

wget -O CFMS.Skinnider.2021.NatureMethods.zip https://zenodo.org/record/4245282/files/Interactomes.zip
unzip -l CFMS.Skinnider.2021.NatureMethods.zip
unzip -p CFMS.Skinnider.2021.NatureMethods.zip Interactomes/Human/CF-MS-interactome.tsv > CFMS.Skinnider.2021.NatureMethods.tsv

wget -O PCP.SILAM.Skinnider.2021.Cell.xlsx https://ars.els-cdn.com/content/image/1-s2.0-S0092867421007042-mmc2.xlsx

wget -O depmap.CRISPR.csv --no-check-certificate https://ndownloader.figshare.com/files/27902226
wget -O depmap.ann.csv --no-check-certificate https://ndownloader.figshare.com/files/27902376

wget -O hippie.txt http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_v2_1.txt

wget -O CorumCore.txt.zip --no-check-certificate https://mips.helmholtz-muenchen.de/corum/download/coreComplexes.txt.zip

wget -O CorumAll.txt.zip --no-check-certificate https://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip