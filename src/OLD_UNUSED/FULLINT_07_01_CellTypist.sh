conda activate celltypist

indir="$GFS/PROJECTS/TfCf/Data/ECCITE4_WT/outs/filtered_feature_bc_matrix/"
outdir="$GFS/PROJECTS/TfCf/Analysis/Celltypist/"

zcat "${indir}/features.tsv.gz" > "${outdir}/features.tsv"
zcat "${indir}/barcodes.tsv.gz" > "${outdir}/barcodes.tsv"
zcat "${indir}/matrix.mtx.gz" > "${outdir}/matrix.mtx"

mkdir $outdir
celltypist \
    -i  "${outdir}/matrix.mtx" \
    -gf "${outdir}/features.tsv" \
    -cf "${outdir}/barcodes.tsv" \
    -m  "Immune_Bonemarrow_High.pkl" \
    -o "${outdir}"

rm "${outdir}/matrix.mtx"
rm "${outdir}/features.tsv"
rm "${outdir}/barcodes.tsv"

conda deactivate