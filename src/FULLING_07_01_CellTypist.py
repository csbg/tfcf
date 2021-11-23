import celltypist
import os
from celltypist import models

models.models_description()

models.download_models(model = ['Immune_Bonemarrow_High.pkl', 'Immune_Bonemarrow_Low.pkl'])


models.models_description(on_the_fly = True)


#Select the model from the above list. If the `model` argument is not provided, will default to `Immune_All_Low.pkl`.
model = models.Model.load(model = 'Immune_Bonemarrow_High.pkl')
#Examine cell types contained in the model.
model.cell_types
#Examine genes/features contained in the model.
model.features
#The stochastic gradient descent logistic regression classifier within the model.
model.classifier
#The standard scaler within the model (used to scale the input query data).
model.scaler
#The model information.
model.description


input_file = celltypist.samples.get_sample_csv()


indir=os.path.join(os.environ.get('GFS'), "PROJECTS/TfCf/Data/")
os.listdir(indir)
infile=os.path.join(indir, "ECCITE5", "outs", "filtered_feature_bc_matrix", "matrix.mtx.gz")
input_file = celltypist.samples.get_sample_csv(infile)


conda activate celltypist

indir="$GFS/PROJECTS/TfCf/Data/ECCITE4_WT/outs/filtered_feature_bc_matrix/"
outdir="$GFS/PROJECTS/TfCf/Analysis/Celltypist/"

zcat "${indir}/features.tsv.gz" > "${outdir}/features.tsv"
zcat "${indir}/barcodes.tsv.gz" > "${outdir}/barcodes.tsv"

mkdir $outdir
celltypist \
    -i  "${indir}/matrix.mtx.gz" \
    -gf "${outdir}/features.tsv" \
    -cf "${outdir}/barcodes.tsv" \
    -m  "Immune_Bonemarrow_High.pkl" \
    -o "$outdir"