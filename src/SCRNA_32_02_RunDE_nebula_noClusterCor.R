source("src/00_init.R")

require(nebula)
require(doMC)
source("src/FUNC_Monocle_PLUS.R")

basedir <- "SCRNA_32_DE_Nebula_simple/"
outB <- dirout(basedir)


# Sample annotation -------------------------------------------------------
SANN <- fread(PATHS$SCRNA$ANN)

TISSUES <- PATHS$SCRNA$MONOCLE.NAMES

# load datasets -----------------------------------------------------------
mobjs <- list()
for(tissuex in TISSUES){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}


# Which analysis to take celltypes from? -----------------------------------
ANALYSIS <- "monocle.singleR"

# Analysis of most/all cells in vivo ---------------------------------------------------------------
cells <- list(
  in.vivo = fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("Guides_Fisher_Mixscape_basic_in.vivo_noMixscape_Cells.tsv"))$rn,
  ex.vivo = fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("Guides_Fisher_Mixscape_basic_ex.vivo_noMixscape_Cells.tsv"))$rn,
  leukemia = fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("Guides_Fisher_Mixscape_basic_leukemia_noMixscape_Cells.tsv"))$rn
)
ct.use <- ""
(tissue.name <- "in.vivo")
#for(tissue.name in "leukemia"){
#for(tissue.name in TISSUES){
for(tissue.name in "in.vivo"){
  
  (timex <- mobjs[[tissue.name]]$timepoint[1])
  for(timex in "28d"){ #unique(mobjs[[tissue.name]]$timepoint)){
  #for(timex in unique(mobjs[[tissue.name]]$timepoint)){

    # Monocle object
    monocle.obj <- mobjs[[tissue.name]]
    monocle.obj <- monocle.obj[, monocle.obj$timepoint == timex]
    monocle.obj <- monocle.obj[, monocle.obj$sample != "WT-LSK_OP0_NM_7d_1"]
    monocle.obj <- monocle.obj[, colnames(monocle.obj) %in% cells[[tissue.name]]]
    if(ncol(monocle.obj) == 0) next
  
    # Filter only celltypes of interest
    ann <- fread(dirout_load(paste0("SCRNA_20_Summary/",tissue.name, "_", ANALYSIS))("Annotation.tsv"))
    ann <- ann[rn %in% colnames(monocle.obj)]
    sort(unique(ann$Clusters))
    # for(x in c("B.cell", "CLP", "T.cell")){ann <- ann[!grepl(x, Clusters)]}
    # sort(unique(ann$Clusters))
  
    # Assign clusters to use as covariate (to avoid seeing shifts in populations but changes within populations)
    monocle.obj <- monocle.obj[, ann$rn]
    monocle.obj$clusterDE <- ann$Clusters
    
    out <- dirout(paste0(basedir, tissue.name, "_", timex, "_", ct.use))
  
    source("src/SCRNA_32_01_DE_FUNC_nebula_noClusterCor.R")
  }
}
