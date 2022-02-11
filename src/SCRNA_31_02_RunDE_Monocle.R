source("src/00_init.R")

require(nebula)
source("src/FUNC_Monocle_PLUS.R")

basedir <- "SCRNA_31_DE_Monocle/"
outB <- dirout(basedir)


# Sample annotation -------------------------------------------------------
SANN <- fread(PATHS$SCRNA$ANN)


# load datasets -----------------------------------------------------------
mobjs <- list()
for(tissuex in PATHS$SCRNA$MONOCLE.NAMES){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}


# Which analysis to take celltypes from? -----------------------------------
ANALYSIS <- "monocle.singleR"

# Normal Monocle analysis ---------------------------------------------------------------
ct.use <- "myeloid"
(tissue.name <- names(mobjs)[2])
for(tissue.name in names(mobjs)){

  # Monocle object
  monocle.obj <- mobjs[[tissue.name]]
  monocle.obj <- monocle.obj[, monocle.obj$timepoint == "14d"]
  monocle.obj <- monocle.obj[, monocle.obj$sample != "WT-LSK_OP0_NM_7d_1"]

  # Filter only celltypes of interest
  ann <- fread(dirout_load(paste0("SCRNA_20_Summary/",tissue.name, "_", ANALYSIS))("Annotation.tsv"))
  ann <- ann[rn %in% colnames(monocle.obj)]
  sort(unique(ann$Clusters))
  for(x in c("Ery", "B.cell", "CLP", "MEP", "T.cell")){ann <- ann[!grepl(x, Clusters)]}
  sort(unique(ann$Clusters))
  
  # Assign clusters to use as covariate (to avoid seeing shifts in populations but changes within populations)
  monocle.obj <- monocle.obj[, ann$rn]
  monocle.obj$clusterDE <- ann$Clusters
  
  out <- dirout(paste0(basedir, tissue.name, "_", ct.use))

  source("src/SCRNA_31_01_DE_FUNC_Monocle.R")
}


