source("src/00_init.R")

require(nebula)
source("src/FUNC_Monocle_PLUS.R")

basedir <- "SCRNA_30_DE/"
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
ANALYSIS <- "monocle"

# Normal Monocle analysis ---------------------------------------------------------------
ct.use <- "myeloid"
(tissue.name <- names(mobjs)[2])
for(tissue.name in names(mobjs)){

  # Monocle object
  monocle.obj <- mobjs[[tissue.name]]
  if(tissue.name == "in.vivo") monocle.obj <- monocle.obj[, monocle.obj$timepoint == "14d"]
  
  # Assign clusters to use as covariate (to avoid seeing shifts in populations but changes within populations)
  monocle.obj$clusterDE <- getCL(monocle.obj)
  
  # Filter only celltypes of interest
  ann <- fread(dirout_load(paste0("SCRNA_20_Summary/",tissue.name, "_", ANALYSIS))("Annotation.tsv"))
  sort(unique(ann$Clusters))
  for(x in c("Ery", "B.cell", "CLP", "Ery", "T.cell")){ann <- ann[!grepl(x, Clusters)]}
  sort(unique(ann$Clusters))
  monocle.obj <- monocle.obj[, colnames(monocle.obj) %in% ann$rn]
  
  out <- dirout(paste0(basedir, tissue.name, "_", ct.use))

  source("src/SCRNA_30_01_DE_FUNC.R")
}


