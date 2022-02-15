source("src/00_init.R")

basedir <- "SCRNA_20_Summary/"
outB <- dirout(basedir)

source("src/FUNC_Monocle_PLUS.R")


# load datasets -----------------------------------------------------------
mobjs <- list()
for(tissuex in PATHS$SCRNA$MONOCLE.NAMES){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}


# Normal Monocle analysis ---------------------------------------------------------------
tissue.name <- names(mobjs)[2]
for(tissue.name in names(mobjs)){

  analysis.name <- "monocle"

  monocle.obj <- mobjs[[tissue.name]]

  monocle.obj$clusters.final <-  getCL(monocle.obj)

  out <- dirout(paste0("SCRNA_20_Summary/", tissue.name, "_", analysis.name))

  source("src/SCRNA_20_01_SummaryFUNC.R")
}

# Monocle analysis with predicted celltypes ---------------------------------------------------------------
(tissue.name <- names(mobjs)[1])
for(tissue.name in names(mobjs)){
  
  analysis.name <- "monocle.singleR"
  
  monocle.obj <- mobjs[[tissue.name]]
  
  singleR.cell.types <- readRDS(dirout_load("SCRNA_06_02_MergeMarkers")("CellTypes_",tissue.name,".RDS"))
  monocle.obj$clusters.final <-  setNames(singleR.cell.types[match(colnames(monocle.obj), cellname),]$labels, colnames(monocle.obj))
  
  out <- dirout(paste0("SCRNA_20_Summary/", tissue.name, "_", analysis.name))
  
  source("src/SCRNA_20_01_SummaryFUNC.R")
}


# Projection to in vivo ---------------------------------------------------------------
proj <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjVivo.RDS"))
umapMT <- as.matrix(proj[,c("UMAP_1", "UMAP_2"), with=F])
colnames(umapMT) <- NULL
row.names(umapMT) <- proj$rn
# Load cell types
cts <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjIzzo_celltypes.RDS"))
for(tissue.name in names(mobjs)){

  analysis.name <- "projectToInvivo"

  monocle.obj <- mobjs[[tissue.name]]
  monocle.obj$clusters.final <-  setNames(cts[match(colnames(monocle.obj), rn),]$functional.cluster, colnames(monocle.obj))
  reducedDims(monocle.obj)$UMAP <- umapMT[colnames(monocle.obj),]
  out <- dirout(paste0("SCRNA_20_Summary/", tissue.name, "_", analysis.name))

  source("src/SCRNA_20_01_SummaryFUNC.R")
}


# Projection to izzo ---------------------------------------------------------------
# Load UMAP
proj <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjIzzo.RDS"))
umapMT <- as.matrix(proj[,c("UMAP_1", "UMAP_2"), with=F])
colnames(umapMT) <- NULL
row.names(umapMT) <- proj$rn
# Load cell types
cts <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjIzzo_celltypes.RDS"))
for(tissue.name in names(mobjs)){

  analysis.name <- "projectToIzzo"

  monocle.obj <- mobjs[[tissue.name]]
  monocle.obj$clusters.final <-  setNames(cts[match(colnames(monocle.obj), rn),]$functional.cluster, colnames(monocle.obj))
  reducedDims(monocle.obj)$UMAP <- umapMT[colnames(monocle.obj),]
  out <- dirout(paste0("SCRNA_20_Summary/", tissue.name, "_", analysis.name))

  source("src/SCRNA_20_01_SummaryFUNC.R")
}

