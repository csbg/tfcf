source("src/00_init.R")
out <- dirout("SCRNA_10_collect_UMAPs/")

require(pdist)
require(doMC)
registerDoMC(cores=10)

# Original UMAPs ----------------------------------------------------------
mobjs <- list()
for(tissuex in PATHS$SCRNA$MONOCLE.NAMES){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}

res <- data.table()
tx <- names(mobjs)[1]
for(tx in names(mobjs)){
  monocle.obj <- mobjs[[tx]]
  dDT.umap <- cbind(
    data.table(reducedDims(monocle.obj)$UMAP, keep.rownames = TRUE),
    data.table(
      sample=monocle.obj$sample,
      tissue=monocle.obj$tissue
    ))
  res <- rbind(res, dDT.umap)
}
colnames(res)[colnames(res) %in% c("V1", "V2")] <- c("UMAP_1", "UMAP_2")
saveRDS(res, out("ProjMonocle.RDS"))

# celltypes from singleR after manual curation -------------------------------------------------------
ff <- list.files(dirout_load("SCRNA_06_02_MergeMarkers")(""), pattern="CellTypes_*", full.names = TRUE)
ff <- ff[grepl(".RDS$", ff)]
singleR.cell.types <- lapply(ff, readRDS)
leukemia.cells <- singleR.cell.types[[which(grepl("leukemia", ff))]]$cellname
singleR.cell.types <- do.call(rbind, singleR.cell.types)
singleR.cell.types[, conf := tuning_scores_first - tuning_scores_second]
singleR.cell.types <- setNames(
  singleR.cell.types[, c("cellname", "sample", "labels", "conf"),with=F],
  c("rn", "sample", "functional.cluster", "functional.cluster.conf")
  )
singleR.cell.types[rn %in% leukemia.cells & functional.cluster %in% c("GMP", "GMP (early)", "HSC", "LSC"), functional.cluster := "LSC"]
saveRDS(singleR.cell.types, out("ProjMonocle_celltypes.RDS"))


# Projection_Invivo -------------------------------------------------------
ff <- list.files(dirout_load("SCRNA_08_01_ProjectionInvivo")(""), pattern="Output_", full.names = TRUE)
dL <- lapply(ff, fread)
dDT <- rbindlist(dL, fill=TRUE)
dDT.umap <- dDT[,c("rn", "sample", "UMAP_1", "UMAP_2", "tissue"), with=F]
saveRDS(dDT.umap, out("ProjVivo.RDS"))

dDT.ct <- dDT[,c("rn", "sample", "functional.cluster", "functional.cluster.conf"), with=F]
saveRDS(dDT.ct, out("ProjVivo_celltypes.RDS"))

# Crossprojection_Invivo -------------------------------------------------------
ff <- list.files(dirout_load("SCRNA_08_01_ProjectionInvivo")(""), pattern="OutputCrossprojection_", full.names = TRUE)
dL <- lapply(ff, fread)
dDT <- rbindlist(dL, fill=TRUE)
dDT.umap <- dDT[,c("rn", "sample", "UMAP_1", "UMAP_2", "tissue"), with=F]
dDT.umap.invivo <- readRDS(out("ProjMonocle.RDS"))
dDT.umap <- rbind(dDT.umap, dDT.umap.invivo[tissue == "in.vivo"])
saveRDS(dDT.umap, out("ProjVivoX.RDS"))


# Projection Izzo WT1 -----------------------------------------------------
ff <- list.files(dirout_load("SCRNA_08_03_ProjectionIzzo_separate/Izzo_WT1/")(""), pattern="Output_", full.names = TRUE)
dL <- lapply(ff, fread)
dDT <- rbindlist(dL, fill=TRUE)
dDT.umap <- dDT[,c("rn", "sample", "UMAP_1", "UMAP_2", "tissue"), with=F]
saveRDS(dDT.umap, out("ProjIzzo.RDS"))

dDT.ct <- dDT[,c("rn", "sample", "functional.cluster", "functional.cluster.conf"), with=F]
saveRDS(dDT.ct, out("ProjIzzo_celltypes.RDS"))
