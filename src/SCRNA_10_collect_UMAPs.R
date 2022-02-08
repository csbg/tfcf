source("src/00_init.R")
out <- dirout("SCRNA_10_collect_UMAPs/")


# Projection_Invivo -------------------------------------------------------
ff <- list.files(dirout_load("SCRNA_08_01_ProjectionInvivo")(""), pattern="Output_", full.names = TRUE)
dL <- lapply(ff, fread)
dDT <- rbindlist(dL, fill=TRUE)
dDT.umap <- dDT[,c("rn", "sample", "UMAP_1", "UMAP_2", "tissue"), with=F]
saveRDS(dDT.umap, out("ProjVivo.RDS"))

dDT.ct <- dDT[,c("rn", "sample", "functional.cluster", "functional.cluster.conf"), with=F]
saveRDS(dDT.ct, out("ProjVivo_celltypes.RDS"))


# Projection Izzo WT1 -----------------------------------------------------
ff <- list.files(dirout_load("SCRNA_08_03_ProjectionIzzo_separate/Izzo_WT1/")(""), pattern="Output_", full.names = TRUE)
dL <- lapply(ff, fread)
dDT <- rbindlist(dL, fill=TRUE)
dDT.umap <- dDT[,c("rn", "sample", "UMAP_1", "UMAP_2", "tissue"), with=F]
saveRDS(dDT.umap, out("ProjIzzo.RDS"))

dDT.ct <- dDT[,c("rn", "sample", "functional.cluster", "functional.cluster.conf"), with=F]
saveRDS(dDT.ct, out("ProjIzzo_celltypes.RDS"))


# Transfer to original Izzo UMAP
#(load(dirout_load("SCRNA_05_01_SingleR")("izzo.RData")))
