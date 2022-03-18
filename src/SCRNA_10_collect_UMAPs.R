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

# celltypes from singleR after manual curation
ff <- list.files(dirout_load("SCRNA_06_02_MergeMarkers")(""), pattern="CellTypes_*", full.names = TRUE)
ff <- ff[grepl(".RDS$", ff)]
singleR.cell.types <- do.call(rbind, lapply(ff, readRDS))
singleR.cell.types[, conf := tuning_scores_first - tuning_scores_second]
singleR.cell.types <- setNames(
  singleR.cell.types[, c("cellname", "sample", "labels", "conf"),with=F],
  c("rn", "sample", "functional.cluster", "functional.cluster.conf")
  )
saveRDS(singleR.cell.types, out("ProjMonocle_celltypes.RDS"))


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

# 
# # Cross project to original in vivo ---------------------------------------
# cp.vivo <- readRDS(out("ProjVivo.RDS"))
# cp.original <- readRDS(out("ProjMonocle.RDS"))
# cp.original <- cp.original[tissue == "in.vivo"]
# umap.to.matrix <- function(x){
#   m <- as.matrix(x[,c("UMAP_1", "UMAP_2"),with=F])
#   row.names(m) <- x$rn
#   m
# }
# tx <- "ex.vivo"
# for(tx in c("ex.vivo", "leukemia")){
#   
#   # For each cell identify the clostest k cells
#   qey <- umap.to.matrix(cp.vivo[tissue == tx])
#   ref <- umap.to.matrix(cp.vivo[tissue == "in.vivo"])
#   
#   qey <- qey[1:1e3,]
#   final.distance <- foreach(x = split(row.names(qey), rep(LETTERS, each=1e3)[1:nrow(qey)])) %dopar% {
#     mypdist  <- pdist(qey[x,], ref)
#     d <- do.call(rbind, split(mypdist@dist, rep(row.names(qey[x,]), each=nrow(ref))))
#     colnames(d) <- row.names(ref)
#     d
#   }
#   
#   
#   # Calculate mean original for closest cells
#   
# }
# 
# 
# # mat1 <- data.frame(x=sample(1:10000,5), 
# #                    y=sample(1:10000,5), 
# #                    z=sample(1:10000,5))
# # row.names(mat1) <- LETTERS[1:5]
# # mat2 <- data.frame(x=sample(1:100,3), 
# #                    y=sample(1:100,3), 
# #                    z=sample(1:1000,3))
# # row.names(mat2) <- LETTERS[6:8]
# # 
# # dist(rbind(mat1, mat2))
# # 
# # euclidean_distance <- function(p,q){
# #   sqrt(sum((p - q)^2))
# # }
# # install.packages("pdist")
# # require(pdist)
# # mypdist  <- pdist(mat1, mat2)
# # 
# # d <- do.call(rbind, split(mypdist@dist, rep(row.names(mat1), each=nrow(mat2))))
# # colnames(d) <- row.names(mat2)
