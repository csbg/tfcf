source("src/00_init.R")
out <- dirout("SCRNA_10_collect_UMAPs/")

require(pdist)
require(doMC)
source("src/FUNC_Monocle_PLUS.R")
registerDoMC(cores=10)

# Original UMAPs ----------------------------------------------------------
mobjs <- list()
for(tissuex in PATHS$SCRNA$MONOCLE.NAMES){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}


# Duplet scors ------------------------------------------------------------
ff <- list.files(dirout_load("SCRNA_03_01_Duplets")(""), pattern="Duplet_Scores_.*.tsv", full.names = TRUE)
ds <- lapply(ff, fread)
names(ds) <- basename(ff)
ds <- rbindlist(ds, idcol = "sample")
ds[, sample := gsub("Duplet_Scores_(.+).tsv", "\\1", sample)]
ds[, rn := paste0(rn, "_", sample)]
ds <- ds[,c("rn", "dupletScore")]


# Projection from Monocle (original) --------------------------------------
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
stopifnot(!any(is.na(res$UMAP_1)))

# Duplets
res <- merge(res, ds, by="rn", all=TRUE)[!is.na(UMAP_1)]
stopifnot(nrow(res[is.na(dupletScore)]) == 0)
table(res$dupletScore > 0.9)
ggplot(res[!is.na(tissue)], aes(x=UMAP_1, y=UMAP_2)) + 
  scale_fill_hexbin() +
  facet_wrap(~tissue, ncol=3) +
  theme_bw(12) +
  stat_summary_hex(aes(z=dupletScore),fun=mean, bins=100)
ggsave(out("ProjMonocle_Duplets.pdf"), w=12,h=4)
table(res$tissue)
#res <- res[dupletScore < 0.9 | tissue != "in.vivo"]
res <- res[dupletScore < 0.9]
table(res$tissue)

# Save
saveRDS(res, out("ProjMonocle.RDS"))

# Clusters
mnam <- names(mobjs)[1]
dDT.ct <- list()
for(mnam in names(mobjs)){
  monocle.obj <- mobjs[[mnam]]
  dDT.ct[[mnam]] <- data.table(data.frame(colData(monocle.obj)@listData), keep.rownames = TRUE)[,c("rn", "sample")][match(colnames(monocle.obj), rn)]
  dDT.ct[[mnam]]$functional.cluster <- getCL(monocle.obj)
  dDT.ct[[mnam]]$functional.cluster.conf <- 1
}
dDT.ct <- rbindlist(dDT.ct)
dDT.ct <- dDT.ct[match(res$rn, rn)]
stopifnot(!any(is.na(dDT.ct$sample)))
saveRDS(dDT.ct, out("ProjMonocle_Clusters.RDS"))




# celltypes from singleR after further curation -------------------------------------------------------
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
singleR.cell.types <- singleR.cell.types[match(res$rn, rn)]
stopifnot(!any(is.na(singleR.cell.types$sample)))
saveRDS(singleR.cell.types, out("ProjMonocle_celltypes.RDS"))


# Projection_Invivo -------------------------------------------------------
ff <- list.files(dirout_load("SCRNA_08_01_ProjectionInvivo")(""), pattern="Output_", full.names = TRUE)
dL <- lapply(ff, fread)
dDT <- rbindlist(dL, fill=TRUE)
dDT.umap <- dDT[,c("rn", "sample", "UMAP_1", "UMAP_2", "tissue"), with=F]
dDT.umap <- merge(dDT.umap, ds, by="rn", all=TRUE)[!is.na(UMAP_1)]
stopifnot(nrow(dDT.umap[is.na(dupletScore)]) == 0)
table(dDT.umap$tissue)
#dDT.umap <- dDT.umap[dupletScore < 0.9 | tissue != "in.vivo"]
dDT.umap <- dDT.umap[dupletScore < 0.9]
table(dDT.umap$tissue)
saveRDS(dDT.umap, out("ProjVivo.RDS"))

dDT.ct <- dDT[,c("rn", "sample", "functional.cluster", "functional.cluster.conf"), with=F]
dDT.ct <- dDT.ct[match(dDT.umap$rn, rn)]
stopifnot(!any(is.na(dDT.ct$sample)))
saveRDS(dDT.ct, out("ProjVivo_celltypes.RDS"))

# Crossprojection_Invivo -------------------------------------------------------
ff <- list.files(dirout_load("SCRNA_08_01_ProjectionInvivo")(""), pattern="OutputCrossprojection_", full.names = TRUE)
dL <- lapply(ff, fread)
dDT <- rbindlist(dL, fill=TRUE)
dDT.umap <- dDT[,c("rn", "sample", "UMAP_1", "UMAP_2", "tissue"), with=F]
dDT.umap.invivo <- readRDS(out("ProjMonocle.RDS"))[,-"dupletScore", with=F]
dDT.umap <- rbind(dDT.umap, dDT.umap.invivo[tissue == "in.vivo"])
dDT.umap <- merge(dDT.umap, ds, by="rn", all=TRUE)[!is.na(UMAP_1)]
stopifnot(nrow(dDT.umap[is.na(dupletScore)]) == 0)
table(dDT.umap$tissue)
#dDT.umap <- dDT.umap[dupletScore < 0.9 | tissue != "in.vivo"]
dDT.umap <- dDT.umap[dupletScore < 0.9]
table(dDT.umap$tissue)
saveRDS(dDT.umap, out("ProjVivoX.RDS"))


# Projection Izzo WT1 -----------------------------------------------------
ff <- list.files(dirout_load("SCRNA_08_03_ProjectionIzzo_separate/Izzo_WT1/")(""), pattern="Output_", full.names = TRUE)
dL <- lapply(ff, fread)
dDT <- rbindlist(dL, fill=TRUE)
dDT.umap <- dDT[,c("rn", "sample", "UMAP_1", "UMAP_2", "tissue"), with=F]
dDT.umap <- merge(dDT.umap, ds, by="rn", all=TRUE)[!is.na(UMAP_1)]
stopifnot(nrow(dDT.umap[is.na(dupletScore) & sample != "Izzo_WT1"]) == 0)
table(dDT.umap$tissue)
# dDT.umap <- dDT.umap[dupletScore < 0.9 | tissue != "in.vivo" | sample == "Izzo_WT1"]
dDT.umap <- dDT.umap[dupletScore < 0.9 | sample == "Izzo_WT1"]
table(dDT.umap$tissue)
saveRDS(dDT.umap, out("ProjIzzo.RDS"))
dDT.ct <- dDT[,c("rn", "sample", "functional.cluster", "functional.cluster.conf"), with=F]
dDT.ct <- dDT.ct[match(dDT.umap$rn, rn)]
stopifnot(!any(is.na(dDT.ct$sample)))
saveRDS(dDT.ct, out("ProjIzzo_celltypes.RDS"))



# # Export TSV for Jake -----------------------------------------------------
# dDT <- readRDS(out("ProjMonocle_celltypes.RDS"))
# write.tsv(dDT, out("ProjMonocle_celltypes.tsv"))
# table(dDT$sample)
