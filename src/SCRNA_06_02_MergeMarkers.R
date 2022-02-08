source("src/00_init.R")
out <- dirout("SCRNA_06_02_MergeMarkers")



# load datasets -----------------------------------------------------------
mobjs <- list()
for(tissuex in PATHS$SCRNA$MONOCLE.NAMES){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}


# load singleR ------------------------------------------------------------
all.samples <- gsub("md5sum_(.+).txt", "\\1", list.files(dirout_load("SCRNA_05_01_SingleR")(""), pattern="md5sum"))

sx <- all.samples[1]
singleR.res <- list()
for(sx in all.samples){
  ff <- list.files(dirout_load("SCRNA_05_01_SingleR")(sx), pattern = "cell_types_.*.csv", full.names = TRUE)
  singleR.res[[sx]] <- setNames(lapply(ff, fread), gsub("cell_types_(.+).csv", "\\1", basename(ff)))
}
singleR.res <- rbindlist(lapply(singleR.res, function(xx) rbindlist(xx, fill=TRUE,idcol = "db")), fill=TRUE, idcol="sample")
singleR.res[, cellname := paste0(cell, "_", sample)]


# Load marker signatures --------------------------------------------------
tx <- names(mobjs)[1]
marker.signatures <- list()
for(tx in names(mobjs)){
  ff <- list.files(dirout_load("SCRNA_06_01_Markers")(tx), pattern="Signatures_Larry.csv", full.names = TRUE)
  marker.signatures[[tx]] <- as.matrix(read.csv(ff))
}



# merge results -----------------------------------------------------------
tx <- names(mobjs)[1]
tx <- "in.vivo"
for(tx in names(mobjs)){
  
  monocle.obj <- mobjs[[tx]]
  
  # load relevant singleR dataset and assign clusters
  sx <- copy(singleR.res[db %in% c("izzo_label.main", "marrow10x_label.main")])
  sx$Cluster <- as.character(monocle.obj@clusters$UMAP$clusters[sx$cellname])
  sx <- sx[!is.na(Cluster)]
  
  # get clsuters where tabula muris (10x) predictions are useful (manually defined)
  sxn <- sx[,.N, by=c("Cluster", "db", "labels")]
  sxn[,total := sum(N), by=c("db", "Cluster")]
  sxn[,fraction := N/total]
  cl.tm <- sxn[db == "marrow10x_label.main" & labels %in% c("Granulocyte progenitors", "Granulocytes", "Immature B cells")][fraction > 0.8]$Cluster
  
  # Combine both predictions
  xDT <- rbind(
    sx[Cluster %in% cl.tm][db == "marrow10x_label.main"],
    sx[!Cluster %in% cl.tm][db == "izzo_label.main"]
  )
  xDT[, db := "merged"]
  
  # clean up labels
  xDT[labels == "IMP", labels := "Mono"]
  xDT[labels == "Neu", labels := "CMP"]
  
  # Split Erys into Ery and MEPs
  xDT.ery <- xDT[grepl("Ery", labels, ignore.case = TRUE)]
  sigs <- marker.signatures[[tx]][xDT.ery$cellname,c("MEP", "EryA")]
  sigs.mep <- sigs[,"MEP"] - sigs[,"EryA"] > 0
  xDT.ery[sigs.mep, labels := "MEP"]
  xDT.ery[!sigs.mep, labels := "Ery"]
  xDT <- rbind(
    xDT[!cellname %in% xDT.ery$cellname],
    xDT.ery
  )
  
  # Majority vote by cluster
  xDT.vote <- xDT[,.N, by=c("Cluster", "labels")]
  xDT.vote[, sum := sum(N), by="Cluster"]
  xDT.vote[, fraction := N/sum]
  ggplot(xDT.vote, aes(x=Cluster, y=fraction, color=labels, shape=labels)) + geom_point() + theme_bw(12) +
    scale_shape_manual(values=rep(c(1,16,2,18,3,4), 20))
  ggsave(out("CellTypes_", tx, ".pdf"), w=10, h=5)
  
  # Figure out which cells to reassign (< 10% in a cluster)
  xDT.majority <- xDT.vote[order(fraction, decreasing = TRUE)][,head(.SD, n=1), by="Cluster"][, c("Cluster", "labels")]
  xDT.keep <- xDT.vote[fraction >= 0.1][, c("Cluster", "labels"), with=F]
  xDT.reassign <- xDT.vote[fraction < 0.1][, c("Cluster", "labels"), with=F]
  
  # finalize
  xDT.reassign <- merge(xDT.reassign, xDT, by=c("Cluster", "labels"))[,-"labels"]
  xDT.reassign <- merge(xDT.reassign, xDT.majority, by=c("Cluster"))
  xDT.keep <- merge(xDT.keep, xDT, by=c("Cluster", "labels"))
  xDT.final <- rbind(xDT.keep, xDT.reassign)
  stopifnot(nrow(xDT.final) == nrow(xDT))
  
  # Majority vote by cluster
  xDT.vote <- xDT.final[,.N, by=c("Cluster", "labels")]
  xDT.vote[, sum := sum(N), by="Cluster"]
  xDT.vote[, fraction := N/sum]
  ggplot(xDT.vote, aes(x=Cluster, y=fraction, color=labels, shape=labels)) + geom_point() + theme_bw(12) +
    scale_shape_manual(values=rep(c(1,16,2,18,3,4), 20)) +
    ylim(0,1)
  ggsave(out("CellTypes_", tx, "_final.pdf"), w=10, h=5)
  
  # check that each cell only occurs once
  stopifnot(!any(duplicated(xDT.final$cellname)))
  
  saveRDS(xDT.final, out("CellTypes_", tx, ".RDS"))
}
