source("src/00_init.R")
out <- dirout("SCRNA_06_02_MergeMarkers")

source("src/FUNC_Monocle_PLUS.R")

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
(tx <- names(mobjs)[1])
tx <- "leukemia"
goi <- fread("metadata/markers.csv", header = F)$V1
for(tx in names(mobjs)){
  
  monocle.obj <- mobjs[[tx]]

  # marker plot
  stopifnot(all(goi %in% row.names(monocle.obj)))
  p <- plot_cells_umap_hex_NF(monocle.obj, scale=TRUE, genes = goi, ncol=7)
  ggsave(out("CellTypes_", tx, "_Markers_UMAP_hex_scale.pdf"), w=20,h=20, plot=p)
  
  # load relevant singleR dataset and assign clusters
  sx <- copy(singleR.res[db %in% c("izzo_label.main", "marrow10x_label.main")])
  sx$Cluster <- as.character(monocle.obj@clusters$UMAP$clusters[sx$cellname])
  sx <- sx[!is.na(Cluster)]
  
  # get clusters where tabula muris (10x) predictions are useful (manually defined)
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
  xDT[labels == "Neu", labels := "GMP"]
  xDT[labels %in% c("Eo", "Ba", "E/B"), labels := "Eo/Ba"]
  xDT[, labels := gsub("progenitors?", "P", labels, ignore.case = TRUE)]
  xDT[, labels := gsub("granulocytes?", "Gran.", labels, ignore.case = TRUE)]
  xDT[, labels := gsub("immature", "Imm.", labels, ignore.case = TRUE)]
  xDT[, labels := gsub("B.cells?", "B-cell", labels, ignore.case = TRUE)]
  
  # Add expression of key genes
  mt <- log1p(NF_TPM_Matrix(monocle.obj, goi))
  xDT <- cbind(xDT, data.table(scale(t(mt[,xDT$cellname]))))
  
  
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
  xDT.majvote <- rbind(xDT.keep, xDT.reassign)
  stopifnot(nrow(xDT.majvote) == nrow(xDT))
  
  # Cell types in leukemia
  if(tx == "leukemia"){
    # Reassign GMPs with high Gata exprsesion to Eo/Ba
    gmp.reassign <- xDT.majvote[labels=="GMP"][,.(Gata2 = mean(Gata2), Gata1 = mean(Gata1)), by="Cluster"][Gata2 > 1 & Gata1 > 1]$Cluster
    xDT.majvote[labels=="GMP" & Cluster %in% gmp.reassign, labels := "Eo/Ba"]
    
    # identify leukemic stem cells
    pDT <- xDT.majvote[,.(Bcat1 = mean(Bcat1), Cd34 = mean(Cd34), Ctsg = mean(Ctsg)), by="Cluster"][order(Bcat1)]
    ggplot(pDT, aes(x=Cd34, y=Ctsg, color=Bcat1, label=Cluster)) + 
      theme_bw(12) +
      geom_point(size=3) +
      geom_text(color="black", aes(x=Cd34 + 0.05)) + 
      scale_color_gradient2()
    ggsave(out("Celltypes_LSCs.pdf"), w=6,h=5)
    xDT.majvote[Cluster %in% pDT[Bcat1 > 0.5 & Cd34 > -0.25 & Ctsg > -0.25]$Cluster, labels := "LSC"]
  }
  
  # Reassign MEPs
  # 1 by Gata expression
  gata2.cluster <- xDT.majvote[labels == "MEP"][,.(Gata1 = mean(Gata1), Gata2 = mean(Gata2)), by="Cluster"][Gata2 > Gata1]$Cluster
  xDT.majvote[Cluster %in% gata2.cluster & labels == "MEP", labels := "MEP (early)"]
  # 2 by cell cycle
  stopifnot(all(xDT.majvote$cellname == colnames(monocle.obj)[match(xDT.majvote$cellname, colnames(monocle.obj))]))
  xDT.majvote$Phase <- monocle.obj$Phase[match(xDT.majvote$cellname, colnames(monocle.obj))]
  ccDT <- xDT.majvote[labels == "MEP"][,.N,by=c("Cluster", "Phase")]
  ccDT[, sum := sum(N), by=c("Cluster")]
  ccDT[, frac := N/sum]
  ccDT <- ccDT[frac > 0.95][Phase != "G2M"]
  if(nrow(ccDT) > 0){
    for(i in 1:nrow(ccDT)){
      xDT.majvote[labels == "MEP" & Cluster == ccDT[i]$Cluster, labels := paste0("MEP ", "(", ccDT[i]$Phase, ")")]
    }
  }
  # 3 separate perturbed cluster
  xDT.majvote$guide <- monocle.obj$guide[match(xDT.majvote$cellname, colnames(monocle.obj))]
  rcor.cluster <- xDT.majvote[labels == "MEP"][,.N, by=c("Cluster", "guide")][grepl("Rcor1", guide)][,sum(N), by="Cluster"][order(V1, decreasing = TRUE)]$Cluster[1]
  rcor.samples <- unique(xDT.majvote[grepl("^Rcor", guide)]$sample)
  rcor.cnt <- xDT.majvote[Cluster == rcor.cluster & labels == "MEP"][sample %in% rcor.samples]
  if(nrow(rcor.cnt) > 0 & nrow(rcor.cnt[grepl("^Rcor", guide)]) / nrow(rcor.cnt) > 0.25) xDT.majvote[Cluster == rcor.cluster & labels == "MEP", labels := "MEP (pert.)"]
  
  # Reassign GMPs
  # 1. based on S100a8 expression (in matured ones)
  gmp.mat.clusters <- xDT.majvote[labels=="GMP"][,mean(S100a8), by=c("Cluster")][V1 > 1]$Cluster
  xDT.majvote[labels=="GMP" & Cluster %in% gmp.mat.clusters, labels := "GMP (late)"]
  # 2. based on CD34 / Ctsg expression
  gmp.cd34 <- xDT.majvote[labels == "GMP"][,.(g1=mean(Cd34), g2=mean(Ctsg)), by=c("Cluster")]
  #ggplot(gmp.cd34, aes(x=g1, y=g2)) + geom_point()
  xDT.majvote[labels=="GMP" & Cluster %in% gmp.cd34[g1 - g2 > 0]$Cluster, labels := "GMP (early)"]
  
  # Reassign HSCs to EBMP (more mature)
  ebmp.clusters <- xDT.majvote[labels=="HSC"][,.(mean(Mecom + Hoxa7 + Hoxa9)/3), by=c("Cluster")][V1 < 0.2]$Cluster
  xDT.majvote[labels=="HSC" & Cluster %in% ebmp.clusters, labels := "EBMP"]
  
  # Plot after voting
  xDT.vote <- xDT.majvote[,.N, by=c("Cluster", "labels")]
  xDT.vote[, sum := sum(N), by="Cluster"]
  xDT.vote[, fraction := N/sum]
  ggplot(xDT.vote, aes(x=Cluster, y=fraction, color=labels, shape=labels)) + geom_point() + theme_bw(12) +
    scale_shape_manual(values=rep(c(1,16,2,18,3,4), 20)) +
    ylim(0,1)
  ggsave(out("CellTypes_", tx, "_final.pdf"), w=10, h=5)
  
  # UMAP Plot
  if(!"U1" %in% colnames(xDT.majvote)) xDT.majvote <- cbind(xDT.majvote, setNames(data.table(reducedDims(monocle.obj)$UMAP[xDT.majvote$cellname,]), c("U1", "U2")))
  p <- ggplot(xDT.majvote, aes(x=U1, y=U2, color=labels)) + 
    theme_bw(12) +
    geom_point(size=0.5) +
    geom_text(data=xDT.majvote[,.(U1=median(U1), U2=median(U2)), by="labels"], aes(label=labels), color="black")
  ggsave(out("CellTypes_", tx, "_UMAP.jpg"), w=7,h=5, plot=p)
  
  # check that each cell only occurs once
  stopifnot(!any(duplicated(xDT.majvote$cellname)))
  
  xDT.majvote <- xDT.majvote[,-c("U1", "U2", "guide", "Phase", goi), with=F]
  saveRDS(xDT.majvote, out("CellTypes_", tx, ".RDS"))
}
