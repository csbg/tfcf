source("src/00_init.R")
base.dir <- "FIG_02_scRNA_UMAPs/"
outBase <- dirout(base.dir)

require(ggrepel)
require(WriteXLS)
require(patchwork)


# FUNCTIONS ---------------------------------------------------------------
ds <- function(path){load(path); return(monocle.obj)}

xu <- xlab("UMAP 1")
yu <- ylab("UMAP 2")

# SETTINGS ----------------------------------------------------------------
ff <- list.files(dirout_load("SCRNA_20_Summary")(""))
(TISSUES <- gsub("_.+", "", ff[grepl("_monocle.singleR$", ff)]))
SIGS.USE <- fread("metadata/markers.signatures.use.scRNA.tsv")
SIGS.USE[, sig := paste(DB, Celltype)]
SIGS.USE.DA <- fread("metadata/markers.signatures.use.scRNA2.tsv")
SIGS.USE.DA[, sig := paste(DB, Celltype)]

# Folders -----------------------------------------------------------------
inDir.funcs <- list()
for(tx in TISSUES){inDir.funcs[[tx]] <- dirout_load(paste0("SCRNA_20_Summary/", tx, "_monocle.singleR"))}

source("src/FUNC_Monocle_PLUS.R")

# Load data ---------------------------------------------------------------

# Annotations
SANN <- fread(PATHS$SCRNA$ANN)

# Marker genes
marker.genes <- fread("metadata/markers.csv", header = FALSE)$V1

# Load signatures
marker.signatures <- lapply(TISSUES, function(tissue.name){
  ff <- list.files(dirout_load("SCRNA_06_01_Markers")(tissue.name), pattern="Signatures_", full.names = TRUE)
  names(ff) <- gsub("^Signatures_(.+?).csv$", "\\1", basename(ff))
  ff <- ff[c(unique(SIGS.USE$DB), "BulkOld", "BulkNew")]
  lapply(ff, function(fx) as.matrix(read.csv(fx)))
})
names(marker.signatures) <- TISSUES
marker.signatures <- rbindlist(lapply(marker.signatures, function(ltx){
  rbindlist(lapply(ltx, function(dt){
      dt <- melt(data.table(dt, keep.rownames = TRUE), id.vars = "rn")
    }), idcol = "DB")
  }), idcol = "tissue")
marker.signatures.bulk <- marker.signatures[DB == "BulkOld"]
marker.signatures.bulk2 <- marker.signatures[DB == "BulkNew"]
marker.signatures <- merge(marker.signatures, SIGS.USE, by.x=c("DB", "variable"), by.y=c("DB", "Celltype"))

# signature differential analysis
# marker.signatures.DA <- lapply(TISSUES, function(tx) fread(inDir.funcs[[tx]]("SigDA.tsv")))
# names(marker.signatures.DA) <- TISSUES

# Marker overlaps
marker.overlaps <- fread(dirout_load("SCRNA_06_03_MarkerOverlaps")("Enrichments.tsv"))

# Cell annotations
annList <- lapply(names(inDir.funcs), function(inDir.current){
  ann <- fread(inDir.funcs[[inDir.current]]("Annotation.tsv"))
  #ann[, perturbed := !(mixscape_class.global %in% c("NP", "NTC") | is.na(mixscape_class.global))]
  ann[, gene := gsub("_.+", "", guide)]
  ann
  })
annList <- rbindlist(annList, fill=TRUE)
annList <- annList[Clusters != "unclear"]

# numeric clusters
ann.numeric <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjMonocle_Clusters.RDS"))
annList$Cluster.number <- ann.numeric[match(annList$rn, rn)]$functional.cluster

# leukemia update clusters
ann.leukemia.original <- readRDS(dirout_load("SCRNA_06_02_MergeMarkers")("CellTypes_leukemia.RDS"))
lsc.clusters <- ann.leukemia.original[, .N, by=c("Cluster", "labels")][order(N,decreasing = TRUE)][,head(.SD,n=1),by=c("Cluster")][labels == "LSC"]$Cluster
ct.centroids <- annList[tissue == "leukemia" & Clusters %in% c("GMP (late)", "Mono", "Eo/Ba")][, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by="Clusters"]
cl.centroids <- annList[tissue == "leukemia"][, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by="Cluster.number"]
mt <- as.matrix(dist(rbind(
  as.matrix.data.frame(data.frame(ct.centroids[,c("UMAP1", "UMAP2"),with=F], row.names = ct.centroids$Clusters)),
  as.matrix.data.frame(data.frame(cl.centroids[,c("UMAP1", "UMAP2"),with=F], row.names = paste0("cl", cl.centroids$Cluster.number)))
)))[ct.centroids$Clusters, paste0("cl", cl.centroids$Cluster.number)]
cl.celltypes <- setNames(row.names(mt)[apply(mt, 2, which.min)], colnames(mt))
cl.celltypes[paste0("cl", lsc.clusters)] <- "LSC"
annList$Clusters.new <- cl.celltypes[paste0("cl", annList$Cluster.number)]
annList[tissue == "leukemia", Clusters := Clusters.new]
annList$Clusters.new <- NULL

# Update cell types ex vivo
annList[tissue == "ex.vivo" & Clusters == "MkP", Clusters := "MEP (early)"]
annList[tissue == "ex.vivo" & Clusters == "MEP", Clusters := "MEP (early)"]


# Other projections
umap.proj <- list(
  original=readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjMonocle.RDS")),
  izzo = readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjIzzo.RDS")),
  in.vivo = readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjVivo.RDS")),
  in.vivo.X = readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjVivoX.RDS"))
)

# DLA factors to plot
dla.table <- fread("metadata/FIGS_02_CFs.main.txt")
dla.healthy <- list(
  supp=dla.table$CF,
  main=with(dla.table, CF[LargePlot]),
  main.small=with(dla.table, CF[SmallPlot]),
  main.Feb28=fread("metadata/FIGS_02_CFs.main_Feb28.txt", header = FALSE)$V1
)
dla.healthy$all <- setdiff(sort(unique(annList$gene)), "NTC")
#for(x in names(dla.healthy)){ dla.healthy[[x]] <- setdiff(dla.healthy[[x]], "Smarcb1")}
dla.cancer <- list(
  supp=dla.healthy$all,
  main=c("Setd1b","Wdr82","Kmt2d","Kmt2a","Men1","Prmt1","Carm1","Prmt5","Smarcb1","Smarcd2","Smarcd1","Brd9","Pbrm1","Smarca5","Chd4","Rbbp4","Hdac3","Setdb1","Atf7ip","Setdb2","Hdac1","Rcor1","Rcor2","Ezh2","Stag2","Ash1l")
)

# SIMPLE SETUP ENDS HERE ---------------------------------------------------------




# load Monocle Objects
mobjs <- list()
(tissuex <- PATHS$SCRNA$MONOCLE.NAMES[3])
for(tissuex in PATHS$SCRNA$MONOCLE.NAMES){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}
# add gene names
for(tissuex in names(mobjs)){
  mobjs[[tissuex]]$CellType <- annList[match(colnames(mobjs[[tissuex]]), rn)]$Clusters
  fData(mobjs[[tissuex]])$gene_short_name <- row.names(fData(mobjs[[tissuex]]))
}
# Remove duplets
for(tissuex in names(mobjs)){
  mobjs[[tissuex]] <- mobjs[[tissuex]][,colnames(mobjs[[tissuex]]) %in% annList$rn]
}

izzo <- fread(dirout_load("SCRNA_08_03_ProjectionIzzo_separate")("Izzo_WT1/Output_izzo.tsv"))


# SETUP ENDS HERE ---------------------------------------------------------




# GENERAL PLOTS -----------------------------------------------------------


# . Marker overlaps -------------------------------------------------------
pDT <- marker.overlaps[database %in% c("IzzoEtAl", "LarrySelf") | (database == "PanglaoDB" & geneset == "Hematopoietic stem cells")]
pDT <- pDT[list != "FcgR3"]
ggplot(pDT, aes(x=list, y=geneset, color=log2OR_cap, size=pmin(5, -log10(padj)))) + 
  themeNF(rotate = TRUE) +
  scale_color_gradient2(low="blue", high="red") +
  scale_size_continuous(range=c(0,5)) +
  geom_point() +
  facet_grid(database ~ ., space = "free", scales = "free")
ggsaveNF(outBase("MarkerOverlaps.pdf"),w=1, h=2.1)

# . Plot Izzo dataset original ----------------------------------------------
hex.obj <- hexbin::hexbin(x=izzo$UMAP_1, y=izzo$UMAP_2, xbins = 100, IDs=TRUE)
pDT <- cbind(izzo, data.table(hex.x=hex.obj@xcm, hex.y=hex.obj@ycm, hex.cell=hex.obj@cell)[match(hex.obj@cID, hex.cell),])
pDT <- pDT[,.N, by=c("hex.x", "hex.y", "functional.cluster")]
pDT[, sum := sum(N), by=c("hex.x", "hex.y")]
pDT[, frac := N / sum]
pDT <- pDT[frac > 0.25]
pDT.labels <- pDT[frac > 0.5, .(hex.x = median(hex.x), hex.y=median(hex.y)), by=c("functional.cluster")]
ggplot(pDT, aes(x=hex.x, y=hex.y)) +
  themeNF() +
  #geom_hex(fill="lightgrey", bins=100) +
  geom_point(aes(color=functional.cluster, alpha=frac), size=0.5) + 
  geom_text(data=pDT.labels, aes(label=functional.cluster), lineheight = 0.8) +
  #scale_color_manual(values=COLORS.CELLTYPES.scRNA.ainhoa) +
  xu + yu
ggsaveNF(outBase("UMAP_Izzo.pdf"), w=1.5,h=1.5)

ggplot(umap.proj$izzo[tissue != "Izzo_WT1"], aes(x=UMAP_1, UMAP_2)) + 
  themeNF() +
  stat_binhex(aes(fill=log10(..count..)), bins=100) + 
  scale_fill_gradientn(colors = c("lightgrey", "#1f78b4", "#e31a1c", "#ff7f00")) +
  facet_grid(. ~ tissue) +
  xu + yu
ggsaveNF(outBase("UMAP_Izzo_projected.pdf"), w=3,h=1.5)


# . UMAP Projections ---------------------------------------------------------
tx <- "in.vivo"
tx <- "leukemia"
tx <- "ex.vivo"
for(tx in names(inDir.funcs)){
  # out directory
  out <- dirout(paste0(base.dir, "/", tx))
  
  # Prepare data
  ann <- annList[tissue ==  tx]
  
  # Hexpoints
  (x <- names(umap.proj)[1])
  for(x in names(umap.proj)){
    xDT.ref <- if(x %in% c("in.vivo", "in.vivo.X")) umap.proj[[x]][match(annList[tissue ==  "in.vivo"]$rn, rn)]
    xDT <- umap.proj[[x]][match(ann$rn, rn)]
    hex.obj <- hexbin::hexbin(x=xDT$UMAP_1, y=xDT$UMAP_2, xbins = 100, IDs=TRUE)
    pDT <- cbind(ann, data.table(hex.x=hex.obj@xcm, hex.y=hex.obj@ycm, hex.cell=hex.obj@cell)[match(hex.obj@cID, hex.cell),])
    pDT <- pDT[,.(N=.N), by=c("hex.x", "hex.y", "Clusters")]
    pDT <- pDT[Clusters != "unclear"]
    pDT[, sum := sum(N), by=c("hex.x", "hex.y")]
    pDT[, frac := N / sum]
    pDT <- pDT[order(frac, decreasing = TRUE)][,head(.SD, 1), by=c("hex.x", "hex.y")]
    #pDT <- pDT[frac > 0.25]
    pDT.labels <- pDT[, .(hex.x = median(hex.x, na.rm=TRUE), hex.y=median(hex.y, na.rm=TRUE)), by=c("Clusters")]
    pDT[, Clusters := cleanCelltypes(Clusters,twoLines = FALSE)]
    pDT.labels[, Clusters := cleanCelltypes(Clusters,twoLines = TRUE)]
    p <- ggplot(pDT, aes(x=hex.x, y=hex.y)) +
      themeNF() + xu + yu + scale_color_manual(values=COLORS.CELLTYPES.scRNA.ainhoa)
    if(!is.null(xDT.ref)) p <- p + geom_hex(data=xDT.ref, fill="lightgrey", bins=100, aes(x=UMAP_1, y=UMAP_2))
    p <- p +
      geom_point(aes(color=Clusters), size=0.5) + 
      geom_text(data=pDT.labels, aes(label=Clusters), lineheight = 0.8)
    p
    ggsaveNF(out("UMAP_Celltypes_",x,".pdf"), w=1.5,h=1.5, plot=p)
  }
}


# . NTC Clusters depletion ------------------------------------------------------------
tx <- "in.vivo"
for(tx in names(inDir.funcs)){
  # out directory
  out <- dirout(paste0(base.dir, "/", tx))
  
  # Prepare data
  pDT.full <- annList[tissue ==  tx]
  
  # Calculate fraction
  pDT <- pDT.full[!is.na(gene)][, .(sum=.N, sumNTC = sum(!is.na(gene) & gene=="NTC")), by=c("Cluster.number", "Clusters", "timepoint")]
  pDT[, perc := sumNTC / sum * 100]
  pDT[, sumKO := sum-sumNTC]
  pDT[, Clusters := cleanCelltypes(Clusters)]
  ggplot(pDT, aes(x=sumKO+1, y=sumNTC+1, color=Clusters, shape=Clusters)) + 
    geom_point() +
    themeNF() +
    geom_abline() +
    scale_x_log10() +
    scale_y_log10() +
    facet_grid(. ~ timepoint) +
    scale_shape_manual(values=rep(c(1,16,2,18,3,4), 20)) +
    scale_color_manual(values=COLORS.CELLTYPES.scRNA.ainhoa) +
    geom_text_repel(data=pDT[sumKO / 10 > sumNTC], aes(label=Cluster.number), color="black") +
    xlab("Number of cells with KO guide + 1") + ylab("Number of cells with NTC guide + 1")
  ggsaveNF(out("NTC_depleted_clusters.pdf"),w=length(unique(pDT$timepoint)) * 1.5, h=1.5)
  
  # plot on UMAP
  ggplot(pDT.full, aes(x=UMAP1, y=UMAP2)) + 
    themeNF() +
    geom_hex(bins=100) +
    scale_fill_gradient(low="lightgrey", high="#ff7f00") +
    geom_text(data=pDT.full[,.(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by=c("Cluster.number")], aes(label=Cluster.number))
  ggsaveNF(out("NTC_depleted_clusters_UMAP.pdf"),w=1.5,h=1.5)
}


# # . Signatures --------------------------------------------------------------
# tx <- "in.vivo"
# tx <- "leukemia"
# tx <- "ex.vivo"
# for(tx in names(inDir.funcs)){
#   # out directory
#   out <- dirout(paste0(base.dir, "/", tx))
#   
#   # annotation
#   ann <- annList[tissue ==  tx & timepoint != "28d"]
#   
#   # Plot signatures scores in each cluster
#   pDT <- merge(
#     ann[, c("rn", "Clusters"),with=F],
#     marker.signatures[, c("rn", "DB", "value", "FinalName"),with=F],
#     by="rn")
#   pDT <- pDT[, mean(value), by=c("Clusters", "FinalName")]
#   pDT[, V1 := scale(V1), by="FinalName"]
#   pDT[, Clusters := cleanCelltypes(Clusters)]
#   #pDT <- hierarch.ordering(pDT, "Clusters", "FinalName", "V1")
#   pDT <- hierarch.ordering(pDT, "FinalName", "Clusters", "V1")
#   ggplot(pDT, aes(x=FinalName, y=Clusters, fill=V1)) + 
#     geom_tile() +
#     themeNF(rotate = TRUE) +
#     scale_fill_gradient2(low="blue", high="red") +
#     xlab("External gene signatures") +
#     ylab("Cell type\n(this study)")
#   ggsaveNF(out("Supp_Signatures.pdf"), w=1.2,h=length(unique(pDT$Clusters)) * 0.05 + 0.5)
#   
#   # Plot expression of marker genes in each cluster
#   
#   # Plot fraction of predicted cells in each cluster
# }


# # . Signatures DE with CRISPR KOs -------------------------------------------
# tx <-TISSUES[1]
# for(tx in TISSUES){
#   
#   out <- dirout(paste0(base.dir, "/", tx))
#   
#   pDT <- marker.signatures.DA[[tx]]
#   pDT <- merge(pDT, SIGS.USE.DA, by=c("sig"))
#   
#   # Summarize across guides
#   pDT[gene == "Pu.1", gene := "Spi1"]
#   pDT <- pDT[, .(
#     d=mean(d, na.rm=TRUE), 
#     dir=length(unique(sign(d[!is.na(d)])))==1, 
#     padj=sum(padj < 0.01), 
#     N=.N), by=c("sample", "FinalName", "gene")]
#   pDT[dir == FALSE, padj := 0]
#   pDT[dir == FALSE, d := NA]
#   
#   # summarize across samples
#   pDT <- pDT[, .(
#     d=mean(d, na.rm=TRUE), 
#     dir=length(unique(sign(d[!is.na(d)])))==1, 
#     padj=sum(padj), 
#     N=sum(N)), by=c("FinalName", "gene")]
#   pDT[dir == FALSE, padj := 0]
#   pDT[dir == FALSE, d := NA]
#   
#   # setup for plotting and plotting
#   pDT[padj == 0 | is.na(d), d := 0]
#   pDT[, sig.perc := padj / N]
#   pDT[,d_cap := pmin(abs(d), 5) * sign(d)]
#   pDT$FinalName <- factor(pDT$FinalName, levels=rev(SIGS.USE.DA$FinalName))
#   #pDT <- hierarch.ordering(dt = pDT, toOrder = "FinalName", orderBy = "gene", value.var = "d")
#   pDT <- hierarch.ordering(dt = pDT, toOrder = "gene", orderBy = "FinalName", value.var = "d")
#   ggplot(pDT, aes(x=gene,y=FinalName, color=d, size=sig.perc)) + 
#     themeNF(rotate = TRUE) +
#     scale_size_continuous(name="% sign.", range = c(0,4)) +
#     scale_color_gradient2(name="delta",low="#1f78b4", high="#e31a1c") +
#     geom_point() +
#     geom_point(shape=1, color="lightgrey") +
#     ylab("Marker signature") + xlab("")
#   ggsaveNF(out("MarkerSignatures_DA.pdf"),w=2,h=0.7)
#   
#   dla <- fread("metadata/FIGS_Order_Fig2_CFs.tsv")$Factor
#   pDT <- pDT[gene %in% dla]
#   pDT$gene <- factor(as.character(pDT$gene), levels=dla)
#   ggplot(pDT, aes(x=gene,y=FinalName, color=d, size=sig.perc)) + 
#     themeNF(rotate = TRUE) +
#     scale_size_continuous(name="% sign.", range = c(0,4)) +
#     scale_color_gradient2(name="delta",low="#1f78b4", high="#e31a1c") +
#     geom_point() +
#     geom_point(shape=1, color="lightgrey") +
#     ylab("Marker signature") + xlab("")
#   ggsaveNF(out("MarkerSignatures_DA_DLA.pdf"),w=1.5,h=0.7)
#   
#   # export table
#   write.tsv(pDT, out("MarkerSignatures_DA.tsv"))
# }


# . Marker gene expression --------------------------------------------------
pDT <- lapply(mobjs, function(obj) DotPlotData(cds = obj, markers = marker.genes, cols = "CellType"))
pDT <- rbindlist(pDT, idcol="tissue")
pDT <- pDT[tissue %in% c("ex.vivo", "in.vivo")]
pDT[, id := paste("tissue", "CellType")]
pDT <- pDT[Gene %in% pDT[,max(percentage), by="Gene"][V1 > 75]$Gene]
pDT[, scale := scale(mean), by=c("Gene", "tissue")]
pDT <- hierarch.ordering(pDT, "Gene", "id", "mean",aggregate = TRUE)
ggplot(pDT, aes(x=CellType, y=Gene, color=scale, size=percentage)) + 
  geom_point() +
  scale_size_continuous(range=c(0,5), limits = c(0,100)) +
  scale_color_gradient(low="lightgrey", high="blue") +
  facet_grid(. ~ tissue, space="free", scales = "free") +
  themeNF(rotate = TRUE)
ggsaveNF(outBase("Markers_Clusters.pdf"), w=2,h=2.5)


# . Cluster enrichment analyses ---------------------------------------------
tx <- "in.vivo"
inDir <- dirout_load("SCRNA_21_02_ClusterEnrichments_simple")
#for(tx in names(inDir.funcs)){
  # out directory
  out <- dirout(paste0(base.dir, "/", "cluster.enrichments/"))
  
  # Collect enrichment scores
  typex <- "basic_leukemia_noMixscape"
  for(typex in gsub("Guides_Fisher_Mixscape_(.+).pdf", "\\1", list.files(inDir(""), pattern="Guides_Fisher_Mixscape_.*.pdf"))){
    fish.file <- inDir("Guides_Fisher_Mixscape_",typex,".tsv")
    if(!file.exists(fish.file)) next
    
    fish.full <- fread(fish.file)
    fish.full[mixscape_class == "Pu.1", mixscape_class := "Spi1"]
    grep("S", unique(fish.full$mixscape_class), value = TRUE)
    #fish.full <- merge(fish.full, unique(SANN[,c("sample_broad", "timepoint"),with=F]), by.x="sample", by.y="sample_broad")
    timex <- "14d"
    for(timex in c(unique(fish.full$sample))){
      fish <- copy(fish.full)
      #fish <- fish[mixscape_class == "Rbbp4"]
      if(timex != "all") fish <- fish[sample == timex]
      
      # summarize across NTCs
      fish <- fish[, .(
        log2OR=mean(log2OR), 
        dir=length(unique(sign(log2OR[padj < 0.01]))) <= 1, 
        #dir=length(unique(sign(log2OR)))==1, 
        padj=sum(padj < 0.01),
        N=.N,
        Ncells=sum(unique(guide.cells))), by=c("sample", "Clusters", "mixscape_class")]
      fish[dir == FALSE, padj := 0]
      fish[dir == FALSE, log2OR := NA]
      
      # legacy
      fish[, gene := mixscape_class]
      
      # Summarize across guides
      # fish[, gene := gsub("_.+", "", mixscape_class)]
      # fish[gene == "Pu.1", gene := "Spi1"]
      # fish <- fish[, .(
      #   log2OR=mean(log2OR, na.rm=TRUE), 
      #   dir=length(unique(sign(log2OR[!is.na(log2OR)])))==1, 
      #   padj=sum(padj), 
      #   N=sum(N)), by=c("sample", "Clusters", "gene")]
      # fish[dir == FALSE, padj := 0]
      # fish[dir == FALSE, log2OR := NA]
      
      # summarize across samples
      # fish <- fish[, .(
      #   log2OR=mean(log2OR, na.rm=TRUE), 
      #   dir=length(unique(sign(log2OR[!is.na(log2OR)])))==1, 
      #   padj=sum(padj), 
      #   N=sum(N)), by=c("Clusters", "gene")]
      # fish[dir == FALSE, padj := 0]
      # fish[dir == FALSE, log2OR := NA]
      
      # setup for plotting
      fish[padj == 0 | is.na(log2OR), log2OR := 0]
      fish[, sig.perc := padj / N]
      fish[,log2OR_cap := pmin(abs(log2OR), 5) * sign(log2OR)]
      fish <- hierarch.ordering(dt = fish, toOrder = "gene", orderBy = "Clusters", value.var = "log2OR")
      fish[, Clusters := gsub("^Gran$", "Gran.", Clusters)]
      #fish[, Clusters := cleanCelltypes(Clusters)]
      #fish <- hierarch.ordering(dt = fish, toOrder = "Clusters", orderBy = "gene", value.var = "log2OR")
      ggplot(fish, aes(x=gene, y=Clusters, size=sig.perc, color=log2OR_cap)) + 
        themeNF(rotate=TRUE) +
        scale_color_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
        scale_size_continuous(name="% sign.", range = c(0,5)) +
        geom_point() +
        geom_point(shape=1, color="lightgrey") +
        xlab("Gene") + ylab("Cell type")
      ggsaveNF(
        out("Cluster_enrichments_",typex,"_", timex, ".pdf"), 
        w=length(unique(fish$gene))*0.05 + 0.5,
        h=length(unique(fish$Clusters))*0.05 + 0.5)
      write.tsv(fish, out("Cluster_enrichments_",typex,"_", timex,".tsv"))
    }
  }
#}


# . Cell numbers ----------------------------------------------------------
pDT <- annList[tissue != "leukemia"][timepoint != "28d"][,.N, by=c("Clusters", "tissue")]
pDT[, sum := sum(N), by="tissue"]
pDT[, perc := N/sum*100]
pDT[, Clusters := cleanCelltypes(Clusters, reverse = TRUE)]
ggplot(pDT, aes(x=Clusters, y=perc, fill=tissue)) + 
    geom_bar(stat="identity", position=position_dodge2(preserve = "single")) + 
    #scale_y_log10() +
    themeNF(rotate=TRUE) +
    ylab("Percent of cells") + xlab("")
ggsaveNF(outBase("CellCounts.pdf"), w=1.5,h=1, guides = TRUE)

pDT <- annList[tissue == "in.vivo"][timepoint != "28d"][,.N, by=c("Clusters", "tissue")]
pDT[, sum := sum(N), by="tissue"]
pDT[, Clusters := cleanCelltypes(Clusters, reverse = FALSE)]
ggplot(pDT, aes(x=Clusters, y=N)) + 
  geom_bar(stat="identity") + 
  #scale_y_log10() +
  themeNF(rotate=TRUE) +
  ylab("Number of cells") + xlab("") +
  scale_y_log10()
ggsaveNF(outBase("CellCounts_invivo.pdf"), w=1.5,h=1, guides = TRUE)

pDT <- annList[tissue != "leukemia"][timepoint != "28d"][,.N, by=c("gene", "Clusters", "tissue")][!is.na(gene)]
pDT[, Clusters := cleanCelltypes(Clusters, reverse = TRUE)]
ggplot(pDT,aes(x=Clusters, y=N, color=tissue)) + 
  themeNF(rotate=T) +
  geom_boxplot(coef=10^10, position=position_dodge2(preserve = "single")) + 
  ylab("Number of cells") + xlab("")
ggsaveNF(outBase("CellCounts_InVivoSparse.pdf"), w=1.5,h=1, guides = TRUE)



# . Cell numbers per guide --------------------------------------------------
pDT <- annList[tissue != "leukemia"][timepoint != "28d"][,.N, by=c("gene", "tissue")][!is.na(gene)]
pDT <- pDT[gene %in% dla.healthy$all]
pDT$gene <- factor(pDT$gene, levels=dla.healthy$all)
#pDT$gene <- factor(pDT$gene, levels=rev(pDT[,sum(N), by="gene"][order(V1)]$gene))
ggplot(pDT, aes(x=gene,y=N, fill=tissue)) + 
  themeNF(rotate = TRUE) +
  geom_bar(stat="identity", position=position_dodge2(width = 0.9, preserve = "single")) +
  scale_y_log10() +
  xlab("CF") + ylab("Number of cells (log10)")
ggsaveNF(outBase("CellCounts_CFs.pdf"), w=4,h=1, guides = TRUE)


# . NTC distributions to assess clonality effects -------------------------------------------------
pDT <- annList[mixscape_class.global == "NTC", .N, by=c("Clusters", "sample_broad", "CRISPR_Cellranger")]
pDT[, sum := sum(N), by=c("CRISPR_Cellranger", "sample_broad")]
pDT <- pDT[sum > 100]
pDT[, frac := N/sum * 100]
pDT <- merge(pDT, unique(SANN[,c("tissue", "timepoint","sample_broad"),with=F]), by="sample_broad")
pDT <- pDT[tissue != "leukemia"]
pDT[, Clusters := cleanCelltypes(Clusters, drop=TRUE)]
pDT[, guide_i := as.numeric(factor(CRISPR_Cellranger)), by=c("sample_broad")]
ggplot(pDT,  aes(x=Clusters, y=frac, fill=factor(guide_i))) + 
  theme_bw(12) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~sample_broad) +
  xRot()
ggsaveNF(outBase("NTC_clonality.pdf"), w=4, h=2)

pDT <- dcast.data.table(pDT, Clusters + tissue + sample_broad ~ paste0("Guide", guide_i), value.var = "frac")
pDT[, timepoint := gsub("^(.+)_(\\d+d)(\\_\\d)?$", "\\2", sample_broad)]
pDT[, batch := as.numeric(factor(sample_broad)), by=c("tissue", "timepoint")]
ggplot(pDT, aes(x=Guide1, y=Guide2, shape=factor(paste("Batch", batch)), color=factor(paste("Batch", batch)))) + 
  geom_point() +
  themeNF() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ tissue + timepoint, ncol = 2) + 
  scale_shape_manual(values=rep(c(1,16,2,18,3,4), 20))
ggsaveNF(outBase("NTC_clonalityScatter.pdf"), w=3, h=2, guides = TRUE)

# . Plot all guides ---------------------------------------------------------
tx <- "in.vivo"
for(tx in names(inDir.funcs)){
  # out directory
  out <- dirout(paste0(base.dir, "/", tx))
  
  pDT.top <- annList[tissue == tx][timepoint != "28d"]
  #pDT.top <- pDT.top[tissue != "in.vivo" | markers == "lin-"]
  pDT.ntc <- pDT.top[mixscape_class.global == "NTC"]
  pDT.final <- copy(pDT.ntc)
  pDT.final$plot <- "NTC"
  for(x in unique(pDT.top[!is.na(gene)]$gene)){
    # FIX MIXSCAPE
    #pDTg <- rbind(pDT.ntc, pDT.top[grepl(paste0("^", x), guide) & mixscape_class.global == "KO"])
    pDTg <- rbind(pDT.ntc, pDT.top[gene == x])
    pDTg$plot <- x
    pDT.final <- rbind(pDT.final, pDTg)
  }
  pDT.final$plot <- relevel(factor(pDT.final$plot),ref = "NTC")
  
  n.factors <- length(unique(pDT.final$plot))
  
  ggplot(pDT.final, aes(x=UMAP1, y=UMAP2)) + 
    themeNF(grid = FALSE) +
    geom_hex(data=pDT.final[mixscape_class.global == "NTC"], bins=100, fill="lightgrey") +
    geom_hex(data=pDT.final[mixscape_class.global != "NTC" | plot == "NTC"], bins=100) +
    scale_fill_gradientn(colours=c("#1f78b4", "#e31a1c")) +
    facet_wrap(~plot, ncol=6) + 
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(), 
      panel.border = element_blank(),
      strip.background = element_blank(), 
      plot.background = element_blank()) +
    xu + yu
  ggsaveNF(out("UMAP_Guides_all.pdf"), w=4,h=ceiling(n.factors/6)*4/6)
  
  # pDT.final <- cbind(pDT.final, umap.proj[["in.vivo.X"]][match(pDT.final$rn, rn)][,c("UMAP_1", "UMAP_2"),with=F])
  # ggplot(pDT.final, aes(x=UMAP_1, y=UMAP_2)) + 
  #   themeNF() +
  #   geom_hex(data=pDT.final[mixscape_class.global == "NTC"], bins=100, fill="lightgrey") +
  #   geom_hex(data=pDT.final[mixscape_class.global != "NTC" | plot == "NTC"], bins=100) +
  #   scale_fill_gradientn(colours=c("#1f78b4", "#e31a1c")) +
  #   facet_wrap(~plot, ncol=6) + 
  #   xu + yu
  # ggsaveNF(out("UMAP_Guides_all_Crossprojected.pdf"), w=4,h=4)
}



# IN VIVO ---------------------------------------------------------
inDir.current <- "in.vivo"
out <- dirout(paste0(base.dir, "/", inDir.current))
ann <- annList[tissue == inDir.current]
# FIX MIXSCAPE
list.files(dirout_load(base.dir)("cluster.enrichments"))
fish.numeric <- fread(dirout_load(base.dir)("cluster.enrichments/Cluster_enrichments_numeric_in.vivo_noMixscape_14d.tsv"))
fish.enrich.broad <- fread(dirout_load(base.dir)("cluster.enrichments/Cluster_enrichments_DavidSpecial_in.vivo_noMixscape_14d",".tsv"))[Clusters != "unclear"]
fish.EryVsMye <- fread(dirout_load(base.dir)("cluster.enrichments/Cluster_enrichments","_eryVsMye", "_in.vivo", "_noMixscape", "_14d",".tsv"))
fish.early <- fread(dirout_load(base.dir)("cluster.enrichments/Cluster_enrichments","_earlyBranches", "_in.vivo", "_noMixscape", "_14d",".tsv"))
fish.monoVsGran <- fread(dirout_load(base.dir)("cluster.enrichments/Cluster_enrichments","_monoVsGran", "_in.vivo", "_noMixscape", "_14d",".tsv"))
fish.d28 <- fread(dirout_load(base.dir)("cluster.enrichments/Cluster_enrichments_DavidSpecial_in.vivo_noMixscape_28d",".tsv"))[Clusters != "unclear"]

viability <- fread(dirout_load("INT_04_Viability")("scRNA_processed.tsv"))
viability <- viability[variable == "scRNA in.vivo 14d"]

# . supp table ------------------------------------------------------------
exDT <- fish.enrich.broad[,c("Clusters", "gene", "log2OR", "padj", "N", "sig.perc"),with=F]
exDT[, sig.perc := sig.perc * 100]
colnames(exDT) <- c("cell type", "CF-KO", "log2 odds ratio", "significant guides (count)", "total guides (count)", "significant guides (percent)")
WriteXLS(x=exDT, ExcelFileName=out("Supplementary_Table_CellTypes_invivo.xls"), AdjWidth=TRUE, BoldHeaderRow=TRUE, FreezeRow=1, SheetNames="Table")
write.tsv(exDT, out("Supplementary_Table_CellTypes_invivo.tsv"))



# . UMAP of timepoint / markers ---------------------------------------------------------
ggplot(ann[gene == "NTC"], aes(x=UMAP1, y=UMAP2)) + 
  themeNF() +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~markers + timepoint, ncol=2) +
  xu + yu
ggsaveNF(out("UMAP_Groups.pdf"), w=1.4,h=1.4, guides=TRUE)

# . UMAP of samples -------------------------------------------------------
ggplot(ann[gene == "NTC"], aes(x=UMAP1, y=UMAP2)) + 
  themeNF() +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~sample_broad, ncol=3) +
  xu + yu
ggsaveNF(out("UMAP_Samples.pdf"), w=2,h=2)


# . Plot manual enrichments -------------------------------------------------
typex <- "all"
for(typex in names(dla.healthy)){
  dla <- dla.healthy[[typex]]
  h=length(unique(dla)) * 0.07 + 0.1
  
  # Plot enrichments with DLA order
  pDT <- fish.enrich.broad[gene %in% dla]
  pDT$gene <- factor(pDT$gene, levels = rev(dla))
  pDT <- pDT[!(Clusters == "MEP (pert.)" & log2OR < 0)]
  pDT[Clusters == "Mega", Clusters := "MkP"]
  pDT[, Clusters := cleanCelltypes(Clusters, clean=TRUE, reverse=FALSE)]
  (p_broad <- ggplot(pDT, aes(y=gene, x=Clusters, size=sig.perc, color=log2OR_cap)) + 
    themeNF(rotate=T) +
    scale_color_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
    scale_size_continuous(name="% sign.", range = c(0,5)) +
    scale_y_discrete(position="right") +
    geom_point() +
    geom_point(shape=1, color="lightgrey") +
    xlab("Lineage") + ylab(""))
  ggsaveNF(out("ClusterEnrichments_manual_cleaned_",typex,".pdf"), h=h,w=0.9, guides = TRUE)
  
  # Plot mye vs GMP enrichments with DLA order
  pDT <- fish.EryVsMye[Clusters == "GMP"][gene %in% dla]#[log2OR > 0 | log2OR < -2]
  #pDT$gene <- factor(pDT$gene, levels = pDT[order(log2OR)]$gene)
  pDT$gene <- factor(pDT$gene, levels = rev(dla))
  (p_EryVsMye <- ggplot(pDT, aes(x=log2OR, y=gene, alpha=sig.perc, fill=sign(log2OR_cap))) + 
    themeNF(rotate=F) +
    geom_col() +
    scale_fill_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
    #scale_size_continuous(name="% sign.", range = c(0,5), limits = c(0,1)) +
    #geom_point() +
    scale_y_discrete(position="right") +
    ylab("") + xlab("Myeloid enrichment\nvs erythroid (log2OR)") +
    geom_vline(xintercept = 0))
  ggsaveNF(out("ClusterEnrichments_manual_GMPvsMEP_",typex,".pdf"), w=1,h=h, guides = TRUE)
  
  pDT <- fish.early[Clusters == "GMP"][gene %in% dla]#[log2OR > 0 | log2OR < -2]
  #pDT$gene <- factor(pDT$gene, levels = pDT[order(log2OR)]$gene)
  pDT$gene <- factor(pDT$gene, levels = rev(dla))
  (p_Early <- ggplot(pDT, aes(x=log2OR, y=gene, alpha=sig.perc, fill=sign(log2OR_cap))) + 
    themeNF(rotate=F) +
    geom_col() +
    scale_fill_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
    #scale_size_continuous(name="% sign.", range = c(0,5), limits = c(0,1)) +
    #geom_point() +
    scale_y_discrete(position="right") +
    ylab("") + xlab("Early myeloid enrichment\nvs erythroid (log2OR)") +
    geom_vline(xintercept = 0))
  ggsaveNF(out("ClusterEnrichments_manual_EarlyGMPvsMEP_",typex,".pdf"), w=1,h=h, guides = TRUE)
  
  pDT <- fish.monoVsGran[Clusters == "Mono"][gene %in% dla]#[log2OR > 0 | log2OR < -2]
  #pDT$gene <- factor(pDT$gene, levels = pDT[order(log2OR)]$gene)
  pDT$gene <- factor(pDT$gene, levels = rev(dla))
  (p_MonoVsGran <- ggplot(pDT, aes(x=log2OR, y=gene, alpha=sig.perc, fill=sign(log2OR_cap))) + 
    themeNF(rotate=F) +
    geom_col() +
    scale_fill_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
    #scale_size_continuous(name="% sign.", range = c(0,5), limits = c(0,1)) +
    #geom_point() +
    scale_y_discrete(position="right") +
    ylab("") + xlab("Monocyte enrichment\nvs granulocyte (log2OR)") +
    geom_vline(xintercept = 0))
  ggsaveNF(out("ClusterEnrichments_manual_MonoVsGran_",typex,".pdf"), w=1,h=h, guides = TRUE)
  
  # Viability
  pDT <- copy(viability)
  pDT[, rn := gsub("^(.)", "\\U\\1", tolower(pDT$rn), perl=TRUE)]
  pDT <- pDT[rn %in% dla]
  pDT$gene <- factor(pDT$rn, levels = rev(dla))
  (p_viability <- ggplot(pDT, aes(x=value, y=gene, fill=sign(value))) + 
      themeNF(rotate=F) +
      geom_col() +
      scale_fill_gradient2(name="log2 fold change",low="blue", midpoint = 0, high="red") +
      #scale_size_continuous(name="% sign.", range = c(0,5), limits = c(0,1)) +
      #geom_point() +
      scale_y_discrete(position="right") +
      ylab("") + xlab("Viability (log2 fold change)") +
      geom_vline(xintercept = 0))
  ggsaveNF(out("ClusterEnrichments_manual_Viability_",typex,".pdf"), w=1,h=h, guides = TRUE)
  
  p_final <- 
    (p_broad + guides(fill=FALSE, color=FALSE, size=FALSE, alpha=FALSE)) +
    (p_Early + guides(fill=FALSE, color=FALSE, size=FALSE, alpha=FALSE)) +
    (p_EryVsMye+ guides(fill=FALSE, color=FALSE, size=FALSE, alpha=FALSE)) +
    (p_MonoVsGran + guides(fill=FALSE, color=FALSE, size=FALSE, alpha=FALSE)) +
    (p_viability + guides(fill=FALSE, color=FALSE, size=FALSE, alpha=FALSE)) +
    plot_layout(ncol=5)
  
  ggsaveNF(out("ClusterEnrichments_manual_Combined_",typex,".pdf"), w=5,h=h, plot=p_final, guides = TRUE)
}

# . Special clusters -----------------------------------------------------------------------
cx <- 26
for(cx in c(26, 45, 51)){
  (title <- paste0(
    "NTCs: ",
    round(nrow(ann[Cluster.number == cx & gene == "NTC"])/nrow(ann[Cluster.number == cx & gene != "NTC"]) * 100,2),
    "% of ",
    nrow(ann[Cluster.number == cx & !is.na(gene)]),
    " cells"
    ))
  pDT <- fish.numeric[Clusters == paste("cl", cx)]
  pDT[, gene := factor(gene, levels=rev(pDT[order(log2OR)]$gene))]
  pDT <- pDT[sig.perc > 0.5 & log2OR > 0]
  ggplot(pDT, aes(x=gene, y=log2OR)) + 
    geom_col() +
    themeNF(rotate = TRUE) +
    labs(y="log2(OR)", x="") +
    ggtitle(title)
  ggsaveNF(out("ClusterEnrichment_",cx,".pdf"), w=length(unique(pDT$gene)) * 0.05 + 0.5,h=0.7, guides=TRUE)
}


# . Characterize cluster 26 (MEP Pert) --------------------------------------

# enrichments
gsea <- fread(dirout("SCRNA_22_01_SpecialClustersDE")("GSEA_significant.tsv"))
dbs <- c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "TRANSFAC_and_JASPAR_PWMs", "TRRUST_Transcription_Factors_2019", "KEGG_2019_Mouse", "MSigDB_Hallmark_2020","WikiPathways_2019_Mouse")
gsea <- gsea[db %in% dbs][cluster == "Markers_26.tsv"][ES > 0]
gsea[, db := factor(db, levels=dbs)]
gsea[, pathway := factor(pathway, levels = unique(gsea[order(ES)]$pathway))]
ggplot(gsea, aes(x=ES, y=pathway, size=-log10(padj))) + 
  geom_point() + 
  facet_grid(db ~ . , space = "free", scales = "free") +
  themeNF() +
  scale_x_continuous(limits=c(0, NA))
ggsaveNF(out("SpecialCluster_26_enrichments.pdf"), w=2.5,h=1.8, guides=FALSE)

# marker genes
cl <- c(28, 16, 26, 29, 11, 38, 45, 51, 1, 25, 30, 54, 49)
gg <- c("Hoxa9", "Gata2", "Gfi1b", "Hba-a1", "Hba-a2", "Hbb-bs", "Mcub", "Mgst1", "Phf13", "Kcne3", "Hmgn3", "Cebpa", "Cebpe", "Irf8", "Spi1", "Thbs1", "Thbs4", "Pecam1", "Dmwd", "Ddit4", "Itgb7", "Lat2")
names(gg) <- gg
cds <- mobjs[["in.vivo"]]
stopifnot(all(gg %in% row.names(cds)))
cds$Cluster.number <- monocle3::clusters(cds)
pDT <- DotPlotData(cds = cds, markers = gg, cols = "Cluster.number")
pDT[, scale := scale(mean), by=c("Gene")]
pDT$Gene <- factor(pDT$Gene, levels=gg)
pDT[,Cluster.number := as.numeric(Cluster.number)]
pDT <- pDT[Cluster.number %in% cl]
pDT[, Cluster.number.fac := factor(as.character(Cluster.number), levels=rev(as.character(cl)))]
h=length(unique(pDT$Cluster.number)) * 0.08 + 0.2
w=length(unique(pDT$Gene)) * 0.08 + 0.2
ggplot(pDT, aes(y=Cluster.number.fac, x=Gene, color=scale, size=percentage)) + 
  geom_point() +
  #geom_point(color="black", shape=1) + 
  scale_size_continuous(range=c(0,5), limits = c(0,100)) +
  scale_color_gradientn(colours = c("lightgrey", "grey", "#1f78b4", "black")) +
  facet_grid(. ~ ., space="free", scales = "free") +
  themeNF(rotate = TRUE) +
  ylab("Cluster") + xlab("Gene")
ggsaveNF(out("SpecialCluster_26_markers.pdf"), w=w,h=h)


# . Plot displasia (d28) mye vs GMP EARLY enrichments -----------------------------------------------------------------------
gg <- c("Brd9", "Smarcd1", "Smarcd2")
# Plot enrichments with DLA order
pDT <- fish.enrich.broad[gene %in% gg]
pDT$gene <- factor(pDT$gene, levels = rev(gg))
pDT <- pDT[!(Clusters == "MEP (pert.)" & log2OR < 0)]
pDT[Clusters == "Mega", Clusters := "MkP"]
pDT[, Clusters := cleanCelltypes(Clusters, clean=TRUE, reverse=FALSE)]
h=length(unique(gg)) * 0.07 + 0.5
(p <- ggplot(pDT, aes(y=gene, x=Clusters, size=sig.perc, color=log2OR_cap)) + 
    themeNF(rotate=T) +
    scale_color_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
    scale_size_continuous(name="% sign.", range = c(0,5)) +
    scale_y_discrete(position="right") +
    geom_point() +
    geom_point(shape=1, color="lightgrey") +
    xlab("Lineage") + ylab(""))
ggsaveNF(out("ClusterEnrichments_manual_d28.pdf"), h=h,w=1.5, guides = TRUE, plot=p)


# . Plot displasia (d28) guides ---------------------------------------------------------
gg <- c("Brd9")
pDT.top <- ann[timepoint == "28d"][grepl("OP1", sample)][gene %in% c("NTC", gg)]
#pDT.ntc <- pDT.top[mixscape_class.global == "NTC"]
pDT.ntc <- ann[timepoint == "28d"]
pDT.ntc[, type := "Background"]
pDT.final <- data.table()
#for(x in c("Kmt2a","Kmt2d","Smarcd2","Smarcd1","Brd9")){
for(x in unique(pDT.top$gene)){
  pDTg <- rbind(pDT.ntc, pDT.top[gene == x], fill=TRUE)
  pDTg$plot <- x
  pDT.final <- rbind(pDT.final, pDTg, fill=TRUE)
}
#pDT.final$plot <- relevel(factor(pDT.final$plot), ref = "NTC")
ggplot(pDT.final, aes(x=UMAP1, y=UMAP2)) + 
  themeNF() +
  geom_hex(data=pDT.final[type == "Background"], bins=100, fill="lightgrey") +
  geom_hex(data=pDT.final[is.na(type)], bins=100) +
  scale_fill_gradientn(colours=c("#1f78b4", "#e31a1c")) +
  facet_wrap(~plot, ncol=3) +
  xu + yu
ggsaveNF(out("UMAP_Guides_displasia.pdf"), w=2,h=1)


# . Plot displasia (d28) numbers -----------------------------------------------
pDT <- copy(ann[markers=="lin-"][grepl("OP1", sample)][mixscape_class.global == "NTC" | gene %in% c("Kmt2a","Kmt2d","Smarcd2","Smarcd1","Brd9", "NTC")])
pDT[, group := gsub("_.+", "", guide)]
pDT <- pDT[!is.na(group)]
pDT <- pDT[,.N, by=c("timepoint", "group")]
pDT[, sum := sum(N), by="timepoint"]
pDT[, perc := N/sum*100]
#pDT <- pDT[mixscape_class.global != "NP"]
#pDT <- pDT[group %in% pDT[, .N, by="group"][N == 2]$group]

# Bar plot
pDT$group <- factor(pDT$group, levels=unique(rev(pDT[order(timepoint, N)]$group)))
ggplot(pDT, aes(x=group, y=perc, fill=group == "NTC")) + 
  themeNF() +
  scale_fill_manual(values=c("black", "red")) +
  geom_bar(position="dodge", stat='identity') +
  facet_grid(. ~ timepoint, space="free", scales = "free") + 
  guides(fill=FALSE) +
  ylab("Percent of cells with guide") + 
  xRot() + 
  xlab("")
ggsaveNF(out("Displasia_Numbers_bars.pdf"), w=2,h=1, guides = TRUE)



# EX VIVO---------------------------------------------------------

# . load data ---------------------------------------------------------
inDir.current <- "ex.vivo"
out <- dirout(paste0(base.dir, "/", inDir.current))
ann <- annList[tissue == inDir.current]
marker.signatures[DB == "Larry"]
# FIX MIXSCAPE
fish.enrich <- list(
  day7=fread(dirout_load(base.dir)("cluster.enrichments/Cluster_enrichments_basic_ex.vivo_noMixscape_7d.tsv")),
  day9=fread(dirout_load(base.dir)("cluster.enrichments/Cluster_enrichments_basic_ex.vivo_noMixscape_9d.tsv"))
  )
fish.enrich <- rbindlist(fish.enrich, idcol="day")


# . supp table ------------------------------------------------------------
exDT <- fish.enrich[,c("day", "Clusters", "gene", "log2OR", "padj", "N", "sig.perc"),with=F]
exDT[, sig.perc := sig.perc * 100]
colnames(exDT) <- c("timepoint", "cell type", "CF-KO", "log2 odds ratio", "significant guides (count)", "total guides (count)", "significant guides (percent)")
WriteXLS(x=exDT, ExcelFileName=out("Supplementary_Table_CellTypes_exvivo.xls"), AdjWidth=TRUE, BoldHeaderRow=TRUE, FreezeRow=1, SheetNames="Table")
write.tsv(exDT, out("Supplementary_Table_CellTypes_exvivo.tsv"))


# . BULK Signatures ---------------------------------------------------------
pDT <- merge(umap.proj$izzo, marker.signatures.bulk, by="rn")
pDT[, value.norm := scale(value), by="variable"]
ggplot(pDT, aes(x=UMAP_1, y=UMAP_2)) +
  stat_summary_hex(aes(z=pmin(3, value.norm)),fun=mean, bins=100) +
  scale_fill_gradient2(low="blue", midpoint = 0, high="red") +
  themeNF() +
  facet_wrap(~variable, ncol = 3) +
  xu + yu
ggsaveNF(out("BulkSignatures.pdf"), w=2.5,h=1)

pDT <- merge(umap.proj$izzo, marker.signatures.bulk2[variable != "LSK"], by="rn")
pDT[, value.norm := scale(value), by="variable"]
ggplot(pDT, aes(x=UMAP_1, y=UMAP_2)) +
  stat_summary_hex(aes(z=pmin(3, value.norm)),fun=mean, bins=100) +
  scale_fill_gradient2(low="blue", midpoint = 0, high="red") +
  themeNF() +
  facet_wrap(~variable, ncol = 3) +
  xu + yu
ggsaveNF(out("BulkSignatures_Supp.pdf"), w=2.5,h=1)


# . Plot enrichments with good order for day 7
suppx <- "main"
for(suppx in names(dla.healthy)){
  dla <- dla.healthy[[suppx]]
  pDT <- fish.enrich[gene %in% dla]
  pDT <- pDT[Clusters %in% unique(ann$Clusters)]
  pDT <- pDT[(Clusters == "GMP (late)" & day == "day9") | (Clusters != "GMP (late)" & day == "day7")]
  with(pDT, table(Clusters, day))
  #pDT$Clusters <- factor(pDT$Clusters, levels = c("HSC", "EBMP", "GMP", "GMP (late)", "MkP", "Eo/Ba"))
  pDT$gene <- factor(pDT$gene, levels = rev(dla))
  pDT[, celltype := cleanCelltypes(Clusters, clean=TRUE, twoLines = FALSE, order = TRUE, reverse = FALSE)]
  h=length(unique(dla)) * 0.03 + 0.2
  ggplot(pDT, aes(y=gene, x=celltype, size=sig.perc, color=log2OR_cap)) + 
    themeNF(rotate=TRUE) +
    scale_color_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
    scale_size_continuous(name="% sign.", range = c(0,4)) +
    #facet_grid(. ~ day, space = "free", scales = "free") +
    geom_point() +
    geom_point(shape=1, color="lightgrey") +
    scale_y_discrete(position="right") +
    xlab("Gene") + ylab("Cell type") +
    theme(strip.text.y = element_text(angle=0))
  ggsaveNF(out("ClusterEnrichments_manual_combined7and9_",suppx,".pdf"), w=0.7,h=h)
}



# LEUKEMIA---------------------------------------------------------

# . Define clustesr -------------------------------------------------------
clusters.plot.cts <- list()

clusters.plot.cts$main <- list(
  "Gran"=24,
  "Gran Prog."=c(21, 18, 25),
  "Eo/Ba"=11,
  "Ery Prog."=c(20, 5), 
  "MkP"=26
)

clusters.plot.cts$supp <- c(
  list(
    "unclear"=c(22, 8, 7),
    "HSC"=c(10, 17, 9, 23, 12, 14, 15),
    "GMP"=c(19, 16, 3, 13, 6, 1, 2),
    "Mono"=c(28, 27,4)
    ),
  clusters.plot.cts$main
)

stopifnot(all(table(unlist(clusters.plot.cts$supp)) == 1))


# . load data ---------------------------------------------------------
inDir.current <- "leukemia"
out <- dirout(paste0(base.dir, "/", inDir.current))
ann <- annList[tissue == inDir.current]
# Use in vivo projected UMAP?
#ann <- merge(ann[,-c("UMAP1", "UMAP2"),with=F], setNames(umap.proj[["in.vivo.X"]][,c("rn", "UMAP_1", "UMAP_2")], c("rn", "UMAP1", "UMAP2")), by="rn")
abs <- fread(inDir.funcs[[inDir.current]]("Antibodies.tsv"))
abs <- merge(abs[,-c("UMAP1", "UMAP2"),with=F], setNames(umap.proj[["in.vivo.X"]][,c("rn", "UMAP_1", "UMAP_2")], c("rn", "UMAP1", "UMAP2")), by="rn")
abs$Clusters <- ann[match(abs$rn, rn)]$Clusters
abs$Cluster.number <- ann[match(abs$rn, rn)]$Cluster.number
abs.main.dla <- toupper(c("CD34","CD41","CD55","FCeR1","CD11b","Ly6C"))

#clusters.plot <- c(8,16,12,5,11,15,6,7,18)

# update ann with cluster numbers
ann[, Clusters := "other"]
for(xnam in names(clusters.plot.cts$main)){
  ann[Cluster.number %in% clusters.plot.cts$main[[xnam]], Clusters := xnam]
}

abs[, Clusters := "other"]
for(xnam in names(clusters.plot.cts$main)){
  abs[Cluster.number %in% clusters.plot.cts$main[[xnam]], Clusters := xnam]
}


# FIX MIXSCAPE
fish.enrich <- fread(dirout_load(base.dir)("cluster.enrichments/Cluster_enrichments_numeric_leukemia_noMixscape_6d.tsv"))
#fish.enrich$celltype <- unique(ann[,c("Cluster.number", "Clusters"),with=F])[match(fish.enrich$Clusters, paste("cl", Cluster.number))]$Clusters


# . supp table ------------------------------------------------------------
exDT <- fish.enrich[,c("Clusters", "gene", "log2OR", "padj", "N", "sig.perc"),with=F]
exDT[, sig.perc := sig.perc * 100]
colnames(exDT) <- c("cluster", "CF-KO", "log2 odds ratio", "significant guides (count)", "total guides (count)", "significant guides (percent)")
WriteXLS(x=exDT, ExcelFileName=out("Supplementary_Table_Clusters_leukemia.xls"), AdjWidth=TRUE, BoldHeaderRow=TRUE, FreezeRow=1, SheetNames="Table")
write.tsv(exDT, out("Supplementary_Table_Clusters_leukemia.tsv"))

# . Clusters on UMAP --------------------------------------------------------
typex <- "original"
for(typex in c("original", "in.vivo.X")){
  xDT <- umap.proj[[typex]][match(ann$rn, rn)]
  hex.obj <- hexbin::hexbin(x=xDT$UMAP_1, y=xDT$UMAP_2, xbins = 100, IDs=TRUE)
  pDT <- cbind(ann, data.table(hex.x=hex.obj@xcm, hex.y=hex.obj@ycm, hex.cell=hex.obj@cell)[match(hex.obj@cID, hex.cell),])
  pDT <- pDT[,.(N=.N), by=c("hex.x", "hex.y", "Cluster.number", "Clusters")]
  pDT[, sum := sum(N), by=c("hex.x", "hex.y")]
  pDT[, frac := N / sum]
  pDT <- pDT[order(frac, decreasing = TRUE)][,head(.SD, 1), by=c("hex.x", "hex.y")]
  cols <- COLORS.CELLTYPES.scRNA.ainhoa
  cols["other"] <- "grey"
  pDT[Cluster.number %in% clusters.plot.cts$supp$HSC, Clusters := "LSK"]
  pDT.labels <- pDT[, .(hex.x = median(hex.x), hex.y=median(hex.y)), by=c("Cluster.number")]
  #pDT[, Clusters := cleanCelltypes(Clusters)]
  #pDT[!Cluster.number %in% clusters.plot, Clusters := "other"]
  p <- ggplot(pDT, aes(x=hex.x, y=hex.y)) +
    themeNF() + xu + yu +
    geom_point(aes(color=Clusters), size=0.5) + 
    scale_color_manual(values=cols) +
    geom_text(data=pDT.labels, aes(label=Cluster.number), lineheight = 0.8)
  ggsaveNF(out("UMAP_numeric_",typex,"_main.pdf"), w=1.5,h=1.5, plot=p)
}


# . Antibodies on UMAP ----------------------------------------------------
ggplot(abs, aes(x=UMAP1, y=UMAP2)) +
  stat_summary_hex(bins = 100, aes(z=pmin(abs(Signal.norm), 2) * sign(Signal.norm)),fun=mean) +
  scale_fill_gradient2(low="blue", high="red") +
  facet_wrap(~Antibody) +
  themeNF() +
  xu + yu
ggsaveNF(out("Antibodies_UMAP.pdf"), w=2, h=2)


# . Summarize Antibodies in clusters --------------------------------------
pDT <- abs[, .(Signal=mean(Signal.norm, na.rm=TRUE)), by=c("Cluster.number", "Antibody")]
suppx <- "supp"
for(suppx in c("main", "supp")){
  cl.labels <- clusters.plot.cts[[suppx]]
  
  pDTx <- pDT[Cluster.number %in% unlist(cl.labels)]
  pDTx[, Cluster.number := as.numeric(Cluster.number)]
  
  pDTx <- merge(pDTx, melt(cl.labels), by.x="Cluster.number", by.y="value")
  pDTx[, celltype := factor(L1, levels=rev(names(cl.labels)))]
  
  if(suppx == "main"){
    pDTx <- pDTx[Antibody %in% abs.main.dla]
    pDTx$Antibody <- factor(pDTx$Antibody, levels = abs.main.dla)
  }
  pDTx[, Cluster.number.fac := factor(as.character(Cluster.number), levels=rev(as.character(unlist(cl.labels))))]
  
  ggplot(pDTx, aes(y=Cluster.number.fac,x=Antibody, fill=Signal)) +
    themeNF() +
    geom_tile() +
    scale_fill_gradient2(low="blue", high="red") +
    xRot() +
    facet_grid(celltype ~ ., space = "free", scales = "free") +
    ylab("Clusters")
  ggsaveNF(out("Antibodies_Average_",suppx,".pdf"), 
           w=length(unique(pDTx$Antibody)) * 0.05 + 0.3, 
           h=length(unique(pDTx$Cluster.number)) * 0.05 + 0.3)
}


# . Cell cycle ------------------------------------------------------------
# UMAP
ggplot(ann, aes(x=UMAP1, y=UMAP2)) + 
  themeNF() +
  geom_hex(bins=100) +
  scale_fill_hexbin() +
  facet_grid(. ~ Phase) +
  xu + yu
ggsaveNF(out("CellCycle_UMAP.pdf"), w=2,h=1)

# in one UMAP
hex.obj <- hexbin::hexbin(x=ann$UMAP1, y=ann$UMAP2, xbins = 100, IDs=TRUE)
pDT <- cbind(ann, data.table(hex.x=hex.obj@xcm, hex.y=hex.obj@ycm, hex.cell=hex.obj@cell)[match(hex.obj@cID, hex.cell),])
pDT <- pDT[,.N, by=c("hex.x", "hex.y", "Phase")]
pDT[, sum := sum(N), by=c("hex.x", "hex.y")]
pDT[, frac := N/sum]
pDT <- dcast.data.table(pDT, hex.x + hex.y ~ Phase, value.var = "frac")
pDT.cols <- as.matrix(pDT[,c("G1", "G2M", "S")])
pDT.cols[is.na(pDT.cols)] <- 0
pDT$col <- apply(pDT.cols %*% rbind(c(255, 161, 77), c(65,255,78), c(120, 77, 255)), 1, function(x){rgb(x[1], x[2], x[3], alpha=255, maxColorValue = 255)})
pDT$id <- paste0("x", 1:nrow(pDT))
ggplot(pDT, aes(x=hex.x, y=hex.y, color=id)) + 
  themeNF() +
  guides(color=FALSE) +
  geom_point(size=0.1) +
  scale_color_manual(values=setNames(pDT$col, pDT$id)) +
  geom_text(x=-7.5, y=7.5, label="G1", color="#FFA14D") +
  geom_text(x=-7.5, y=6.5, label="G2M", color="#41FF4E") +
  geom_text(x=-7.5, y=5.5, label="S", color="#784DFF") +
  xu + yu
ggsaveNF(out("CellCycle_UMAP_single.pdf"), w=1,h=1)


# Percentage
pDT <- ann[!is.na(gene)][,.N, by=c("Phase", "CRISPR_Cellranger", "timepoint", "markers")]
pDT[,sum := sum(N), by=c("CRISPR_Cellranger", "timepoint", "markers")]
pDT[, perc := N/sum*100]
pDT[, gene := gsub("_.+", "", CRISPR_Cellranger)]
pDT$gene <- factor(pDT$gene, levels=pDT[,mean(perc), by=c("Phase", "gene")][order(V1)][Phase == "G1"]$gene)
ggplot(pDT, aes(x=CRISPR_Cellranger,y=perc,fill=Phase)) + 
  themeNF(rotate = TRUE) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c(G1 = "#b2df8a", G2M = "#6a3d9a", S = "grey")) +
  facet_grid(. ~ gene, scales = "free", space = "free", switch = "x") +
  theme(strip.text.x = element_text(angle=90)) +
  theme(panel.spacing = unit(0.05, "cm")) +
  ylab("Percentage of cells") + xlab("")
ggsaveNF(out("CellCycle_Numbers.pdf"), w=3,h=1.2)


# . manual enrichments ----------------------------------------------------
pDT <- copy(fish.enrich)
#pDT$Complex <- dla[match(pDT$gene, Factor)]$Complex
#pDT[, celltype := cleanCelltypes(celltype, clean=TRUE, twoLines = FALSE, order = TRUE, reverse = TRUE)]
pDT[, log2OR_cap := pmin(2, abs(log2OR)) * sign(log2OR)]
suppx <- "supp"
for(suppx in names(dla.cancer)){
  dla <- dla.cancer[[suppx]]
  cl.labels <- clusters.plot.cts[[suppx]]
  
  pDTx <- pDT[Clusters %in% paste("cl", unlist(cl.labels))]
  pDTx <- pDTx[gene %in% dla]
  pDTx$gene <- factor(pDTx$gene, levels = dla)
  pDTx[,Clusters := as.numeric(gsub("cl ", "", Clusters))]
  
  pDTx <- merge(pDTx, melt(cl.labels), by.x="Clusters", by.y="value")
  pDTx[, celltype := factor(L1, levels=rev(names(cl.labels)))]
  pDTx <- pDTx[!(Clusters == 24 & log2OR < 0)]
  
  pDTx[, Cluster.number.fac := factor(as.character(Clusters), levels=rev(as.character(unlist(cl.labels))))]
  
  w=length(unique(pDTx$gene)) * 0.04 + 1
  h=length(unique(pDTx$Clusters)) * 0.06 + 0.3
  ggplot(pDTx, aes(x=gene, y=Cluster.number.fac, size=sig.perc*100, color=log2OR_cap)) + 
    themeNF(rotate=TRUE) +
    scale_color_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
    scale_size_continuous(name="% sign.", range = c(0,4)) +
    geom_point() +
    scale_y_discrete(position="right") +
    geom_point(shape=1, color="lightgrey") +
    xlab("Gene") + ylab("Cluster / cell type") +
    facet_grid(celltype ~ . , space="free", scales = "free") +
    theme(strip.text.y = element_text(angle=0))
  ggsaveNF(out("ClusterEnrichments_manual_cleaned_",suppx,".pdf"), w=w,h=h)
}


# . LSC marker plot ---------------------------------------------------------
#gg <- c("Bcat1", "Hif1a", "Myc", "Gata1", "Prss34", "Hba-a1", "Irf8", "S100a9", "Ltf")
#gg <- c("Itgam","Ly6c1","Ltf","S100a9","Cd55","Fcer1a","Prss34","Gata1","Hba-a1","Itga2b","Pf4")
gg <- list(
  main=c("Pf4","Itga2b","Hba-a1","Gata1","Cd55","Prss34","S100a9","Itgam","Ltf"),
  supp=c("Bcat1","Meis1","Kit","Sox4","Pf4","Itga2b","Hba-a1","Gata1","Cd55","Prss34","Mpo","S100a9","Itgam","Ltf","Csf1r","F13a1")
)

cds <- mobjs[["leukemia"]]
stopifnot(all(unlist(gg) %in% row.names(cds)))
cds$Cluster.number <- monocle3::clusters(cds)

suppx <- "main"
for(suppx in c("supp", "main")){
  pDT <- DotPlotData(cds = cds, markers = gg[[suppx]], cols = "Cluster.number")
  pDT[, scale := scale(mean), by=c("Gene")]
  pDT$Gene <- factor(pDT$Gene, levels=gg[[suppx]])
  pDT[,Cluster.number := as.numeric(Cluster.number)]
  
  cl.labels <- clusters.plot.cts[[suppx]]
  pDTx <- merge(pDT, melt(cl.labels), by.x="Cluster.number", by.y="value")
  pDTx[, celltype := factor(L1, levels=rev(names(cl.labels)))]
  
  pDTx <- pDTx[Cluster.number %in% unlist(cl.labels)]
  pDTx[, Cluster.number.fac := factor(as.character(Cluster.number), levels=rev(as.character(unlist(cl.labels))))]
  
  h=length(unique(pDTx$Cluster.number)) * 0.08 + 0.2
  w=length(unique(pDTx$Gene)) * 0.08 + 0.2
  ggplot(pDTx, aes(y=Cluster.number.fac, x=Gene, color=scale, size=percentage)) + 
    geom_point() +
    #geom_point(color="black", shape=1) + 
    scale_size_continuous(range=c(0,5), limits = c(0,100)) +
    scale_color_gradientn(colours = c("lightgrey", "grey", "#1f78b4", "black")) +
    facet_grid(celltype ~ ., space="free", scales = "free") +
    themeNF(rotate = TRUE) +
    ylab("Cluster") + xlab("Gene")
  ggsaveNF(out("Markers_Leukemia_",suppx,".pdf"), w=w,h=h)
}



# Larry Barplot -----------------------------------------------------------
#gg <- c("Smarcd2","Smarcd1","Brd9","Kmt2d","Ncoa6","Kmt2a","Men1","Wdr82","Setd1a","Setd1b","Chd4", "Hdac3","Setdb1","Stag2")
gg <- dla.healthy$all
ctx <- c("HSC","GMP","Granulocyte","Mono", "B","MEP","EryA")
refx <- fread(dirout_load("SCRNA_06_01_Markers")("GEO.txt.gz"), skip = 1)
x <- melt(refx[NAME %in% gg], id.vars = c("UNIQUD", "NAME"))
#refx[grepl("Ncoa", NAME)]
x <- x[variable %in% ctx]
x[, celltype := factor(variable, levels=ctx)]     
x[, gene := factor(NAME, levels=gg)]     
x <- x[, .(value = mean(value)),by=c("gene", "celltype")]
ggplot(x, aes(x=celltype,y=value)) + 
  geom_col() +
  facet_wrap(~gene, ncol = 10, scales = "free") +
  themeNF(rotate=TRUE) +
  labs(x="", y="")
ggsaveNF(outBase("Larry_Genes.pdf"), w=4.5,h=6, guides=TRUE)
