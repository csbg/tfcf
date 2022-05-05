source("src/00_init.R")
base.dir <- "FIG_02_scRNA_UMAPs/"
outBase <- dirout(base.dir)


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
marker.signatures.DA <- lapply(TISSUES, function(tx) fread(inDir.funcs[[tx]]("SigDA.tsv")))
names(marker.signatures.DA) <- TISSUES

# Marker overlaps
marker.overlaps <- fread(dirout_load("SCRNA_06_03_MarkerOverlaps")("Enrichments.tsv"))

# Doublet scores?
ff <- list.files(dirout_load("SCRNA_03_01_Duplets")(""), pattern="^Duplet_Scores_.+.tsv$", full.names = TRUE)
names(ff) <- gsub("^Duplet_Scores_(.+).tsv$", "\\1", basename(ff))
ff <- lapply(ff, fread)
duplet.scores <- rbindlist(ff, idcol="sample")
duplet.scores[, rn := paste0(rn, "_", sample)]

# Cell annotations
annList <- lapply(names(inDir.funcs), function(inDir.current){
  ann <- fread(inDir.funcs[[inDir.current]]("Annotation.tsv"))
  ann[, perturbed := !(mixscape_class.global %in% c("NP", "NTC") | is.na(mixscape_class.global))]
  ann[, gene := gsub("_.+", "", guide)]
  ann
  })
annList <- rbindlist(annList, fill=TRUE)
annList$dupletScore <- duplet.scores[match(annList$rn, rn),]$dupletScore


# Other projections
umap.proj <- list(
  original=readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjMonocle.RDS")),
  izzo = readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjIzzo.RDS")),
  in.vivo = readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjVivo.RDS")),
  in.vivo.X = readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjVivoX.RDS"))
)

# # Update ex vivo and leukemia to use in.vivo.X UMAPs
# res <- rbind(
#   annList[tissue == "in.vivo"], 
#   merge(
#     annList[tissue != "in.vivo"][,-c("UMAP1", "UMAP2")], 
#     setNames(umap.proj[["in.vivo.X"]][tissue != "in.vivo"][,c("rn", "UMAP_1", "UMAP_2")], c("rn", "UMAP1", "UMAP2")),
#     by="rn"
#     )
# )
# stopifnot(length(union(res$rn, annList$rn)) == length(intersect(res$rn, annList$rn)))
# annList <- res

# load Monocle Objects
mobjs <- list()
tissuex <- PATHS$SCRNA$MONOCLE.NAMES[1]
for(tissuex in PATHS$SCRNA$MONOCLE.NAMES){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}
for(tissuex in names(mobjs)){
  mobjs[[tissuex]]$CellType <- annList[match(colnames(mobjs[[tissuex]]), rn)]$Clusters
  fData(mobjs[[tissuex]])$gene_short_name <- row.names(fData(mobjs[[tissuex]]))
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


# . Dupletscores ------------------------------------------------------------
ggplot(annList, aes(x=UMAP1, y=UMAP2)) + 
  themeNF() +
  stat_summary_hex(aes(z=dupletScore),fun=mean, bins=100) +
  facet_grid(. ~ tissue)
ggsaveNF(outBase("Duplets.pdf"),w=3,h=10, guides = TRUE)


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
    pDT <- pDT[,.(N=.N, dupletScore=mean(dupletScore)), by=c("hex.x", "hex.y", "Clusters")]
    pDT[, sum := sum(N), by=c("hex.x", "hex.y")]
    pDT[, frac := N / sum]
    pDT <- pDT[order(frac, decreasing = TRUE)][,head(.SD, 1), by=c("hex.x", "hex.y")]
    #pDT <- pDT[frac > 0.25]
    pDT.labels <- pDT[, .(hex.x = median(hex.x), hex.y=median(hex.y)), by=c("Clusters")]
    pDT[, Clusters := cleanCelltypes(Clusters,twoLines = FALSE)]
    pDT.labels[, Clusters := cleanCelltypes(Clusters,twoLines = TRUE)]
    p <- ggplot(pDT, aes(x=hex.x, y=hex.y)) +
      themeNF() + xu + yu + scale_color_manual(values=COLORS.CELLTYPES.scRNA.ainhoa)
    if(!is.null(xDT.ref)) p <- p + geom_hex(data=xDT.ref, fill="lightgrey", bins=100, aes(x=UMAP_1, y=UMAP_2))
    p <- p +
      geom_point(aes(color=Clusters), size=0.5) + 
      geom_text(data=pDT.labels, aes(label=Clusters), lineheight = 0.8)
    ggsaveNF(out("UMAP_Celltypes_",x,".pdf"), w=1.5,h=1.5, plot=p)
  }
}


# . Signatures --------------------------------------------------------------
tx <- "in.vivo"
tx <- "leukemia"
tx <- "ex.vivo"
for(tx in names(inDir.funcs)){
  # out directory
  out <- dirout(paste0(base.dir, "/", tx))
  
  # annotation
  ann <- annList[tissue ==  tx & timepoint != "28d"]
  
  # Plot signatures scores in each cluster
  pDT <- merge(
    ann[, c("rn", "Clusters"),with=F],
    marker.signatures[, c("rn", "DB", "value", "FinalName"),with=F],
    by="rn")
  pDT <- pDT[, mean(value), by=c("Clusters", "FinalName")]
  pDT[, V1 := scale(V1), by="FinalName"]
  pDT[, Clusters := cleanCelltypes(Clusters)]
  #pDT <- hierarch.ordering(pDT, "Clusters", "FinalName", "V1")
  pDT <- hierarch.ordering(pDT, "FinalName", "Clusters", "V1")
  ggplot(pDT, aes(x=FinalName, y=Clusters, fill=V1)) + 
    geom_tile() +
    themeNF(rotate = TRUE) +
    scale_fill_gradient2(low="blue", high="red") +
    xlab("External gene signatures") +
    ylab("Cell type\n(this study)")
  ggsaveNF(out("Supp_Signatures.pdf"), w=1.2,h=length(unique(pDT$Clusters)) * 0.05 + 0.5)
  
  # Plot expression of marker genes in each cluster
  
  # Plot fraction of predicted cells in each cluster
}


# . Signatures DE with CRISPR KOs -------------------------------------------
tx <-TISSUES[1]
for(tx in TISSUES){
  
  out <- dirout(paste0(base.dir, "/", tx))
  
  pDT <- marker.signatures.DA[[tx]]
  pDT <- merge(pDT, SIGS.USE.DA, by=c("sig"))
  
  # Summarize across guides
  pDT[gene == "Pu.1", gene := "Spi1"]
  pDT <- pDT[, .(
    d=mean(d, na.rm=TRUE), 
    dir=length(unique(sign(d[!is.na(d)])))==1, 
    padj=sum(padj < 0.01), 
    N=.N), by=c("sample", "FinalName", "gene")]
  pDT[dir == FALSE, padj := 0]
  pDT[dir == FALSE, d := NA]
  
  # summarize across samples
  pDT <- pDT[, .(
    d=mean(d, na.rm=TRUE), 
    dir=length(unique(sign(d[!is.na(d)])))==1, 
    padj=sum(padj), 
    N=sum(N)), by=c("FinalName", "gene")]
  pDT[dir == FALSE, padj := 0]
  pDT[dir == FALSE, d := NA]
  
  # setup for plotting and plotting
  pDT[padj == 0 | is.na(d), d := 0]
  pDT[, sig.perc := padj / N]
  pDT[,d_cap := pmin(abs(d), 5) * sign(d)]
  pDT$FinalName <- factor(pDT$FinalName, levels=rev(SIGS.USE.DA$FinalName))
  #pDT <- hierarch.ordering(dt = pDT, toOrder = "FinalName", orderBy = "gene", value.var = "d")
  pDT <- hierarch.ordering(dt = pDT, toOrder = "gene", orderBy = "FinalName", value.var = "d")
  ggplot(pDT, aes(x=gene,y=FinalName, color=d, size=sig.perc)) + 
    themeNF(rotate = TRUE) +
    scale_size_continuous(name="% sign.", range = c(0,4)) +
    scale_color_gradient2(name="delta",low="#1f78b4", high="#e31a1c") +
    geom_point() +
    geom_point(shape=1, color="lightgrey") +
    ylab("Marker signature") + xlab("")
  ggsaveNF(out("MarkerSignatures_DA.pdf"),w=2,h=0.7)
  
  dla <- fread("metadata/FIGS_Order_Fig2_CFs.tsv")$Factor
  pDT <- pDT[gene %in% dla]
  pDT$gene <- factor(as.character(pDT$gene), levels=dla)
  ggplot(pDT, aes(x=gene,y=FinalName, color=d, size=sig.perc)) + 
    themeNF(rotate = TRUE) +
    scale_size_continuous(name="% sign.", range = c(0,4)) +
    scale_color_gradient2(name="delta",low="#1f78b4", high="#e31a1c") +
    geom_point() +
    geom_point(shape=1, color="lightgrey") +
    ylab("Marker signature") + xlab("")
  ggsaveNF(out("MarkerSignatures_DA_DLA.pdf"),w=1.5,h=0.7)
  
  # export table
  write.tsv(pDT, out("MarkerSignatures_DA.tsv"))
}


# . Marker gene expression --------------------------------------------------
# for(tx in names(mobjs)){
#   out <- dirout(paste0(base.dir, "/", tx))
#   plot_genes_by_group(mobjs[[tx]], markers = marker.genes, group_cells_by = "CellType") + 
#     scale_size_continuous(range=c(0,5)) +
#     themeNF(rotate = TRUE) +
#     xlab("Cell type")
#   ggsaveNF(out("Markers_Clusters.pdf"), w=1.5,h=2.5)
# }

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
  typex <- "basic_in.vivo_noMixscape"
  for(typex in gsub("Guides_Fisher_Mixscape_(.+).pdf", "\\1", list.files(inDir(""), pattern="Guides_Fisher_Mixscape_.*.pdf"))){
    fish.file <- inDir("Guides_Fisher_Mixscape_",typex,".tsv")
    if(!file.exists(fish.file)) next
    
    fish.full <- fread(fish.file)
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
        padj=sum(padj < 0.01)/.N,
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
pDT <- annList[tissue != "leukemia"][timepoint != "d28"][,.N, by=c("Clusters", "tissue")]
pDT[, sum := sum(N), by="tissue"]
pDT[, perc := N/sum*100]
pDT[, Clusters := cleanCelltypes(Clusters, reverse = TRUE)]
ggplot(pDT, aes(x=Clusters, y=perc, fill=tissue)) + 
    geom_bar(stat="identity", position=position_dodge2(preserve = "single")) + 
    #scale_y_log10() +
    themeNF(rotate=TRUE) +
    ylab("Number of cells") + xlab("")
ggsaveNF(outBase("CellCounts.pdf"), w=1.5,h=1, guides = TRUE)

pDT <- annList[tissue != "leukemia"][timepoint != "d28"][,.N, by=c("gene", "Clusters", "tissue")][!is.na(gene)]
# ggplot(pDT,aes(x=N, group=CRISPR_Cellranger)) + 
#   stat_ecdf() +
#   facet_grid(tissue ~ .)
pDT[, Clusters := cleanCelltypes(Clusters, reverse = TRUE)]
ggplot(pDT,aes(x=Clusters, y=N, color=tissue)) + 
  themeNF(rotate=T) +
  geom_boxplot(coef=10^10, position=position_dodge2(preserve = "single")) + 
  ylab("Number of cells") + xlab("")
ggsaveNF(outBase("CellCounts_InVivoSparse.pdf"), w=1.5,h=1, guides = TRUE)


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
ggplot(pDT, aes(x=Guide1, y=Guide2, color=sample_broad)) + 
  geom_point() +
  themeNF() +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(. ~ tissue)
ggsaveNF(outBase("NTC_clonalityScatter.pdf"), w=2, h=1)

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
  for(x in unique(pDT.top[perturbed == TRUE]$gene)){
    # FIX MIXSCAPE
    #pDTg <- rbind(pDT.ntc, pDT.top[grepl(paste0("^", x), guide) & mixscape_class.global == "KO"])
    pDTg <- rbind(pDT.ntc, pDT.top[grepl(paste0("^", x), guide)])
    pDTg$plot <- x
    pDT.final <- rbind(pDT.final, pDTg)
  }
  pDT.final$plot <- relevel(factor(pDT.final$plot),ref = "NTC")
  ggplot(pDT.final, aes(x=UMAP1, y=UMAP2)) + 
    themeNF() +
    geom_hex(data=pDT.final[mixscape_class.global == "NTC"], bins=100, fill="lightgrey") +
    geom_hex(data=pDT.final[mixscape_class.global != "NTC" | plot == "NTC"], bins=100) +
    scale_fill_gradientn(colours=c("#1f78b4", "#e31a1c")) +
    facet_wrap(~plot, ncol=6) + 
    xu + yu
  ggsaveNF(out("UMAP_Guides_all.pdf"), w=4,h=4)
  
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
fish.bcells <- fread(dirout_load(base.dir)("cluster.enrichments/Cluster_enrichments","_basic", "_in.vivo", "_noMixscape", "_14d",".tsv"))
fish.enrich.broad <- fread(dirout_load(base.dir)("cluster.enrichments/Cluster_enrichments_broad_in.vivo_noMixscape_14d",".tsv"))
fish.EryVsMye <- fread(dirout_load(base.dir)("cluster.enrichments/Cluster_enrichments","_eryVsMye", "_in.vivo", "_noMixscape", "_14d",".tsv"))


# . Plot distribution of one gene -------------------------------------------
pDT <- ann[markers == "lin-"][timepoint == "14d"]
gg <- "Rbbp4"
pDT <- pDT[(gene == gg & perturbed == TRUE) | mixscape_class == "NTC"]
pDT <- pDT[, .N, by=c("gene", "Clusters")]
pDT[, sum := sum(N), by=c("Clusters")]
pDT[, rel2NTCs := N/(sum-N)]
pDT[rel2NTCs == Inf, rel2NTCs := 1]
pDT <- pDT[gene != "NTC"]
ggplot(pDT, aes(x=cleanCelltypes(Clusters, reverse = FALSE), y=rel2NTCs)) + 
  themeNF(rotate=TRUE) +
  geom_bar(stat="identity") +
  ggtitle(gg)


# . NTC Clusters depletion ------------------------------------------------------------
pDT <- copy(ann)
clx <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjMonocle_Clusters.RDS"))
pDT$cluster.numeric <- clx[match(pDT$rn, rn)]$functional.cluster
pDT.full <- copy(pDT)
pDT <- pDT[!is.na(gene)][, .(sum=.N, sumNTC = sum(!is.na(gene) & gene=="NTC")), by=c("cluster.numeric", "Clusters", "timepoint")]
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
  geom_text_repel(data=pDT[sumKO / 10 > sumNTC], aes(label=cluster.numeric), color="black") +
  xlab("Number of cells with KO guide + 1") + ylab("Number of cells with NTC guide + 1")
ggsaveNF(out("NTC_depleted_clusters.pdf"),w=3,h=1.5)

ggplot(pDT.full, aes(x=UMAP1, y=UMAP2)) + 
  themeNF() +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="#ff7f00") +
  geom_text(data=pDT.full[,.(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by=c("cluster.numeric")], aes(label=cluster.numeric))
ggsaveNF(out("NTC_depleted_clusters_UMAP.pdf"),w=1.5,h=1.5)



# . UMAP of timepoint / markers ---------------------------------------------------------
ggplot(ann[perturbed == FALSE], aes(x=UMAP1, y=UMAP2)) + 
  themeNF() +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~markers + timepoint, ncol=2) +
  xu + yu
ggsaveNF(out("UMAP_Groups.pdf"), w=1,h=1)

# . UMAP of samples -------------------------------------------------------
ggplot(ann[perturbed == FALSE], aes(x=UMAP1, y=UMAP2)) + 
  themeNF() +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~sample_broad, ncol=3) +
  xu + yu
ggsaveNF(out("UMAP_Samples.pdf"), w=2,h=2)

# . Plot top guides ---------------------------------------------------------
# pDT.top <- ann[timepoint == "14d"]
# pDT.ntc <- pDT.top[mixscape_class.global == "NTC"]
# pDT.final <- copy(pDT.ntc)
# pDT.final$plot <- "NTC"
# for(x in c("Wdr82", "Rcor1", "Ehmt1", "Setdb1", "Hdac3")){
#   pDTg <- rbind(pDT.ntc, pDT.top[grepl(paste0("^", x), guide) & perturbed == TRUE])
#   pDTg$plot <- x
#   pDT.final <- rbind(pDT.final, pDTg)
# }
# pDT.final$plot <- relevel(factor(pDT.final$plot), ref = "NTC")
# ggplot(pDT.final, aes(x=UMAP1, y=UMAP2)) + 
#   themeNF() +
#   geom_hex(data=pDT.final[mixscape_class.global == "NTC"], bins=100, fill="lightgrey") +
#   geom_hex(data=pDT.final[mixscape_class.global != "NTC" | plot == "NTC"], bins=100) +
#   scale_fill_gradientn(colours=c("#1f78b4", "#e31a1c")) +
#   facet_wrap(~plot, ncol=3) +
#   xu + yu
# ggsaveNF(out("UMAP_Guides.pdf"), w=2,h=2)


# . Plot manual enrichments -------------------------------------------------
dla <- fread("metadata/FIGS_02_CFs.main.txt", header = F)$V1

# B cell enrichments
pDT <- fish.bcells[Clusters == "Imm. B-cell"][abs(log2OR) > 0.5]
pDT$gene <- factor(pDT$gene, levels = pDT[order(log2OR)]$gene)
ggplot(pDT, aes(x=log2OR, y=gene, size=sig.perc, color = log2OR)) + 
  themeNF(12) +
  geom_vline(xintercept = 0) +
  scale_size_continuous(range = c(1,4)) +
  geom_point() + 
  geom_point(shape=1, color="lightgrey") + 
  scale_color_gradient2(low="blue", high="red") +
  ylab("") + xlab("Enrichment in immature B-cells")
ggsaveNF(out("ClusterEnrichments_manual_BcellsOnly.pdf"), w=1,h=1.5)

# Plot enrichments with DLA order
pDT <- fish.enrich.broad[gene %in% dla]
pDT$gene <- factor(pDT$gene, levels = dla)
pDT[, Clusters := cleanCelltypesBroad(Clusters, clean=FALSE)]
ggplot(pDT, aes(x=gene, y=Clusters, size=sig.perc, color=log2OR_cap)) + 
  themeNF(rotate=TRUE) +
  scale_color_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
  scale_size_continuous(name="% sign.", range = c(0,5)) +
  geom_point() +
  geom_point(shape=1, color="lightgrey") +
  xlab("Gene") + ylab("Cell type")
ggsaveNF(out("ClusterEnrichments_manual_cleaned.pdf"), w=1.5,h=0.7)

# Plot my vs GMP enrichments with DLA order
pDT <- fish.EryVsMye[Clusters == "GMP"][gene %in% dla]#[log2OR > 0 | log2OR < -2]
pDT$gene <- factor(pDT$gene, levels = pDT[order(log2OR)]$gene)
ggplot(pDT, aes(x=log2OR, y=gene, size=sig.perc, color=log2OR_cap)) + 
  themeNF(rotate=TRUE) +
  scale_color_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
  scale_size_continuous(name="% sign.", range = c(0,5), limits = c(0,1)) +
  geom_point() +
  ylab("") + xlab("Myeloid enrichment vs Erythroid (log2OR)") +
  geom_vline(xintercept = 0)
ggsaveNF(out("ClusterEnrichments_manual_MEPvsGMP.pdf"), w=1,h=1)

# . Plot displasia guides ---------------------------------------------------------
pDT.top <- ann[timepoint == "28d"][grepl("OP1", sample)][mixscape_class.global == "NTC" | perturbed == TRUE]
pDT.ntc <- pDT.top[mixscape_class.global == "NTC"]
pDT.final <- copy(pDT.ntc)
pDT.final$plot <- "NTC"
for(x in c("Kmt2a","Kmt2d","Smarcd2","Smarcd1","Brd9")){
  pDTg <- rbind(pDT.ntc, pDT.top[grepl(paste0("^", x), guide)])
  pDTg$plot <- x
  pDT.final <- rbind(pDT.final, pDTg)
}
pDT.final$plot <- relevel(factor(pDT.final$plot), ref = "NTC")
ggplot(pDT.final, aes(x=UMAP1, y=UMAP2)) + 
  themeNF() +
  geom_hex(data=pDT.final[mixscape_class.global == "NTC"], bins=100, fill="lightgrey") +
  geom_hex(data=pDT.final[mixscape_class.global != "NTC" | plot == "NTC"], bins=100) +
  scale_fill_gradientn(colours=c("#1f78b4", "#e31a1c")) +
  facet_wrap(~plot, ncol=3) +
  xu + yu
ggsaveNF(out("UMAP_Guides_displasia.pdf"), w=2,h=1.5)

# . Plot displasia numbers -----------------------------------------------
pDT <- copy(ann[markers=="lin-"][grepl("OP1", sample)][mixscape_class.global == "NTC" | perturbed == TRUE][gene %in% c("Kmt2a","Kmt2d","Smarcd2","Smarcd1","Brd9", "NTC")])
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
dla <- fread("metadata/FIGS_Order_Fig2_CFs.tsv")
pDT <- fish.enrich[day == "day7"][gene %in% dla$Factor][Clusters %in% c("GMP", "MkP", "Eo/Ba", "LSK")]
pDT$gene <- factor(pDT$gene, levels = dla$Factor)
#pDT$Complex <- dla[match(pDT$gene, Factor)]$Complex
pDT[, Clusters := cleanCelltypes(Clusters, clean=FALSE, twoLines = FALSE, order = TRUE, reverse = TRUE)]
ggplot(pDT, aes(x=gene, y=Clusters, size=sig.perc, color=log2OR_cap)) + 
  themeNF(rotate=TRUE) +
  scale_color_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
  scale_size_continuous(name="% sign.", range = c(0,4)) +
  geom_point() +
  geom_point(shape=1, color="lightgrey") +
  xlab("Gene") + ylab("Cell type") +
  #facet_grid(Clusters ~ . , space="free", scales = "free") +
  theme(strip.text.y = element_text(angle=0))
ggsaveNF(out("ClusterEnrichments_manual_cleaned_day7.pdf"), w=1.6,h=0.6)

# . Plot enrichments with late GMPs from day 9
dla <- fread("metadata/FIGS_Order_Fig2_CFs.tsv")
pDT <- fish.enrich[day == "day9"][gene %in% dla$Factor][Clusters %in% c("Late GMP")][log2OR != 0]
pDT$gene <- factor(pDT$gene, levels = pDT[order(log2OR)]$gene)
#pDT$Complex <- dla[match(pDT$gene, Factor)]$Complex
pDT[, Clusters := cleanCelltypes(Clusters, clean=FALSE, twoLines = FALSE, order = TRUE, reverse = TRUE)]
ggplot(pDT, aes(y=gene, x=log2OR_cap, size=sig.perc, color=log2OR_cap)) + 
  themeNF(rotate=TRUE) +
  scale_color_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
  scale_size_continuous(name="% sign.", range = c(0,4)) +
  geom_point() +
  geom_point(shape=1, color="lightgrey") +
  xlab("log2 OR") + ylab("Gene") +
  #facet_grid(Clusters ~ . , space="free", scales = "free") +
  theme(strip.text.y = element_text(angle=0))
ggsaveNF(out("ClusterEnrichments_manual_lateGMPs_day9.pdf"), w=1,h=1)


# LEUKEMIA---------------------------------------------------------

# . load data ---------------------------------------------------------
inDir.current <- "leukemia"
out <- dirout(paste0(base.dir, "/", inDir.current))
ann <- annList[tissue == inDir.current]
#ann <- merge(ann[,-c("UMAP1", "UMAP2"),with=F], setNames(umap.proj[["in.vivo.X"]][,c("rn", "UMAP_1", "UMAP_2")], c("rn", "UMAP1", "UMAP2")), by="rn")
abs <- fread(inDir.funcs[[inDir.current]]("Antibodies.tsv"))
abs <- merge(abs[,-c("UMAP1", "UMAP2"),with=F], setNames(umap.proj[["in.vivo.X"]][,c("rn", "UMAP_1", "UMAP_2")], c("rn", "UMAP1", "UMAP2")), by="rn")
# FIX MIXSCAPE
fish.enrich <- fread(dirout_load(base.dir)("cluster.enrichments/Cluster_enrichments_basic_leukemia_noMixscape_6d.tsv"))

# . Antibodies on UMAP ----------------------------------------------------
ggplot(abs, aes(x=UMAP1, y=UMAP2)) +
  stat_summary_hex(bins = 100, aes(z=pmin(abs(Signal.norm), 2) * sign(Signal.norm)),fun=mean) +
  scale_fill_gradient2(low="blue", high="red") +
  facet_wrap(~Antibody) +
  themeNF() +
  xu + yu
ggsaveNF(out("Antibodies_UMAP.pdf"), w=2, h=2)

# . Combine signatures with Antibodies --------------------------------------
pDT <- abs[, .(Signal=mean(Signal.norm, na.rm=TRUE)), by=c("Clusters", "Antibody")]
pDT[, Clusters := cleanCelltypes(Clusters)]
#pDT <- hierarch.ordering(pDT, "Clusters", "Antibody", "Signal")
pDT <- hierarch.ordering(pDT, "Antibody", "Clusters", "Signal")
ggplot(pDT, aes(y=Clusters,x=Antibody, fill=Signal)) +
  themeNF() +
  geom_tile() +
  scale_fill_gradient2(low="blue", high="red") +
  xRot() +
  ylab("Cell types")
ggsaveNF(out("Antibodies_Average.pdf"), w=1, h=0.7)

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
pDT <- ann[perturbed == TRUE | mixscape_class == "NTC"][,.N, by=c("Phase", "CRISPR_Cellranger", "timepoint", "markers")]
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
dla <- fread("metadata/FIGS_Order_Fig2_CFs.tsv")
pDT <- fish.enrich[gene %in% dla$Factor]
pDT$gene <- factor(pDT$gene, levels = dla$Factor)
#pDT$Complex <- dla[match(pDT$gene, Factor)]$Complex
pDT[, Clusters := cleanCelltypes(Clusters, clean=FALSE, twoLines = FALSE, order = TRUE, reverse = TRUE)]
ggplot(pDT, aes(x=gene, y=Clusters, size=sig.perc, color=log2OR_cap)) + 
  themeNF(rotate=TRUE) +
  scale_color_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
  scale_size_continuous(name="% sign.", range = c(0,4)) +
  geom_point() +
  geom_point(shape=1, color="lightgrey") +
  xlab("Gene") + ylab("Cell type") +
  #facet_grid(Clusters ~ . , space="free", scales = "free") +
  theme(strip.text.y = element_text(angle=0))
ggsaveNF(out("ClusterEnrichments_manual_cleaned.pdf"), w=1.6,h=0.6)


# Cancer vs Healthy -------------------------------------------------------
out <- dirout(paste0(base.dir, "/", "CvH"))

pDT <- list(
  leukemia = fread(outBase("leukemia/Cluster_enrichments_basic_all.tsv")),
  normal = fread(outBase("ex.vivo/Cluster_enrichments_basic_all.tsv"))
)
pDT <- rbindlist(pDT, idcol = "tissue")
pDT[, Clusters := cleanCelltypes(Clusters, reverse = FALSE, clean=FALSE)]

ggplot(pDT, aes(x=Clusters, y=tissue, size=sig.perc, color=log2OR_cap)) + 
  themeNF(rotate = TRUE) +
  scale_color_gradient2(name = "log2(OR)", low="blue", midpoint = 0, high="red") +
  scale_size_continuous(name = "% sign", range=c(0,5)) +
  geom_point() +
  geom_point(shape=1, color="lightgrey") +
  facet_grid(gene ~ .) +
  theme(strip.text.y = element_text(angle=0)) +
  xlab("Cell types") + ylab("")
ggsaveNF(out("Cluster_enrichments.pdf"), w=1.5,h=4.5, guides = TRUE)
