source("src/00_init.R")

# Load packages and functions
require(umap)
require(igraph)
require(nebula)
require(fgsea)
library(SingleR)
source("src/FUNC_Monocle_PLUS.R")

# Figure out command line arguments (which tissue to analyze)
args = commandArgs(trailingOnly=TRUE)
# args[1] <- "leukemia"
# args[1] <- "in vivo"
# args[1] <- "in vitro"

# Define output directory based on the tissue
baseDir <- "FULLINT_10_01_BasicAnalysis"
baseDir.add <- if(length(args) == 0){
  "combined"
} else {
  if(!args[1] %in% c("leukemia", "in vitro", "in vivo")){
    stop("Argument ", args[1], " not valid. Expecting 'leukemia', 'in vivo', or 'in vitro'")
  } else {
    args[1]
  }
}
out <- dirout(paste0(baseDir, "_", make.names(baseDir.add), "/"))


# LOAD DATA ---------------------------------------------------------------
# load single cell data
# Subset data and reprocess based on the tissue (command line argument)
if(baseDir.add == "combined"){
  print("Loading full dataset")
  load(PATHS$FULLINT$Monocle)
}else{
  monocle.file <- out("MonocleObject.RData")
  if(file.exists(monocle.file)){
    print("Loading processed partial dataset")
    load(monocle.file)
  } else {
    print("Processing partial dataset for tissue:")
    print(args[1])
    load(PATHS$FULLINT$Monocle)
    monocle.obj <- monocle.obj[,monocle.obj$tissue == args[1]]
    monocle.obj <-
      preprocess_cds(monocle.obj, verbose = TRUE) %>%
      reduce_dimension(preprocess_method = "PCA", verbose = TRUE)
    set.seed(42)
    monocle.obj <-
      align_cds(monocle.obj, alignment_group = "sample", verbose = TRUE) %>%
      reduce_dimension(
        reduction_method = "UMAP",
        preprocess_method = "Aligned",
        verbose = TRUE)
    set.seed(12121)
    monocle.obj = cluster_cells(monocle.obj)
    
    additional.info <- additional.info[unique(monocle.obj$sample)]
    AGG.CSV <- AGG.CSV[sample_id %in% monocle.obj$sample]
    
    save(monocle.obj, additional.info, AGG.CSV, file=monocle.file)
  }
}
fData(monocle.obj)$gene_short_name <- row.names(fData(monocle.obj))


# Markers
marker.genes <- fread("metadata/markers.csv")

# Marker signautres
ff <- list.files(dirout_load("FULLINT_08_01_Markers")(""), pattern="Signatures_")
ff <- dirout_load("FULLINT_08_01_Markers")(ff)
names(ff) <- gsub("^Signatures_(.+?).csv$", "\\1", basename(ff))
ff <- ff[!grepl("Augmented_2021", ff)]
marker.signatures <- lapply(ff, function(fx) as.matrix(read.csv(fx)))

# SingleR
ff <- list.files(dirout_load("FULLINT_05_01_SingleR")(""), pattern = "cell_types_.*.csv", full.names = TRUE)
singleR.res <- setNames(lapply(ff, fread), gsub("cell_types_(.+).csv", "\\1", basename(ff)))

# CytoTRACE
tryCatch({
  (load(dirout_load("FULLINT_06_01_CytoTRACE")("CytoTRACE.RData")))
}, error=function(e){
  message("CytoTRACE import failed")
})

# ChIP Targets
(load(dirout_load("CHIP_20_01_Peaks_julen")("ChIP.Targets.RData")))

# ANNOTATION ------------------------------------------------------
sann <- fread("metadata/annotation.tsv", sep="\t")

# Collect ANNOTATION --------------------------------------------------------------
ann <- data.table(data.frame(colData(monocle.obj)@listData), keep.rownames = TRUE)
ann$MonocleClusters <- as.character(monocle.obj@clusters$UMAP$clusters[ann$rn])
umap <- setNames(data.table(reducedDims(monocle.obj)$UMAP, keep.rownames = TRUE), c("rn", "UMAP1", "UMAP2"))
ann <- merge(ann, umap, by="rn", all=TRUE)
if("cytoRes" %in% ls()) ann$CytoTRACE <- cytoRes$CytoTRACE[ann$rn]
#ann$tissue <- sann[match(gsub("_.+", "", ann$sample), sample),]$tissue


# Final clusters ----------------------------------------------------------
ann[, Clusters := MonocleClusters]

# define files ------------------------------------------------------------
neb.file <- out("DEG_Results_nebula.RData")


# SETUP ENDS HERE ---------------------------------------------------------





# Output annotation -------------------------------------------------------
write.tsv(ann, out("Annotation.tsv"))


# ADDITIONAL QC  --------------------------------------------------------
qcm <- "nCount_RNA"
for(qcm in c("percent.mt", "nFeature_RNA", "nCount_RNA")){
  print(qcm)
  ann$measure <- ann[[qcm]]
  p <- ggplot(ann, aes(y=measure + 0.1, x=factor(Clusters))) + 
    geom_violin(color=NA, fill="lightblue") + 
    geom_boxplot(fill=NA, coef=Inf) +
    scale_y_log10() +
    theme_bw(12) +
    ylab(qcm) +
    ggtitle(qcm)
  ggsave(out("QC_", qcm, "_Clusters.pdf"), w=5,h=4, plot=p)
  
  ggplot(ann, aes(x=UMAP1, y=UMAP2)) + 
    theme_bw(12) +
    stat_summary_hex(bins=100, aes(z=measure),fun=mean) +
    scale_fill_hexbin() +
    #scale_fill_gradient(low="lightgrey", high="blue") +
    ggtitle(qcm)
  ggsave(out("QC_", qcm, "_UMAP.pdf"), w=5,h=4)
  }
ann$measure <- NULL
#ann[,cluster.qual.keep :=TRUE]
#ann[Clusters %in% ann[,median(percent.mt), by="Clusters"][V1 < 1]$Clusters, cluster.qual.keep := FALSE]


# # CELLRANGER -------------------------------------------
# if("AGG.CSV" %in% ls()){
#   ann.exp <- merge(ann, AGG.CSV[,c("sample_id", "i"),with=F], by.x="sample", by.y="sample_id", all.x=TRUE)
#   stopifnot(!any(is.na(ann.exp$i)))
#   ann.exp[,Barcode := paste0(gsub("-.+$", "", rn), "-", i)]
#   write.table(ann.exp[,c("Barcode", "UMAP1", "UMAP2"),with=F], file=out("Cellranger_UMAP.csv"), sep=",", col.names = c("Barcode", "UMAP-1", "UMAP-2"), quote=F, row.names = F)
#   write.table(ann.exp[,c("Barcode", "mixscape_class"),with=F], file=out("Cellranger_MIXSCAPE.csv"), sep=",", col.names = c("Barcode", "MIXSCAPE"), quote=F, row.names = F)
#   write.table(ann.exp[,c("Barcode", "Clusters"),with=F], file=out("Cellranger_Clusters.csv"), sep=",", col.names = c("Barcode", "Clusters_Seurat"), quote=F, row.names = F)
# }



# TRANSFER / PREDICT CLUSTERS IN FULL DATASET ----------------------------------------
transf.file <- out("TransferClusters.tsv")
if(baseDir.add == "in vivo"){
  if(file.exists(transf.file)){
    
  } else {
    print("Transferring cell clusters to full datatset")
    fullDS <- function(){load(PATHS$FULLINT$Monocle); return(monocle.obj)}
    monocle.full <- fullDS()
    
    ref=SummarizedExperiment(
      assays = SimpleList(logcounts=SCRNA.TPXToLog(SCRNA.RawToTPX(counts(monocle.obj), 1e6))), 
      colData=DataFrame(label.main=ann$Clusters, row.names = ann$rn))
    
    res <- SingleR(
      test = counts(monocle.full),
      ref = ref,
      labels = ann$Clusters)
    
    write.tsv(data.table(labels=res$labels, res@rownames), transf.file)
  }
}



# SAMPLES -----------------------------------------------------------------
ggplot(ann, aes(x=sample)) +
  theme_bw() +
  xRot() +
  geom_bar()
ggsave(out("Samples_Numbers.pdf"), w=0.5 * length(unique(ann$sample)) + 1, h=6)

ggplot(ann, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~sample, ncol=5)
ggsave(out("Samples_UMAP.pdf"), w=5*2+2,h=ceiling(length(unique(ann$sample))/5) * 2 + 1)

pDT <- ann[,.N, by=c("sample", "Clusters")]
pDT[,sumS := sum(N), by="sample"]
pDT[,sumC := sum(N), by="Clusters"]
pDT[,percentS := N/sumS*100]
pDT[,percentC := N/sumC*100]
ggplot(pDT, aes(y=sample,x=factor(as.numeric(Clusters)), size=percentS, color=percentC)) +
  scale_size_continuous(name="% of sample") + 
  scale_color_gradient(name="% of cluster", low="black", high="red") + 
  theme_bw() + 
  geom_point()
ggsave(out("Samples_Clusters.pdf"), w=8,h=length(unique(pDT$sample)) * 0.3+1)


# TISSUES -----------------------------------------------------------------
ggplot(ann, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~tissue, ncol=3)
ggsave(out("Tissues_UMAP.pdf"), w=3*2+2,h=3 + 1)

pDT <- ann[,.N, by=c("tissue", "Clusters")]
pDT[,sumS := sum(N), by="tissue"]
pDT[,sumC := sum(N), by="Clusters"]
pDT[,percentS := N/sumS*100]
pDT[,percentC := N/sumC*100]
ggplot(pDT, aes(y=tissue,x=factor(as.numeric(Clusters)), size=percentS, color=percentC)) +
  scale_size_continuous(name="% of tissue") + 
  scale_color_gradient(name="% of cluster", low="black", high="red") + 
  theme_bw() + 
  geom_point()
ggsave(out("Tissues_Clusters.pdf"), w=8,h=length(unique(pDT$tissue)) * 1+1)


# Cell cycle --------------------------------------------------------------
ggplot(ann, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_hexbin() +
  facet_grid(. ~ Phase)
ggsave(out("CellCycle_UMAP.pdf"), w=16,h=5)



# CLUSTERS ----------------------------------------------------
# UMAP
ggplot(ann, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  geom_label(data=ann[,.(UMAP1=median(UMAP1), UMAP2=median(UMAP2)), by="Clusters"], aes(label=Clusters), fill="#ffffffaa")
ggsave(out("Clusters_UMAP.pdf"), w=6,h=5)

plot_cells(monocle.obj,
           reduction_method="UMAP",
           show_trajectory_graph=FALSE, 
           color_cells_by="cluster")
ggsave(out("Clusters_UMAP.points.jpg"), w=5,h=5)


# Cell type MARKERS  ------------------------------------------------------

# Markers
plot_genes_by_group(monocle.obj, markers = marker.genes$Name, group_cells_by = "cluster") + scale_size_continuous(range=c(0,5))
ggsave(out("Markers_Clusters.pdf"), w=6,h=8)

p <- plot_cells(monocle.obj,
           reduction_method="UMAP",
           genes = sort(marker.genes$Name),
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE
           )
ggsave(out("Markers_UMAP.jpg"), w=30,h=30)

# gg.order <- t(as.matrix(counts(monocle.obj[marker.genes$Name,])))
# gg.order <- corS(mm)
# gg.order <- hclust(as.dist(1-cMT))
# gg.order <- gg.order$labels[gg.order$order]

p <- plot_cells_umap_hex_NF(monocle.obj,genes = sort(marker.genes$Name))
ggsave(out("Markers_UMAP_hex.pdf"), w=30,h=30, plot=p)

p <- plot_cells_umap_hex_NF(monocle.obj, scale=TRUE, genes = sort(marker.genes$Name))
ggsave(out("Markers_UMAP_hex_scale.pdf"), w=30,h=30, plot=p)


# CELLTYPES SingleR -------------------------------------------------------
singleR.sum <- data.table()
srx <- names(singleR.res)[1]
for(srx in names(singleR.res)){
  print(srx)
  singleR.resX <- singleR.res[[srx]]
  
  # Predicted cells on UMAP
  pDT.pc <- merge(
    melt(singleR.resX[,c("cell", "labels"), with=F], id.vars = "cell"),
    ann[,c("UMAP1", "UMAP2", "rn", "Clusters"),with=F],
    by.x="cell", by.y="rn")
  pDT.pc <- pDT.pc[value %in% pDT.pc[,.N, by="value"][N > 10]$value]
  p <- ggplot(pDT.pc, aes(x=UMAP1, y=UMAP2)) +
    theme_bw(12) +
    geom_hex(bins=100) + 
    scale_fill_hexbin() +
    facet_wrap(~value, ncol=5) +
    theme_bw(12) +
    ggtitle(srx)
  ggsave(
    out("SingleR_", srx, "_UMAP_Predicted",".pdf"), 
    w=5*3+2, 
    h=ceiling(length(unique(pDT.pc$value))/5)*3+1, 
    limitsize = FALSE,
    plot=p)
  
  # Clusters - Predictions
  pDT.ann <- merge(singleR.resX[,c("labels", "cell")],ann, by.x="cell", by.y="rn")
  pDT.ann <- pDT.ann[,.N, by=c("Clusters", "labels", "tissue")]
  pDT.ann[,sum := sum(N), by=c("Clusters", "tissue")]
  pDT.ann[,percent := N/sum*100]
  ggplot(pDT.ann, aes(x=factor(as.numeric(Clusters)), y=labels, fill=percent)) + 
    theme_bw(12) +
    geom_tile() +
    facet_grid(. ~ tissue) +
    scale_fill_gradient(limits=c(0,100), low="white", high="red") +
    ggtitle(srx)
  ggsave(out("SingleR_", srx, "_Clusters_", "PercPredicted", ".pdf"),          
         h=length(unique(pDT.ann$labels)) * 0.3+1,
         w=length(unique(pDT.ann$Clusters)) * 0.3 * 3+2,
         limitsize = FALSE)
  
  pDT.ann$dataset <- srx
  singleR.sum <- rbind(singleR.sum, pDT.ann)
}
singleR.sum[,id := paste(dataset, labels, tissue)]
pDT <- merge(singleR.sum[,max(percent), by=c("id")][V1 > 20][,c("id")], singleR.sum, by=c("id"))
pDT <- hierarch.ordering(pDT, toOrder = "Clusters", orderBy = "labels", value.var = "percent", aggregate = TRUE)
pDT <- hierarch.ordering(pDT, toOrder = "labels", orderBy = "Clusters", value.var = "percent", aggregate = TRUE)
ggplot(pDT, aes(y=Clusters, x=labels, fill=percent)) + 
  theme_bw(12) + 
  geom_tile() +
  facet_grid(tissue ~ gsub("_", "\n", dataset), scales = "free", space = "free") + 
  scale_fill_gradient(limits=c(0,100), low="white", high="red") +
  xRot()
ggsave(out("SingleR_0_Clusters_", "PercPredicted", ".pdf"),          
       w=nrow(pDT[,.N, by=c("dataset", "labels")]) * 0.2+2,
       h=length(unique(pDT$Clusters)) * 0.2 * 3+1,
       limitsize = FALSE)



# CellTypes from Marker signatures --------------------------------------------------
mnam <- names(marker.signatures)[4]
for(mnam in names(marker.signatures)){
  mx <- marker.signatures[[mnam]]
  
  pDT <- merge(ann[,c("rn", "UMAP1", "UMAP2")], melt(data.table(mx, keep.rownames = TRUE), id.vars = "rn"), by="rn")
  pDT[, value.norm := scale(value), by="variable"]
  # filter what to show?
  # pDT[value.norm > 2][,.(UMAP1 = sd(UMAP1), UMAP2 = sd(UMAP2)), by="variable"][,.(mean(UMAP1, UMAP2)), by="variable"][order(V1)]
  
  ggplot(pDT, aes(x=UMAP1, y=UMAP2)) +
    stat_summary_hex(aes(z=value),fun=mean, bins=100) +
    #scale_fill_gradient2(low="blue", midpoint = 0, high="red") +
    scale_fill_hexbin() +
    theme_bw(12) +
    facet_wrap(~variable) +
    ggtitle("Marker Signatures - Larry et al, Science")
  ggsave(out("Markers_Signatures_",mnam,"_UMAP_raw.pdf"), w=12,h=12)
  
  ggplot(pDT, aes(x=UMAP1, y=UMAP2)) +
    stat_summary_hex(aes(z=pmin(3, value.norm)),fun=mean, bins=100) +
    scale_fill_gradient2(low="blue", midpoint = 0, high="red") +
    #scale_fill_hexbin() +
    theme_bw(12) +
    facet_wrap(~variable) +
    ggtitle("Marker Signatures - Larry et al, Science")
  ggsave(out("Markers_Signatures_",mnam,"_UMAP_scaled.pdf"), w=12,h=12)
  
  cleanDev(); pdf(out("Markers_Signatures_",mnam,"_Clusters.pdf"),w=8,h=6)
  pheatmap(sapply(with(ann, split(rn, Clusters)), function(cx) colMeans(mx[cx,])))
  dev.off()
}



# CYTOTRACE ---------------------------------------------------------------
if("cytoRes" %in% ls()){
  # Umap
  ggplot(ann, aes(x=UMAP1, y=UMAP2)) +
    stat_summary_hex(aes(z=CytoTRACE),fun=mean, bins=100) +
    #scale_fill_gradient2(low="blue", midpoint = 0.5, high="red") +
    scale_fill_hexbin() +
    theme_bw(12) +
    ggtitle("CytoTRACE - (1: less diff; 0: more diff)")
  ggsave(out("CytoTRACE_UMAP.pdf"), w=5,h=4)

  # Clusters
  ggplot(ann, aes(x=factor(as.numeric(Clusters)), y=CytoTRACE)) +
    geom_violin(color=NA, fill="lightblue") +
    geom_boxplot(fill=NA, coef=Inf) +
    theme_bw(12) +
    xRot() +
    ggtitle("CytoTRACE - (1: less diff; 0: more diff)")
  ggsave(out("CytoTRACE_Clusters.pdf"), w=7,h=4)
}



# ANTIBODIES --------------------------------------------------------------
if("CITESEQ2" %in% ann$sample){
  abMT <- additional.info$CITESEQ2$`Antibody Capture`
  abMT <- SCRNA.TPXToLog(SCRNA.RawToTPX(abMT, scale.factor = 1e6))
  ann.c1 <- ann[sample == "CITESEQ2"]
  
  # UMAP
  res <- data.table()
  abx <- row.names(abMT)[1]
  for(abx in row.names(abMT)){
    pDT <- copy(ann.c1)
    pDT$Signal <- abMT[abx, gsub("_CITESEQ2", "", ann.c1$rn)]
    pDT$Antibody <- abx
    res <- rbind(res, pDT)
  }
  res[,Signal.norm := scale(Signal), by="Antibody"]
  ggplot(res, aes(x=UMAP1, y=UMAP2)) +
    stat_summary_hex(bins = 100, aes(z=Signal.norm),fun=mean) +
    scale_fill_gradient2(low="blue", high="red") +
    facet_wrap(~Antibody) +
    theme_bw(12)
  ggsave(out("Antibodies_UMAP.pdf"), w=12+2, h=9+1)
  
  # Correlations
  cMT <- corS(t(as.matrix(abMT)), use="pairwise.complete.obs")
  diag(cMT) <- NA
  dist <- as.dist(1-cMT)
  cleanDev(); pdf(out("Antibodies_Correlation.pdf"), w=5,h=4)
  pheatmap(cMT,
           clustering_distance_rows = dist,
           clustering_distance_cols = dist,
           breaks=seq(-1,1, 0.01), color=COLORS.HM.FUNC(200),
  )
  dev.off()
  
  # Percentile plots
  probx=0.9
  res[,percentile := quantile(Signal.norm, probs=probx, na.rm=TRUE), by="Antibody"]
  resN <- res[, sum(Signal.norm > percentile), by=c("Clusters", "Antibody")]
  resN[,clSize := sum(V1), by="Clusters"]
  stopifnot(all(resN[, length(unique(clSize)), by="Clusters"]$V1 == 1))
  resN[,percentage := V1/clSize*100]
  resN <- hierarch.ordering(resN, toOrder = "Clusters", orderBy = "Antibody", value.var = "percentage")
  resN <- hierarch.ordering(resN, toOrder = "Antibody", orderBy = "Clusters", value.var = "percentage")
  resN[is.na(percentage), percentage := 0]
  ggplot(resN, aes(x=Clusters, y=Antibody, fill=percentage)) +
    geom_tile() +
    ggtitle("Percent of cell in Clusters\nthat are above 90th percentile of antibody signal") +
    scale_fill_gradient(low="white", high="blue")
  ggsave(out("Antibodies_Percentile.pdf"), w=5, h=4)
}


# GUIDES ------------------------------------------------------

# . Number of guides --------------------------------------------------------
# Number of cells per guide / mixscape
ggplot(ann[!is.na(mixscape_class.global)], aes(x=guide, fill=mixscape_class.global)) + 
  geom_bar(position="dodge") + 
  facet_wrap(~sample, scales = "free") +
  theme_bw(12) + xRot() +
  ylab("Cells")
ggsave(out("Guides_Counts.pdf"), w=12,h=10)
# % of cells assigned
ggplot(ann[,sum(!is.na(guide))/.N*100, by="sample"], aes(x=sample, y=V1)) + 
  geom_bar(stat="identity") + 
  theme_bw(12) + xRot() +
  ylab("Percent of cells with assigned guide (NP or KO)")
ggsave(out("Guides_Counts_PercAssigned.pdf"), w=6,h=4)

# . Guides in UMAP ----------------------------------------------------------
sx <- "ECCITE1"
for(sx in unique(ann[!is.na(mixscape_class)]$sample)){
  pDT <- ann[sample == sx][!is.na(mixscape_class.global)]
  if(nrow(pDT) == 0) next
  grps <- length(unique(pDT$mixscape_class))
  ggplot(pDT, aes(x=UMAP1, y=UMAP2)) + 
    geom_hex(bins=50) +
    scale_fill_gradient(low="lightgrey", high="blue") +
    facet_wrap(~mixscape_class + mixscape_class.global, ncol = 7) +
    theme_bw(12)
  ggsave(out("Guides_UMAP_", sx, ".pdf"), w=7*4,h=ceiling(grps/7)*4+1)
}


# . Guides per cluster - MIXSCAPE ------------------------------------------------------
res <- data.table()
pDT1 <- ann[!is.na(mixscape_class)][mixscape_class.global != "NP"]
for(sx in unique(pDT1$sample)){
  pDT2 <- pDT1[sample == sx]
  for(gx in unique(pDT2[mixscape_class != "NTC"]$mixscape_class)){
    for(cx in unique(pDT2$Clusters)){
      #message(gx, "-", cx)
      mx <- as.matrix(with(pDT2[mixscape_class %in% c(gx, "NTC")], table(Clusters == cx, mixscape_class == gx)))
      #print(mx)
      if(dim(mx) == c(2,2)){
        fish <- fisher.test(mx)
        res <- rbind(res, data.table(Clusters=cx, mixscape_class=gx, p=fish$p.value, OR=fish$estimate, sample=sx, total.cells=sum(mx)))
      }
    }
  }
}
res[,padj := p.adjust(p, method="BH")]
res[, log2OR := log2(pmin(5, OR + min(res[OR != 0]$OR)))]
res[,grp := paste(mixscape_class, sample)]
res <- hierarch.ordering(res, toOrder = "grp", orderBy = "Clusters", value.var = "log2OR")
#res <- hierarch.ordering(res, toOrder = "Clusters", orderBy = "grp", value.var = "log2OR")
ggplot(res, aes(
  x=Clusters,
  y=sample, 
  color=log2OR, 
  size=pmin(-log10(padj), 5))) + 
  geom_point(shape=16) +
  scale_color_gradient2(name="log2OR", low="blue", high="red") +
  scale_size_continuous(name="padj") + 
  facet_grid(mixscape_class ~ ., space = "free", scales = "free") +
  theme_bw(12) +
  theme(strip.text.y = element_text(angle=0))
ggsave(out("Guides_Fisher_Mixscape.pdf"), w=10, h=length(unique(res$grp)) * 0.25 + 1, limitsize = FALSE)
write.tsv(res[,-c("grp"), with=F], out("Guides_Fisher_Mixscape.tsv"))


# . Guides per cluster ------------------------------------------------------
res <- data.table()
pDT1 <- ann[!is.na(guide)]
for(sx in unique(pDT1$sample)){
  pDT2 <- pDT1[sample == sx]
  for(gx in unique(pDT2[mixscape_class != "NTC"]$guide)){
    for(cx in unique(pDT2$Clusters)){
      #message(gx, "-", cx)
      mx <- as.matrix(with(pDT2[guide %in% c(gx, "NTC")], table(Clusters == cx, guide == gx)))
      #print(mx)
      if(dim(mx) == c(2,2)){
        fish <- fisher.test(mx)
        res <- rbind(res, data.table(Clusters=cx, guide=gx, p=fish$p.value, OR=fish$estimate, sample=sx, total.cells=sum(mx)))
      }
    }
  }
}
res[,padj := p.adjust(p, method="BH")]
res[, log2OR := log2(pmin(5, OR + min(res[OR != 0]$OR)))]
res[,grp := paste(guide, sample)]
res <- hierarch.ordering(res, toOrder = "grp", orderBy = "Clusters", value.var = "log2OR")
#res <- hierarch.ordering(res, toOrder = "Clusters", orderBy = "grp", value.var = "log2OR")
ggplot(res, aes(
  x=Clusters,
  y=sample,
  color=log2OR, 
  size=pmin(-log10(padj), 5))) + 
  geom_point(shape=16) +
  scale_color_gradient2(name="log2OR", low="blue", high="red") +
  scale_size_continuous(name="padj") + 
  facet_grid(guide ~ ., space = "free", scales = "free") +
  theme_bw(12) +
  theme(strip.text.y = element_text(angle=0))
ggsave(out("Guides_Fisher_noMixscape.pdf"), w=10, h=length(unique(res$grp)) * 0.25 + 1, limitsize = FALSE)
write.tsv(res[,-c("grp"), with=F], out("Guides_Fisher_noMixscape.tsv"))



# Guides Signature differential analysis ----------------------------------

# Calculate stats
mMT <- cbind(marker.signatures$Larry, marker.signatures$PanglaoDB[,c("Basophils", "Megakaryocytes")])
resSDA <- data.table()
guides <- unique(ann[mixscape_class.global == "KO"][,.N, by="guide"][N > 10]$guide)
sigx <- colnames(mMT)[1]
for(sigx in colnames(mMT)){
  xNTC <- mMT[ann[mixscape_class.global == "NTC"]$rn, sigx]
  (guidex <- ann$guide[1])
  for(guidex in guides){
    x <- mMT[ann[guide == guidex & mixscape_class.global == "KO"]$rn, sigx]
    resSDA <- rbind(resSDA, data.table(
      guide=guidex, 
      sig=sigx, 
      p=wilcox.test(x, xNTC)$p.value, 
      d=mean(x) - mean(xNTC)
      ))
  }
}
resSDA[, padj := p.adjust(p, method="BH")]

# Plot stats
ggplot(resSDA, aes(x=guide,y=sig, color=d, size=pmin(5, -log10(padj)))) + 
  theme_bw(12) +
  geom_point() +
  scale_color_gradient2(low="blue", high="red") +
  xRot()
ggsave(out("SigDA_Stats.pdf"), w=8,h=5)

# Guides in UMAP
pDT <- ann[mixscape_class.global %in% c("KO", "NTC") & guide %in% c(guides, "NTC")]
grps <- length(unique(pDT$guide))
ggplot(pDT, aes(x=UMAP1, y=UMAP2)) + 
  geom_hex(bins=50) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~guide, ncol = 7) +
  theme_bw(12)
ggsave(out("SigDA_UMAP.pdf"), w=7*4,h=ceiling(grps/7)*4+1)

# follow up specific guide
# xDT <- copy(ann)[guide %in% c("Kmt2d_4G", "NTC") & mixscape_class.global %in% c("KO", "NTC")]
# xDT$Sig <- marker.signatures.use[xDT$rn,"Granulocyte"]
# ggplot(xDT, aes(x=UMAP1, y=UMAP2, color=Sig)) + 
#   geom_point(size=0.2) +
#   scale_color_gradient(high="red") +
#   scale_fill_gradient(low="lightgrey", high="blue") +
#   facet_wrap(~guide, ncol = 7) +
#   theme_bw(12)


# DE for GUIDES -----------------------------------------
obj.de <- monocle.obj

# Keep only KO and NTCs
obj.de <- obj.de[,obj.de$mixscape_class.global %in% c("KO", "NTC")]
stopifnot(all(names(monocle3::clusters(obj.de)) == colnames(obj.de)))
obj.de <- obj.de[,getCL(obj.de) %in% names(which(table(getCL(obj.de)) > 30))]

# Remove lowly expressed genes
obj.de <- obj.de[Matrix::rowSums(counts(obj.de)) > 20,]

# Order by sample
obj.de <- obj.de[,order(obj.de$sample)]
# ggplot(data.table(sample=obj.de$sample, i=1:ncol(obj.de)), aes(x=i, y=sample)) + geom_point()

# prepare annotation for DE
obj.de.ann <- data.frame(
  row.names=colnames(obj.de),
  GuideDE=gsub("_.+$", "", obj.de$guide),
  ClusterDE=getCL(obj.de),
  tissueDE=obj.de$tissue
)
obj.de.ann <- filter(obj.de.ann, !GuideDE %in% c("Pu.1", "Spi1"))
obj.de.ann <- mutate(obj.de.ann, GuideDE = gsub("^Men$", "Men1", GuideDE))
sort(unique(obj.de.ann$GuideDE))
obj.de.ann$tissueDE <- factor(obj.de.ann$tissueDE, levels=c("in vitro", "leukemia", "in vivo"))
obj.de.ann$GuideDE <- factor(obj.de.ann$GuideDE, levels=c("NTC", setdiff(unique(obj.de.ann$GuideDE), "NTC")))
obj.de <- obj.de[, row.names(obj.de.ann)]
stopifnot(row.names(obj.de.ann) == colnames(obj.de))
str(obj.de.ann)
write.tsv(data.table(obj.de.ann, keep.rownames = TRUE), out("DEG_Annnotation.tsv"))

# Run nebula
if(file.exists(neb.file)){
  (load(neb.file))
} else {
  # Model matrix
  if(baseDir.add != "combined"){ # for the combined one we look at tissue differences
    mm <- model.matrix(data=obj.de.ann, ~ GuideDE + ClusterDE)
    colnames(mm) <- gsub("^(GuideDE.+)$", paste0("\\1_", baseDir.add), colnames(mm))
  } else {
    mm <- model.matrix(data=obj.de.ann, ~ GuideDE * tissueDE + ClusterDE)
    x <- data.table(obj.de.ann)
    for(i in 1:ncol(x)) x[[i]] <- as.character(x[[i]])
    gx <- x[GuideDE != "NTC"]$GuideDE[1]
    
    # fix coefficients
    for(gx in unique(obj.de.ann$GuideDE)){
      tx <- unique(x[GuideDE == gx]$tissueDE)
      if("in vitro" %in% tx){
        # Rename the main guide effect to "in vitro
        colnames(mm) <- gsub(paste0("^(","GuideDE",gx,")$"), "\\1_in vitro", colnames(mm))
      } else{
        # Remove current main guide effect (for in vitro)
        mm <- mm[,-which(colnames(mm) == paste0("GuideDE",gx))]
        # add main effects for leukemia and in vivo
        if("leukemia" %in% tx) colnames(mm)[which(colnames(mm) == paste0("GuideDE", gx, ":", "tissueDE", "leukemia"))] <- paste0("GuideDE",gx,"_leukemia")
        if("in vivo" %in% tx)  colnames(mm)[which(colnames(mm) == paste0("GuideDE", gx, ":", "tissueDE", "in vivo"))] <- paste0("GuideDE",gx,"_in vivo")
      }
    }
  }
  
  mm <- mm[,colSums(mm) != 0]
  x <- as.matrix(unique(data.table(mm)))
  #pheatmap(x[,sort(colnames(x))], cluster_cols = F)
  
  # run nebula
  nebRes <- nebula(
    count = counts(obj.de),
    id = obj.de$sample,
    pred = mm
  )
  
  # export results
  res <- data.table()
  for(cx in colnames(mm)){
    res <- rbind(res, data.table(
      term=cx,
      p_value=nebRes$summary[,paste("p", cx, sep="_")],
      se=nebRes$summary[,paste("se", cx, sep="_")],
      estimate=nebRes$summary[,paste("logFC", cx, sep="_")],
      gene_id=nebRes$summary$gene,
      convergence=nebRes$convergence
    ))
  }
  res[,q_value := p.adjust(p_value, method="BH")]
  res[, estimate := estimate / log(2)] # convert to log2 FC
  save(res, file=neb.file)
}


#  . Export results -------------------------------------------------------
resGuides <- res[grepl("GuideDE", term)][convergence >= -15]
table(res$convergence)
resGuides[!grepl("tissueDE", term), tissue := gsub("^.+?_", "", term)]
resGuides[grepl("tissueDE", term), tissue := gsub("^.+\\:tissueDE", "", term)]
resGuides[, interaction := grepl("tissueDE", term)]
resGuides[, guide := gsub("^GuideDE", "", term)]
resGuides[, guide := gsub("\\:tissueDE.+$", "", guide)]
resGuides[, guide := gsub("_.+$", "", guide)]
resGuides <- resGuides[!is.na(estimate)]
resGuides[, estimate_raw := estimate]
resGuides[, estimate := ifelse(p_value > 0.9, 0, estimate)]
write.tsv(resGuides[q_value < 1][,-"term",with=F], file=out("DEG_Results_Export.tsv"))
write.tsv(resGuides, file=out("DEG_Results_all.tsv"))

#  . Vulcano / p-val distribution -----------------------------------------
ggplot(resGuides, aes(x=estimate, y=-log10(p_value))) + 
  theme_bw(12) +
  geom_hex(bins=100) +
  facet_wrap(~gsub("_", "\n", gsub("\\:tissueDE", "\nInteraction: ", term)), scales = "free", ncol = 5)
ggsave(out("DEG_Vulcano.pdf"), w=15,h=15)

ggplot(resGuides, aes(x=p_value)) + 
  theme_bw(12) +
  geom_histogram() +
  facet_wrap(~gsub("_", "\n", gsub("\\:tissueDE", "\nInteraction: ", term)), scales = "free", ncol = 5)
ggsave(out("DEG_PVal_histogram.pdf"), w=15,h=15)


if(length(unique(resGuides$term)) < 2){
  stop("Not proceeding as only one or no gene produced DEG results. Guides: ", paste(unique(resGuides$term), collapse=";"))
}

#  . Plot top genes -------------------------------------------------------
# Estimate and P-value
gg <- resGuides[q_value < 0.05][order(-abs(estimate))][,head(.SD,n=10), by="term"]$gene_id
pDT <- resGuides[gene_id %in% gg]
pDT <- hierarch.ordering(pDT, toOrder="gene_id", orderBy = "term", value.var = "estimate")
pDT <- hierarch.ordering(pDT, toOrder="term", orderBy = "gene_id", value.var = "estimate")
ggplot(pDT, aes(x=guide, y=gene_id, size=pmin(5, -log10(q_value)), color=sign(estimate) * pmin(5,abs(estimate)))) +
  scale_size_continuous(name="padj", range=c(0,5)) +
  scale_color_gradient2(name="delta", low="blue", high="red") +
  facet_grid(. ~ paste("Interaction:", interaction) + tissue, space = "free", scales = "free") + 
  theme_bw(12) +
  geom_point() +
  xRot()
ggsave(out("DEG_examples.pdf"), w=10,h=30)



# UMAP and correlations of DEG --------------------------------------------
resGuides <- fread(out("DEG_Results_all.tsv"))


# . LogFC MATRIX ------------------------------------------------------------
resGuides.I <- merge(resGuides[interaction==TRUE], resGuides[interaction==FALSE, c("guide", "estimate", "gene_id"),with=F], by=c("guide", "gene_id"), all.x=TRUE, allow.cartesian=FALSE)
stopifnot(sum(is.na(resGuides.I$estimate.x)) == 0)
stopifnot(sum(is.na(resGuides.I$estimate.y)) == 0)
resGuides.I[,estimate := estimate.y + estimate.x]
resGuides.I <- rbind(
  resGuides.I[,c("guide", "estimate", "gene_id", "tissue"),with=F], 
  resGuides[interaction==FALSE,c("guide", "estimate", "gene_id", "tissue"),with=F],
  fill=TRUE)
resGuides.I[,id := paste(guide, tissue)]
umapMT <- toMT(resGuides.I, row = "gene_id", col = "id", val = "estimate")
write.table(umapMT, sep=",", row.names = TRUE, col.names = TRUE, quote = FALSE, file=out("DEG_Results_logFCMT.csv"))


# . CF ChIP-seq targets -----------------------------------------------------

# Look for ChIP targets based on logFC
chip.distributions <- data.table()
for(fx in names(chip.targets)){
  x <- copy(resGuides.I)
  x[,ChIP := fx]
  x[,ChIP.target := gene_id %in% chip.targets[[fx]]]
  chip.distributions <- rbind(chip.distributions, x)
}
chip.distributions[,ChIP.guide := gsub("^.+_", "", ChIP)]
pDT <- chip.distributions[ChIP.guide == guide]
pDT[,id := paste("Targets:\n", sub(" ", "\n", id))]
pDT[,ChIP := paste("ChIP:\n", sub("_", "\n", ChIP))]
p <- ggplot(pDT, aes(x=ChIP, color=ChIP.target, y=estimate)) +
  theme_bw(12) +
  geom_violin() +
  facet_wrap(~ id, scales = "free") +
  xRot()
ggsave(out("DEG_ChIP.pdf"), w=15,h=15, plot=p)

# Calculate enrichment
idx <- "Kmt2a in vitro"
chip.gsea.file <- out("DEG_ChIP_GSEA.tsv")
if(file.exists(chip.gsea.file)){
  chip.gsea.res <- fread(chip.gsea.file)
} else {
  chip.gsea.res <- data.table()
  for(idx in unique(resGuides.I$id)){
    rnaDT <- resGuides.I[id == idx]
    quantile(rnaDT$estimate)
    chipL <- chip.targets[grepl(rnaDT$guide[1], names(chip.targets))]
    res <- fgsea(chipL, stats=setNames(rnaDT$estimate, rnaDT$gene_id))
    res$list <- idx
    chip.gsea.res <- rbind(chip.gsea.res, res)
  }
  chip.gsea.res$leadingEdge <- NULL #sapply(chip.gsea.res$leadingEdge, function(x) paste(x, collapse = ","))
  write.tsv(chip.gsea.res, chip.gsea.file)
}
save(chip.gsea.res, resGuides.I, chip.targets, file=out("DEG_ChIP_GSEA.RData"))
ggplot(chip.gsea.res, aes(x=pathway, y=list, size=pmin(5, -log10(padj)), color=NES))+ 
  theme_bw(12) +
  scale_color_gradient2(low="blue", high="red") +
  geom_point() + 
  xRot()
ggsave(out("DEG_ChIP_GSEA.pdf"), w=9,h=6)


# Correlation
chip.fc <- copy(chip.targets.fc)
chip.fc[,ChIP.guide := gsub("^.+_", "", ChIP)]
chip.fc <- merge(resGuides.I, chip.fc, by.x=c("guide", "gene_id"), by.y=c("ChIP.guide", "Gene"))
chip.fc[,id := paste("Targets:\n", sub(" ", "\n", id))]
chip.fc[,ChIP := paste("ChIP:\n", sub("_", "\n", ChIP))]
chip.fc$FC <- unlist(chip.fc$FC)
p <- ggplot(chip.fc, aes(y=FC,x=estimate)) + 
  theme_bw() +
  stat_binhex(aes(fill=log10(..count..))) +
  facet_wrap(~ id + ChIP, scales = "free")
ggsave(out("DEG_ChIP_Correlation.pdf"), w=30,h=30, plot=p)


# . COR of DEG --------------------------------------------------------------
cMT <- corS(umapMT)
gn <- ncol(umapMT)
dd <- as.dist(1-cMT)
diag(cMT) <- NA
cleanDev(); pdf(out("DEG_CorHM.pdf"),w=gn/6+2, h=gn/6+1.5)
pheatmap(cMT,clustering_distance_rows = dd, clustering_distance_cols = dd)
dev.off()

cleanDev(); pdf(out("DEG_CorHM_Colors.pdf"),w=gn/6+2, h=gn/6+1.5)
pheatmap(cMT,
         clustering_distance_rows = dd, clustering_distance_cols = dd,
         breaks=seq(-1,1, 0.01),
         color=colorRampPalette(c("#6a3d9a", "#a6cee3", "white", "#fdbf6f", "#e31a1c"))(200))
dev.off()



# . UMAP of DEG -----------------------------------------
umap.log2FC.cutoff <- 3
for(umap.type in c("all", "top")){
  
  umap.type.name <- paste0(umap.type, ".genes")
  
  set.seed(1212)
  umap.file <- out("RegulatoryMap_UMAP_",umap.type.name,".tsv")
  if(file.exists(umap.file)){
    umap <- fread(umap.file)
  } else {
    gg <- if(umap.type == "top") unique(resGuides[q_value < 0.05 & abs(estimate) > 1]$gene_id) else row.names(umapMT)
    mt <- umapMT[gg,]
    mt[mt > umap.log2FC.cutoff] <- umap.log2FC.cutoff
    mt[mt < -umap.log2FC.cutoff] <- -umap.log2FC.cutoff
    umap.res <- umap(mt)
    umap <- data.table(umap.res$layout, keep.rownames = TRUE)
    umap <- setNames(umap, c("Gene", "UMAP1", "UMAP2"))
    ggplot(umap, aes(x=UMAP1, y=UMAP2)) + geom_hex(bins=100) + theme_bw(12)
    ggsave(out("RegulatoryMap_UMAP_",umap.type.name,".pdf"), w=6,h=5)
    
    # Cluster
    set.seed(1212)
    idx <- umap.res$knn$indexes
    g <- do.call(rbind, apply(idx[, 2:ncol(idx)], 2, function(col){data.table(row.names(idx)[col], row.names(idx)[idx[,1]])}))
    (g <- graph.edgelist(as.matrix(g),directed=FALSE))
    set.seed(1234)
    cl <- cluster_walktrap(g)
    clx <- setNames(cl$membership, V(g)$name)
    umap$Cluster <- clx[umap$Gene]
    
    ggplot(umap, aes(x=UMAP1, y=UMAP2, color=factor(Cluster))) + 
      geom_point() + 
      theme_bw(12) +
      geom_label(data=umap[, .(UMAP1=median(UMAP1), UMAP2=median(UMAP2)), by="Cluster"], aes(label=Cluster))
    ggsave(out("RegulatoryMap_UMAP_",umap.type.name,"_Clusters.pdf"), w=6,h=5)
    
    # Export annotation
    write.tsv(umap, umap.file)
  }
  
  # Prepare for plotting
  pUMAP.de <- merge(umap, setNames(melt(data.table(umapMT, keep.rownames = TRUE), id.vars = "rn"), c("gene_id", "term", "estimate")), by.x="Gene", by.y="gene_id")
  dim.umap1 <- floor(max(abs(pUMAP.de$UMAP1))) + 0.5
  dim.umap2 <- floor(max(abs(pUMAP.de$UMAP2))) + 0.5
  pUMAP.de[, guide := gsub(" .+", "", term)]
  pUMAP.de[, tissue := gsub(".+? ", "", term)]
  tn <- length(unique(pUMAP.de$guide))
  pUMAP.de[, estimate_cap := pmin(umap.log2FC.cutoff, abs(estimate)) * sign(estimate)]
  
  # Plot estimates on UMAP
  ggplot(pUMAP.de, aes(x=UMAP1, y=UMAP2)) +
    stat_summary_hex(
      aes(z=estimate_cap),
      fun=mean) +
    scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") +
    facet_grid(guide~tissue) + theme_bw(12) +
    xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
    xlim(-dim.umap1,dim.umap1) + ylim(-dim.umap2,dim.umap2)
  ggsave(out("RegulatoryMap_UMAP_",umap.type.name,"_Values.pdf"), w=2*3+2,h=tn * 2 + 1, limitsize=FALSE)
  
  # values by cluster
  pDT <- pUMAP.de[, mean(estimate_cap), by=c("Cluster", "term")]
  pDT <- hierarch.ordering(pDT, toOrder = "Cluster", orderBy = "term", value.var = "V1")
  pDT <- hierarch.ordering(pDT, orderBy = "Cluster", toOrder = "term", value.var = "V1")
  ggplot(pDT, aes(x=factor(Cluster), y=term, fill=V1)) + 
    theme_bw(12) + 
    geom_tile() +
    scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") +
    xlab("Gene modules (Gene-UMAP Clusters)")
  ggsave(out("RegulatoryMap_UMAP_",umap.type.name,"_ClusterValues.pdf"), w=10,h=tn * 0.2 + 1, limitsize=FALSE)
  
  
  # CF targets on UMAP
  # umap <- fread(out("RegulatoryMap_UMAP.tsv"))
  pDT <- data.table()
  for(cf in names(chip.targets)){
    x <- copy(umap)
    x[,Target := Gene %in% chip.targets[[cf]] + 0]
    x$CF <- cf
    pDT <- rbind(pDT, x)
  }
  pDT[, tissue := gsub("_.+$", "", CF)] 
  pDT[, factor := gsub("^.+?_", "", CF)]
  # dimensions for umap
  dim.umap1 <- floor(max(abs(pDT$UMAP1))) + 0.5
  dim.umap2 <- floor(max(abs(pDT$UMAP2))) + 0.5
  # UMAP percentage
  hex.percent <- function(x){sum(x)/length(x) * 100}
  ggplot(pDT, aes(x=UMAP1, y=UMAP2)) +
    stat_summary_hex(
      aes(z=Target),
      fun=hex.percent) +
    scale_fill_gradient(high="#e31a1c", low="#ffffff") +
    facet_grid(factor~tissue) + theme_bw(12) +
    xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
    xlim(-dim.umap1,dim.umap1) + ylim(-dim.umap2,dim.umap2)
  ggsave(out("RegulatoryMap_UMAP_",umap.type.name,"_CFtargets.pdf"), 
         w=length(unique(pDT$tissue))*1+2,
         h=length(unique(pDT$factor))*1+1,
         limitsize=FALSE)
  
  # UMAP normalized
  df <- pDT[CF == pDT$CF[1]]
  makeHexData <- function(df) {
    h <- hexbin::hexbin(df$UMAP1, df$UMAP2, IDs = TRUE)
    ret <- data.frame(hexbin::hcell2xy(h),
                      percent = tapply(df$Target, h@cID, FUN = function(z) sum(z)/length(z)*100),
                      cid = h@cell)
    ret <- mutate(ret, norm = percent / max(ret$percent))
    ret
  }
  pDT.hex <- do.call(rbind, lapply(unique(pDT$CF), function(cfx){
    data.table(makeHexData(pDT[CF == cfx]), CF=cfx)
  }))
  pDT.hex[, tissue := gsub("_.+$", "", CF)] 
  pDT.hex[, factor := gsub("^.+?_", "", CF)]
  ggplot(pDT.hex) +
    geom_hex(aes(x = x, y = y, fill = norm), stat = "identity") +
    scale_fill_gradient(high="#e31a1c", low="#ffffff") +
    facet_grid(factor~tissue) + theme_bw(12) +
    xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
    xlim(-dim.umap1,dim.umap1) + ylim(-dim.umap2,dim.umap2)
  ggsave(out("RegulatoryMap_UMAP_",umap.type.name,"_CFtargets_normalized.pdf"), 
         w=length(unique(pDT.hex$tissue))*1+2,
         h=length(unique(pDT.hex$factor))*1+1,
         limitsize=FALSE)
}

