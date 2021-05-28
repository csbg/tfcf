source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("INT_03_SeuratIntegration_CellCycle/")

require(Seurat)

# Read cellranger analysis results --------------------------------------------
clusters.c1 <- fread(PATHS$CITESEQ1_CLEAN$DATA$clusters)
clusters.c1$dataset <- "CITESEQ1"
clusters.e1 <- fread(PATHS$ECCITE1$DATA$clusters)
clusters.e1$dataset <- "ECCITE1"
clusters <- rbind(clusters.c1, clusters.e1)

# Read guides merge with clusters
guides <- fread(PATHS$ECCITE1$DATA$guides)
names(guides) <- gsub("^(.)", "\\U\\1", names(guides), perl = T)
guides$dataset <- "ECCITE1"
clusters <- merge(clusters, guides, by=c("Barcode", "dataset"), all=TRUE)

# Read data  and Seurat integration --------------------------------------------
seurat.file <- out("SeuratObject.RData")
if(!file.exists(seurat.file)){
  data.c1 <- Read10X_h5(PATHS$CITESEQ1_CLEAN$DATA$matrix)
  data.e1 <- Read10X_h5(PATHS$ECCITE1$DATA$matrix)
  
  # 1
  #seurat.c1 <- CreateSeuratObject(counts = data.c1$"Gene Expression"[1:4000,1:1000], project = "CITESEQ1", min.cells = 5)
  seurat.c1 <- CreateSeuratObject(counts = data.c1$"Gene Expression", project = "CITESEQ1", min.cells = 5)
  seurat.c1 <- subset(seurat.c1, subset = nFeature_RNA > 500)
  seurat.c1 <- NormalizeData(seurat.c1, verbose = FALSE)
  seurat.c1 <- FindVariableFeatures(seurat.c1, selection.method = "vst", nfeatures = 2000)
  
  # 2
  #seurat.e1 <- CreateSeuratObject(counts = data.e1[1:4000,1:1000], project = "ECCITE1", min.cells = 5)
  seurat.e1 <- CreateSeuratObject(counts = data.e1, project = "ECCITE1", min.cells = 5)
  seurat.e1 <- subset(seurat.e1, subset = nFeature_RNA > 500)
  seurat.e1 <- NormalizeData(seurat.e1, verbose = FALSE)
  seurat.e1 <- FindVariableFeatures(seurat.e1, selection.method = "vst", nfeatures = 2000)
  
  # Integrate
  anchors <- FindIntegrationAnchors(object.list = list(seurat.c1, seurat.e1), dims = 1:20)
  sobj <- IntegrateData(anchorset = anchors, dims = 1:20)
  DefaultAssay(sobj) <- "integrated"
  sobj@meta.data$dataset <- sobj@meta.data$orig.ident
  
  # Cell Cycle scoring
  sobj <- CellCycleScoring(sobj, 
                           s.features = cc.genes$s.genes, 
                           g2m.features = cc.genes$g2m.genes, 
                           set.ident = TRUE)
  
  # Process full dataset
  sobj <- ScaleData(sobj, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
  sobj <- RunPCA(sobj, npcs = 30, verbose = FALSE)
  sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:20)
  sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:20)
  sobj <- FindClusters(sobj, resolution = 0.5)
  
  # Store full dataset
  sobj2 <- sobj
  sobj <- SCRNA.DietSeurat(sobj)
  save(sobj, file=seurat.file)
} else {
  print("Loading Seurat file")
  load(seurat.file)
  sobj <- SCRNA.UndietSeurat(sobj)
}


# Collect and export metadata --------------------------------------------------------------

ann <- data.table(sobj[["umap"]]@cell.embeddings, keep.rownames = TRUE)
ann <- merge(ann, data.table(sobj@meta.data, keep.rownames = TRUE), by="rn")
ann[, Barcode := gsub("_\\d+$", "", rn)]
ann <- merge(ann, clusters, by=c("Barcode", "dataset"))
colnames(ann) <- gsub("_", ".", colnames(ann))
ann[,SCluster := as.numeric(seurat.clusters) + 1]
ann$seurat.clusters <- NULL
ann$integrated.snn.res.0.5 <- NULL
write.tsv(ann, out("Metadata.tsv"))

stopifnot(nrow(ann[dataset == "CITESEQ1" & !is.na(Guide)]) == 0)



# SETUP ENDS HERE ---------------------------------------------------------



# PLOT DATASET AND CELL CYCLE -------------------------------------------------------------------
ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_point(aes(color=dataset), alpha=1, size=0.1) +
  theme_bw(12)
ggsave(out("UMAP_datasets.pdf"), w=6, h=5)

ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_point(aes(color=Phase), alpha=1, size=0.1) +
  theme_bw(12)
ggsave(out("UMAP_CellCycle.pdf"), w=6, h=5)


# PLOT ORIGINAL CLUSTERS ------------------------------------------------
ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_point(aes(color=factor(Cluster))) +
  facet_grid(. ~ dataset) +
  theme_bw(12) +
  geom_label(data=ann[,.(UMAP.1=median(UMAP.1), UMAP.2=median(UMAP.2)), by=c("Cluster", "dataset")], aes(label=Cluster))
ggsave(out("UMAP_Clusters_points.pdf"), w=11, h=5)

ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_hex() +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_grid(. ~ dataset) +
  theme_bw(12) +
  geom_label(data=ann[,.(UMAP.1=median(UMAP.1), UMAP.2=median(UMAP.2)), by=c("Cluster", "dataset")], aes(label=Cluster))
ggsave(out("UMAP_Clusters.pdf"), w=11, h=5)


# PLOT INTEGRATED CLUSTERS ------------------------------------------------
ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_point(aes(color=factor(SCluster))) +
  facet_grid(. ~ dataset) +
  theme_bw(12) +
  geom_label(data=ann[,.(UMAP.1=median(UMAP.1), UMAP.2=median(UMAP.2)), by=c("SCluster", "dataset")], aes(label=SCluster))
ggsave(out("UMAP_SClusters_points.pdf"), w=11, h=5)

ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_hex() +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_grid(. ~ dataset) +
  theme_bw(12) +
  geom_label(data=ann[,.(UMAP.1=median(UMAP.1), UMAP.2=median(UMAP.2)), by=c("SCluster", "dataset")], aes(label=SCluster))
ggsave(out("UMAP_SClusters.pdf"), w=11, h=5)


# PLOT GUIDES ------------------------------------------------
pDT <- ann[!is.na(Guide) & !grepl("None", Guide) & !grepl(" ", Guide)]
ggplot(pDT, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_hex() +
  scale_fill_gradient(low="lightgrey", high="blue") +
  #geom_point(data=, aes(color=Guide)) +
  facet_wrap(~ Guide) +
  theme_bw(12)
ggsave(out("UMAP_Guides.pdf"), w=10, h=6)

# GUIDE ENRICHMENT
res <- data.table()
pDT <- ann[!is.na(Guide)][!grepl(" ", Guide)]
for(gx in unique(pDT[Guide != "None"]$Guide)){
  for(cx in unique(pDT$SCluster)){
    fish <- fisher.test(as.matrix(with(ann[Guide %in% c(gx, "None")], table(SCluster == cx, Guide == gx))))
    res <- rbind(res, data.table(Cluster=cx, Guide=gx, p=fish$p.value, OR=fish$estimate))
  }
}
res[,padj := p.adjust(p, method="BH")]
res[padj < 0.05]
ggplot(res, aes(
  x=factor(Cluster),
  y=Guide, 
  color=log2(OR + min(res[OR != 0]$OR)), 
  size=pmin(-log10(padj), 5))
) + 
  geom_point(shape=16) +
  scale_color_gradient2(name="log2OR", low="blue", high="red") +
  scale_size_continuous(name="padj") + 
  theme_bw(12)
ggsave(out("Guides_Fisher.pdf"), w=5, h=4)


# ANTIBODIES --------------------------------------------------------------
abMT <- data.c1$`Antibody Capture`
abMT <- SCRNA.TPXToLog(SCRNA.RawToTPX(abMT, scale.factor = 1e6))
ann.c1 <- ann[dataset == "CITESEQ1"]

# Plot signal
res <- data.table()
abx <- row.names(abMT)[1]
for(abx in row.names(abMT)){
  pDT <- copy(ann.c1)
  pDT$Signal <- abMT[abx,ann.c1$Barcode]
  pDT$Antibody <- abx
  res <- rbind(res, pDT)
}
res[,Signal.norm := scale(Signal), by="Antibody"]
ggplot(res, aes(x=UMAP.1, y=UMAP.2)) + 
  stat_summary_hex(aes(z=Signal.norm),fun=mean) +
  scale_fill_gradient2(low="blue", high="red") +
  facet_wrap(~Antibody) +
  theme_bw(12)
ggsave(out("Antibodies_UMAP.pdf"), w=12+2, h=9+1)

cMT <- corS(t(as.matrix(abMT)))
diag(cMT) <- NA
dist <- as.dist(1-cMT)
cleanDev(); pdf(out("Antibodies_Correlation.pdf"), w=5,h=4)
pheatmap(cMT,
         clustering_distance_rows = dist,
         clustering_distance_cols = dist,
         breaks=seq(-1,1, 0.01), color=COLORS.HM.FUNC(200),
)
dev.off()

probx=0.9
res[,percentile := quantile(Signal.norm, probs=probx), by="Antibody"]
resN <- res[, sum(Signal.norm > percentile), by=c("SCluster", "Antibody")]
resN[,clSize := sum(V1), by="SCluster"]
stopifnot(all(resN[, length(unique(clSize)), by="SCluster"]$V1 == 1))
resN[,percentage := V1/clSize*100]
resN <- hierarch.ordering(resN, toOrder = "SCluster", orderBy = "Antibody", value.var = "percentage")
resN <- hierarch.ordering(resN, toOrder = "Antibody", orderBy = "SCluster", value.var = "percentage")
ggplot(resN, aes(x=SCluster, y=Antibody, fill=percentage)) + 
  geom_tile() + 
  ggtitle("Percent of cell in SCluster\nthat are above 90th percentile of antibody signal") + 
  scale_fill_gradient(low="white", high="blue")
ggsave(out("Antibodies_Percentile.pdf"), w=5, h=4)



# MIXSCAPE -----------------------------------------------------
sobj2 <- sobj
x <- unique(gsub(".+_(\\d+)$", "\\1", data.table(sobj2@meta.data, keep.rownames = T)[dataset == "ECCITE1"]$rn))
stopifnot(length(x) == 1)
g2 <- guides[!grepl(" ", Guide)][Guide != "None"]
sobj2@meta.data$guide <- g2[match(row.names(sobj2@meta.data), paste0(g2$Barcode, "_", x))]$Guide
with(data.table(sobj2@meta.data), table(guide, dataset))
with(g2, table(Guide))
stopifnot(all(data.table(sobj2@meta.data, keep.rownames = TRUE)[guide == "NTC"][order(rn)]$rn == paste0(g2[Guide == "NTC"][order(Barcode)]$Barcode, "_",x)))
sobj2 <- sobj2[,!is.na(sobj2@meta.data$guide)]
with(data.table(sobj2@meta.data), table(guide, dataset))
eccite <- CalcPerturbSig(
  object = sobj2, assay = "RNA", slot = "data", 
  gd.class ="guide",nt.cell.class = "NTC",
  reduction = "pca",ndims = 20,num.neighbors = 20,
  new.assay.name = "PRTB")
eccite <- ScaleData(object = eccite, assay = "PRTB", do.scale = F, do.center = T)
eccite <- RunMixscape(
  object = eccite,assay = "PRTB",slot = "scale.data",
  labels = "guide",nt.class.name = "NTC",
  min.de.genes = 5,
  iter.num = 10,
  de.assay = "RNA",
  verbose = TRUE,
  prtb.type = "KO")

table(eccite$mixscape_class.global, eccite$guide)
prop.table(table(eccite$mixscape_class.global, eccite$guide),2)

Idents(eccite) <- "mixscape_class"
res <- data.table()
for(guidex in unique(grep(" KO$", eccite@meta.data$mixscape_class, value = TRUE))){
  print(guidex)
  markers <- FindMarkers(
    eccite, ident.1 = guidex, ident.2 = "NTC", logfc.threshold = log(2),
    assay = "RNA",only.pos = FALSE, test.use="negbinom")
  res <- rbind(res, data.table(markers, guide=guidex, keep.rownames = TRUE))
}




# Create Dotplot ----------------------------------------------------------
gg <- unique(c(
  res[order(avg_log2FC, decreasing = TRUE)][,head(.SD,n=5), by="guide"]$rn,
  res[order(avg_log2FC, decreasing = FALSE)][,head(.SD,n=5), by="guide"]$rn))

eccite@meta.data$mix_cluster <- with(eccite@meta.data, paste(gsub(" ", "_", mixscape_class), seurat_clusters))
Idents(eccite) <- "mix_cluster"
avExp <- AverageExpression(eccite, assays="RNA")
perExp <- sapply(
  with(data.table(eccite@meta.data, keep.rownames = TRUE), split(rn, mix_cluster)),
  function(cells){
    rowSums(eccite@assays$RNA@counts[,cells, drop = FALSE] != 0) / length(cells)
  }) * 100

dotplotDT <- merge(
  melt(data.table(avExp$RNA, keep.rownames = TRUE), id.vars = "rn", value.name = "average"),
  melt(data.table(perExp, keep.rownames = TRUE), id.vars = "rn", value.name = "percentage"),
  by=c("rn", "variable"))

dotplotDT[, scaleExp := scale(average), by=c("rn")]

ggplot(dotplotDT[rn %in% gg], aes(x=as.character(variable), y=rn, size=percentage, color=scaleExp)) + 
  geom_point() + 
  scale_color_gradient2(low='blue', high = "red") +
  scale_size_continuous(limits=c(0, 100)) +
  theme_bw() +
  facet_grid(. ~ variable) + 
  xRot()

matrix <- as.matrix(colMeans(exp  > 0))*100

DotPlot(eccite, features = gg, assay="RNA") + xRot() + coord_flip()
