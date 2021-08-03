source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("INT_03_SeuratIntegration_CellCycle/")

require(Seurat)

# Read cellranger analysis results --------------------------------------------
# clusters.c1 <- fread(PATHS$CITESEQ1_CLEAN$DATA$clusters)
# clusters.c1$dataset <- "CITESEQ1"
# clusters.e1 <- fread(PATHS$ECCITE1$DATA$clusters)
# clusters.e1$dataset <- "ECCITE1"
# clusters <- rbind(clusters.c1, clusters.e1)

# Read guides merge with clusters
guides <- fread(PATHS$ECCITE1$DATA$guides)
names(guides) <- gsub("^(.)", "\\U\\1", names(guides), perl = T)
guides$dataset <- "ECCITE1"
#clusters <- merge(clusters, guides, by=c("Barcode", "dataset"), all=TRUE)

# Read data  and Seurat integration --------------------------------------------
seurat.file <- out("SeuratObject.RData")
if(!file.exists(seurat.file)){
  
  seurat.list <- list()
  additional.info <- list()
  
  for(dsx in c("CITESEQ1", "CITESEQ2", "ECCITE1", "ECCITE2")){
    print("---------------")
    print(dsx)
    print("---------------")
    if(!file.exists(PATHS[[dsx]]$DATA$matrix)){
      message("File for ", dsx, " not found")
      next
    }
    data <- Read10X_h5(PATHS[[dsx]]$DATA$matrix)
    data.gx <- NA
    if(class(data) == "dgCMatrix"){
      data.gx <- data
    }
    if(is.list(data) && "Gene Expression" %in% names(data)){
      data.gx <- data[["Gene Expression"]]
      additional.info[[dsx]] <- data[names(data) != "Gene Expression"]
    }
    if(class(data.gx) != "dgCMatrix"){
      message("Data for ", dsx, " not in the right format")
    }
    seurat.obj <- CreateSeuratObject(counts = data.gx, project = dsx, min.cells = 5)
    seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")
    seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
    seurat.obj <- NormalizeData(seurat.obj, verbose = FALSE)
    seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
    seurat.list[[dsx]] <- seurat.obj
  }
  
  # Integrate
  anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:20)
  rm(list=c("seurat.list", "data.gx", "data"))
  sobj <- IntegrateData(anchorset = anchors, dims = 1:20)
  rm(list="anchors")
  
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
  
  
  # add ECCITE1 guides
  
  read.csv(PATHS$ECCITE1$DATA$guides)
  matrix_dir = paste(Sys.getenv("DATA"), "ECCITE1_citeseq_combined//umi_count/", sep="/")
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  colnames(mat) = paste0(barcode.names$V1, "-1")
  rownames(mat) = gsub("\\-.+$", "", feature.names$V1)
  mat <- mat[row.names(mat) != "unmapped",]
  additional.info[["ECCITE1"]] <- list("CRISPR Guide Capture" = as(mat, "dgCMatrix"))
  
  # Store full dataset
  sobj <- SCRNA.DietSeurat(sobj)
  save(sobj, additional.info, file=seurat.file)
} else {
  print("Loading Seurat file")
  load(seurat.file)
}
sobj <- SCRNA.UndietSeurat(sobj)


# Collect metadata --------------------------------------------------------------
ann <- data.table(sobj[["umap"]]@cell.embeddings, keep.rownames = TRUE)
ann <- merge(ann, data.table(sobj@meta.data, keep.rownames = TRUE), by="rn")
ann[, Barcode := gsub("_\\d+$", "", rn)]
#ann <- merge(ann, clusters, by=c("Barcode", "dataset"))
colnames(ann) <- gsub("_", ".", colnames(ann))
ann[,cluster := as.numeric(as.character(seurat.clusters))]
ann$seurat.clusters <- NULL
ann$integrated.snn.res.0.5 <- NULL
ann$orig.ident <- NULL
ann$old.ident <- NULL
ann.orig <- copy(ann)


# Add extra meta data (antibodies, guides) --------------------------------
cutoff <- 5
dsx <- names(additional.info)[1]
for(dsx in names(additional.info)){
  print(dsx)
  tx <- names(additional.info[[dsx]])[1]
  for(tx in names(additional.info[[dsx]])){
    txn <- make.names(tx)
    print(tx)
    x <- additional.info[[dsx]][[tx]]
    x <- as.matrix(x[,apply(x, 2, max) >= cutoff]) # only keep cells with any guide counted above cutoff
    ii <- apply(x, 2, function(col) which(col >= cutoff)) # for each cell (column) get the rows that are above the cutoff
    labelsx <- sapply(ii, function(x) paste(sort(names(x)), collapse=","))
    labelsx <- setNames(data.table(data.frame(labelsx), keep.rownames = TRUE), c("Barcode", txn))
    prep.ann <- if(txn %in% colnames(ann)) ann[, -txn,with=F] else ann
    ann <- rbind(ann[dataset != dsx], merge(prep.ann[dataset == dsx], labelsx, by="Barcode", all.x=TRUE), fill=TRUE)
  }
}
ann[, Guide := CRISPR.Guide.Capture]
ann[, Guide := gsub("_.+$", "", Guide)]
ann[grepl(",", Guide), Guide := NA]

write.tsv(ann, out("Metadata.tsv"))

for(dsx in unique(ann$dataset)){
  if(grepl("CITESEQ", dsx)){
    stopifnot(nrow(ann[dataset == dsx & !is.na(CRISPR.Guide.Capture)]) == 0)
    stopifnot(nrow(ann[dataset == dsx & !is.na(Antibody.Capture)]) != 0)
  } else {
    stopifnot(nrow(ann[dataset == dsx & !is.na(CRISPR.Guide.Capture)]) != 0)
    stopifnot(nrow(ann[dataset == dsx & !is.na(Antibody.Capture)]) == 0)
  }
}




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

# 
# # PLOT ORIGINAL CLUSTERS ------------------------------------------------
# ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
#   geom_point(aes(color=factor(Cluster))) +
#   facet_grid(. ~ dataset) +
#   theme_bw(12) +
#   geom_label(data=ann[,.(UMAP.1=median(UMAP.1), UMAP.2=median(UMAP.2)), by=c("Cluster", "dataset")], aes(label=Cluster))
# ggsave(out("UMAP_Clusters_points.pdf"), w=11, h=5)
# 
# ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
#   geom_hex() +
#   scale_fill_gradient(low="lightgrey", high="blue") +
#   facet_grid(. ~ dataset) +
#   theme_bw(12) +
#   geom_label(data=ann[,.(UMAP.1=median(UMAP.1), UMAP.2=median(UMAP.2)), by=c("Cluster", "dataset")], aes(label=Cluster))
# ggsave(out("UMAP_Clusters.pdf"), w=11, h=5)


# PLOT INTEGRATED CLUSTERS ------------------------------------------------
ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_point(aes(color=factor(cluster))) +
  facet_wrap( ~ dataset, ncol = 2) +
  theme_bw(12) +
  geom_label(data=ann[,.(UMAP.1=median(UMAP.1), UMAP.2=median(UMAP.2)), by=c("cluster", "dataset")], aes(label=cluster))
ggsave(out("UMAP_clusters_points.pdf"), w=11, h=10)

ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_hex() +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap( ~ dataset, ncol = 2) +
  theme_bw(12) +
  geom_label(data=ann[,.(UMAP.1=median(UMAP.1), UMAP.2=median(UMAP.2)), by=c("cluster", "dataset")], aes(label=cluster))
ggsave(out("UMAP_cluster.pdf"), w=11, h=10)


# PLOT GUIDES ------------------------------------------------
ggplot(ann[!is.na(Guide)], aes(x=UMAP.1, y=UMAP.2)) + 
  geom_hex() +
  scale_fill_gradient(low="lightgrey", high="blue") +
  #geom_point(data=, aes(color=Guide)) +
  facet_wrap(~ gsub("_.+", "", Guide)) +
  theme_bw(12)
ggsave(out("UMAP_Guides.pdf"), w=7, h=5)

# GUIDE ENRICHMENT
res <- data.table()
pDT <- ann[!is.na(Guide)][!grepl(" ", Guide)]
for(gx in unique(pDT[Guide != "NTC"]$Guide)){
  for(cx in unique(pDT$cluster)){
    mx <- as.matrix(with(ann[Guide %in% c(gx, "NTC")], table(cluster == cx, Guide == gx)))
    if(dim(mx) == c(2,2)){
      fish <- fisher.test(mx)
      res <- rbind(res, data.table(Cluster=cx, Guide=gx, p=fish$p.value, OR=fish$estimate))
    }
  }
}
res[,padj := p.adjust(p, method="BH")]
res[padj < 0.05]
res[, log2OR := log2(pmin(5, OR + min(res[OR != 0]$OR)))]
res <- hierarch.ordering(res, toOrder = "Guide", orderBy = "Cluster", value.var = "log2OR")
res <- hierarch.ordering(res, toOrder = "Cluster", orderBy = "Guide", value.var = "log2OR")
ggplot(res, aes(
  x=Cluster,
  y=Guide, 
  color=log2OR, 
  size=pmin(-log10(padj), 5))) + 
  geom_point(shape=16) +
  scale_color_gradient2(name="log2OR", low="blue", high="red") +
  scale_size_continuous(name="padj") + 
  theme_bw(12)
ggsave(out("Guides_Fisher.pdf"), w=5, h=4)


# ANTIBODIES --------------------------------------------------------------
abMT <- additional.info$CITESEQ2$`Antibody Capture`
abMT <- SCRNA.TPXToLog(SCRNA.RawToTPX(abMT, scale.factor = 1e6))
ann.c1 <- ann[dataset == "CITESEQ2"]

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

probx=0.9
res[,percentile := quantile(Signal.norm, probs=probx, na.rm=TRUE), by="Antibody"]
resN <- res[, sum(Signal.norm > percentile), by=c("cluster", "Antibody")]
resN[,clSize := sum(V1), by="cluster"]
stopifnot(all(resN[, length(unique(clSize)), by="cluster"]$V1 == 1))
resN[,percentage := V1/clSize*100]
resN <- hierarch.ordering(resN, toOrder = "cluster", orderBy = "Antibody", value.var = "percentage")
resN <- hierarch.ordering(resN, toOrder = "Antibody", orderBy = "cluster", value.var = "percentage")
resN[is.na(percentage), percentage := 0]
ggplot(resN, aes(x=cluster, y=Antibody, fill=percentage)) + 
  geom_tile() + 
  ggtitle("Percent of cell in cluster\nthat are above 90th percentile of antibody signal") + 
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

markers.all <- res


# Create Dotplot ----------------------------------------------------------
gg <- unique(c(
  res[order(avg_log2FC, decreasing = TRUE)][,head(.SD,n=5), by="guide"]$rn,
  res[order(avg_log2FC, decreasing = FALSE)][,head(.SD,n=5), by="guide"]$rn))

eccite@meta.data$mix_cluster <- with(eccite@meta.data, paste(gsub(" ", "_", mixscape_class), seurat_clusters))
Idents(eccite) <- "mixscape_class"
avExp <- AverageExpression(eccite, assays="RNA")
perExp <- sapply(
  with(data.table(eccite@meta.data, keep.rownames = TRUE), split(rn, mixscape_class)),
  function(cells){
    rowSums(eccite@assays$RNA@counts[,cells, drop = FALSE] != 0) / length(cells)
  }) * 100

dotplotDT <- merge(
  melt(data.table(avExp$RNA, keep.rownames = TRUE), id.vars = "rn", value.name = "average",variable.factor = FALSE),
  melt(data.table(perExp, keep.rownames = TRUE), id.vars = "rn", value.name = "percentage",variable.factor = FALSE),
  by=c("rn", "variable"))

dotplotDT <- merge(dotplotDT, dotplotDT[variable == "NTC"][,c("rn", "average")], by="rn", suffixes = c("", "_NTC"))
dotplotDT[, scaleExp := average - average_NTC]
dotplotDT[, scaleExp := scale(scaleExp, center = FALSE), by="rn"]

#dotplotDT <- cbind(dotplotDT, setNames(data.table(do.call(rbind, strsplit(dotplotDT$variable, " "))), c("guide", "cluster")))
markers.all.plot <- markers.all
markers.all.plot[, percentage := pct.1 * 100]
markers.all.plot[, variable := guide]
ggplot(dotplotDT[rn %in% gg], aes(x=variable, y=rn, size=percentage, color=scaleExp)) + 
  geom_point() + 
  geom_point(data=markers.all.plot[rn %in% gg & p_val_adj < 0.05], shape=1, color="black") + 
  scale_color_gradient2(low='blue', high = "red") +
  scale_size_continuous(limits=c(0, 100)) +
  theme_bw() +
  xRot()
ggsave(out("Mixscape_Guides.pdf"), w=8,h=8)



ggplot(markers.all[rn %in% gg], aes(x=guide, y=rn, size=percentage, color=scaleExp)) + 
  geom_point() + 
  scale_color_gradient2(low='blue', high = "red") +
  scale_size_continuous(limits=c(0, 100)) +
  theme_bw() +
  xRot()

DoHeatmap(size = 3,
  eccite, features = gg, group.by = "mix_cluster", assay = "RNA", slot="data",
  cells=data.table(eccite@meta.data, keep.rownames = TRUE)[guide %in% c("NTC", "Pu.1")]$rn)
ggsave(out("Mixscape_HM.jpeg"), w=20,h=6, dpi = 100)

# ggplot(dotplotDT[rn %in% gg], aes(x=factor(as.numeric(cluster)), y=rn, size=percentage, color=scaleExp)) + 
#   geom_point() + 
#   scale_color_gradient2(low='blue', high = "red") +
#   scale_size_continuous(limits=c(0, 100)) +
#   theme_bw() +
#   facet_grid(. ~ guide) + 
#   xRot()
# 
# 
# matrix <- as.matrix(colMeans(exp  > 0))*100

DotPlot(eccite, features = gg, assay="RNA") + xRot() + coord_flip()
ggsave(out("Mixscape_Guides_Seurat.pdf"), w=8,h=8)




# Second approach, by cluster ---------------------------------------------
Idents(eccite) <- "mixscape_class"
cons.markers <- FindConservedMarkers(min.cells.group=0,
  object = eccite, ident.1="Pu.1 KO", ident.2 = "NTC", 
  assay = "RNA", grouping.var = "seurat_clusters")

x <- data.table(cons.markers, keep.rownames = TRUE)
res <- data.table()
i <- 0
for(i in unique(eccite@meta.data$seurat_clusters)){
  ret <- data.table(x[,c("rn", grep(paste0("^", i), colnames(x), value=TRUE)),with=F], cluster=i)
  colnames(ret) <- gsub("^\\d+_", "", colnames(ret))
  res <- rbind(res, ret, fill=TRUE)
}

ggplot(res, aes(x=factor(as.numeric(cluster)), y=rn, color=avg_log2FC, size=pmin(5, -log10(p_val_adj)))) + 
  geom_point() + scale_color_gradient2(low="blue", high="red") +
  theme_bw(12)
ggsave(out("Mixscape_byCluster.pdf"), w=8,h=8)
