source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("INT_02_SeuratIntegration/")

require(Seurat)

# Integrating the two datasets --------------------------------------------

# Read data
data.c1 <- Read10X_h5(PATHS$CITESEQ1_CLEAN$DATA$matrix)
data.e1 <- Read10X_h5(PATHS$ECCITE1$DATA$matrix)

# Read clusters
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

# Create Seurat objects
# 1
seurat.c1 <- CreateSeuratObject(counts = data.c1$"Gene Expression", project = "CITESEQ1", min.cells = 5)
seurat.c1 <- subset(seurat.c1, subset = nFeature_RNA > 500)
seurat.c1 <- NormalizeData(seurat.c1, verbose = FALSE)
seurat.c1 <- FindVariableFeatures(seurat.c1, selection.method = "vst", nfeatures = 2000)
# 2
seurat.e1 <- CreateSeuratObject(counts = data.e1, project = "ECCITE1", min.cells = 5)
seurat.e1 <- subset(seurat.e1, subset = nFeature_RNA > 500)
seurat.e1 <- NormalizeData(seurat.e1, verbose = FALSE)
seurat.e1 <- FindVariableFeatures(seurat.e1, selection.method = "vst", nfeatures = 2000)

# Integrate
anchors <- FindIntegrationAnchors(object.list = list(seurat.c1, seurat.e1), dims = 1:20)
sobj <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(sobj) <- "integrated"
sobj@meta.data$dataset <- sobj@meta.data$orig.ident

# Process full dataset
sobj <- ScaleData(sobj, verbose = FALSE)
sobj <- RunPCA(sobj, npcs = 30, verbose = FALSE)
sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:20)
sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:20)
sobj <- FindClusters(sobj, resolution = 0.5)

# Store full dataset
#sobj <- DietSeurat(sobj)
save(sobj, file=out("SeuratObject.RData"))

# Collect information
ann <- data.table(sobj[["umap"]]@cell.embeddings, keep.rownames = TRUE)
ann <- merge(ann, data.table(sobj@meta.data, keep.rownames = TRUE), by="rn")
ann[, Barcode := gsub("_\\d+$", "", rn)]
ann <- merge(ann, clusters, by=c("Barcode", "dataset"))
colnames(ann) <- gsub("_", ".", colnames(ann))
write.tsv(ann, out("Metadata.tsv"))

# Plots -------------------------------------------------------------------
ggplot(ann, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_point(aes(color=dataset), alpha=1, size=0.1) +
  theme_bw(12)
ggsave(out("UMAP_datasets.pdf"), w=6, h=5)

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

pDT <- ann[!is.na(Guide) & !grepl("None", Guide) & !grepl(" ", Guide)]
ggplot(pDT, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_hex() +
  scale_fill_gradient(low="lightgrey", high="blue") +
  #geom_point(data=, aes(color=Guide)) +
  facet_wrap(~ Guide) +
  theme_bw(12)
ggsave(out("UMAP_Guides.pdf"), w=10, h=6)

ann[dataset == "CITESEQ1" & !is.na(Guide)]



# TESTING IF OUR FUNCTION WORKS
# c1.dat.scrna <- SCRNA.read_10Xh5.610(PATHS$CITESEQ1_CLEAN$DATA$matrix) 
# c1.dat.scrna.ex <- c1.dat.scrna$matrix[c1.dat.scrna$features[feature_type == "Gene Expression"]$id,]
# gg.unique <- c1.dat.scrna$features[feature_type == "Gene Expression"][,.N, by="name"][order(N)][N == 1]
# gg.unique <- c1.dat.scrna$features[name %in% gg.unique$name]
# c1.dat.scrna.ex.test <- c1.dat.scrna.ex[gg.unique$id,]
# row.names(c1.dat.scrna.ex.test) <- gg.unique$name
# # Now see if the values are the same
# bar <- colnames(c1.dat.seurat$"Gene Expression")
# gg <- gg.unique$name[1:1000]
# stopifnot(sum(as.matrix(c1.dat.seurat$"Gene Expression"[gg,bar])) > 0)
# stopifnot(all(as.matrix(c1.dat.seurat$"Gene Expression"[gg,bar]) == as.matrix(c1.dat.scrna.ex.test[gg,bar])))

