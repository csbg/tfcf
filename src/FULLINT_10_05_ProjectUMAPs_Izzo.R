source("src/00_init.R")

require(ProjecTILs)
require(umap)

base.dir <- "FULLINT_10_05_ProjectUMAPs_Izzo/"
out <- dirout(base.dir)


# Folders -----------------------------------------------------------------
inDir.funcs <- list(
  "in.vivo"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vivo"),
  "in.vitro"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vitro"),
  "leukemia"=dirout_load("FULLINT_10_01_BasicAnalysis_leukemia")
)


# Read in datasets -------------------
mobjs <- list()
for(tissuex in c("in.vitro", "leukemia", "in.vivo")){
  (load(inDir.funcs[[tissuex]]("MonocleObject.RData")))
  mobjs[[tissuex]] <- monocle.obj
}


# Function to transform monocle3 to Seurat objects ----------------------------------
x <- mobjs$in.vivo
as.Seurat.NF <- function(x){
  logcounts(x) <- counts(x)
  x <- as.Seurat(x)
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e6)
  x <- RenameAssays(x, RNA = "integrated")
  x
}


# Load izzo dataset -------------------------------------------------------
(load(dirout_load("FULLINT_05_01_SingleR")("izzo.RData")))
izzoCDS <- new_cell_data_set(expression_data = izzoMT, cell_metadata = data.frame(row.names=colnames(izzoMT), tissue=rep("Izzo", ncol(izzoMT))))

# Prepare reference -------------------------------------------------------
ref.file <- out("reference.rds")
if(file.exists(ref.file)){
  ref <- readRDS(ref.file)
} else {
  ref <- as.Seurat.NF(izzoCDS)
  ref <- FindVariableFeatures(ref)
  
  # PCA
  set.seed(1234)
  which.assay="integrated"
  varfeat <- ref@assays[[which.assay]]@var.features
  refdata <- data.frame(t(ref@assays[[which.assay]]@data[varfeat,]))
  refdata <- refdata[, sort(colnames(refdata))]
  ref.pca <- prcomp(refdata, rank. = 50, scale. = TRUE, center = TRUE, retx=TRUE)
  
  # UMAP
  seed=1234
  n.neighbors=30
  min.dist=0.3
  metric="cosine"
  ndim=10
  umap.config <- umap.defaults
  umap.config$n_neighbors = n.neighbors
  umap.config$min_dist = min.dist
  umap.config$metric = metric
  umap.config$n_components = 2
  umap.config$random_state = seed
  umap.config$transform_state = seed
  ref.umap <- umap(ref.pca$x[,1:ndim], config=umap.config)
  colnames(ref.umap$layout) <- c("UMAP_1","UMAP_2")
  
  # add to object
  
  # create a list with required components
  ref <- ScaleData(ref)
  ref <- RunPCA(ref)
  ref <- RunUMAP(ref, reduction = "pca", dims = 1:5)
  ref@reductions$PCA <- ref@reductions$pca
  ref@reductions$UMAP <- ref@reductions$umap
  ref@reductions$UMAP@cell.embeddings <- ref.umap$layout
  ref@reductions$PCA@cell.embeddings <- ref.pca$x
  ref@reductions$PCA@feature.loadings <- ref.pca$rotation
  colnames(ref@reductions$PCA@cell.embeddings) <- gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$x), perl=TRUE)
  colnames(ref@reductions$PCA@feature.loadings) <- gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$rotation), perl=TRUE)
  #Store the complete PCA and UMAP object in @misc
  ref@misc$pca_object <- ref.pca
  ref@misc$umap_object <- ref.umap
  ref@misc$projecTILs="Izzo"
  
  # Add labels
  ref <- AddMetaData(ref, as.factor(gsub("\\-?\\d$", "", izzo.ann[match(colnames(ref), V1),]$clusterName)), col.name = "functional.cluster")
  
  # Export table
  write.tsv(
    merge(
      data.table(ref@meta.data, keep.rownames = TRUE),
      data.table(ref@reductions$UMAP@cell.embeddings, keep.rownames = TRUE),
      by="rn"
    ), out("Output_izzo",".tsv"))
  
  # Save
  saveRDS(ref, ref.file)
}


# Run projection and cell type prediction ----------------------------------------------------------
ref.use <- ref
ref.use@reductions$umap <- ref.use@reductions$UMAP
ref.use@reductions$pca <- ref.use@reductions$PCA

tx <- "leukemia"
for(tx in names(mobjs)){
  
  query <- as.Seurat.NF(mobjs[[tx]])
  
  sx <- query$sample[1]
  for(sx in unique(query$sample)){
    
    # Only cells from one sample
    query.use <- subset(query, cells=colnames(query)[query$sample == sx])
    query.use@reductions$umap <- query.use@reductions$UMAP
    query.use@reductions$pca <- query.use@reductions$PCA
    
    # Make projection
    proj <- make.projection(
      query = query.use,
      ref = ref.use,
      filter.cells = FALSE,
      fast.mode = TRUE,
      seurat.k.filter=100
    )
    
    # Prediction
    pred <- cellstate.predict(ref=ref.use, query=proj)
    
    # Export results
    write.tsv(
      merge(
        data.table(pred@meta.data, keep.rownames = TRUE),
        data.table(proj@reductions$umap@cell.embeddings, keep.rownames = TRUE),
        by="rn"
      ), out("Output_", tx, "_",sx,".tsv"))
  }
}




# Summarize results -------------------------------------------------------
ff <- list.files(out(""), pattern="Output_")
ff <- lapply(ff, function(fx) fread(out(fx)))
pDT <- rbindlist(ff, fill=TRUE)

# Hex plot
ggplot(pDT, aes(x=UMAP_1, y=UMAP_2)) + 
  theme_bw(12) +
  geom_hex(data=pDT[tissue == "Izzo"], bins=100) + 
  geom_density_2d(data=pDT[tissue != "Izzo"], aes(color=tissue))
ggsave(out("Hexplot.pdf"), w=6,h=5)

# Hexplot by tissue
ggplot(pDT, aes(x=UMAP_1, y=UMAP_2)) + 
  theme_bw(12) +
  stat_binhex(aes(fill=log10(..count..)), bins=100) + 
  facet_grid(. ~ tissue)
ggsave(out("Hexplot.byTissue.pdf"), w=21,h=5)

# Celltypes
ggplot(pDT, aes(x=UMAP_1, y=UMAP_2, color=functional.cluster)) + 
  theme_bw(12) +
  scale_color_brewer(palette = "Paired") +
  geom_point(shape=1, alpha=0.3)
ggsave(out("Celltypes.png"), w=7,h=5)

# Numbers of cell types
ggplot(pDT, aes(x=functional.cluster)) + 
  theme_bw(12) +
  geom_bar() + 
  facet_grid(. ~ tissue) +
  xRot() +
  scale_y_log10()
ggsave(out("Celltypes.numbers.pdf"), w=12,h=4)


# Antibodies
inDir.current <- "leukemia"
abs <- fread(inDir.funcs[[inDir.current]]("Antibodies.tsv"))
absDT <- merge(pDT, abs[,c("Antibody", "Signal.norm", "rn"), with=F], by="rn")
ggplot(absDT, aes(x=UMAP_1, y=UMAP_2)) +
  stat_summary_hex(bins = 100, aes(z=pmin(abs(Signal.norm), 2) * sign(Signal.norm)),fun=mean) +
  scale_fill_gradient2(low="blue", high="red") +
  facet_wrap(~Antibody) +
  theme_bw(12)
ggsave(out("Antibodies_UMAP.pdf"), w=12+2, h=9+1)