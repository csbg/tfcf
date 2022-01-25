source("src/00_init.R")

require(ProjecTILs)
require(umap)

base.dir <- "FULLINT_10_04_ProjectUMAPs/"
out <- dirout(base.dir)


# Folders -----------------------------------------------------------------
inDir.funcs <- list(
  "in.vivo"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vivo"),
  "in.vitro"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vitro"),
  "leukemia"=dirout_load("FULLINT_10_01_BasicAnalysis_leukemia")
)


# Read in vivo data and perform differnetial expression -------------------
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


# Prepare reference -------------------------------------------------------
ref.file <- out("reference.rds")
if(file.exists(ref.file)){
  ref <- readRDS(ref.file)
} else {
  ref <- as.Seurat.NF(mobjs$in.vivo)
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
  ref@reductions$UMAP@cell.embeddings <- ref.umap$layout
  ref@reductions$PCA@cell.embeddings <- ref.pca$x
  ref@reductions$PCA@feature.loadings <- ref.pca$rotation
  colnames(ref@reductions$PCA@cell.embeddings) <- gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$x), perl=TRUE)
  colnames(ref@reductions$PCA@feature.loadings) <- gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$rotation), perl=TRUE)
  #Store the complete PCA and UMAP object in @misc
  ref@misc$pca_object <- ref.pca
  ref@misc$umap_object <- ref.umap
  ref@misc$projecTILs="in vivo"
  
  # Add labels
  cell.types <- fread(dirout_load("FULLINT_05_01_SingleR")("cell_types_marrow10x_label.main.csv"))
  ref <- AddMetaData(ref, as.factor(cell.types[match(colnames(ref), cell),]$labels), col.name = "functional.cluster")
  
  # Save
  saveRDS(ref, ref.file)
}


# Run projection and cell type prediction ----------------------------------------------------------
ref.use <- ref
ref.use@reductions$umap <- ref.use@reductions$UMAP
ref.use@reductions$pca <- ref.use@reductions$PCA

tx <- "leukemia"
for(tx in c("in.vitro", "leukemia")){
  
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


