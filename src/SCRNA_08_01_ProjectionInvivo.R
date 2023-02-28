source("src/00_init.R")

require(ProjecTILs)
require(umap)

base.dir <- "SCRNA_08_01_ProjectionInvivo/"
out <- dirout(base.dir)

source("src/FUNC_ProjecTILs_PLUS.R")


# Annotation --------------------------------------------------------------
SANN <- fread(PATHS$SCRNA$ANN)

# Read in vivo data and perform differnetial expression -------------------
mobjs <- list()
for(tissuex in PATHS$SCRNA$MONOCLE.NAMES){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}


# singleR cell types ------------------------------------------------------
singleR.cell.types <- readRDS(dirout_load("SCRNA_06_02_MergeMarkers")("CellTypes_in.vivo.RDS"))


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

# the following code and code at the end of this section enables loading of previously stored data, this is dangerous if data gets updated
# ref.file <- out("reference.rds")
# if(file.exists(ref.file)){
#   ref <- readRDS(ref.file)
# } else {

ref.monocle <- mobjs$in.vivo
ref <- as.Seurat.NF(ref.monocle)
ref.umap.original <- reducedDims(ref.monocle)$UMAP
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
ref@misc$projecTILs="in.vivo"

# Add labels
stopifnot(all(colnames(ref) %in% singleR.cell.types$cellname))
ref <- AddMetaData(ref, as.factor(singleR.cell.types[match(colnames(ref), cellname),]$labels), col.name = "functional.cluster")

# Export table
write.tsv(
  merge(
    data.table(ref@meta.data, keep.rownames = TRUE),
    data.table(ref@reductions$UMAP@cell.embeddings, keep.rownames = TRUE),
    by="rn"
  ), out("Output_in.vivo",".tsv"))
  
#   # Save
#   saveRDS(ref, ref.file)
# }


# Run projection and cell type prediction ----------------------------------------------------------
ref.use <- ref
ref.use@reductions$umap <- ref.use@reductions$UMAP
ref.use@reductions$pca <- ref.use@reductions$PCA

tx <- "leukemia"
for(tx in setdiff(PATHS$SCRNA$MONOCLE.NAMES, "in.vivo")){
  
  query <- as.Seurat.NF(mobjs[[tx]])
  
  sx <- query$sample[14]
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
    
    # Cross-map UMAP
    proj.umap.original <- ref.umap.predict(ref=ref.use, query=proj, ref.umap = ref.umap.original)
    
    
    # Export results
    # Normal
    write.tsv(
      merge(
        data.table(pred@meta.data, keep.rownames = TRUE),
        data.table(proj@reductions$umap@cell.embeddings, keep.rownames = TRUE),
        by="rn"
      ), out("Output_", tx, "_",sx,".tsv"))
    
    # Cross-projected to original UMAP
    write.tsv(
      merge(
        data.table(pred@meta.data, keep.rownames = TRUE),
        data.table(proj.umap.original, keep.rownames = TRUE),
        by="rn"
      ), out("OutputCrossprojection_", tx, "_",sx,".tsv"))
    
  }
}




# Summarize results -------------------------------------------------------
ff <- list.files(out(""), pattern="Output_")
ff <- lapply(ff, function(fx) fread(out(fx)))

ff2 <- list.files(out(""), pattern="OutputCrossprojection_")
ff2 <- lapply(ff2, function(fx) fread(out(fx)))

pDT.list <- list(
  original=rbindlist(ff, fill=TRUE),
  crossproject=rbindlist(ff2, fill=TRUE)
  )

for(xx in names(pDT.list)){
  pDT <- pDT.list[[xx]]
  
  # Hex plot
  ggplot(pDT, aes(x=UMAP_1, y=UMAP_2)) + 
    theme_bw(12) +
    geom_hex(data=pDT[tissue == "in.vivo"], bins=100) + 
    geom_density_2d(data=pDT[tissue != "in.vivo"], aes(color=tissue))
  ggsave(out(xx, "_Hexplot.pdf"), w=6,h=5)
  
  # Hexplot by tissue
  ggplot(pDT, aes(x=UMAP_1, y=UMAP_2)) + 
    theme_bw(12) +
    stat_binhex(aes(fill=log10(..count..)), bins=100) + 
    facet_grid(. ~ tissue)
  ggsave(out(xx, "_Hexplot.byTissue.pdf"), w=16,h=5)
  
  # Celltypes
  ggplot(pDT, aes(x=UMAP_1, y=UMAP_2, color=functional.cluster)) + 
    theme_bw(12) +
    scale_color_brewer(palette = "Paired") +
    geom_point(shape=1, alpha=0.3)
  ggsave(out(xx, "_Celltypes.png"), w=7,h=5)
  
  # Numbers of cell types
  ggplot(pDT, aes(x=functional.cluster)) + 
    theme_bw(12) +
    geom_bar() + 
    facet_grid(. ~ tissue) +
    xRot() +
    scale_y_log10()
  ggsave(out(xx, "_Celltypes.numbers.pdf"), w=9,h=4)
  
  
  # Compare to SingleR ------------------------------------------------------
  xDT <- merge(pDT, singleR.cell.types, by.x="rn", by.y="cellname")
  
  if(nrow(xDT) == 0) next
  
  jDT <- data.table()
  for(tx in unique(xDT$tissue)){
    jMT <- jaccard.twolists(
      l1=with(xDT[tissue == tx], split(rn, functional.cluster)),
      l2=with(xDT[tissue == tx], split(rn, labels))
    )
    
    jDT <- rbind(jDT, data.table(melt(data.table(jMT, keep.rownames = TRUE), id.vars = "rn"), tissue=tx))
  }
  ggplot(jDT, aes(x=rn, y=variable, fill=value)) +
    theme_bw(12) +
    geom_tile() +
    facet_grid(. ~ tissue) +
    xRot() +
    scale_fill_gradient(low="white", high="blue")
  ggsave(out(xx, "_ComparisonToSingleR.pdf"), w=13,h=5)
}

