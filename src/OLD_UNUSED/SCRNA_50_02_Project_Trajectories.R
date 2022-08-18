source("src/00_init.R")

require(ProjecTILs)
require(umap)

base.dir <- "SCRNA_50_02_ProjectionTrajectories/"
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
#singleR.cell.types <- readRDS(dirout_load("SCRNA_06_02_MergeMarkers")("CellTypes_in.vivo.RDS"))
pseudotime <- fread(dirout_load("SCRNA_50_01_Trajectories")("Values.tsv"))


# Function to transform monocle3 to Seurat objects ----------------------------------
x <- mobjs$in.vivo
as.Seurat.NF <- function(x){
  logcounts(x) <- counts(x)
  x <- as.Seurat(x)
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e6)
  x <- RenameAssays(x, RNA = "integrated")
  x
}


# Prep
ref.monocle <- mobjs$in.vivo
ref.monocle$pseudotime <- pseudotime[match(colnames(ref.monocle), rn)]$traj
ref.monocle$pseudotime.celltype <- pseudotime[match(colnames(ref.monocle), rn)]$celltype
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
# stopifnot(all(colnames(ref) %in% singleR.cell.types$cellname))
# ref <- AddMetaData(ref, as.factor(singleR.cell.types[match(colnames(ref), cellname),]$labels), col.name = "functional.cluster")

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
    
    # Cross-map UMAP
    proj.trajectory <- ref.trajectory.predict(ref=ref.use, query=proj)
    
    # Cross-projected to original UMAP
    write.tsv(
      merge(
        proj.trajectory,
        data.table(proj@reductions$umap@cell.embeddings, keep.rownames = TRUE),
        by="rn"
      ), out("Output_", tx, "_",sx,".tsv"))
  }
}






# test statistcs ---------------------------------------------------------------
ann1 <- fread(dirout_load("SCRNA_20_Summary/ex.vivo_monocle.singleR")("Annotation.tsv"))
ann2 <- fread(dirout_load("SCRNA_20_Summary/leukemia_monocle.singleR")("Annotation.tsv"))
ann <- rbind(ann1, ann2, fill=TRUE)

ff <- list.files(out(""), pattern="Output.*.tsv")
ff <- ff[!grepl("in.vivo", ff)]
names(ff) <- gsub("Output_leukemia_(.+).tsv", "\\1", ff)
pDT <- rbindlist(lapply(ff, function(fx) fread(out(fx))), idcol = "sample")
pDT[, traj := pseudotime]
pDT$pseudotime <- NULL
write.tsv(pDT, out("Values.tsv"))
pDT <- merge(pDT, ann, by="rn")[!is.na(mixscape_class.global)]
pDT[, gene := gsub("_.+$", "", CRISPR_Cellranger)]
pDT[, celltype := ct]
pDT[, traj.scale := scale(traj), by="celltype"]
pDT[, id := paste0(tissue, "_", timepoint)]

(idx <- unique(pDT$id)[1])
for(idx in unique(pDT$id)){
  outS <- dirout(paste0(base.dir, idx))
  pDTx <- pDT[id == idx]

  # . test ------------------------------------------------------------------
  typex <- "Ery"
  gx <- "Rcor1"
  res <- data.table()
  for(typex in unique(pDT$celltype)){
    pDT1 <- pDTx[celltype == typex]
    for(gx in unique(pDT1[mixscape_class.global != "NTC"]$gene)){
      x1 <- pDT1[gene == gx]$traj
      x2 <- pDT1[gene == "NTC"]$traj
      if(length(x1) > 10 & length(x2) > 10){
        res <- rbind(res, data.table(
          p.wx=wilcox.test(x1, x2)$p.value,
          p.ks=ks.test(x1, x2)$p.value,
          d=median(x1) - median(x2),
          type=typex,
          gene=gx
        ))
      }
    }
  }
  res[, padj.wx := p.adjust(p.wx, method="BH")]
  res[, padj.ks := p.adjust(p.ks, method="BH")]
  write.tsv(res, outS("Statistics.tsv"))
  
  
  # . load ------------------------------------------------------------------
  #res <- fread(outS("Statistics.tsv"))
  
  # . plot stats ------------------------------------------------------------
  ggplot(res, aes(y=-log10(p.wx+1e-10), x=-log10(p.ks+1e-10), color=d)) + 
    geom_point() +
    geom_point(shape=1, color="black") +
    theme_bw(12)+ 
    geom_text_repel(aes(label=paste(gene)), color="black")+
    scale_color_gradient2(low="blue", high="red") +
    facet_grid(. ~ type)
  ggsave(outS("Statistics_Comparison.pdf"),w=20,h=6)
  
  # xDT <- melt(res, id.vars = c("type", "gene"))
  # xDT[, measurement := gsub("\\..+$", "", variable)]
  # xDT[, type := gsub("^.+?\\.", "", variable)]
  ggplot(res, aes(y=gene, x=type, size=pmin(5, -log10(padj.wx)), color=d)) + 
    geom_point() +
    theme_bw(12)+ 
    scale_color_gradient2(low="blue", high="red") +
    xRot()
  ggsave(outS("Statistics.pdf"), w=4,h=10)
  
  # . UMAP ------------------------------------------------------------------
  pDT.UMAP <- pDTx[abs(traj.scale) < 3]
  ggplot(pDT.UMAP, aes(x=UMAP1, y=UMAP2)) + 
    stat_summary_hex(bins = 100, aes(z=traj.scale),fun=mean) +
    theme_bw(12) +
    scale_fill_gradientn(colors=c("lightgrey", "blue", "purple", "red", "orange")) 
  ggsave(outS("UMAP.pdf"), w=5,h=5)
  

  # . celltypes -------------------------------------------------------------
  ggplot(pDTx, aes(x=UMAP1, y=UMAP2, color=celltype)) + 
    geom_point() +
    theme_bw(12)
  ggsave(outS("Celltypes_UMAP.jpg"), w=5,h=5)
  
  ggplot(pDTx, aes(x=UMAP1, y=UMAP2)) + 
    geom_hex() +
    facet_grid(. ~ celltype) + 
    theme_bw(12)
  ggsave(outS("Celltypes_UMAP_hex.pdf"), w=20,h=5)
  
  ggplot(pDTx, aes(x=celltype)) + 
    geom_bar() +
    scale_y_log10() +
    theme_bw(12)
  ggsave(outS("Celltypes_Nubmers.pdf"), w=5,h=5)
  
  # . plot distributions ----------------------------------------------------
  pDT.distr <- copy(pDTx)
  ggplot(pDT.distr, aes(x=traj, y=traj.scale)) + geom_hex() + facet_wrap(~celltype, scales = "free")
  pDT.distr <- pDT.distr[abs(traj.scale) < 3]
  pDT.sum <- pDT.distr[, .(traj = median(traj), q1 = quantile(traj, 0.25), q2 = quantile(traj, 0.75)),by=c("gene", "celltype")]
  pDT.stats <- copy(res)
  pDT.stats[, celltype := type]
  pDT.stats[, type := "not.sig"]
  pDT.stats[padj.ks < 0.1, type := "sig.low"]
  pDT.stats[padj.ks < 0.01, type := "sig.high"]
  pDT.distr <- merge(pDT.distr, pDT.stats[,c("gene", "celltype", "type"),with=F], by=c("gene", "celltype"), all.x=TRUE)
  pDT.distr[gene == "NTC", type := "NTC"]
  ggplot(pDT.distr, aes(y=gene, x=traj)) + 
    geom_violin(color=NA, aes(fill=type), scale="width") + 
    scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
    geom_point(data=pDT.sum, color="black") + 
    geom_errorbarh(data=pDT.sum, color="black", aes(xmin=q1, xmax=q2), height = .2) + 
    theme_bw(12) + 
    facet_grid(. ~ celltype, scale="free_x") +
    xRot()
  ggsave(outS("Distribution.pdf"), w=20,h=10)
  
  ggplot(pDT.distr[gene %in% c("Brd9", "Smarcd2", "Smarcd1", 'NTC')][celltype %in% "Mye"],
         aes(x=traj, color=gene)) + 
    theme_bw(12) +
    geom_density()
  ggsave(outS("Distribution_Brd9.pdf"), w=5,h=4)
  
  
  # . scatterplot -----------------------------------------------------------
  pDT2 <- pDTx[celltype %in% c("Ery", "Mye")][, median(traj), by=c("celltype", "gene")]
  pDT2[, V1 := scale(V1), by="celltype"]
  pDT2 <- dcast.data.table(pDT2, gene ~ celltype, value.var = "V1")
  ggplot(pDT2, aes(x=Ery, y=Mye)) + 
    geom_point() + 
    theme_bw(12) +
    geom_abline() + 
    geom_text_repel(aes(label=gene))
  ggsave(outS("Scatter_EryVsMye.pdf"),w=8,h=8)
}
