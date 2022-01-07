source("src/00_init.R")

base.dir <- "FULLINT_10_03_NewUMAPs/"
out <- dirout(base.dir)


# Folders -----------------------------------------------------------------
list.files(dirout_load("")(""))
inDir.funcs <- list(
  "in.vivo"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vivo"),
  "in.vitro"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vitro"),
  "leukemia"=dirout_load("FULLINT_10_01_BasicAnalysis_leukemia")
)


# Read in vivo data and perform differnetial expression -------------------
(load(inDir.funcs[["in.vivo"]]("MonocleObject.RData")))
monocle.obj.invivo <- monocle.obj
in.vivo.markers <- top_markers(
  monocle.obj,
  group_cells_by = "cluster",
  genes_to_test_per_group = 200,
  reduction_method = "UMAP",
  marker_sig_test = TRUE,
  reference_cells = NULL,
  speedglm.maxiter = 25,
  cores = 10,
  verbose = FALSE
)
in.vivo.markers <- data.table(in.vivo.markers)

# Plot
gg <- unique(in.vivo.markers[fraction_expressing >= 0.10][order(pseudo_R2)][, head(.SD, n=3), by=c("cell_group")]$gene_id)
fData(monocle.obj)$gene_short_name <- row.names(fData(monocle.obj))
plot_genes_by_group(monocle.obj,markers = gg,
                    ordering_type="cluster_row_col",
                    max.size=3)
ggsave(out("Invivo_Markers.pdf"), w=12, h=12)

# Store
write.tsv(in.vivo.markers, out("Invivo_Markers.tsv"))


# Create UMAP in vitro and leukemia based on in vivo DE genes -------------
gg <- unique(in.vivo.markers[marker_test_q_value < 0.05]$gene_id)
tissuex <- "in.vitro"
for(tissuex in c("in.vitro", "leukemia")){
  (load(inDir.funcs[[tissuex]]("MonocleObject.RData")))
  
  monocle.obj <- monocle.obj[gg,]
  
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
  
  monocle.obj = cluster_cells(monocle.obj)
  plot_cells(monocle.obj)
  ggsave(out("NewUMAP_", tissuex, ".jpg"), w=6,h=6)
  
  umap <- setNames(data.table(reducedDims(monocle.obj)$UMAP, keep.rownames = TRUE), c("rn", "UMAP1", "UMAP2"))
  write.tsv(umap, out("NewUMAP_", tissuex, ".tsv"))
}

