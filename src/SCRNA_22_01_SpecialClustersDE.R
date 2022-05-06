source("src/00_init.R")
out <- dirout("SCRNA_22_01_SpecialClustersDE/")

tissuex <- "in.vivo"

# Load object
(load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
mobj <- monocle.obj
fData(mobj)$gene_short_name <- row.names(fData(mobj))
ann <- data.table(data.frame(colData(mobj)), keep.rownames = TRUE)


# annotate cell types
ann.cts <- readRDS(dirout_load("SCRNA_06_02_MergeMarkers")("CellTypes_", tissuex,".RDS"))
mobj$celltype <- ann.cts[match(colnames(mobj), cellname),]$labels

# Clusters
cl.test <- fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("ClustersRemoved_in.vivo.tsv"))$Cluster.number

clx <- cl.test[1]
for(clx in cl.test){
  message(clx)
  cellx <- names(sort(table(mobj$celltype[clusters(mobj) %in% clx]), decreasing = TRUE))[1]
  print(cellx)
  mobj.test <- mobj[,mobj$celltype == cellx]
  mobj.test$de.group <- factor(ifelse(clusters(mobj.test) == clx, "x", "ref"), levels=c("ref", "x"))
  res <- fit_models(mobj.test, model_formula_str = "~de.group")
  res <- coefficient_table(res)
  res <- data.table(res)[status == "OK"][term == "de.groupx"][,-c("model", "model_summary", "term", "std_err", "test_val", "normalized_effect", "model_component", "status")]
  write.tsv(res, out("Markers_",clx,".tsv"))
}


ff <- list.files(out(""), full.names = TRUE, pattern = "Markers_\\d+.tsv$")
names(ff) <- basename(ff)
pDT <- rbindlist(lapply(ff, fread), idcol = "cluster", fill = TRUE)

cl <- pDT$cluster[1]
for(cl in unique(pDT$cluster)){
  gg <- pDT[cluster == cl][q_value < 0.05][order(-estimate)][1:100]$gene_short_name
  if(length(gg) > 2)
  p <- plot_genes_by_group(mobj, markers = gg, group_cells_by = "cluster")
  ggsave(out(cl,".pdf"), w=15,h=20)
}