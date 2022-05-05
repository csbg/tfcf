source("src/00_init.R")
out <- dirout("SCRNA_22_01_SpecialClustersDE/")

tissuex <- "in.vivo"

# Load object
(load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
mobj <- monocle.obj
fData(mobj)$gene_short_name <- row.names(fData(mobj))
ann <- data.table(data.frame(colData(mobj)), keep.rownames = TRUE)


# annotate cell types
cts <- readRDS(dirout_load("SCRNA_06_02_MergeMarkers")("CellTypes_", tissuex,".RDS"))
mobj$celltype <- cts[match(colnames(mobj), cellname),]$labels

# Clusters
cl.test <- fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("ClustersRemoved_in.vivo.tsv"))$Cluster.number

clx <- cl.test[1]
for(clx in cl.test){
  ctx <- names(sort(table(ann[clusters(mobj) %in% cl.test]$celltype), decreasing = TRUE))[1]
  mobj.test <- mobj[,mobj$celltype == ctx]
  mobj.test$de.group <- factor(ifelse(clusters(mobj.test) == clx, "x", "ref"), levels=c("ref", "x"))
  res <- fit_models(mobj.test, model_formula_str = "~de.group")
  res <- coefficient_table(res)
  res <- data.table(res)[status == "OK"][term == "de.groupx"][,-c("model", "model_summary", "term", "std_err", "test_val", "normalized_effect", "model_component", "status")]
  res <- res[q_value < 0.1]
  write.tsv(markers, out("Markers_",clx,".tsv"))
}
