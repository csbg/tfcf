source("src/00_init.R")
out <- dirout("SCRNA_22_01_SpecialClustersDE/")

tissuex <- "in.vivo"

require(fgsea)
(load(PATHS$RESOURCES$Enrichr.mouse))

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



# Make plots --------------------------------------------------------------
ff <- list.files(out(""), full.names = TRUE, pattern = "Markers_\\d+.tsv$")
names(ff) <- basename(ff)
pDT <- rbindlist(lapply(ff, fread), idcol = "cluster", fill = TRUE)
pDT[q_value > 0.75, estimate := 0]

cl <- pDT$cluster[1]
for(cl in unique(pDT$cluster)){
  gg <- pDT[cluster == cl][q_value < 0.05][order(-estimate)][1:100]$gene_short_name
  if(length(gg) > 2)
  p <- plot_genes_by_group(mobj, markers = gg, group_cells_by = "cluster")
  ggsave(out(cl,".pdf"), w=15,h=20)
}



# Vulcano plot ------------------------------------------------------------
ggplot(pDT, aes(x=estimate, y=-log10(p_value))) + geom_hex() + facet_grid(. ~ cluster)
ggsave(out("Vulcano.pdf"), w = 20,h=4)

ggplot(pDT, aes(x=p_value)) + geom_histogram() + facet_grid(. ~ cluster)
ggsave(out("Pvalue_Distribution.pdf"), w = 20,h=4)


# Enrichments -------------------------------------------------------------
gsea.res <- data.table() 
de.grp <- pDT$cluster[1]
for(tx in unique(pDT$cluster)){
    for(dbx in names(enr.terms)){
      gsea.res <- rbind(gsea.res, data.table(fgsea(
        pathways=enr.terms[[dbx]], 
        stats=with(pDT[cluster == tx], setNames(estimate, nm=gene_id))), 
        cluster=tx,
        db=dbx))
  }
}
saveRDS(gsea.res, file=out("FGSEA.RDS"))


# Enrichment Plots -------------------------------------------------------------------
if(!"gsea.res" %in% ls()) gsea.res <- readRDS(out("FGSEA.RDS"))

# cleanup / export results
gsea.res[is.nan(NES), NES := 0]
gsea.res.export <- gsea.res[padj < 0.05][,-c("log2err", "NES", "size", "pval"),with=F]
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge, function(vec) paste(vec[1:10], collapse = ","))
write.tsv(gsea.res.export, out("GSEA_significant",".tsv"))

# Prepare for plotting
for(dbx in unique(gsea.res$db)){
  pDT <- gsea.res[db == dbx]
  pw.display <- unique(pDT[padj < 0.1][order(NES)][, tail(.SD, n=5), by=c("cluster")]$pathway)
  pDT <- pDT[pathway %in% pw.display]
  if(nrow(pDT) < 1) next
  try({
    pDT <- hierarch.ordering(pDT, "cluster", "pathway", "NES", TRUE)
    pDT <- hierarch.ordering(pDT, "pathway", "cluster", "NES", TRUE)
  }, silent=TRUE)
  ggplot(pDT, aes(x=cluster, y=pathway, color=NES, size=pmin(5, -log10(padj)))) +
    geom_point() + scale_color_gradient2(low="blue", mid="white", high="red") +
    geom_point(data=pDT[padj < 0.05], shape=1, color="black") +
    scale_size_continuous(range=c(0,5), limits = c(0,5)) +
    theme_bw(12) +
    xRot() +
    theme(strip.text.y=element_text(angle=0))
  ggsave(out("GSEA_plot_",dbx,".pdf"), w=8,h=length(unique(pDT$pathway)) * 0.2 + 3, limitsize = FALSE)
}

