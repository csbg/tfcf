source("src/00_init.R")
out <- dirout("SCRNA_07_01_MarkerGenes/")

tissuex <- "in.vivo"
(load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
mobj <- monocle.obj
fData(mobj)$gene_short_name <- row.names(fData(mobj))
cts <- readRDS(dirout_load("SCRNA_06_02_MergeMarkers")("CellTypes_", tissuex,".RDS"))
mobj$celltype <- cts[match(colnames(mobj), cellname),]$labels
table(mobj$celltype)

markers <- top_markers(mobj, group_cells_by = "celltype", cores=ceiling(length(unique(mobj$celltype))/2))
write.tsv(markers, out("Markers.tsv"))

str(markers.split <- with(markers, split(gene_id, cell_group)))

# Marker gene expression --------------------------------------------------
for(ct in names(markers.split)){
  plot_genes_by_group(mobj, markers = markers.split[[ct]], group_cells_by = "celltype")
  ggsave(out("Markers_",make.names(ct),".pdf"), w=10,h=8)
}
