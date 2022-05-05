source("src/00_init.R")
out <- dirout("SCRNA_07_02_MarkerGenes_Numbers/")

tissuex <- "in.vivo"
(load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
mobj <- monocle.obj
fData(mobj)$gene_short_name <- row.names(fData(mobj))

markers <- top_markers(mobj, group_cells_by = "cluster", cores=10)
write.tsv(markers, out("Markers.tsv"))
#markers <- fread(out("Markers.tsv"))

str(markers.split <- with(markers, split(gene_id, cell_group)))

# Marker gene expression --------------------------------------------------
for(ct in names(markers.split)){
  plot_genes_by_group(mobj, markers = markers.split[[ct]], group_cells_by = "cluster")
  ggsave(out("Markers_",make.names(ct),".pdf"), w=10,h=8)
}
