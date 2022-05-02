source("src/00_init.R")
out <- dirout("SCRNA_50_01_Trajectories/")


# Load annotation ---------------------------------------------------------
(load(PATHS$SCRNA$MONOCLE.DIR("in.vivo")))
# ann <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjMonocle_celltypes.RDS"))
# ann <- ann[rn %in% colnames(monocle.obj)]
ann <- fread(dirout_load("SCRNA_20_Summary/in.vivo_monocle.singleR")("Annotation.tsv"))
celltypes <- fread('metadata/FIGS_celltypes.tsv')
celltypes <- celltypes[Name %in% ann$Clusters]

#cds <- monocle.obj[,sample(colnames(monocle.obj), 500)]

typex <- "Mye"
for(typex in unique(celltypes[Type != "HSC"]$Type)){
  print(typex)
  cells <- ann[Clusters %in% celltypes[Type %in% c(typex, "HSC")]$Name]$rn
  #cells <- sample(cells, 100)
  cds <- monocle.obj[,cells]
  
  # recluster (not used but necessary)
  cds <- monocle3::cluster_cells(cds)
  
  # learn graph
  cds <- learn_graph(
    cds,
    verbose = TRUE,
    use_partition = FALSE,
    close_loop = TRUE)
  
  # order cells
  cds <- order_cells(cds, root_cells=ann[rn %in% colnames(cds)][functional.cluster == "HSC"]$rn)
  
  # Make plot
  plot_cells(cds,color_cells_by = "pseudotime")
  ggsave(out("Pseudotime_",typex,".jpg"), w=10,h=10)
  
  # export table
  write.tsv(data.table(data.frame(pseudotime(cds)), keep.rownames = TRUE), out("Pseudotime_",typex,".tsv"))
}



# test statistcs ---------------------------------------------------------------
ff <- list.files(out(""), pattern="Pseudotime.*.tsv")
names(ff) <- gsub("Pseudotime_(.+).tsv", "\\1", ff)
pDT <- rbindlist(lapply(ff, function(fx) fread(out(fx))), idcol = "celltype")
pDT[, traj := pseudotime.cds.]
pDT$pseudotime.cds. <- NULL
write.tsv(pDT, out("Values.tsv"))
pDT <- merge(pDT, ann[Clusters != "HSC"], by="rn")[!is.na(mixscape_class.global)]
pDT[, gene := gsub("_.+$", "", CRISPR_Cellranger)]


typex <- "Ery"
gx <- "Rcor1"
res <- data.table()
for(typex in unique(pDT$celltype)){
  pDT1 <- pDT[celltype == typex]
  for(gx in unique(pDT1[mixscape_class.global != "NTC"]$gene)){
    x1 <- pDT1[gene == gx]$traj
    x2 <- pDT1[gene == "NTC"]$traj
    res <- rbind(res, data.table(
      p=wilcox.test(x1, x2)$p.value,
      d=median(x1) - median(x2),
      type=typex,
      gene=gx
    ))
  }
}
res[, padj := p.adjust(p, method="BH")]
write.tsv(res, out("Wilcoxon.tsv"))

ggplot(res, aes(y=gene, x=type, size=pmin(5, -log10(padj)), color=d)) + 
  geom_point() +
  theme_bw(12)+ 
  scale_color_gradient2(low="blue", high="red") +
  xRot()
ggsave(out("Statistics.pdf"), w=4,h=10)

ggplot(pDT, aes(y=gene, x=traj)) + 
  geom_violin(color=NA, fill="blue") + 
  geom_boxplot(coef=1e10, fill=NA, color="black") +
  theme_bw(12) + 
  facet_grid(. ~ celltype, scale="free_x") +
  xRot()
ggsave(out("Distribution.pdf"), w=20,h=10)

