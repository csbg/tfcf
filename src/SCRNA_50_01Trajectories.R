source("src/00_init.R")
out <- dirout("SCRNA_50_01_Trajectories/")

require(ggrepel)

ann <- fread(dirout_load("SCRNA_20_Summary/in.vivo_monocle.singleR")("Annotation.tsv"))

# Load annotation ---------------------------------------------------------
(load(PATHS$SCRNA$MONOCLE.DIR("in.vivo")))
# ann <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjMonocle_celltypes.RDS"))
# ann <- ann[rn %in% colnames(monocle.obj)]
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
pDT[, traj.scale := scale(traj), by="celltype"]
pDT <- pDT[timepoint != "28d"]

# . test ------------------------------------------------------------------
typex <- "Ery"
gx <- "Rcor1"
res <- data.table()
for(typex in unique(pDT$celltype)){
  pDT1 <- pDT[celltype == typex]
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
write.tsv(res, out("Statistics.tsv"))


# . load ------------------------------------------------------------------
res <- fread(out("Statistics.tsv"))


# . plot stats ------------------------------------------------------------
ggplot(res, aes(y=-log10(p.wx+1e-10), x=-log10(p.ks+1e-10), color=d)) + 
  geom_point() +
  geom_point(shape=1, color="black") +
  theme_bw(12)+ 
  geom_text_repel(aes(label=paste(gene)), color="black")+
  scale_color_gradient2(low="blue", high="red") +
  facet_grid(. ~ type)
ggsave(out("Statistics_Comparison.pdf"),w=20,h=6)

# xDT <- melt(res, id.vars = c("type", "gene"))
# xDT[, measurement := gsub("\\..+$", "", variable)]
# xDT[, type := gsub("^.+?\\.", "", variable)]
ggplot(res, aes(y=gene, x=type, size=pmin(5, -log10(padj.wx)), color=d)) + 
  geom_point() +
  theme_bw(12)+ 
  scale_color_gradient2(low="blue", high="red") +
  xRot()
ggsave(out("Statistics.pdf"), w=4,h=10)



# . UMAP ------------------------------------------------------------------
pDT.UMAP <- pDT[abs(traj.scale) < 3]
ggplot(pDT.UMAP, aes(x=UMAP1, y=UMAP2)) + 
  stat_summary_hex(bins = 100, aes(z=traj.scale),fun=mean) +
  theme_bw(12) +
  scale_fill_gradientn(colors=c("lightgrey", "blue", "purple", "red", "orange")) 
ggsave(out("UMAP.pdf"), w=5,h=5)


# . plot distributions ----------------------------------------------------
pDT.distr <- copy(pDT)
ggplot(pDT.distr, aes(x=traj, y=traj.scale)) + geom_hex() + facet_wrap(~celltype, scales = "free")
pDT.distr <- pDT.distr[abs(traj.scale) < 3]
pDT.sum <- pDT.distr[, .(traj = median(traj), q1 = quantile(traj, 0.25), q2 = quantile(traj, 0.75)),by=c("gene", "celltype")]
pDT.stats <- copy(res)
pDT.stats[, celltype := type]
pDT.stats[, type := "not.sig"]
pDT.stats[padj.ks < 0.1, type := "sig.low"]
pDT.stats[padj.ks < 0.01, type := "sig.high"]
pDT.distr <- merge(pDT.distr, pDT.stats[,c("gene", "celltype", "type"),with=F], by=c("gene", "celltype"), all.x=TRUE)
pDT.distr[is.na(type), type := "NTC"]
ggplot(pDT.distr, aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.sum, color="black") + 
  geom_errorbarh(data=pDT.sum, color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ celltype, scale="free_x") +
  xRot()
ggsave(out("Distribution.pdf"), w=20,h=10)

ggplot(pDT.distr[gene %in% c("Brd9", "Smarcd2", "Smarcd1", 'NTC')][celltype %in% "Mye"],
       aes(x=traj, color=gene)) + 
  theme_bw(12) +
  geom_density()
ggsave(out("Distribution_Brd9.pdf"), w=5,h=4)


# . scatterplot -----------------------------------------------------------
pDT2 <- pDT[celltype %in% c("Ery", "Mye")][, median(traj), by=c("celltype", "gene")]
pDT2[, V1 := scale(V1), by="celltype"]
pDT2 <- dcast.data.table(pDT2, gene ~ celltype, value.var = "V1")
ggplot(pDT2, aes(x=Ery, y=Mye)) + 
  geom_point() + 
  theme_bw(12) +
  geom_abline() + 
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_EryVsMye.pdf"),w=8,h=8)
