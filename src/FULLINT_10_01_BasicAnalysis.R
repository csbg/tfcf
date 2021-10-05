source("src/00_init.R")
out <- dirout("FULLINT_10_01_BasicAnalysis/")


(load(PATHS$FULLINT$Monocle))


marker.genes <- fread("metadata/markers.csv")





# Calculate clusters ------------------------------------------------------
monocle.obj = cluster_cells(monocle.obj, resolution=1e-5)


# Collect annotation --------------------------------------------------------------
ann <- data.table(data.frame(colData(monocle.obj)@listData), keep.rownames = TRUE)
ann$Clusters <- as.character(monocle.obj@clusters$UMAP$clusters[ann$rn])
umap <- setNames(data.table(reducedDims(monocle.obj)$UMAP, keep.rownames = TRUE), c("rn", "UMAP1", "UMAP2"))
ann <- merge(ann, umap, by="rn", all=TRUE)



# ADDITIONAL QC to remove bad clusters --------------------------------------------------------
for(qcm in c("percent.mt", "nFeature_RNA", "nCount_RNA")){
  print(qcm)
  ann$measure <- ann[[qcm]]
  p <- ggplot(ann, aes(y=measure + 0.1, x=factor(Clusters))) + 
    geom_violin(color=NA, fill="lightblue") + 
    geom_boxplot(fill=NA, coef=Inf) +
    scale_y_log10() +
    theme_bw(12) +
    ylab(qcm)
  ggtitle(qcm)
  ggsave(out("QC_", qcm, ".pdf"), w=5,h=4, plot=p)
}
ann[,cluster.qual.keep :=TRUE]
ann[Clusters %in% ann[,median(percent.mt), by="Clusters"][V1 < 1]$Clusters, cluster.qual.keep := FALSE]



# CELLRANGER -------------------------------------------
if("AGG.CSV" %in% ls()){
  ann.exp <- fread(out("Metadata_3_Mixscape.tsv"))
  ann.exp <- merge(ann.exp, AGG.CSV[,c("sample_id", "i"),with=F], by.x="dataset", by.y="sample_id", all.x=TRUE)
  stopifnot(!any(is.na(ann.exp$i)))
  ann.exp[,Barcode := paste0(gsub("-.+$", "", Barcode), "-", i)]
  write.table(ann.exp[cluster.qual.keep == TRUE][,c("Barcode", "UMAP.1", "UMAP.2"),with=F], file=out("Cellranger_UMAP.csv"), sep=",", col.names = c("Barcode", "UMAP-1", "UMAP-2"), quote=F, row.names = F)
  write.table(ann.exp[cluster.qual.keep == TRUE][,c("Barcode", "mixscape_class"),with=F], file=out("Cellranger_MIXSCAPE.csv"), sep=",", col.names = c("Barcode", "MIXSCAPE"), quote=F, row.names = F)
  write.table(ann.exp[cluster.qual.keep == TRUE][,c("Barcode", "cluster"),with=F], file=out("Cellranger_Clusters.csv"), sep=",", col.names = c("Barcode", "Clusters_Seurat"), quote=F, row.names = F)
}

# Analysis of clusters ----------------------------------------------------




# ANTIBODIES --------------------------------------------------------------
abMT <- additional.info$CITESEQ2
abMT <- SCRNA.TPXToLog(SCRNA.RawToTPX(abMT, scale.factor = 1e6))
ann.c1 <- ann[dataset == "CITESEQ2"]

# Plot signal
res <- data.table()
abx <- row.names(abMT)[1]
for(abx in row.names(abMT)){
  pDT <- copy(ann.c1)
  pDT$Signal <- abMT[abx,ann.c1$Barcode]
  pDT$Antibody <- abx
  res <- rbind(res, pDT)
}
res[,Signal.norm := scale(Signal), by="Antibody"]
ggplot(res, aes(x=UMAP.1, y=UMAP.2)) + 
  stat_summary_hex(aes(z=Signal.norm),fun=mean) +
  scale_fill_gradient2(low="blue", high="red") +
  facet_wrap(~Antibody) +
  theme_bw(12)
ggsave(out("Antibodies_UMAP.pdf"), w=12+2, h=9+1)

cMT <- corS(t(as.matrix(abMT)), use="pairwise.complete.obs")
diag(cMT) <- NA
dist <- as.dist(1-cMT)
cleanDev(); pdf(out("Antibodies_Correlation.pdf"), w=5,h=4)
pheatmap(cMT,
         clustering_distance_rows = dist,
         clustering_distance_cols = dist,
         breaks=seq(-1,1, 0.01), color=COLORS.HM.FUNC(200),
)
dev.off()

probx=0.9
res[,percentile := quantile(Signal.norm, probs=probx, na.rm=TRUE), by="Antibody"]
resN <- res[, sum(Signal.norm > percentile), by=c("cluster", "Antibody")]
resN[,clSize := sum(V1), by="cluster"]
stopifnot(all(resN[, length(unique(clSize)), by="cluster"]$V1 == 1))
resN[,percentage := V1/clSize*100]
resN <- hierarch.ordering(resN, toOrder = "cluster", orderBy = "Antibody", value.var = "percentage")
resN <- hierarch.ordering(resN, toOrder = "Antibody", orderBy = "cluster", value.var = "percentage")
resN[is.na(percentage), percentage := 0]
ggplot(resN, aes(x=cluster, y=Antibody, fill=percentage)) + 
  geom_tile() + 
  ggtitle("Percent of cell in cluster\nthat are above 90th percentile of antibody signal") + 
  scale_fill_gradient(low="white", high="blue")
ggsave(out("Antibodies_Percentile.pdf"), w=5, h=4)



# GUIDES ------------------------------------------------------

# Number of guides --------------------------------------------------------
# Number of cells per guide / mixscape
ggplot(ann[!is.na(mixscape_class.global)], aes(x=guide, fill=mixscape_class.global)) + 
  geom_bar(position="dodge") + 
  facet_wrap(~sample, scales = "free") +
  theme_bw(12) + xRot() +
  ylab("Cells")
ggsave(out("Guides_Counts.pdf"), w=12,h=10)
# % of cells assigned
ggplot(ann[,sum(!is.na(guide))/.N*100, by="sample"], aes(x=sample, y=V1)) + 
  geom_bar(stat="identity") + 
  theme_bw(12) + xRot() +
  ylab("Percent of cells with assigned guide (NP or KO)")
ggsave(out("Guides_CountsNA.pdf"), w=6,h=4)

# Guides in UMAP ----------------------------------------------------------
sx <- "ECCITE1"
for(sx in unique(ann[!is.na(mixscape_class)]$sample)){
  pDT <- ann[sample == sx][!is.na(mixscape_class.global)]
  grps <- length(unique(pDT$mixscape_class))
  ggplot(pDT, aes(x=UMAP1, y=UMAP2)) + 
    geom_hex() +
    scale_fill_gradient(low="lightgrey", high="blue") +
    facet_wrap(~mixscape_class + mixscape_class.global, ncol = 7) +
    theme_bw(12)
  ggsave(out("Guides_UMAP_", sx, ".pdf"), w=7*4,h=ceiling(grps/7)*4+1)
}

# Guides per cluster ------------------------------------------------------
res <- data.table()
pDT1 <- ann[!is.na(mixscape_class)][mixscape_class.global != "NP"]
for(sx in unique(pDT1$sample)){
  pDT2 <- pDT1[sample == sx]
  for(gx in unique(pDT2[mixscape_class != "NTC"]$mixscape_class)){
    for(cx in unique(pDT2$Clusters)){
      message(gx, "-", cx)
      mx <- as.matrix(with(ann[sample == sx][mixscape_class %in% c(gx, "NTC")], table(Clusters == cx, mixscape_class == gx)))
      print(mx)
      if(dim(mx) == c(2,2)){
        fish <- fisher.test(mx)
        res <- rbind(res, data.table(Clusters=cx, mixscape_class=gx, p=fish$p.value, OR=fish$estimate, sample=sx))
      }
    }
  }
}
res[,padj := p.adjust(p, method="BH")]
res[, log2OR := log2(pmin(5, OR + min(res[OR != 0]$OR)))]
res[,grp := paste(mixscape_class, sample)]
res <- hierarch.ordering(res, toOrder = "grp", orderBy = "Clusters", value.var = "log2OR")
res <- hierarch.ordering(res, toOrder = "Clusters", orderBy = "grp", value.var = "log2OR")
ggplot(res, aes(
  x=Clusters,
  y=mixscape_class, 
  color=log2OR, 
  size=pmin(-log10(padj), 5))) + 
  geom_point(shape=16) +
  scale_color_gradient2(name="log2OR", low="blue", high="red") +
  scale_size_continuous(name="padj") + 
  facet_grid(sample ~ ., space = "free", scales = "free") +
  theme_bw(12) +
  theme(strip.text.y = element_text(angle=0))
ggsave(out("Guides_Fisher.pdf"), w=10, h=10)