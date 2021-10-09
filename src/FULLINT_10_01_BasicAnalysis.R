source("src/00_init.R")
out <- dirout("FULLINT_10_01_BasicAnalysis/")

require(umap)
require(igraph)


# LOAD DATA ---------------------------------------------------------------
# single cell
(load(PATHS$FULLINT$Monocle))
fData(monocle.obj)$gene_short_name <- row.names(fData(monocle.obj))

# Markers
marker.genes <- fread("metadata/markers.csv")

# SingleR
ff <- list.files(dirout_load("FULLINT_05_01_SingleR")(""), pattern = "cell_types_.*.csv", full.names = TRUE)
singleR.res <- setNames(lapply(ff, fread), gsub("cell_types_(.+).csv", "\\1", basename(ff)))

# CytoTRACE


# CLUSTERING ------------------------------------------------------
set.seed(12121)
monocle.obj = cluster_cells(monocle.obj, resolution=1e-5)


# Collect ANNOTATION --------------------------------------------------------------
ann <- data.table(data.frame(colData(monocle.obj)@listData), keep.rownames = TRUE)
ann$Clusters <- as.character(monocle.obj@clusters$UMAP$clusters[ann$rn])
umap <- setNames(data.table(reducedDims(monocle.obj)$UMAP, keep.rownames = TRUE), c("rn", "UMAP1", "UMAP2"))
ann <- merge(ann, umap, by="rn", all=TRUE)



# SETUP ENDS HERE ---------------------------------------------------------




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
  ann.exp <- merge(ann, AGG.CSV[,c("sample_id", "i"),with=F], by.x="sample", by.y="sample_id", all.x=TRUE)
  stopifnot(!any(is.na(ann.exp$i)))
  ann.exp[,Barcode := paste0(gsub("-.+$", "", rn), "-", i)]
  write.table(ann.exp[cluster.qual.keep == TRUE][,c("Barcode", "UMAP1", "UMAP2"),with=F], file=out("Cellranger_UMAP.csv"), sep=",", col.names = c("Barcode", "UMAP-1", "UMAP-2"), quote=F, row.names = F)
  write.table(ann.exp[cluster.qual.keep == TRUE][,c("Barcode", "mixscape_class"),with=F], file=out("Cellranger_MIXSCAPE.csv"), sep=",", col.names = c("Barcode", "MIXSCAPE"), quote=F, row.names = F)
  write.table(ann.exp[cluster.qual.keep == TRUE][,c("Barcode", "Clusters"),with=F], file=out("Cellranger_Clusters.csv"), sep=",", col.names = c("Barcode", "Clusters_Seurat"), quote=F, row.names = F)
}



# SAMPLES -----------------------------------------------------------------
ggplot(ann, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw(12) +
  geom_hex() +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~sample, ncol=5)
ggsave(out("Samples_UMAP.pdf"), w=5*2+2,h=ceiling(length(unique(ann$sample))/5) * 2 + 1)

pDT <- ann[,.N, by=c("sample", "Clusters")]
pDT[,sumS := sum(N), by="sample"]
pDT[,sumC := sum(N), by="Clusters"]
pDT[,percentS := N/sumS*100]
pDT[,percentC := N/sumC*100]
ggplot(pDT, aes(y=sample,x=factor(as.numeric(Clusters)), size=percentS, color=percentC)) +
  scale_size_continuous(name="% of sample") + 
  scale_color_gradient(name="% of cluster", low="black", high="red") + 
  theme_bw() + 
  geom_point()
ggsave(out("Samples_Clusters.pdf"), w=8,h=length(unique(pDT$sample)) * 0.3+1)


# CLUSTERS ----------------------------------------------------
# UMAP
ggplot(ann, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw(12) +
  geom_hex() +
  scale_fill_gradient(low="lightgrey", high="blue") +
  geom_label(data=ann[,.(UMAP1=median(UMAP1), UMAP2=median(UMAP2)), by="Clusters"], aes(label=Clusters), fill="#ffffffaa")
ggsave(out("Clusters_UMAP.pdf"), w=6,h=5)

# Markers
plot_genes_by_group(monocle.obj, markers = marker.genes$Name, group_cells_by = "cluster") + scale_size_continuous(range=c(0,5))
ggsave(out("Clusters_Markers.pdf"), w=6,h=8)




# CELLTYPES SingleR -------------------------------------------------------
srx <- names(singleR.res)[1]
for(srx in names(singleR.res)){
  print(srx)
  #   x <- singleR.res[[srx]]
  #   x[3:5]
  #   x <- merge(
  #     melt(x[,c("cell", grep("^score", colnames(x), value = TRUE)), with=F], id.vars = "cell"),
  #     ann[,c("UMAP1", "UMAP2", "rn"),with=F],
  #     by.x="cell", by.y="rn")
  #   print(quantile(x[,median(value), by="variable"]$V1))
  # }
  
  singleR.resX <- singleR.res[[srx]]
  #x[3:5]
  pDT <- merge(
    melt(singleR.resX[,c("cell", grep("^score", colnames(singleR.resX), value = TRUE)), with=F], id.vars = "cell"),
    ann[,c("UMAP1", "UMAP2", "rn", "Clusters"),with=F],
    by.x="cell", by.y="rn")
  
  # Violin plots
  # ggplot(pDT, aes(x=variable, y=value)) + 
  #   theme_bw(12) + 
  #   geom_violin() + 
  #   xRot()
  # ggsave(out("SingleR_", srx, "_Violin.pdf"), w=12, h=7)
  
  # Clusters - Scores
  # pDT2 <- pDT[,.(mean=mean(value)), by=c("Clusters", "variable")]
  # pDT2[,scaleMean := scale(mean), by=c("Clusters")]
  # for(mx in c("mean")){
  #   pDT3 <- copy(pDT2)
  #   #pDT3 <- hierarch.ordering(pDT3, toOrder = "Clusters", orderBy = "variable", value.var = mx)
  #   pDT3 <- hierarch.ordering(pDT3, toOrder = "variable", orderBy = "Clusters", value.var = mx)
  #   ggplot(pDT3, aes_string(x="Clusters", y="variable", fill=mx)) + 
  #     theme_bw(12) +
  #     geom_tile() +
  #     xRot() +
  #     scale_fill_gradient(low="white", high="blue")
  #   ggsave(out("SingleR_", srx, "_Clusters_", mx, ".pdf"), 
  #          h=length(unique(pDT3$variable)) * 0.3+1,
  #          w=length(unique(pDT3$Clusters)) * 0.3+2
  #          )
  # }
  
  # UMAP
  # n <- length(unique(pDT$variable))
  # ggplot(pDT, aes(x=UMAP1, y=UMAP2)) + 
  #   theme_bw(12) + 
  #   stat_summary_hex(aes(z=value),fun=mean) +
  #   scale_fill_gradient(low="white", high="blue") +
  #   facet_wrap(~variable, ncol=5) +
  #   theme_bw(12)
  # ggsave(out("SingleR_", srx, ".pdf"), w=5*2+2, h=ceiling(n/5)*2+1, limitsize = FALSE)
  
  # Clusters - Predictions
  pDT.ann <- merge(singleR.resX[,c("pruned_labels", "cell")],ann, by.x="cell", by.y="rn")
  pDT.ann <- pDT.ann[,.N, by=c("Clusters", "pruned_labels")]
  pDT.ann[,sum := sum(N), by="Clusters"]
  pDT.ann[,percent := N/sum*100]
  ggplot(pDT.ann, aes(x=factor(as.numeric(Clusters)), y=pruned_labels, fill=percent)) + 
    geom_tile() +
    scale_fill_gradient(low="white", high="red")
  ggsave(out("SingleR_", srx, "_Clusters_", "PercPredicted", ".pdf"),          
         h=length(unique(pDT.ann$pruned_labels)) * 0.3+1,
         w=length(unique(pDT.ann$Clusters)) * 0.3+2)
  
}


# ANTIBODIES --------------------------------------------------------------
abMT <- additional.info$CITESEQ2$`Antibody Capture`
abMT <- SCRNA.TPXToLog(SCRNA.RawToTPX(abMT, scale.factor = 1e6))
ann.c1 <- ann[sample == "CITESEQ2"]

# UMAP
res <- data.table()
abx <- row.names(abMT)[1]
for(abx in row.names(abMT)){
  pDT <- copy(ann.c1)
  pDT$Signal <- abMT[abx, gsub("_CITESEQ2", "", ann.c1$rn)]
  pDT$Antibody <- abx
  res <- rbind(res, pDT)
}
res[,Signal.norm := scale(Signal), by="Antibody"]
ggplot(res, aes(x=UMAP1, y=UMAP2)) + 
  stat_summary_hex(aes(z=Signal.norm),fun=mean) +
  scale_fill_gradient2(low="blue", high="red") +
  facet_wrap(~Antibody) +
  theme_bw(12)
ggsave(out("Antibodies_UMAP.pdf"), w=12+2, h=9+1)

# Correlations
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

# Percentile plots
probx=0.9
res[,percentile := quantile(Signal.norm, probs=probx, na.rm=TRUE), by="Antibody"]
resN <- res[, sum(Signal.norm > percentile), by=c("Clusters", "Antibody")]
resN[,clSize := sum(V1), by="Clusters"]
stopifnot(all(resN[, length(unique(clSize)), by="Clusters"]$V1 == 1))
resN[,percentage := V1/clSize*100]
resN <- hierarch.ordering(resN, toOrder = "Clusters", orderBy = "Antibody", value.var = "percentage")
resN <- hierarch.ordering(resN, toOrder = "Antibody", orderBy = "Clusters", value.var = "percentage")
resN[is.na(percentage), percentage := 0]
ggplot(resN, aes(x=Clusters, y=Antibody, fill=percentage)) + 
  geom_tile() + 
  ggtitle("Percent of cell in Clusters\nthat are above 90th percentile of antibody signal") + 
  scale_fill_gradient(low="white", high="blue")
ggsave(out("Antibodies_Percentile.pdf"), w=5, h=4)



# GUIDES ------------------------------------------------------

# . Number of guides --------------------------------------------------------
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
ggsave(out("Guides_Counts_PercAssigned.pdf"), w=6,h=4)

# . Guides in UMAP ----------------------------------------------------------
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

# . Guides per cluster ------------------------------------------------------
res <- data.table()
pDT1 <- ann[!is.na(mixscape_class)][mixscape_class.global != "NP"]
for(sx in unique(pDT1$sample)){
  pDT2 <- pDT1[sample == sx]
  for(gx in unique(pDT2[mixscape_class != "NTC"]$mixscape_class)){
    for(cx in unique(pDT2$Clusters)){
      message(gx, "-", cx)
      mx <- as.matrix(with(pDT2[mixscape_class %in% c(gx, "NTC")], table(Clusters == cx, mixscape_class == gx)))
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



# DE for GUIDES -----------------------------------------
obj.de <- monocle.obj
obj.de <- obj.de[,colData(obj.de)$mixscape_class.global %in% c("KO", "NTC")]
obj.de <- obj.de[rowSums(counts(obj.de)) > 20,]
colData(obj.de)$ClusterDE <- obj.de@clusters$UMAP$clusters[row.names(colData(obj.de))]
colData(obj.de)$GuideDE <- gsub("_.+$", "", colData(obj.de)$guide)
colData(obj.de)$GuideDE <- factor(colData(obj.de)$GuideDE, levels=c("NTC", setdiff(unique(colData(obj.de)$GuideDE), "NTC")))

de.file <- out("DEG_Results.RData")
if(file.exists(de.file)){
  load(de.file)
} else {
  x <- fit_models(obj.de, model_formula_str = "~GuideDE + ClusterDE + sample", verbose=TRUE, cores=10)
  res <- data.table(coefficient_table(x), keep.rownames = TRUE)
  res <- res[,-c("model", "model_summary", "rn", "gene_short_name", "model_component"),with=F]
  save(res, file=de.file)
}

#  . Export results -------------------------------------------------------
resGuides <- res[grepl("GuideDE", term)]
resGuides[,term := gsub("GuideDE", "", term)]
write.tsv(resGuides, file=out("DEG_Results.tsv"))


#  . Plot top genes -------------------------------------------------------

# Estimate and P-value
gg <- resGuides[q_value < 0.05][order(-abs(estimate))][,head(.SD,n=10), by="term"]$gene_id
pDT <- resGuides[gene_id %in% gg]
pDT <- hierarch.ordering(pDT, toOrder="gene_id", orderBy = "term", value.var = "estimate")
pDT <- hierarch.ordering(pDT, toOrder="term", orderBy = "gene_id", value.var = "estimate")
ggplot(pDT, aes(x=term, y=gene_id, size=pmin(5, -log10(q_value)), color=sign(estimate) * pmin(5,abs(estimate)))) + 
  scale_size_continuous(name="padj", range=c(0,5)) + 
  scale_color_gradient2(name="delta", low="blue", high="red") +
  theme_bw(12) +
  geom_point() +
  xRot()
ggsave(out("DEG_examples.pdf"), w=6,h=15)

# Dotplots
ll <- with(resGuides[q_value < 0.05][order(-abs(estimate))][,head(.SD,n=10), by="term"], split(gene_id, term))
lx <- "Chd4"
for(lx in names(ll)){
  pDT <- DotPlotData(obj.de, markers=ll[[lx]], cols=c("GuideDE", "ClusterDE", "sample"))
  #pDT[GuideDE == lx & Gene == "Fyb2"]
  #0.075 * 40
  ggplot(pDT, aes(x=ClusterDE, y=sample, color=mean, size=percentage)) + 
    theme_bw(12) +
    geom_point() + 
    scale_size_continuous(range=c(0,5)) +
    scale_color_gradient(low="black", high="red") +
    facet_grid(Gene ~ GuideDE) +
    xRot() +
    ggtitle(lx)
  ggsave(out("DEG_examples_Dotplot_",lx,".pdf"), w=30,h=length(unique(pDT$Gene)) * length(unique(pDT$sample)) * 0.15 + 1)
}

# Checking one specific gene
colData(obj.de)$GroupDE <- paste(colData(obj.de)$GuideDE, colData(obj.de)$sample, colData(obj.de)$ClusterDE)
samples.ko <- ann[grepl("Chd4", guide) & mixscape_class.global != "NP"]$rn
samples.NTC <- ann[mixscape_class.global == "NTC"]$rn
plot_genes_violin(cds_subset = obj.de["Fyb2",c(samples.ko, samples.NTC)], group_cells_by = "GroupDE", normalize = TRUE, log_scale = FALSE, pseudocount = 0) + xRot()
ggsave(out("DEG_Fyb2.pdf"), w=12,h=4)
# counts(obj.de["Fyb2", ann[sample == "ECCITE4_Cas9" & grepl("Chd4", guide) & mixscape_class.global != "NP" & Clusters == 2]$rn])
# pDT[GuideDE == lx & Gene == "Fyb2"]



# COR of DEG --------------------------------------------------------------
umapMT <- toMT(dt=resGuides, row = "gene_id", col = "term", val = "estimate")
colnames(umapMT) <- gsub("GuideDE", "", colnames(umapMT))
cMT <- corS(umapMT)
gn <- ncol(umapMT)
diag(cMT) <- NA
cleanDev(); pdf(out("DEG_CorHM.pdf"),w=gn/6+2, h=gn/6+1.5)
pheatmap(cMT)#, breaks=seq(-1,1, 0.01), color=colorRampPalette(c("#6a3d9a", "#a6cee3", "white", "#fdbf6f", "#e31a1c"))(200))
dev.off()



# UMAP of DEG -----------------------------------------

# UMAP
umapMT <- toMT(dt=resGuides, row = "gene_id", col = "term", val = "estimate")
colnames(umapMT) <- gsub("GuideDE", "", colnames(umapMT))

cMT <- corS(umapMT)
diag(cMT) <- NA
pheatmap(cMT)

set.seed(1212)
umap.res <- umap(umapMT)
umap <- data.table(umap.res$layout, keep.rownames = TRUE)
ggplot(umap, aes(x=V1, y=V2)) + geom_hex() + theme_bw(12)
ggsave(out("RegulatoryMap_UMAP.pdf"), w=6,h=5)

# Cluster
set.seed(1212)
idx <- umap.res$knn$indexes
g <- do.call(rbind, apply(idx[, 2:ncol(idx)], 2, function(col){data.table(row.names(idx)[col], row.names(idx)[idx[,1]])}))
(g <- graph.edgelist(as.matrix(g),directed=FALSE))
cl <- cluster_walktrap(g)
clx <- setNames(cl$membership, V(g)$name)
umap$Cluster <- clx[umap$rn]
ggplot(umap, aes(x=V1, y=V2, color=factor(Cluster))) + geom_point() + theme_bw(12)
ggsave(out("RegulatoryMap_UMAP_Clusters.pdf"), w=6,h=5)

# Export annotation
umap <- setNames(umap, c("Gene", "UMAP1", "UMAP2", "Cluster"))
write.tsv(umap, out("RegulatoryMap_UMAP.tsv"))

# Plot estimates on UMAP
pUMAP.de <- merge(umap, resGuides[,c("gene_id", "term", "estimate")], by.x="Gene", by.y="gene_id")
summary.function <- function(x){ret <- mean(x);return(min(5, abs(ret)) * sign(ret))}
dim.umap1 <- floor(max(abs(pUMAP.de$UMAP1))) + 0.5
dim.umap2 <- floor(max(abs(pUMAP.de$UMAP2))) + 0.5
tn <- length(unique(pUMAP.de$term))
ggplot(pUMAP.de, aes(x=UMAP1, y=UMAP2)) + 
  stat_summary_hex(
    aes(z=estimate),
    fun=summary.function) +
  scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") + 
  facet_wrap(~gsub("GuideDE", "", term), ncol = 5) + theme_bw(12) +
  xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
  xlim(-dim.umap1,dim.umap1) + ylim(-dim.umap2,dim.umap2)
ggsave(out("RegulatoryMap_UMAP_Values.pdf"), w=11,h=ceiling(tn/5) * 2 + 1)
