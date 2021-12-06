source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("CITESEQ1_03_Analysis/")

data <- SCRNA.read_10Xh5.610(PATHS$CITESEQ1_CLEAN$DATA$matrix)
umap <- fread(PATHS$CITESEQ1_CLEAN$DATA$umap, check.names = TRUE)
clusters <- fread(PATHS$CITESEQ1_CLEAN$DATA$clusters, check.names = TRUE)
de.genes <- fread(PATHS$CITESEQ1_CLEAN$DATA$de, check.names = TRUE)

aDT <- merge(umap, clusters, by="Barcode")
stopifnot(all(aDT$Barcode == colnames(data$matrix)))


exMT <- data$matrix[data$features[feature_type == "Gene Expression"]$id,]
stopifnot(all(diff(exMT@p)[1:10] == apply(exMT[,1:10], 2, function(col) sum(col!=0))))
aDT$nrGenes <- diff(exMT@p)
exMT <- SCRNA.TPXToLog(SCRNA.RawToTPX(exMT, scale.factor = 1e6))

abMT <- data$matrix[data$features[feature_type == "Antibody Capture"]$id,]
abMT <- SCRNA.TPXToLog(SCRNA.RawToTPX(abMT, scale.factor = 1e6))

gDT <- data$features



res <- data.table()
for(idx in unique(gDT[feature_type == "Antibody Capture"]$id)){
  pDT <- copy(aDT)
  pDT$Signal <- abMT[idx,]
  pDT$Antibody <- idx
  res <- rbind(res, pDT)
}

res[,Signal.norm := scale(Signal), by="Antibody"]

ggplot(res, aes(x=UMAP.1, y=UMAP.2)) + 
  stat_summary_hex(aes(z=Signal.norm),fun=mean) +
  scale_fill_gradient2(low="blue", high="red") +
  facet_wrap(~Antibody) +
  theme_bw(12)
ggsave(out("Antibodies_UMAP.pdf"), w=12+2, h=9+1)

cMT <- corS(t(as.matrix(abMT)))
diag(cMT) <- NA
dist <- as.dist(1-cMT)
cleanDev(); pdf(out("Antibodies_Correlation.pdf"), w=5,h=4)
pheatmap(cMT,
         clustering_distance_rows = dist,
         clustering_distance_cols = dist,
         breaks=seq(-1,1, 0.01), color=COLORS.HM.FUNC(200),
         )
dev.off()

ggplot(res, aes(x=factor(Cluster), y=Signal)) + 
  geom_violin(fill="lightblue") +
  facet_grid(Antibody ~ ., scales = "free") +
  theme_bw(12)
ggsave(out("Antibodies_Clusters.pdf"), w=6, h=16)

probx=0.9
res[,percentile := quantile(Signal.norm, probs=probx), by="Antibody"]
resN <- res[, sum(Signal.norm > percentile), by=c("Cluster", "Antibody")]
resN[,clSize := sum(V1), by="Cluster"]
stopifnot(all(resN[, length(unique(clSize)), by="Cluster"]$V1 == 1))
resN[,percentage := V1/clSize*100]
resN <- hierarch.ordering(resN, toOrder = "Cluster", orderBy = "Antibody", value.var = "percentage")
resN <- hierarch.ordering(resN, toOrder = "Antibody", orderBy = "Cluster", value.var = "percentage")
ggplot(resN, aes(x=Cluster, y=Antibody, fill=percentage)) + 
  geom_tile() + 
  ggtitle("Percent of cell in cluster\nthat are above 90th percentile of antibody signal") + 
  scale_fill_gradient(low="white", high="blue")
ggsave(out("Antibodies_Percentile.pdf"), w=5, h=4)


res <- data.table()
for(idx in unique(c("Elane", "Itgam", "F13a1", "Gata1", "Itga2b"))){
  pDT <- copy(aDT)
  pDT$Signal <- exMT[gDT[name  == idx]$id,]
  pDT$Antibody <- idx
  res <- rbind(res, pDT)
}

res[,Signal.norm := scale(Signal), by="Antibody"]

ggplot(res, aes(x=UMAP.1, y=UMAP.2)) + 
  stat_summary_hex(aes(z=Signal.norm),fun=mean) +
  scale_fill_gradient2(low="blue", high="red") +
  facet_wrap(~Antibody) +
  theme_bw(12)
ggsave(out("Expression_UMAP.pdf"), w=12+2, h=9+1)

ggplot(res, aes(x=factor(Cluster), y=Signal)) + 
  geom_violin(fill="lightblue") +
  facet_grid(Antibody ~ ., scales = "free") +
  theme_bw(12)
ggsave(out("Expression_Clusters.pdf"), w=6, h=16)

ggplot(aDT, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_hex() +
  scale_fill_gradient(low="lightgrey", high="blue") +
  theme_bw(12) +
  geom_label(data=aDT[,.(UMAP.1=median(UMAP.1), UMAP.2=median(UMAP.2)), by="Cluster"], aes(label=Cluster)) + 
  ggtitle("Number of cells")
ggsave(out("Clusters_UMAP.pdf"), w=5, h=4)

ggplot(aDT, aes(x=factor(Cluster), y=nrGenes)) + 
  geom_violin(fill="lightblue") +
  theme_bw(12)
ggsave(out("NrGenes_Clusters.pdf"), w=5, h=4)

ggplot(aDT, aes(x=nrGenes)) + 
  geom_density() +
  theme_bw(12)
ggsave(out("NrGenes_Total.pdf"), w=5, h=4)