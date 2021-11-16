source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("CITESEQ1_01_InitialAnalysis/")

data <- SCRNA.read_10Xh5.610(PATHS$CITESEQ1$DATA$matrix)
umap <- fread(PATHS$CITESEQ1$DATA$umap, check.names = TRUE)
clusters <- fread(PATHS$CITESEQ1$DATA$clusters, check.names = TRUE)
de.genes <- fread(PATHS$CITESEQ1$DATA$de, check.names = TRUE)

aDT <- merge(umap, clusters, by="Barcode")
stopifnot(all(aDT$Barcode == colnames(data$matrix)))


exMT <- data$matrix[data$features[feature_type == "Gene Expression"]$id,]
aDT$LibrarySize <- Matrix::colSums(exMT)
aDT$MitochondriaPercent <- Matrix::colSums(exMT[data$features[grepl("^MT-", name, ignore.case = TRUE)]$id,])
aDT[,MitochondriaPercent := MitochondriaPercent/LibrarySize * 100]
exMT <- SCRNA.TPXToLog(SCRNA.RawToTPX(exMT, scale.factor = 1e6))

abMT <- data$matrix[data$features[feature_type == "Antibody Capture"]$id,]
aDT$LibrarySizeAb <- Matrix::colSums(abMT)
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
ggsave(out("UMAP_Signals.pdf"), w=12+2, h=9+1)

ggplot(res, aes(x=factor(Cluster), y=Signal)) + 
  geom_violin(fill="lightblue") +
  facet_grid(Antibody ~ ., scales = "free") +
  theme_bw(12)
ggsave(out("Cluster_Signals.pdf"), w=6, h=16)


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
ggsave(out("UMAP_ExSignals.pdf"), w=12+2, h=9+1)

ggplot(res, aes(x=factor(Cluster), y=Signal)) + 
  geom_violin(fill="lightblue") +
  facet_grid(Antibody ~ ., scales = "free") +
  theme_bw(12)
ggsave(out("Cluster_ExSignals.pdf"), w=6, h=16)


ggplot(aDT, aes(x=UMAP.1, y=UMAP.2)) + 
  stat_summary_hex(aes(z=LibrarySize),fun=mean) +
  scale_fill_gradient(low="blue", high="red") +
  theme_bw(12) +
  ggtitle("Library size (genes)")
ggsave(out("UMAP_LibrarySizeGeneEx.pdf"), w=5, h=4)

ggplot(aDT, aes(x=factor(Cluster), y=LibrarySize)) + 
  geom_violin() + 
  scale_fill_gradient(low="blue", high="red") +
  scale_y_log10() + 
  theme_bw(12) +
  ggtitle("Library size (genes)")
ggsave(out("Cluster_LibrarySizeGeneEx.pdf"), w=5, h=4)

ggplot(aDT, aes(x=UMAP.1, y=UMAP.2)) + 
  stat_summary_hex(aes(z=LibrarySizeAb),fun=mean) +
  scale_fill_gradient(low="blue", high="red") +
  theme_bw(12) +
  ggtitle("Library size (antibodies)")
ggsave(out("UMAP_LibrarySizeAbs.pdf"), w=5, h=4)

ggplot(aDT, aes(x=factor(Cluster), y=LibrarySizeAb)) + 
  geom_violin() + 
  scale_fill_gradient(low="blue", high="red") +
  scale_y_log10() + 
  theme_bw(12) +
  ggtitle("Library size (antibodies)")
ggsave(out("Cluster_LibrarySizeGeneAbs.pdf"), w=5, h=4)

ggplot(aDT, aes(x=UMAP.1, y=UMAP.2)) + 
  stat_summary_hex(aes(z=MitochondriaPercent),fun=mean) +
  scale_fill_gradient(low="blue", high="red") +
  theme_bw(12) +
  ggtitle("Mitochondrial reads (all cells)")
ggsave(out("UMAP_Mitochondria.pdf"), w=5, h=4)

ggplot(aDT, aes(x=factor(Cluster), y=MitochondriaPercent+0.1)) + 
  geom_violin() + 
  theme_bw(12) +
  scale_y_log10() + 
  ggtitle("Mitochondrial reads (all cells)")
ggsave(out("Cluster_Mitochondria.pdf"), w=5, h=4)

ggplot(aDT[MitochondriaPercent < 10], aes(x=UMAP.1, y=UMAP.2)) + 
  stat_summary_hex(aes(z=MitochondriaPercent),fun=mean) +
  scale_fill_gradient(low="blue", high="red") +
  theme_bw(12) +
  ggtitle("Mitochondrial reads\n(only cells with < 10% are plotted)")
ggsave(out("UMAP_Mitochondria2.pdf"), w=5, h=4)

ggplot(aDT, aes(x=UMAP.1, y=UMAP.2)) + 
  geom_hex() +
  scale_fill_gradient(low="lightgrey", high="blue") +
  theme_bw(12) +
  geom_label(data=aDT[,.(UMAP.1=median(UMAP.1), UMAP.2=median(UMAP.2)), by="Cluster"], aes(label=Cluster)) + 
  ggtitle("Number of cells")
ggsave(out("UMAP_Cluster.pdf"), w=5, h=4)


sum(aDT[,.N, by="Cluster"][Cluster %in% 1:3]$N)/nrow(aDT)

cells.keep <- aDT[LibrarySize > 1000 & LibrarySizeAb > 1000 & MitochondriaPercent < 10]
write.table(cells.keep, out("Cells_Keep.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE)

pDT <- merge(cells.keep[,.N, by="Cluster"], aDT[,.N, by="Cluster"], by="Cluster", all=TRUE, suffixes = c("_keep", "_orig"))
pDT <- melt(pDT, id.vars = "Cluster")
ggplot(pDT, aes(x=factor(Cluster), y=value, fill=variable)) + geom_bar(stat="identity", position="dodge")
ggsave(out("Removal_Numbers.pdf"), w=6,h=4)
