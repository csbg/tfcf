source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("CITESEQ1_01_InitialAnalysis/")

data <- SCRNA.read_10Xh5.610(PATHS$CITESEQ1$DATA$matrix)
umap <- fread(PATHS$CITESEQ1$DATA$umap, check.names = TRUE)
clusters <- fread(PATHS$CITESEQ1$DATA$clusters, check.names = TRUE)

aDT <- merge(umap, clusters, by="Barcode")

dMT <- data$matrix
dMT.raw <- dMT
dMT <- SCRNA.TPXToLog(SCRNA.RawToTPX(dMT, scale.factor = 1e6))
gDT <- data$features

stopifnot(all(aDT$Barcode == colnames(dMT)))
aDT$LibrarySize <- Matrix::colSums(dMT.raw)



res <- data.table()
for(idx in unique(gDT[feature_type == "Antibody Capture"]$id)){
  pDT <- copy(aDT)
  pDT$Signal <- dMT[idx,]
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


ggplot(aDT, aes(x=UMAP.1, y=UMAP.2)) + 
  stat_summary_hex(aes(z=LibrarySize),fun=mean) +
  scale_fill_gradient(low="blue", high="red") +
  theme_bw(12)
ggsave(out("UMAP_LibrarySize.pdf"), w=5, h=4)

