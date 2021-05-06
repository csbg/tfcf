source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("ECCITE1_03_ClusterEnrichment/")

clusters <- fread(PATHS$ECCITE1$DATA$clusters)
guides <- fread(PATHS$ECCITE1$DATA$guides)
names(guides) <- gsub("^(.)", "\\U\\1", names(guides), perl = T)


dat <- merge(clusters, guides, by=c("Barcode"))


res <- data.table()
for(gx in unique(dat[!grepl(" ", Guide)]$Guide)){
  for(cx in unique(dat$Cluster)){
    fish <- fisher.test(as.matrix(with(dat, table(Cluster == cx, Guide == gx))))
    res <- rbind(res, data.table(Cluster=cx, Guide=gx, p=fish$p.value, OR=fish$estimate))
  }
}

res[,padj := p.adjust(p, method="BH")]
res[padj < 0.05]

ggplot(res, aes(
      x=factor(Cluster),
      y=Guide, 
      color=log2(OR + min(res[OR != 0]$OR)), 
      size=pmin(-log10(padj), 5))
      ) + 
  geom_point(shape=16) +
  scale_color_gradient2(name="log2OR", low="blue", high="red") +
  scale_size_continuous(name="padj") + 
  theme_bw(12)
ggsave(out("Fisher.pdf"), w=5, h=5)