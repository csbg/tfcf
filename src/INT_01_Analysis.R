source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("INT_01_Analysis/")

clusters.c1 <- fread(PATHS$CITESEQ1_CLEAN$DATA$clusters, check.names = TRUE)
clusters.e1 <- fread(PATHS$ECCITE1$DATA$clusters, check.names = TRUE)


data.c1 <- SCRNA.read_10Xh5.610(PATHS$CITESEQ1_CLEAN$DATA$matrix)
mt.c1 <- data.c1$matrix[data.c1$features[feature_type == "Gene Expression"]$id,]
mt.c1 <- SCRNA.TPXToLog(SCRNA.RawToTPX(mt.c1, scale.factor = 1e6))

data.e1 <- SCRNA.read_10Xh5.610(PATHS$ECCITE1$DATA$matrix)
mt.e1 <- SCRNA.TPXToLog(SCRNA.RawToTPX(data.e1$matrix, scale.factor = 1e6))

cl.mean.c1 <- sapply(with(clusters.c1, split(Barcode, Cluster)), function(cxb){Matrix::rowMeans(mt.c1[,cxb])})
colnames(cl.mean.c1) <- paste("CITE", colnames(cl.mean.c1))

cl.mean.e1 <- sapply(with(clusters.e1, split(Barcode, Cluster)), function(cxb){Matrix::rowMeans(mt.e1[,cxb])})
colnames(cl.mean.e1) <- paste("ECCITE", colnames(cl.mean.e1))

stopifnot(all(row.names(cl.mean.c1) == row.names(cl.mean.e1)))


# Cluster correlations ----------------------------------------------------
keep <- rowSums(cl.mean.c1) != 0 & rowSums(cl.mean.e1) != 0
cMT <- corS(cl.mean.c1[keep, ], cl.mean.e1[keep, ], use="pairwise.complete.obs")
cleanDev(); pdf(out("Cluster_Correlation.pdf"), w=6,h=5)
pheatmap(cMT)
dev.off()

cMT <- corS(cbind(cl.mean.c1[keep, ], cl.mean.e1[keep, ]))
dd <- as.dist(1-cMT)
cleanDev(); pdf(out("Cluster_Correlation_All.pdf"), w=8,h=8)
pheatmap(cMT, clustering_distance_rows = dd, clustering_distance_cols = dd)
dev.off()

x <- data.table(cbind(cl.mean.c1[keep, ], cl.mean.e1[keep, ]),keep.rownames = TRUE, check.names = TRUE)
x <- melt(x, id.vars = c("rn", grep("ECCITE", colnames(x), value = TRUE)))
x <- melt(x, id.vars = c("rn", "variable", "value"))
ggplot(x, aes(x=value, y=value.1)) + 
  stat_binhex(aes(fill=log10(..count..))) + 
  facet_grid(variable ~ variable.1) + 
  theme_bw(12)
ggsave(out("Cluster_Correlation_hex.pdf"), width=20, height=20)



# Individual cells --------------------------------------------------------
gg <- names(which(rowSums(cl.mean.c1) != 0))
scCor <- corS(as.matrix(mt.e1[gg,]), cl.mean.c1[gg,])
pDT <- merge(melt(
    data.table(scCor, keep.rownames = TRUE), 
    id.vars = "rn", variable.name = "Cluster"
    ), 
  clusters.e1, 
  by.x="rn", by.y="Barcode", suffixes = c("_c1", "_e1")
  )
ggplot(pDT[, .(mean = mean(value), sd=sd(value)), by=c("Cluster_e1", "Cluster_c1")], 
       aes(x=factor(Cluster_e1), y=factor(Cluster_c1), color=mean, size=1/sd)) + 
  geom_point() +
  scale_color_distiller() + 
  theme_bw(12)

ggplot(pDT[, .(mean = mean(value), sd=sd(value)), by=c("Cluster_e1", "Cluster_c1")], 
       aes(x=factor(Cluster_e1), y=factor(Cluster_c1), fill=mean)) + 
  geom_tile() +
  theme_bw(12)


ggplot(pDT, aes(x=factor(Cluster_e1), y=value)) + 
  geom_violin() +
  facet_grid(Cluster_c1 ~ .) +
  theme_bw(12)
ggsave(out("Violin.pdf"), w=15,h=15)
