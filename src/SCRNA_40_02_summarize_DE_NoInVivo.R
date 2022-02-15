source("src/00_init.R")

base.dir <- "SCRNA_40_02_DE_summary_NoInVivo/"
out <- dirout(base.dir)


require(ggrepel)
require(umap)
require(igraph)


# FUNCTIONS ---------------------------------------------------------------
ds <- function(path){load(path); return(monocle.obj)}


# EnrichR -----------------------------------------------------------------
(load(PATHS$RESOURCES$Enrichr.mouse))


# Load data ---------------------------------------------------------------
ff <- list.files(dirout_load("SCRNA_30_DE_Nebula")(""), full.names = TRUE)
names(ff) <- basename(ff)
ff <- ff[!grepl("in.vivo", ff)]
ff <- lapply(ff, function(x) paste0(x, "/DEG_Results_all.tsv"))
ff <- lapply(ff, fread)
DE.RES <- rbindlist(ff, idcol = "tissue")
stopifnot(all(grepl("^GuideDE", DE.RES$term)))
write.tsv(DE.RES[, -c("se", "convergence", "estimate_raw", "term"), with=F], out("DEG_Statistics.tsv"))
write.tsv(DE.RES[, -c("se", "convergence", "estimate_raw", "term"), with=F][q_value < 0.05], out("DEG_Statistics_significant.tsv"))

DE.RES[, id := paste(guide, make.names(tissue))]
FC.MT <- toMT(dt = DE.RES, col = "id", row = "gene_id", val = "estimate")
write.csv(round(FC.MT, 3), out("FoldChanges.csv"))


# Correlations ------------------------------------------------------------
cMT <- corS(FC.MT[apply(is.na(FC.MT), 1, sum) == 0,], use="pairwise.complete.obs")
gn <- ncol(FC.MT)
dd <- as.dist(1-cMT)

write.csv(cMT, out("DEG_Cor.csv"))

diag(cMT) <- NA
cleanDev(); pdf(out("DEG_CorHM.pdf"),w=gn/6+2, h=gn/6+1.5)
pheatmap(cMT,clustering_distance_rows = dd, clustering_distance_cols = dd)
dev.off()

cleanDev(); pdf(out("DEG_CorHM_Colors.pdf"),w=gn/6+2, h=gn/6+1.5)
lim <- max(abs(cMT), na.rm=TRUE)
pheatmap(cMT,
         clustering_distance_rows = dd, clustering_distance_cols = dd,
         breaks=seq(-lim,lim, length.out = 200),
         color=colorRampPalette(c("#6a3d9a", "#a6cee3", "white", "#fdbf6f", "#e31a1c"))(200))
dev.off()

# pDT <- melt(data.table(cMT, keep.rownames = T), id.vars = "rn")
# pDT[, gene1 := gsub(" .+", "", rn)]
# pDT[, gene2 := gsub(" .+", "", variable)]
# ggplot(pDT, aes(x=rn, y=variable, fill=value)) + 
#   theme_bw() +
#   geom_tile() + 
#   facet_grid(gene2 ~ gene1, scales = "free", space = "free")
set.seed(1212)
mds.res <- data.table(cmdscale(d=as.dist(1-cMT), k=2), keep.rownames=TRUE)
mds.res <- cbind(mds.res, setNames(data.table(do.call(rbind, strsplit(mds.res$rn, " "))), c("gene", "tissue")))
ggplot(mds.res, aes(x=V1, y=V2, color=tissue, label=gene)) + 
  geom_point(size=1) + 
  ggrepel::geom_text_repel()+#color="black") + 
  theme_bw() +
  xlab("MDS dimension 1") +
  ylab("MDS dimension 2")
ggsave(out("Correlation_Heatmap_MDS.pdf"), w=8, h=7)




# UMAP of genes -----------------------------------------------------------
umap.log2FC.cutoff <- 3
umapMT <- FC.MT
umap.type.name <- "all"
gg <- if(umap.type.name == "top") unique(DE.RES[q_value < 0.05 & abs(estimate) > 1]$gene_id) else row.names(umapMT)
mt <- umapMT[gg,]
mt[mt > umap.log2FC.cutoff] <- umap.log2FC.cutoff
mt[mt < -umap.log2FC.cutoff] <- -umap.log2FC.cutoff
mt[is.na(mt)] <- 0
set.seed(1212)
umap.res <- umap(mt)
umap <- data.table(umap.res$layout, keep.rownames = TRUE)
umap <- setNames(umap, c("Gene", "UMAP1", "UMAP2"))
ggplot(umap, aes(x=UMAP1, y=UMAP2)) + geom_hex(bins=100) + theme_bw(12)
ggsave(out("RegulatoryMap_UMAP_",umap.type.name,".pdf"), w=6,h=5)

# Cluster
idx <- umap.res$knn$indexes
g <- do.call(rbind, apply(idx[, 2:ncol(idx)], 2, function(col){data.table(row.names(idx)[col], row.names(idx)[idx[,1]])}))
(g <- graph.edgelist(as.matrix(g),directed=FALSE))
set.seed(1234)
cl <- cluster_walktrap(g)
clx <- setNames(cl$membership, V(g)$name)
umap$Cluster <- clx[umap$Gene]
# POInts
ggplot(umap, aes(x=UMAP1, y=UMAP2, color=factor(Cluster))) + 
  geom_point() + 
  theme_bw(12) +
  geom_label(data=umap[, .(UMAP1=median(UMAP1), UMAP2=median(UMAP2)), by="Cluster"], aes(label=Cluster))
ggsave(out("RegulatoryMap_UMAP_",umap.type.name,"_Clusters.pdf"), w=6,h=5)
#HEX
ggplot(umap, aes(x=UMAP1, y=UMAP2)) + 
  geom_hex(bins=100) +
  theme_bw(12) +
  geom_label(data=umap[, .(UMAP1=median(UMAP1), UMAP2=median(UMAP2)), by="Cluster"], aes(label=Cluster))
ggsave(out("RegulatoryMap_UMAP_",umap.type.name,"_Clusters_hex.pdf"), w=6,h=5)

# Export annotation
write.tsv(umap, out("RegulatoryMap_UMAP_",umap.type.name,".tsv"))


# Plot logFC on UMAP ---------------------------------------------------
pUMAP.de <- merge(umap, setNames(melt(data.table(umapMT, keep.rownames = TRUE), id.vars = "rn"), c("gene_id", "term", "estimate")), by.x="Gene", by.y="gene_id")
dim.umap1 <- floor(max(abs(pUMAP.de$UMAP1))) + 0.5
dim.umap2 <- floor(max(abs(pUMAP.de$UMAP2))) + 0.5
pUMAP.de[, guide := gsub(" .+", "", term)]
pUMAP.de[, tissue := gsub(".+? ", "", term)]
tn <- length(unique(pUMAP.de$guide))
pUMAP.de[, estimate_cap := pmin(umap.log2FC.cutoff, abs(estimate)) * sign(estimate)]

# Plot estimates on UMAP
ggplot(pUMAP.de, aes(x=UMAP1, y=UMAP2)) +
  stat_summary_hex(
    aes(z=estimate_cap),
    fun=function(x){mean(x, na.rm=TRUE)}) +
  scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") +
  facet_grid(guide~tissue) + theme_bw(12) +
  xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
  xlim(-dim.umap1,dim.umap1) + ylim(-dim.umap2,dim.umap2)
ggsave(out("RegulatoryMap_UMAP_",umap.type.name,"_Values.pdf"), w=2*3+2,h=tn * 2 + 1, limitsize=FALSE)

# values by cluster
pDT <- pUMAP.de[, mean(estimate_cap, na.rm=TRUE), by=c("Cluster", "term")]
pDT[, gene := gsub(" .+", "", term)]
pDT[, tissue := gsub(".+? ", "", term)]
pDT <- hierarch.ordering(pDT, toOrder = "Cluster", orderBy = "term", value.var = "V1")
pDT <- hierarch.ordering(pDT, orderBy = "Cluster", toOrder = "gene", value.var = "V1", aggregate = TRUE)
ggplot(pDT, aes(y=factor(Cluster), x=tissue, fill=V1)) + 
  theme_bw(12) + 
  geom_tile() +
  scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") +
  facet_grid(. ~ gene, scales = "free", space = "free", switch = "x") +
  xlab("Gene modules (Gene-UMAP Clusters)") +
  theme(strip.text.x = element_text(angle=90)) +
  xRot()
ggsave(out("RegulatoryMap_UMAP_",umap.type.name,"_ClusterValues.pdf"), w=15,h=5, limitsize=FALSE)


# Enrichment in gene clusters ---------------------------------------------
fish.genes <- with(umap, split(Gene, Cluster))
fish.res <- fisher.test.enrichment(
  geneSets = enr.terms,
  gene.list = fish.genes,
  cores=10,
  bg=unique(umap$Gene))
fish.res[, padj := p.adjust(pval, method="BH")]
write.tsv(fish.res, out("UMAP_GSEA.tsv"))

# Plots
fish.res[, db := database]
fish.res[, pathway := geneset]
fish.res[, log2OR := log2(pmin(10, pmax(1/10, oddsRatio)))]
fish.res[db == "ChromatinFactors", pathway := gsub("^(.+?)_(.+)$", "\\2 \\1", pathway)]
dbx <- fish.res$db[1]
for(dbx in unique(fish.res$db)){
  pDT <- fish.res[db == dbx]
  pwx <- unique(pDT[padj < 0.1][order(-log2OR)][,head(.SD,n=3), by="list"]$pathway)
  pDT <- pDT[pathway %in% pwx]
  ggplot(pDT, aes(y=factor(as.numeric(list)), x=pathway, size=pmin(5, -log10(padj)), color=log2OR)) +
    theme_bw(12) +
    #facet_grid(. ~ guide, space = "free", scales = "free") +
    scale_color_gradient2(low="blue", high="red") +
    geom_point() +
    xRot()
  ggsave(out("UMAP_GSEA_", dbx, ".pdf"), w=10,h=6)
}