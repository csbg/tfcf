source("src/00_init.R")
base.dir <- "FIG_03_scRNA_DE/"
out <- dirout(base.dir)

require(ggrepel)


# Read data ---------------------------------------------------------------
inDir <- dirout_load("FULLINT_10_02_DE")

cMT <- as.matrix(read.csv(inDir("DEG_Cor.csv"), row.names = 1))
fcMT <- as.matrix(read.csv(inDir("FoldChanges.csv"), row.names=1))
umapDT <- fread(inDir("RegulatoryMap_UMAP_","all",".tsv"))
gseaDT <- fread(inDir("UMAP_GSEA.tsv"))

umapDT.dim1 <- floor(max(abs(umapDT$UMAP1))) + 0.5
umapDT.dim2 <- floor(max(abs(umapDT$UMAP2))) + 0.5


# Link UMAP and DE --------------------------------------------------------
umap.log2FC.cutoff <- 3
pUMAP.de <- merge(umapDT, setNames(melt(data.table(fcMT, keep.rownames = TRUE), id.vars = "rn"), c("gene_id", "term", "estimate")), by.x="Gene", by.y="gene_id")
pUMAP.de[, guide := sub("\\..+$", "", term)]
pUMAP.de[, tissue := sub("^.+?\\.", "", term)]
pUMAP.de[, estimate_cap := pmin(umap.log2FC.cutoff, abs(estimate)) * sign(estimate)]




# Plot values on UMAP -----------------------------------------------------
# tn <- length(unique(pUMAP.de$guide))

mean.umap <- function(x){mean(x, na.rm=TRUE)}

gg <- c("Kmt2a", "Kmt2d", "Men1", "Rbbp4", "Setdb1", "Smarcd2", "Wdr82")
ggplot(pUMAP.de[guide %in% gg], aes(x=UMAP1, y=UMAP2)) +
  themeNF() + 
  stat_summary_hex(
    aes(z=estimate_cap),
    fun=mean.umap) +
  scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") +
  facet_grid(tissue~guide) +
  xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
  xlim(-umapDT.dim1,umapDT.dim1) + ylim(-umapDT.dim2,umapDT.dim2)
ggsaveNF(out("UMAP_Values_selection.pdf"), w=2,h=1)

ggplot(pUMAP.de, aes(x=UMAP1, y=UMAP2)) +
  themeNF() + 
  stat_summary_hex(
    aes(z=estimate_cap),
    fun=mean.umap) +
  scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") +
  facet_grid(guide ~ tissue) +
  xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
  xlim(-umapDT.dim1,umapDT.dim1) + ylim(-umapDT.dim2,umapDT.dim2)
ggsaveNF(out("UMAP_Values.pdf"), w=2,h=8)


# Plot values by UMAP cluster -----------------------------------------------------
pDT <- pUMAP.de[, mean(estimate_cap, na.rm=TRUE), by=c("Cluster", "term", "guide", "tissue")]
pDT <- hierarch.ordering(pDT, toOrder = "Cluster", orderBy = "term", value.var = "V1")
pDT <- hierarch.ordering(pDT, orderBy = "Cluster", toOrder = "guide", value.var = "V1", aggregate = TRUE)
ggplot(pDT, aes(y=factor(Cluster), x=tissue, fill=V1)) + 
  themeNF() + 
  geom_tile() +
  scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") +
  facet_grid(. ~ guide, scales = "free", space = "free", switch = "x") +
  ylab("Gene modules") +
  theme(strip.text.x = element_text(angle=90)) +
  theme(panel.spacing = unit(0.01, "cm")) +
  xlab("") +
  xRot()
ggsaveNF(out("UMAP_ClusterLogFC.pdf"), w=4,h=1)



# GSEA --------------------------------------------------------------------
fish.res <- copy(gseaDT)
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

