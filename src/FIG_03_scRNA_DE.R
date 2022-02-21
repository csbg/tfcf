source("src/00_init.R")
base.dir <- "FIG_03_scRNA_DE/"
out <- dirout(base.dir)

require(ggrepel)


# Read data ---------------------------------------------------------------
inDir <- dirout_load("SCRNA_40_01_DE_summary")

cMT <- as.matrix(read.csv(inDir("DEG_Cor.csv"), row.names = 1))
fcMT <- as.matrix(read.csv(inDir("FoldChanges.csv"), row.names=1))
umapDT <- fread(inDir("RegulatoryMap_UMAP_","all",".tsv"))
gseaDT <- fread(inDir("UMAP_GSEA.tsv"))
deDT <- fread(inDir("DEG_Statistics.tsv"))

umapDT.dim1 <- floor(max(abs(umapDT$UMAP1))) + 0.5
umapDT.dim2 <- floor(max(abs(umapDT$UMAP2))) + 0.5


# Link UMAP and DE --------------------------------------------------------
umap.log2FC.cutoff <- 3
pUMAP.de <- merge(umapDT, setNames(melt(data.table(fcMT, keep.rownames = TRUE), id.vars = "rn"), c("gene_id", "term", "estimate")), by.x="Gene", by.y="gene_id")
pUMAP.de[, guide := sub("\\..+$", "", term)]
pUMAP.de[, tissue := sub("^.+?\\.", "", term)]
pUMAP.de[, estimate_cap := pmin(umap.log2FC.cutoff, abs(estimate)) * sign(estimate)]


GOI <- c("Kmt2a", "Kmt2d", "Men1", "Rbbp4", "Setdb1", "Smarcd2", "Wdr82")

# SETUP ENDS HERE ---------------------------------------------------------





# Correlation analyses ----------------------------------------------------
pDT <- melt(data.table(cMT, keep.rownames = TRUE), id.vars = "rn")
pDT[, guide1 := sub(" .+$", "", rn)]
pDT[, tissue1 := sub("^.+? ", "", rn)]
pDT[, guide2 := sub("\\..+$", "", variable)]
pDT[, tissue2 := sub("^.+?\\.", "", variable)]
pDT <- hierarch.ordering(pDT, toOrder = "guide1", orderBy = "variable", value.var = "value", aggregate = TRUE)
pDT <- hierarch.ordering(pDT, toOrder = "guide2", orderBy = "rn", value.var = "value", aggregate = TRUE)

# Overall correlations
ggplot(pDT, aes(x=guide1, y=guide2, fill=value)) + 
  themeNF() +
  geom_tile() + 
  facet_grid(tissue2 ~ tissue1, space = "free", scales = "free") +
  scale_fill_gradient2(low="blue", high="red") +
  xlab("") + ylab("") +
  xRot()
ggsaveNF(out("CorrelationHM.pdf"), w=3,h=3)

# Normal vs cancer correlations
pDT[,variable := as.character(variable)]
pDT[make.names(rn) == make.names(variable), value := NA]
pDT.cvh <- pDT[!grepl("in.vivo", tissue1) & !grepl("in.vivo", tissue2)]
ggplot(pDT.cvh, 
       aes(
         x=paste(guide1, tissue1), 
         y=paste(guide2, tissue2), 
         fill=value)) + 
  themeNF() +
  geom_tile() + 
  facet_grid(guide2 ~ guide1, space = "free", scales = "free") +
  scale_fill_gradient2(low="blue", high="red") +
  xlab("") + ylab("") +
  theme(strip.text.x = element_text(angle=90)) +
  theme(strip.text.y = element_text(angle=0)) +
  theme(panel.spacing = unit(0, "lines")) +
  xRot()
ggsaveNF(out("Correlation_CVH_HM.pdf"), w=4,h=4)

# Rank factors by their normal vs healthy correlation
pDT.sum <- pDT.cvh[guide1 == guide2 & tissue1 == "ex.vivo_myeloid" & tissue1 != tissue2]
pDT.sum$guide <- factor(pDT.sum$guide1, levels = pDT.sum[order(-value)]$guide1)
ggplot(pDT.sum, aes(x=guide, y=value)) + 
  theme_bw() +
  geom_bar(stat="identity") + 
  ylab("Correlation normal vs cancer") +
  xRot()
ggsaveNF(out("Correlation_CVH_bars.pdf"), w=2,h=1, guides = TRUE)

# MDS of correlations
set.seed(1212)
mds.res <- data.table(cmdscale(d=as.dist(1-cMT), k=2), keep.rownames=TRUE)
mds.res <- cbind(mds.res, setNames(data.table(do.call(rbind, strsplit(mds.res$rn, " "))), c("gene", "tissue")))
ggplot(mds.res, aes(x=V1, y=V2, color=tissue, label=gene)) + 
  themeNF() +
  geom_point(size=1) + 
  ggrepel::geom_text_repel()+#color="black") + 
  xlab("MDS dimension 1") +
  ylab("MDS dimension 2")
ggsaveNF(out("CorrelationHM_MDS.pdf"), w=2, h=2)



# Scatterplot ----------------------------------------------------
pDT <- data.table(fcMT, keep.rownames = TRUE)
pDT <- melt(pDT, id.vars = "rn")
pDT[, guide := gsub("\\..+$", "", variable)]
pDT[, tissue := sub("^.+?\\.", "", variable)]
pDT <- dcast.data.table(pDT[!grepl("in.vivo", tissue)], guide + rn ~ tissue, value.var = "value")
ggplot(pDT, aes(x=ex.vivo_myeloid, y=leukemia_myeloid)) + 
  theme_bw(12) +
  stat_binhex(aes(fill=log10(..count..))) +
  facet_wrap(~guide, scales = "free")
ggsave(out("Correlation_CvH_Scatterplots.pdf"), w=20,h=20)


# Vulcano plots -----------------------------------------------------------
ggplot(deDT[guide %in% GOI], aes(x=estimate, y=pmin(30, -log10(p_value)))) + 
  theme_bw(12) +
  facet_grid(guide ~ tissue, scale="free_y") +
  stat_binhex(aes(fill=log10(..count..)))
ggsave(out("Vulcano.pdf"), w=12,h=16)



# Regulatory Map UMAP -----------------------------------------------------
ggplot(umapDT, aes(x=UMAP1, y=UMAP2, color=factor(Cluster))) + 
  geom_point(shape=1) + 
  theme_bw(12) +
  geom_label(data=umapDT[, .(UMAP1=median(UMAP1), UMAP2=median(UMAP2)), by="Cluster"], aes(label=Cluster))
ggsave(out("UMAP","_Clusters.pdf"), w=6,h=5)

#HEX
ggplot(umapDT, aes(x=UMAP1, y=UMAP2)) + 
  geom_hex(bins=100) +
  theme_bw(12) +
  geom_label(data=umapDT[, .(UMAP1=median(UMAP1), UMAP2=median(UMAP2)), by="Cluster"], aes(label=Cluster))
ggsave(out("UMAP","_Clusters_hex.pdf"), w=6,h=5)



# Plot values on UMAP -----------------------------------------------------
# tn <- length(unique(pUMAP.de$guide))
mean.umap <- function(x){mean(x, na.rm=TRUE)}

ggplot(pUMAP.de[guide %in% GOI], aes(x=UMAP1, y=UMAP2)) +
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

