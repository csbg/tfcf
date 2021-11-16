source("src/00_init.R")
out <- dirout("FULLINT_10_02_DEG")

require(monocle3)
require(fgsea)
source("src/FUNC_Monocle_PLUS.R")
require(doMC)
require(foreach)
require(doParallel)

# Load full Monocle object and differential expression analysis -----------
(load(PATHS$FULLINT$Monocle))
(load(PATHS$FULLINT$DEG))
res.clean <- fread(PATHS$FULLINT$DEG.clean)
ann <- fread(PATHS$FULLINT$DEG.ann)
logFCMT <- as.matrix(read.csv(PATHS$FULLINT$DEG.logFCMT))
cds <- monocle.obj[,ann$rn]
rm(list="monocle.obj")

load(PATHS$CHIP$Targets)
load(PATHS$RESOURCES$Enrichr.mouse)
enr.terms$ChromatinFactors <- chip.targets

umap <- fread(PATHS$FULLINT$DEG.UMAP)




# Plots some genes (only in vitro) -------------------------------------------
res.ex <- res.clean[interaction == FALSE & tissue == "in vitro" & q_value < 0.05 & estimate > 0]
res.ex[,.N, by="guide"]

genex <- "Kmt2d"

# Plot analysis for one guide
gg <- res.ex[guide == genex]$gene_id
annx <- ann[GuideDE %in% c("NTC", genex) & tissueDE == "in vitro"]

# Prepare exprsesion matrix
eMT <- NF_TPM_Matrix(cds, gg)
eMT <- eMT[,annx$rn]
eMT <- SCRNA.TPXToLog(eMT)
eMT <- t(scale(t(as.matrix(eMT))))
eMT[eMT > 5] <- 5
eMT[eMT < -5] <- -5

# column annotation
annC <- data.frame(
  row.names = annx$rn,
  guide=annx$GuideDE,
  cluster=annx$ClusterDE
)

# Heatmap
cleanDev(); pdf(out("HM.pdf"),w=10,h=10)
pheatmap(
  cluster_cols = FALSE,
  annotation_col = annC,
  eMT[,annx[order(GuideDE, ClusterDE)]$rn], 
  show_rownames = FALSE, 
  show_colnames = FALSE
  )
dev.off()

# Violin plots
str(colData(cds))
if(length(gg) > 40) gg <- sample(gg, 40)
pObj <- cds[gg,annx$rn]
pObj$cluster <- annx$ClusterDE
pObj$group <- paste(pObj$guide, pObj$tissue, paste0("c", pObj$cluster))
plot_genes_violin(
  pObj,
  log_scale = FALSE,
  label_by_short_name=FALSE,
  group_cells_by="group", 
  ncol = 10) + 
  scale_y_continuous(trans = "log1p") +
  xRot()
ggsave(out("Violins.pdf"), w=30,h=30)



# Recalculate interaction coefficients -------------------------------------------
res.invitro <- res.clean[interaction == FALSE & tissue == "in vitro"]
res.interaction.leukemia <- res.clean[interaction == TRUE & tissue == "leukemia"]
res.leukemia <- fread(dirout_load("FULLINT_10_01_BasicAnalysis_leukemia")("DEG_Results.tsv"))
gg <- intersect(intersect(res.invitro$guide, res.interaction.leukemia$guide), res.leukemia$guide)

pDT <- merge(res.invitro, res.leukemia, by=c("guide", "gene_id"), suffixes=c("_vitro", "_leuk"))
pDT <- merge(pDT, res.interaction.leukemia, by=c("guide", "gene_id"))
ggplot(pDT, aes(x=estimate_leuk, y=estimate_vitro + estimate)) + 
  theme_bw(12) +
  facet_wrap(~guide) +
  geom_hex() +
  xlab("Measured estimate (leukemia)") +
  ylab("Calculated estimate (vitro + interaction)")
ggsave(out("InteractionCoefficients.pdf"), w=10,h=10)

ggplot(pDT[,cor(estimate_leuk, estimate_vitro + estimate), by="guide"],
       aes(x=guide, y=V1)) +
  theme_bw(12) +
  geom_bar(stat='identity') +
  xRot()
ggsave(out("InteractionCoefficients_Corr.pdf"), w=6,h=4)




# Signal in gene clustes --------------------------------------------------
# gg <- umap[Cluster ==1]$Gene
# umap.logFC <- sapply(with(umap, split(Gene, Cluster)), function(gg){colMeans(logFCMT[gg,])})
# pheatmap(umap.logFC)


# Enrichment in gene clusters ---------------------------------------------
fish.genes <- with(umap, split(Gene, Cluster))
fish.res <- fisher.test.enrichment(geneSets = enr.terms, gene.list = fish.genes, cores=10)
fish.res[, padj := p.adjust(pval, method="BH")]
write.tsv(fish.res, out("UMAP_GSEA.tsv"))

# Plots
fish.res[, db := database]
fish.res[, pathway := geneset]
fish.res[, log2OR := log2(pmin(10, pmax(1/10, oddsRatio)))]
fish.res[db == "ChromatinFactors", pathway := gsub("^(.+?)_(.+)$", "\\2 \\1", pathway)]
dbx <- fish.res$db[1]
for(dbx in unique(fgsea$db)){
  pDT <- fish.res[db == dbx]
  pDT <- pDT[pathway %in% unique(pDT[padj < 0.05][order(-abs(log2OR))]$pathway)[1:30]]
  ggplot(pDT, aes(y=factor(as.numeric(list)), x=pathway, size=pmin(5, -log10(padj)), color=log2OR)) + 
    theme_bw(12) +
    #facet_grid(. ~ guide, space = "free", scales = "free") +
    scale_color_gradient2(low="blue", high="red") +
    geom_point() +
    xRot()
  ggsave(out("UMAP_GSEA_", dbx, ".pdf"), w=10,h=6)
}



# Enrichments -------------------------------------------------------------
fgsea.file <- out("GSEA.tsv")
if(file.exists(fgsea.file)){
  fgsea <- fread(fgsea.file)
} else {
  fgsea.hitlists <- setNames(lapply(colnames(logFCMT), function(cx){
    setNames(logFCMT[,cx], row.names(logFCMT))
  }), colnames(logFCMT))
  dbx <- names(enr.terms)[1]
  hl <- names(fgsea.hitlists)[1]
  registerDoMC(cores=8)
  fgsea <- foreach(hl = names(fgsea.hitlists)) %dopar% {
    ret <- data.table()
    for(dbx in names(enr.terms)){
      fgsea.res <- fgsea(
        pathways=enr.terms[[dbx]],
        stats=fgsea.hitlists[[hl]],
        BPPARAM=BiocParallel::SnowParam(1, "SOCK")
      )
      fgsea.res$leadingEdge <- sapply(fgsea.res$leadingEdge, function(gx) paste(gx, collapse = ","))
      fgsea.res$hitlist <- hl
      fgsea.res$db <- dbx
      ret <- rbind(ret, fgsea.res)
    }
    return(ret)
  }
  fgsea <- do.call(rbind, fgsea)
  write.tsv(fgsea, fgsea.file)
}

# Plots
fgsea[db == "ChromatinFactors", pathway := gsub("^(.+?)_(.+)$", "\\2 \\1", pathway)]
dbx <- fgsea$db[1]
for(dbx in unique(fgsea$db)){
  pDT <- fgsea[db == dbx]
  pDT <- pDT[pathway %in% unique(pDT[padj < 0.05][order(-abs(NES))]$pathway)[1:30]]
  pDT[, guide := gsub("\\..+$", "", hitlist)]
  ggplot(pDT, aes(x=hitlist, y=pathway, size=pmin(5, -log10(padj)), color=NES)) + 
    theme_bw(12) +
    facet_grid(. ~ guide, space = "free", scales = "free") +
    scale_color_gradient2(low="blue", high="red") +
    geom_point() +
    xRot()
  ggsave(out("GSEA_", dbx, ".pdf"), w=15,h=10)
}