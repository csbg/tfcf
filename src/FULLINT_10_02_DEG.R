source("src/00_init.R")

base.dir <- "FULLINT_10_02_DE/"
out <- dirout(base.dir)


require(ggrepel)


# FUNCTIONS ---------------------------------------------------------------
ds <- function(path){load(path); return(monocle.obj)}


# EnrichR -----------------------------------------------------------------
(load(PATHS$RESOURCES$Enrichr.mouse))


# Folders -----------------------------------------------------------------
list.files(dirout_load("")(""))
inDir.funcs <- list(
  "in.vivo"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vivo"),
  "in.vitro"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vitro"),
  "leukemia"=dirout_load("FULLINT_10_01_BasicAnalysis_leukemia")
)
inDir.full <- dirout_load("FULLINT_10_01_BasicAnalysis_combined")



# Load data ---------------------------------------------------------------
ff <- lapply(inDir.funcs, function(x) x("DEG_Results_all.tsv"))
ff <- lapply(ff, fread)
DE.RES <- do.call(rbind, ff)
write.tsv(DE.RES[, -c("se", "convergence", "interaction", "estimate_raw"), with=F], out("DEG_Statistics.tsv"))

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


# Compare logFC to UMAP ---------------------------------------------------
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



# out <- dirout("FULLINT_10_02_DEG")
# 
# require(monocle3)
# require(fgsea)
# source("src/FUNC_Monocle_PLUS.R")
# require(doMC)
# require(foreach)
# require(doParallel)
# 
# # Load full Monocle object and differential expression analysis -----------
# (load(PATHS$FULLINT$Monocle))
# (load(PATHS$FULLINT$DEG))
# res.clean <- fread(PATHS$FULLINT$DEG.clean)
# res.clean.t <- list(
#   combined = copy(res.clean),
#   leukemia = fread(PATHS$FULLINT$DEG.clean.leukemia),
#   invitro = fread(PATHS$FULLINT$DEG.clean.invitro),
#   invivo = fread(PATHS$FULLINT$DEG.clean.invivo)
# )
# res.clean.t <- lapply(res.clean.t, function(dt){dt[,direction := ifelse(estimate > 0, "up", "down")]; dt})
# ann <- fread(PATHS$FULLINT$DEG.ann)
# logFCMT <- as.matrix(read.csv(PATHS$FULLINT$DEG.logFCMT))
# cds <- monocle.obj[,ann$rn]
# rm(list="monocle.obj")
# 
# load(PATHS$CHIP$Targets)
# load(PATHS$RESOURCES$Enrichr.mouse)
# enr.terms$ChromatinFactors <- chip.targets
# 
# umap <- fread(PATHS$FULLINT$DEG.UMAP)
# 
# 
# 
# # Plots some genes (only in vitro) -------------------------------------------
# res.ex <- res.clean[interaction == FALSE & tissue == "in vitro" & q_value < 0.05 & estimate > 1]
# res.ex[,.N, by="guide"]
# 
# genex <- "Kmt2d"
# 
# # Plot analysis for one guide
# gg <- res.ex[guide == genex]$gene_id
# annx <- ann[GuideDE %in% c("NTC", genex) & tissueDE == "in vitro"]
# 
# # Prepare exprsesion matrix
# eMT <- NF_TPM_Matrix(cds, gg)
# eMT <- eMT[,annx$rn]
# eMT <- SCRNA.TPXToLog(eMT)
# eMT <- t(scale(t(as.matrix(eMT))))
# eMT[eMT > 5] <- 5
# eMT[eMT < -5] <- -5
# 
# # column annotation
# annC <- data.frame(
#   row.names = annx$rn,
#   guide=annx$GuideDE,
#   cluster=annx$ClusterDE
# )
# 
# # Heatmap
# cleanDev(); pdf(out("Example_HM.pdf"),w=10,h=10)
# pheatmap(
#   cluster_cols = FALSE,
#   annotation_col = annC,
#   eMT[,annx[order(GuideDE, ClusterDE)]$rn], 
#   show_rownames = FALSE, 
#   show_colnames = FALSE
#   )
# dev.off()
# 
# # Violin plots
# str(colData(cds))
# if(length(gg) > 40) gg <- sample(gg, 40)
# pObj <- cds[gg,annx$rn]
# pObj$cluster <- annx$ClusterDE
# pObj$group <- paste(pObj$guide, pObj$tissue, paste0("c", pObj$cluster))
# plot_genes_violin(
#   pObj,
#   log_scale = FALSE,
#   label_by_short_name=FALSE,
#   group_cells_by="group", 
#   ncol = 10) + 
#   scale_y_continuous(trans = "log1p") +
#   xRot()
# ggsave(out("Example_Violins.pdf"), w=30,h=30)
# 
# 
# 
# # Recalculate interaction coefficients -------------------------------------------
# res.invitro <- res.clean[interaction == FALSE & tissue == "in vitro"]
# res.interaction.leukemia <- res.clean[interaction == TRUE & tissue == "leukemia"]
# res.leukemia <- fread(dirout_load("FULLINT_10_01_BasicAnalysis_leukemia")("DEG_Results_all.tsv"))
# gg <- intersect(intersect(res.invitro$guide, res.interaction.leukemia$guide), res.leukemia$guide)
# 
# pDT <- merge(res.invitro, res.leukemia, by=c("guide", "gene_id"), suffixes=c("_vitro", "_leuk"))
# pDT <- merge(pDT, res.interaction.leukemia, by=c("guide", "gene_id"))
# ggplot(pDT, aes(x=estimate_leuk, y=estimate_vitro + estimate)) + 
#   theme_bw(12) +
#   facet_wrap(~guide) +
#   geom_hex() +
#   xlab("Measured estimate (leukemia)") +
#   ylab("Calculated estimate (vitro + interaction)")
# ggsave(out("InteractionCoefficients.pdf"), w=10,h=10)
# 
# ggplot(pDT[,cor(estimate_leuk, estimate_vitro + estimate), by="guide"],
#        aes(x=guide, y=V1)) +
#   theme_bw(12) +
#   geom_bar(stat='identity') +
#   xRot()
# ggsave(out("InteractionCoefficients_Corr.pdf"), w=6,h=4)
# 
# 
# 
# 
# # Number of genes ---------------------------------------------------------
# cnts <- lapply(res.clean.t, function(x) x[q_value < 0.05 & abs(estimate) > 1][,.N, by=c("guide", "tissue", "interaction", "direction")][order(N)])
# pDT <- do.call(rbind, cnts)
# pDT[direction == "down", N := -N]
# ggplot(pDT, aes(x=tissue, y=N, fill=interaction)) +
#   geom_bar(stat="identity", position = "dodge") +
#   facet_grid(. ~ guide) +
#   xRot()
# 
# 
# 
# # Plot top genes ----------------------------------------------------------
# pDT <- do.call(rbind, res.clean.t)
# str(gg <- unique(pDT[q_value < 0.01][order(-abs(estimate))][,head(.SD, 5), by=c("guide", "tissue")]$gene_id))
# pDT <- pDT[gene_id %in% gg]
# pDT[, term := paste(guide, tissue, interaction)]
# pDT <- hierarch.ordering(pDT, toOrder="gene_id", orderBy = "term", value.var = "estimate", aggregate = TRUE)
# pDT[, int := ifelse(interaction == TRUE, "interaction", "basic")]
# pDT$UMAP.Cluster <- umap[match(pDT$gene_id, Gene)]$Cluster
# ggplot(pDT, aes(x=paste(tissue), y=gene_id, size=pmin(5, -log10(q_value)), color=sign(estimate) * pmin(5,abs(estimate)))) +
#   scale_size_continuous(name="padj", range=c(0,5)) +
#   scale_color_gradient2(name="delta", low="blue", high="red") +
#   facet_grid(UMAP.Cluster ~ guide + int, space = "free", scales = "free") + 
#   theme_bw(12) +
#   geom_point() +
#   xRot()
# ggsave(out("TopGenes.pdf"), w=25,h=30)
# 
# # Signal in gene clustes --------------------------------------------------
# # gg <- umap[Cluster ==1]$Gene
# # umap.logFC <- sapply(with(umap, split(Gene, Cluster)), function(gg){colMeans(logFCMT[gg,])})
# # pheatmap(umap.logFC)
# 
# 
# # Enrichment in gene clusters ---------------------------------------------
# fish.genes <- with(umap, split(Gene, Cluster))
# fish.res <- fisher.test.enrichment(
#   geneSets = enr.terms, 
#   gene.list = fish.genes, 
#   cores=10,
#   bg=unique(res.clean$gene_id))
# fish.res[, padj := p.adjust(pval, method="BH")]
# write.tsv(fish.res, out("UMAP_GSEA.tsv"))
# 
# # Plots
# fish.res[, db := database]
# fish.res[, pathway := geneset]
# fish.res[, log2OR := log2(pmin(10, pmax(1/10, oddsRatio)))]
# fish.res[db == "ChromatinFactors", pathway := gsub("^(.+?)_(.+)$", "\\2 \\1", pathway)]
# dbx <- fish.res$db[1]
# for(dbx in unique(fish.res$db)){
#   pDT <- fish.res[db == dbx]
#   pwx <- unique(pDT[padj < 0.1][order(-log2OR)][,head(.SD,n=3), by="list"]$pathway)
#   pDT <- pDT[pathway %in% pwx]
#   ggplot(pDT, aes(y=factor(as.numeric(list)), x=pathway, size=pmin(5, -log10(padj)), color=log2OR)) + 
#     theme_bw(12) +
#     #facet_grid(. ~ guide, space = "free", scales = "free") +
#     scale_color_gradient2(low="blue", high="red") +
#     geom_point() +
#     xRot()
#   ggsave(out("UMAP_GSEA_", dbx, ".pdf"), w=10,h=6)
# }
# 
# 
# 
# # Enrichments -------------------------------------------------------------
# fgsea.file <- out("GSEA.tsv")
# if(file.exists(fgsea.file)){
#   fgsea <- fread(fgsea.file)
# } else {
#   fgsea.hitlists <- setNames(lapply(colnames(logFCMT), function(cx){
#     setNames(logFCMT[,cx], row.names(logFCMT))
#   }), colnames(logFCMT))
#   dbx <- names(enr.terms)[1]
#   hl <- names(fgsea.hitlists)[1]
#   registerDoMC(cores=8)
#   fgsea <- foreach(hl = names(fgsea.hitlists)) %dopar% {
#     ret <- data.table()
#     for(dbx in names(enr.terms)){
#       fgsea.res <- fgsea(
#         pathways=enr.terms[[dbx]],
#         stats=fgsea.hitlists[[hl]],
#         BPPARAM=BiocParallel::SnowParam(1, "SOCK")
#       )
#       fgsea.res$leadingEdge <- sapply(fgsea.res$leadingEdge, function(gx) paste(gx, collapse = ","))
#       fgsea.res$hitlist <- hl
#       fgsea.res$db <- dbx
#       ret <- rbind(ret, fgsea.res)
#     }
#     return(ret)
#   }
#   fgsea <- do.call(rbind, fgsea)
#   write.tsv(fgsea, fgsea.file)
# }
# 
# # Check against other analysis
# # (load(dirout_load("FULLINT_10_01_BasicAnalysis_combined")("DEG_ChIP_GSEA.RData")))
# # fgsea.hitlists
# # gplots::venn(list(enr.terms$ChromatinFactors$Mye_Kmt2d, chip.targets$Mye_Kmt2d))
# # x <- resGuides.I[id == "Kmt2d leukemia"]
# # plot(x$estimate, fgsea.hitlists$Kmt2d.leukemia[x$gene_id])
# # fgsea[hitlist == "Kmt2d.leukemia" & pathway == "Mye_Kmt2d"][,1:7]
# # chip.gsea.res[list == "Kmt2d leukemia" & pathway == "Mye_Kmt2d"]
# 
# # Plots
# fgsea[db == "ChromatinFactors", pathway := gsub("^(.+?)_(.+)$", "\\2 \\1", pathway)]
# dbx <- fgsea$db[1]
# for(dbx in unique(fgsea$db)){
#   pDT <- fgsea[db == dbx]
#   pDT <- pDT[pathway %in% unique(pDT[padj < 0.05][order(-abs(NES))]$pathway)[1:30]]
#   pDT[, guide := gsub("\\..+$", "", hitlist)]
#   ggplot(pDT, aes(x=hitlist, y=pathway, size=pmin(5, -log10(padj)), color=NES)) + 
#     theme_bw(12) +
#     facet_grid(. ~ guide, space = "free", scales = "free") +
#     scale_color_gradient2(low="blue", high="red") +
#     geom_point() +
#     xRot()
#   ggsave(out("GSEA_", dbx, ".pdf"), w=15,h=10)
# }
# 
# 
# 
# 
# # Enrichment based on significant genes in each guide? -----------------------------------
# 
# # Plot Values on UMAP (logFC or % of genes) -------------------------------
# 
# # # Plot estimates on UMAP
# # pUMAP.de <- merge(umap, setNames(melt(data.table(umapMT, keep.rownames = TRUE), id.vars = "rn"), c("gene_id", "term", "estimate")), by.x="Gene", by.y="gene_id")
# # summary.function <- function(x){ret <- mean(x);return(min(5, abs(ret)) * sign(ret))}
# # dim.umap1 <- floor(max(abs(pUMAP.de$UMAP1))) + 0.5
# # dim.umap2 <- floor(max(abs(pUMAP.de$UMAP2))) + 0.5
# # pUMAP.de[, guide := gsub(" .+", "", term)]
# # pUMAP.de[, tissue := gsub(".+? ", "", term)]
# # tn <- length(unique(pUMAP.de$guide))
# # ggplot(pUMAP.de, aes(x=UMAP1, y=UMAP2)) +
# #   stat_summary_hex(
# #     aes(z=estimate),
# #     fun=summary.function) +
# #   scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") +
# #   facet_grid(guide~tissue) + theme_bw(12) +
# #   xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
# #   xlim(-dim.umap1,dim.umap1) + ylim(-dim.umap2,dim.umap2)
# # ggsave(out("RegulatoryMap_UMAP_",umap.type.name,"_Values.pdf"), w=2*3+2,h=tn * 2 + 1)
# # 
# # # values by cluster
# # pDT <- pUMAP.de[, mean(estimate), by=c("Cluster", "term")]
# # pDT <- hierarch.ordering(pDT, toOrder = "Cluster", orderBy = "term", value.var = "V1")
# # pDT <- hierarch.ordering(pDT, orderBy = "Cluster", toOrder = "term", value.var = "V1")
# # ggplot(pDT, aes(x=factor(Cluster), y=term, fill=V1)) + 
# #   theme_bw(12) + 
# #   geom_tile() +
# #   scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") +
# #   xlab("Gene modules (Gene-UMAP Clusters)")
# # ggsave(out("RegulatoryMap_UMAP_",umap.type.name,"_ClusterValues.pdf"), w=10,h=tn * 0.2 + 1)