source("src/00_init.R")
base.dir <- "FIG_03_scRNA_DE/"
out <- dirout(base.dir)

require(ggrepel)


# Read data ---------------------------------------------------------------
inDir <- dirout_load("SCRNA_40_01_DE_summary")

annList <- lapply(PATHS$SCRNA$MONOCLE.NAMES, function(tx){
  ann <- fread(dirout_load(paste0("SCRNA_20_Summary/", tx, "_monocle.singleR"))("Annotation.tsv"))
  ann[, gene := gsub("_.+", "", guide)]
  ann
})
annList <- rbindlist(annList, fill=TRUE)

cMT <- as.matrix(read.csv(inDir("DEG_Cor.csv"), row.names = 1))
fcMT <- as.matrix(read.csv(inDir("FoldChanges.csv"), row.names=1))
umapDT <- fread(inDir("RegulatoryMap_UMAP_","all",".tsv"))
umapDT.dim1 <- floor(max(abs(umapDT$UMAP1))) + 0.5
umapDT.dim2 <- floor(max(abs(umapDT$UMAP2))) + 0.5
gseaDT <- fread(inDir("UMAP_GSEA.tsv"))
deDT <- fread(inDir("DEG_Statistics.tsv"))

atacDT.lsk <- fread(paste0(gsub("Data", "CollaborationData", PATHS$LOCATIONS$DATA), "LSK_narrowPeak_consensusPeaks.boolean.annotatePeaks.extended.fc.txt"))
atacDT <- data.table::melt(
  atacDT.lsk,
  id.vars = c("chr", "start", "end", "Distance to TSS", "Gene Name"),
  measure.vars = grep("log2FC$", colnames(atacDT.mye), value = TRUE)
  )
atacDT[, variable := gsub("LSK-(.+?)_.+", "\\L\\1", variable, perl=TRUE)]
atacDT[, variable := gsub("^(.)", "\\U\\1", variable, perl=TRUE)]


# Link UMAP and DE --------------------------------------------------------
umap.log2FC.cutoff <- 3
pUMAP.de <- merge(umapDT, setNames(data.table::melt(data.table(fcMT, keep.rownames = TRUE), id.vars = "rn"), c("gene_id", "term", "estimate")), by.x="Gene", by.y="gene_id")
pUMAP.de[, guide := sub("\\..+$", "", term)]
pUMAP.de[, tissue := sub("^.+?\\.", "", term)]
pUMAP.de[, estimate_cap := pmin(umap.log2FC.cutoff, abs(estimate)) * sign(estimate)]


GOI <- c("Kmt2a", "Kmt2d", "Men1", "Rbbp4", "Setdb1", "Smarcd2", "Wdr82", "Cbx3", "Hdac1", "Crebbp")

# SETUP ENDS HERE ---------------------------------------------------------




# Plot specific genes -----------------------------------------------------
tx <- deDT$tissue[1]
for(tx in unique(deDT$tissue)){
  pDT <- deDT[tissue == tx]
  pDT.l <- list(
    markers = pDT[grep("Mpo|Elane|Ctsg|S100a9|S100a8|Mtf2|Car1|Car2|Epor|Pklr|Pf4|Cpa3|Hba|Hbb|Tfrc|Cpox|Mt2|Trem2|Mefv|F13a1|Ly6c2|Fcgr3|Cd34|Kit|Itgam",gene_id),],
    #tfs = pDT[grep("Cited2|Satb1|Hoxa9|Hoxa7|Hoxb4|Fli1|Erg|Spi1|Cebpa|Cebpb|Cebpe|Mef2c|Klf4|Irf8|Gfi1|Gata1|Gata2|Klf1|Klf9|Klf3|Nfia|Nfix|Nfe2|Sox6|Mafk|Gfi1b|Tal1|Cited4|Lmo1|Lmo4",gene_id),]
    tfs = pDT[grep("Spi1|Cebpa|Cebpb|Cebpe|Klf4|Gfi1|Irf8|Mef2c|Gata2|Gata1|Klf1|Zfpm1|Mafk|Sox6|Klf9|Lmo4|Nfia|Nfix|Nfe2|Nfe2l2|Tal1|Cited4",gene_id),]
    )
  for(lnam in names(pDT.l)){
    xDT <- pDT.l[[lnam]]
    xDT <- hierarch.ordering(xDT, "guide", "gene_id", value.var = "estimate")
    xDT <- hierarch.ordering(xDT, "gene_id", "guide", value.var = "estimate")
    ggplot(xDT, aes(x=guide, y=gene_id, color=estimate, size=pmin(5, -log10(q_value)))) + 
      themeNF(rotate=TRUE) +
      geom_point() +
      scale_color_gradient2(name="log2FC", low="blue", high="red") +
      scale_size_continuous(name="-log10(padj)", range = c(0,5)) +
      ggtitle(tx) +
      ylab("Gene expression") + xlab("CRISPR targets")
    ggsaveNF(out("GeneDE_", tx, "_", lnam, ".pdf"),w=3,h=2, guides = TRUE)
  }
}




# ATAC-seq comparison -----------------------------------------------------
atacDT[order(abs(get("Distance to TSS")), decreasing = FALSE)][, head(.SD, n=1), by=c("Gene Name")]
pDT <- merge(
  deDT[tissue == "ex.vivo_myeloid"][, c("gene_id", "guide", "estimate"), with=FALSE],
  atacDT[order(abs(get("Distance to TSS")), decreasing = FALSE)][, head(.SD, n=1), by=c("Gene Name", "variable")][,c("Gene Name", "variable", "value")],
  by.x=c("gene_id", "guide"), by.y=c("Gene Name", "variable"))
pCor <- pDT[, cor(estimate, value, use="pairwise.complete.obs", method="spearman"), by=c("guide")]
ggplot(pDT, aes(x=estimate, y=value)) + 
  themeNF() +
  stat_binhex(aes(fill=log10(..count..))) + 
  #ggtitle(paste(round(pt$estimate, 3))) +
  geom_text(data=pCor, aes(label=paste0("R=", round(V1,3))), x=-Inf,y=Inf, hjust=0, vjust=1)+
  facet_grid(. ~ guide) +
  xlab("RNA-seq logFC KO vs WT") +
  ylab("ATAC-seq logFC KO vs WT")
ggsaveNF(out("ATAC_correlation.pdf"), w=4,h=1)



# Correlation in BAF ------------------------------------------------------
pDT <- data.table::melt(data.table(cMT, keep.rownames = TRUE), id.vars = "rn")
pDT <- pDT[grepl("ex.vivo", rn) | grepl("in.vivo", rn)]
pDT <- pDT[grepl("ex.vivo", variable) | grepl("in.vivo", variable)]
pDT <- pDT[grepl("_myeloid", rn) & grepl("_myeloid", variable)]
pDT <- pDT[grepl("Brd9|Kmt2d|Kmt2a|Smarcb1", rn) & grepl("Brd9|Kmt2d|Kmt2a|Smarcb1", variable)]
pDT[, rn := gsub("_myeloid", "", rn)]
pDT[, variable := gsub("_myeloid", "", variable)]
pDT[, variable := gsub("\\.", " ", variable)]
pDT[, rn := gsub("\\.", " ", rn)]
pDT <- pDT[rn != variable]
pDT <- hierarch.ordering(pDT, toOrder = "rn", orderBy = "variable", value.var = "value")
pDT <- hierarch.ordering(pDT, toOrder = "variable", orderBy = "rn", value.var = "value")
ggplot(pDT, aes(x=rn, y=variable, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="blue", high="red") +
  themeNF(rotate=TRUE)
ggsaveNF(out("Phenocopy_BAF.pdf"), w=1,h=1)


# Number of up- and downregualted genes ----------------------------------------------------
pDT.up_down <- copy(deDT[tissue == "leukemia_myeloid"])
pDT.up_down[, direction := ifelse(estimate > 0, "up", "down")]
pDT.up_down <- pDT.up_down[q_value < 0.05][abs(estimate) > 0.5][,length(unique(gene_id)), by=c("guide", "direction")]
#pDT[direction == "down", V1 := -V1]
ggplot(pDT.up_down, aes(x=guide, y=V1, fill=direction)) +
  themeNF() +
  geom_bar(stat="identity", position="dodge") +
  scale_y_log10() +
  coord_flip()
ggsaveNF(out("NumberOfGenes.pdf"), w=1,h=2)


# Correlation Cancer vs Healthy ----------------------------------------------------
pDT <- data.table::melt(data.table(cMT, keep.rownames = TRUE), id.vars = "rn")
pDT[, guide1 := sub(" .+$", "", rn)]
pDT[, tissue1 := sub("^.+? ", "", rn)]
pDT[, guide2 := sub("\\..+$", "", variable)]
pDT[, tissue2 := sub("^.+?\\.", "", variable)]
pDT <- hierarch.ordering(pDT, toOrder = "guide1", orderBy = "variable", value.var = "value", aggregate = TRUE)
pDT <- hierarch.ordering(pDT, toOrder = "guide2", orderBy = "rn", value.var = "value", aggregate = TRUE)

# Overall correlations
# ggplot(pDT, aes(x=guide1, y=guide2, fill=value)) + 
#   themeNF() +
#   geom_tile() + 
#   facet_grid(tissue2 ~ tissue1, space = "free", scales = "free") +
#   scale_fill_gradient2(low="blue", high="red") +
#   xlab("") + ylab("") +
#   xRot()
# ggsaveNF(out("CorrelationHM.pdf"), w=3,h=3)

# Cancer vs Healthy correlations
pDT[,variable := as.character(variable)]
pDT[make.names(rn) == make.names(variable), value := NA]
pDT.cvh <- pDT[!grepl("in.vivo", tissue1) & !grepl("in.vivo", tissue2)]
double.guides <- pDT.cvh[guide1 == guide2][!is.na(value)]$guide1
pDT.cvh <- pDT.cvh[guide1 %in% double.guides]
pDT.cvh <- pDT.cvh[guide2 %in% double.guides]
pDT.cvh[,tissue1 := gsub("_.+","", tissue1)]
pDT.cvh[,tissue2 := gsub("_.+","", tissue2)]

# Cancer vs Healthy same factor Barplot
pDT.sum <- pDT.cvh[guide1 == guide2 & tissue1 == "ex.vivo" & tissue1 != tissue2]
pDT.sum$guide <- factor(pDT.sum$guide1, levels = pDT.sum[order(value)]$guide1)
ggplot(pDT.sum, aes(x=guide, y=value)) + 
  themeNF() +
  geom_bar(stat="identity") + 
  ylab("Correlation normal vs cancer") +
  coord_flip() +
  ylab("")
ggsaveNF(out("Correlation_CVH_bars.pdf"), w=1,h=2, guides = TRUE)



pDT.sum.control <- merge(
  pDT.sum,
  dcast.data.table(annList[mixscape_class.global == "KO", .N, by=c("gene", "tissue")][tissue != "in.vivo"], gene ~ tissue, value.var="N"),
  by.x="guide1", by.y="gene")
ggplot(pDT.sum.control,
       aes(x=ex.vivo, y=leukemia, label=guide1)) +
  geom_point(aes(color=value), size=3) +
  geom_text_repel() +
  themeNF(rotate=TRUE) +
  scale_color_gradient2(low="blue", high="red") +
  scale_y_log10(limits=c(10,1e4)) +
  scale_x_log10(limits=c(10,1e4)) +
  geom_abline()
ggsaveNF(out("Correlation_CVH_bars_GuideCellNumbers.pdf"), w=2,h=2)


# Cancer vs Healthy Heatmap
# ggplot(pDT.cvh,
#        aes(
#          x=paste(tissue1),
#          y=paste(tissue2),
#          fill=value)) +
#   themeNF() +
#   geom_tile() +
#   facet_grid(guide2 ~ guide1, space = "free", scales = "free", switch = "both") +
#   scale_fill_gradient2(low="blue", high="red") +
#   xlab("") + ylab("") +
#   theme(strip.text.x = element_text(angle=90)) +
#   theme(strip.text.y.left = element_text(angle=0)) +
#   theme(panel.spacing = unit(0, "lines")) +
#   xRot()
# ggsaveNF(out("Correlation_CVH_HM.pdf"), w=3,h=3)


# Cancer vs Healthy Heatmap of changes phenocopies
pDT.cvh.diff <- dcast.data.table(pDT.cvh[tissue1 == tissue2], guide1 + guide2 ~ tissue1, value.var = "value")
pDT.cvh.diff[, diff := leukemia - ex.vivo]
pDT.cvh.diff <- hierarch.ordering(pDT.cvh.diff, "guide1", "guide2", "diff")
pDT.cvh.diff <- hierarch.ordering(pDT.cvh.diff, "guide2", "guide1", "diff")
pDT.cvh.diff[, diff.grp := paste("sign", ifelse(abs(diff) > 0.1, sign(diff), 0))]
ggplot(pDT.cvh.diff, aes(x=guide1, y=guide2, fill=ex.vivo, size=abs(diff), color=diff.grp)) + 
  geom_point(shape=21) +
  scale_color_manual(name="Change\ndirection", values=c("blue", "grey", "red", "grey")) +
  scale_fill_gradient2(name="Correlation\nex vivo", low="blue", high="red") +
  scale_size_continuous(name="Change\nmagnitude") +
  themeNF(rotate = TRUE) +
  ggtitle("Change in correlation leukemia vs ex vivo") +
  xlab("CRISPR targets") + ylab("CRISPR targets")
ggsaveNF(out("Correlation_CVH_HM_Details_HM.pdf"), w=2,h=2)


# Cancer vs Healthy changes phenocopies control for number of cells
xDT <- dcast.data.table(annList[mixscape_class.global == "KO", .N, by=c("gene", "tissue")][tissue != "in.vivo"], gene ~ tissue, value.var="N")
xDT2 <- pDT.cvh.diff[,.(meandiff = mean(diff, na.rm=TRUE)), by=c("guide1")]
xDT2[, gene := as.character(guide1)]
xDT <- merge.data.table(xDT, xDT2, by="gene")
ggplot(xDT,
       aes(x=ex.vivo, y=leukemia, label=gene)) +
  geom_point(aes(color=meandiff), size=3) +
  geom_text_repel() +
  themeNF(rotate=TRUE) +
  scale_color_gradient2(low="blue", high="red") +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline()
ggsaveNF(out("Correlation_CVH_HM_Details_HM_GuideCellNumbers.pdf"), w=2,h=2)

# Cancer vs Healthy of changes phenocopies Heatmap examples
pDT.cvh.diff.2 <- copy(pDT.cvh.diff[sign(ex.vivo) != sign(leukemia)][order(-abs(diff))][1:50])
pDT.cvh.diff.2 <- pDT.cvh.diff.2[t(apply(as.matrix(pDT.cvh.diff.2[,c("guide1", "guide2")]), 1, order))[,1] == 1]
pDT.cvh.diff.2[, pair := paste0(guide1, "-", guide2)]
pDT.cvh.diff.2$pair <- factor(pDT.cvh.diff.2$pair, levels=pDT.cvh.diff.2[order(diff)]$pair)
pDT.cvh.diff.2 <- melt(pDT.cvh.diff.2, id.vars = c("pair"), measure.vars=c("ex.vivo", "leukemia", "diff"))
ggplot(pDT.cvh.diff.2, aes(x=variable, y=pair, fill=value)) + 
  themeNF(rotate=TRUE) +
  geom_tile() +
  scale_fill_gradient2(low="blue", high="red") +
  xlab("Correlation")
ggsaveNF(out("Correlation_CVH_HM_Details_HM_Examples.pdf"), w=1,h=2)

# Cancer vs Healthy of changes phenocopies Scatterplot
pDT.cvh.diff.scatter <- pDT.cvh.diff[t(apply(as.matrix(pDT.cvh.diff[,c("guide1", "guide2")]), 1, order))[,1] == 1]
gg <- unique(rbind(
  pDT.cvh.diff.scatter[ex.vivo > 0.4 | ex.vivo < -0.15],
  pDT.cvh.diff.scatter[sign(leukemia) != sign(ex.vivo)][order(abs(diff), decreasing = TRUE)][1:20]
))
# ggplot(pDT.cvh.diff, aes(x=ex.vivo, y=diff)) + 
#   geom_point() +
#   geom_text_repel(data=gg, aes(label=paste0(guide1, "-", guide2))) +
#   themeNF() +
#   xlab("Correlation in ex vivo") +
#   ylab("Difference in correlation\nleukemia vs ex vivo")
# ggsaveNF(out("Correlation_CVH_HM_Details.pdf"), w=2,h=2)
pmin <- min(c(pDT.cvh.diff.scatter$leukemia, pDT.cvh.diff$ex.vivo), na.rm = T)
pmax <- max(c(pDT.cvh.diff.scatter$leukemia, pDT.cvh.diff$ex.vivo), na.rm = T)
p <- ggplot(pDT.cvh.diff.scatter, aes(x=ex.vivo, y=leukemia)) + 
  geom_point(color="lightblue", alpha=0.5) +
  geom_text_repel(data=gg, aes(label=paste0(guide1, "-", guide2))) +
  themeNF() +
  xlab("Correlation in ex vivo") +
  ylab("Correlation in leukemia") + 
  geom_abline() +
  xlim(pmin, pmax) + ylim(pmin,pmax)
ggsaveNF(out("Correlation_CVH_HM_Scatterplot.pdf"), w=2,h=2, plot=p)


# Scatterplot ----------------------------------------------------
pDT <- data.table(fcMT, keep.rownames = TRUE)
pDT <- data.table::melt(pDT, id.vars = "rn")
pDT[, guide := gsub("\\..+$", "", variable)]
pDT[, tissue := sub("^.+?\\.", "", variable)]
pDT <- dcast.data.table(pDT[!grepl("in.vivo", tissue)], guide + rn ~ tissue, value.var = "value")
ggplot(pDT, aes(x=ex.vivo_myeloid, y=leukemia_myeloid)) + 
  themeNF() +
  stat_binhex(aes(fill=log10(..count..))) +
  facet_wrap(~guide, scales = "free", ncol = 6) +
  gradient_fill_gradient2(low="grey", high="blue")
ggsaveNF(out("Correlation_CvH_Scatterplots.pdf"), w=4,h=4)


# MDS of correlations ----------------------------------------------------
set.seed(1212)
for(tx in unique(deDT$tissue)){
  cMTx <- cMT[grepl(tx, row.names(cMT)),grepl(tx, colnames(cMT))]
  mds.res <- data.table(cmdscale(d=as.dist(1-cMTx), k=2), keep.rownames=TRUE)
  mds.res <- cbind(mds.res, setNames(data.table(do.call(rbind, strsplit(mds.res$rn, " "))), c("gene", "tissue")))
  ggplot(mds.res, aes(x=V1, y=V2, label=gene)) +
    themeNF() +
    geom_point(size=1) +
    ggrepel::geom_text_repel()+#color="black") +
    xlab("MDS dimension 1") +
    ylab("MDS dimension 2")
  ggsaveNF(out("CorrelationHM_MDS_",tx,".pdf"), w=2, h=2)
}

# Vulcano plots -----------------------------------------------------------
ggplot(deDT[guide %in% GOI], aes(x=estimate, y=pmin(30, -log10(p_value)))) + 
  themeNF() +
  facet_grid(guide ~ tissue, scale="free_y") +
  stat_binhex(aes(fill=log10(..count..)))
ggsaveNF(out("Vulcano.pdf"), w=4,h=4)

dir.create(out("Vulcanos/"))
guidex <- deDT$guide[1]
for(guidex in unique(deDT$guide)){
  pDT <- deDT[guide == guidex]
  pDT.examples <- rbind(
    pDT[order(estimate)][,head(.SD,n=5), by=c("tissue")],
    pDT[order(-estimate)][,head(.SD,n=5), by=c("tissue")],
    pDT[estimate > 0][order(p_value)][,head(.SD,n=5), by=c("tissue")],
    pDT[estimate < 0][order(p_value)][,head(.SD,n=5), by=c("tissue")]
  )
  ggplot(pDT, aes(x=estimate, y=pmin(30, -log10(p_value)))) + 
    facet_grid(guide ~ tissue, scale="free_y") +
    stat_binhex(aes(fill=log10(..count..))) +
    scale_fill_gradient(low="lightgrey", high="blue")+ 
    geom_text_repel(data=pDT.examples, aes(label=gene_id)) +
    themeNF()
  ggsaveNF(out("Vulcanos/",guidex,".pdf"), w=length(unique(pDT$tissue))*2,h=2)
}



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

ggplot(pUMAP.de[guide %in% GOI][!grepl("in.vivo", tissue)], aes(x=UMAP1, y=UMAP2)) +
  themeNF() + 
  stat_summary_hex(
    aes(z=estimate_cap),
    fun=mean.umap) +
  scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") +
  facet_grid(gsub("_", "\n", tissue)~guide) +
  xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
  xlim(-umapDT.dim1,umapDT.dim1) + ylim(-umapDT.dim2,umapDT.dim2)
ggsaveNF(out("UMAP_Values_selection.pdf"), w=3,h=0.75)


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
ggsaveNF(out("UMAP_ClusterLogFC.pdf"), w=6,h=1)


pDT <- pUMAP.de[, mean(estimate_cap, na.rm=TRUE), by=c("Cluster", "term", "guide", "tissue")]
pDT <- pDT[!grepl("in.vivo", tissue)]
pDT <- hierarch.ordering(pDT, toOrder = "Cluster", orderBy = "term", value.var = "V1")
pDT <- hierarch.ordering(pDT, orderBy = "Cluster", toOrder = "guide", value.var = "V1", aggregate = TRUE)
pDT[, tissue := gsub("_.+", "", tissue)]
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
ggsaveNF(out("UMAP_ClusterLogFC_CvH.pdf"), w=4,h=1)


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

