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
stop("Better define what to import? The next commented command takes forever")
# deDT <- fread(inDir("DEG_Statistics.tsv"))
deDT$use <- FALSE
deDT[tissue %in% grep("_s", unique(deDT$tissue), value=TRUE), use := TRUE]

# Reformat correlations
pDT <- cMT[grepl("_s$", rownames(cMT)), grepl("_s$", colnames(cMT))]
pDT <- data.table::melt(data.table(pDT, keep.rownames = TRUE), id.vars = "rn")
pDT[, rn := gsub(" ", ".", rn)]
pDT <- cbind(pDT, setNames(data.table(do.call(rbind, strsplit(gsub("_everything.+$", "", sub("\\.", "_", pDT$rn)), "_"))), paste0(c("gene", "tissue", "time"), "_x")))
pDT <- cbind(pDT, setNames(data.table(do.call(rbind, strsplit(gsub("_everything.+$", "", sub("\\.", "_", pDT$variable)), "_"))), paste0(c("gene", "tissue", "time"), "_y")))
pDT <- pDT[, -c("rn", "variable"), with=F]
cDT <- pDT

# ATAC
atacDT.lsk <- fread(paste0(gsub("Data", "CollaborationData", PATHS$LOCATIONS$DATA), "LSK_narrowPeak_consensusPeaks.boolean.annotatePeaks.extended.fc.txt"))
atacDT <- data.table::melt(
  atacDT.lsk,
  id.vars = c("chr", "start", "end", "Distance to TSS", "Gene Name"),
  measure.vars = grep("log2FC$", colnames(atacDT.lsk), value = TRUE)
  )
atacDT[, variable := gsub("LSK-(.+?)_.+", "\\L\\1", variable, perl=TRUE)]
atacDT[, variable := gsub("^(.)", "\\U\\1", variable, perl=TRUE)]


# Link UMAP and DE --------------------------------------------------------
umap.log2FC.cutoff <- 3
pUMAP.de <- merge(umapDT, setNames(data.table::melt(data.table(fcMT, keep.rownames = TRUE), id.vars = "rn"), c("gene_id", "term", "estimate")), by.x="Gene", by.y="gene_id")
pUMAP.de[, guide := sub("\\..+$", "", term)]
pUMAP.de[, tissue := sub("^.+?\\.", "", term)]
pUMAP.de[, estimate_cap := pmin(umap.log2FC.cutoff, abs(estimate)) * sign(estimate)]


# Define interesting gene sets --------------------------------------------
GOI <- list(
  BAF = c("Kmt2a", "Men1", "Brd9", "Smarcd1", "Smarcb1","Smarcd2", "Kmt2d"),
  Fig4 = c("")
)

GOI.targets <- list(
  # markers = list(markers=c("Mpo","Elane","Ctsg","S100a9","S100a8","Mtf2","Car1","Car2","Epor","Pklr","Pf4","Cpa3","Hba","Hbb","Tfrc","Cpox","Mt2","Trem2","Mefv","F13a1","Ly6c2","Fcgr3","Cd34","Kit","Itgam")),
  # tfs = list(tfs=c("Spi1","Cebpa","Cebpb","Cebpe","Klf4","Gfi1","Irf8","Mef2c","Gata2","Gata1","Klf1","Gfi1b", "Zfpm1","Mafk",
  #         "Sox6","Klf9","Lmo4","Nfia","Nfix","Nfe2","Nfe2l2","Tal1","Cited4")),
  # markers2 = list(
  #   Prog=c("Cd34","Kit","Meis1","Hoxa5","Hoxa7","Hoxa9"),
  #   Mye=c("Elane","Mpo","Csf1r","Spi1","Cebpa","Cebpb","Cebpe","Gfi1","Irf8"," Mef2c"),
  #   "Mega-ery"=c("Car1","Gypa","Hbb-bt","Hbb-1","Hba-a1","Gata2","Gata1","Gfi1b","Klf1","Sox6","Klf9","Nfe2l2","Cited4","Tal1"),
  #   "Mega-only"=c("Itga2b","Cavin2","Vwf","Gp1ba","Gp1bb","Nfix","Nfia"),
  #   Baso=c("Cpa3","Prss34","Fcer1a"),
  #   "B-cell"=c("Il7r","Cd19","Pax5","Ebf1"))
  markers2 = with(fread("metadata/FIGS_02_DE_Genes.tsv", check.names = TRUE), split(Gene.list, Programme))
)
GOI.targets$LSC <- list(LSC=c("Hif1a", "Myc", "Bcat1", grep("Hox", unique(deDT$gene_id), value=TRUE)))



dla.vulcano.genes <- fread("metadata/FIGS_VulcanoGenes.tsv", header = FALSE)
dla.vulcano.genes <- setdiff(unique(do.call(c, dla.vulcano.genes)), "")

# dla.factors <- fread("metadata/FIGS_Order_Fig3_v2.tsv")
# dla.factors <- lapply(as.list(dla.factors), function(ll) setdiff(ll, ""))
dla.factors <- list(
  supp=fread("metadata/FIGS_02_CFs.supp.txt")$Factor,
  main=fread("metadata/FIGS_02_CFs.main.txt")$Factor
)

# SETUP ENDS HERE ---------------------------------------------------------




# Plot specific genes -----------------------------------------------------
tx <- deDT[use == TRUE]$tissue[1]
for(tx in unique(deDT[use == TRUE]$tissue)){
  pDT <- deDT[tissue == tx]
  lnam <- names(GOI.targets)[2]
  for(lnam in names(GOI.targets)){
    goi <- rbindlist(lapply(GOI.targets[[lnam]], data.table), idcol = "celltype")
    goi[! V1 %in% pDT$gene_id]
    xDT <- merge(pDT, goi, by.x="gene_id", by.y="V1")
    xDT$gene_id <- factor(xDT$gene_id, levels = goi$V1)
    xDT$celltype <- factor(xDT$celltype, levels = unique(goi$celltype))
    #xDT$guide <- factor(xDT$guide, levels = GOI$BAF)
    xDT <- hierarch.ordering(xDT, "guide", "gene_id", value.var = "estimate")
    # xDT <- hierarch.ordering(xDT, "gene_id", "guide", value.var = "estimate")
    h=length(unique(xDT$gene_id)) * 0.05 + length(unique(xDT$celltype)) * 0.05 + 0.5
    w=length(unique(xDT$guide)) * 0.07 + 0.7
    ggplot(xDT, aes(x=guide, y=gene_id, color=estimate, size=pmin(5, -log10(q_value)))) + 
      themeNF(rotate=TRUE) +
      geom_point() +
      scale_color_gradient2(name="log2FC", low="blue", high="red") +
      scale_size_continuous(name="-log10(padj)", range = c(0,5)) +
      ggtitle(tx) +
      facet_grid(celltype ~ ., space = "free", scales = "free") +
      ylab("Gene expression") + xlab("CRISPR targets")
    ggsaveNF(out("GeneDE_", tx, "_", lnam, "_allGuides.pdf"),w=w,h=h, guides = TRUE)
    
    (dla.nam <- names(dla.factors)[3])
    for(dla.nam in names(dla.factors)){
      dla <- dla.factors[[dla.nam]]
      xDT2 <- xDT[guide %in% dla]
      if(nrow(xDT2) == 0) next
      xDT2$guide <- factor(xDT2$guide, levels=dla)
      w=length(unique(xDT2$guide)) * 0.07 + 0.7
      ggplot(xDT2, aes(x=guide, y=gene_id, color=estimate, size=pmin(5, -log10(q_value)))) + 
        themeNF(rotate=TRUE) +
        geom_point() +
        scale_color_gradient2(name="log2FC", low="blue", high="red") +
        scale_size_continuous(name="-log10(padj)", range = c(0,5)) +
        ggtitle(tx) +
        facet_grid(celltype ~ ., space = "free", scales = "free") +
        ylab("Gene expression") + xlab("CRISPR targets")
      ggsaveNF(out("GeneDE_", tx, "_", lnam, "_", dla.nam,".pdf"),w=w,h=h, guides = TRUE)}
  }
}


# ATAC-seq comparison -----------------------------------------------------
atacDT[order(abs(get("Distance to TSS")), decreasing = FALSE)][, head(.SD, n=1), by=c("Gene Name")]
pDT <- merge(
  deDT[tissue == "ex.vivo_7d_everything_s"][, c("gene_id", "guide", "estimate"), with=FALSE],
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


# Vulcano plots -----------------------------------------------------------
dir.create(out("Vulcanos/"))
tissuex <- deDT$tissue[1]
for(tissuex in unique(deDT$tissue)){
  pDT <- deDT[tissue == tissuex]
  # pDT.examples <- rbind(
  #   pDT[order(estimate)][,head(.SD,n=5), by=c("tissue")],
  #   pDT[order(-estimate)][,head(.SD,n=5), by=c("tissue")],
  #   pDT[estimate > 0][order(p_value)][,head(.SD,n=5), by=c("tissue")],
  #   pDT[estimate < 0][order(p_value)][,head(.SD,n=5), by=c("tissue")]
  # )
  pDT.examples <- pDT[gene_id %in% dla.vulcano.genes]
  cnt <- length(unique(pDT$guide))
  ncol=5
  p <- ggplot(pDT, aes(x=estimate, y=pmin(30, -log10(p_value)))) + 
    facet_wrap(~ guide, scale="free", ncol = ncol) +
    stat_binhex(aes(fill=log10(..count..))) +
    scale_fill_gradient(low="lightgrey", high="blue")+ 
    geom_text_repel(data=pDT.examples, aes(label=gene_id)) +
    themeNF()
  ggsaveNF(out("Vulcanos/Vulcano_",tissuex,".pdf"), h=ceiling(cnt/ncol),w=ncol, limitsize=FALSE, plot=p)
}



# Correlation Heatmaps ----------------------------------------------------
pDT <- cDT[tissue_x == tissue_y & time_x == time_y]
pDT <- pDT[time_x %in% c("7d", "14d")]
pDT <- pDT[gene_x != gene_y]
pDT <- hierarch.ordering(pDT, "gene_x", "gene_y", "value", TRUE)
pDT$gene_y <- factor(pDT$gene_y, levels = levels(pDT$gene_x))
p <- ggplot(pDT, aes(x=gene_x, y=gene_y, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low="blue", high="red")+
  facet_grid(. ~ tissue_x) +
  themeNF(rotate=TRUE) +
  xlab("") + ylab("")
ggsaveNF(out("ExpressionCorrelation.pdf"), w=4,h=2, plot=p)

pDT <- pDT[gene_x %in% dla.factors$`CRISPR list-Main` & gene_y %in% dla.factors$`CRISPR list-Main`]
p <- ggplot(pDT, aes(x=gene_x, y=gene_y, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low="blue", high="red")+
  facet_grid(. ~ tissue_x) +
  themeNF(rotate=TRUE) +
  xlab("") + ylab("")
ggsaveNF(out("ExpressionCorrelation_Figure2.pdf"), w=2,h=1, plot=p)


# Repressors ------------------------------------------------------------

# Correlation HM
gg <- c("Setdb1", "Atf7ip", "Chd4", "Rbbp4", "Hdac3", "Hdac1")
pDT <- cDT[tissue_x == "in.vivo" & time_x == "14d"]
pDT <- pDT[time_y == time_x & tissue_x == tissue_y]
pDT <- pDT[gene_x %in% gg & gene_y %in% gg]
pDT <- pDT[gene_x != gene_y]
pDT <- hierarch.ordering(pDT, "gene_x", "gene_y", "value", TRUE)
pDT$gene_y <- factor(pDT$gene_y, levels = levels(pDT$gene_x))
p <- ggplot(pDT, aes(x=gene_x, y=gene_y, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low="blue", high="red", limits=c(-1,1))+
  facet_grid(. ~ tissue_x) +
  themeNF(rotate=TRUE) +
  xlab("") + ylab("")
ggsaveNF(out("Repressors_ExpressionCorrelation.pdf"), w=0.8,h=0.8, plot=p)

# Genes
pDT <- deDT[tissue == "in.vivo_14d_everything_s"]
pDT <- pDT[guide %in% gg]
pDT[, direction := ifelse(estimate < 0, "down", "up")]
xDT <- pDT[gene_id %in% pDT[q_value < 0.05 & abs(estimate) > 1][, .N, by=c("gene_id", "direction")][direction == "up"][!grepl("Rik$", gene_id)][order(N, decreasing = TRUE)][1:50]$gene_id]
xDT <- hierarch.ordering(xDT, "guide", "gene_id", value.var = "estimate")
xDT <- hierarch.ordering(xDT, "gene_id", "guide", value.var = "estimate")
# xDT <- hierarch.ordering(xDT, "gene_id", "guide", value.var = "estimate")
h=length(unique(xDT$gene_id)) * 0.05 + 0.5
w=length(unique(xDT$guide)) * 0.07 + 0.7
ggplot(xDT, aes(x=guide, y=gene_id, color=estimate, size=pmin(5, -log10(q_value)))) + 
  themeNF(rotate=TRUE) +
  geom_point() +
  scale_color_gradient2(name="log2FC", low="blue", high="red") +
  scale_size_continuous(name="-log10(padj)", range = c(0,5)) +
  #facet_grid(celltype ~ ., space = "free", scales = "free") +
  ylab("Gene expression") + xlab("CRISPR targets")
ggsaveNF(out("Repressors_CommonGenes.pdf"),w=w,h=h, guides = TRUE)


# Correlation in BAF ------------------------------------------------------
# pDT <- data.table::melt(data.table(cMT, keep.rownames = TRUE), id.vars = "rn")
# pDT <- pDT[(grepl("ex.vivo", rn) & grepl("ex.vivo", variable)) | (grepl("in.vivo", rn) & grepl("in.vivo", variable))]
# pDT <- pDT[grepl("_myeloid", rn) & grepl("_myeloid", variable)]
# pDT[, tissue := ifelse(grepl("in.vivo", rn), "in vivo", "ex vivo")]
# pDT[, variable := gsub("\\.", " ", variable)]
# pDT[, variable := gsub(" .+$", "", variable)]
# pDT[, rn := gsub("\\.", " ", rn)]
# pDT[, rn := gsub(" .+$", "", rn)]
# pDT <- cDT[tissue_x %in% c("ex.vivo", "in.vivo") & tissue_y %in% c("ex.vivo", "in.vivo")][gene_x %in% GOI$BAF & gene_y %in% GOI$BAF]
# #pDT <- pDT[rn != variable]
# pDT$gene_x <- factor(pDT$rn, levels = GOI$BAF)
# pDT$gene_y <- factor(pDT$variable, levels = GOI$BAF)
# # pDT <- hierarch.ordering(pDT, toOrder = "rn", orderBy = "variable", value.var = "value")
# # pDT <- hierarch.ordering(pDT, toOrder = "variable", orderBy = "rn", value.var = "value")
# ggplot(pDT, aes(x=gene_x, y=gene_y, fill=value)) + 
#   geom_tile() + 
#   facet_grid(. ~ tissue) +
#   scale_fill_gradient2(low="blue", high="red") +
#   themeNF(rotate=TRUE) +
#   xlab("") + ylab("")
# ggsaveNF(out("Phenocopy_BAF.pdf"), w=2,h=1)


# Correlation Cancer vs Healthy ----------------------------------------------------
# pDT <- data.table::melt(data.table(cMT, keep.rownames = TRUE), id.vars = "rn")
# pDT[, guide1 := sub(" .+$", "", rn)]
# pDT[, tissue1 := sub("^.+? ", "", rn)]
# pDT[, guide2 := sub("\\..+$", "", variable)]
# pDT[, tissue2 := sub("^.+?\\.", "", variable)]
# pDT <- hierarch.ordering(pDT, toOrder = "guide1", orderBy = "variable", value.var = "value", aggregate = TRUE)
# pDT <- hierarch.ordering(pDT, toOrder = "guide2", orderBy = "rn", value.var = "value", aggregate = TRUE)
# 
# # Overall correlations
# # ggplot(pDT, aes(x=guide1, y=guide2, fill=value)) + 
# #   themeNF() +
# #   geom_tile() + 
# #   facet_grid(tissue2 ~ tissue1, space = "free", scales = "free") +
# #   scale_fill_gradient2(low="blue", high="red") +
# #   xlab("") + ylab("") +
# #   xRot()
# # ggsaveNF(out("CorrelationHM.pdf"), w=3,h=3)
# 
# # Cancer vs Healthy correlations
# pDT[,variable := as.character(variable)]
# pDT[make.names(rn) == make.names(variable), value := NA]
# tissues.use <- c("ex.vivo_noMixscape", "leukemia_noMixscape")
# pDT.cvh <- pDT[tissue1 %in% tissues.use & tissue2 %in% tissues.use]
# double.guides <- pDT.cvh[guide1 == guide2][!is.na(value)]$guide1
# pDT.cvh <- pDT.cvh[guide1 %in% double.guides]
# pDT.cvh <- pDT.cvh[guide2 %in% double.guides]
# pDT.cvh[,tissue1 := gsub("_.+","", tissue1)]
# pDT.cvh[,tissue2 := gsub("_.+","", tissue2)]
# 
# # Cancer vs Healthy same factor Barplot
# pDT.sum <- pDT.cvh[guide1 == guide2 & tissue1 == "ex.vivo" & tissue1 != tissue2]
# pDT.sum$guide <- factor(pDT.sum$guide1, levels = pDT.sum[order(value)]$guide1)
# pDT.sum <- pDT.sum[guide != "Cbx3"]
# ggplot(pDT.sum, aes(x=guide, y=value)) + 
#   themeNF() +
#   geom_bar(stat="identity") + 
#   ylab("Correlation normal vs cancer") +
#   coord_flip() +
#   facet_grid(. ~ "x") +
#   ylab("")
# ggsaveNF(out("Correlation_CVH_bars.pdf"), w=1,h=2, guides = TRUE)
# 
# 
# 
# # Number of up- and downregualted genes
# pDT.up_down <- copy(deDT[tissue %in% tissues.use])
# pDT.up_down[, tissue := gsub("_noMixscape", "", tissue)]
# pDT.up_down[, direction := ifelse(estimate > 0, "up", "down")]
# pDT.up_down <- pDT.up_down[q_value < 0.05][abs(estimate) > 0.5][,length(unique(gene_id)), by=c("guide", "direction", "tissue")]
# pDT.up_down[direction == "down", V1 := -V1]
# pDT.up_down <- pDT.up_down[guide %in% pDT.sum$guide]
# pDT.up_down$guide <- factor(pDT.up_down$guide, levels=levels(pDT.sum$guide))
# #pDT[direction == "down", V1 := -V1]
# p <- ggplot(pDT.up_down, aes(x=guide, y=V1, fill=direction)) +
#   themeNF() +
#   geom_bar(stat="identity", position="stack") +
#   coord_flip() +
#   facet_grid(. ~ tissue)
# ggsaveNF(out("NumberOfGenes.pdf"), w=1,h=2, plot=p)
# #ggsaveNF(out("NumberOfGenes_log.pdf"), w=1,h=2, plot=p + scale_y_log10())
# 
# 
# 
# # check relationship of correlation cvh to number of cells
# pDT.sum.control <- merge(
#   pDT.sum,
#   dcast.data.table(annList[mixscape_class.global == "KO", .N, by=c("gene", "tissue")][tissue != "in.vivo"], gene ~ tissue, value.var="N"),
#   by.x="guide1", by.y="gene")
# ggplot(pDT.sum.control,
#        aes(x=ex.vivo, y=leukemia, label=guide1)) +
#   geom_point(aes(color=value), size=3) +
#   geom_text_repel() +
#   themeNF(rotate=TRUE) +
#   scale_color_gradient2(low="blue", high="red") +
#   scale_y_log10(limits=c(10,1e4)) +
#   scale_x_log10(limits=c(10,1e4)) +
#   geom_abline()
# ggsaveNF(out("Correlation_CVH_bars_GuideCellNumbers.pdf"), w=2,h=2)
# 
# 
# # Cancer vs Healthy Heatmap
# # ggplot(pDT.cvh,
# #        aes(
# #          x=paste(tissue1),
# #          y=paste(tissue2),
# #          fill=value)) +
# #   themeNF() +
# #   geom_tile() +
# #   facet_grid(guide2 ~ guide1, space = "free", scales = "free", switch = "both") +
# #   scale_fill_gradient2(low="blue", high="red") +
# #   xlab("") + ylab("") +
# #   theme(strip.text.x = element_text(angle=90)) +
# #   theme(strip.text.y.left = element_text(angle=0)) +
# #   theme(panel.spacing = unit(0, "lines")) +
# #   xRot()
# # ggsaveNF(out("Correlation_CVH_HM.pdf"), w=3,h=3)
# 
# 
# # Cancer vs Healthy Heatmap of changes phenocopies
# pDT.cvh.diff <- dcast.data.table(pDT.cvh[tissue1 == tissue2], guide1 + guide2 ~ tissue1, value.var = "value")
# pDT.cvh.diff[, diff := leukemia - ex.vivo]
# pDT.cvh.diff <- hierarch.ordering(pDT.cvh.diff, "guide1", "guide2", "diff")
# pDT.cvh.diff <- hierarch.ordering(pDT.cvh.diff, "guide2", "guide1", "diff")
# pDT.cvh.diff[, diff.grp := paste("sign", ifelse(abs(diff) > 0.1, sign(diff), 0))]
# ggplot(pDT.cvh.diff, aes(x=guide1, y=guide2, fill=ex.vivo, size=abs(diff), color=diff.grp)) + 
#   geom_point(shape=21) +
#   scale_color_manual(name="Change\ndirection", values=c("blue", "grey", "red", "grey")) +
#   scale_fill_gradient2(name="Correlation\nex vivo", low="blue", high="red") +
#   scale_size_continuous(name="Change\nmagnitude") +
#   themeNF(rotate = TRUE) +
#   ggtitle("Change in correlation leukemia vs ex vivo") +
#   xlab("CRISPR targets") + ylab("CRISPR targets")
# ggsaveNF(out("Correlation_CVH_HM_Details_HM.pdf"), w=2,h=2)
# 
# 
# # Cancer vs Healthy changes phenocopies control for number of cells
# xDT <- dcast.data.table(annList[mixscape_class.global == "KO", .N, by=c("gene", "tissue")][tissue != "in.vivo"], gene ~ tissue, value.var="N")
# xDT2 <- pDT.cvh.diff[,.(meandiff = mean(diff, na.rm=TRUE)), by=c("guide1")]
# xDT2[, gene := as.character(guide1)]
# xDT <- merge.data.table(xDT, xDT2, by="gene")
# ggplot(xDT,
#        aes(x=ex.vivo, y=leukemia, label=gene)) +
#   geom_point(aes(color=meandiff), size=3) +
#   geom_text_repel() +
#   themeNF(rotate=TRUE) +
#   scale_color_gradient2(low="blue", high="red") +
#   scale_y_log10() +
#   scale_x_log10() +
#   geom_abline()
# ggsaveNF(out("Correlation_CVH_HM_Details_HM_GuideCellNumbers.pdf"), w=2,h=2)
# 
# # Cancer vs Healthy of changes phenocopies Heatmap examples
# pDT.cvh.diff.2 <- copy(pDT.cvh.diff[sign(ex.vivo) != sign(leukemia)][order(-abs(diff))][1:50])
# pDT.cvh.diff.2 <- pDT.cvh.diff.2[t(apply(as.matrix(pDT.cvh.diff.2[,c("guide1", "guide2")]), 1, order))[,1] == 1]
# pDT.cvh.diff.2[, pair := paste0(guide1, "-", guide2)]
# pDT.cvh.diff.2$pair <- factor(pDT.cvh.diff.2$pair, levels=pDT.cvh.diff.2[order(diff)]$pair)
# pDT.cvh.diff.2 <- melt(pDT.cvh.diff.2, id.vars = c("pair"), measure.vars=c("ex.vivo", "leukemia", "diff"))
# ggplot(pDT.cvh.diff.2, aes(x=variable, y=pair, fill=value)) + 
#   themeNF(rotate=TRUE) +
#   geom_tile() +
#   scale_fill_gradient2(low="blue", high="red") +
#   xlab("Correlation")
# ggsaveNF(out("Correlation_CVH_HM_Details_HM_Examples.pdf"), w=1,h=2)
# 
# # Cancer vs Healthy of changes phenocopies Scatterplot
# pDT.cvh.diff.scatter <- pDT.cvh.diff[t(apply(as.matrix(pDT.cvh.diff[,c("guide1", "guide2")]), 1, order))[,1] == 1]
# gg <- unique(rbind(
#   pDT.cvh.diff.scatter[ex.vivo > 0.4 | ex.vivo < -0.15],
#   pDT.cvh.diff.scatter[sign(leukemia) != sign(ex.vivo)][order(abs(diff), decreasing = TRUE)][1:20]
# ))
# # ggplot(pDT.cvh.diff, aes(x=ex.vivo, y=diff)) + 
# #   geom_point() +
# #   geom_text_repel(data=gg, aes(label=paste0(guide1, "-", guide2))) +
# #   themeNF() +
# #   xlab("Correlation in ex vivo") +
# #   ylab("Difference in correlation\nleukemia vs ex vivo")
# # ggsaveNF(out("Correlation_CVH_HM_Details.pdf"), w=2,h=2)
# pmin <- min(c(pDT.cvh.diff.scatter$leukemia, pDT.cvh.diff$ex.vivo), na.rm = T)
# pmax <- max(c(pDT.cvh.diff.scatter$leukemia, pDT.cvh.diff$ex.vivo), na.rm = T)
# p <- ggplot(pDT.cvh.diff.scatter, aes(x=ex.vivo, y=leukemia)) + 
#   geom_point(color="lightblue", alpha=0.5) +
#   geom_text_repel(data=gg, aes(label=paste0(guide1, "-", guide2))) +
#   themeNF() +
#   xlab("Correlation in ex vivo") +
#   ylab("Correlation in leukemia") + 
#   geom_abline() +
#   xlim(pmin, pmax) + ylim(pmin,pmax)
# ggsaveNF(out("Correlation_CVH_HM_Scatterplot.pdf"), w=2,h=2, plot=p)
# 
# 
# # Scatterplot ----------------------------------------------------
# pDT <- data.table(fcMT, keep.rownames = TRUE)
# pDT <- data.table::melt(pDT, id.vars = "rn")
# pDT[, guide := gsub("\\..+$", "", variable)]
# pDT[, tissue := sub("^.+?\\.", "", variable)]
# pDT <- dcast.data.table(pDT, guide + rn ~ tissue, value.var = "value")
# ggplot(pDT, aes(x=ex.vivo_noMixscape, y=leukemia_noMixscape)) + 
#   themeNF() +
#   stat_binhex(aes(fill=log10(..count..))) +
#   facet_wrap(~guide, scales = "free", ncol = 6)
#   #scale_fill_gradient2(low="grey", high="blue")
# ggsaveNF(out("Correlation_CvH_Scatterplots.pdf"), w=4,h=4)


# # MDS of correlations ----------------------------------------------------
# set.seed(1212)
# for(tx in unique(deDT$tissue)){
#   cMTx <- cMT[grepl(tx, row.names(cMT)),grepl(tx, colnames(cMT))]
#   mds.res <- data.table(cmdscale(d=as.dist(1-cMTx), k=2), keep.rownames=TRUE)
#   mds.res <- cbind(mds.res, setNames(data.table(do.call(rbind, strsplit(mds.res$rn, " "))), c("gene", "tissue")))
#   ggplot(mds.res, aes(x=V1, y=V2, label=gene)) +
#     themeNF() +
#     geom_point(size=1) +
#     ggrepel::geom_text_repel()+#color="black") +
#     xlab("MDS dimension 1") +
#     ylab("MDS dimension 2")
#   ggsaveNF(out("CorrelationHM_MDS_",tx,".pdf"), w=2, h=2)
}


# # Regulatory Map UMAP -----------------------------------------------------
# ggplot(umapDT, aes(x=UMAP1, y=UMAP2, color=factor(Cluster))) + 
#   geom_point(shape=1) + 
#   theme_bw(12) +
#   geom_label(data=umapDT[, .(UMAP1=median(UMAP1), UMAP2=median(UMAP2)), by="Cluster"], aes(label=Cluster))
# ggsave(out("UMAP","_Clusters.pdf"), w=6,h=5)
# 
# #HEX
# ggplot(umapDT, aes(x=UMAP1, y=UMAP2)) + 
#   geom_hex(bins=100) +
#   theme_bw(12) +
#   geom_label(data=umapDT[, .(UMAP1=median(UMAP1), UMAP2=median(UMAP2)), by="Cluster"], aes(label=Cluster))
# ggsave(out("UMAP","_Clusters_hex.pdf"), w=6,h=5)
# 
# 
# 
# # Plot values on UMAP -----------------------------------------------------
# # tn <- length(unique(pUMAP.de$guide))
# mean.umap <- function(x){mean(x, na.rm=TRUE)}
# 
# ggplot(pUMAP.de, aes(x=UMAP1, y=UMAP2)) +
#   themeNF() + 
#   stat_summary_hex(
#     aes(z=estimate_cap),
#     fun=mean.umap) +
#   scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") +
#   facet_grid(guide ~ tissue) +
#   xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
#   xlim(-umapDT.dim1,umapDT.dim1) + ylim(-umapDT.dim2,umapDT.dim2)
# ggsaveNF(out("UMAP_Values.pdf"), w=2,h=8)
# 
# 
# # Plot values by UMAP cluster -----------------------------------------------------
# pDT <- pUMAP.de[, mean(estimate_cap, na.rm=TRUE), by=c("Cluster", "term", "guide", "tissue")]
# pDT <- hierarch.ordering(pDT, toOrder = "Cluster", orderBy = "term", value.var = "V1")
# pDT <- hierarch.ordering(pDT, orderBy = "Cluster", toOrder = "guide", value.var = "V1", aggregate = TRUE)
# ggplot(pDT, aes(y=factor(Cluster), x=tissue, fill=V1)) + 
#   themeNF() + 
#   geom_tile() +
#   scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") +
#   facet_grid(. ~ guide, scales = "free", space = "free", switch = "x") +
#   ylab("Gene modules") +
#   theme(strip.text.x = element_text(angle=90)) +
#   theme(panel.spacing = unit(0.01, "cm")) +
#   xlab("") +
#   xRot()
# ggsaveNF(out("UMAP_ClusterLogFC.pdf"), w=6,h=1)
# 
# 
# pDT <- pUMAP.de[, mean(estimate_cap, na.rm=TRUE), by=c("Cluster", "term", "guide", "tissue")]
# pDT <- pDT[!grepl("in.vivo", tissue)]
# pDT <- hierarch.ordering(pDT, toOrder = "Cluster", orderBy = "term", value.var = "V1")
# pDT <- hierarch.ordering(pDT, orderBy = "Cluster", toOrder = "guide", value.var = "V1", aggregate = TRUE)
# pDT[, tissue := gsub("_.+", "", tissue)]
# ggplot(pDT, aes(y=factor(Cluster), x=tissue, fill=V1)) + 
#   themeNF() + 
#   geom_tile() +
#   scale_fill_gradient2(high="#e31a1c",mid="#ffffff", low="#1f78b4") +
#   facet_grid(. ~ guide, scales = "free", space = "free", switch = "x") +
#   ylab("Gene modules") +
#   theme(strip.text.x = element_text(angle=90)) +
#   theme(panel.spacing = unit(0.01, "cm")) +
#   xlab("") +
#   xRot()
# ggsaveNF(out("UMAP_ClusterLogFC_CvH.pdf"), w=4,h=1)
# 
# 
# # GSEA --------------------------------------------------------------------
# fish.res <- copy(gseaDT)
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

