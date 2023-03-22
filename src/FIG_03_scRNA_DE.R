source("src/00_init.R")
base.dir <- "FIG_03_scRNA_DE/"
out <- dirout(base.dir)

require(ggrepel)
require(WriteXLS)


# Read data ---------------------------------------------------------------
inDir <- dirout_load("SCRNA_40_01_DE_summary")

annList <- lapply(PATHS$SCRNA$MONOCLE.NAMES, function(tx){
  ann <- fread(dirout_load(paste0("SCRNA_20_Summary/", tx, "_monocle.singleR"))("Annotation.tsv"))
  ann[, gene := gsub("_.+", "", guide)]
  ann
})
annList <- rbindlist(annList, fill=TRUE)
annList <- annList[Clusters != "unclear"]

cMT <- as.matrix(read.csv(inDir("DEG_Cor.csv"), row.names = 1))
fcMT <- as.matrix(read.csv(inDir("FoldChanges.csv"), row.names=1))
# umapDT <- fread(inDir("RegulatoryMap_UMAP_","all",".tsv"))
# umapDT.dim1 <- floor(max(abs(umapDT$UMAP1))) + 0.5
# umapDT.dim2 <- floor(max(abs(umapDT$UMAP2))) + 0.5
# gseaDT.umap <- fread(inDir("UMAP_GSEA.tsv"))
gseaDT <- readRDS(dirout_load("SCRNA_41_01_GSEA")("FGSEA.RDS"))
deDT <- readRDS(inDir("DEG_Statistics_simple.RDS"))
deDT$use <- TRUE

deDT2 <- fread(inDir("DEG_Statistics_significant.tsv"))
#deDT[tissue %in% grep("_s", unique(deDT$tissue), value=TRUE), use := TRUE]

# Reformat correlations
pDT <- cMT #[grepl("_s$", rownames(cMT)), grepl("_s$", colnames(cMT))]
pDT <- data.table::melt(data.table(pDT, keep.rownames = TRUE), id.vars = "rn")
pDT[, rn := gsub(" ", ".", rn)]
pDT <- cbind(pDT, setNames(data.table(do.call(rbind, strsplit(gsub("_everything.+$", "", sub("\\.", "_", pDT$rn)), "_"))), paste0(c("gene", "tissue", "time", "type"), "_x")))
pDT <- cbind(pDT, setNames(data.table(do.call(rbind, strsplit(gsub("_everything.+$", "", sub("\\.", "_", pDT$variable)), "_"))), paste0(c("gene", "tissue", "time", "type"), "_y")))
pDT <- pDT[, -c("rn", "variable"), with=F]
cDT <- pDT

# ATAC
# atacDT.lsk <- fread(paste0(gsub("Data", "CollaborationData", PATHS$LOCATIONS$DATA), "LSK_narrowPeak_consensusPeaks.boolean.annotatePeaks.extended.fc.txt"))
# atacDT <- data.table::melt(
#   atacDT.lsk,
#   id.vars = c("chr", "start", "end", "Distance to TSS", "Gene Name"),
#   measure.vars = grep("log2FC$", colnames(atacDT.lsk), value = TRUE)
#   )
# atacDT[, variable := gsub("LSK-(.+?)_.+", "\\L\\1", variable, perl=TRUE)]
# atacDT[, variable := gsub("^(.)", "\\U\\1", variable, perl=TRUE)]


# ChIP --------------------------------------------------------------------
# chip.targets <- fread("metadata/FIGS_06_ChIPtargetsJulen.txt")

# Link UMAP and DE --------------------------------------------------------
# umap.log2FC.cutoff <- 3
# pUMAP.de <- merge(umapDT, setNames(data.table::melt(data.table(fcMT, keep.rownames = TRUE), id.vars = "rn"), c("gene_id", "term", "estimate")), by.x="Gene", by.y="gene_id")
# pUMAP.de[, guide := sub("\\..+$", "", term)]
# pUMAP.de[, tissue := sub("^.+?\\.", "", term)]
# pUMAP.de[, estimate_cap := pmin(umap.log2FC.cutoff, abs(estimate)) * sign(estimate)]


# Define interesting gene sets --------------------------------------------
GOI <- list(
  BAF = c("Kmt2a", "Men1", "Brd9", "Smarcd1", "Smarcb1","Smarcd2", "Kmt2d"),
  Fig4 = c("")
)

GOI.targets <- fread("metadata/FIGS_02_DE_Genes.tsv", check.names = TRUE)
GOI.targets.TFs <- c("Runx1","Runx2","Cebpa","Fos","Jun","Nfil3","Hlf","Atf4","Atf1","Atf3","Elf4","Spib","Spi1","Gabpa","Irf1","Irf8","Gata1","Gata2","Nfia","Nfix","Nrf1","Klf1","Klf6","Klf9","Klf10")

dla.vulcano.genes <- fread("metadata/FIGS_VulcanoGenes.tsv", header = FALSE)
dla.vulcano.genes <- setdiff(unique(do.call(c, dla.vulcano.genes)), "")

# Factors
dla.table <- fread("metadata/FIGS_02_CFs.main.txt")
dla.healthy <- list(
  supp=dla.table$CF,
  main=with(dla.table, CF[LargePlot]),
  main.small=with(dla.table, CF[SmallPlot]),
  main.Feb28=fread("metadata/FIGS_02_CFs.main_Feb28.txt", header = FALSE)$V1
)
dla.healthy$all <- setdiff(sort(unique(annList$gene)), "NTC")
for(x in names(dla.healthy)){ dla.healthy[[x]] <- setdiff(dla.healthy[[x]], "Smarcb1")}
dla.factors <- dla.healthy

# SETUP ENDS HERE ---------------------------------------------------------



# . Supplementary tables ----------------------------------------------------

# DE Genes
for(tx in unique(deDT2$tissue)){
  pDT <- deDT2[tissue == tx][q_value < 0.01][abs(estimate) > 1][order(guide, q_value, -estimate)][,c("guide", "gene_id", "estimate", "q_value"),with=F]
  pDT[, q_value := format(q_value, scientific=TRUE, digits=3)]
  for(cx in 1:ncol(pDT)){if(is.numeric(pDT[[cx]])) pDT[[cx]] <- round(pDT[[cx]],3)}
  colnames(pDT) <- c("CF KO", "gene", "log2 fold change", "adjusted p-value")
  WriteXLS(x=pDT, ExcelFileName=out("Supplementary_Table_DE_",tx,".xls"), AdjWidth=TRUE, BoldHeaderRow=TRUE, FreezeRow=1, SheetNames="Table")
  write.tsv(pDT, out("Supplementary_Table_DE_",tx,".tsv"))
}

# GSEA
pDT <- gseaDT[tissue == "in.vivo_14d_noClusters"][padj < 0.05][order(grp, padj, -NES)][,c("grp", "pathway", "db", "NES", "padj"),with=F]
pDT[, padj := format(padj, scientific=TRUE, digits=3)]
for(cx in 1:ncol(pDT)){if(is.numeric(pDT[[cx]])) pDT[[cx]] <- round(pDT[[cx]],3)}
colnames(pDT) <- c("CF KO", "gene set", "database", "log2 fold change", "adjusted p-value")
WriteXLS(x=pDT, ExcelFileName=out("Supplementary_Table_GSEA.xls"), AdjWidth=TRUE, BoldHeaderRow=TRUE, FreezeRow=1, SheetNames="Table")
write.tsv(pDT, out("Supplementary_Table_GSEA.tsv"))



# Plot specific genes -----------------------------------------------------
(tx <- deDT[use == TRUE]$tissue[1])
for(tx in unique(deDT[use == TRUE]$tissue)){
  pDT <- deDT[tissue == tx]
  (lnam <- "x")
  
  xDT <- merge(pDT, GOI.targets, by.x="gene_id", by.y="Gene.list")
  xDT$gene_id <- factor(xDT$gene_id, levels = GOI.targets$Gene.list)
  xDT$celltype <- factor(xDT$Programme, levels = unique(GOI.targets$Programme))
  #xDT$guide <- factor(xDT$guide, levels = GOI$BAF)
  xDT <- hierarch.ordering(xDT, "guide", "gene_id", value.var = "estimate")
  # xDT <- hierarch.ordering(xDT, "gene_id", "guide", value.var = "estimate")
  w=length(unique(xDT$gene_id)) * 0.05 + length(unique(xDT$celltype)) * 0.05 + 0.5
  h=length(unique(xDT$guide)) * 0.05 + 0.5
  xDT[abs(estimate) > 5, estimate := pmin(abs(estimate),5) * sign(estimate)]
  ggplot(xDT, aes(y=guide, x=gene_id, color=estimate, size=pmin(5, -log10(q_value)))) + 
    themeNF(rotate=TRUE) +
    geom_point() +
    scale_color_gradient2(name="log2FC", low="blue", high="red") +
    scale_size_continuous(name="-log10(padj)", range = c(0,5)) +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(vjust=1, hjust=0)) +
    ggtitle(tx) +
    facet_grid(. ~ celltype, space = "free", scales = "free") +
    xlab("Gene expression") + ylab("CRISPR targets")
  ggsaveNF(out("GeneDE_", tx, "_", lnam, "_allGuides.pdf"),w=w,h=h, guides = TRUE)
  
  (dla.nam <- names(dla.factors)[3])
  for(dla.nam in names(dla.factors)){
    dla <- dla.factors[[dla.nam]]
    xDT2 <- xDT[guide %in% dla]
    if(nrow(xDT2) == 0) next
    xDT2$guide <- factor(xDT2$guide, levels=rev(dla))
    h=length(unique(xDT2$guide)) * 0.05 + 0.5
    ggplot(xDT2, aes(y=guide, x=gene_id, color=estimate, size=pmin(5, -log10(q_value)))) + 
      themeNF(rotate=TRUE) +
      geom_point() +
      scale_color_gradient2(name="log2FC", low="blue", high="red") +
      scale_size_continuous(name="-log10(padj)", range = c(0,5)) +
      scale_x_discrete(position = "top") +
      theme(axis.text.x = element_text(vjust=1, hjust=0)) +
      ggtitle(tx) +
      facet_grid(. ~ celltype, space = "free", scales = "free") +
      xlab("Gene expression") + ylab("CRISPR targets")
    ggsaveNF(out("GeneDE_", tx, "_", lnam, "_", dla.nam,".pdf"),w=w,h=h, guides = TRUE)}
}


# Specific TF genes -------------------------------------------------------
# Plot specific genes -----------------------------------------------------
(tx <- deDT2$tissue[1])
for(tx in unique(deDT2[grepl("ex.vivo", tissue)]$tissue)){
  pDT <- deDT2[tissue == tx]
  (lnam <- "x")
  
  xDT <- pDT[gene_id %in% GOI.targets.TFs]
  xDT$gene_id <- factor(xDT$gene_id, levels = GOI.targets.TFs)
  xDT <- hierarch.ordering(xDT, "guide", "gene_id", value.var = "estimate")
  w=length(unique(xDT$gene_id)) * 0.05 + 1
  h=length(unique(xDT$guide)) * 0.05 + 0.5
  xDT[abs(estimate) > 5, estimate := pmin(abs(estimate),5) * sign(estimate)]
  ggplot(xDT, aes(y=guide, x=gene_id, color=estimate, size=pmin(5, -log10(q_value)))) + 
    themeNF(rotate=TRUE) +
    geom_point() +
    scale_color_gradient2(name="log2FC", low="blue", high="red") +
    scale_size_continuous(name="-log10(padj)", range = c(0,5)) +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(vjust=1, hjust=0)) +
    ggtitle(tx) +
    facet_grid(. ~ ., space = "free", scales = "free") +
    xlab("Gene expression") + ylab("CRISPR targets")
  ggsaveNF(out("TF_DE_", tx, "_", lnam, "_allGuides.pdf"),w=w,h=h, guides = TRUE)
  
  (dla.nam <- names(dla.factors)[3])
  for(dla.nam in names(dla.factors)){
    dla <- dla.factors[[dla.nam]]
    xDT2 <- xDT[guide %in% dla]
    if(nrow(xDT2) == 0) next
    xDT2$guide <- factor(xDT2$guide, levels=rev(dla))
    h=length(unique(xDT2$guide)) * 0.05 + 0.5
    ggplot(xDT2, aes(y=guide, x=gene_id, color=estimate, size=pmin(5, -log10(q_value)))) + 
      themeNF(rotate=TRUE) +
      geom_point() +
      scale_color_gradient2(name="log2FC", low="blue", high="red") +
      scale_size_continuous(name="-log10(padj)", range = c(0,5)) +
      scale_x_discrete(position = "top") +
      theme(axis.text.x = element_text(vjust=1, hjust=0)) +
      ggtitle(tx) +
      facet_grid(. ~ ., space = "free", scales = "free") +
      xlab("Gene expression") + ylab("CRISPR targets")
    ggsaveNF(out("TF_DE_", tx, "_", lnam, "_", dla.nam,".pdf"),w=w,h=h, guides = TRUE)}
}



# GSEA ------------------------------------------------------

# . in vivo 14d -----------------------------------------------------------
dla <- fread("metadata/FIGS_02_GSEA.txt")
pDT <- gseaDT[tissue == "in.vivo_14d_noClusters"]
pDT.sel <- pDT[grp %in% dla.factors$supp]
pval.cutoff <- 1e-4
i <- 2
xDT <- data.table()
for(i in 1:nrow(dla)){
  ret <- pDT[pathway == dla[i]$pathway & db == dla[i]$db]
  if(nrow(ret) == 0){
    pw <- dla[i]$pathway
    pw <- gsub(" ", ".", pw)
    pw <- gsub("([\\(\\)])", "\\\\\\1", pw)
    ret <- pDT[grepl(pw, pathway, ignore.case = TRUE) & db == dla[i]$db]
  }
  xDT <- rbind(xDT, ret)
}
xDT <- hierarch.ordering(xDT, "pathway", "grp", value.var = "NES", aggregate = TRUE)
xDT <- xDT[grp %in% dla.factors$supp]
xDT$grp <- factor(xDT$grp, levels=dla.factors$supp)
h=length(unique(xDT$pathway)) * 0.07 + 0.5
w=length(unique(xDT$grp)) * 0.07 + 2
xDT$db <- factor(gsub("\\_", "\n", xDT$db), levels = unique(gsub("\\_", "\n", dla$db)))
ggplot(xDT, aes(x=grp, y=pathway, color=NES, size=pmin(5, -log10(padj)))) + 
  themeNF(rotate=TRUE) +
  geom_point() +
  scale_color_gradient2(name="NES", low="blue", high="red") +
  scale_size_continuous(name="-log10(padj)", range = c(0,5)) +
  facet_grid(db ~ ., space = "free", scales = "free") +
  theme(strip.text.y = element_text(angle=0)) +
  ylab("") + xlab("")
ggsaveNF(out("GSEA.pdf"),w=w,h=h, guides = TRUE)



# . apoptosis -------------------------------------------------------------
pDT <- gseaDT[tissue == "in.vivo_14d_noClusters"]
pDT <- pDT[grepl("apop", pathway, ignore.case = TRUE)]
xDT <- pDT
# pval.cutoff <- 1e-4
# i <- 2
# xDT <- data.table()
# for(i in 1:nrow(dla)){
#   ret <- pDT[pathway == dla[i]$pathway & db == dla[i]$db]
#   if(nrow(ret) == 0){
#     pw <- dla[i]$pathway
#     pw <- gsub(" ", ".", pw)
#     pw <- gsub("([\\(\\)])", "\\\\\\1", pw)
#     ret <- pDT[grepl(pw, pathway, ignore.case = TRUE) & db == dla[i]$db]
#   }
#   xDT <- rbind(xDT, ret)
# }
xDT <- hierarch.ordering(xDT, "pathway", "grp", value.var = "NES", aggregate = TRUE)
h=length(unique(xDT$pathway)) * 0.07 + 0.5
w=length(unique(xDT$grp)) * 0.07 + 2
xDT$db <- factor(gsub("\\_", "\n", xDT$db), levels = unique(gsub("\\_", "\n", xDT$db)))
ggplot(xDT, aes(x=grp, y=pathway, color=NES, size=pmin(5, -log10(padj)))) + 
  themeNF(rotate=TRUE) +
  geom_point() +
  scale_color_gradient2(name="NES", low="blue", high="red") +
  scale_size_continuous(name="-log10(padj)", range = c(0,5)) +
  facet_grid(db ~ ., space = "free", scales = "free") +
  theme(strip.text.y = element_text(angle=0)) +
  ylab("") + xlab("")
ggsaveNF(out("GSEA_apoptosis.pdf"),w=5,h=1, guides = TRUE)

# correlation when using clusters -----------------------------------------
pDT <- cDT[time_x == "14d" & time_y == time_x]
pDT <- pDT[tissue_x == "in.vivo" & tissue_x == tissue_y]
pDT <- rbindlist(list(
  "same approach_different CF"=pDT[type_x == type_y & gene_x != gene_y],
  "same CF_different approach"=pDT[type_x != type_y & gene_x == gene_y]
), idcol = "type")
ggplot(pDT, aes(x=gsub("_", "\n", type), y=value)) + 
  geom_violin(fill="lightblue", color="black") +
  geom_boxplot(width=.1, coef=Inf) +
  themeNF() +
  labs(x="", y="Correlation")
ggsaveNF(out("Comparison_withOrWithoutClusters.pdf"),w=2,h=1, guides = TRUE)  


# Examples
uDT <- deDT[tissue %in% c("in.vivo_14d_noClusters", "in.vivo_14d_useClusters")]
uDT[, type := gsub(".+_", "", tissue)]
uDT[, time := "14d"]
uDT[, issue := "in.vivo"]
uDT[, gene := guide]

gg <- c("Kmt2d", "Ash1l", "Wdr82")
gx <- gg[1]
xDT <- rbind()
for(gx in gg){
  xDT <- rbind(xDT, merge(
    uDT[guide == gx & tissue == "in.vivo_14d_noClusters"],
    uDT[guide %in% gg],
    by=c("gene_id"),
    suffixes=c("_x", "_y")
    ))
}

xDT[, same_type := ifelse(type_x == type_y, "same approach", "different approach")]
xDT[, same_guide := ifelse(guide_x == guide_y, "same CF", "different CF")]

cors <- cDT %>%
  filter(gene_x %in% gg & gene_y %in% gg) %>%
  filter(type_x == "noClusters") %>%
  filter(tissue_x == "in.vivo" & tissue_y == "in.vivo") %>%
  filter(time_x == "14d" & time_y == "14d")

rename.type <- function(x){
  x[x == "useClusters"] <- "cluster covariate"
  x[x == "noClusters"] <- "cluster-agnostic"
  factor(x, levels=c("cluster-agnostic", "cluster covariate"))
}

ggplot(xDT) + 
  stat_binhex(aes(fill=log10(..count..), x=estimate_x, y=estimate_y)) +
  facet_grid(gene_x + rename.type(type_x) ~ gene_y + rename.type(type_y), switch="y") +
  themeNF() +
  geom_text(data=cors, aes(label=paste("R =", round(value, 3))), x=-4, y=7, hjust=0,vjust=1) +
  scale_fill_gradientn(name="Nr. genes\n(log10)", colours=c("#a6cee3", "#fdbf6f", "#ff7f00")) + 
  labs(
    x="log fold change",
    y="log fold change"
    )
ggsaveNF(out("Comparison_withOrWithoutClusters_Examples.pdf"),w=3,h=1.5, guides = TRUE)  

