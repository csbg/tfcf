source("src/00_init.R")
base.dir <- "FIG_02_scRNA_UMAPs/"
outBase <- dirout(base.dir)


# FUNCTIONS ---------------------------------------------------------------
ds <- function(path){load(path); return(monocle.obj)}


# SETTINGS ----------------------------------------------------------------
ff <- list.files(dirout_load("SCRNA_20_Summary")(""))
TISSUES <- gsub("_.+", "", ff[grepl("_monocle$", ff)])
SIGS.USE <- fread("metadata/markers.signatures.use.scRNA.tsv")
SIGS.USE[, sig := paste(DB, Celltype)]

# Folders -----------------------------------------------------------------
inDir.funcs <- list()
for(tx in TISSUES){inDir.funcs[[tx]] <- dirout_load(paste0("SCRNA_20_Summary/", tx, "_monocle.singleR"))}


# Load data ---------------------------------------------------------------

# Annotations
SANN <- fread(PATHS$SCRNA$ANN)

# Marker genes
marker.genes <- fread("metadata/markers.csv")

# Load signatures
marker.signatures <- lapply(TISSUES, function(tissue.name){
  ff <- list.files(dirout_load("SCRNA_06_01_Markers")(tissue.name), pattern="Signatures_", full.names = TRUE)
  names(ff) <- gsub("^Signatures_(.+?).csv$", "\\1", basename(ff))
  ff <- ff[unique(SIGS.USE$DB)]
  lapply(ff, function(fx) as.matrix(read.csv(fx)))
})
names(marker.signatures) <- TISSUES
marker.signatures <- rbindlist(lapply(marker.signatures, function(ltx){
  rbindlist(lapply(ltx, function(dt){
      dt <- melt(data.table(dt, keep.rownames = TRUE), id.vars = "rn")
    }), idcol = "DB")
  }), idcol = "tissue")
marker.signatures <- merge(marker.signatures, SIGS.USE, by.x=c("DB", "variable"), by.y=c("DB", "Celltype"))

# Cell annotations
annList <- lapply(names(inDir.funcs), function(inDir.current){
  ann <- fread(inDir.funcs[[inDir.current]]("Annotation.tsv"))
  ann[, perturbed := !(mixscape_class.global %in% c("NP", "NTC") | is.na(mixscape_class.global))]
  ann[, gene := gsub("_.+", "", guide)]
  ann
  })
annList <- rbindlist(annList, fill=TRUE)

# Other projections
umap.proj <- list(
  original=readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjMonocle.RDS")),
  izzo = readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjIzzo.RDS")),
  in.vivo = readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjVivo.RDS"))
)

# load Monocle Objects
mobjs <- list()
tissuex <- PATHS$SCRNA$MONOCLE.NAMES[1]
for(tissuex in PATHS$SCRNA$MONOCLE.NAMES){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}
for(tissuex in names(mobjs)){
  mobjs[[tissuex]]$CellType <- annList[match(colnames(mobjs[[tissuex]]), rn)]$Clusters
  fData(mobjs[[tissuex]])$gene_short_name <- row.names(fData(mobjs[[tissuex]]))
}



# SETUP ENDS HERE ---------------------------------------------------------




# UMAP Projections ---------------------------------------------------------
tx <- "in.vivo"
tx <- "leukemia"
tx <- "ex.vitro"
for(tx in names(inDir.funcs)){
  # out directory
  out <- dirout(paste0(base.dir, "/", tx))
  
  # Prepare data
  ann <- annList[tissue ==  tx]
  
  # Hexpoints
  for(x in names(umap.proj)){
    xDT <- umap.proj[[x]][match(ann$rn, rn)]
    hex.obj <- hexbin::hexbin(x=xDT$UMAP_1, y=xDT$UMAP_2, xbins = 100, IDs=TRUE)
    pDT <- cbind(ann, data.table(hex.x=hex.obj@xcm, hex.y=hex.obj@ycm, hex.cell=hex.obj@cell)[match(hex.obj@cID, hex.cell),])
    pDT <- pDT[,.N, by=c("hex.x", "hex.y", "Clusters")]
    pDT[, sum := sum(N), by=c("hex.x", "hex.y")]
    pDT[, frac := N / sum]
    pDT <- pDT[frac > 0.25]
    ggplot(pDT, aes(x=hex.x, y=hex.y)) +
      theme_bw(12) +
      #geom_hex(fill="lightgrey", bins=100) +
      geom_point(aes(color=Clusters, shape=Clusters), size=1) + 
      scale_shape_manual(values=rep(c(1,16,2,18,3,4), 20))
    ggsave(out("UMAP_Celltypes_",x,".pdf"), w=6,h=4)
  }
} 


# Signatures --------------------------------------------------------------
tx <- "in.vivo"
tx <- "leukemia"
tx <- "ex.vitro"
for(tx in names(inDir.funcs)){
  # out directory
  out <- dirout(paste0(base.dir, "/", tx))
  
  # annotation
  ann <- annList[tissue ==  tx & timepoint != "28d"]
  
  # Plot signatures scores in each cluster
  pDT <- merge(
    ann[, c("rn", "Clusters"),with=F],
    marker.signatures[, c("rn", "DB", "value", "FinalName"),with=F],
    by="rn")
  pDT <- pDT[, mean(value), by=c("Clusters", "FinalName")]
  pDT[, V1 := scale(V1), by="FinalName"]
  pDT <- hierarch.ordering(pDT, "Clusters", "FinalName", "V1")
  pDT <- hierarch.ordering(pDT, "FinalName", "Clusters", "V1")
  ggplot(pDT, aes(x=FinalName, y=Clusters, fill=V1)) + 
    geom_tile() +
    scale_fill_gradient2(low="blue", high="red") +
    theme_bw(12) +
    xRot()
  ggsave(out("Supp_Signatures.pdf"), w=6,h=5)
  
  # Plot expression of marker genes in each cluster
  
  # Plot fraction of predicted cells in each cluster
}


# Marker gene expression --------------------------------------------------
for(tx in names(mobjs)){
  out <- dirout(paste0(base.dir, "/", tx))
  plot_genes_by_group(mobjs[[tx]], markers = marker.genes$Name, group_cells_by = "CellType") + scale_size_continuous(range=c(0,5))
  ggsave(out("Markers_Clusters.pdf"), w=6,h=8)
}


# Cluster enrichment analyses ---------------------------------------------
tx <- "in.vivo"
for(tx in names(inDir.funcs)){
  # out directory
  out <- dirout(paste0(base.dir, "/", tx))
  
  # Collect enrichment scores
  typex <- "noBcells"
  for(typex in gsub(".+_(.+).pdf", "\\1", list.files(inDir.funcs[[tx]](""), pattern="Guides_Fisher_Mixscape_.*.pdf"))){
    fish.file <- inDir.funcs[[tx]]("Guides_Fisher_Mixscape_",typex,".tsv")
    if(!file.exists(fish.file)) next
    
    fish.full <- fread(fish.file)
    fish.full <- merge(fish.full, unique(SANN[,c("sample_broad", "timepoint"),with=F]), by.x="sample", by.y="sample_broad")
    for(timex in c(unique(fish.full$timepoint), "all")){
      fish <- copy(fish.full)
      if(timex != "all") fish <- fish[timepoint == timex]
      
      # summarize across NTCs
      fish <- fish[, .(
        log2OR=mean(log2OR), 
        dir=length(unique(sign(log2OR)))==1, 
        padj=sum(padj < 0.01), 
        N=.N), by=c("sample", "Clusters", "mixscape_class")]
      fish[dir == FALSE, padj := 0]
      fish[dir == FALSE, log2OR := NA]
      
      # Summarize across guides
      fish[, gene := gsub("_.+", "", mixscape_class)]
      fish <- fish[, .(
        log2OR=mean(log2OR, na.rm=TRUE), 
        dir=length(unique(sign(log2OR[!is.na(log2OR)])))==1, 
        padj=sum(padj), 
        N=sum(N)), by=c("sample", "Clusters", "gene")]
      fish[dir == FALSE, padj := 0]
      fish[dir == FALSE, log2OR := NA]
      
      # summarize across samples
      fish <- fish[, .(
        log2OR=mean(log2OR, na.rm=TRUE), 
        dir=length(unique(sign(log2OR[!is.na(log2OR)])))==1, 
        padj=sum(padj), 
        N=sum(N)), by=c("Clusters", "gene")]
      fish[dir == FALSE, padj := 0]
      fish[dir == FALSE, log2OR := NA]
      
      # setup for plotting
      fish[padj == 0 | is.na(log2OR), log2OR := 0]
      fish[, sig.perc := padj / N]
      fish[,log2OR_cap := pmin(abs(log2OR), 5) * sign(log2OR)]
      fish <- hierarch.ordering(dt = fish, toOrder = "gene", orderBy = "Clusters", value.var = "log2OR")
      fish <- hierarch.ordering(dt = fish, toOrder = "Clusters", orderBy = "gene", value.var = "log2OR")
      ggplot(fish, aes(x=gene, y=Clusters, size=sig.perc, color=log2OR_cap)) + 
        theme_bw() +
        scale_color_gradient2(low="blue", midpoint = 0, high="red") +
        geom_point() +
        geom_point(shape=1, color="lightgrey") +
        xRot()
      ggsave(out("Cluster_enrichments_",typex,"_", timex, ".pdf"), w=length(unique(fish$gene))*0.2 + 2,h=length(unique(fish$Clusters))*0.2+1)
      write.tsv(fish, out("Cluster_enrichments_",typex,"_", timex,".tsv"))
    }
  }
}



# IN VIVO ---------------------------------------------------------
inDir.current <- "in.vivo"
out <- dirout(paste0(base.dir, "/", inDir.current))
ann <- annList[tissue == inDir.current]
fish.bcells <- fread(out("Cluster_enrichments","_basic", "_14d",".tsv"))

# . UMAP of timepoint / markers ---------------------------------------------------------
ggplot(ann[perturbed == FALSE], aes(x=UMAP1, y=UMAP2)) + 
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~markers + timepoint, ncol=3)
ggsave(out("UMAP_Groups.pdf"), w=7,h=5)

# . UMAP of samples -------------------------------------------------------
ggplot(ann[perturbed == FALSE], aes(x=UMAP1, y=UMAP2)) + 
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~sample, ncol=3)
ggsave(out("UMAP_Samples.pdf"), w=7,h=5)

# . Plot top guides ---------------------------------------------------------
pDT.top <- ann[timepoint == "14d" & markers == "lin-"]
pDT.ntc <- pDT.top[mixscape_class.global == "NTC"]
pDT.final <- copy(pDT.ntc)
pDT.final$plot <- "NTC"
for(x in c("Wdr82", "Rcor1", "Ehmt1", "Men1", "Kmt2a")){
  pDTg <- rbind(pDT.ntc, pDT.top[grepl(paste0("^", x), guide) & perturbed == TRUE])
  pDTg$plot <- x
  pDT.final <- rbind(pDT.final, pDTg)
}
pDT.final$plot <- relevel(factor(pDT.final$plot), ref = "NTC")
ggplot(pDT.final, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw() +
  geom_hex(data=pDT.final[mixscape_class.global == "NTC"], bins=100, fill="lightgrey") +
  geom_hex(data=pDT.final[mixscape_class.global != "NTC" | plot == "NTC"], bins=100) +
  scale_fill_gradientn(colours=c("#6a3d9a", "#e31a1c", "#ff7f00", "#ffff99")) +
  facet_wrap(~plot, ncol=3)
ggsave(out("UMAP_Guides.pdf"), w=9,h=7)


# . Plot b cell enrichments -------------------------------------------------
pDT <- fish.bcells[Clusters == "Imm. B-cell"]
pDT$gene <- factor(pDT$gene, levels = pDT[order(log2OR)]$gene)
ggplot(pDT, aes(x=log2OR, y=gene, size=sig.perc, color = log2OR)) + 
  theme_bw(12) +
  geom_vline(xintercept = 0) +
  scale_size_continuous(range = c(1,4)) +
  geom_point() + 
  geom_point(shape=1, color="lightgrey") + 
  scale_color_gradient2(low="blue", high="red") +
  ylab("") + xlab("Enrichment in immature B-cells")
ggsave(out("Cluster_enrichments_BcellsOnly.pdf"), w=5,h=4)

# . NTC distributions to assess clonality effects -------------------------------------------------
pDT <- ann[mixscape_class.global == "NTC", .N, by=c("Clusters", "sample_broad", "CRISPR_Cellranger")]
pDT[, sum := sum(N), by=c("CRISPR_Cellranger", "sample_broad")]
pDT[, frac := N/sum * 100]
ggplot(pDT,  aes(x=Clusters, y=frac, fill=CRISPR_Cellranger)) + 
  theme_bw(12) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~sample_broad) +
  xRot()
ggsave(out("NTC_clonality.pdf"), w=12, h=6)
write.tsv(pDT[, sum(N), by=c("CRISPR_Cellranger", "sample_broad")][order(sample_broad)], out("NTC_clonality.tsv"))


# . Plot displasia guides ---------------------------------------------------------
pDT.top <- ann[timepoint == "28d"][markers=="lin-"]
pDT.ntc <- pDT.top[mixscape_class.global == "NTC"]
pDT.final <- copy(pDT.ntc)
pDT.final$plot <- "NTC"
for(x in c("Brd9", "Gltscr1", "Phf10")){
  pDTg <- rbind(pDT.ntc, pDT.top[grepl(paste0("^", x), guide) & perturbed == TRUE])
  pDTg$plot <- x
  pDT.final <- rbind(pDT.final, pDTg)
}
pDT.final$plot <- relevel(factor(pDT.final$plot), ref = "NTC")
ggplot(pDT.final, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw() +
  geom_hex(data=pDT.final[mixscape_class.global == "NTC"], bins=100, fill="lightgrey") +
  geom_hex(data=pDT.final[mixscape_class.global != "NTC" | plot == "NTC"], bins=100) +
  scale_fill_gradientn(colours=c("#6a3d9a", "#e31a1c", "#ff7f00", "#ffff99")) +
  facet_wrap(~plot, ncol=3)
ggsave(out("UMAP_Guides_displasia.pdf"), w=9,h=7)

# . Plot displasia numbers -----------------------------------------------
pDT <- copy(ann[markers=="lin-"])
pDT[, group := gsub("_.+", "", guide)]
pDT <- pDT[!is.na(group)]
pDT <- pDT[,.N, by=c("timepoint", "mixscape_class.global", "group")]
pDT[, sum := sum(N), by="timepoint"]
pDT[, perc := N/sum*100]
pDT <- pDT[mixscape_class.global != "NP"]
pDT <- pDT[group %in% pDT[, .N, by="group"][N == 2]$group]

# Bar plot
pDT$group <- factor(pDT$group, levels=pDT[timepoint == "28d"][order(N)]$group)
ggplot(pDT, aes(x=group, y=perc, fill=group == "NTC")) + 
  theme_bw() +
  scale_fill_manual(values=c("black", "red")) +
  geom_bar(position="dodge", stat='identity') +
  facet_grid(. ~ timepoint, space="free", scales = "free") + 
  guides(fill=FALSE) +
  ylab("Percent of cells with guide") + 
  xRot()
ggsave(out("Displasia_Numbers_bars.pdf"), w=5,h=3)



# LEUKEMIA---------------------------------------------------------

# . load data ---------------------------------------------------------
inDir.current <- "leukemia"
out <- dirout(paste0(base.dir, "/", inDir.current))
ann <- annList[tissue == inDir.current]
abs <- fread(inDir.funcs[[inDir.current]]("Antibodies.tsv"))

# . Antibodies on UMAP ----------------------------------------------------
ggplot(abs, aes(x=UMAP1, y=UMAP2)) +
  stat_summary_hex(bins = 100, aes(z=pmin(abs(Signal.norm), 2) * sign(Signal.norm)),fun=mean) +
  scale_fill_gradient2(low="blue", high="red") +
  facet_wrap(~Antibody) +
  theme_bw(12)
ggsave(out("Antibodies_UMAP.pdf"), w=12+2, h=9+1)

# . Combine signatures with Antibodies --------------------------------------
pDT <- abs[, .(Signal=mean(Signal.norm, na.rm=TRUE)), by=c("Clusters", "Antibody")]
pDT <- hierarch.ordering(pDT, "Clusters", "Antibody", "Signal")
pDT <- hierarch.ordering(pDT, "Antibody", "Clusters", "Signal")
ggplot(pDT, aes(x=Clusters,y=Antibody, fill=Signal)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", high="red") +
  xRot()
ggsave(out("Antibodies_Average.pdf"), w=5, h=4)

# . Cell cycle ------------------------------------------------------------
# UMAP
ggplot(ann, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_hexbin() +
  facet_grid(. ~ Phase)
ggsave(out("CellCycle_UMAP.pdf"), w=16,h=5)

# Percentage
pDT <- ann[perturbed == TRUE | mixscape_class == "NTC"][,.N, by=c("Phase", "CRISPR_Cellranger", "timepoint", "markers")]
pDT[,sum := sum(N), by=c("CRISPR_Cellranger", "timepoint", "markers")]
pDT[, perc := N/sum*100]
pDT[, gene := gsub("_.+", "", CRISPR_Cellranger)]
pDT$gene <- factor(pDT$gene, levels=pDT[,mean(perc), by=c("Phase", "gene")][order(V1)][Phase == "G1"]$gene)
ggplot(pDT, aes(x=CRISPR_Cellranger,y=perc,fill=Phase)) + 
  theme_bw(12) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c(G1 = "#b2df8a", G2M = "#6a3d9a", S = "grey")) +
  facet_grid(timepoint + markers ~ gene, scales = "free", space = "free", switch = "x") +
  theme(strip.text.x = element_text(angle=90)) +
  xRot()
ggsave(out("CellCycle_Numbers.pdf"), w=15,h=5)


# . Plot top guides ---------------------------------------------------------
pDT.top <- copy(ann)
pDT.ntc <- pDT.top[mixscape_class.global == "NTC"]
pDT.final <- copy(pDT.ntc)
pDT.final$plot <- "NTC"
for(x in unique(pDT.top[perturbed == TRUE]$gene)){
  pDTg <- rbind(pDT.ntc, pDT.top[grepl(paste0("^", x), guide) & mixscape_class.global == "KO"])
  pDTg$plot <- x
  pDT.final <- rbind(pDT.final, pDTg)
}
pDT.final$plot <- relevel(factor(pDT.final$plot),ref = "NTC")
ggplot(pDT.final, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw() +
  geom_hex(data=pDT.final[mixscape_class.global == "NTC"], bins=100, fill="lightgrey") +
  geom_hex(data=pDT.final[mixscape_class.global != "NTC" | plot == "NTC"], bins=100) +
  scale_fill_gradientn(colours=c("#6a3d9a", "#e31a1c", "#ff7f00", "#ffff99")) +
  facet_wrap(~plot, ncol=3)
ggsave(out("UMAP_Guides.pdf"), w=9,h=15)



# Cancer vs Healthy -------------------------------------------------------
out <- dirout(paste0(base.dir, "/", "CvH"))

pDT <- list(
  leukemia = fread(outBase("leukemia/Cluster_enrichments_basic_all.tsv")),
  normal = fread(outBase("ex.vivo/Cluster_enrichments_basic_all.tsv"))
)
pDT <- rbindlist(pDT, idcol = "tissue")

ggplot(pDT, aes(x=Clusters, y=tissue, size=sig.perc, color=log2OR_cap)) + 
  theme_bw() +
  scale_color_gradient2(low="blue", midpoint = 0, high="red") +
  geom_point() +
  geom_point(shape=1, color="lightgrey") +
  facet_grid(gene ~ .) +
  xRot() +
  theme(strip.text.y = element_text(angle=0))
ggsave(out("Cluster_enrichments.pdf"), w=4,h=16)
