source("src/00_init.R")
base.dir <- "FIG_02_scRNA_healthy/"
outBase <- dirout(base.dir)


# FUNCTIONS ---------------------------------------------------------------
ds <- function(path){load(path); return(monocle.obj)}


# Folders -----------------------------------------------------------------
list.files(dirout_load("")(""))
inDir.funcs <- list(
  "in.vivo"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vivo"),
  "in.vitro"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vitro"),
  "leukemia"=dirout_load("FULLINT_10_01_BasicAnalysis_leukemia")
)
inDir.full <- dirout_load("FULLINT_10_01_BasicAnalysis_combined")



# Load data ---------------------------------------------------------------
marker.signatures.use <- fread("metadata/markers.signatures.use.scRNA.tsv")
marker.signatures.use[, sig := paste(DB, Celltype)]

# Load signatures
ff <- list.files(dirout_load("FULLINT_08_01_Markers")(""), pattern="Signatures_")
ff <- dirout_load("FULLINT_08_01_Markers")(ff)
names(ff) <- gsub("^Signatures_(.+?).csv$", "\\1", basename(ff))
ff <- ff[!grepl("Augmented_2021", ff)]
marker.signatures <- lapply(ff, function(fx) as.matrix(read.csv(fx)))
# Prep data
markS <- lapply(unique(marker.signatures.use$DB), function(x) data.table(melt(data.table(marker.signatures[[x]], keep.rownames = TRUE), id.vars = "rn"), DB=x))
markS <- merge(do.call(rbind, markS), marker.signatures.use, by.x=c("variable", "DB"), by.y=c("Celltype", "DB"))
# Cell annotations
annList <- lapply(names(inDir.funcs), function(inDir.current){
  ann <- fread(inDir.funcs[[inDir.current]]("Annotation.tsv"))
  ann[, perturbed := !(mixscape_class.global %in% c("NP", "NTC") | is.na(mixscape_class.global))]
  ann
  })
names(annList) <- names(inDir.funcs)


# Marker signatures ---------------------------------------------------------
inDir.current <- "in.vivo"
inDir.current <- "leukemia"
inDir.current <- "in.vitro"
for(inDir.current in names(inDir.funcs)){
  # out directory
  out <- dirout(paste0(base.dir, "/", inDir.current))
  
  # Prepare data
  ann <- annList[[inDir.current]]
  sda <- fread(inDir.funcs[[inDir.current]]("SigDA.tsv"))
  sda <- merge(sda, unique(ann[,c("timepoint", "markers", "sample")]), c("sample"))
  
  # Prepare plots
  pDT <- merge(ann[perturbed == FALSE][,c("rn", "UMAP1", "UMAP2"),with=F], markS, by="rn")
  pDT[, value.norm := scale(value), by="FinalName"]
  
  # Nomalized signatures split by cell type
  ggplot(pDT, aes(x=UMAP1, y=UMAP2)) +
    stat_summary_hex(aes(z=pmin(3, value.norm)),fun=mean, bins=100) +
    scale_fill_gradient2(low="blue", midpoint = 0, high="red") +
    #scale_fill_hexbin() +
    theme_bw(12) +
    facet_wrap(~FinalName) +
    ggtitle(paste("Marker Signatures"))
  ggsave(out("UMAP_MarkerSignatures.pdf"), w=12,h=12)
  
  # Points (overplotted)
  pDTx <- pDT[value.norm > 2][,paste(FinalName, collapse="xxx"), by=c("rn")][!grepl("xxx", V1)]
  pDTx <- merge(ann, pDTx, by="rn", all=TRUE)
  ggplot(pDTx[!is.na(V1)], aes(x=UMAP1, y=UMAP2, color=V1)) +
    geom_hex(data=pDTx[is.na(V1)], fill="lightgrey", bins=100) +
    geom_point(size=0.2) +
    scale_color_manual(values=COLORS.CELLTYPES.scRNA) +
    theme_bw(12) +
    ggtitle(paste("Marker Signatures"))
  ggsave(out("UMAP_MarkerSignatures_Simple_points.pdf"), w=5,h=4)
  
  # Hexpoints (best?)
  hex.obj <- hexbin::hexbin(x=pDTx$UMAP1, y=pDTx$UMAP2, xbins = 100, IDs=TRUE)
  pDTh <- cbind(pDTx, data.table(hex.x=hex.obj@xcm, hex.y=hex.obj@ycm, hex.cell=hex.obj@cell)[match(hex.obj@cID, hex.cell),])[,.N, by=c("hex.x", "hex.y", "V1")][!is.na(V1)]
  pDTh[, rank := rank(-N), by=c("hex.x", "hex.y")]
  pDTh[rank == 1]
  ggplot(pDTh[rank == 1]) +
    theme_bw(12) +
    geom_hex(data=pDTx, aes(x=UMAP1, y=UMAP2), fill="lightgrey", bins=100) +
    geom_point(aes(x=hex.x, y=hex.y, color=V1), size=0.2) + 
    scale_color_manual(values = COLORS.CELLTYPES.scRNA) +
    ggtitle(paste("Marker Signatures"))
  ggsave(out("UMAP_MarkerSignatures_Simple_hexpoints.pdf"), w=6,h=4)
  
  # Hex but with factor (problematic)
  hex.fun.x <- function(x){
    if(sum(x != "NA") >= 1) names(sort(table(setdiff(x, "NA")), decreasing=TRUE))[1] else NA
  }
  pDTx[, finalLabel := ifelse(is.na(V1), "NA", V1)]
  pDTx$finalLabel <- factor(pDTx$finalLabel)
  ggplot(pDTx, aes(x=UMAP1, y=UMAP2, z=finalLabel)) +
    stat_summary_hex(fun=hex.fun.x, bins=100) +
    theme_bw(12) +
    scale_fill_manual(values = COLORS.CELLTYPES.scRNA) +
    ggtitle(paste("Marker Signatures"))
  ggsave(out("UMAP_MarkerSignatures_Simple.pdf"), w=5,h=4)
  
  
  # Plot enrichment scores
  pDT <- sda[, .(
    d=mean(d), 
    sig.pos=sum(padj < 0.05 & d > 0),
    sig.neg=sum(padj < 0.05 & d < 0),
    total = .N
  ), by=c("guide", "gene", "timepoint", "markers", "sig")]
  pDT <- merge(pDT, marker.signatures.use, by="sig")
  pDT[, sig.perc := ifelse(d > 0, sig.pos, sig.neg)/total * 100]
  pDT <- hierarch.ordering(pDT, toOrder = "gene", orderBy = "FinalName", value.var = "d", aggregate = TRUE)
  pDT <- hierarch.ordering(pDT, toOrder = "FinalName", orderBy = "guide", value.var = "d", aggregate = TRUE)
  
  if(inDir.current == "in.vivo") pDT <- pDT[timepoint == "14d" & markers == "Lin-"]
  ggplot(pDT, aes(x=FinalName, y=guide, size=sig.perc, color=d)) + 
    theme_bw() +
    scale_color_gradient2(low="blue", midpoint = 0, high="red") +
    geom_point() +
    xRot() + 
    facet_grid(gene ~ timepoint + markers, scale="free", space="free") +
    theme(strip.text.y = element_text(angle=0))
  ggsave(out("MarkerSignatures_DA.pdf"), w=12,h=15)
}



# IN VIVO ---------------------------------------------------------
inDir.current <- "in.vivo"
out <- dirout(paste0(base.dir, "/", inDir.current))
ann <- fread(inDir.funcs[[inDir.current]]("Annotation.tsv"))


# . UMAP of timepoint / markers ---------------------------------------------------------
ggplot(ann[perturbed == FALSE][markers != "ckit+"], aes(x=UMAP1, y=UMAP2)) + 
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~markers + timepoint, ncol=3)
ggsave(out("UMAP_Groups.pdf"), w=7,h=5)

# . UMAP of samples -------------------------------------------------------
ggplot(ann[perturbed == FALSE][markers != "ckit+"], aes(x=UMAP1, y=UMAP2)) + 
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~sample, ncol=3)
ggsave(out("UMAP_Samples.pdf"), w=7,h=5)

# . Plot top guides ---------------------------------------------------------
pDT.top <- ann[timepoint == "14d" & markers == "Lin-"]
pDT.ntc <- pDT.top[mixscape_class.global == "NTC"]
pDT.final <- copy(pDT.ntc)
pDT.final$plot <- "NTC"
for(x in c("Wdr82", "Rcor1", "Ehmt1", "Men1", "Kmt2a")){
  pDTg <- rbind(pDT.ntc, pDT.top[grepl(paste0("^", x), guide) & mixscape_class.global == "KO"])
  pDTg$plot <- x
  pDT.final <- rbind(pDT.final, pDTg)
}
pDT.final$plot <- factor(pDT.final$plot, levels=c("NTC", levels(pDT$gene)))
ggplot(pDT.final, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw() +
  geom_hex(data=pDT.final[mixscape_class.global == "NTC"], bins=100, fill="lightgrey") +
  geom_hex(data=pDT.final[mixscape_class.global != "NTC" | plot == "NTC"], bins=100) +
  scale_fill_gradientn(colours=c("#6a3d9a", "#e31a1c", "#ff7f00", "#ffff99")) +
  facet_wrap(~plot, ncol=3)
ggsave(out("UMAP_Guides.pdf"), w=9,h=7)

# . Plot displasia guides ---------------------------------------------------------
pDT.top <- ann[timepoint == "28d"]
pDT.ntc <- pDT.top[mixscape_class.global == "NTC"]
pDT.final <- copy(pDT.ntc)
pDT.final$plot <- "NTC"
for(x in c("Brd9", "Gltscr1", "Phf10")){
  pDTg <- rbind(pDT.ntc, pDT.top[grepl(paste0("^", x), guide) & mixscape_class.global == "KO"])
  pDTg$plot <- x
  pDT.final <- rbind(pDT.final, pDTg)
}
pDT.final$plot <- factor(pDT.final$plot, levels=c("NTC", levels(pDT$gene)))
ggplot(pDT.final, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw() +
  geom_hex(data=pDT.final[mixscape_class.global == "NTC"], bins=100, fill="lightgrey") +
  geom_hex(data=pDT.final[mixscape_class.global != "NTC" | plot == "NTC"], bins=100) +
  scale_fill_gradientn(colours=c("#6a3d9a", "#e31a1c", "#ff7f00", "#ffff99")) +
  facet_wrap(~plot, ncol=3)
ggsave(out("UMAP_Guides_displasia.pdf"), w=9,h=7)

# . Plot displasia numberes -----------------------------------------------
pDT <- copy(ann[markers=="Lin-"])
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
ann <- annList[[inDir.current]]

# . Combine signatures with Antibodies --------------------------------------
pDT <- merge(ann[perturbed == FALSE][,c("rn", "UMAP1", "UMAP2"),with=F], markS, by="rn")
pDT[, value.norm := scale(value), by="FinalName"]
abs <- fread(inDir.funcs[[inDir.current]]("Antibodies.tsv"))
absDT <- merge(
  pDT[,c("rn", "sig", "value.norm"), with=F],
  abs[,c("rn", "Antibody", "Signal.norm"), with=F],
  by="rn", allow.cartesian=TRUE)
absC <- absDT[, cor(value.norm, Signal.norm, use="pairwise.complete.obs"), by=c("sig", "Antibody")][order(V1)]
absC <- hierarch.ordering(absC, toOrder = "sig", orderBy = "Antibody", value.var = "V1")
absC <- hierarch.ordering(absC, orderBy = "sig", toOrder = "Antibody", value.var = "V1")
ggplot(absC,
       aes(x=sig, y=Antibody, fill=V1)) + 
  theme_bw(12) +
  geom_tile() + 
  scale_fill_gradient2(low="blue", high="red") +
  xRot()
ggsave(out("Antibodies_signatures_correlation.pdf"), w=5,h=5)


# . Cell cycle ------------------------------------------------------------
# UMAP
ggplot(ann, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_hexbin() +
  facet_grid(. ~ Phase)
ggsave(out("CellCycle_UMAP.pdf"), w=16,h=5)

# Percentage
pDT <- ann[perturbed == TRUE | mixscape_class == "NTC"][,.N, by=c("Phase", "guide", "timepoint", "markers")]
pDT[,sum := sum(N), by=c("guide", "timepoint", "markers")]
pDT[, perc := N/sum*100]
pDT[, gene := gsub("_.+", "", guide)]
pDT$gene <- factor(pDT$gene, levels=pDT[,mean(perc), by=c("Phase", "gene")][order(V1)][Phase == "G1"]$gene)
ggplot(pDT, aes(x=guide,y=perc,fill=Phase)) + 
  theme_bw(12) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c(G1 = "#b2df8a", G2M = "#6a3d9a", S = "grey")) +
  facet_grid(timepoint + markers ~ gene, scales = "free", space = "free", switch = "x") +
  theme(strip.text.x = element_text(angle=90)) +
  xRot()
ggsave(out("CellCycle_Numbers.pdf"), w=10,h=5)