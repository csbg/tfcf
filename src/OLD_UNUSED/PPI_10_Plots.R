source("src/00_init.R")
out <- dirout("PPI_10_01_Plots")

# Read interaction data ---------------------------------------------------
(load(dirout_load("PPI_01_CollectData")("GOI.Interactions.RData")))
interactions[db == "DepMap", dataset := gsub("^(.)", "\\U\\1", gsub("_", " ", dataset), perl=TRUE)]
interactions[,.N, by=c("db", "dataset", "organism")]
table(interactions$db)
interactions[db %in% c("AML", "DepMap"), type := "Genetic interaction"]
interactions[db %in% c("Manual", "Corum"), type := "Complexes"]
interactions[grepl("distance", dataset, ignore.case = TRUE), type := "Indirect PPI"]
interactions[is.na(type), type := "PPI"]
interactions[,.N, by=c("type", "db", "dataset", "organism")][order(type)]


# List of TF-CF pairs -----------------------------------------------------
tfcf <- fread(dirout_load("EXT_01_GetAnnotationsDavid/")("Manual.TFCF.mapping.tsv"))
tfs <- tfcf[Category == "TF"]$GENE
cfs <- tfcf[Category == "CF"]$GENE
tfcf.pairs <- unique(interactions[,c("A","B", "pair")])[(A %in% tfs & B %in% cfs) | (A %in% cfs & B %in% tfs)]$pair


# Our groups --------------------------------------------------------------
our.grps <- fread(dirout_load("FIG_02_POOLED")("SimpleHM.tsv"))
our.grps.pairs <- with(our.grps, split(Gene, Cluster))
our.grps.pairs <- lapply(our.grps.pairs, unique)
our.grps.pairs <- lapply(our.grps.pairs, function(vc){
  res <- data.table()
  for(i in vc){for(j in vc){res <- rbind(res, data.table(A=i, B=j))}}
  res
})
our.grps.pairs <- do.call(rbind, our.grps.pairs)
our.grps.pairs <- with(unique(data.table(t(apply(as.matrix(our.grps.pairs[A != B]), 1, sort)))), paste(V1, "-", V2))


# DepMap analysis ---------------------------------------------------------
depmap.file <- out("DepMap.summaries.tsv")
if(!file.exists(depmap.file)){
  dpm <- toMT(interactions[db == "DepMap"], row = "pair", col = "dataset", val = "Score")
  
  # Get those specific to one tissue (blood)
  dpm.blood <- dpm[,"Blood"]
  dpm.means <- apply(dpm[,!colnames(dpm) %in% c("All", "Blood")], 1, mean)
  dpm.sds <- apply(dpm[,!colnames(dpm) %in% c("All", "Blood")], 1, sd)
  dpm.blood.z <- (dpm.blood-dpm.means)/dpm.sds
  tail(sort(dpm.blood.z))
  plot(density(dpm[1,]))
  
  # Get those broadly regulated
  dpm.p <- apply(dpm[,!colnames(dpm) %in% c("All")], 1, function(row) wilcox.test(row, alternative = "greater")$p.value)
  dpm.padj <- p.adjust(dpm.p, method = "BH")
  dpm.median <- apply(dpm[,!colnames(dpm) %in% c("All")], 1, median)
  
  # export values
  dpm.res <- data.table(pair=row.names(dpm), means=dpm.means, sds=dpm.sds, blood.z=dpm.blood.z, blood=dpm.blood, pvals=dpm.p, padj=dpm.padj, median=dpm.median)
  write.tsv(dpm.res, depmap.file)
} else {
  dpm.res <- fread(depmap.file)
}

# SETUP ENDS HERE ---------------------------------------------------------





# Define first groups - for plot ------------------------------------------
pairs.wide <- list(
  DepMap.blood = dpm.res[blood.z > qnorm(p = 0.99) & blood > 0.2 & means > -0.2]$pair,
  DepMap.broad = dpm.res[padj < 0.05]$pair,
  AML.screen = interactions[db == "AML"][Score > 0.7]$pair,
  Complex.known = unique(interactions[type == "Complexes"]$pair),
  TFCF.pair = tfcf.pairs,
  this.study=our.grps.pairs
)
# Define ppis that are not in complexes
for(ppix in c("Indirect PPI","PPI")){
  pairs.wide[[make.names(ppix)]] <-setdiff(interactions[type %in% c(ppix)][Score > 0.7]$pair, pairs.wide$Complex.known)
}


# Define groups of pairs --------------------------------------------------
pair.small <- list(
  TFCF.pair = pairs.wide$TFCF.pair,
  Complex.known = pairs.wide$Complex.known,
  Complex.extensions = unique(c(pairs.wide$PPI, pairs.wide$Indirect.PPI)),
  Screens = unique(c(pairs.wide$DepMap.blood, pairs.wide$DepMap.broad, pairs.wide$AML.screen)),
  This.study=pairs.wide$this.study
)
stopifnot(length(intersect(pair.small$Complex.extension, pair.small$Complex.known))==0)
pair.groups2 <- pair.small
names(pair.groups2) <- gsub("\\.", "\n", names(pair.groups2))
cleanDev(); pdf(out("Venn.pdf"),w=5,h=5)
gplots::venn(pair.groups2)
dev.off()

all.pairs <- unique(interactions$pair)
fisher.test(as.matrix(table(all.pairs %in% pair.small$Complex.known, all.pairs %in% pair.small$Screens)))



# Define interesting pairs ------------------------------------------------
pairs.interesting <- list(
  Complexes = with(pair.small, intersect(Screens, intersect(Complex.known, This.study))),
  Complex.extensions = with(pair.small, setdiff(Complex.extensions, Complex.known))
)


# Plots pairs -------------------------------------------------------------

pDT <- do.call(rbind, list(
  interactions[pair %in% pair.groups$Complex.extension & pair %in% c(pair.groups$Depmap.global, pair.groups$Depmap.blood)],
  interactions[pair %in% pair.groups$Depmap.blood & pair %in% c(pair.groups$Complex.known)]
))
# pDT <- rbind(
#   pDT[pair %in% pair.groups$TFCF.pair],
#   pDT[!pair %in% pair.groups$TFCF.pair]
# )
pDT <- pDT[pair %in% sample(pDT$pair, 50)]
length(unique(pDT$pair))

height <- length(unique(pDT$pair)) * 0.2 + 2

pDT.ppi <- pDT[db != "DepMap"][dataset != "Core"]
pDT.ppi[db %in% c("Corum", "Manual"), db2 := "Complex"]
pDT.ppi[grepl("_distances", dataset, ignore.case = TRUE), db2 := "Indirect PPI"]
pDT.ppi[grepl("_distances", dataset, ignore.case = TRUE), dataset := gsub("_distances", "", dataset, ignore.case = TRUE)]
pDT.ppi[is.na(db2), db2 := "Direct PPI"]
pDT.ppi$pair <- factor(pDT.ppi$pair, levels=unique(pDT$pair))
p.ppi <- ggplot(pDT.ppi, aes(y=pair, x=paste0(dataset, " (", organism, ") "), color=Score, size=Score)) +
  geom_point() +
  scale_size_continuous(range=c(0,5),limits = c(0,1)) +
  theme_bw(12) +
  facet_grid(. ~ db2, switch = "x", space = "free", scales = "free") +
  scale_color_gradient(low="white", high="#6a3d9a", limits=c(0,1)) +
  scale_y_discrete(drop=FALSE) +
  ylab("")
ggsave(out("Plot.PPI.pdf"), w=6,  h=height, plot=p.ppi + xRot())

pDT.dpm <- pDT[db == "DepMap"]
pDT.dpm$pair <- factor(pDT.dpm$pair, levels=unique(pDT$pair))
p.dpm <- ggplot(pDT.dpm, aes(y=pair, x=paste0(dataset, " (", organism, ") "), color=Score, size=abs(Score))) +
  scale_y_discrete(drop=FALSE) +
  geom_point() +
  scale_size_continuous(range=c(0,5)) +
  scale_color_gradient2(name="Correlation", low="blue", high="red", limits=c(-1,1)) +
  theme_bw(12) +
  facet_grid(. ~ db, switch = "x", space = "free", scales = "free") +
  ylab("")
ggsave(out("Plot.DepMap.pdf"), w=6,  h=height, plot=p.dpm + xRot())

pDT.group <- data.table(pair=pDT$pair, sapply(pair.groups, function(pg) pDT$pair %in% pg) + 0)
pDT.group <- melt(pDT.group, id.vars = "pair")
pDT.group[,db := "Group"]
pDT.group$pair <- factor(pDT.group$pair, levels=unique(pDT$pair))
pDT.group$variable <- factor(gsub("\\.", " ", pDT.group$variable), levels=gsub("\\.", " ", names(pair.groups)))
p.group <- ggplot(pDT.group, aes(y=pair, x=variable, fill=paste0("x", value))) +
  scale_y_discrete(drop=FALSE) +
  geom_tile() +
  scale_fill_manual(values=c(x1="#33a02c", x0="white")) +
  theme_bw(12) +
  facet_grid(. ~ db, switch = "x", space = "free", scales = "free") +
  guides(color=FALSE, fill=FALSE) +
  ylab("")
ggsave(out("Plot.Group.pdf"), w=3,  h=height, plot=p.group + xRot())

p <- gridExtra::grid.arrange(
  p.group +
    theme(
      plot.margin=unit(c(0,0,0,0), units="mm"),
      strip.background.y=element_blank()
      ) + guides(color=FALSE, fill=FALSE, size=FALSE),
  p.ppi +
    theme(
      plot.margin=unit(c(0,0,0,0), units="mm"),
      strip.background.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y = element_blank()
      ) + guides(color=FALSE, fill=FALSE, size=FALSE),
  p.dpm +
    theme(
      plot.margin=unit(c(0,0,0,0), units="mm"),
      strip.background.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y = element_blank()
      ) + guides(color=FALSE, fill=FALSE, size=FALSE),
  nrow=1, ncol=3, widths=c(1.5, 1,3))
ggsave(out("Plot.combined.pdf"), h=height-2, w=8, plot=p)

