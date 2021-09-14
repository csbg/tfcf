source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("PPI_10_02_Plots")

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
pairsDT <- data.table(pair=unique(interactions$pair))
pairsDT[,this.study := pair %in% our.grps.pairs]
pairsDT[,TFCF.pair := pair %in% tfcf.pairs]
pairsDT[,Complex.known := pair %in% unique(interactions[type == "Complexes"]$pair)]
for(ppix in c("Indirect PPI","PPI")){
  pairsDT[[make.names(ppix)]] <- pairsDT$pair %in% setdiff(interactions[type %in% c(ppix)][Score > 0.7]$pair, pairsDT[Complex.known == TRUE]$pair)
}
pairsDT[, Complex.extensions := PPI | Indirect.PPI]
pairsDT[,AML.screen := pair %in% interactions[db == "AML"][Score > 0.7]$pair]
pairsDT[,DepMap.blood := pair %in% dpm.res[blood.z > qnorm(p = 0.99) & blood > 0.2 & means > -0.2]$pair]
pairsDT[,DepMap.broad := pair %in% dpm.res[padj < 0.05]$pair]
pairsDT[, Screens := DepMap.blood | DepMap.broad | AML.screen]


fisher.test(as.matrix(with(pairsDT, table(Complex.known, Screens))))


# Define interesting pairs ------------------------------------------------
pairs.interesting <- rbind(
  pairsDT[Complex.known & this.study & Screens],
  pairsDT[this.study & !Complex.known & (Complex.extensions | Screens)],
  pairsDT[!Complex.known & TFCF.pair & Complex.extensions & Screens]
)


# PLOT --------------------------------------------------------------------
pairs.interesting <- pairs.interesting[order(this.study, Complex.known, TFCF.pair, Complex.extensions, AML.screen, DepMap.blood, DepMap.broad)]
pDT <- interactions[pair %in% pairs.interesting$pair]
height <- length(unique(pDT$pair)) * 0.2 + 2

pDT.group <- pairs.interesting[,-c("Screens")]
for(cx in colnames(pDT.group)){
  if(is.logical(pDT.group[[cx]])){pDT.group[[cx]] <- pDT.group[[cx]] + 0}
}
pDT.group <- melt(pDT.group, id.vars = "pair")
pDT.group[,db := "Group"]
pDT.group$pair <- factor(pDT.group$pair, levels=unique(pairs.interesting$pair))
str(o <- intersect(colnames(pairs.interesting), pDT.group$variable))
clean.grps <- function(x){gsub("\\.", " ", x)}
pDT.group$variable <- factor(clean.grps(pDT.group$variable), levels=clean.grps(o))
p.group <- ggplot(pDT.group, aes(y=pair, x=variable, fill=paste0("x", value))) +
  scale_y_discrete(drop=FALSE) +
  geom_tile() +
  scale_fill_manual(values=c(x1="#33a02c", x0="white")) +
  theme_bw(12) +
  facet_grid(. ~ db, switch = "x", space = "free", scales = "free") +
  guides(color=FALSE, fill=FALSE) +
  ylab("")
ggsave(out("Plot.Group.pdf"), w=3,  h=height, plot=p.group + xRot())

pDT.ppi <- pDT[type %in% c("PPI", "Indirect PPI", "Complexes")][db != "Corum" | grepl("Core", dataset)]
pDT.ppi[, type := gsub("Indirect PPI", "PPI (ind)", type)]
pDT.ppi$pair <- factor(pDT.ppi$pair, levels=unique(pairs.interesting$pair))
p.ppi <- ggplot(pDT.ppi, aes(y=pair, x=paste0(dataset, " (", organism, ") "), color=Score, size=Score)) +
  geom_point() +
  scale_size_continuous(range=c(0,5),limits = c(0,1)) +
  theme_bw(12) +
  facet_grid(. ~ type, switch = "x", space = "free", scales = "free") +
  scale_color_gradient(low="white", high="#6a3d9a", limits=c(0,1)) +
  scale_y_discrete(drop=FALSE) +
  ylab("")
ggsave(out("Plot.PPI.pdf"), w=6,  h=height, plot=p.ppi + xRot())

pDT.dpm <- pDT[db == "DepMap" | db == "AML"]
pDT.dpm$pair <- factor(pDT.dpm$pair, levels=unique(pairs.interesting$pair))
p.dpm <- ggplot(pDT.dpm, aes(y=pair, x=paste0(dataset, " (", organism, ") "), color=Score, size=abs(Score))) +
  scale_y_discrete(drop=FALSE) +
  geom_point() +
  scale_size_continuous(range=c(0,5)) +
  scale_color_gradient2(name="Correlation", low="blue", high="red", limits=c(-1,1)) +
  theme_bw(12) +
  facet_grid(. ~ db, switch = "x", space = "free", scales = "free") +
  ylab("")
ggsave(out("Plot.DepMap.pdf"), w=7,  h=height, plot=p.dpm + xRot())

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
  nrow=1, ncol=3, widths=c(3, 2,4))
ggsave(out("Plot.combined.pdf"), h=height-2, w=9, plot=p)

