source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("PPI_02_FirstPlots")

(load(dirout_load("PPI_01_CollectData")("GOI.Interactions.RData")))


interactions[db == "DepMap", dataset := gsub("^(.)", "\\U\\1", gsub("_", " ", dataset), perl=TRUE)]
interactions[,.N, by=c("db", "dataset", "organism")]



# DepMap analysis ---------------------------------------------------------
dpm <- toMT(interactions[db == "DepMap"], row = "pair", col = "dataset", val = "Score")

# Get those specific to one tissue (blood)
dpm.means <- apply(dpm[,!colnames(dpm) %in% c("All", "Blood")], 1, mean)
dpm.sds <- apply(dpm[,!colnames(dpm) %in% c("All", "Blood")], 1, sd)
dpm.blood.z <- (dpm[,"Blood"]-dpm.means)/dpm.sds
tail(sort(dpm.blood.z))
plot(density(dpm[1,]))

# Get those broadly regulated
dpm.p <- apply(dpm[,!colnames(dpm) %in% c("All")], 1, function(row) wilcox.test(row, alternative = "greater")$p.value)
dpm.padj <- p.adjust(dpm.p, method = "BH")
dpm.median <- apply(dpm[,!colnames(dpm) %in% c("All")], 1, median)
#plot(dpm.median, -log10(dpm.p))
dpm.res <- data.table(pair=row.names(dpm), means=dpm.means, sds=dpm.sds, blood.z=dpm.blood.z, pvals=dpm.p, padj=dpm.padj, median=dpm.median)
write.tsv(dpm.res, out("DepMap.summaries.tsv"))

# Plot vulcano plot
ggplot(dpm.res, aes(x=median, y=-log10(pvals))) + 
  geom_hex() +
  theme_bw(12)
ggsave(out("DepMap_Vulcano.pdf"), w=5,h=4)


# Define groups of pairs --------------------------------------------------
pair.groups <- list(
  Complex.extension = unique(interactions[!db %in% c("Corum", "DepMap")][!pair %in% unique(interactions[db %in% c("Corum")]$pair)][Score > 0.7]$pair),
  Complex.known = unique(interactions[db %in% c("Corum")]$pair),
  Depmap.global = dpm.res[padj < 0.05]$pair,
  Depmap.blood = names(which(dpm.blood.z > qnorm(p = 0.99)))
)
pair.groups2 <- pair.groups
names(pair.groups2) <- gsub("\\.", "\n", names(pair.groups2))
cleanDev(); pdf(out("Venn.pdf"),w=5,h=5)
gplots::venn(pair.groups2)
dev.off()



# # Plots pairs -------------------------------------------------------------
pair.groups
pDT <- do.call(rbind, list(
  interactions[pair %in% pair.groups$complex.extension & pair %in% c(pair.groups$depmap.global, pair.groups$depmap.blood)],
  interactions[pair %in% pair.groups$depmap.blood & pair %in% c(pair.groups$complex.known)]
))
length(unique(pDT$pair))

pDT.ppi <- pDT[db != "DepMap"][dataset != "Core"]
pDT.ppi[,db2 := ifelse(db == "Corum", "Complex", "Other")]
pDT.ppi$pair <- factor(pDT.ppi$pair, levels=unique(pDT$pair))
p.ppi <- ggplot(pDT.ppi, aes(y=pair, x=paste0(dataset, " (", organism, ") "), color=Score, size=Score)) +
  geom_point() +
  scale_size_continuous(range=c(0,5),limits = c(0,1)) +
  theme_bw(12) +
  facet_grid(. ~ db2, switch = "x", space = "free", scales = "free") +
  scale_color_gradient(low="white", high="#6a3d9a", limits=c(0,1)) +
  scale_y_discrete(drop=FALSE) +
  ylab("")
ggsave(out("Plot.PPI.pdf"), w=6,  h=15, plot=p.ppi + xRot())

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
ggsave(out("Plot.DepMap.pdf"), w=6,  h=15, plot=p.dpm + xRot())


pDT.group <- data.table(pair=pDT$pair, sapply(pair.groups, function(pg) pDT$pair %in% pg) + 0)
pDT.group <- melt(pDT.group, id.vars = "pair")
pDT.group[,db := "Group"]
pDT.group$pair <- factor(pDT.group$pair, levels=unique(pDT$pair))
p.group <- ggplot(pDT.group, aes(y=pair, x=gsub("\\.", " ", variable), fill=paste0("x", value))) +
  scale_y_discrete(drop=FALSE) +
  geom_tile() +
  scale_fill_manual(values=c(x1="#33a02c", x0="white")) +
  theme_bw(12) +
  facet_grid(. ~ db, switch = "x", space = "free", scales = "free") +
  guides(color=FALSE, fill=FALSE) +
  ylab("")
ggsave(out("Plot.Group.pdf"), w=3,  h=15, plot=p.group + xRot())

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
ggsave(out("Plot.combined.pdf"), h=10, w=8, plot=p)

