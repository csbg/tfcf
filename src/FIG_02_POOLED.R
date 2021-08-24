source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("FIG_02_POOLED/")

require(latex2exp)
require(ggrepel)
require(igraph)
require(ggtext)

# Load data ---------------------------------------------------------------
m <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix))
ann <- fread(PATHS$POOLED$DATA$annotation)
stopifnot(all(ann$sample == colnames(m)))


# Load analysis results ---------------------------------------------------
inDir <- dirout_load("POOLED_10_03_IndividualAnalysis_NormFactors_Controls")
list.files(inDir(""), pattern=".tsv")
RESULTS.wt <- fread(inDir("Results_Pvalues.tsv"))[Library != "A"]
RESULTS.wt.agg <- fread(inDir("Aggregated.tsv"))
RESULTS.wt.mds <- fread(inDir("Correlation_hits_MDS.tsv"))


# Hopefully Cool plot ---------------------------------------------------------------
# Add Main branch
cleanComparisons2 <- function(x){
  orderX <- c("Main branch", "LSK.CKIT", "GMP.MEP", "MYE.UND", "GMPcd11.DN")
  x <- factor(x, levels=orderX)
}

# Comparison groups
compDT <- unique(RESULTS.wt[Library != "A"][,c("Population1", "Population2", "Comparison"), with=F])
compDT[Comparison %in% c("GMP.LSK", "MEP.LSK","UND.MEP", "MYE.GMP"), Comparison.Group := "Main branch"]
compDT[is.na(Comparison.Group), Comparison.Group := Comparison]
compDT2 <- unique(melt(compDT[,-"Comparison",with=F], id.vars = "Comparison.Group")[,-"variable"])

# Aggregated data within population
agDT <- copy(RESULTS.wt.agg)
mainBranchOrdering <- c("Und", "MEP", "LSKd7", "GMP", "Mye")
agDT$Population <- factor(agDT$variable, levels=c(mainBranchOrdering, setdiff(agDT$variable, mainBranchOrdering)))
agDT <- merge(agDT, compDT2, by.y="value", by.x="variable", allow.cartesian=TRUE)
agDT$rn <- factor(agDT$rn, levels=RESULTS.wt.mds$rn[hclust(dist(as.matrix(data.frame(RESULTS.wt.mds[,c("V1", "V2"), with=F]))))$order])

# Comparisons between populations
pDT.stats <- RESULTS.wt[hit == TRUE][Library != "A"]
pDT.stats <- merge(pDT.stats, unique(compDT[,c("Comparison", "Comparison.Group")]), by="Comparison")
# Filter 2:  50 % of guides being significant
pDT.stats <- pDT.stats[, .(mean(z), length(unique(Guide[padj < 0.05])),n=length(unique(Guide))), by=c("Gene", "Comparison",  "Comparison.Group", "Population1", "Population2")]
pDT.stats[,percSig := V2/n*100]
pDT.stats <- pDT.stats[percSig > 50]

# Plot
ggplot(agDT[rn %in% pDT.stats$Gene], aes(x=Population, y=rn)) +
  theme_bw(12) +
  geom_point(aes(fill=log2FC), shape=21, color="white", size=5) +
  facet_grid(. ~ cleanComparisons2(Comparison.Group), scales = "free", space = "free") +
  geom_segment(data=pDT.stats, aes(xend=Population1, x=Population2, y=Gene, yend=Gene, color=V1), arrow=arrow(type="closed", length = unit(0.3, "cm"))) +
  scale_fill_gradient2(name=TeX(r'($\\overset{\Delta_{Cas9-WT}}{(dots)}$)')) +
  #geom_point(aes(fill=log2FC), shape=21, color="white", size=2) +
  scale_color_gradient2(name=TeX(r'($\\overset{\Delta_{Populations}}{(arrows)}$)')) +
  xRot()
ggsave(out("Aggregated_Edges.pdf"), w=8,h=15)



# Selected comparisons (David) --------------------------------------------
pDT.stats <- RESULTS.wt[Gene %in% RESULTS.wt[hit == TRUE]$Gene]
pDT.stats <- pDT.stats[, .(mean(z), length(unique(Guide[padj < 0.05])),n=length(unique(Guide))), by=c("Gene", "Comparison", "Population1", "Population2")]
pDT.stats[,percSig := V2/n*100]
pDT.stats <- pDT.stats[Gene %in% pDT.stats[percSig > 50]$Gene]
pDT.stats <- pDT.stats[Comparison %in% COMPARISONS.USE]

mx <- toMT(pDT.stats, row = "Gene", col = "Comparison", val = "V1")
mx[is.na(mx)] <- 0
pDT.stats$Cluster <- cutree(hclust(dist(mx)), k = 7)[pDT.stats$Gene]

ggplot(pDT.stats, aes(y=cleanComparisons(Comparison, ggtext = TRUE), x=Gene)) +
  theme_bw(12) + 
  geom_point(aes(color=V1), size=4) +
  scale_color_gradient2(name=TeX(r'($\\overset{\Delta_{Cas9-WT}}$)')) + #, low="#e31a1c", high="#1f78b4") +
  facet_grid(. ~ Cluster, scales = "free", space = "free") +
  xRot() + 
  ylab("") +
  theme(axis.text.y = element_markdown())
ggsave(out("SimpleHM.pdf"), w=15,h=2.5)


# Vulcano plots -----------------------------------------------------------
ggplot(RESULTS.wt[Comparison %in% COMPARISONS.USE], aes(x=z, y=-log10(p))) + 
  geom_hex() +
  facet_grid(. ~ cleanComparisons(Comparison)) +
  theme_bw(12)
ggsave(out("Vulcano.pdf"), w=15,h=2)



# Comparison of two main populations --------------------------------------
cx <- c("GMP.MEP", "LSK.CKIT")
pDT.full <- RESULTS.wt[Comparison %in% cx][Genotype == "Cas9"]
pDT.sig <- pDT.full[hit == TRUE][,paste(sort(unique(cleanComparisons(Comparison))), collapse = ","), by="Gene"]
pDT.sig[grepl(",", V1),V1 := "Both"]
pDT.long <- pDT.full[,.(mean(z)), by=c("Gene", "Comparison", "Genotype")]
pDT <- dcast.data.table(pDT.long, Gene + Genotype ~ Comparison, value.var = "V1")
pDT <- merge(pDT, pDT.sig, by="Gene", all.x=TRUE)
pDT[, sig := V1]
pDT[is.na(sig), sig := "None"]
pDT$sig <- factor(pDT$sig, levels=c("Both", cleanComparisons(cx, order = FALSE), "None"))

xx <- 4

formatArrows <- function(x){
  paste0(gsub("^(.+?)\\.(.+)$", "\\1", x),"  ", r'($\leftarrow$)',"  .  ",r'($\rightarrow$)',"  ", gsub("^(.+?)\\.(.+)$", "\\2", x))
}

dlim <- max(abs(c(pDT[[cx[1]]], pDT[[cx[2]]])))
ggplot(pDT, aes_string(x=cx[1], y=cx[2], color="sig")) + 
  #geom_rect(xmin=-xx, ymin=-xx, ymax=xx, xmax=xx, color=NA, fill="#eeeeee60") +
  geom_hline(yintercept = 0, color="lightgrey", alpha=0.5) +
  geom_vline(xintercept = 0, color="lightgrey", alpha=0.5) +
  # geom_hex() + 
  # scale_fill_gradient(low="lightgrey", high="lightblue") + 
  # geom_point(data=pDT[Gene %in% pDT.full[hit==TRUE]$Gene], shape=1) +
  geom_text_repel(data=pDT[Gene %in% intersect(pDT.full[hit==TRUE]$Gene, pDT.long[abs(V1) > xx]$Gene)], aes(label=Gene)) +
  geom_point(alpha=0.5) +
  scale_color_manual(name="Significant in", values=c("#e31a1c", "#1f78b4", "#6a3d9a", "grey")) +
  #facet_grid(. ~ Genotype) +
  theme_bw() +
  xlab(TeX(formatArrows(cx[1]))) +
  ylab(TeX(formatArrows(cx[2]))) +
  ylim(-dlim, dlim) + xlim(-dlim,dlim)
ggsave(out("Comparison_Scatter.pdf"), w=5.5,h=4.5)





# Plot network ------------------------------------------------------------
outGraphs <- dirout("FIG_02_POOLED/Graphs/")
genex <- "Kmt2d"
for(genex in unique(RESULTS.wt[hit == TRUE]$Gene)){
  COLORS.graph <- c("#fb9a99", "lightgrey", "#a6cee3")
  
  agx <- RESULTS.wt.agg[rn == genex]
  agx <- unique(agx[,c("variable", "log2FC"), with=F])
  
  statx <- RESULTS.wt[Gene == genex][Genotype == "Cas9"]
  statx[, keep := sum(padj < 0.05) >= 2 & (all(sign(Score[padj < 0.05]) > 0) | all(sign(Score[padj < 0.05]) < 0)), by=c("Gene", "Analysis")]
  statx <- statx[, .(
    z=mean(z), 
    up=length(unique(Guide[padj < 0.05 & z > 0])),
    dn=length(unique(Guide[padj < 0.05 & z < 0])),
    n=length(unique(Guide))), 
    by=c("Gene", "Comparison")]
  
  
  el <- data.table(do.call(rbind, COMPARISONS), keep.rownames = TRUE)
  #statx <- statx[match(row.names(el), Comparison)][!is.na(Gene)]
  statx <- merge(statx, el, by.x="Comparison", by.y="rn")
  statx[Comparison == "LSK.CKIT", V1 := "LSKd7"]
  statx[Comparison == "LSK.CKIT", V2 := "LSKd7"]
  # statx[Comparison == "GMPcd11.DN", V2 := "MEP"]
  # statx[Comparison == "GMPcd11.DN", V1 := "Und"]
  # statx  <- statx[Comparison != "UND.MEP"]
  
  g <- graph.edgelist(as.matrix(statx[,c("V2", "V1")]))
  #V(g)$log2FC <- agx[match(V(g)$name, variable),]$log2FC
  E(g)$z <- statx$z
  E(g)$up <- statx$up
  E(g)$dn <- statx$dn
  E(g)$n <- statx$n
  E(g)$sig.label <- with(statx, paste(up, dn, n, sep="/"))
  E(g)$cnt <- ifelse(E(g)$z > 0, E(g)$up, E(g)$dn)
  E(g)$perc <- round(E(g)$cnt/E(g)$n) * 100
  
  g <- delete.vertices(g, is.na(V(g)$log2FC))
  #V(g)$color <- mapNumericToColors(V(g)$log2FC, cols = COLORS.graph)
  V(g)$color <- "lightgrey"
  V(g)$frame.color <- NA
  V(g)$label.color <- "black"
  E(g)$color <- mapNumericToColors(E(g)$z, cols = COLORS.graph)
  E(g)$width <- 2+E(g)$perc/100*3
  E(g)$arrow.width <- 2
  E(g)$arrow.size <- 1
  #E(g)$label <- paste0("", round(E(g)$z, 1),"\n", "(",E(g)$sig.label, ")")
  #E(g)$label.color <- "black"
  
  
  layout <- list(
    # "cKit" = c(3,3),
    # "LSKd9" = c(3,2),
    "GMP.CD11bGr1" = c(3,1),
    "GMP.DN" = c(3,0),
    "LSKd7" = c(0,3),
    "GMP" = c(-1,1.5),
    "MEP" = c(1,1.5),
    "Mye" = c(-1,0),
    "Und" = c(1,0)
    )
  
  cleanDev(); pdf(outGraphs("Graph_", genex, ".pdf"), w=5,h=5)
  plot.igraph(g, layout=do.call(rbind, layout)[V(g)$name,], main=genex)
  dev.off()
}




# Heatmaps for graphs -----------------------------------------------------
res.stats <- RESULTS.wt[Gene %in% RESULTS.wt[hit == TRUE & !grepl("NonTargeting", Gene)]$Gene]

statx <- res.stats[, .(
  z=mean(z), 
  up=length(unique(Guide[padj < 0.05 & z > 0])),
  dn=length(unique(Guide[padj < 0.05 & z < 0])),
  n=length(unique(Guide))), 
  by=c("Gene", "Comparison")]

bg <- data.table()
for(gx in unique(statx$Gene)){
  for(cx in unique(statx$Comparison)){
    bg <- rbind(bg, data.table(Gene = gx, Comparison = cx, z=-5:5))
  }
}


ggplot(statx, aes(x=Gene)) + 
  facet_grid(cleanComparisons(Comparison) ~ .) +
  geom_tile(data=bg, aes(fill=z, y=factor(z))) + 
  scale_fill_gradient2(high="#a6cee3", low="#fb9a99") + 
  geom_point(aes(y=factor(round(z))), size=3, color="black", shape=18) +
  theme_bw(12) + 
  xRot()
ggsave(out("GraphHM.pdf"),w=17, h=10)



# TODO --> example of how p-values were calculated