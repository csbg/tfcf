source("src/00_init.R")
base.dir <- "FIG_01_POOLED_vsWT/"
out <- dirout(base.dir)

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
RESULTS.wt <- fread(inDir("Results_Pvalues.tsv"))[Library != "A"][Comparison %in% COMPARISONS.healthy]
RESULTS.wt.agg <- fread(inDir("Aggregated.tsv"))


# Gene annotations --------------------------------------------------------
ANN.genes <- fread("metadata/TFCF_Annotations.tsv")
ANN.genes[,Complex_simple := make.names(COMPLEX)]
ANN.genes[,Complex_simple := gsub("\\..*$", "", Complex_simple)]
ANN.genes[Complex_simple == "X", Complex_simple := NA]
ANN.genes[,rank := rank(Complex_simple), by="GENE"]
ANN.genes <- ANN.genes[rank == 1][,-"rank"]
stopifnot(all(ANN.genes[,.N, by="GENE"]$N == 1))
unique(ANN.genes$Complex_simple)


# Define comparison groups ------------------------------------------------
# Add Main branch
cleanComparisons2 <- function(x){
  orderX <- c("Main branch", "LSK.CKIT", "GMP.MEP", "MYE.UND", "GMPcd11.DN")
  x <- factor(x, levels=orderX)
}
# Define groups
compDT <- unique(RESULTS.wt[Library != "A"][,c("Population1", "Population2", "Comparison"), with=F])
compDT[Comparison %in% c("GMP.LSK", "MEP.LSK","UND.MEP", "MYE.GMP"), Comparison.Group := "Main branch"]
compDT[is.na(Comparison.Group), Comparison.Group := Comparison]
FIGS.COMPARISONS <- copy(compDT)


# Summarize results of differential analysis for each gene ----------------
pDT.stats <- merge(RESULTS.wt[Genotype == "Cas9"], unique(FIGS.COMPARISONS[,c("Comparison", "Comparison.Group")]), by="Comparison", all.x=TRUE)
pDT.stats <- pDT.stats[, .(
  z=mean(z), 
  sig.up=length(unique(Guide[padj < 0.05 & z > 0])),
  sig.down=length(unique(Guide[padj < 0.05 & z < 0])),
  n=length(unique(Guide))
  ), by=c("Gene", "Comparison",  "Comparison.Group", "Population1", "Population2")]
pDT.stats[, sig := ifelse(z > 0, sig.up, sig.down)]
pDT.stats[,percSig := sig/n*100]
stopifnot(all(pDT.stats$percSig >= 0))
stopifnot(all(pDT.stats$percSig <= 100))
pDT.stats[, hit := percSig >= 50]
pDT.stats[, complete.screen := all(COMPARISONS.healthy %in% Comparison), by="Gene"]
RESULTS.wt.agg.gene <- copy(pDT.stats)

# Did this keep all genes?
stopifnot(length(unique(RESULTS.wt.agg.gene$Gene)) == length(unique(RESULTS.wt$Gene)))


# Example calculation -----------------------------------------------------
SCORES.EXAMPLE.comp <- "GMP.MEP"
SCORES.EXAMPLE.lib <- "As"
SCORES.EXAMPLE <- fread(inDir("Results_Scores.tsv"))
SCORES.EXAMPLE <- SCORES.EXAMPLE[Comparison == SCORES.EXAMPLE.comp][Library == SCORES.EXAMPLE.lib]
SCORES.EXAMPLE.bg <- fread(inDir("Results_SummaryStats.tsv"))[Comparison == SCORES.EXAMPLE.comp][Library == SCORES.EXAMPLE.lib]



# Data from replicates ----------------------------------------------------
replicationDir <- dirout_load("POOLED_11_03_Replicates_NormFactorsFromControls")

list.files(replicationDir(""))
# Scores
REPLICATES.SCORES <- fread(replicationDir("Results.tsv"))
REPLICATES.SCORES <- REPLICATES.SCORES[grepl("vs", coef)]
REPLICATES.SCORES[, Comparison := gsub("vs", ".", coef)]
# unique(REPLICATES.SCORES[!toupper(Comparison) %in% toupper(names(COMPARISONS))]$Comparison)
REPLICATES.SCORES.EXAMPLE <- REPLICATES.SCORES[coef == gsub("\\.", "vs", SCORES.EXAMPLE.comp)][grepl(SCORES.EXAMPLE.lib, analysis)]

# Values
REPLICATES.VALUES <- fread(replicationDir("analysis_", REPLICATES.SCORES.EXAMPLE$analysis[1], "/", "Data_DateRemoved.csv"))
rns <- REPLICATES.VALUES$rn
REPLICATES.VALUES <- as.matrix(REPLICATES.VALUES[, -"rn"])
row.names(REPLICATES.VALUES) <- rns
REPLICATES.VALUES <- REPLICATES.VALUES[unique(REPLICATES.SCORES.EXAMPLE[significant_TF == TRUE]$rn),ann[Library == SCORES.EXAMPLE.lib][Population %in% COMPARISONS[[SCORES.EXAMPLE.comp]]]$sample]
REPLICATES.VALUES <- melt(data.table(REPLICATES.VALUES, keep.rownames = TRUE), id.vars = "rn")
REPLICATES.VALUES <- merge(ann, REPLICATES.VALUES, by.x="sample", by.y="variable")
REPLICATES.VALUES[, Gene := gsub("_.+", "", rn)]
REPLICATES.VALUES[, sample_i := rank(sample), by=c("Gene", "rn", "Genotype")]
REPLICATES.VALUES[, guide_i := rank(rn), by=c("Gene", "sample")]


# SETUP ENDS HERE ---------------------------------------------------------



# Check inconsistencies ---------------------------------------------------
pDT <- merge(RESULTS.wt[Genotype == "Cas9"], RESULTS.wt.agg.gene[sig.down != 0 & sig.up != 0][,c("Gene", "Comparison"),with=F],by=c("Gene", "Comparison"))
ggplot(pDT, aes(x=Comparison, y=paste(Guide, Library), color=z, size=pmin(5, -log10(padj)))) + 
  theme_bw(12) +
  geom_point() +
  scale_color_gradient2() +
  facet_grid(Gene ~ ., scales = "free", space = "free")
ggsave(out("Check_inconsistencies.pdf"), w=10,h=29)



# Method validation with replicates ---------------------------------------
cols.merge <- c("Comparison", "Guide", "Gene")
cols.vals <- c("Ratio.log2", "p")
repDT <- copy(REPLICATES.SCORES)
repDT[,Ratio.log2 := logFC]
repDT[,p := P.Value]
repDT[,Comparison := toupper(Comparison)]
repDT[, Guide := rn]
repDT[, Gene := gsub("_.+$", "", Guide)]
wtDT <- copy(RESULTS.wt)
wtDT[,Comparison := toupper(Comparison)]
wtDT[,Ratio.log2 := Score]
pDT <- merge.data.table(repDT[,c(cols.merge, cols.vals),with=F], wtDT[Genotype == "Cas9"][,c(cols.merge, cols.vals, "Library"),with=F], by=cols.merge, suffixes = c("_rep", "_wt"))
pDT[,cor(Ratio.log2_rep, Ratio.log2_wt), by=c("Comparison")]
pDT[,corS(Ratio.log2_rep, Ratio.log2_wt), by=c("Comparison")]
pDT[,cor(-log10(p_rep), -log10(p_wt)), by=c("Comparison")]
pDT[,corS(-log10(p_rep), -log10(p_wt)), by=c("Comparison")]

# Ratio log2
ggplot(pDT, aes(x=Ratio.log2_wt, y=Ratio.log2_rep)) + 
  geom_vline(xintercept = 0, color="blue") +
  geom_hline(yintercept = 0, color="blue") +
  theme_bw(12) +
  geom_point(shape=1) + 
  facet_wrap(~ Comparison + Library, scales = "free")
ggsave(out("Methods_comparison_WTvsREP_log2Ratio.pdf"), w=10,h=10)

# log10 p-value
ggplot(pDT, aes(x=sign(Ratio.log2_wt) * -log10(p_wt), y=sign(Ratio.log2_rep) * -log10(p_rep))) + 
  geom_vline(xintercept = 0, color="blue") +
  geom_hline(yintercept = 0, color="blue") +
  theme_bw(12) +
  geom_point(shape=1) + 
  facet_wrap(~ Comparison + Library, scales = "free")
ggsave(out("Methods_comparison_WTvsREP_log10p.pdf"), w=10,h=10)



# Number of genes ---------------------------------------------------------
ggplot(RESULTS.wt[Genotype != "WT" & !grepl("NonT", Gene),length(unique(Guide)), by=c("Gene", "Comparison")], aes(x=Comparison, fill=factor(V1))) + 
  theme_bw(12) +
  geom_bar() +
  xRot()
ggsave(out("Numbers.pdf"), w=5,h=5)


# Example calculation and scores -----------------------------------------------------
ggplot(SCORES.EXAMPLE, aes(x=Score)) + 
  geom_density(data=data.table(Score=rnorm(10000, mean=SCORES.EXAMPLE.bg$mean, sd=SCORES.EXAMPLE.bg$sd)), fill="#b2df8a", color=NA) +
  scale_color_manual(values=COLOR.Genotypes) +
  geom_density(aes(color=Genotype)) +
  theme_bw(12)
ggsave(out("Scores_Example_Calculation.pdf"), w=5,h=3)

# Scores
str(gg <- unique(RESULTS.wt[Comparison == SCORES.EXAMPLE.comp][Library == SCORES.EXAMPLE.lib][hit == TRUE]$Gene))
pDT <- SCORES.EXAMPLE
pDT[,Gene := gsub("_.+$", "", Guide)]
pDT[, label := ifelse(Genotype == "WT", "WT", Gene)]
pDT$label <- factor(pDT$label, levels = unique(c("WT", pDT[,median(Score), by=c("label")][order(V1)]$label)))
ggplot(pDT[Gene %in% gg | grepl("^NonT", label) | Genotype == "WT"], aes(x=label, y=Score, color=grepl("^NonT", label))) + 
  geom_boxplot(color="lightgrey") + 
  scale_color_manual(values=c("black", "red")) + 
  theme_bw(12) +
  geom_hline(color="lightblue", yintercept = 0) +
  geom_jitter(height=0, width=0.2, shape=1) +
  xRot()
ggsave(out("Scores_Example_Scores.pdf"), w=10,h=4)


# Distribution
compx <- "MYE.UND"
pDT <- RESULTS.wt[Comparison == compx][!grepl("NonT", Gene)]
pDT.median <- pDT[,.(z=mean(z), sd=sd(z), n=.N), by=c("Gene")]
pDT.median[, se := sd/sqrt(n)]
pDT.median[Gene %in% RESULTS.wt.agg.gene[Comparison == compx][hit == TRUE]$Gene, hit := TRUE]
pDT.median[is.na(hit), hit := FALSE]
pDT$Gene <- factor(pDT$Gene, levels=c(pDT.median[order(z)]$Gene))
ggplot(pDT, aes(x=Gene, y=z)) + 
  theme_bw(12) +
  geom_point(shape=1) +
  geom_point(data=pDT.median, aes(color=hit), shape=16) +
  xRot()
ggsave(out("Scores_Example_Distribution.pdf"), w=70,h=5, limitsize = FALSE)

# Mean and SE
pDT.median$Gene <- factor(pDT.median$Gene, levels=levels(pDT$Gene))
ggplot(pDT.median, aes(x=Gene, y=z)) + 
  theme_bw(12) +
  geom_point(aes(color=hit)) +
  geom_errorbar(aes(ymax=z-sd, ymin=z+sd)) +
  xRot()
ggsave(out("Scores_Example_Distribution_SD.pdf"), w=70,h=5, limitsize = FALSE)

ggplot(pDT.median, aes(x=Gene, y=z)) + 
  theme_bw(12) +
  geom_point() +
  geom_point(data=pDT.median[hit == TRUE], color="red") +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  xlab("Target genes")
ggsave(out("Scores_Example_Distribution_Simple.pdf"), w=5,h=4)


# Validation with replicates ----------------------------------------------
ggplot(REPLICATES.VALUES, aes(x=paste(Genotype, Population), shape=factor(sample_i), y=value, color=factor(guide_i))) + 
  theme_bw(12) +
  geom_point() +
  facet_grid(~ Gene) +
  xRot()
ggsave(out("Valiation_Replicates.pdf"), w=12,h=4)


# Networks in one dimension ---------------------------------------------------------------
compDT2 <- unique(melt(FIGS.COMPARISONS[,-"Comparison",with=F], id.vars = "Comparison.Group")[,-"variable"])

# Aggregated data annotated with comparisons
agDT <- copy(RESULTS.wt.agg)
mainBranchOrdering <- c("Und", "MEP", "LSKd7", "GMP", "Mye")
agDT$Population <- factor(agDT$variable, levels=c(mainBranchOrdering, setdiff(agDT$variable, mainBranchOrdering)))
agDT <- merge(agDT, compDT2, by.y="value", by.x="variable", allow.cartesian=TRUE)
agDT[,id := paste(rn, Population,Comparison.Group)]
agDT <- hierarch.ordering(agDT, toOrder = "rn", orderBy = "id", value.var = "log2FC")
agDT <- merge(agDT, ANN.genes, by.x="rn", by.y="GENE", all.x=TRUE)
#agDT$rn <- factor(agDT$rn, levels=RESULTS.wt.mds$rn[hclust(dist(as.matrix(data.frame(RESULTS.wt.mds[,c("V1", "V2"), with=F]))))$order])

# Only show genes where 50 % of guides are significant
pDT.stats <- RESULTS.wt.agg.gene[hit == TRUE][complete.screen == TRUE]

# Annotate clusters
pDT.stats <- merge(pDT.stats, ANN.genes[,c("GENE", "Complex_simple"),with=F], by.x="Gene", by.y="GENE", all.x=TRUE)

# Plot
ggplot(agDT[rn %in% pDT.stats$Gene], aes(x=Population, y=rn)) +
  theme_bw(12) +
  geom_point(aes(fill=log2FC), shape=21, color="white", size=5) +
  facet_grid(Complex_simple ~ cleanComparisons2(Comparison.Group), scales = "free", space = "free") +
  geom_segment(data=pDT.stats, aes(xend=Population1, x=Population2, y=Gene, yend=Gene, color=z), arrow=arrow(type="closed", length = unit(0.3, "cm"))) +
  scale_fill_gradient2(name=TeX(r'($\\overset{\Delta_{Cas9-WT}}{(dots)}$)')) +
  #geom_point(aes(fill=log2FC), shape=21, color="white", size=2) +
  scale_color_gradient2(name=TeX(r'($\\overset{\Delta_{Populations}}{(arrows)}$)')) +
  theme(strip.text.y = element_text(angle=0)) + 
  xRot()
ggsave(out("Aggregated_Edges.pdf"), w=8,h=18)


# Selected comparisons (David) --------------------------------------------
pDT.stats <- copy(RESULTS.wt.agg.gene)
pDT.stats <- unique(pDT.stats[,-"Comparison.Group"])
pDT.stats <- pDT.stats[Gene %in% pDT.stats[hit == TRUE][complete.screen == TRUE]$Gene]
pDT.stats <- pDT.stats[Comparison %in% COMPARISONS.healthy]
pDT.stats <- merge(pDT.stats, ANN.genes, by.x="Gene", by.y="GENE", all.x=TRUE)
# mx <- toMT(pDT.stats, row = "Gene", col = "Comparison", val = "z")
# mx[is.na(mx)] <- 0
# pDT.stats$Cluster <- cutree(hclust(dist(mx)), k = 7)[pDT.stats$Gene]
pDT.stats <- hierarch.ordering(pDT.stats, toOrder = "Gene", orderBy = "Comparison", value.var="z")
write.tsv(pDT.stats, out("SimpleHM.tsv"))

# Plot
ggplot(pDT.stats, aes(x=cleanComparisons(Comparison, ggtext = TRUE), y=Gene)) +
  theme_bw(12) + 
  geom_point(aes(color=z, size=percSig)) +
  scale_color_gradient2(name=TeX(r'($\\overset{\Delta_{Cas9-WT}}$)')) + #, low="#e31a1c", high="#1f78b4") +
  facet_grid(Complex_simple ~ ., scales = "free", space = "free") +
  theme(strip.text.y = element_text(angle=0)) +
  xRot() + 
  ylab("") +
  theme(axis.text.x = element_markdown())
ggsave(out("SimpleHM.pdf"), w=5,h=16)


# Vulcano plots -----------------------------------------------------------
ggplot(RESULTS.wt[Comparison %in% COMPARISONS.healthy], aes(x=z, y=-log10(p))) + 
  geom_hex() +
  facet_grid(. ~ cleanComparisons(Comparison)) +
  theme_bw(12)
ggsave(out("Vulcano.pdf"), w=15,h=2)


# Comparison of two main populations --------------------------------------
cx <- c("GMP.MEP", "LSK.CKIT")
pDT <- RESULTS.wt.agg.gene[Comparison %in% cx]
pDT.sig <- pDT[hit == TRUE][,paste(sort(unique(cleanComparisons(Comparison))), collapse = ","), by="Gene"]
pDT.sig[grepl(",", V1),V1 := "Both"]
pDT <- dcast.data.table(pDT, Gene ~ Comparison, value.var = "z")
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
  geom_text_repel(data=pDT[sig != "None"], aes(label=Gene)) +
  geom_point(alpha=0.5) +
  scale_color_manual(name="Significant in", values=c("#e31a1c", "#1f78b4", "#6a3d9a", "grey")) +
  #facet_grid(. ~ Genotype) +
  theme_bw() +
  xlab(TeX(formatArrows(cx[1]))) +
  ylab(TeX(formatArrows(cx[2]))) +
  ylim(-dlim, dlim) + xlim(-dlim,dlim)
ggsave(out("Comparison_Scatter.pdf"), w=5.5,h=4.5)


# Plot networks ------------------------------------------------------------
outGraphs <- dirout(paste0(base.dir, "/Graphs/"))
genex <- "Kmt2d"
for(genex in unique(RESULTS.wt.agg.gene[hit == TRUE][complete.screen == TRUE]$Gene)){
  COLORS.graph <- c("#fb9a99", "lightgrey", "#a6cee3")
  
  agx <- RESULTS.wt.agg[rn == genex]
  agx <- unique(agx[,c("variable", "log2FC"), with=F])
  
  statx <- RESULTS.wt.agg.gene[Gene == genex]
  
  el <- data.table(do.call(rbind, COMPARISONS), keep.rownames = TRUE)
  #statx <- statx[match(row.names(el), Comparison)][!is.na(Gene)]
  statx <- merge(statx, el, by.x="Comparison", by.y="rn")
  statx[Comparison == "LSK.CKIT", V1 := "LSKd7"]
  statx[Comparison == "LSK.CKIT", V2 := "LSKd7"]
  # statx[Comparison == "GMPcd11.DN", V2 := "MEP"]
  # statx[Comparison == "GMPcd11.DN", V1 := "Und"]
  # statx  <- statx[Comparison != "UND.MEP"]
  
  g <- graph.edgelist(as.matrix(statx[,c("V2", "V1")]))
  V(g)$log2FC <- agx[match(V(g)$name, variable),]$log2FC
  E(g)$z <- statx$z
  E(g)$up <- statx$sig.up
  E(g)$dn <- statx$sig.down
  E(g)$n <- statx$n
  E(g)$sig.label <- with(statx, paste(sig.up, sig.down, n, sep="/"))
  E(g)$cnt <- ifelse(E(g)$z > 0, E(g)$up, E(g)$dn)
  E(g)$perc <- round(E(g)$cnt/E(g)$n) * 100
  
  g <- delete.vertices(g, is.na(V(g)$log2FC))
  V(g)$color <- mapNumericToColors(V(g)$log2FC, cols = COLORS.graph)
  #V(g)$color <- "lightgrey"
  V(g)$frame.color <- NA
  V(g)$label.color <- "black"
  E(g)$color <- mapNumericToColors(E(g)$z, cols = COLORS.graph)
  E(g)$width <- 2+E(g)$perc/100*3
  E(g)$arrow.width <- 2
  E(g)$arrow.size <- 1
  E(g)$label <- paste0("", round(E(g)$z, 1),"\n", "(",E(g)$sig.label, ")")
  E(g)$label.color <- "black"
  
  
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



# Heatmaps for networks -----------------------------------------------------
res.stats <- RESULTS.wt.agg.gene[hit == TRUE][complete.screen == TRUE]
res.stats[abs(z) >= 5, z := 5 * sign(z)]

# Generate background
bg <- data.table()
for(gx in unique(res.stats$Gene)){
  for(cx in unique(res.stats$Comparison)){
    bg <- rbind(bg, data.table(Gene = gx, Comparison = cx, z=-5:5))
  }
}

# Plotting
ggplot(res.stats, aes(x=Gene)) + 
  facet_grid(cleanComparisons(Comparison) ~ .) +
  geom_tile(data=bg, aes(fill=z, y=factor(z))) + 
  scale_fill_gradient2(high="#a6cee3", low="#fb9a99") + 
  geom_point(aes(y=factor(round(z))), size=3, color="black", shape=18) +
  theme_bw(12) + 
  xRot()
ggsave(out("GraphHM.pdf"),w=17, h=10)



# TODO --> example of how p-values were calculated