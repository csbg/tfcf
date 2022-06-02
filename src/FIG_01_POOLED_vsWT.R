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
RESULTS.wt.agg <- fread(inDir("ComparisonToWT.tsv"))


# Gene annotations --------------------------------------------------------
ANN.genes <- fread("metadata/TFCF_Annotations_v2.tsv", check.names = TRUE)
ANN.genes[,Complex_simple := BROAD.COMPLEX]
ANN.genes <- unique(ANN.genes[!is.na(Complex_simple) & Complex_simple != ""][,c("GENE", "Complex_simple"), with=F])
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
pDT.stats <- merge(RESULTS.wt, unique(FIGS.COMPARISONS[,c("Comparison", "Comparison.Group")]), by="Comparison", all.x=TRUE)
pDT.stats <- pDT.stats[, .(
  z=mean(z), 
  sig.up=length(unique(Guide[padj < 0.05 & z > 0])),
  sig.down=length(unique(Guide[padj < 0.05 & z < 0])),
  n=length(unique(Guide))
  ), by=c("Gene", "Comparison",  "Comparison.Group", "Population1", "Population2", "Genotype")]
pDT.stats[, sig := ifelse(z > 0, sig.up, sig.down)]
pDT.stats[,percSig := sig/n*100]
stopifnot(all(pDT.stats$percSig >= 0))
stopifnot(all(pDT.stats$percSig <= 100))
pDT.stats[, hit := percSig >= 50 & n > 1]
pDT.stats[, complete.screen := all(COMPARISONS.healthy[1:4] %in% Comparison), by=c("Gene", "Genotype")]
RESULTS.wt.agg.gene <- pDT.stats[Genotype == "Cas9"]
RESULTS.wt.agg.gene.wt <- pDT.stats[Genotype == "WT"]

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
REPLICATES.SCORES.EXAMPLE <- REPLICATES.SCORES[coef == gsub("\\.", "vs", SCORES.EXAMPLE.comp)]#[grepl(SCORES.EXAMPLE.lib, analysis)]

# Values
ff <- lapply(unique(REPLICATES.SCORES.EXAMPLE$analysis), function(fx) fread(replicationDir("analysis_", fx, "/", "Data_DateRemoved.csv")))
ff <- lapply(ff, function(fx){ melt(data.table(fx), id.vars="rn")})
ff <- do.call(rbind, ff)
REPLICATES.VALUES <- toMT(ff, row = "rn", col = "variable", val = "value")
guidesx <- unique(REPLICATES.SCORES.EXAMPLE[gene %in% REPLICATES.SCORES.EXAMPLE[, sum(significant_TF) / .N, by="gene"][order(V1, decreasing=TRUE)][1:6]$gene]$rn)
#guidesx <- guidesx[!grepl("Spi1", guidesx)]
samplesx <- ann[Population %in% COMPARISONS[[SCORES.EXAMPLE.comp]]]$sample
REPLICATES.VALUES <- REPLICATES.VALUES[guidesx,intersect(samplesx, colnames(REPLICATES.VALUES))]
REPLICATES.VALUES <- melt(data.table(REPLICATES.VALUES, keep.rownames = TRUE), id.vars = "rn")
REPLICATES.VALUES <- merge(ann, REPLICATES.VALUES, by.x="sample", by.y="variable")
REPLICATES.VALUES[, Gene := gsub("_.+", "", rn)]
REPLICATES.VALUES <- REPLICATES.VALUES[!is.na(value)]
REPLICATES.VALUES[, sample_i := as.numeric(factor(sample)), by=c("Gene", "rn", "Genotype")]
REPLICATES.VALUES[, guide_i := as.numeric(factor(rn)), by=c("Gene")]

# # Values
# REPLICATES.VALUES <- fread(replicationDir("analysis_", REPLICATES.SCORES.EXAMPLE$analysis[1], "/", "Data_DateRemoved.csv"))
# rns <- REPLICATES.VALUES$rn
# REPLICATES.VALUES <- as.matrix(REPLICATES.VALUES[, -"rn"])
# row.names(REPLICATES.VALUES) <- rns
# guidesx <- unique(REPLICATES.SCORES.EXAMPLE[gene %in% REPLICATES.SCORES.EXAMPLE[significant_TF == TRUE][,.N, by="gene"][N >= 2]$gene]$rn)
# guidesx <- guidesx[!grepl("Spi1", guidesx)]
# REPLICATES.VALUES <- REPLICATES.VALUES[guidesx,ann[Library == SCORES.EXAMPLE.lib][Population %in% COMPARISONS[[SCORES.EXAMPLE.comp]]]$sample]
# REPLICATES.VALUES <- melt(data.table(REPLICATES.VALUES, keep.rownames = TRUE), id.vars = "rn")
# REPLICATES.VALUES <- merge(ann, REPLICATES.VALUES, by.x="sample", by.y="variable")
# REPLICATES.VALUES[, Gene := gsub("_.+", "", rn)]
# REPLICATES.VALUES[, sample_i := rank(sample), by=c("Gene", "rn", "Genotype")]
# REPLICATES.VALUES[, guide_i := rank(rn), by=c("Gene")]



# SETUP ENDS HERE ---------------------------------------------------------



# Check inconsistencies ---------------------------------------------------
pDT <- merge(RESULTS.wt[Genotype == "Cas9"], RESULTS.wt.agg.gene[sig.down != 0 & sig.up != 0][,c("Gene", "Comparison"),with=F],by=c("Gene", "Comparison"))
ggplot(pDT, aes(x="x", y=paste(Guide, Library), color=z, size=pmin(5, -log10(padj)))) + 
  theme_bw(12) +
  geom_point() +
  scale_color_gradient2() +
  facet_wrap(~ Gene + Comparison, scales = "free") +
  #facet_grid(Gene ~ ., scales = "free", space = "free")
  xRot()
ggsave(out("Check_inconsistencies.pdf"), w=20,h=29)

inc <- RESULTS.wt[, var(z, na.rm = TRUE), by=c("Gene", "Comparison")][order(V1, decreasing = TRUE)][1:20]
pDT <- merge(RESULTS.wt[Genotype == "Cas9"], inc, by=c("Gene", "Comparison"))
p <- ggplot(pDT, aes(x=Library, y=Guide, color=z, size=pmin(5, -log10(padj)))) + 
  theme_bw(12) +
  geom_point() +
  scale_color_gradient2() +
  facet_grid(Gene ~ Comparison, scales = "free", space = "free") +
  xRot()
ggsave(out("Check_inconsistencies2.pdf"), w=10,h=25, plot=p)


# Method validation comparing analysis vs wildtype to analysis with replicates ---------------------------------------
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
pDT[, Comparison := cleanComparisons(Comparison)]

# Ratio log2
ggplot(pDT[!Library %in% c("P1", "R1")], aes(x=Ratio.log2_wt, y=Ratio.log2_rep)) + 
  geom_vline(xintercept = 0, color="blue") +
  geom_hline(yintercept = 0, color="blue") +
  themeNF() +
  xlab(TeX("Comparison to wildtype (log_{2}FC)")) + 
  ylab(TeX("Replicate analysis (log_{2}FC)")) + 
  geom_point(shape=1, alpha=0.3) + 
  facet_grid(Comparison ~ Library, scales = "free")
ggsaveNF(out("Methods_comparison_WTvsREP_log2Ratio.pdf"), w=2,h=2, guides = TRUE)

# log10 p-value
# ggplot(pDT, aes(x=sign(Ratio.log2_wt) * -log10(p_wt), y=sign(Ratio.log2_rep) * -log10(p_rep))) + 
#   geom_vline(xintercept = 0, color="blue") +
#   geom_hline(yintercept = 0, color="blue") +
#   theme_bw(12) +
#   geom_point(shape=1) + 
#   facet_wrap(~ Comparison + Library, scales = "free")
# ggsave(out("Methods_comparison_WTvsREP_log10p.pdf"), w=10,h=10)


# Number of genes ---------------------------------------------------------
pDT <- RESULTS.wt[Genotype != "WT" & !grepl("NTC", Gene),length(unique(Guide)), by=c("Gene", "Comparison")]
pDT[,Comparison := cleanComparisons(Comparison)]
ggplot(pDT, aes(x=Comparison, fill=factor(V1))) + 
  themeNF() +
  geom_bar() +
  xRot() +
  xlab("") + ylab("Number of genes")
ggsaveNF(out("Numbers.pdf"), w=0.8, h=1)


# Example calculation and scores -----------------------------------------------------
ggplot(SCORES.EXAMPLE[!(Genotype == "Cas9" & grepl("^NTC", Guide))], aes(x=Score)) + 
  geom_density(data=data.table(Score=rnorm(10000, mean=SCORES.EXAMPLE.bg$mean, sd=SCORES.EXAMPLE.bg$sd)), fill="#b2df8a", color=NA) +
  scale_color_manual(values=COLOR.Genotypes) +
  geom_density(aes(color=Genotype)) +
  themeNF() +
  xlab(TeX(paste(cleanComparisons(SCORES.EXAMPLE[1]$Comparison), "log_{2}FC"))) +
  ylab("Density")
ggsaveNF(out("Scores_Example_Calculation.pdf"), w=1.2)

# Scores
str(gg <- unique(RESULTS.wt[Comparison == SCORES.EXAMPLE.comp][Library == SCORES.EXAMPLE.lib][hit == TRUE]$Gene))
pDT <- SCORES.EXAMPLE
pDT[,Gene := gsub("_.+$", "", Guide)]
pDT[, label := ifelse(Genotype == "WT", "WT", Gene)]
pDT[grepl("^NTC", Gene), label := "NTC"]
pDT$label <- factor(pDT$label, levels = unique(c("WT", pDT[,median(Score), by=c("label")][order(V1)]$label)))
pDT <- pDT[Gene %in% gg | grepl("^NTC", label) | Genotype == "WT"]
ggplot(pDT, aes(x=label, y=Score, color=grepl("^NTC", label))) + 
  geom_boxplot(color="lightgrey") + 
  scale_color_manual(values=c("black", "red")) + 
  themeNF() +
  geom_hline(color="lightblue", yintercept = 0) +
  geom_jitter(height=0, width=0.2, shape=1) +
  xRot() +
  ylab(TeX(paste(cleanComparisons(pDT[1]$Comparison), "log_{2}FC"))) +
  xlab("")
ggsaveNF(out("Scores_Example_Scores.pdf"), w=2, h=1)


# Distribution
compx <- "MYE.UND"
pDT <- RESULTS.wt[Genotype == "Cas9"][Comparison == compx][!grepl("NTC", Gene)]
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
ggplot(pDT.median, aes(y=Gene, x=z)) + 
  themeNF() +
  geom_errorbarh(aes(xmax=z-sd, xmin=z+sd), size=0.5) +
  geom_point(aes(color=hit), size=0.5) +
  scale_color_manual(values=c("grey", "#e31a1c")) +
  theme(axis.text.y = element_text(size=2)) +
  theme(panel.grid = element_blank()) +
  xlab("z-Scores")
ggsaveNF(out("Scores_Example_Distribution_SD.pdf"), w=1,h=5, limitsize = FALSE)

ggplot(pDT.median, aes(x=Gene, y=z)) + 
  themeNF() +
  geom_hline(yintercept = 0) +
  geom_point(color="grey") +
  geom_point(data=pDT.median[hit == TRUE], color="#e31a1c", shape=1) +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  xlab("Target genes") +
  ylab("z-Scores") +
  expand_limits(x= c(-10, length(unique(pDT.median$Gene)) + 10))
ggsaveNF(out("Scores_Example_Distribution_Simple.pdf"), w=1,h=1)


# Result validation with replicates ----------------------------------------------
pDT <- dcast.data.table(REPLICATES.VALUES, guide_i + sample_i + Population + Gene ~ Genotype, value.var = "value")
pDT[, diff := Cas9 - WT]
ggplot(pDT, aes(x=factor(guide_i), y=diff, fill=paste(Population, sample_i))) + 
  #geom_point() +
  geom_bar(stat="identity", position="dodge") +
  themeNF() +
  facet_grid(.~Gene, scales = "free", space = "free") +
  #geom_vline(xintercept = 2.5) +
  geom_hline(yintercept = 0) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  #coord_flip() +
  #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("Guides") + ylab(TeX("Change over WT (log_{2}FC)")) +
  xRot()
ggsaveNF(out("Validation_Replicates.pdf"), w=2,h=1)


# Networks in one dimension ---------------------------------------------------------------
compDT2 <- unique(melt(FIGS.COMPARISONS[,-"Comparison",with=F], id.vars = "Comparison.Group")[,-"variable"])

# Aggregated data annotated with comparisons
agDT <- copy(RESULTS.wt.agg)
mainBranchOrdering <- c("Und", "MEP", "LSKd7", "GMP", "Mye")
agDT$Population <- factor(agDT$Population, levels=c(mainBranchOrdering, setdiff(agDT$Population, mainBranchOrdering)))
agDT <- merge(agDT, compDT2, by.x="Population", by.y="value", allow.cartesian=TRUE)
agDT[,id := paste(Gene, Population,Comparison.Group)]
agDT <- hierarch.ordering(agDT, toOrder = "Gene", orderBy = "id", value.var = "log2FC")
agDT <- merge(agDT, ANN.genes, by.x="Gene", by.y="GENE", all.x=TRUE)
#agDT$rn <- factor(agDT$rn, levels=RESULTS.wt.mds$rn[hclust(dist(as.matrix(data.frame(RESULTS.wt.mds[,c("V1", "V2"), with=F]))))$order])

# Comparisons between popuulations
pDT.stats <- RESULTS.wt.agg.gene[hit == TRUE][complete.screen == TRUE][!grepl("^NTC", Gene)]
pDT.stats <- merge(pDT.stats, ANN.genes[,c("GENE", "Complex_simple"),with=F], by.x="Gene", by.y="GENE", all.x=TRUE)

# Plot
ggplot(agDT[Gene %in% pDT.stats$Gene], aes(x=Population, y=Gene)) +
  themeNF() +
  geom_point(aes(fill=pmin(2, abs(log2FC)) * sign(log2FC)), shape=21, color="white", size=5) +
  facet_grid(Complex_simple ~ cleanComparisons(Comparison.Group), scales = "free", space = "free") +
  geom_segment(data=pDT.stats, aes(xend=Population1, x=Population2, y=Gene, yend=Gene, color=pmin(5, abs(z)) * sign(z)), arrow=arrow(type="closed", length = unit(0.3, "cm"))) +
  scale_fill_gradient2(name=TeX(r'($\\overset{\Delta_{Cas9-WT}}{(dots)}$)')) +
  #geom_point(aes(fill=log2FC), shape=21, color="white", size=2) +
  scale_color_gradient2(name=TeX(r'($\\overset{\Delta_{Populations}}{(arrows)}$)')) +
  theme(strip.text.y = element_text(angle=0)) +
  xRot()
ggsaveNF(out("Aggregated_Edges.pdf"), w=3,h=6, guide=TRUE)



# Selected comparisons (David) --------------------------------------------
dla.list <- list(
  main = fread("metadata/FIGS_Order_Fig1E.tsv"),
  supp = fread("metadata/FIGS_Order_Fig1E_supp.tsv"),
  all = data.table(Factor=sort(unique(RESULTS.wt.agg.gene[hit == TRUE]$Gene)), Complex="NA")
)
(dla.nam <- names(dla.list)[1])
for(dla.nam in names(dla.list)){
  dla <- dla.list[[dla.nam]]
  if(is.null(dla$Heatmap)) dla$Heatmap <- 1
  pDT.stats <- copy(RESULTS.wt.agg.gene)
  pDT.stats <- unique(pDT.stats[,-"Comparison.Group"])
  pDT.stats <- pDT.stats[Gene %in% dla$Factor]
  pDT.stats$Gene <- factor(pDT.stats$Gene, levels=dla$Factor)
  pDT.stats$ComplexDLA <- factor(dla[match(pDT.stats$Gene, Factor)]$Complex, levels=unique(dla$Complex))
  # pDT.stats <- pDT.stats[Gene %in% c(pDT.stats[hit == TRUE][complete.screen == TRUE]$Gene, "Smarcd1", "Ezh2", "Rcor2")]
  # pDT.stats <- pDT.stats[!Gene %in% "Ctcf"]
  pDT.stats <- pDT.stats[Comparison %in% COMPARISONS.healthy]
  pDT.stats <- merge(pDT.stats, ANN.genes, by.x="Gene", by.y="GENE", all.x=TRUE)
  #pDT.stats <- hierarch.ordering(pDT.stats, toOrder = "Gene", orderBy = "Comparison", value.var="z")
  write.tsv(pDT.stats, out("SimpleHM_",dla.nam,".tsv"))
  
  # Plot
  pDT.stats[, z.cap := pmin(5, abs(z)) * sign(z)]
  p <- ggplot(pDT.stats[Genotype == "Cas9"], aes(
    y=cleanComparisons(Comparison, ggtext = TRUE, reverse = TRUE, colors=c("0000B3", "A60000")), 
    x=Gene)
  ) +
    themeNF() + 
    geom_point(aes(fill=z.cap, size=percSig), shape=21, color="lightgrey") +
    facet_grid(. ~ ComplexDLA, scales = "free", space = "free") +
    theme(strip.text.x = element_text(angle=90)) +
    scale_size_continuous(range = c(2,4)) +
    xRot() + 
    ylab("") +
    theme(axis.text.y = element_markdown()) +
    xlab("") +
    theme(panel.spacing = unit(0.01, "cm"))
  w=length(unique(pDT.stats$Gene)) * 0.06 + 0.6
  ggsaveNF(out("SimpleHM_",dla.nam,"_RedBlue.pdf"), w=w,h=1, limitsize = FALSE,
           plot = p + scale_fill_gradient2(
             name=TeX(r'($\\overset{\Delta_{Cas9-WT}}$)'), 
             low="#0000B3", high="#A60000"))
  ggsaveNF(out("SimpleHM_",dla.nam,"_GreenPurple.pdf"), w=w,h=1, limitsize = FALSE,
           plot = p + scale_fill_gradient2(
             name=TeX(r'($\\overset{\Delta_{Cas9-WT}}$)'), 
             high="#0D5E04", low="#5700C2"))
  
  # cleanDev(); pdf(out("SimpleHM_Dendrogram.pdf"), w=15,h=5)
  # plot(hclust(dist(toMT(pDT.stats[Genotype == "Cas9"], row = "Gene", col = "Comparison", val = "z"))))
  # dev.off()
  # 
  # require(umap)
  # umObj <- toMT(pDT.stats[Genotype == "Cas9"], row = "Gene", col = "Comparison", val = "z")
  # umObj[is.na(umObj)] <- 0
  # umObj <- umap(umObj)
  # umap <- data.table(umObj$layout, keep.rownames = TRUE)
  # umap <- setNames(umap, c("Gene", "UMAP1", "UMAP2"))
  # ggplot(umap, aes(x=UMAP1, y=UMAP2)) + 
  #   themeNF() +
  #   geom_point(color="#1f78b4") +
  #   geom_text_repel(aes(label = Gene))
  # ggsaveNF(out("SimpleHM_UMAP.pdf"), w=2,h=2, guides = TRUE)
}


# Vulcano plots -----------------------------------------------------------
# ggplot(RESULTS.wt[Comparison %in% COMPARISONS.healthy], aes(x=z, y=-log10(p))) + 
#   geom_hex() +
#   facet_grid(. ~ cleanComparisons(Comparison)) +
#   theme_bw(12)
# ggsave(out("Vulcano.pdf"), w=15,h=2)


# Comparison of two main populations (2D - Scatterplot) --------------------------------------
cap <- 5

xxx <- list(
  c("GMP.MEP", "CKIT.LSK"),
  c("MYE.UND", "GMPcd11.DN")
)

for(cx in xxx){
  # Prepare data
  pDT <- rbind(RESULTS.wt.agg.gene.wt, RESULTS.wt.agg.gene)
  pDT <- pDT[Comparison %in% cx]
  pDT.sig <- pDT[Genotype == "Cas9"][hit == TRUE][,paste(sort(unique(cleanComparisons(Comparison))), collapse = ","), by="Gene"]
  pDT.sig[grepl(",", V1),V1 := "Both"]
  pDT <- dcast.data.table(pDT, Gene + Genotype ~ Comparison, value.var = "z")
  pDT <- merge(pDT, pDT.sig, by="Gene", all.x=TRUE)
  pDT[, sig := V1]
  pDT[is.na(sig), sig := "None"]
  pDT$sig <- factor(pDT$sig, levels=c("Both", cleanComparisons(cx, order = FALSE), "None"))
  pDT[, colorCode := abs(get(cx[1])) + abs(get(cx[2])) > 5]
  
  # Define dimensions
  pDT$dim1 <- pDT[[cx[1]]]
  pDT$dim2 <- pDT[[cx[2]]]
  pDT[,dim1 := pmin(cap, abs(dim1)) * sign(dim1)]
  pDT[,dim2 := pmin(cap, abs(dim2)) * sign(dim2)]
  
  # Axis labels with axes
  formatArrows <- function(x){
    paste0(gsub("^(.+?)\\.(.+)$", "\\1", x),"  ", r'($\leftarrow$)',"  .  ",r'($\rightarrow$)',"  ", gsub("^(.+?)\\.(.+)$", "\\2", x))
  }
  
  # Define complexes
  pDT$Complex <- ANN.genes[match(pDT$Gene, GENE)]$Complex_simple
  pDT[is.na(Complex), Complex := "None"]
  
  write.tsv(pDT, out("Comparison_2D_Scatter_", paste(cx, collapse = "vs"),".tsv"))
  
  dlim <- cap
  ggplot(pDT, aes(x=dim1, y=dim2)) + 
    themeNF() +
    geom_hline(yintercept = 0, color="lightgrey", alpha=0.5) +
    geom_vline(xintercept = 0, color="lightgrey", alpha=0.5) +
    stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = NA, bins = 10) +
    #scale_fill_distiller(palette = "Blues", direction = 1) +
    scale_fill_gradientn(colours=c("white", "#a6cee3", "#fdbf6f")) +
    geom_point(data=pDT[Genotype == "Cas9" & grepl("NTC", Gene)], size=1, shape=3, color="#1f78b4") +
    geom_point(data=pDT[Genotype == "Cas9" & !grepl("NTC", Gene) & sig != "None"], alpha=1, color="#e31a1c", size=2, shape=4) +
    geom_point(data=pDT[Genotype == "Cas9" & !grepl("NTC", Gene) & sig == "None"], alpha=0.5, color="#000000", size=1, shape=16) +
    geom_text_repel(data=pDT[Genotype == "Cas9" & !grepl("NTC", Gene) & colorCode == TRUE], aes(label=Gene)) +
    xlab(TeX(formatArrows(cx[1]))) +
    ylab(TeX(formatArrows(cx[2]))) +
    ylim(-dlim, dlim) + xlim(-dlim,dlim)
  ggsaveNF(out("Comparison_2D_Scatter_", paste(cx, collapse = "vs"),".pdf"), w=2,h=2)
}


# Plot networks ------------------------------------------------------------
outGraphs <- dirout(paste0(base.dir, "/Graphs/"))
genex <- "Kmt2d"
for(genex in unique(RESULTS.wt.agg.gene[hit == TRUE][complete.screen == TRUE]$Gene)){
  COLORS.graph <- c("#fb9a99", "lightgrey", "#a6cee3")
  
  agx <- RESULTS.wt.agg[Gene == genex]
  agx <- unique(agx[,c("Population", "log2FC"), with=F])
  
  statx <- RESULTS.wt.agg.gene[Gene == genex]
  
  el <- data.table(do.call(rbind, COMPARISONS), keep.rownames = TRUE)
  #statx <- statx[match(row.names(el), Comparison)][!is.na(Gene)]
  statx <- merge(statx, el, by.x="Comparison", by.y="rn")
  statx[Comparison == "CKIT.LSK", V1 := "LSKd7"]
  statx[Comparison == "CKIT.LSK", V2 := "LSKd7"]
  # statx[Comparison == "GMPcd11.DN", V2 := "MEP"]
  # statx[Comparison == "GMPcd11.DN", V1 := "Und"]
  # statx  <- statx[Comparison != "UND.MEP"]
  
  g <- graph.edgelist(as.matrix(statx[,c("V2", "V1")]))
  V(g)$log2FC <- agx[match(V(g)$name, Population),]$log2FC
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

