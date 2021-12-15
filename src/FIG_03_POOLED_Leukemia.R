source("src/00_init.R")
out <- dirout("FIG_01_POOLED_vsControls/")

require(latex2exp)
require(ggrepel)
require(igraph)
require(ggtext)




# Load data ---------------------------------------------------------------

# Individual sampels and guides
m <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix))
ann <- fread(PATHS$POOLED$DATA$annotation)

# Samples of libraries aggregated 
mA <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix.aggregated))
annA <- fread(PATHS$POOLED$DATA$annotation.aggregated)
# Aggregate by guide
mA <- t(sapply(split(row.names(mA) , gsub("_.+$", "", row.names(mA))), function(gx){
  colSums(mA[gx,], na.rm = TRUE)
}))
# cpm transformed
mA.cpm <- t(t(mA) / colSums(mA, na.rm = TRUE)) * 1e6

# Values across libraries aggregated
grps <- unique(gsub("_[A-Za-z0-9]+$", "", colnames(mA)))
mAA <-sapply(grps, function(grpx){
  message(grpx)
  colx <- grep(paste0(grpx, "_"), colnames(mA.cpm), value=TRUE)
  print(paste(colx, collapse = " - "))
  rowMeans(mA.cpm[,colx], na.rm = TRUE)
})

# Load analysis results ---------------------------------------------------
RESULTS.ctr.pval <- fread(dirout_load("POOLED_12_01_Controls")("Results_Pvalues.tsv"))
RESULTS.ctr.pval[, Library := gsub("^.+_", "", Group)]
RESULTS.ctr.pval[, Genotype := gsub("_.+$", "", Group)]

RESULTS.wt.pval <- fread(dirout_load("POOLED_10_03_IndividualAnalysis_NormFactors_Controls")("Results_Pvalues.tsv"))
RESULTS.wt.pval[, Ratio.log2 := Score]


# Aggregate results across guides
pDT.stats <- RESULTS.ctr.pval[, .(
  Ratio.log2=mean(Ratio.log2), 
  sig.up=length(unique(Guide[padj < 0.05 & Ratio.log2 > 0])),
  sig.down=length(unique(Guide[padj < 0.05 & Ratio.log2 < 0])),
  n=length(unique(Guide))
), by=c("Gene", "Comparison")]
pDT.stats[, sig := ifelse(Ratio.log2 > 0, sig.up, sig.down)]
pDT.stats[,percSig := sig/n*100]
stopifnot(all(pDT.stats$percSig >= 0))
stopifnot(all(pDT.stats$percSig <= 100))
pDT.stats[, hit := percSig >= 50]
RESULTS.AGG <- copy(pDT.stats)


# Define comparisons that are to be combined ------------------------------
COMBINATIONS <- list(
  healthy=list(m="GMP.LSK", e="MEP.LSK", i="GMP.MEP"),
  #h2=list(m="", e="", i="MYE.UND"),
  leukemia=list(m="DMMye.LSC", e="DMEry.LSC", i="DMMye.Ery")
)
SAMPLES <- list(
  healthy=list(m="Cas9_GMP", e="Cas9_MEP", i="Cas9_LSKd7"),
  leukemia=list(m="Cas9_DM.Mye", e="Cas9_DM.Ery", i="Cas9_DM.LSC")
)


# SETUP ENDS HERE ---------------------------------------------------------




# Compare methods for p-values --------------------------------------------
cols.merge <- c("Genotype", "Library", "Comparison", "Guide", "Gene")
cols.vals <- c("Ratio.log2", "z", "p", "padj")
pDT <- merge.data.table(RESULTS.ctr.pval[,c(cols.merge, cols.vals),with=F], RESULTS.wt.pval[,c(cols.merge, cols.vals),with=F], by=cols.merge, suffixes = c("_ctr", "_wt"))
pDT[,cor(z_wt, z_ctr), by=c("Genotype", "Library", "Comparison")]
pDT[,corS(z_wt, z_ctr), by=c("Genotype", "Library", "Comparison")]
pDT[,cor(Ratio.log2_wt, Ratio.log2_ctr), by=c("Genotype", "Library", "Comparison")]
pDT[,corS(Ratio.log2_wt, Ratio.log2_ctr), by=c("Genotype", "Library", "Comparison")]
pDT[,cor(-log10(p_wt), -log10(p_ctr)), by=c("Genotype", "Library", "Comparison")]
pDT[,corS(-log10(p_wt), -log10(p_ctr)), by=c("Genotype", "Library", "Comparison")]

# Ratio log2
ggplot(pDT[Genotype == "Cas9"], aes(x=Ratio.log2_wt, y=Ratio.log2_ctr)) + 
  theme_bw(12) +
  geom_point(shape=1) + 
  facet_wrap(~ Comparison + Library, scales = "free")
ggsave(out("Methods_comparison_WTvsCTRL_log2Ratio.pdf"), w=10,h=10)

# log10 p-value
ggplot(pDT[Genotype == "Cas9"], aes(x=sign(Ratio.log2_wt) * -log10(p_wt), y=sign(Ratio.log2_ctr) * -log10(p_ctr))) + 
  theme_bw(12) +
  geom_point(shape=1) + 
  facet_wrap(~ Comparison + Library, scales = "free")
ggsave(out("Methods_comparison_WTvsCTRL_log10p.pdf"), w=10,h=10)



# Scatterplot with arrows -------------------------------------------------------------

# Combine logFCs and p-values in one table
xDT_r <- dcast.data.table(
  data = RESULTS.AGG[Comparison %in% do.call(c, lapply(COMBINATIONS, unlist))], 
  formula = Gene ~ Comparison, 
  value.var = "Ratio.log2")
xDT_p <- dcast.data.table(
  data = RESULTS.AGG[Comparison %in% do.call(c, lapply(COMBINATIONS, unlist))], 
  formula = Gene ~ Comparison, 
  value.var = "percSig")
id.cols <- setdiff(colnames(xDT_r), unique(RESULTS.ctr.pval$Comparison))
xDT <- merge(xDT_r, xDT_p, by=id.cols, all=TRUE, suffixes=c("_r", "_p"))
stopifnot(all(duplicated(xDT$Gene)==FALSE))

# Triangle can be calculated from the other comparisons 1
ggplot(xDT, aes(x=DMEry.LSC_r, y=DMMye.LSC_r, color=DMMye.Ery_r)) + 
  theme_bw() +
  geom_point() + 
  scale_color_gradient2()
ggsave(out("ContraTriangle1.pdf"), w=5,h=4)

# Triangle can be calculated from the other comparisons 2 
ggplot(xDT, aes(x=DMMye.LSC_r - DMEry.LSC_r, y=DMMye.Ery_r)) + 
  ggtitle(paste("R=", round(with(xDT, cor(DMMye.LSC_r - DMEry.LSC_r, DMMye.Ery_r, use="pairwise.complete.obs")),3))) +
  theme_bw() +
  geom_point()
ggsave(out("ContraTriangle2.pdf"), w=5,h=4)

# Individual plots and prepare data for combined plot below
pDT <- data.table()
ci <- "leukemia"
for(ci in names(COMBINATIONS)){
  m <- COMBINATIONS[[ci]]$m
  e <- COMBINATIONS[[ci]]$e
  i <- COMBINATIONS[[ci]]$i
  ggplot(xDT, aes_string(x=paste0(m, "_r"), y=paste0(e, "_r"))) + 
    theme_bw() +
    geom_point(alpha=0.3) + 
    geom_abline() +
    xlab(cleanComparisons(m)) +
    ylab(cleanComparisons(e)) +
    ggtitle(ci)
  ggsave(out("Scatterplot_", ci, ".pdf"), w=5,h=4)
  
  pDT <- rbind(pDT, data.table(tissue=ci, setNames(
    xDT[,c("Gene", paste0(c(m, e, i), "_r"), paste0(c(m, e, i), "_p")), with=F],
    c("Gene", paste0(c("m", "e", "i"), "_r"), paste0(c("m", "e", "i"), "_p"))
  )))
}

# Combined plot

lDT <- xDT[
  sign(GMP.LSK_r) != sign(DMMye.LSC_r) & Gene %in% RESULTS.AGG[Comparison %in% unlist(COMBINATIONS$healthy)]$Gene |
  sign(MEP.LSK_r) != sign(DMEry.LSC_r) & Gene %in% RESULTS.AGG[Comparison %in% unlist(COMBINATIONS$leukemia)]$Gene
]
lDT <- unique(lDT)

# lDT.labels <- lDT[,.(m_r=pmean(DMMye.LSC_r, GMP.LSK_r), e_r=pmean(DMEry.LSC_r, MEP.LSK_r)), by="Gene"]
ggplot(pDT, aes(x=m_r, y=e_r, color=tissue)) + 
  ggplot2::theme_classic() +
  geom_abline(color="black", size=0.5, linetype=3) +
  geom_hline(yintercept = 0, color="black", size=0.5, linetype=3) +
  geom_vline(xintercept = 0, color="black", size=0.5, linetype=3) +
  geom_point(aes(size=pmax(m_p, e_p), alpha=pmax(m_p, e_p))) + 
  geom_segment(
    color="#666666", size=0.3, data=lDT,  alpha=1,
    arrow=arrow(length = unit(0.25, "cm")),
    aes(x=GMP.LSK_r, y=MEP.LSK_r, xend=DMMye.LSC_r, yend=DMEry.LSC_r)) +
  geom_text_repel(data=pDT[Gene %in% lDT$Gene & pmax(m_p, e_p) > 0.5], aes(label=Gene), color="black") +
  scale_color_manual(values=c(leukemia="#fb9a99", healthy="#a6cee3")) +
  scale_size_continuous(limits=c(0,100), name="% sig. gRNAs") +
  #geom_label(data=pDT[Gene %in% lDT$Gene & pmax(m_p, e_p) > 0.5], aes(label=Gene), alpha=0.3) +
  # geom_text(data=lDT.labels,  aes(label=Gene), color="black") +
  ggplot2::xlab("Average log2 ratio: GMP/LSK (across guides)") +
  ggplot2::ylab("Average log2 ratio: MEP/LSK (across guides)") +
  guides(alpha=FALSE) +
  ggtitle("Comparison of KO effects between in vitro and leukemia")
ggsave(out("Scatterplot_comb.pdf"), w=10,h=8)




# Triangle / ternary plots ------------------------------------------------
colnames(mA)
require(ggtern) # Has to stay here, otherwise breaks everything else - especially ggplot
scale.hexgradient <- scale_fill_gradientn(colours=c("#a6cee3", "#fdbf6f", "#ff7f00"))  

tx <- "leukemia"
for(tx in names(SAMPLES)){
  mx <- mAA[,unlist(SAMPLES[[tx]])]
  mx <- mx[rowSums(mx) > 0,]
  colnames(mx) <- gsub("Cas9_(DM\\.)?", "", colnames(mx))
  cx <- colnames(mx)
  mx2 <- apply(mx, 1, max) - mx + apply(mx, 1, min)
  
  pDT <- data.table(data.frame(mx2), keep.rownames = TRUE)
  gg <- unique(RESULTS.AGG[Comparison %in% COMBINATIONS[[tx]]][hit == TRUE]$Gene)
  ggtern(pDT, aes_string(x=cx[1], y=cx[2], z=cx[3])) + 
    theme_bw() +
    #geom_point(shape=1, color="#33a02c") +
    geom_hex_tern() +
    scale.hexgradient +
    geom_text(data=pDT[rn %in% gg], aes(label=rn), color="black", size=2)
  ggsave(out("Ternary_", tx, ".pdf"),w=8,h=8)
}
