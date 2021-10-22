source("src/00_init.R")
out <- dirout("FIG_03_POOLED_Leukemia/")

require(latex2exp)
require(ggrepel)
require(igraph)
require(ggtext)

# Load data ---------------------------------------------------------------
m <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix))
ann <- fread(PATHS$POOLED$DATA$annotation)
stopifnot(all(ann$sample == colnames(m)))


# Load analysis results ---------------------------------------------------
RESULTS.ctr.pval <- fread(dirout_load("POOLED_12_01_Controls")("Results_Pvalues.tsv"))



# Define comparisons that are to be combined ------------------------------
COMBINATIONS <- list(
  healthy=list(m="GMP.LSK", e="MEP.LSK", i="GMP.MEP"),
  #h2=list(m="", e="", i="MYE.UND"),
  leukemia=list(m="DMMye.LSC", e="DMEry.LSC", i="DMMye.Ery")
)


# SETUP ENDS HERE ---------------------------------------------------------




# Analysis per guide -------------------------------------------------------
xDT_r <- dcast.data.table(
  data = RESULTS.ctr.pval[Comparison %in% do.call(c, lapply(COMBINATIONS, unlist))], 
  formula = Group + Guide + Gene ~ Comparison, 
  value.var = "Ratio.log2")

xDT_p <- dcast.data.table(
  data = RESULTS.ctr.pval[Comparison %in% do.call(c, lapply(COMBINATIONS, unlist))], 
  formula = Group + Guide + Gene ~ Comparison, 
  value.var = "padj")

id.cols <- setdiff(colnames(xDT_r), unique(RESULTS.ctr.pval$Comparison))

xDT <- merge(xDT_r, xDT_p, by=id.cols, all=TRUE, suffixes=c("_r", "_p"))

stopifnot(all(duplicated(with(xDT, paste(Group, Guide)))==FALSE))

ggplot(xDT, aes(x=DMEry.LSC, y=DMMye.LSC, color=DMMye.Ery)) + 
  theme_bw() +
  geom_point() + 
  scale_color_gradient2()
ggsave("ContraTriangle1.pdf", w=5,h=4)

ggplot(xDT, aes(x=DMMye.LSC - DMEry.LSC, y=DMMye.Ery)) + 
  ggtitle(paste("R=", round(with(x, cor(DMMye.LSC - DMEry.LSC, DMMye.Ery)),3))) +
  theme_bw() +
  geom_point()
ggsave("ContraTriangle2.pdf", w=5,h=4)


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
}



lDT <- xDT[
  sign(GMP.LSK_r) != sign(DMMye.LSC_r) & (GMP.LSK_p < 0.05 | DMMye.LSC_p < 0.05) |
  sign(MEP.LSK_r) != sign(DMEry.LSC_r) & (MEP.LSK_p < 0.05 | DMEry.LSC_p < 0.05)
]
ggplot(xDT, aes(x=GMP.LSK_r, y=MEP.LSK_r)) + 
  theme_bw() +
  geom_point(color="grey", alpha=0.3) + 
  geom_abline() +
  geom_segment(data=lDT, aes(xend=DMMye.LSC_r, yend=DMEry.LSC_r)) +
  ggtitle(ci)
ggsave(out("Scatterplot_", ci, ".pdf"), w=5,h=4)




# Analysis per gene -------------------------------------------------------

RESULTS.AGG <- RESULTS.ctr.pval[!grepl("^WT_", Group),.(
  Ratio.log2 = mean(Ratio.log2), 
  N=.N, 
  Nsig.pos=sum(padj < 0.05 * sign(Ratio.log2 > 0)),
  Nsig.neg=sum(padj < 0.05 * sign(Ratio.log2 < 0))
  ), by=c("Comparison", "Gene")]
RESULTS.AGG[, Frag.sig := pmax(Nsig.pos, Nsig.neg)/N]
RESULTS.AGG[Nsig.pos > 0 & Nsig.neg > 0, Frag.sig := 0]


xDT_r <- dcast.data.table(
  data = RESULTS.AGG[Comparison %in% do.call(c, lapply(COMBINATIONS, unlist))], 
  formula = Gene ~ Comparison, 
  value.var = "Ratio.log2")

xDT_p <- dcast.data.table(
  data = RESULTS.AGG[Comparison %in% do.call(c, lapply(COMBINATIONS, unlist))], 
  formula = Gene ~ Comparison, 
  value.var = "Frag.sig")

id.cols <- setdiff(colnames(xDT_r), unique(RESULTS.ctr.pval$Comparison))

xDT <- merge(xDT_r, xDT_p, by=id.cols, all=TRUE, suffixes=c("_r", "_p"))

stopifnot(all(duplicated(xDT$Gene)==FALSE))

ggplot(xDT, aes(x=DMEry.LSC_r, y=DMMye.LSC_r, color=DMMye.Ery_r)) + 
  theme_bw() +
  geom_point() + 
  scale_color_gradient2()
ggsave("ContraTriangle1.pdf", w=5,h=4)

ggplot(xDT, aes(x=DMMye.LSC_r - DMEry.LSC_r, y=DMMye.Ery_r)) + 
  ggtitle(paste("R=", round(with(xDT, cor(DMMye.LSC_r - DMEry.LSC_r, DMMye.Ery_r, use="pairwise.complete.obs")),3))) +
  theme_bw() +
  geom_point()
ggsave("ContraTriangle2.pdf", w=5,h=4)


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



lDT <- xDT[
  sign(GMP.LSK_r) != sign(DMMye.LSC_r) & (GMP.LSK_p > 0.5 | DMMye.LSC_p > 0.5) |
    sign(MEP.LSK_r) != sign(DMEry.LSC_r) & (MEP.LSK_p > 0.5 | DMEry.LSC_p > 0.5)
]

# lDT.labels <- lDT[,.(m_r=pmean(DMMye.LSC_r, GMP.LSK_r), e_r=pmean(DMEry.LSC_r, MEP.LSK_r)), by="Gene"]

ggplot(pDT, aes(x=m_r, y=e_r, color=tissue)) + 
  theme_classic() +
  geom_abline(color="black", size=0.5, linetype=3) +
  geom_hline(yintercept = 0, color="black", size=0.5, linetype=3) +
  geom_vline(xintercept = 0, color="black", size=0.5, linetype=3) +
  geom_point(aes(size=pmax(m_p, e_p)*100, alpha=pmax(m_p, e_p))) + 
  geom_segment(
    color="#666666", size=0.3, data=lDT,  alpha=1,
    arrow=arrow(length = unit(0.25, "cm")),
    aes(x=GMP.LSK_r, y=MEP.LSK_r, xend=DMMye.LSC_r, yend=DMEry.LSC_r)) +
  geom_text_repel(data=pDT[Gene %in% lDT$Gene & pmax(m_p, e_p) > 0.5], aes(label=Gene), color="black") +
  scale_color_manual(values=c(leukemia="#fb9a99", healthy="#a6cee3")) +
  scale_size_continuous(limits=c(0,100), name="% sig. gRNAs") +
  #geom_label(data=pDT[Gene %in% lDT$Gene & pmax(m_p, e_p) > 0.5], aes(label=Gene), alpha=0.3) +
  # geom_text(data=lDT.labels,  aes(label=Gene), color="black") +
  xlab("Average log2 ratio: GMP/LSK (across guides)") +
  ylab("Average log2 ratio: MEP/LSK (across guides)") +
  guides(alpha=FALSE) +
  ggtitle("Comparison of KO effects between in vitro and leukemia")
ggsave(out("Scatterplot_comb.pdf"), w=10,h=8)
