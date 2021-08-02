source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
baseDir <- "POOLED_09_CleanData/"
out <- dirout(baseDir)

require(ineq)

m <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix))
ann <- fread(PATHS$POOLED$DATA$annotation)
stopifnot(all(ann$sample == colnames(m)))

str(m)
str(ann)

qc <- copy(ann)
qc$reads <- colSums(m, na.rm=TRUE)

cMT <- corS(m, use="pairwise.complete.obs")
diag(cMT) <-NA

i <- 16
for(i in 1:nrow(qc)){
  c.rep <- mean(cMT[qc[i]$sample, unique(merge(qc[i], qc, by=c("Genotype", "Population", "Library", "System"))$sample.y)], na.rm=TRUE)
  c.pop <- mean(cMT[qc[i]$sample, unique(merge(qc[i], qc, by=c("Genotype", "Population", "System"))$sample.y)], na.rm=TRUE)
  c.all <- mean(cMT[qc[i]$sample, ], na.rm=TRUE)
  qc[i, cor.rep := c.rep]
  qc[i, cor.pop := c.pop]
  qc[i, cor.all := c.all]
}

qc$cor.rep
qc$cor.pop
qc$cor.all

?Gini
qc$Gini <- apply(m, 2, Gini, corr=T)

pDT <- melt.data.table(qc, measure.vars = c("Gini", "reads", grep("^cor.rep", colnames(qc), value=TRUE)), id.vars=c("sample", "Genotype", "Population", "Library", "System"))
ggplot(pDT, aes(x=sample, y=value)) + geom_bar(stat="identity") + 
  facet_grid(variable ~ Genotype + Population, space = "free_x", scales = "free") +
  theme_bw(12) +
  xRot()
ggsave(out("QC.pdf"),w=20,h=20)


