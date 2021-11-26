source("src/00_init.R")
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
ggsave(out("QC.pdf"),w=30,h=20)



# Aggregate data ----------------------------------------------------------
ann.agg <- copy(ann)
ann.agg[grepl("^LSK", System) & Population == "LSK", Population := System]
ann.agg[, id := paste(Genotype, Population, Library, sep="_")]
ann.agg[, group := paste(Genotype, Library, sep="_")]
m2 <- sapply(with(ann.agg, split(sample, id)), function(sx){
  apply(m[,sx, drop=F], 1, function(row){
    if(all(is.na(row))){
      NA
    } else {
      sum(row, na.rm=TRUE)
    }
  })
})
write.table(m2, quote = F, sep = ",", row.names = TRUE, col.names = TRUE, file = out("Matrix_aggregated.csv"))
write.tsv(ann.agg, out("Annotation_aggregated.tsv"))

