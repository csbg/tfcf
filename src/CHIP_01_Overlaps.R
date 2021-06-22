source(paste0(Sys.getenv("CODE"), "src/00_init.R"))

require(GenomicRanges)



# Settings ----------------------------------------------------------------
out <- dirout("CHIP_01_Overlaps/")
cutoff <- -log10(0.01)


# Functions ---------------------------------------------------------------
convertSampleNames <- function(x){
  ret <- data.table(
    sample=x,
    group=gsub("^(.+?)_.+$", "\\1", x),
    number=gsub("^.+_(S\\d+)$", "\\1", x),
    factor=toupper(gsub("\\-.+$", "", gsub("\\d+K", "", gsub("^(.+?)_(.+?)_.+$", "\\2", x))))
  )
  ret[number==sample, number := "SX"]
  ret[factor==sample, number := "SX"]
  return(ret)
}


# Load data ---------------------------------------------------------------
ff <- list.files(paste0(Sys.getenv("DATA"), "ChIP_Peaks/"), full.names = TRUE, pattern=".narrowPeak")
names(ff) <- gsub("_peaks.narrowPeak", "", basename(ff))
dat <- lapply(ff, fread)
sigPeaks <- lapply(dat, function(dt){return(dt[V9 > cutoff])})
lapply(sigPeaks, nrow)
sapply(dat, nrow)
sum(sapply(dat, nrow))


# Plot number of peaks ----------------------------------------------------
pDT <- data.table(peaks=sapply(dat, nrow), sample=names(dat))
xx <- convertSampleNames(pDT$sample)
pDT <- merge(pDT, xx, by="sample")
ggplot(pDT, aes(x=sample, y=peaks, fill=group)) + geom_bar(stat="identity") + theme_bw() + 
  xRot() + 
  scale_y_log10() + 
  facet_grid(.~factor, scales = "free", space = "free")
ggsave(out("NumberPeaks.pdf"), w=14,h=6)




# merge peaks -------------------------------------------------------------
peaks <- lapply(dat, function(x){
  gr <- x[,1:3,with=F]
  colnames(gr) <- c("CHR","START","END")
  gr[, START := as.numeric(START)]
  gr[, END := as.numeric(END)]
  as(gr, "GRanges")
})

peaks <- peaks[!sapply(peaks, is.null)]
peaks.consensus <- GRangesList(peaks)
peaks.consensus <- unlist(peaks.consensus)
peaks.consensus <- reduce(peaks.consensus)
hist(width(peaks.consensus))
quantile(width(peaks.consensus))

ggplot(data.table(width(peaks.consensus)), aes(x=V1)) + stat_ecdf() + 
  scale_x_log10(limits=c(1,NA)) + theme_bw(12) + xlab("width")
ggsave(out("Peak_widths.pdf"),w=5, h=5)

olMT <- sapply(peaks, function(peak) countOverlaps(peaks.consensus, peak))

cMT <- cor(olMT)
diag(cMT) <- NA
dd <- as.dist(1-abs(cMT))

cleanDev(); pdf(out("Correlation.pdf"),w=12,h=12)
pheatmap(cMT, clustering_distance_rows = dd, clustering_distance_cols = dd)
dev.off()


i1 <- 1
i2 <- 2
res <- data.table()
for(i1 in 1:(ncol(olMT)-1)){
  message(i1)
  for(i2 in (i1+1):ncol(olMT)){
    cmt <- table(olMT[,i1] > 0, olMT[,i2] > 0)
    if(!all(dim(cmt) == c(2,2))) next
    fish <- fisher.test(cmt)
    res <- rbind(res, data.table(s1=colnames(olMT)[i1], s2=colnames(olMT)[i2], p=fish$p.value, OR=fish$estimate))
    res <- rbind(res, data.table(s2=colnames(olMT)[i1], s1=colnames(olMT)[i2], p=fish$p.value, OR=fish$estimate))
  }
}

res[, padj := p.adjust(p, method="BH")]
res[, OR_cap := pmin(OR, max(res[OR != Inf]$OR))]
quantile(res$OR)

write.tsv(res, out("Enrichments.tsv"))

res <- fread(out("Enrichments.tsv"))

cleanDev(); pdf(out("Enrichment_OR.pdf"),w=12,h=12)
pheatmap(log2(toMT(res, row="s1", col="s2", val = "OR_cap")))
dev.off()

# fish.ordering <- function(dt, toOrder, orderBy, value.var, aggregate = FALSE){
#   orMT <- as.matrix(dcast.data.table(dt, get(orderBy) ~ get(toOrder), value.var=value.var)[,-"orderBy",with=F])
#   orMT[is.na(orMT)] <- 1
#   orMT[orMT == Inf] <- max(orMT[orMT != Inf])
#   orMT <- log2(orMT)
#   
#   dd <- as.dist(max(orMT[orMT != Inf]) - abs(orMT))
#   pheatmap(orMT, clustering_distance_rows = dd, clustering_distance_cols = dd)
#   pheatmap(orMT)
#   
#   hclustObj <- hclust(as.dist(max(orMT[orMT != Inf]) - abs(orMT)))
#   dt[[toOrder]] <- factor(dt[[toOrder]], levels=hclustObj$labels[hclustObj$order])
#   return(dt)
# }
# 
# HM.COLORS.FUNC <- colorRampPalette(c("#6a3d9a", "#a6cee3", "white", "#fdbf6f", "#e31a1c"))
# pheatmap(orMT, breaks=seq(-10,10, 0.1), color=HM.COLORS.FUNC(201))
# 
# res <- hierarch.ordering(dt=res, toOrder="s1", orderBy="s2", value.var="OR_cap")
# res <- hierarch.ordering(res, "s2", "s1", "OR_cap")
# 
# 
# ggplot(res, aes(x=s1, y=s2, color=log2(OR_cap), size=pmin(5, -log10(padj)))) + 
#   geom_point() +
#   scale_color_gradient2(low="blue", high='red') +
#   theme_bw() + 
#   xRot()
# ggsave(out("Enrichment.pdf"), w=15,h=15)

       