source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
baseDir <- "POOLED_10_IndividualAnalysis/"



# Load data ---------------------------------------------------------------
m <- as.matrix(read.csv(PATHS$DATA$matrix))
str(m)
ann <- fread(PATHS$DATA$annotation)
ann[,Date := paste0(Date, "_",Date2)]
stopifnot(all(ann$sample == colnames(m)))


res <- data.table()

datex <- "Nov2019"
libx <- "A"
for(datex in unique(ann$Date)){
  for(libx in unique(ann[Date == datex]$Library)){
    
    xAnn <- ann[Date == datex & Library == libx]
    xMT <- m[,xAnn$sample]
    #table(apply(!is.na(m[grepl("NonTargetingControl", row.names(m)),xAnn$sample]), 1, sum))
    xMT <- xMT[apply(!is.na(xMT), 1, sum) > ncol(xMT) * 0.8,]
    xMT <- voom(xMT, plot=FALSE)$E
    
    popx <- "Mye"
    for(popx in unique(xAnn$Population)){
      pAnn <- xAnn[Population == popx]
      stopifnot(nrow(pAnn) == 2)
      stopifnot(c("Cas9", "WT") %in% pAnn$Genotype)
      res <- rbind(res, data.table(
        Population = popx, 
        Score=xMT[,pAnn[Genotype == "Cas9"]$sample] - xMT[,pAnn[Genotype == "WT"]$sample], 
        Guide=row.names(xMT),
        Date = datex,
        Library=libx
      ))
    }
  }
}


ggplot(res, aes(x=Score, color=Population)) + 
  geom_density() +
  scale_color_manual(values=RColorBrewer::brewer.pal(name = "Dark2", n=length(unique(res$Population)))) +
  facet_wrap(~Library + Date, scales = "free")

ggplot(res, aes(x=Score, color=Population)) + 
  stat_ecdf() +
  scale_color_manual(values=RColorBrewer::brewer.pal(name = "Dark2", n=length(unique(res$Population)))) +
  facet_wrap(~Library + Date, scales = "free")
