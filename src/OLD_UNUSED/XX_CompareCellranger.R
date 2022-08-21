source("src/00_init.R")

require(plyr)

out <- dirout("XX_CompareCellranger")


ff <- list.files(PATHS$LOCATIONS$DATA, pattern="metrics_summary.csv", recursive = TRUE, full.names = TRUE)
ff <- c(ff, list.files(gsub("Data", "Data2", PATHS$LOCATIONS$DATA), pattern="metrics_summary.csv", recursive = TRUE, full.names = TRUE))
ff <- ff[!grepl("OLD", ff)]
res <- rbindlist(lapply(ff, fread), fill=TRUE)
res$sample <- ff

res <- cbind(res, do.call(rbind, strsplit(res$sample, "/"))[,c(6,8)])

colnames(res)

pDT <- melt(res[,c("Estimated Number of Cells", "Mean Reads per Cell", "Fraction Reads in Cells", "V1", "V2")], id.vars = c("V1", "V2"))
pDT[, value := gsub("\\%", "", gsub(",", "", value))]
pDT[, value := as.numeric(value)]
pDT <- dcast.data.table(pDT, V2 + variable ~ V1, value.var = "value")
pDT <- pDT[!is.na(Data) & !is.na(Data2)]

ggplot(pDT, aes(x=Data, y=Data2, color=V2, shape=V2)) + 
  theme_bw(12) +
  geom_abline(color="lightgrey") +
  geom_point() +
  facet_wrap(~variable, scales = "free") +
  scale_shape_manual(values=rep(c(1,16,2,18,3,4), 20)) +
  xlab("Nik") + ylab("Julen")
ggsave(out("CellrangerCompare.pdf"),w=15,h=5)
