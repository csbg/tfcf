source("src/00_init.R")
base.dir <- "FIG_80_STables/"
out <- dirout(base.dir)

require(WriteXLS)


# Table single-cell stats -------------------------------------------------
ff <- list.files(paste0(PATHS$LOCATIONS$DATA), pattern="metrics_summary.csv", recursive = TRUE, full.names = TRUE)
names(ff) <- basename(gsub("\\/outs.+$", "", ff))
pDT <- rbindlist(lapply(ff, fread), use.names = TRUE, fill=TRUE, idcol="sample")
pDT <- merge(SANN[, c("sample_new", "tissue", "markers","timepoint", "sample_found")], pDT, by.x="sample_found", by.y="sample")[, -"sample_found",with=F]
names(pDT) <- gsub("^(.)", "\\U\\1", names(pDT), perl=TRUE)
WriteXLS(x=pDT, ExcelFileName=out("Supplementary_Table_SingleCellStats.xls"), AdjWidth=TRUE, BoldHeaderRow=TRUE, FreezeRow=1, SheetNames="Table")
write.tsv(pDT, out("Supplementary_Table_SingleCellStats.tsv"))
