source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("INT_02_CellrangerStatistics/")


ff <- list.files(paste(Sys.getenv("DATA"), sep='/'), pattern = "metrics_summary.csv", recursive = TRUE)
ff <- setNames(lapply(ff, function(fx){fread(paste(Sys.getenv("DATA"), fx, sep='/'))}), gsub("outs/.+$", "", ff))

res <- do.call(rbind, c(ff,fill=TRUE))
res$Dataset <- gsub("\\/$", "", names(ff))

write.tsv(res, out("CellrangerStats.tsv"))
