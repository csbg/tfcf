source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("FIG_02_POOLED/")

# Load data ---------------------------------------------------------------
m <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix))
ann <- fread(PATHS$POOLED$DATA$annotation)
stopifnot(all(ann$sample == colnames(m)))



list.files(dirout_load("POOLED_10_01_IndividualAnalysis")(""))
