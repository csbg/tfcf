source("src/00_init.R")
out <- dirout("FULLINT_06_01_CytoTRACE/")

(load(PATHS$FULLINT$Monocle))

# INSTALLATION ------------------------------------------------------------
# Acutally used
# download.file("https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz", destfile = out("CytoTRACE_0.3.3.tar.gz"))
# renv::install("bioc::sva")
# renv::install("ccaPP")
# renv::install(out("CytoTRACE_0.3.3.tar.gz"), type = "source")

require(CytoTRACE)

counts <- counts(monocle.obj)
counts <- counts[Matrix::rowSums(counts) > 20,]
str(as.matrix(counts))
cytoRes <- CytoTRACE(as.matrix(counts), batch = colData(monocle.obj)[,"sample"], ncores = 10)
cytoRes$exprMatrix <- NULL
save(cytoRes, file=out("CytoTRACE.RData"))
