source("src/00_init.R")
out <- dirout("FULLINT_06_01_CytoTRACE/")

(load(PATHS$FULLINT$Monocle))

# INSTALLATION ------------------------------------------------------------
# Acutally used
# download.file("https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz", destfile = out("CytoTRACE_0.3.3.tar.gz"))
# renv::install("bioc::sva")
# renv::install("ccaPP")
# renv::install(out("CytoTRACE_0.3.3.tar.gz"), type = "source")

require(sva)
require(ccaPP)
require(CytoTRACE)
require(scran)

monocle.obj.original <- monocle.obj

sx <- monocle.obj$sample[1]
cytoResDT <- data.table()
for(sx in unique(monocle.obj.original$sample)){
  print("-------")
  print(sx)
  monocle.obj <- monocle.obj.original[,monocle.obj.original$sample == sx]
  
  # Get most variable genes
  vg <- monocle.obj %>%
    scran::modelGeneVar(assay.type = "counts") %>% 
    scran::getTopHVGs(n = 10000)
  
  # Calculate rowsums
  counts <- counts(monocle.obj)
  rowsumsMT <- Matrix::rowSums(counts)
  
  # Get top 10 000 genes (most variable first, then by rowsums)
  gg <- unique(c(vg, names(sort(rowsumsMT, decreasing = TRUE))))[1:1e4]
  
  counts <- counts[gg,]
  counts <- counts[Matrix::rowSums(counts) > 20,]
  # str(as.matrix(counts))
  #batch <- colData(monocle.obj)[,"sample"]
  #names(batch) <- row.names(colData(monocle.obj))
  data <- as.matrix(counts)
  #stopifnot(all(!is.na(batch)))
  #stopifnot(ncol(data) == length(batch))
  #stopifnot(all(colnames(data) == names(batch)))
  stopifnot(!any(duplicated(row.names(data))))
  stopifnot(!any(duplicated(colnames(data))))
  #stopifnot(!any(duplicated(names(batch))))
  
  cytoRes <- CytoTRACE(mat = data, ncores = 10)
  cytoRes$exprMatrix <- NULL
  
  ret <- data.table(rn=names(cytoRes$CytoTRACE))
  cytoGenes <- list()
  for(mx in names(cytoRes)){
    if(length(cytoRes[[mx]]) != nrow(ret)){
      # cytoGenes[[sx]] <- cytoRes$cytoGenes
      next
    }
    print(mx)
    stopifnot(names(cytoRes[[mx]]) == ret$rn)
    ret[[mx]] <- cytoRes[[mx]]
  }
  cytoResDT <- rbind(cytoResDT, ret)

}
save(cytoResDT, file=out("CytoTRACE.RData"))


