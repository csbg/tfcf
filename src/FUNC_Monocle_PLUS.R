
DotPlotData <- function(cds, markers, cols, pseudocount = 1, scale_max = 3, scale_min = -3, lower_threshold = 0){
  exprs_mat <- t(as.matrix(exprs(cds)[markers, ]))
  exprs_mat <- reshape2::melt(exprs_mat)
  colnames(exprs_mat) <- c("Cell", "Gene", "Expression")
  exprs_mat$Gene <- as.character(exprs_mat$Gene)
  
  exprs_mat <- data.table(cbind(exprs_mat, data.frame(colData(cds)@listData)[as.character(exprs_mat$Cell),]))
  for(colx in cols){
    exprs_mat <- exprs_mat[!is.na(get(colx))]
  }
  
  exprs_mat <- exprs_mat[, .(
    mean = mean(log(Expression + pseudocount)), 
    percentage = sum(Expression > lower_threshold)/length(Expression)*100,
    N=.N
  ), by=c(cols, "Gene")]
  exprs_mat[mean < scale_min, mean := scale_min]
  exprs_mat[mean > scale_max, mean := scale_max]
  return(exprs_mat)
}


getCL <- function(obj){
  as.character(monocle3::clusters(obj))
}

NF_TPM_Matrix <- function(cds, genes){
  cds_exprs <- SingleCellExperiment::counts(cds)[genes, , drop = FALSE]
  Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))
}


plot_cells_umap_hex_NF <- function(cds, genes, ncol=5, nbins=50, scale=FALSE){
  genes <- genes[genes %in% row.names(cds)]
  #genes <- sample(genes, 5)
  cds_exprs <- NF_TPM_Matrix(cds, genes)
  # cds_exprs <- SingleCellExperiment::counts(cds)[genes, , drop = FALSE]
  # cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))
  pDT <- data.table()
  for(gx in genes){
    pDT <- rbind(pDT, data.table(reducedDims(cds)[["UMAP"]], Expression=cds_exprs[gx,], Gene=gx))
  }
  func_summary <- function(x){
    log(mean(x)+1)
  }
  pDT[,scaleE := scale(Expression), by="Gene"]
  pDT[,scaleE := min(3, abs(scaleE)) * sign(scaleE)]
  
  p <- ggplot(pDT, aes(x=V1, y=V2)) +
    scale_fill_hexbin() +
    facet_wrap(~Gene, ncol=ncol) +
    theme_bw(12)
  if(scale){
    p <- p + stat_summary_hex(aes(z=scaleE),fun=mean, bins=nbins)
  } else {
    p <- p + stat_summary_hex(aes(z=Expression),fun=func_summary, bins=nbins)
  }
  return(p)
}