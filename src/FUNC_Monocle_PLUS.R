
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



