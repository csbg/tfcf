ref.umap.predict <- function (ref, query, ref.umap, reduction = "pca", ndim = 10, k = 20) {
  require(Seurat)
  stopifnot(all(colnames(ref) == row.names(ref.umap)))
  tdim <- dim(ref@reductions[[reduction]]@cell.embeddings)[2]
  if (ndim > tdim) {
    warning(sprintf("Number of dimensions ndim=%i is larger than the dimensions in reduction %s - Using only first %i dimensions", 
                    ndim, reduction, tdim))
    ndim = tdim
  }
  ref.space <- ref@reductions[[reduction]]@cell.embeddings[, 1:ndim]
  query.space <- query@reductions[[reduction]]@cell.embeddings[, 1:ndim]
  nn.method = "rann"
  nn.ranked <- Seurat:::NNHelper(data = ref.space, query = query.space, k = k, method = nn.method)
  
  umap.proj <- list()
  for (r in 1:dim(nn.ranked@nn.idx)[1]) {
    top.k <- nn.ranked@nn.idx[r, ]
    umap.proj[[r]] <- apply(ref.umap[top.k,], 2, mean)
  }
  umap.proj <- do.call(rbind, umap.proj)
  row.names(umap.proj) <- row.names(query.space)
  colnames(umap.proj) <- c("UMAP_1", "UMAP_2")
  return(umap.proj)
}



ref.trajectory.predict <- function (ref, query, reduction = "pca", ndim = 10, k = 20) {
  require(Seurat)
  tdim <- dim(ref@reductions[[reduction]]@cell.embeddings)[2]
  if (ndim > tdim) {
    warning(sprintf("Number of dimensions ndim=%i is larger than the dimensions in reduction %s - Using only first %i dimensions", 
                    ndim, reduction, tdim))
    ndim = tdim
  }
  ref.space <- ref@reductions[[reduction]]@cell.embeddings[, 1:ndim]
  query.space <- query@reductions[[reduction]]@cell.embeddings[, 1:ndim]
  nn.method = "rann"
  nn.ranked <- Seurat:::NNHelper(data = ref.space, query = query.space, k = k, method = nn.method)
  
  umap.proj <- list()
  for (r in 1:dim(nn.ranked@nn.idx)[1]) {
    top.k <- nn.ranked@nn.idx[r, ]
    cts <- ref$pseudotime.celltype[top.k]
    scores <- sort(table(cts), decreasing=TRUE)/k
    umap.proj[[r]] <- data.table(ct=names(scores)[1], conf=scores[1], pseudotime=mean(ref$pseudotime[top.k][cts == names(scores)[1]]))
  }
  umap.proj <- do.call(rbind, umap.proj)
  umap.proj$rn <- row.names(query.space)
  return(umap.proj)
}