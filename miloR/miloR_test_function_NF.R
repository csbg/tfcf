# x=traj_milo
# design=~miloR_gene
# design.df=traj_design
# fdr.weighting = c("k-distance")
# #fdr.weighting = c("k-distance", "neighbour-distance", "max", "graph-overlap", "none")
# min.mean = 0
# model.contrasts = NULL
# robust = TRUE
# reduced.dim = "PCA"
# norm.method = "TMM"

testNhoods_NF <- function(x,design, design.df,fdr.weighting = c("k-distance", "neighbour-distance", "max", "graph-overlap", "none"), min.mean = 0,model.contrasts = NULL, robust = TRUE, reduced.dim = "PCA",norm.method = "TMM"){
  
  if (!is.null(model.contrasts)) stop("This adapted function does not work with model contrasts any more")
    
  if (is(design, "formula")) {
    model <- model.matrix(design, data = design.df)
    rownames(model) <- rownames(design.df)
  } else if (is(design, "matrix")) {
    model <- design
    if (nrow(model) != nrow(design.df)) {
      stop("Design matrix and model matrix are not the same dimensionality")
    }
    if (any(rownames(model) != rownames(design.df))) {
      warning("Design matrix and model matrix dimnames are not the same")
      check.names <- any(rownames(model) %in% rownames(design.df))
      if (isTRUE(check.names)) {
        rownames(model) <- rownames(design.df)
      } else {
        stop("Design matrix and model matrix rownames are not a subset")
      }
    }
  }
  if (!is(x, "Milo")) {
    stop("Unrecognised input type - must be of class Milo")
  } else if (miloR:::.check_empty(x, "nhoodCounts")) {
    stop("Neighbourhood counts missing - please run countCells first")
  }
  if (!any(norm.method %in% c("TMM", "logMS", "RLE"))) {
    stop("Normalisation method ", norm.method, " not recognised. Must be either TMM, RLE or logMS")
  }
  if (!reduced.dim %in% reducedDimNames(x)) {
    stop(reduced.dim, " is not found in reducedDimNames. Avaiable options are ", 
         paste(reducedDimNames(x), collapse = ","))
  }
  subset.counts <- FALSE
  if (ncol(nhoodCounts(x)) != nrow(model)) {
    if (all(rownames(model) %in% colnames(nhoodCounts(x)))) {
      message("Design matrix is a strict subset of the nhood counts")
      subset.counts <- TRUE
    }else {
      stop("Design matrix (", nrow(model), ") and nhood counts (", 
           ncol(nhoodCounts(x)), ") are not the same dimension")
    }
  }
  if (min.mean > 0) {
    if (isTRUE(subset.counts)) {
      keep.nh <- rowMeans(nhoodCounts(x)[, rownames(model)]) >= 
        min.mean
    }else {
      keep.nh <- rowMeans(nhoodCounts(x)) >= min.mean
    }
   } else {
    if (isTRUE(subset.counts)) {
      keep.nh <- rep(TRUE, nrow(nhoodCounts(x)[, rownames(model)]))
    }else {
      keep.nh <- rep(TRUE, nrow(nhoodCounts(x)))
    }
  }
  if (isTRUE(subset.counts)) {
    keep.samps <- intersect(rownames(model), colnames(nhoodCounts(x)[keep.nh, 
    ]))
  }else {
    keep.samps <- colnames(nhoodCounts(x)[keep.nh, ])
  }
  if (any(colnames(nhoodCounts(x)[keep.nh, keep.samps]) != 
          rownames(model)) & !any(colnames(nhoodCounts(x)[keep.nh, 
                                                          keep.samps]) %in% rownames(model))) {
    stop("Sample names in design matrix and nhood counts are not matched.\n             Set appropriate rownames in design matrix.")
  }else if (any(colnames(nhoodCounts(x)[keep.nh, keep.samps]) != 
               rownames(model)) & any(colnames(nhoodCounts(x)[keep.nh, 
                                                              keep.samps]) %in% rownames(model))) {
    warning("Sample names in design matrix and nhood counts are not matched. Reordering")
    model <- model[colnames(nhoodCounts(x)[keep.nh, keep.samps]), 
    ]
  }
  if (length(norm.method) > 1) {
    message("Using TMM normalisation")
    dge <- DGEList(counts = nhoodCounts(x)[keep.nh, keep.samps], 
                   lib.size = colSums(nhoodCounts(x)[keep.nh, keep.samps]))
    dge <- calcNormFactors(dge, method = "TMM")
  }else if (norm.method %in% c("TMM")) {
    message("Using TMM normalisation")
    dge <- DGEList(counts = nhoodCounts(x)[keep.nh, keep.samps], 
                   lib.size = colSums(nhoodCounts(x)[keep.nh, keep.samps]))
    dge <- calcNormFactors(dge, method = "TMM")
  }else if (norm.method %in% c("RLE")) {
    message("Using RLE normalisation")
    dge <- DGEList(counts = nhoodCounts(x)[keep.nh, keep.samps], 
                   lib.size = colSums(nhoodCounts(x)[keep.nh, keep.samps]))
    dge <- calcNormFactors(dge, method = "RLE")
  }else if (norm.method %in% c("logMS")) {
    message("Using logMS normalisation")
    dge <- DGEList(counts = nhoodCounts(x)[keep.nh, keep.samps], 
                   lib.size = colSums(nhoodCounts(x)[keep.nh, keep.samps]))
  }
  dge <- estimateDisp(dge, model)
  fit <- glmQLFit(dge, model, robust = robust)
  if (!is.null(model.contrasts)) {
    stop("This adapted function does not work with model contrasts any more")
    # mod.constrast <- makeContrasts(contrasts = model.contrasts, 
    #                                levels = model)
    # res <- as.data.frame(topTags(glmQLFTest(fit, contrast = mod.constrast), 
    #                              sort.by = "none", n = Inf))
  }else {
    return.list <- list()
    for(coefx in colnames(model)){
      #n.coef <- ncol(model)
      res <- as.data.frame(topTags(glmQLFTest(fit, coef = coefx), 
                                   sort.by = "none", n = Inf))
      
      res$Nhood <- as.numeric(rownames(res))
      message("Performing spatial FDR correction with", fdr.weighting[1], 
              " weighting")
      mod.spatialfdr <- graphSpatialFDR(x.nhoods = nhoods(x), graph = graph(x), 
                                        weighting = fdr.weighting, k = x@.k, pvalues = res[order(res$Nhood), 
                                        ]$PValue, indices = nhoodIndex(x), distances = nhoodDistances(x), 
                                        reduced.dimensions = reducedDim(x, reduced.dim))
      res$SpatialFDR[order(res$Nhood)] <- mod.spatialfdr
      return.list[[coefx]] <- res
    }
    return(bind_rows(return.list, .id="coef"))
  }
}

