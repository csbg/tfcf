source("src/00_init.R")
out <- dirout("SCRNA_99_02_TestDEG_AlmostAllGenes/")

require(Seurat)
require(tidyverse)
require(limma)
require(edgeR)


# Load data ---------------------------------------------------------------
inFile <- out("TCR_Data.RDS")
if(file.exists(inFile)){
  sobj <- readRDS(inFile)
} else{
  files <- c("TCR_ClassLabels.csv", "TCR_Data.h5")
  for(fx in files){
    download.file(paste0("http://kpnn.computational-epigenetics.org/", fx), fx)
  }
  x <- Seurat::Read10X_h5("TCR_Data.h5")
  labels <- read_csv("TCR_ClassLabels.csv")
  x <- x[,labels$barcode]
  saveRDS(x, inFile)
  sobj <- x
}



# Make monocle object -----------------------------------------------------
mobj <- new_cell_data_set(expression_data = sobj)
colData(mobj)$sample <- gsub("^[ACTG]+\\-", "s", colnames(mobj))



# Differential expression -------------------------------------------------
steps <- 10
(nrs <- floor(10**(seq(1, log10(ncol(mobj)), length.out = steps))))
i <- 1000
j <- 1
for(i in nrs){
  message("Run ",i)
  
  # 10 rounds
  for(j in 1:10){
    print(paste("replicate", j))
    
    if(i == ncol(mobj) & j > 1) next
    
    label <- paste0(i, "r", j)
    
    # cells to use
    cells <- sample(colnames(mobj), i)
    
    # Model matrix
    mm <- model.matrix(~colData(mobj)[cells, "sample"])
    colnames(mm)[2] <- "sample"
    
    # dataset
    dx <- counts(mobj)[, cells]
    gfilter <- rowSums(dx) > 0
    #gfilter <- filterByExpr(dx, design=mm, min.prop = 0.1, min.count=1)
    table(gfilter)
    
    # Limma voom
    print("Limma")
    s <- Sys.time()
    vObj <- voom(dx[gfilter,], design = mm, plot=TRUE)
    fit <- lmFit(vObj, design=mm)
    fit <- eBayes(fit)
    vRes <- data.table(topTable(fit, coef="sample", number = Inf), keep.rownames = TRUE)
    print(s - Sys.time())
    write.tsv(vRes, out("Limma_", label, ".tsv"))
    
    # Monocle
    print("Monocle")
    s <- Sys.time()
    mfit <- fit_models(mobj[gfilter,cells], model_formula_str = "~sample")
    mRes <- data.table(coefficient_table(mfit) %>% filter(term != "(Intercept)"), keep.rownames = TRUE)
    print(s - Sys.time())
    mRes$model <- NULL
    mRes$model_summary <- NULL
    write.tsv(mRes, out("Monocle_", label, ".tsv"))
    
    # Venn diagram of tested genes
    pdf(out("Venn_", label, ".pdf"))
    gplots::venn(list(monocle=mRes$gene_id, limma=vRes$rn))
    dev.off()
    
    # Match results
    mRes <- mRes[match(vRes$rn, gene_id)]
    
    # Export result
    x <- rbind(data.table(
      monocle=mRes$estimate,
      limma=vRes$logFC,
      type="estimate"),
      data.table(
        monocle=mRes$normalized_effect,
        limma=vRes$logFC,
        type="norm_effect"
      ),
      data.table(
        monocle=-log10(mRes$q_value) * sign(mRes$estimate),
        limma=-log10(vRes$adj.P.Val) * sign(vRes$logFC),
        type="padj"
      ))
    
    x$cells <- i
    x$replicate <- j
    
    write.tsv(x, out("Result_", label, ".tsv"))
  }
}

# Plot --------------------------------------------------------------------

# Load in data
ff <- list.files(out(""), pattern="Result.*.tsv", full.names = TRUE)
ff <- lapply(ff, fread)
ff <- rbindlist(ff)

# Plot correlations
tx <- "estimate"
for(tx in unique(ff$type)){
  ggplot(ff[type == tx], aes(x=monocle, y=limma)) +
    geom_hline(yintercept = 0, color="grey") +
    geom_vline(xintercept = 0, color="grey") +
    geom_abline(color="grey") +
    geom_hex() +
    facet_grid(factor(cells)~replicate, scales = "free") +
    theme_bw() +
    ggtitle(tx)
  ggsave(out("Hexbin_", tx, ".pdf"), w=20,h=20)
}

# Plot correlation
ff[,.(cor = cor(monocle, limma, m="s")), by=c("cells", "replicate", "type")] %>%
  ggplot(aes(x=factor(as.numeric(replicate)),y=factor(as.numeric(cells)), fill=cor)) + 
  geom_tile() + 
  facet_grid(. ~ type) + 
  scale_fill_gradient(low="white", high="red") +
  theme_bw()
ggsave(out("Hexbin_Correlations", ".pdf"), w=10,h=5)

# Check whether the discrepency is in the lowly expressed genes?
ff <- list.files(out(""), pattern="Monocle_\\d+r1.tsv")
xx <- gsub("Monocle_(\\d+)r1.tsv", "\\1", ff)
m <- lapply(xx, function(x){
  ret <- fread(out("Monocle_",x,"r1.tsv"))
  ret$cells <- x
  ret
})
m <- rbindlist(m)
v <- lapply(xx, function(x){
  ret <- fread(out("Limma_",x,"r1.tsv"))
  ret$cells <- x
  ret
})
v <- rbindlist(v)
m <- m[order(cells, gene_id)]
v <- v[order(cells, rn)]
stopifnot(all(m$gene_id == v$rn))
stopifnot(all(m$cells == v$cells))

m$diff <- v$logFC - m$normalized_effect
ggplot(m, aes(x=num_cells_expressed, y=diff)) + 
  geom_hex() +
  theme_bw() +
  facet_grid(status ~ factor(as.numeric(cells)), scales = "free_x") +
  scale_x_log10() +
  ylab("Difference limma logFC vs monocle norm effect")
ggsave(out("NrGenes_Differences.pdf"), w=20,h=6)  


ggplot(m, aes(x=num_cells_expressed, y=-log10(q_value))) + 
  geom_point(alpha=0.2) +
  theme_bw() +
  facet_grid(status ~ factor(as.numeric(cells)), scales = "free_x") +
  scale_x_log10()
ggsave(out("NrGenes_padj.jpg"), w=20,h=6)
