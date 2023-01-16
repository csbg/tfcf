source("src/00_init.R")
base.dir <- "FIG_01_POOLED_Replicates/"
out <- dirout(base.dir)

require(latex2exp)
require(ggrepel)
require(igraph)
require(ggtext)
require(WriteXLS)
require(ComplexHeatmap)


# Load data ---------------------------------------------------------------
m <- as.matrix(read.csv(PATHS$POOLED$DATA$matrix))
ann <- fread(PATHS$POOLED$DATA$annotation)
stopifnot(all(ann$sample == colnames(m)))


# Gene annotations --------------------------------------------------------
ANN.genes <- fread("metadata/TFCF_Annotations_v2.tsv", check.names = TRUE)
ANN.genes[,Complex_simple := BROAD.COMPLEX]
ANN.genes <- unique(ANN.genes[!is.na(Complex_simple) & Complex_simple != ""][,c("GENE", "Complex_simple"), with=F])
stopifnot(all(ANN.genes[,.N, by="GENE"]$N == 1))
unique(ANN.genes$Complex_simple)


# ggpairs hexbin function -------------------------------------------------

ggpairs_hex <- function(df, hexbins = 10) {
  # REF: https://stackoverflow.com/questions/20872133/using-stat-binhex-with-ggpairs
  p <- ggpairs(df, lower="blank")
  seq <- 1:ncol(df)
  for (x in seq)
    for (y in seq) 
      if (y>x) 
        p <- putPlot(p, ggplot(df, aes_string(x=names(df)[x],y=names(df)[y])) + stat_binhex(bins = hexbins), y,x)
  
  return(p)
}


# SETUP ENDS HERE ---------------------------------------------------------




# Overview plot -----------------------------------------------------------
libx <- "As"
for(libx in unique(ann$Library)){
  annX <- ann[Library == libx]
  annX <- annX[Genotype == "Cas9" & !grepl("DM", sample)]
  annX <- annX[Population %in% annX[,.N, by="Population"][N>1]$Population]
  annX[Population == "LSK", Population := System]
  annX[, sx := as.numeric(factor(sample)), by=c("Population")]
  annX[, name_clean := paste0(Population, "-s", sx)]
  if(nrow(annX) <= 4) next
  cMT <- corS(m[,annX$sample], use="p")
  diag(cMT) <- NA
  hm <- Heatmap(cMT, 
                row_labels = annX$name_clean, 
                column_labels = annX$name_clean,
                column_split = annX$Population,
                row_split = annX$Population,
                col=RColorBrewer::brewer.pal(8, "OrRd"),
                heatmap_legend_param = list(title = "Spearman\ncorrelation"))
  pdf(out("Heatmap", libx, ".pdf"), 
      w=nrow(annX) * 0.2 + 2,
      h=nrow(annX) * 0.2 + 1)
  draw(hm)
  dev.off()
}


# Dedicated plot ----------------------------------------------------------
libx <- "Br"
for(libx in c("As", "Br")){
for(ctx in c("cKit", "GMP")){
  annX <- ann[Library == libx]
  annX[Population == "LSK", Population := System]
  annX <- annX[Genotype == "Cas9"]
  annX <- if(ctx == "cKit") annX[Population %in% c("LSKd9", "cKit")] else annX[Population %in% c("GMP", "MEP")]
  annX[Population == "LSKd9", Population := "LSK"]
  annX[, sx := as.numeric(factor(sample)), by=c("Population")]
  annX[, name_clean := paste0(Population, "_s", sx)]
  annX <- annX[order(name_clean)]
  mat <- m[,annX$sample]
  colnames(mat) <- annX$name_clean
  cleanDev()
  pm <- ggpairs(data.frame(mat), 
                lower = list(continuous = wrap("points", shape=1)),
                #upper = list(continuous = "density"),
                upper = list(continuous = wrap("cor", method = "spearman", color="black")),
                diag = list(continuous = "barDiag")) +
    scale_x_log10() +
    scale_y_log10() +
    themeNF(rotate = TRUE)
  
  pdf(out("Details", libx, "_", ctx, "_Pairs.pdf"), 
      w=nrow(annX) * 1 + 2,
      h=nrow(annX) * 1 + 1)
  print(pm)
  dev.off()
  
  cMT <- corS(mat, use="p")
  diag(cMT) <- NA
  hm <- Heatmap(cMT, 
                row_labels = annX$name_clean, 
                column_labels = annX$name_clean,
                column_split = annX$Population,
                row_split = annX$Population,
                col=RColorBrewer::brewer.pal(8, "OrRd"),
                heatmap_legend_param = list(title = "Spearman\ncorrelation"))
  pdf(out("Details", libx, "_", ctx, "_HM.pdf"), 
      w=nrow(annX) * 0.5 + 2,
      h=nrow(annX) * 0.5 + 1)
  draw(hm)
  dev.off()
}
}

