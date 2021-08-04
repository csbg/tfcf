source(paste0(Sys.getenv("CODE"), "src/00_init.R"))
out <- dirout("PPI_02_FirstPlots")

(load(dirout_load("PPI_01_CollectData")("GOI.Interactions.RData")))


interactions[db == "DepMap", dataset := gsub("^(.)", "\\U\\1", gsub("_", " ", dataset), perl=TRUE)]
interactions[,.N, by=c("db", "dataset", "organism")]



# DepMap analysis ---------------------------------------------------------
dpm <- toMT(interactions[db == "DepMap"], row = "pair", col = "dataset", val = "Score")
dpm[1,]

# Define groups of pairs --------------------------------------------------
pair.groups <- list(
  complex.extension = unique(interactions[!db %in% c("Corum", "DepMap")][!pair %in% unique(interactions[db %in% c("Corum")]$pair)][Score > 0.7]$pair),
  complexes = unique(interactions[db %in% c("Corum")]$pair),
  depmap.global = "",
  depmap.blood = ""
)



# # Identify interesting pairs ----------------------------------------------
# mds <- fread(dirout_load("POOLED_10_IndividualAnalysis")("Correlation_","hits","_MDS.tsv"))
# mds.dist <- as.matrix(mds[,2:3])
# row.names(mds.dist) <- mds$rn
# mds.dist <- as.matrix(dist(mds.dist))
# diag(mds.dist) <- max(mds.dist)
# ii <- unique(data.table(t(apply(which(mds.dist <= sort(mds.dist[upper.tri(mds.dist)], decreasing = FALSE)[200], arr.ind = TRUE), 1, sort))))
# pairs <- data.table(A=mds$rn[ii$V1], B=mds$rn[ii$V2])
# gg <- unique(c(pairs$A, pairs$B))
# 
# all.interactions.small <- all.interactions[Mouse_A %in% gg & Mouse_B %in% gg]
# 
# # Plots pairs -------------------------------------------------------------
# res <- data.table()
# for(pi in 1:nrow(pairs)){
#   pa=pairs[pi]$A
#   pb=pairs[pi]$B
#   res <- rbind(res, data.table(all.interactions.small[(Mouse_A == pa & Mouse_B == pb) | (Mouse_A == pb & Mouse_B == pa)], pair=paste(pa, pb)))
# }
# res <- res[!is.na(dataset)]
# 
# p.ppi <- ggplot(res, aes(y=pair, x=paste0(dataset, " (", organism, ") "), color=Score)) +
#   geom_point() +
#   theme_bw(12) + 
#   facet_grid(. ~ db, switch = "x", space = "free", scales = "free") +
#   scale_color_gradient(low="blue", high="red", limits=c(0,1))
# ggsave(out("PPI.result.pdf"), w=6,  h=15, plot=p.ppi + xRot())
# 
# 
# 
# p <- gridExtra::grid.arrange(
#   p.ppi,
#   p.dm, 
#   nrow=1, ncol=2, widths=c(3,4))
# ggsave(out("Results_combined.pdf"), h=10, w=9, plot=p)
