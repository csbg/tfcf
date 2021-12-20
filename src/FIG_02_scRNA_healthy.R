source("src/00_init.R")
base.dir <- "FIG_02_scRNA_healthy/"
out <- dirout(base.dir)


# FUNCTIONS ---------------------------------------------------------------
ds <- function(path){load(path); return(monocle.obj)}



# Folders -----------------------------------------------------------------
list.files(dirout_load("")(""))
inDir.funcs <- list(
  "in.vivo"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vivo"),
  "in.vitro"=dirout_load("FULLINT_10_01_BasicAnalysis_in.vitro"),
  "leukemia"=dirout_load("FULLINT_10_01_BasicAnalysis_leukemia")
)
inDir.full <- dirout_load("FULLINT_10_01_BasicAnalysis_combined")



# Load data ---------------------------------------------------------------
marker.signatures.use <- fread("metadata/markers.signatures.use.scRNA.tsv")
marker.signatures.use[, sig := paste(DB, Celltype)]

# Load signatures
ff <- list.files(dirout_load("FULLINT_08_01_Markers")(""), pattern="Signatures_")
ff <- dirout_load("FULLINT_08_01_Markers")(ff)
names(ff) <- gsub("^Signatures_(.+?).csv$", "\\1", basename(ff))
ff <- ff[!grepl("Augmented_2021", ff)]
marker.signatures <- lapply(ff, function(fx) as.matrix(read.csv(fx)))
# Prep data
markS <- lapply(unique(marker.signatures.use$DB), function(x) data.table(melt(data.table(marker.signatures[[x]], keep.rownames = TRUE), id.vars = "rn"), DB=x))
markS <- merge(do.call(rbind, markS), marker.signatures.use, by.x=c("variable", "DB"), by.y=c("Celltype", "DB"))



# Specific tissue ---------------------------------------------------------
inDir.current <- "in.vivo"

#dsx <- ds(inDir.funcs[[inDir.current]]("MonocleObject.RData"))
ann <- fread(inDir.funcs[[inDir.current]]("Annotation.tsv"))
sda <- fread(inDir.funcs[[inDir.current]]("SigDA.tsv"))
sda <- merge(sda, unique(ann[,c("timepoint", "markers", "sample")]), c("sample"))

ann[,.N, by=c("sample", "markers", "timepoint")][order(markers, timepoint)]
ann[, perturbed := !(mixscape_class.global %in% c("NP", "NTC") | is.na(mixscape_class.global))]

# UMAP of timepoint / markers
ggplot(ann[perturbed == FALSE], aes(x=UMAP1, y=UMAP2)) + 
  theme_bw(12) +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~markers + timepoint, ncol=3)
ggsave(out("UMAP_Samples.pdf"), w=7,h=5)

# Plot top signatures
pDT <- merge(ann[perturbed == FALSE][,c("rn", "UMAP1", "UMAP2"),with=F], markS, by="rn")
pDT[, value.norm := scale(value), by="FinalName"]
ggplot(pDT, aes(x=UMAP1, y=UMAP2)) +
  stat_summary_hex(aes(z=pmin(3, value.norm)),fun=mean, bins=100) +
  scale_fill_gradient2(low="blue", midpoint = 0, high="red") +
  #scale_fill_hexbin() +
  theme_bw(12) +
  facet_wrap(~FinalName) +
  ggtitle(paste("Marker Signatures"))
ggsave(out("UMAP_MarkerSignatures.pdf"), w=12,h=12)

# Plot enrichment scores
pDT <- sda[, .(
  d=mean(d), 
  sig.pos=sum(padj < 0.05 & d > 0),
  sig.neg=sum(padj < 0.05 & d < 0),
  total = .N
  ), by=c("guide", "gene", "timepoint", "markers", "sig")]
pDT <- merge(pDT, marker.signatures.use, by="sig")
pDT[, sig.perc := ifelse(d > 0, sig.pos, sig.neg)/total * 100]
ggplot(pDT[timepoint == "14d" & markers == "Lin-"], aes(y=FinalName, x=guide, size=sig.perc, color=d)) + 
  theme_bw() +
  scale_color_gradient2(low="blue", midpoint = 0, high="red") +
  geom_point() +
  xRot() + 
  facet_grid(. ~ gene, scale="free", space="free")
ggsave(out("MarkerSignatures.pdf"), w=12,h=6)

# Plot top guides
pDT <- ann[timepoint == "14d" & markers == "Lin-"]
pDT2 <- pDT[mixscape_class.global == "NTC"]
for(x in c("Wdr82", "Rcor1")){
  pDT2 <- rbind(pDT2, pDT[grepl(paste0("^", x), guide) & mixscape_class.global == "KO"])
}
ggplot(pDT2, aes(x=UMAP1, y=UMAP2)) + 
  theme_bw() +
  geom_hex(bins=100) +
  scale_fill_gradient(low="lightgrey", high="blue") +
  facet_wrap(~guide + timepoint, ncol=3)
ggsave(out("UMAP_Guides.pdf"), w=9,h=10)





