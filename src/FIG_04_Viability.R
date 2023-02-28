source("src/00_init.R")
base.dir <- "FIG_04_Viability/"
out <- dirout(base.dir)

require(readxl)
require(depmap)
require(ggrepel)

# LOAD DATA ---------------------------------------------------------------

# Comparison to wildtype --------------------------------------------------
inDir <- dirout_load("POOLED_10_03_IndividualAnalysis_NormFactors_Controls")
RESULTS.wt.agg <- fread(inDir("ComparisonToWT.tsv"))
RESULTS.wt.agg.clean <- RESULTS.wt.agg %>%
  filter(!grepl("^GMP\\.", Population) & !Population %in% c("Mye", "Und"))

# Ainhoa's table ----------------------------------------------
ainhoa.cnts <- readxl::read_xlsx("metadata/Viability_Ainhoa_from_scRNA.xlsx")
ainhoa.cnts <- data.table(ainhoa.cnts)
ainhoa.cnts[, sgRNA_ID := toupper(gsub("NonTargetingControlGuideForMouse", "NTC", sgRNA_ID))]

# scData ------------------------------------------------------------------
scRNA.file <- out("scRNA.RDS")
if(file.exists(scRNA.file)){
  scRNA.cnts <- readRDS(scRNA.file)
} else {
  ff <- list.files(dirout_load("SCRNA_20_Summary")(""))
  (TISSUES <- gsub("_.+", "", ff[grepl("_monocle.singleR$", ff)]))
  inDir.funcs <- list()
  for(tx in TISSUES){inDir.funcs[[tx]] <- dirout_load(paste0("SCRNA_20_Summary/", tx, "_monocle.singleR"))}
  annList <- lapply(names(inDir.funcs), function(inDir.current){
    ann <- fread(inDir.funcs[[inDir.current]]("Annotation.tsv"))
    ann[, gene := gsub("_.+", "", guide)]
    ann
  })
  annList <- rbindlist(annList, fill=TRUE)
  scRNA.cnts <- annList[, length(unique(rn)), by=c("CRISPR_Cellranger", "tissue", "timepoint", "sample", "sample_broad")]
  scRNA.cnts <- scRNA.cnts[!is.na(CRISPR_Cellranger)]
  saveRDS(scRNA.cnts, scRNA.file)
}
scRNA.cnts <- scRNA.cnts[, .(cnt=sum(V1)), by=c("CRISPR_Cellranger", "tissue", "timepoint")]
scRNA.cnts <- dcast.data.table(scRNA.cnts, CRISPR_Cellranger ~ paste("scRNA", tissue, timepoint), value.var = "cnt")
scRNA.cnts[, CRISPR_Cellranger := toupper(CRISPR_Cellranger)]


# Sabatini lab ------------------------------------------------------------
sabatini <- readxl::read_xlsx("metadata/Viability_Sabatini_AML.xlsx", skip=1)
sabatini <- data.table(sabatini)



# Which genes to show -----------------------------------------------------
(GENES <- unique(gsub("_.+$", "", scRNA.cnts$CRISPR_Cellranger)))



# Depmap ------------------------------------------------------------------
depmap.file <- out("Depmap.RDS")
if(file.exists(depmap.file)){
  depmap <- readRDS(depmap.file)
} else {
  crispr <-  depmap::depmap_crispr()
  depmap.metdata <- depmap::depmap_metadata()
  depmap.ann <- depmap.metdata %>% 
    filter(lineage_subtype == "AML") %>%
    select(depmap_id, cell_line_name)
  crispr <- crispr %>%
    filter(gene_name %in% GENES)
  depmap <- merge(crispr, depmap.ann, by="depmap_id")
  saveRDS(data.table(depmap), depmap.file)
}


# SETUP ENDS HERE ---------------------------------------------------------



# Our counts scRNA --------------------------------------------------------------

# # Plot number of cells
# pDT <- copy(scRNA.cnts)
# pDT[, gene := gsub("_.+$", "", CRISPR_Cellranger)]
# pDT <- pDT[, .(V1 = sum(V1)), by=c("gene", "timepoint", "tissue")]
# ggplot(pDT[gene != "NTC"], aes(x=gene, y=V1 + 1)) + 
#   geom_col() + 
#   facet_grid(tissue + timepoint ~ .) + 
#   themeNF(rotate = TRUE) +
#   scale_y_log10()
# ggsave(out("scRNA_Guide_counts.pdf"), w=15, h=8)

# calculate CPMs
our.cnts <- merge(ainhoa.cnts[,-"oPool_ID",with=F], scRNA.cnts, by.x="sgRNA_ID", by.y="CRISPR_Cellranger", all=TRUE)
cMT <- as.matrix(our.cnts[, -"sgRNA_ID"])
row.names(cMT) <- our.cnts$sgRNA_ID
cMT[is.na(cMT)] <- 0
fMT <- t(t(cMT) / colSums(cMT) * 1e6)
stopifnot(all(colSums(fMT) == 1e6))

# Plot of CPMs
pDT <- melt(data.table(fMT, keep.rownames = TRUE), id.vars = "rn")
pDT[, gene := gsub("_.+$", "", rn)]
pDT[, analysis := ifelse(grepl("^scRNA", variable), "scRNA", "Ainhoa")]
pDT <- pDT[variable %in% pDT[,sum(value != 0), by="variable"][V1 > 1]$variable]
ggplot(pDT, aes(x=variable, y=rn, fill=value + 1)) + 
  geom_tile() + 
  scale_fill_gradientn(colours = c("white", "yellow", "red", "blue"), trans="log10") + 
  facet_grid(gene ~ analysis, space="free", scale="free") + 
  theme_bw() +
  xRot()
ggsave(out("scRNA_counts.pdf"), w=8,h=40, limitsize = FALSE)


# logFCs 
# aggregate by gene
gMTx <- t(sapply(split(row.names(fMT), gsub("_.+$", "", row.names(fMT))), function(rows) colSums(fMT[rows,, drop=F])))
quantile(gMTx)
#for(typex in c("all", "cleaned")){
typex <- "cleaned"
gMT <- gMTx
if(typex == "cleaned") gMT[gMT == 0] <- NA
# Normalize by NTC cpm (should be the same across all)
gMT <- t(t(gMT) * (mean(gMT["NTC",]) / gMT["NTC",]))
# divide by d0 counts
gMT <- log2(gMT + 0.1) - log2(gMT[,"Read_counts_d0"] + 0.1)


# Calculate values
pDT <- melt(data.table(gMT, keep.rownames = TRUE), id.vars = "rn")
pDT[, analysis := ifelse(grepl("^scRNA", variable), "scRNA", "Ainhoa")]
pDT <- pDT[variable %in% pDT[,sum(value != 0, na.rm=TRUE), by="variable"][V1 > 1]$variable]
pDT[, value.cap := sign(value) * pmin(5, abs(value))]
write.tsv(pDT, out("scRNA_processed.tsv"))

# Plot
pDT <- pDT[analysis == "scRNA"]
pDT <- pDT[variable != "scRNA in.vivo 28d"]
pDT <- pDT[rn != "NTC"]
pDT <- pDT[!is.na(value)]
pDT$rn <- factor(pDT$rn, levels=pDT[variable == "scRNA in.vivo 14d"][order(value)]$rn)
p <- ggplot(pDT, aes(x=variable, y=rn)) + 
  scale_fill_gradient2(name=expression(log[2](FC)), low="#1f78b4", high="#e31a1c") + 
  #facet_grid(. ~ analysis, space="free", scale="free") + 
  labs(x="Population", y="CRISPR Target") +
  themeNF(rotate=TRUE)
ggsave(out("scRNA_counts_summarized_",typex,"_capped.pdf"), w=3,h=8, limitsize = FALSE,
       plot=p + geom_tile(aes(fill=value.cap)))
ggsave(out("scRNA_counts_summarized_",typex,"_raw.pdf"), w=3,h=8, limitsize = FALSE,
       plot=p + geom_tile(aes(fill=value)))

ORDER.GUIDES <- levels(pDT$rn)


# Our counts vs wildtype --------------------------------------------------------------
pDT <- RESULTS.wt.agg[toupper(Gene) %in% GENES]
ggplot(pDT, aes(x=Population, y=Gene, fill=log2FC)) + 
  geom_tile() + 
  scale_fill_gradient2(low="blue", high="red") + 
  theme_bw() +
  xRot()
ggsave(out("POOLED_counts.pdf"), w=5,h=8, limitsize = FALSE)


pDT <- RESULTS.wt.agg.clean[toupper(Gene) %in% ORDER.GUIDES]
pDT$Gene <- factor(toupper(pDT$Gene), levels=ORDER.GUIDES)
ggplot(pDT, aes(x=Population, y=Gene, fill=log2FC)) + 
    geom_tile() + 
    scale_fill_gradient2(name=expression(log[2](FC)), low="#1f78b4", high="#e31a1c") + 
    themeNF(rotate=TRUE) +
  labs(x="Population", y="CRISPR Target")
ggsave(out("POOLED_counts_cleaned.pdf"), w=3,h=8, limitsize = FALSE)


# Sabitini ----------------------------------------------------------------
pDT <- melt(sabatini[,-"sgRNAs included",with=F][Gene %in% toupper(GENES)], id="Gene")
ggplot(pDT, aes(x=variable, y=Gene, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="blue", high="red") + 
  theme_bw() +
  xRot()
ggsave(out("Sabatini_counts.pdf"), w=5,h=8, limitsize = FALSE)


# Depmap ------------------------------------------------------------------
pDT <- depmap %>% select(gene_name, dependency, cell_line_name)
pDT <- data.table(pDT)[gene_name %in% ORDER.GUIDES]
pDT$gene_name <- factor(toupper(pDT$gene_name), levels=ORDER.GUIDES)
ggplot(pDT, aes(x=cell_line_name, y=gene_name, fill=dependency)) + 
  geom_tile() + 
  scale_fill_gradient2(name="DepMap\nDependency", low="#1f78b4", high="#e31a1c") + 
  themeNF(rotate=TRUE) +
  labs(x="Cell line", y="CRISPR Target")
ggsave(out("Depmap_logFC.pdf"), w=5,h=8, limitsize = FALSE)



# Summarized form ---------------------------------------------------------
x1 <- depmap %>%
  mutate(Gene = gene_name) %>%
  group_by(Gene) %>%
  summarize(depmap=min(dependency))

x2 <- RESULTS.wt.agg.clean %>%
  mutate(Gene = toupper(Gene)) %>%
  group_by(Gene) %>%
  summarize(pooled=min(log2FC))

pDT <- inner_join(x1, x2, by="Gene")
cor(pDT$depmap, pDT$pooled)
cor(pDT$depmap, pDT$pooled, method="spearman")


pDT2 <- data.table(pDT)[order(-abs(depmap - pooled))][1:10]
ggplot(pDT, aes(x=pooled, y=depmap)) + 
  theme_bw() +
  geom_hline(yintercept = 0, color="grey") +
  geom_vline(xintercept = 0, color="grey") +
  geom_abline(color="#1f78b4") +
  geom_point() +
  geom_point(data=pDT2, color="red") +
  geom_text_repel(data=pDT2, aes(label=Gene), color="red") +
  labs(x="log fold change (this study)", y="DepMap dependency")
ggsave(out("Diff_PooledVsDepMap.pdf"),w=5,h=4)  
  
